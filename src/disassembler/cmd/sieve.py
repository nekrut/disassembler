"""sieve: thin near-duplicate haplotypes with the k-mer diversity sieve.

Wraps `extern/kmer-diversity-sieve/l1_pipeline_serial`. The sieve is a MinHash/LSH
lineage-level reducer: on defaults tuned to gene scale (tau=0.005, sketch=256) it
collapses a cohort of ~identical lineages while preserving topology; the sieve
binary's own defaults (tau=0.15, sketch=64) target archive-scale collapse and reduce
a gene-scale cohort to ~1 representative. See docs/sieve.md and
docs/notes/sieve_calibration.md.

Both the CLI and the composite pipelines call `sieve_run(...)`, which removes exact
duplicates first and reports what the binary's completeness gate will drop.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import signal
import subprocess
import tempfile
import time
import sys
from pathlib import Path
from typing import NamedTuple

from disassembler.resolve import resolve_sieve, not_found_message, EX_TOOL_NOT_FOUND
from disassembler.textio import open_text
from disassembler.validators import fasta_record_lengths, fraction, positive_int

from disassembler import __version__
from disassembler import summary as summary_io

# The binary's completeness gate (`--max-n-frac`) counts every character other than A/C/G/T,
# '-' included. Walker rows carry no gaps, so there 5% is a sequence-quality gate. A LexicMap
# row is on the reference grid with '-' wherever no HSP reached, so the same limit becomes a
# coverage gate; from-lexicmap derives its value from the cohort's own coverage gate instead.
DEFAULT_MAX_GAP_FRAC = 0.05


def auto_min_len(ref: Path) -> int | None:
    """95% of the shortest record in `ref`, or None (with the reason printed) if it has none.

    The shortest, not the first: with several references (a multi-anchor panel), each
    accession's primary is reconstructed against whichever one won, so a floor sized on
    record 1 could drop every reconstruction won by a shorter reference.
    """
    records = [r for r in fasta_record_lengths(ref) if r[1] > 0]
    if not records:
        print(f"sieve: {ref}: no record with sequence; cannot derive --min-len.",
              file=sys.stderr)
        return None
    rid, length = min(records, key=lambda r: r[1])
    if len(records) > 1:
        print(f"[sieve] {ref} has {len(records)} records; sizing --min-len on the shortest "
              f"({rid}, {length} nt).", flush=True)
    # Reconstructions are typically a few bases shorter than the reference (the walk may
    # not reach the exact termini). Defaulting to the exact query length would then drop
    # every haplotype, so the auto-default takes a 5% haircut.
    return int(length * 0.95)


# PM-01: `--warmup` is read with `atoi` and the admission test is `n_mine > warmup`, so a
# value nothing can exceed disarms the gate. After 2,000 records the binary used to reject a
# record whose count of never-seen k-mers exceeded `mean + outlier_sd*sqrt(var) + 0.5`, and a
# rejected record is never sketched, promoted or counted, so the rest of its clade was
# rejected with it: in a cohort of near-identical unique haplotypes the cut-off is about
# 1/(3*sqrt(n)) novel divergence, 0.5% at n = 5,000, whatever the gene length. A divergent
# lineage placed after record 2,000 vanished from the MSA and tree at exit 0, and which
# lineages survived depended on file order. `--outlier-sd 1e9` is not the fix: it still
# rejects when the variance is 0.
ADMISSION_OFF_WARMUP = 2 ** 31 - 1
_REJECT_RE = re.compile(r"REJECT \(error\)\s*:\s*(\d+)")

# PM-02: the promoted set's LSH index is a 2^20-slot open-addressing table that is never
# grown, and lookup stops only at an empty slot, so once 2^20 distinct band keys are stored
# the next lookup of an absent key loops forever at 100% CPU. Each promoted record inserts up
# to 16 band keys; 65,536 records with fully distinct bands is the worst case, and closer
# cohorts fill at roughly 150k to 430k promoted. This is upstream's code
# (github.com/stevenweaver/kmer-diversity-sieve) and is not patched here; the wrapper refuses
# rather than hanging. See the backlog entry for the upstream fix (grow at half load).
LSH_TABLE_LIMIT = 65536

# PM-08: the per-record k-mer buffer is a fixed 1<<16 entries and is never grown, so k-mers
# past position ~65,550 are dropped without a message and two records identical over their
# first 65.5 kb are judged near-duplicates. Also upstream's; warned about here.
SKETCH_LENGTH_LIMIT = 65550


def report_rc(rc: int, what: str, prog: str) -> int:
    """Turn a `subprocess` return code into a message and a shell-style exit status.

    PM-19. A signal death is negative here, and it was printed and returned as is, so the
    shell saw 252 for SIGILL and 250 for SIGABRT and the message read "exited -4". The sieve
    Makefile builds with `-march=native`, so a binary moved to another CPU dies with SIGILL,
    which is exactly the case the old message hid.
    """
    if rc < 0:
        try:
            sig = signal.Signals(-rc).name
        except ValueError:
            sig = f"signal {-rc}"
        hint = (" (built with -march=native on another CPU? rebuild it here with "
                "`make serial`)" if -rc == signal.SIGILL and prog == "sieve" else "")
        print(f"{prog}: {what} was killed by {sig}{hint}", file=sys.stderr)
        return 128 - rc
    if rc:
        print(f"{prog}: {what} exited {rc}", file=sys.stderr)
    return rc


class DedupReport(NamedTuple):
    """What one dedup pass saw. The counts used to be two ints and a second read of the file."""
    n_in: int
    n_out: int
    n_empty: int          # header-only records, which used to be counted as duplicates
    headers: list         # the original header line of record i, which the binary never sees
    short: int            # below --min-len, which the binary's completeness gate will drop
    gappy: int            # over --max-gap-frac
    longest: int = 0      # the longest record, for the sketch-length limit (PM-08)


def byte_exact_dedup(src: Path, dst: Path, min_len: int = 0,
                     max_gap_frac: float = 1.0) -> DedupReport:
    """Write to `dst` the first record of each distinct sequence in `src`, in order.

    Returns `DedupReport`. Sequences are compared upper-cased with whitespace removed, so
    line wrapping does not keep two copies apart. `src` may be gzipped; `dst` is always
    plain, which is how the sieve binary, which reads only plain FASTA, gets a gzipped input.

    A later identical copy has novelty 0 and is never promoted, so dropping it loses nothing.
    What duplicates do change is admission: each copy adds a zero to the binary's running
    singleton-count statistics. That gate is now switched off (PM-01), so dedup is a speed-up
    and a memory saving, not a correctness fix.

    **The records are renumbered `>r0`, `>r1`, ... and the headers returned separately**
    (PM-04, PM-13). The binary folds a month bucket parsed from `|YYYY-MM` in the header into
    every LSH band key, so a record was compared only with promoted records of the same
    month: twelve copies of one sequence with `|2021-01-15` ... `|2021-12-15` headers all
    survived, while the same twelve with plain ids collapsed to one. It also selects what to
    emit by the first 255 characters of the header, so two records sharing a header were both
    emitted when only one had been promoted, and a long `from-lexicmap` header lost its tail.
    Neutral ids remove both, and the wrapper puts the original headers back.

    Sequences are keyed by a 16-byte digest rather than held whole (PM-24): the set used
    about 1 GB at 250,000 gene-length rows, transient but before the binary even starts.
    """
    seen: set[bytes] = set()
    headers: list[str] = []
    n_in = n_out = n_empty = short = gappy = longest = 0
    header: str | None = None
    seq_lines: list[str] = []

    def _emit(fout) -> None:
        nonlocal n_out, n_empty, short, gappy, longest
        if header is None:
            return
        # PM-07: the binary strips only \n and \r and counts every other byte as an N, while
        # the wrapper's own predictor stripped all whitespace, so a file with trailing spaces
        # was dropped by the completeness gate while the wrapper predicted no loss at all.
        # Normalise once, here, so both see the same thing.
        seq = "".join("".join(seq_lines).split())
        if not seq:
            n_empty += 1
            return
        key = hashlib.blake2b(seq.upper().encode(), digest_size=16).digest()
        if key in seen:
            return
        seen.add(key)
        # The gate counts measured in the same pass that writes the file (PM-24).
        if len(seq) > longest:
            longest = len(seq)
        if len(seq) < min_len:
            short += 1
        elif 1.0 - sum(map(seq.count, "ACGTacgt")) / len(seq) > max_gap_frac:
            gappy += 1
        fout.write(f">r{len(headers)}\n{seq}\n")
        headers.append(header)
        n_out += 1

    with open_text(src) as fin, dst.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                _emit(fout)
                header = line
                seq_lines = []
                n_in += 1
            elif header is None:
                if line.strip():
                    raise ValueError(f"{src}: text before the first '>' header "
                                     f"({line.strip()[:40]!r})")
            else:
                seq_lines.append(line)
        _emit(fout)
    return DedupReport(n_in, n_out, n_empty, headers, short, gappy, longest)


def completeness_failures(path: Path, min_len: int, max_gap_frac: float) -> tuple[int, int]:
    """Records the binary's completeness gate will drop: (shorter than min_len, too gappy).

    Mirrors l1_pipeline's check so the wrapper can say why records were dropped (the binary
    reports one LOW-QUAL count and calls every non-ACGT character "N"). Length counts every
    character, '-' included; the gap fraction is the share that is not A/C/G/T in either case.
    """
    short = gappy = 0

    def _check(seq: str) -> None:
        nonlocal short, gappy
        if len(seq) < min_len:
            short += 1
        elif seq and 1.0 - sum(map(seq.count, "ACGTacgt")) / len(seq) > max_gap_frac:
            gappy += 1

    seq_lines: list[str] = []
    in_record = False
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if in_record:
                    _check("".join(seq_lines))
                in_record = True
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if in_record:
            _check("".join(seq_lines))
    return short, gappy


def sieve_run(
    input_fasta: Path,
    output_fasta: Path,
    *,
    ref: Path | None = None,
    min_len: int | None = None,
    max_gap_frac: float = DEFAULT_MAX_GAP_FRAC,
    threshold: float = 0.005,
    sketch: int = 256,
    sieve_bin: str | None = None,
) -> int:
    """Run the sieve on `input_fasta`, concatenate rank shards to `output_fasta`.

    If `min_len` is None: derived as 95% of the shortest record in `ref`. One of
    `min_len` or `ref` must be given. Exact duplicates are removed before the binary runs
    (see `byte_exact_dedup`); `max_gap_frac` is passed to it as `--max-n-frac`.

    Returns process exit code (0 on success).
    """
    resolved = resolve_sieve(sieve_bin)
    if not resolved:
        print(not_found_message("sieve", "sieve"), file=sys.stderr)
        return EX_TOOL_NOT_FOUND

    if min_len is None:
        if ref is None:
            print("sieve: pass --min-len, or --ref to derive it from the query.",
                  file=sys.stderr)
            return 2
        try:
            min_len = auto_min_len(ref)
        except (OSError, ValueError) as e:
            print(f"sieve: cannot read --ref: {e}", file=sys.stderr)
            return 1
        if min_len is None:
            return 2

    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)
    try:
        output_fasta.parent.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        # PM-20: a path component that is a regular file raises FileExistsError here, which
        # used to leave a traceback.
        print(f"sieve: cannot create the output directory "
              f"{output_fasta.parent}: {e}", file=sys.stderr)
        return 1

    emit_prefix = output_fasta.with_suffix("")

    started = time.time()

    def emit_summary(status, n_in=0, n_unique=0, n_out=0, short=0, gappy=0, error=None,
                     records=()):
        """The summary, on every exit path including the failures (docs/summaries.md).

        Before 0.2.0 this stage wrote no machine-readable output at all: the dedup counts,
        the completeness-gate drops and the promotion rate went to stdout and nowhere else,
        so a wrapping workflow could not tell a sieve that admitted everything from one that
        dropped most of its input.
        """
        try:
            summary_io.write(
                output_fasta.with_suffix(".summary.json"),
                tool=summary_io.tool("sieve", __version__, external=[
                    {"name": "l1_pipeline_serial", "version": None, "path": str(resolved)}]),
                inputs={"references": None,
                        "source": summary_io.source("file", str(input_fasta), n_in)},
                settings={"min_len": min_len, "max_gap_frac": max_gap_frac,
                          "threshold": threshold, "sketch": sketch},
                run=summary_io.run_block(status, started, time.time() - started,
                                         n_records=n_out, n_ok=n_out, error=error),
                discards={"records": {"exact_duplicates": n_in - n_unique,
                                      "shorter_than_min_len": short,
                                      "too_gappy": gappy,
                                      "not_promoted": max(0, n_unique - short - gappy - n_out)}},
                cohort={"n_in": n_in, "n_unique": n_unique, "n_promoted": n_out,
                        "admission_rate": round(n_out / n_in, 6) if n_in else 0.0},
                records=list(records),
            )
        except OSError:
            pass

    # PM-05, PM-27: the binary runs in a private directory with neutral names. Its shards
    # used to be written beside the user's output and found with
    # `glob(f"{stem}.rank*.fasta")`, unescaped, which deleted `cohort.ranked_by_abundance
    # .fasta` and `cohort.rank1_manual.fasta` before the input was even read, merged any
    # look-alike into the output, missed the real shard when the stem held `[` or `]`, and
    # let two concurrent runs in one directory sweep each other. The scratch names were also
    # `.tmp_<user name>.<pid>`, predictable enough to plant a symlink at.
    try:
        work = Path(tempfile.mkdtemp(prefix=".sieve_", dir=output_fasta.parent))
    except OSError as e:
        # PM-20: an unwritable or full output directory is the output's problem, not the
        # input's, and this is the first thing that touches it.
        print(f"sieve: cannot write beside the output {output_fasta}: {e}", file=sys.stderr)
        emit_summary("failed", error=f"cannot write beside the output: {e}")
        return 1
    deduped = work / "in.fasta"
    emit_prefix = work / "promoted"
    try:
        try:
            rep = byte_exact_dedup(input_fasta, deduped, min_len, max_gap_frac)
        except ValueError as e:
            print(f"sieve: cannot read the input FASTA: {e}", file=sys.stderr)
            emit_summary("failed", error=f"cannot read the input FASTA: {e}")
            return 1
        except OSError as e:
            # PM-20: `byte_exact_dedup` both reads the input and writes the scratch copy, so
            # an unwritable or full output directory was reported as "cannot read the input
            # FASTA", which sends the reader to the wrong file.
            if e.filename and e.filename != str(input_fasta):
                print(f"sieve: cannot write beside the output {output_fasta}: {e}",
                      file=sys.stderr)
                emit_summary("failed", error=f"cannot write beside the output: {e}")
            else:
                print(f"sieve: cannot read the input FASTA: {e}", file=sys.stderr)
                emit_summary("failed", error=f"cannot read the input FASTA: {e}")
            return 1
        n_in, n_unique, short, gappy = rep.n_in, rep.n_out, rep.short, rep.gappy
        if n_unique == 0:
            print(f"sieve: {input_fasta} has no sequences.", file=sys.stderr)
            emit_summary("failed", n_in=n_in, error="the input FASTA has no sequences")
            return 1
        # PM-15: a header-only record is not an exact duplicate, and counting it as one made
        # "4 -> 2 records (2 exact duplicates removed)" of two real and two empty records.
        empty_note = f", {rep.n_empty} empty records dropped" if rep.n_empty else ""
        print(f"[sieve] dedup: {n_in} -> {n_unique} records "
              f"({n_in - n_unique - rep.n_empty} exact duplicates removed, first copy kept"
              f"{empty_note})", flush=True)
        if short or gappy:
            print(f"[sieve] completeness gate drops {short + gappy} of {n_unique}: "
                  f"{short} shorter than {min_len} nt (--min-len), {gappy} with more than "
                  f"{max_gap_frac:.0%} gaps or other non-ACGT characters (--max-gap-frac)",
                  flush=True)

        # PM-08: the binary's per-record k-mer buffer is a fixed 1<<16 entries and is never
        # grown, so k-mers past position ~65,550 are dropped silently and two records
        # identical over their first 65.5 kb are judged near-duplicates: one is dropped at
        # exit 0 with nothing on stderr. Upstream's code, so this warns rather than patching.
        if rep.longest > SKETCH_LENGTH_LIMIT:
            print(f"[warn] the longest record is {rep.longest:,} nt and the sieve binary "
                  f"sketches only the first {SKETCH_LENGTH_LIMIT:,} of each record "
                  f"(upstream bug, kmer-diversity-sieve): records identical over that "
                  f"prefix are judged near-duplicates whatever follows it.", file=sys.stderr)

        # PM-02: refuse rather than hand the binary an input it will hang on. Only admitted
        # records reach the LSH table.
        admitted = n_unique - short - gappy
        if admitted > LSH_TABLE_LIMIT:
            msg = (f"{admitted:,} records would be admitted, and the sieve binary's LSH "
                   f"index is a fixed 2^20-slot table that is never grown: past "
                   f"{LSH_TABLE_LIMIT:,} admitted records it can loop forever instead of "
                   f"finishing (upstream bug, kmer-diversity-sieve). Raise --min-len or "
                   f"--max-gap-frac to admit fewer, or sieve a subset. Chunking and "
                   f"re-sieving is not a way round it: the union can exceed the limit again.")
            print(f"sieve: {msg}", file=sys.stderr)
            emit_summary("failed", n_in=n_in, n_unique=n_unique, short=short, gappy=gappy,
                         error=msg)
            return 1

        cmd = [
            resolved,
            "--input", str(deduped),
            "--min-len", str(min_len),
            "--max-n-frac", str(max_gap_frac),
            "--threshold", str(threshold),
            "--warmup", str(ADMISSION_OFF_WARMUP),
            "--sketch", str(sketch),
            "--emit-promoted", str(emit_prefix),
        ]
        print(f"[sieve] {' '.join(cmd)}", flush=True)
        # PM-01: the binary reports its rejections on stdout and nothing read them.
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
        sys.stdout.write(proc.stdout)
        sys.stdout.flush()
        rc = report_rc(proc.returncode, "kmer-diversity-sieve", "sieve")
        if rc != 0:
            emit_summary("failed", n_in=n_in, n_unique=n_unique, short=short, gappy=gappy,
                         error=f"kmer-diversity-sieve exited {proc.returncode}")
            return rc
        m = _REJECT_RE.search(proc.stdout)
        if m is None:
            print("[warn] the sieve binary did not report a REJECT count; it may be a "
                  "version that still has the admission gate armed.", file=sys.stderr)
        elif int(m.group(1)) > 0:
            msg = (f"the sieve binary rejected {m.group(1)} record(s) on admission, which "
                   f"should be impossible with the gate disarmed. The binary is probably "
                   f"older than the wrapper; rebuild it.")
            print(f"sieve: {msg}", file=sys.stderr)
            emit_summary("failed", n_in=n_in, n_unique=n_unique, short=short, gappy=gappy,
                         error=msg)
            return 1

        # PM-14: the output used to be opened for writing before anything was counted, so a
        # run that promoted nothing left a 0-byte file over the previous result. Staged, and
        # renamed in only when there is something to publish. PM-13: the binary's own output
        # headers are the truncated neutral ids, so the saved originals go back here.
        shards = sorted(work.glob("promoted.rank*.fasta"),
                        key=lambda q: int(q.name[len("promoted.rank"):-len(".fasta")]))
        staged = work / "out.fasta"
        n_out = 0
        with staged.open("w") as fh:
            for rf in shards:
                with rf.open() as rin:
                    for line in rin:
                        if line.startswith(">"):
                            n_out += 1
                            idx = line[1:].split()[0]
                            line = rep.headers[int(idx[1:])] if idx[:1] == "r" else line
                        fh.write(line)
        if n_out:
            os.replace(staged, output_fasta)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    # One record per promoted haplotype: its id and the hash of its sequence, which is how
    # a summary names a sequence without carrying it (docs/summaries.md).
    promoted = []
    if n_out:
        ident, seq = None, []
        with output_fasta.open() as fh:
            for line in fh:
                if line.startswith(">"):
                    if ident is not None:
                        promoted.append({"id": ident,
                                         "seq_sha256": summary_io.seq_sha256("".join(seq)),
                                         "length": len("".join(seq))})
                    ident, seq = line[1:].split()[0], []
                else:
                    seq.append(line.strip())
        if ident is not None:
            promoted.append({"id": ident, "seq_sha256": summary_io.seq_sha256("".join(seq)),
                             "length": len("".join(seq))})
    print(f"[sieve] {n_in} -> {n_out} haplotypes "
          f"(tau={threshold}, sketch={sketch}) -> {output_fasta}", flush=True)
    emit_summary("ok" if n_out else "failed", n_in=n_in, n_unique=n_unique, n_out=n_out,
                 short=short, gappy=gappy, records=promoted,
                 error=None if n_out else "promoted 0 haplotypes")
    if n_out == 0:
        # PM-14: `is_novel` on an empty promoted set returns 1, so zero promoted cannot be a
        # --threshold problem: every record failed the completeness gate or carries no
        # A/C/G/T 15-mer. The old hint led with the one cause it could not be.
        print(f"sieve: promoted 0 of {n_unique} unique records; {output_fasta} is unchanged. "
              f"{short} were shorter than --min-len {min_len} and {gappy} were over "
              f"--max-gap-frac {max_gap_frac:.0%}"
              + ("." if short or gappy else
                 ", so the records carry no A/C/G/T 15-mer the sieve can sketch."),
              file=sys.stderr)
        return 1
    return 0


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "sieve",
        help="Thin near-duplicate haplotypes with the k-mer diversity sieve.",
        description=(
            "MinHash/LSH lineage-level reducer. Use only for LARGE cross-accession "
            "cohorts of near-duplicate lineages, NOT to preserve intra-host SNVs "
            "(the walker's --emit-vcf / --emit-snv-table do that). Exact duplicates are "
            "removed first, keeping the first copy."
        ),
    )
    p.add_argument("-i", "--input", required=True, help="Input FASTA")
    p.add_argument("-o", "--output", required=True, help="Output (promoted) FASTA")
    p.add_argument("--min-len", type=positive_int, default=None,
                   help="Minimum sequence length gate. Default: auto = 95%% of the "
                        "shortest record in --ref.")
    p.add_argument("--ref", default=None,
                   help="Reference CDS FASTA, used only to auto-derive --min-len "
                        "when --min-len is not given.")
    p.add_argument("--max-gap-frac", type=fraction, default=DEFAULT_MAX_GAP_FRAC,
                   help="Completeness gate: drop a sequence whose share of gaps ('-') or "
                        "other non-ACGT characters exceeds this (default 0.05). Walker rows "
                        "have no gaps. For a LexicMap cohort, whose rows carry '-' where no "
                        "HSP reached, use 1 - the cohort's coverage gate + 0.05, as "
                        "from-lexicmap does.")
    p.add_argument("--threshold", type=fraction, default=0.005,
                   help="LSH promotion threshold tau (default 0.005, gene-scale). "
                        "The sieve binary's own default of 0.15 targets archive-scale.")
    p.add_argument("--sketch", type=positive_int, default=256,
                   help="MinHash sketch size m (default 256, gene-scale). The "
                        "estimator cannot resolve novelty finer than 1/m.")
    p.add_argument("--sieve-bin", default=None,
                   help="Path to l1_pipeline_serial (auto-discovered otherwise).")
    p.set_defaults(func=run)


def run(args) -> int:
    return sieve_run(
        Path(args.input),
        Path(args.output),
        ref=Path(args.ref) if args.ref else None,
        min_len=args.min_len,
        max_gap_frac=args.max_gap_frac,
        threshold=args.threshold,
        sketch=args.sketch,
        sieve_bin=args.sieve_bin,
    )
