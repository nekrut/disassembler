"""sieve: thin near-duplicate haplotypes with the k-mer diversity sieve.

Wraps `extern/kmer-diversity-sieve/l1_pipeline_serial`. The sieve is a MinHash/LSH
lineage-level reducer: on defaults tuned to gene scale (tau=0.005, sketch=256) it
collapses a cohort of ~identical lineages while preserving topology; the sieve
binary's own defaults (tau=0.15, sketch=64) target archive-scale collapse and reduce
a gene-scale cohort to ~1 representative. See docs/sieve.md and
docs/issue_sieve_calibration.md.

Both the CLI and the composite pipelines call `sieve_run(...)`, which removes exact
duplicates first and reports what the binary's completeness gate will drop.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

from disassembler.resolve import resolve_sieve
from disassembler.validators import fasta_record_lengths, fraction, positive_int

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


def byte_exact_dedup(src: Path, dst: Path) -> tuple[int, int]:
    """Write to `dst` the first record of each distinct sequence in `src`, in order.

    Returns (n_in, n_out). Headers are kept verbatim from the first copy; sequences are
    compared upper-cased with line breaks removed, so line wrapping does not keep two copies
    apart. Empty records are dropped.

    A later identical copy has novelty 0 and is never promoted, so dropping it loses nothing.
    What duplicates do change is admission: each copy adds a zero to the binary's running
    singleton-count statistics, and once its 2000-record warmup is over the shrunken mean and
    SD make it reject genuinely new sequences as outliers ("REJECT (error)"). On four C. auris
    cohorts (2026-10-30) that rejected 26 to 273 sequences, and none once deduplicated.
    Below the warmup the two inputs give byte-identical output
    (tests/test_dedup_sieve_equivalence.py).
    """
    seen: set[str] = set()
    n_in = 0
    n_out = 0
    header: str | None = None
    seq_lines: list[str] = []

    def _emit(fout) -> None:
        nonlocal n_out
        if header is None:
            return
        seq = "".join(seq_lines).replace("\n", "").replace("\r", "").upper()
        if seq and seq not in seen:
            seen.add(seq)
            fout.write(header)
            fout.writelines(seq_lines)
            n_out += 1

    with src.open() as fin, dst.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                _emit(fout)
                header = line
                seq_lines = []
                n_in += 1
            else:
                seq_lines.append(line)
        _emit(fout)
    return n_in, n_out


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
        print(
            "sieve: kmer-diversity-sieve not found. Build it "
            "(`cd extern/kmer-diversity-sieve && make serial`), put it on PATH, "
            "or pass --sieve-bin.",
            file=sys.stderr,
        )
        return 127

    if min_len is None:
        if ref is None:
            print("sieve: pass --min-len, or --ref to derive it from the query.",
                  file=sys.stderr)
            return 2
        min_len = auto_min_len(ref)
        if min_len is None:
            return 2

    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    emit_prefix = output_fasta.with_suffix("")
    # Clear any stale shards from a previous run so the concat glob can't pick them up.
    for stale in output_fasta.parent.glob(f"{emit_prefix.name}.rank*.fasta"):
        stale.unlink()

    # The binary runs on a deduplicated copy, written beside the output and removed after.
    deduped = output_fasta.with_name(f".tmp_{output_fasta.stem}.deduped.{os.getpid()}.fasta")
    try:
        n_in, n_unique = byte_exact_dedup(input_fasta, deduped)
        if n_unique == 0:
            print(f"sieve: {input_fasta} has no sequences.", file=sys.stderr)
            return 1
        print(f"[sieve] dedup: {n_in} -> {n_unique} records "
              f"({n_in - n_unique} exact duplicates removed, first copy kept)", flush=True)
        short, gappy = completeness_failures(deduped, min_len, max_gap_frac)
        if short or gappy:
            print(f"[sieve] completeness gate drops {short + gappy} of {n_unique}: "
                  f"{short} shorter than {min_len} nt (--min-len), {gappy} with more than "
                  f"{max_gap_frac:.0%} gaps or other non-ACGT characters (--max-gap-frac)",
                  flush=True)

        cmd = [
            resolved,
            "--input", str(deduped),
            "--min-len", str(min_len),
            "--max-n-frac", str(max_gap_frac),
            "--threshold", str(threshold),
            "--warmup", "2000",
            "--sketch", str(sketch),
            "--emit-promoted", str(emit_prefix),
        ]
        print(f"[sieve] {' '.join(cmd)}", flush=True)
        rc = subprocess.run(cmd).returncode
    finally:
        deduped.unlink(missing_ok=True)
    if rc != 0:
        print(f"sieve: kmer-diversity-sieve exited {rc}", file=sys.stderr)
        return rc

    rank_files = sorted(output_fasta.parent.glob(f"{emit_prefix.name}.rank*.fasta"))
    with output_fasta.open("w") as fh:
        for rf in rank_files:
            fh.write(rf.read_text())

    n_out = sum(1 for line in output_fasta.open() if line.startswith(">"))
    print(f"[sieve] {n_in} -> {n_out} haplotypes "
          f"(tau={threshold}, sketch={sketch}) -> {output_fasta}", flush=True)
    if n_out == 0:
        print("sieve: promoted 0 haplotypes; try a lower --threshold, a lower --min-len "
              "(reconstructions may be shorter than the reference), or a higher "
              "--max-gap-frac (gapped rows, such as LexicMap partials).",
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
