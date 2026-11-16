"""msa: reference-guided codon MSA + FastTree phylogeny from reconstructed haplotypes.

Aligns each haplotype to the reference CDS with MAFFT `--keeplength --addfragments`
so reference coordinates 1..L are preserved and the reading frame stays invariant,
then builds a FastTree GTR tree from the alignment.

MAFFT never reverse-complements, and `--keeplength` makes any input fit the reference's
length, so a row on the other strand used to align as noise and exit 0. Each row is
oriented to the reference by shared k-mers first, and the alignment is checked for
identity to the reference after.
"""

from __future__ import annotations

import argparse
import functools
import tempfile
import os
import re
import shutil
import stat
import statistics
import subprocess
import time
import sys
from pathlib import Path

from disassembler.textio import open_text
from disassembler.validators import fraction, positive_int

from disassembler import __version__
from disassembler import summary as summary_io

def find_bin(name: str, *aliases: str) -> str:
    """Resolve a tool on PATH.

    `aliases` are alternative spellings of the same tool, tried in order. bioconda ships
    FastTree under both `fasttree` and `FastTree` depending on the build.

    Raises rather than returning the bare name so an unresolvable tool fails here with
    a message naming it, instead of surfacing as a FileNotFoundError from subprocess
    after the expensive stages have already run. There is no fallback location: the one
    there was (a `codon` env under ~/micromamba) belonged to one old setup.
    """
    names = (name, *aliases)
    for n in names:
        found = shutil.which(n)
        if found:
            return found
    spellings = " or ".join(repr(n) for n in names)
    raise FileNotFoundError(
        f"{spellings} not found on PATH. Activate the `disassembler` env "
        f"(environment.yml pins mafft and fasttree), or put {name} on PATH."
    )


def clean_id(s: str) -> str:
    """A record id reduced to `[A-Za-z0-9_]`, **not truncated** (PM-26).

    It used to be cut to 40 characters, which is exactly where a GISAID or NCBI header
    carries the part that tells two records apart: two
    `hCoV-19/USA/CA-CDC-QDX12345678/2021|EPI_ISL_100000{1,2}|x` records became
    `hCoV_19_USA_CA_CDC_QDX12345678_2021_EPI_` and `..._EPI__1`, with the accessions gone,
    so a tree tip could be traced back only by file order. FastTree, Newick and HyPhy all
    accept long names.
    """
    return re.sub(r"[^A-Za-z0-9_]", "_", s).strip("_") or "seq"


def assign_ids(raw_ids: list) -> tuple:
    """`(ids, changed)`: a unique alignment id per raw id, and whether any differs from it.

    Two passes (PM-26): an id that is already clean, unique and not `Reference` keeps
    itself, whatever the order it appears in. Only the rest are renamed, and a suffix never
    takes a name pass 1 reserved. Before this, an id that was already clean could be renamed
    while another record took its name.

    The suffix search remembers where each base stopped (PM-22). It used to restart at 1 for
    every record, so n ids sharing a prefix cost about n squared probes: 10,000 such ids
    spent 31.8 s in renaming alone, against 0.05 s now, for the same names.
    """
    cleaned = [clean_id(r) for r in raw_ids]
    used = {"Reference"}
    out = [None] * len(raw_ids)
    # Pass 1 claims a name only for the record that already HAS it. Counting how many ids
    # clean to the same name is the wrong test: `SRR1.1` and `SRR1_1` both clean to
    # `SRR1_1`, and it is the record actually called `SRR1_1` that should keep it. Two
    # records with the identical raw id are indistinguishable, so the first taking it is
    # not an order dependence that means anything.
    for i, (raw, c) in enumerate(zip(raw_ids, cleaned)):
        if c == raw and c not in used:
            out[i] = c
            used.add(c)
    next_idx = {}
    for i, c in enumerate(cleaned):
        if out[i] is not None:
            continue
        candidate = c
        if candidate in used:
            idx = next_idx.get(c, 1)
            candidate = f"{c}_{idx}"
            while candidate in used:
                idx += 1
                candidate = f"{c}_{idx}"
            next_idx[c] = idx
        used.add(candidate)
        out[i] = candidate
    return out, any(a != b for a, b in zip(out, raw_ids))


def sanitize_id(s: str, used: set) -> str:
    """Kept for callers outside this module; `assign_ids` is what `run` uses."""
    candidate = base = clean_id(s)
    idx = 1
    while candidate in used:
        candidate = f"{base}_{idx}"
        idx += 1
    used.add(candidate)
    return candidate


BAD_MAFFT = {(7, 525): "segfaults on the --addfragments path msa uses; fixed in 7.526"}


def check_mafft(mafft_bin: str) -> None:
    """Refuse a MAFFT release known to break the path this stage uses (PM-17).

    `docs/packaging.md` has recorded since 0.1.0 that 7.525 segfaults on `--addfragments`,
    and the published recipe leaves MAFFT unpinned, but nothing ever ran `mafft --version`:
    an environment that resolved 7.525 failed with a raw `MAFFT failed:` log and no hint,
    and under `reconstruct` only after the whole walk. Unparsable output passes, so stubs
    and older releases keep working.
    """
    try:
        # stdin closed and a short timeout: a tool asked for its version must never be able
        # to sit waiting for input. A stub `mafft` that only `cat`s does exactly that.
        r = subprocess.run([mafft_bin, "--version"], capture_output=True, text=True,
                           stdin=subprocess.DEVNULL, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return
    m = re.search(r"v(\d+)\.(\d+)", (r.stderr or "") + (r.stdout or ""))
    if m and (int(m[1]), int(m[2])) in BAD_MAFFT:
        raise RuntimeError(f"{mafft_bin} is MAFFT {m[1]}.{m[2]}, which "
                           f"{BAD_MAFFT[(int(m[1]), int(m[2]))]}. Install mafft>=7.526.")


@functools.lru_cache(maxsize=8)
def tool_version(binary: str, *argv: str) -> str | None:
    """The first non-empty line a tool prints when asked its version (PM-16).

    stdin is closed and the timeout is short, because a tool that does not recognise the
    flag may fall back to reading a sequence from standard input and wait there; cached,
    because the answer cannot change inside one run.
    """
    try:
        r = subprocess.run([binary, *argv], capture_output=True, text=True,
                           stdin=subprocess.DEVNULL, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    for line in ((r.stdout or "") + "\n" + (r.stderr or "")).splitlines():
        if line.strip():
            return line.strip()
    return None


def is_special(path: Path) -> bool:
    """Is this output something a rename would destroy: a symlink, device, FIFO or socket?"""
    try:
        return path.is_symlink() or not stat.S_ISREG(os.stat(path).st_mode)
    except (FileNotFoundError, NotADirectoryError):
        return False


def publish(tmp: Path, dest: Path) -> None:
    """Move `tmp` onto `dest`, by copy when `dest` is not a plain file (PM-09).

    A rename swaps a directory entry, so `--out-tree /dev/null` turned the null device into
    a regular file holding a Newick tree (as root, which is the default in a container), and
    a FIFO output was replaced so its reader blocked forever. Where the target's directory
    cannot take new files at all, `/dev` for a normal user or `/dev/fd` for bash `>(...)`,
    the temp `open` raised an uncaught exception after MAFFT had already run.
    """
    if is_special(dest):
        with open(tmp, "rb") as src, open(dest, "wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        os.replace(tmp, dest)


STOP_CODONS = ("TAA", "TAG", "TGA")


def codon_qc(row: str, input_acgt: int) -> dict:
    """What `--keeplength` did to one row, which nothing used to look at (PM-03).

    MAFFT is run with `--auto --keeplength --addfragments`, which deletes every base a row
    has beyond the reference and puts gaps where nucleotide scoring wants them, not on codon
    boundaries. The only check afterwards was `len(row) != ref_len`, which `--keeplength`
    guarantees, printed under the banner "Validating codon lengths...". So a +1 insertion row
    came out identical to the reference, a -1 deletion shifted the rest of the row into the
    reference frame, in-frame insertions vanished from the alignment and the tree, and an
    internal stop passed: the run ended "Codon alignment complete", exit 0. The alignment is
    not changed here; it is described.
    """
    aligned_acgt = real_bases(row)
    internal = [len(m.group(0)) for m in re.finditer(r"(?<=[^-])-+(?=[^-])", row)]
    return {
        "trimmed": input_acgt - aligned_acgt,
        "internal_gap_runs": internal,
        "split_codons": [c // 3 + 1 for c in range(0, len(row), 3)
                         if 0 < row[c:c + 3].count("-") < 3],
        "internal_stops": [c // 3 + 1 for c in range(0, len(row) - 3, 3)
                           if row[c:c + 3].upper() in STOP_CODONS],
        "frameshift": (input_acgt - aligned_acgt) % 3 != 0 or any(r % 3 for r in internal),
    }


def real_bases(seq: str) -> int:
    """A, C, G and T in either case; gaps, N and other IUPAC codes do not count."""
    return sum(map(seq.count, "ACGTacgt"))


# Orientation is decided on shared k-mers. Measured on four real cohorts (HIV-1 pol, a
# SARS-CoV-2 spike walker run, LexicMap C. auris FKS1 and H. pylori vacA), every row shared
# at least 896 15-mers with its reference's strand and at most 23 with the other.
ORIENT_K = 15
# Every 4th k-mer of a row is plenty at those counts, and keeps a 250,000-row cohort cheap.
ORIENT_STEP = 4
_COMPLEMENT = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def revcomp(seq: str) -> str:
    """Reverse complement, IUPAC codes included; gaps stay gaps."""
    return seq.translate(_COMPLEMENT)[::-1]


def kmer_set(seq: str, step: int = 1, k: int = None) -> set:
    """Upper-case k-mers of `seq` made only of A/C/G/T, taking every `step`-th start."""
    seq = seq.upper()
    k = ORIENT_K if k is None else k
    acgt = set("ACGT")
    return {seq[i:i + k] for i in range(0, len(seq) - k + 1, step)
            if acgt.issuperset(seq[i:i + k])}


ORIENT_FALLBACK_K = 11


def orientation(seq: str, fwd: set, rev: set, fallback: tuple = None) -> int:
    """+1 if `seq` looks like the reference, -1 like its reverse complement, 0 undecided.

    PM-12. The exact 15-mers sampled every 4th position find nothing in common with either
    strand once a row is about 30% divergent, and a tie was written out as given: MAFFT does
    not reverse-complement, so a divergent minus-strand row aligned at about 0.43 identity
    and sat in the alignment and the tree as one noise taxon, with a warning that named the
    wrong cause ("check they are the same gene"). On a tie, every 11-mer is compared against
    the reference's, and the call is made only with a clear margin, so a row that really is
    undecidable still comes back 0.
    """
    kmers = kmer_set(seq, ORIENT_STEP)
    n_fwd, n_rev = len(kmers & fwd), len(kmers & rev)
    if n_fwd == n_rev and fallback is not None:
        small = kmer_set(seq, 1, ORIENT_FALLBACK_K)
        a, b = len(small & fallback[0]), len(small & fallback[1])
        if max(a, b) < 2 * min(a, b) + 3:
            return 0
        return (a > b) - (b > a)
    return (n_fwd > n_rev) - (n_rev > n_fwd)


ACGT_SET = frozenset("ACGT")


def identity_to_reference(ref: str, row: str) -> float:
    """Matches over the alignment columns where both the reference and the row have A/C/G/T.

    PM-23: the reference used to be upper-cased again for every row, and the loop tested
    membership in a string literal per column. Same answer, one pass, no per-column set-up.
    """
    ref_u, row_u = ref.upper(), row.upper()
    compared = matches = 0
    acgt = ACGT_SET
    for a, b in zip(ref_u, row_u):
        if a in acgt and b in acgt:
            compared += 1
            if a == b:
                matches += 1
    return matches / compared if compared else 0.0


def identity_ok(aligned_recs, min_identity: float) -> bool:
    """Print identity to the reference (the first record); False if the median is too low.

    `--keeplength` fits any input to the reference's length, so the length check above can
    not catch an alignment that is noise. Identity can: a strand mismatch measured 0.40-0.50.
    """
    ref = str(aligned_recs[0].seq)
    # PM-23: the reference is upper-cased once, not once per row. Rows that share a
    # sequence (PM-10's expansion) are measured once and the value reused, which is what
    # keeps the median weighted by copy number without doing the work N times.
    ref_u = ref.upper()
    cache = {}
    ids = []
    for r in aligned_recs[1:]:
        s = str(r.seq)
        if s not in cache:
            cache[s] = identity_to_reference(ref_u, s)
        ids.append((r.id, cache[s]))
    if not ids or min_identity <= 0:
        return True
    median = statistics.median(v for _, v in ids)
    low = [(i, v) for i, v in ids if v < min_identity]
    print(f"[*] Identity to the reference: median {median:.3f}, lowest "
          f"{min(v for _, v in ids):.3f}")
    if low:
        shown = ", ".join(f"{i} ({v:.2f})" for i, v in sorted(low, key=lambda x: x[1])[:5])
        print(f"[warn] {len(low)} of {len(ids)} sequences match the reference at under "
              f"{min_identity:.2f} identity: {shown}"
              f"{f' and {len(low) - 5} more' if len(low) > 5 else ''}. "
              f"Check they are the same gene.", file=sys.stderr)
    if median < min_identity:
        print(f"Error: median identity to the reference is {median:.3f}, under "
              f"--min-identity {min_identity:.2f}: the alignment is mostly noise, so it was "
              f"not written. Check that --ref is the right gene and an in-frame CDS.",
              file=sys.stderr)
        return False
    return True


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "msa",
        help="Codon MSA + FastTree GTR phylogeny from reconstructed haplotypes.",
        description=(
            "MAFFT --keeplength --addfragments against the reference CDS, then FastTree GTR. "
            "Intended for walker output, which is per-accession reconstructions in their own "
            "coordinates. MAFFT does not reverse-complement, so each row is oriented to the "
            "reference by shared 15-mers first, and the alignment is checked for identity to "
            "the reference after. lexicmap output is already on the reference grid; under "
            "lexicmap --require-full-length its only gaps are in-frame deletions, and rows "
            "without one come through MAFFT unchanged. See docs/codon_msa.md."
        ),
    )
    p.add_argument("--ref", "-r", required=True, help="Reference in-frame CDS FASTA")
    p.add_argument(
        "--input", "-i", required=True,
        help="Haplotypes FASTA (walker output, or a lexicmap clean cohort FASTA). Rows on "
             "the reference's other strand are reverse-complemented. See the subcommand "
             "description.",
    )
    p.add_argument("--out-aln", "-a", required=True, help="Output codon MSA FASTA")
    p.add_argument("--out-tree", "-t", default=None, help="Output Newick tree")
    p.add_argument("--no-tree", action="store_true",
                   help="Write only the alignment. FastTree is then neither looked up nor "
                        "run, so an alignment can be produced on a machine that does not "
                        "have it. Exactly one of --out-tree and --no-tree is required.")
    p.add_argument("--mafft-threads", type=positive_int, default=8, help="MAFFT threads")
    p.add_argument("--msa-min-cov", type=fraction, default=0.70,
                   help="Minimum coverage to enter the MSA: a sequence's A/C/G/T bases, as a "
                        "fraction of the reference length (default 0.70). Gaps ('-') and N "
                        "do not count, so a gapped lexicmap row is measured by what it covers.")
    p.add_argument("--min-identity", type=fraction, default=0.60,
                   help="Identity check after alignment: warn about every sequence that "
                        "matches the reference at less than this, and fail without writing "
                        "the alignment if the median does (default 0.60). Same-gene cohorts "
                        "measure 0.88 or more; a sequence on the wrong strand, or of another "
                        "gene, aligns at 0.40-0.50. 0 disables the check.")
    p.set_defaults(func=run)


def _write_summary(args, state, status: str) -> None:
    """The summary, whatever happened (`docs/summaries.md`).

    This stage produced no machine-readable output at all before 0.2.0, which is awkward for
    the one stage that can fail on a tool version: MAFFT 7.525 mishandles `--addfragments`,
    and the only trace of how many rows survived the coverage gate, how many were reverse
    complemented and how many could not be oriented was three lines on stdout.
    """
    try:
        summary_io.write(
            Path(args.out_aln).with_suffix(".summary.json"),
            tool=summary_io.tool("msa", __version__, external=[
                {"name": "mafft", "version": state.get("mafft"), "path": state.get("mafft_bin")},
                {"name": "fasttree", "version": None, "path": state.get("fasttree_bin")},
            ]),
            inputs={
                "references": ([{"id": state["ref_id"], "length": state["ref_len"],
                                 "sha256": state["ref_sha256"]}]
                               if state.get("ref_id") else None),
                "source": summary_io.source("file", str(args.input), state["n_in"]),
            },
            settings={"msa_min_cov": args.msa_min_cov, "min_identity": args.min_identity,
                      "mafft_threads": getattr(args, "mafft_threads", None)},
            run=summary_io.run_block(status, state["started"], time.time() - state["started"],
                                     n_records=state["taxa"], n_ok=state["taxa"],
                                     error=state.get("error")),
            # Rows lost between the input and the alignment. `below_coverage` is the gate;
            # `unoriented` could not be placed on either strand.
            discards={"rows": {
                "below_coverage": max(0, state["n_in"] - state["retained"]),
                "unoriented": state.get("unoriented", 0),
                "lost_in_alignment": max(0, state["retained"] - state["taxa"]),
            }},
            cohort={"n_in": state["n_in"], "retained": state["retained"],
                    "taxa": state["taxa"], "reverse_complemented": state.get("flipped", 0),
                    "reference_length": state["ref_len"]},
            records=state.get("records", []),
        )
    except OSError:
        # The alignment's own directory may be what failed; the caller already said so.
        pass


def run(args) -> int:
    """Align and build the tree, then write a summary on every path."""
    state = {"started": time.time(), "n_in": 0, "retained": 0, "flipped": 0, "unoriented": 0,
             "ref_id": None, "ref_len": 0, "ref_sha256": None, "taxa": 0, "error": None,
             "records": []}
    try:
        rc = _run(args, state)
    except BaseException as e:
        state["error"] = f"{type(e).__name__}: {e}"
        _write_summary(args, state, "failed")
        raise
    if rc != 0 and state["error"] is None:
        state["error"] = f"exited {rc}; see stderr"
    _write_summary(args, state, "ok" if rc == 0 else "failed")
    return rc


def _run(args, state) -> int:
    from Bio import SeqIO  # imported lazily so `disassembler --help` does not need biopython

    # PM-18: `--out-tree` used to be required and FastTree looked up before any work, so an
    # alignment-only run needed FastTree installed even though it was never going to be run.
    want_tree = not getattr(args, "no_tree", False)
    if want_tree and not args.out_tree:
        print("Error: pass --out-tree <path>, or --no-tree to write only the alignment.",
              file=sys.stderr)
        return 2
    if not want_tree and args.out_tree:
        print("Error: --no-tree and --out-tree are mutually exclusive.", file=sys.stderr)
        return 2

    ref_p = Path(args.ref)
    in_p = Path(args.input)
    out_aln = Path(args.out_aln)
    out_tree = Path(args.out_tree) if want_tree else None
    out_aln.parent.mkdir(parents=True, exist_ok=True)
    if out_tree is not None:
        out_tree.parent.mkdir(parents=True, exist_ok=True)

    fasttree_bin = None
    try:
        mafft_bin = find_bin("mafft")
        if want_tree:
            fasttree_bin = find_bin("fasttree", "FastTree")
        state["mafft_bin"] = str(mafft_bin)
        state["fasttree_bin"] = str(fasttree_bin) if fasttree_bin else None
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # PM-17: refuse a MAFFT release known to break --addfragments, before any work.
    try:
        check_mafft(mafft_bin)
    except RuntimeError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # PM-16: nothing recorded which MAFFT or FastTree produced the result, although the
    # published recipe leaves MAFFT unpinned, so two installs could give different files
    # with nothing in either output naming the tools.
    state["mafft_version"] = tool_version(mafft_bin, "--version")
    state["fasttree_version"] = tool_version(fasttree_bin, "-expert") if fasttree_bin else None
    print(f"[*] tools: disassembler {__version__}; "
          f"MAFFT {state['mafft_version'] or 'unknown'}"
          + (f"; FastTree {state['fasttree_version'] or 'unknown'}" if fasttree_bin else ""))

    # A missing, zstd-compressed or binary input ends in one line, not a traceback.
    try:
        with open_text(ref_p) as fh:
            ref_records = list(SeqIO.parse(fh, "fasta"))
        with open_text(in_p) as fh:
            # PM-23: `SimpleFastaParser` gives (title, sequence) and builds no SeqRecord,
            # which is almost all of the 222 MB a 20,000-row cohort cost here. The rows are
            # iterated exactly once either way.
            from Bio.SeqIO.FastaIO import SimpleFastaParser
            input_records = [(title.split(None, 1)[0], seq)
                             for title, seq in SimpleFastaParser(fh)]
    except (OSError, ValueError) as e:
        print(f"Error: cannot read the reference or input FASTA: {e}", file=sys.stderr)
        return 1
    if not ref_records:
        print(f"Error: no records found in reference FASTA {ref_p}", file=sys.stderr)
        return 1
    ref_rec = ref_records[0]
    ref_len = len(ref_rec.seq)
    state["ref_id"] = ref_rec.id
    state["ref_len"] = ref_len
    state["ref_sha256"] = summary_io.seq_sha256(str(ref_rec.seq))
    state["n_in"] = len(input_records)
    if len(ref_records) > 1:
        print(f"[warn] {ref_p}: {len(ref_records)} records; aligning against the first "
              f"({ref_rec.id}, {ref_len} nt). The MSA takes one reference: pass a "
              f"single-record FASTA to choose another.", file=sys.stderr)
    if ref_len % 3 != 0:
        print(f"[warn] Reference CDS length {ref_len} is not divisible by 3!", file=sys.stderr)

    if not input_records:
        print(f"Error: no records found in input FASTA {in_p}", file=sys.stderr)
        return 1

    print(
        f"[*] Processing {len(input_records)} reconstructed sequences against reference "
        f"{ref_rec.id} ({ref_len} bp)..."
    )

    ref_seq = str(ref_rec.seq).upper()
    fwd_kmers, rev_kmers = kmer_set(ref_seq), kmer_set(revcomp(ref_seq))
    # The shorter-k sets, built once, for rows whose 15-mers tie (PM-12).
    fallback = (kmer_set(ref_seq, 1, ORIENT_FALLBACK_K),
                kmer_set(revcomp(ref_seq), 1, ORIENT_FALLBACK_K))

    # PM-06, PM-09, PM-27: one private directory with fixed names, beside the output so
    # `os.replace` stays on one filesystem. The three scratch files used to be named after
    # the user's own files (`.tmp_<name>.<pid>`) in the output directory, where ordinary
    # choices collide: `-i raw/cohort.fasta -a aligned/cohort.fasta` truncated MAFFT's input
    # with its output, and `-r pol.fasta -i ref_pol.fasta` made the reference overwrite the
    # cleaned input so MAFFT aligned the reference to itself and a two-row alignment
    # replaced any earlier good one. Predictable names in a shared directory are also
    # plantable: `open(..., "w")` follows a symlink and truncates its target, and the rename
    # then moved the link itself into place as the alignment.
    tmp_dir = Path(tempfile.mkdtemp(prefix=".tmp_msa_", dir=out_aln.parent))
    clean_input_p = tmp_dir / "input.fasta"
    clean_ref_p = tmp_dir / "reference.fasta"
    tmp_aln = tmp_dir / "aligned.fasta"
    tmp_tree = tmp_dir / "tree.nwk"
    try:
        retained = flipped = unoriented = 0
        rows = []            # (alignment id, representative id) in input order
        rep_of = {}          # oriented sequence -> the id MAFFT is given for it
        input_acgt = {}      # alignment id -> A/C/G/T in the row as it went in
        kept_raw, kept_seq, kept_flipped = [], [], []
        for raw_id, raw_seq in input_records:
            seq = raw_seq.upper()
            # Count bases, not characters: a lexicmap row is reference-length with '-'
            # where no HSP reached, so its length passed this gate whatever its coverage.
            if real_bases(seq) < ref_len * args.msa_min_cov:
                continue
            strand = orientation(seq, fwd_kmers, rev_kmers, fallback)
            if strand < 0:
                seq = revcomp(seq)
                flipped += 1
            elif strand == 0:
                unoriented += 1
            kept_raw.append(raw_id)
            kept_seq.append(seq)
            kept_flipped.append(strand < 0)
            retained += 1

        # PM-26: ids are assigned over the surviving rows at once, so one that is already
        # clean and unique keeps itself whatever order it arrived in.
        ids, ids_changed = assign_ids(kept_raw)
        with open(clean_input_p, "w") as fh:
            for clean, seq in zip(ids, kept_seq):
                input_acgt[clean] = real_bases(seq)
                # PM-10: under --addfragments each row is aligned to the reference alone, so
                # every copy of a sequence gets the same row. Sending one copy and expanding
                # afterwards gives a byte-identical alignment: 600 rows drawn from 30
                # variants took 129.7 s and the 30 distinct ones 8.4 s.
                rep = rep_of.setdefault(seq, clean)
                rows.append((clean, rep))
                if rep == clean:
                    fh.write(f">{clean}\n{seq}\n")

        state["retained"] = retained
        if retained == 0:
            state["error"] = "no sequences passed quality/coverage filtering"
            print("Error: no sequences passed quality/coverage filtering", file=sys.stderr)
            return 1

        print(
            f"[*] {retained}/{len(input_records)} sequences passed coverage filter "
            f"(>= {args.msa_min_cov*100:.0f}% of the reference in A/C/G/T bases)"
        )
        if len(rep_of) != retained:
            print(f"[*] {len(rep_of)} distinct of {retained} rows go to MAFFT; the copies "
                  f"are expanded back afterwards")
        if flipped:
            print(f"[*] {flipped} of {retained} sequences were on the reference's other "
                  f"strand and were reverse-complemented (MAFFT does not do this itself)")
        if flipped and flipped == retained:
            print(f"[warn] every sequence was on the other strand from {ref_rec.id}. That is "
                  f"expected when the reference is a CDS annotated on the minus strand; if the "
                  f"reference itself is the reversed one, its codon frame is wrong.",
                  file=sys.stderr)
        if unoriented:
            print(f"[warn] {unoriented} sequence(s) share as many {ORIENT_K}-mers with the "
                  f"reference's reverse complement as with the reference, so the strand could "
                  f"not be decided; aligned as given, and possibly the other way round.",
                  file=sys.stderr)

        # PM-26: the mapping, whenever a name changed or a row was flipped. Truncation is
        # gone, but renaming is not reversible either and nothing recorded it.
        if ids_changed or flipped:
            id_map = out_aln.with_name(out_aln.name + ".ids.tsv")
            tmp_map = tmp_dir / "ids.tsv"
            with open(tmp_map, "w") as fh:
                fh.write("aligned_id\toriginal_header\treverse_complemented\n")
                for clean, raw, rev in zip(ids, kept_raw, kept_flipped):
                    fh.write(f"{clean}\t{raw}\t{'yes' if rev else 'no'}\n")
            publish(tmp_map, id_map)
            print(f"[*] id map -> {id_map}")

        with open(clean_ref_p, "w") as fh:
            fh.write(f">Reference\n{str(ref_rec.seq).upper()}\n")

        print("[*] Running MAFFT --keeplength --addfragments...")
        mafft_cmd = [
            mafft_bin,
            "--auto",
            "--keeplength",
            # PM-11: without this MAFFT lower-cases its output, so the codon MSA was the one
            # lower-case file in the pipeline while the docs claimed rows came through
            # unchanged, and soft-masking tools saw a fully masked alignment.
            "--preservecase",
            "--thread", str(args.mafft_threads),
            "--addfragments", str(clean_input_p),
            str(clean_ref_p),
        ]

        with open(tmp_aln, "w") as fh:
            res = subprocess.run(mafft_cmd, stdout=fh, stderr=subprocess.PIPE, text=True)

        if res.returncode != 0:
            print(f"MAFFT failed: {res.stderr}", file=sys.stderr)
            return 1

        # Checked before the rename, so a failed check never replaces a good alignment.
        distinct = list(SeqIO.parse(tmp_aln, "fasta"))
        # PM-06: nothing used to check that MAFFT returned the records it was given, so a
        # run that aligned the reference to itself reported "MAFFT aligned 2 records" and
        # published a two-row alignment over a good one.
        by_id = {r.id: r for r in distinct}
        if (not distinct or distinct[0].id != "Reference"
                or len(distinct) != len(rep_of) + 1
                or any(rep not in by_id for _c, rep in rows)):
            print(f"Error: MAFFT returned {len(distinct)} records, expected the reference "
                  f"plus {len(rep_of)} sequences", file=sys.stderr)
            return 1

        # PM-10: expand the representatives back to every id, in input order. A new record
        # sharing the representative's `Seq` object, not a deep copy of it: `Seq` is
        # immutable, and copying one per row would put the whole alignment in memory twice.
        from Bio.SeqRecord import SeqRecord
        aligned_recs = [distinct[0]]
        for clean, rep in rows:
            aligned_recs.append(SeqRecord(by_id[rep].seq, id=clean, name=clean,
                                          description=""))

        print(f"[*] MAFFT aligned {len(rep_of)} distinct rows into {len(aligned_recs)} "
              f"records. Checking what --keeplength did to them...")
        for r in aligned_recs:
            if len(r.seq) != ref_len:
                print(
                    f"Error: aligned record {r.id} length {len(r.seq)} != ref length {ref_len}",
                    file=sys.stderr,
                )
                return 1

        # PM-25: coverage is measured again on the ALIGNED row. The only gate used to run on
        # the unaligned one, so a contig carrying part of the gene plus flanking sequence
        # passed on its raw base count and `--keeplength` then trimmed the flank away,
        # leaving a 20-25% fragment in the alignment and the tree with no message, although
        # the docs say such a row is dropped before alignment. A small tolerance, because
        # `--keeplength` also deletes real insertions, so a row just over the gate can fall
        # under it legitimately.
        floor = ref_len * max(0.0, args.msa_min_cov - 0.02)
        dropped = [r.id for r in aligned_recs[1:] if real_bases(str(r.seq)) < floor]
        if dropped:
            for rid in dropped[:10]:
                rec = next(r for r in aligned_recs if r.id == rid)
                print(f"[warn] {rid} covers {real_bases(str(rec.seq))/ref_len:.2f} of the "
                      f"reference after alignment; dropped", file=sys.stderr)
            if len(dropped) > 10:
                print(f"[warn] ... and {len(dropped) - 10} more", file=sys.stderr)
            aligned_recs = [aligned_recs[0]] + [r for r in aligned_recs[1:]
                                                if r.id not in set(dropped)]
            if len(aligned_recs) == 1:
                state["error"] = "every row fell under the coverage gate after alignment"
                print("Error: every row fell under the coverage gate after alignment",
                      file=sys.stderr)
                return 1

        # PM-03: say what --keeplength did. Warnings only; the alignment is unchanged.
        qc_rows, flags = [], {"frameshift": [], "split_codons": [], "internal_stops": [],
                              "trimmed": []}
        for r in aligned_recs[1:]:
            q = codon_qc(str(r.seq), input_acgt.get(r.id, 0))
            qc_rows.append((r.id, q))
            if q["frameshift"]:
                flags["frameshift"].append(r.id)
            if q["split_codons"]:
                flags["split_codons"].append(r.id)
            if q["internal_stops"]:
                flags["internal_stops"].append(r.id)
            if q["trimmed"]:
                flags["trimmed"].append(r.id)
        wording = {
            "frameshift": "carry an indel that is not a multiple of 3",
            "split_codons": "have a gap that splits a codon",
            "internal_stops": "have an internal stop codon",
            "trimmed": "lost bases to --keeplength (an insertion, or overhang past the "
                       "reference ends)",
        }
        for key, ids_flagged in flags.items():
            if ids_flagged:
                shown = ", ".join(ids_flagged[:5])
                more = f" and {len(ids_flagged) - 5} more" if len(ids_flagged) > 5 else ""
                print(f"[warn] {len(ids_flagged)} row(s) {wording[key]}: {shown}{more}",
                      file=sys.stderr)
        if any(flags.values()):
            qc_path = out_aln.with_name(out_aln.name + ".qc.tsv")
            tmp_qc = tmp_dir / "qc.tsv"
            with open(tmp_qc, "w") as fh:
                fh.write("id\ttrimmed_bases\tinternal_gap_runs\tsplit_codons\t"
                         "internal_stops\tframeshift\n")
                for rid, q in qc_rows:
                    fh.write(f"{rid}\t{q['trimmed']}\t"
                             f"{','.join(map(str, q['internal_gap_runs'])) or 'None'}\t"
                             f"{','.join(map(str, q['split_codons'])) or 'None'}\t"
                             f"{','.join(map(str, q['internal_stops'])) or 'None'}\t"
                             f"{'yes' if q['frameshift'] else 'no'}\n")
            publish(tmp_qc, qc_path)
            print(f"[*] codon QC -> {qc_path}")

        if not identity_ok(aligned_recs, args.min_identity):
            return 1
        # Rewrite only when the records are not what MAFFT wrote: expanded copies, or rows
        # dropped by the second coverage gate.
        if len(aligned_recs) != len(distinct):
            SeqIO.write(aligned_recs, tmp_aln, "fasta")
        publish(tmp_aln, out_aln)

        if not want_tree:
            print("[*] No tree requested.")
            state["taxa"] = len(aligned_recs)
            return 0

        # Written to a temp and published on success: opening out_tree directly truncated it
        # before FastTree ran, so a failure left a 0-byte .nwk that looks like a real output
        # to anything downstream (and destroyed a previously good tree). FastTree is fed from
        # the temp alignment, because a special `--out-aln` cannot be read back (PM-09).
        print(f"[*] Generating FastTree GTR phylogeny -> {out_tree}...")
        src_aln = tmp_aln if tmp_aln.is_file() else out_aln
        with open(src_aln, "r") as infh, open(tmp_tree, "w") as outfh:
            ft_res = subprocess.run(
                [fasttree_bin, "-nt", "-gtr", "-quiet"],
                stdin=infh, stdout=outfh, stderr=subprocess.PIPE, text=True,
            )
        if ft_res.returncode != 0:
            print(f"FastTree failed: {ft_res.stderr}", file=sys.stderr)
            return 1
        publish(tmp_tree, out_tree)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


    state["taxa"] = len(aligned_recs)
    state["flipped"] = flipped
    state["records"] = [
        {"id": r.id, "length": len(r.seq),
         "seq_sha256": summary_io.seq_sha256(str(r.seq))}
        for r in aligned_recs
    ]
    print(
        f"\n[✓] Codon alignment complete: {len(aligned_recs)} taxa x {ref_len} nt "
        f"({ref_len // 3} codons)"
    )
    print(f"    Alignment: {out_aln}")
    print(f"    Phylogeny: {out_tree}")
    return 0
