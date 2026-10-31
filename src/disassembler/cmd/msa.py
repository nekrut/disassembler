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
import os
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

from disassembler.textio import open_text
from disassembler.validators import fraction, positive_int

MAMBA_PREFIX = Path(os.environ.get("MAMBA_ROOT_PREFIX", str(Path.home() / "micromamba")))
CODON_ENV_BIN = MAMBA_PREFIX / "envs" / "codon" / "bin"


def find_bin(name: str, *aliases: str) -> str:
    """Resolve a tool to an executable path. PATH wins; the `codon` env is a fallback.

    `aliases` are alternative spellings of the same tool, tried in order. bioconda ships
    FastTree under both `fasttree` and `FastTree` depending on the build.

    Raises rather than returning the bare name so an unresolvable tool fails here with
    a message naming it, instead of surfacing as a FileNotFoundError from subprocess
    after the expensive stages have already run.
    """
    names = (name, *aliases)
    for n in names:
        found = shutil.which(n)
        if found:
            return found
    for n in names:
        fallback = shutil.which(str(CODON_ENV_BIN / n))
        if fallback:
            print(
                f"[*] {n}: not on PATH, falling back to {fallback} "
                f"(environment.yml pins do not describe this binary)",
                file=sys.stderr,
            )
            return fallback
    spellings = " or ".join(repr(n) for n in names)
    raise FileNotFoundError(
        f"{spellings} not found on PATH or in {CODON_ENV_BIN}. "
        f"Activate the `disassembler` env (environment.yml pins mafft and fasttree), "
        f"or put {name} on PATH."
    )


def sanitize_id(s: str, used: set) -> str:
    clean = re.sub(r"[^A-Za-z0-9_]", "_", s).strip("_")
    if not clean:
        clean = "seq"
    base = clean[:40]
    candidate = base
    idx = 1
    while candidate in used:
        candidate = f"{base}_{idx}"
        idx += 1
    used.add(candidate)
    return candidate


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


def kmer_set(seq: str, step: int = 1) -> set:
    """Upper-case k-mers of `seq` made only of A/C/G/T, taking every `step`-th start."""
    seq = seq.upper()
    acgt = set("ACGT")
    return {seq[i:i + ORIENT_K] for i in range(0, len(seq) - ORIENT_K + 1, step)
            if set(seq[i:i + ORIENT_K]) <= acgt}


def orientation(seq: str, fwd: set, rev: set) -> int:
    """+1 if `seq` shares more k-mers with the reference than with its reverse complement,
    -1 if fewer, 0 if as many (including none at all)."""
    kmers = kmer_set(seq, ORIENT_STEP)
    n_fwd, n_rev = len(kmers & fwd), len(kmers & rev)
    return (n_fwd > n_rev) - (n_rev > n_fwd)


def identity_to_reference(ref: str, row: str) -> float:
    """Matches over the alignment columns where both the reference and the row have A/C/G/T."""
    compared = matches = 0
    for a, b in zip(ref.upper(), row.upper()):
        if a in "ACGT" and b in "ACGT":
            compared += 1
            matches += a == b
    return matches / compared if compared else 0.0


def identity_ok(aligned_recs, min_identity: float) -> bool:
    """Print identity to the reference (the first record); False if the median is too low.

    `--keeplength` fits any input to the reference's length, so the length check above can
    not catch an alignment that is noise. Identity can: a strand mismatch measured 0.40-0.50.
    """
    ref = str(aligned_recs[0].seq)
    ids = [(r.id, identity_to_reference(ref, str(r.seq))) for r in aligned_recs[1:]]
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
    p.add_argument("--out-tree", "-t", required=True, help="Output Newick tree")
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


def run(args) -> int:
    from Bio import SeqIO  # imported lazily so `disassembler --help` does not need biopython

    ref_p = Path(args.ref)
    in_p = Path(args.input)
    out_aln = Path(args.out_aln)
    out_tree = Path(args.out_tree)
    out_aln.parent.mkdir(parents=True, exist_ok=True)
    out_tree.parent.mkdir(parents=True, exist_ok=True)

    try:
        mafft_bin = find_bin("mafft")
        fasttree_bin = find_bin("fasttree", "FastTree")
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # A missing, zstd-compressed or binary input ends in one line, not a traceback.
    try:
        with open_text(ref_p) as fh:
            ref_records = list(SeqIO.parse(fh, "fasta"))
        with open_text(in_p) as fh:
            input_records = list(SeqIO.parse(fh, "fasta"))
    except (OSError, ValueError) as e:
        print(f"Error: cannot read the reference or input FASTA: {e}", file=sys.stderr)
        return 1
    if not ref_records:
        print(f"Error: no records found in reference FASTA {ref_p}", file=sys.stderr)
        return 1
    ref_rec = ref_records[0]
    ref_len = len(ref_rec.seq)
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

    used_ids = {"Reference"}
    clean_input_p = out_aln.parent / f".tmp_{in_p.name}.{os.getpid()}"
    retained = flipped = unoriented = 0
    with open(clean_input_p, "w") as fh:
        for r in input_records:
            seq = str(r.seq).upper()
            # Count bases, not characters: a lexicmap row is reference-length with '-' where
            # no HSP reached, so its length passed this gate whatever its coverage.
            if real_bases(seq) < ref_len * args.msa_min_cov:
                continue
            strand = orientation(seq, fwd_kmers, rev_kmers)
            if strand < 0:
                seq = revcomp(seq)
                flipped += 1
            elif strand == 0:
                unoriented += 1
            clean_id = sanitize_id(r.id, used_ids)
            fh.write(f">{clean_id}\n{seq}\n")
            retained += 1

    if retained == 0:
        clean_input_p.unlink(missing_ok=True)
        print("Error: no sequences passed quality/coverage filtering", file=sys.stderr)
        return 1

    print(
        f"[*] {retained}/{len(input_records)} sequences passed coverage filter "
        f"(>= {args.msa_min_cov*100:.0f}% of the reference in A/C/G/T bases)"
    )
    if flipped:
        print(f"[*] {flipped} of {retained} sequences were on the reference's other strand and "
              f"were reverse-complemented (MAFFT does not do this itself)")
    if flipped and flipped == retained:
        print(f"[warn] every sequence was on the other strand from {ref_rec.id}. That is "
              f"expected when the reference is a CDS annotated on the minus strand; if the "
              f"reference itself is the reversed one, its codon frame is wrong.",
              file=sys.stderr)
    if unoriented:
        print(f"[warn] {unoriented} sequence(s) share as many {ORIENT_K}-mers with the "
              f"reference's reverse complement as with the reference (usually none with "
              f"either); aligned as given.", file=sys.stderr)

    clean_ref_p = out_aln.parent / f".tmp_ref_{ref_p.name}.{os.getpid()}"
    with open(clean_ref_p, "w") as fh:
        fh.write(f">Reference\n{str(ref_rec.seq).upper()}\n")

    print("[*] Running MAFFT --keeplength --addfragments...")
    mafft_cmd = [
        mafft_bin,
        "--auto",
        "--keeplength",
        "--thread", str(args.mafft_threads),
        "--addfragments", str(clean_input_p),
        str(clean_ref_p),
    ]

    # Same temp-then-rename as the tree below: opening out_aln directly truncated it
    # before MAFFT ran, so a MAFFT failure replaced a previously good alignment with
    # an empty file.
    tmp_aln = out_aln.with_name(f".tmp_{out_aln.name}.{os.getpid()}")
    try:
        with open(tmp_aln, "w") as fh:
            res = subprocess.run(mafft_cmd, stdout=fh, stderr=subprocess.PIPE, text=True)

        clean_input_p.unlink(missing_ok=True)
        clean_ref_p.unlink(missing_ok=True)

        if res.returncode != 0:
            print(f"MAFFT failed: {res.stderr}", file=sys.stderr)
            return 1

        # Checked before the rename, so a failed check never replaces a good alignment.
        aligned_recs = list(SeqIO.parse(tmp_aln, "fasta"))
        print(f"[✓] MAFFT aligned {len(aligned_recs)} records. Validating codon lengths...")
        for r in aligned_recs:
            if len(r.seq) != ref_len:
                print(
                    f"Error: aligned record {r.id} length {len(r.seq)} != ref length {ref_len}",
                    file=sys.stderr,
                )
                return 1
        if not identity_ok(aligned_recs, args.min_identity):
            return 1
        os.replace(tmp_aln, out_aln)
    finally:
        Path(tmp_aln).unlink(missing_ok=True)

    # Written to a temp and renamed on success: opening out_tree directly truncated it
    # before FastTree ran, so a failure left a 0-byte .nwk that looks like a real output
    # to anything downstream (and destroyed a previously good tree).
    print(f"[*] Generating FastTree GTR phylogeny -> {out_tree}...")
    tmp_tree = out_tree.with_name(f".tmp_{out_tree.name}.{os.getpid()}")
    try:
        with open(out_aln, "r") as infh, open(tmp_tree, "w") as outfh:
            ft_res = subprocess.run(
                [fasttree_bin, "-nt", "-gtr", "-quiet"],
                stdin=infh, stdout=outfh, stderr=subprocess.PIPE, text=True,
            )
        if ft_res.returncode != 0:
            print(f"FastTree failed: {ft_res.stderr}", file=sys.stderr)
            return 1
        os.replace(tmp_tree, out_tree)
    finally:
        Path(tmp_tree).unlink(missing_ok=True)

    print(
        f"\n[✓] Codon alignment complete: {len(aligned_recs)} taxa x {ref_len} nt "
        f"({ref_len // 3} codons)"
    )
    print(f"    Alignment: {out_aln}")
    print(f"    Phylogeny: {out_tree}")
    return 0
