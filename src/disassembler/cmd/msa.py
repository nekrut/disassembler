"""msa: reference-guided codon MSA + FastTree phylogeny from reconstructed haplotypes.

Aligns each haplotype to the reference CDS with MAFFT `--keeplength --addfragments`
so reference coordinates 1..L are preserved and the reading frame stays invariant,
then builds a FastTree GTR tree from the alignment.

Ported verbatim from python/build_walker_codon_msa.py; behaviour is unchanged.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

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


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "msa",
        help="Codon MSA + FastTree GTR phylogeny from reconstructed haplotypes.",
        description=(
            "MAFFT --keeplength --addfragments against the reference CDS, then FastTree GTR."
        ),
    )
    p.add_argument("--ref", "-r", required=True, help="Reference in-frame CDS FASTA")
    p.add_argument(
        "--input", "-i", required=True,
        help="Reconstructed haplotypes FASTA (walker output)",
    )
    p.add_argument("--out-aln", "-a", required=True, help="Output codon MSA FASTA")
    p.add_argument("--out-tree", "-t", required=True, help="Output Newick tree")
    p.add_argument("--mafft-threads", type=int, default=8, help="MAFFT threads")
    p.add_argument("--msa-min-cov", type=float, default=0.70,
                   help="Minimum sequence coverage fraction (fraction of reference length "
                        "a sequence must reach to enter the MSA)")
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

    ref_records = list(SeqIO.parse(ref_p, "fasta"))
    if not ref_records:
        print(f"Error: no records found in reference FASTA {ref_p}", file=sys.stderr)
        return 1
    ref_rec = ref_records[0]
    ref_len = len(ref_rec.seq)
    if ref_len % 3 != 0:
        print(f"[warn] Reference CDS length {ref_len} is not divisible by 3!", file=sys.stderr)

    input_records = list(SeqIO.parse(in_p, "fasta"))
    if not input_records:
        print(f"Error: no records found in input FASTA {in_p}", file=sys.stderr)
        return 1

    print(
        f"[*] Processing {len(input_records)} reconstructed sequences against reference "
        f"{ref_rec.id} ({ref_len} bp)..."
    )

    used_ids = {"Reference"}
    clean_input_p = out_aln.parent / f".tmp_{in_p.name}.{os.getpid()}"
    retained = 0
    with open(clean_input_p, "w") as fh:
        for r in input_records:
            if len(r.seq) < ref_len * args.msa_min_cov:
                continue
            clean_id = sanitize_id(r.id, used_ids)
            fh.write(f">{clean_id}\n{str(r.seq).upper()}\n")
            retained += 1

    if retained == 0:
        clean_input_p.unlink(missing_ok=True)
        print("Error: no sequences passed quality/coverage filtering", file=sys.stderr)
        return 1

    print(
        f"[*] {retained}/{len(input_records)} sequences passed coverage filter "
        f"(>= {args.msa_min_cov*100:.0f}%)"
    )

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
        os.replace(tmp_aln, out_aln)
    finally:
        Path(tmp_aln).unlink(missing_ok=True)

    aligned_recs = list(SeqIO.parse(out_aln, "fasta"))
    print(f"[✓] MAFFT aligned {len(aligned_recs)} records. Validating codon lengths...")
    for r in aligned_recs:
        if len(r.seq) != ref_len:
            print(
                f"Error: aligned record {r.id} length {len(r.seq)} != ref length {ref_len}",
                file=sys.stderr,
            )
            return 1

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
