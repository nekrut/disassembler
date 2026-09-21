#!/usr/bin/env python3
"""
build_walker_codon_msa.py: Build reference-guided codon alignment and FastTree phylogeny
from Logan path-walker reconstructed haplotypes.

Features:
1. Aligns reconstructed haplotypes to reference CDS using MAFFT --keeplength --addfragments.
   This guarantees that reference coordinates 1..L are strictly preserved and insertions
   relative to reference are trimmed, ensuring an invariant reading frame.
2. Sanitizes headers to alphanumeric/underscores for downstream HyphAeon / HYPHY compatibility.
3. Ensures reference sequence is sequence #1.
4. Generates FastTree GTR maximum-likelihood phylogeny.
5. Emits alignment QC stats (number of sequences, length, codon count, pairwise identity).
"""

import sys
import os
import re
import argparse
import subprocess
import shutil
from pathlib import Path
from Bio import SeqIO

MAMBA_PREFIX = Path(os.environ.get("MAMBA_ROOT_PREFIX", str(Path.home() / "micromamba")))
CODON_ENV_BIN = MAMBA_PREFIX / "envs" / "codon" / "bin"


def find_bin(name: str, *aliases: str) -> str:
    """Resolve a tool to an executable path. PATH wins; the `codon` env is a fallback.

    `aliases` are alternative spellings of the same tool, tried in order. bioconda ships
    FastTree under both `fasttree` and `FastTree` depending on the build, and resolving now
    rather than at exec time means a wrong guess would block a run that used to work.

    This used to check $MAMBA_ROOT_PREFIX/envs/codon/bin *before* PATH, a leftover from the
    upstream layout. On any machine that happens to have a micromamba env named `codon`, MAFFT
    and FastTree came from there while the caller had preflighted the ones on PATH, so
    environment.yml's pins stopped describing what actually ran and could differ silently from
    the versions the golden fixtures were captured with. The fallback is kept, because it is
    what makes the tool work on a machine where these live only in that env, but it now
    announces itself.

    Raises rather than returning the bare name: an unresolvable tool has to fail here, with a
    message naming it, instead of surfacing as a FileNotFoundError traceback from subprocess
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
            print(f"[*] {n}: not on PATH, falling back to {fallback} "
                  f"(environment.yml pins do not describe this binary)", file=sys.stderr)
            return fallback
    spellings = " or ".join(repr(n) for n in names)
    raise FileNotFoundError(
        f"{spellings} not found on PATH or in {CODON_ENV_BIN}. "
        f"Activate the `disassembler` env (environment.yml pins mafft and fasttree), "
        f"or put {name} on PATH."
    )

def sanitize_id(s: str, used: set) -> str:
    clean = re.sub(r'[^A-Za-z0-9_]', '_', s).strip('_')
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

def main():
    parser = argparse.ArgumentParser(description="Build codon alignment and phylogeny from reconstructed haplotypes")
    parser.add_argument("--ref", "-r", required=True, help="Reference in-frame CDS FASTA")
    parser.add_argument("--input", "-i", required=True, help="Reconstructed haplotypes FASTA from logan_path_walker")
    parser.add_argument("--out-aln", "-a", required=True, help="Output codon MSA FASTA")
    parser.add_argument("--out-tree", "-t", required=True, help="Output Newick tree")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads for MAFFT")
    parser.add_argument("--min-cov", type=float, default=0.70, help="Minimum sequence coverage fraction")
    args = parser.parse_args()

    ref_p = Path(args.ref)
    in_p = Path(args.input)
    out_aln = Path(args.out_aln)
    out_tree = Path(args.out_tree)
    out_aln.parent.mkdir(parents=True, exist_ok=True)
    out_tree.parent.mkdir(parents=True, exist_ok=True)

    # Resolve both tools before doing anything expensive. FastTree used to be resolved after
    # the alignment, so a machine with mafft but no fasttree burned the whole walk and MSA
    # before dying on an unhandled traceback.
    try:
        mafft_bin = find_bin("mafft")
        fasttree_bin = find_bin("fasttree", "FastTree")
    except FileNotFoundError as e:
        sys.exit(f"Error: {e}")

    ref_records = list(SeqIO.parse(ref_p, "fasta"))
    if not ref_records:
        sys.exit(f"Error: no records found in reference FASTA {ref_p}")
    ref_rec = ref_records[0]
    ref_len = len(ref_rec.seq)
    if ref_len % 3 != 0:
        print(f"[warn] Reference CDS length {ref_len} is not divisible by 3!", file=sys.stderr)

    input_records = list(SeqIO.parse(in_p, "fasta"))
    if not input_records:
        sys.exit(f"Error: no records found in input FASTA {in_p}")

    print(f"[*] Processing {len(input_records)} reconstructed sequences against reference {ref_rec.id} ({ref_len} bp)...")

    # Filter and sanitize input sequences
    used_ids = { "Reference" }
    clean_input_p = out_aln.parent / f".tmp_{in_p.name}.{os.getpid()}"
    retained = 0
    with open(clean_input_p, "w") as fh:
        for r in input_records:
            # Check length fraction
            if len(r.seq) < ref_len * args.min_cov:
                continue
            clean_id = sanitize_id(r.id, used_ids)
            fh.write(f">{clean_id}\n{str(r.seq).upper()}\n")
            retained += 1

    if retained == 0:
        clean_input_p.unlink(missing_ok=True)
        sys.exit("Error: no sequences passed quality/coverage filtering")

    print(f"[*] {retained}/{len(input_records)} sequences passed coverage filter (>= {args.min_cov*100:.0f}%)")

    # Create temporary single-record reference FASTA
    clean_ref_p = out_aln.parent / f".tmp_ref_{ref_p.name}.{os.getpid()}"
    with open(clean_ref_p, "w") as fh:
        fh.write(f">Reference\n{str(ref_rec.seq).upper()}\n")

    print(f"[*] Running MAFFT --keeplength --addfragments...")
    mafft_cmd = [
        mafft_bin,
        "--auto",
        "--keeplength",
        "--thread", str(args.threads),
        "--addfragments", str(clean_input_p),
        str(clean_ref_p)
    ]
    
    # Same temp-then-rename as the tree below: opening out_aln directly truncated it before
    # MAFFT ran, so a MAFFT failure replaced a previously good alignment with an empty file.
    tmp_aln = out_aln.with_name(f".tmp_{out_aln.name}.{os.getpid()}")
    try:
        with open(tmp_aln, "w") as fh:
            res = subprocess.run(mafft_cmd, stdout=fh, stderr=subprocess.PIPE, text=True)

        clean_input_p.unlink(missing_ok=True)
        clean_ref_p.unlink(missing_ok=True)

        if res.returncode != 0:
            sys.exit(f"MAFFT failed: {res.stderr}")
        os.replace(tmp_aln, out_aln)
    finally:
        Path(tmp_aln).unlink(missing_ok=True)

    # Validate alignment
    aligned_recs = list(SeqIO.parse(out_aln, "fasta"))
    print(f"[✓] MAFFT aligned {len(aligned_recs)} records. Validating codon lengths...")
    for r in aligned_recs:
        if len(r.seq) != ref_len:
            sys.exit(f"Error: aligned record {r.id} length {len(r.seq)} != ref length {ref_len}")

    # Build FastTree GTR phylogeny.
    #
    # Written to a temp and renamed on success: opening out_tree directly truncated it before
    # FastTree ran, so a failure left a 0-byte .nwk that looks like a real output to anything
    # downstream (and destroyed a previously good tree).
    print(f"[*] Generating FastTree GTR phylogeny -> {out_tree}...")
    tmp_tree = out_tree.with_name(f".tmp_{out_tree.name}.{os.getpid()}")
    try:
        with open(out_aln, "r") as infh, open(tmp_tree, "w") as outfh:
            ft_res = subprocess.run(
                [fasttree_bin, "-nt", "-gtr", "-quiet"],
                stdin=infh, stdout=outfh, stderr=subprocess.PIPE, text=True
            )
        if ft_res.returncode != 0:
            sys.exit(f"FastTree failed: {ft_res.stderr}")
        os.replace(tmp_tree, out_tree)
    finally:
        Path(tmp_tree).unlink(missing_ok=True)

    print(f"\n[✓] Codon alignment complete: {len(aligned_recs)} taxa x {ref_len} nt ({ref_len//3} codons)")
    print(f"    Alignment: {out_aln}")
    print(f"    Phylogeny: {out_tree}")

if __name__ == "__main__":
    main()
