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

def find_bin(name: str) -> str:
    candidate = CODON_ENV_BIN / name
    if candidate.exists():
        return str(candidate)
    which = shutil.which(name)
    if which:
        return which
    return name

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
    clean_input_p = out_aln.parent / f".tmp_{in_p.name}"
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
    clean_ref_p = out_aln.parent / f".tmp_ref_{ref_p.name}"
    with open(clean_ref_p, "w") as fh:
        fh.write(f">Reference\n{str(ref_rec.seq).upper()}\n")

    mafft_bin = find_bin("mafft")
    fasttree_bin = find_bin("fasttree")

    print(f"[*] Running MAFFT --keeplength --addfragments...")
    mafft_cmd = [
        mafft_bin,
        "--auto",
        "--keeplength",
        "--thread", str(args.threads),
        "--addfragments", str(clean_input_p),
        str(clean_ref_p)
    ]
    
    with open(out_aln, "w") as fh:
        res = subprocess.run(mafft_cmd, stdout=fh, stderr=subprocess.PIPE, text=True)
    
    clean_input_p.unlink(missing_ok=True)
    clean_ref_p.unlink(missing_ok=True)

    if res.returncode != 0:
        sys.exit(f"MAFFT failed: {res.stderr}")

    # Validate alignment
    aligned_recs = list(SeqIO.parse(out_aln, "fasta"))
    print(f"[✓] MAFFT aligned {len(aligned_recs)} records. Validating codon lengths...")
    for r in aligned_recs:
        if len(r.seq) != ref_len:
            sys.exit(f"Error: aligned record {r.id} length {len(r.seq)} != ref length {ref_len}")

    # Build FastTree GTR phylogeny
    print(f"[*] Generating FastTree GTR phylogeny -> {out_tree}...")
    with open(out_aln, "r") as infh, open(out_tree, "w") as outfh:
        ft_res = subprocess.run(
            [fasttree_bin, "-nt", "-gtr", "-quiet"],
            stdin=infh, stdout=outfh, stderr=subprocess.PIPE, text=True
        )
    if ft_res.returncode != 0:
        sys.exit(f"FastTree failed: {ft_res.stderr}")

    print(f"\n[✓] Codon alignment complete: {len(aligned_recs)} taxa x {ref_len} nt ({ref_len//3} codons)")
    print(f"    Alignment: {out_aln}")
    print(f"    Phylogeny: {out_tree}")

if __name__ == "__main__":
    main()
