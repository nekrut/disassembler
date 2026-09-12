#!/usr/bin/env python3
"""disassemble: accessions + a query gene, to a phylogeny, in one command.

Chains the two halves of this repo end to end:

  1. `walker/` (logan-walker)          stream Logan unitigs, reconstruct the coding haplotype
  2. `build_walker_codon_msa.py`       reference-guided in-frame codon MSA + FastTree GTR tree

so you go straight from an SRA accession list to `reconstructed.fasta`, a codon alignment and a
Newick tree without running the two stages by hand.

Needs, all on PATH (the `disassembler` conda env from environment.yml provides them):
  - GraphAligner            (the walker shells out to it)
  - MAFFT, FastTree         (the codon MSA / tree step)
  - biopython               (importable by the interpreter running this script)

Examples
--------
  # one accession
  python/disassemble.py -a SRR6251967 -q queries/markers/zika_ns5.fasta -o out/zika

  # a cohort, with the intra-host SNP table, 8 workers
  python/disassemble.py -l accessions.txt -q queries/markers/zika_ns5.fasta -o out/zika \\
      --emit-snv-table -w 8

Anything after a literal `--` is passed straight through to logan-walker, e.g.
  python/disassemble.py -a SRR8209082 -q q.fasta -o out -- --prune-frac 0.0 --no-prefer-majority
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_WALKER = ROOT / "walker" / "target" / "release" / "logan-walker"
CODON_MSA = ROOT / "python" / "build_walker_codon_msa.py"


def die(msg: str) -> None:
    sys.exit(f"disassemble: {msg}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="End-to-end: accessions + query gene -> haplotypes -> codon MSA -> FastTree.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Pass extra logan-walker flags after a literal `--`.",
    )
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("-a", "--accession", help="a single SRA accession")
    src.add_argument("-l", "--accession-list", help="file with one accession per line")
    ap.add_argument("-q", "--query", required=True, help="reference in-frame CDS FASTA")
    ap.add_argument("-o", "--out-dir", required=True, help="output directory")
    ap.add_argument("-w", "--workers", type=int, default=8, help="walker compute workers (default 8)")
    ap.add_argument("--fetch-workers", type=int, default=8, help="walker S3 fetch workers (default 8)")
    ap.add_argument("--threads", type=int, default=8, help="MAFFT threads (default 8)")
    ap.add_argument("--min-cov", type=float, default=0.70, help="min coverage for a sequence to enter the MSA (default 0.70)")
    ap.add_argument("--emit-vcf", action="store_true", help="also write per-accession intra-host SNP VCFs")
    ap.add_argument("--emit-snv-table", action="store_true", help="also write the long snvs.tsv cohort table")
    ap.add_argument("--walker-bin", default=str(DEFAULT_WALKER), help="path to the logan-walker binary")
    ap.add_argument("--skip-tree", action="store_true", help="stop after reconstruction; do not build the MSA/tree")
    ap.add_argument("walker_args", nargs=argparse.REMAINDER,
                    help="extra args passed to logan-walker after `--`")
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    query = Path(args.query)
    if not query.is_file():
        die(f"query not found: {query}")

    walker = args.walker_bin
    if not (Path(walker).is_file() or shutil.which(walker)):
        die(f"logan-walker not found at {walker}. Build it: (cd walker && cargo build --release), "
            f"or pass --walker-bin.")
    if not args.skip_tree and shutil.which("mafft") is None:
        die("mafft not on PATH (needed for the codon MSA). Activate the `disassembler` env, or --skip-tree.")

    recon = out / "reconstructed.fasta"

    # ---- stage 1: reconstruct ------------------------------------------------
    cmd = [walker, "-q", str(query), "-o", str(out), "-f", str(recon),
           "-w", str(args.workers), "--fetch-workers", str(args.fetch_workers)]
    if args.accession:
        cmd += ["-a", args.accession]
    else:
        cmd += ["-l", args.accession_list]
    if args.emit_vcf:
        cmd += ["--emit-vcf"]
    if args.emit_snv_table:
        cmd += ["--emit-snv-table"]
    # drop a leading "--" separator that argparse.REMAINDER keeps
    extra = args.walker_args[1:] if args.walker_args[:1] == ["--"] else args.walker_args
    cmd += extra

    print(f"[disassemble] 1/2 reconstruct: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)

    n_recon = sum(1 for line in recon.open() if line.startswith(">")) if recon.is_file() else 0
    print(f"[disassemble] reconstructed {n_recon} haplotype(s) -> {recon}", flush=True)
    if n_recon == 0:
        die("no haplotypes reconstructed; nothing to align.")
    if args.skip_tree:
        print("[disassemble] --skip-tree: done after reconstruction.")
        return

    # ---- stage 2: codon MSA + FastTree --------------------------------------
    stem = query.stem
    aln = out / f"{stem}_codon_msa.fasta"
    tree = out / f"{stem}_fasttree.nwk"
    cmd2 = [sys.executable, str(CODON_MSA), "--ref", str(query), "--input", str(recon),
            "--out-aln", str(aln), "--out-tree", str(tree),
            "--threads", str(args.threads), "--min-cov", str(args.min_cov)]
    print(f"[disassemble] 2/2 codon MSA + tree: {' '.join(cmd2)}", flush=True)
    subprocess.run(cmd2, check=True)

    print(f"[disassemble] done:\n"
          f"  haplotypes : {recon}\n"
          f"  codon MSA  : {aln}\n"
          f"  phylogeny  : {tree}")


if __name__ == "__main__":
    main()
