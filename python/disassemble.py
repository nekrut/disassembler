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
DEFAULT_SIEVE = ROOT / "extern" / "kmer-diversity-sieve" / "l1_pipeline_serial"
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
    ap.add_argument("--no-emit-alts", action="store_true", help="do not write alt_haplotypes.fasta")
    ap.add_argument("--prune-frac", type=float, default=0.01, help="unitig depth prune fraction for reference bias removal (default 0.01)")
    ap.add_argument("--max-size-mb", type=int, default=1024, help="maximum compressed SRA size in MB (default 1024)")
    ap.add_argument("--no-prefer-majority", action="store_true", help="leave primary on reference allele instead of majority allele")
    ap.add_argument("--min-abund", type=float, default=5.0, help="minimum ka:f abundance for a variant (default 5.0)")
    ap.add_argument("--min-vaf", type=float, default=0.10, help="minimum VAF for an alternative branch (default 0.10)")
    ap.add_argument("--walker-bin", default=str(DEFAULT_WALKER), help="path to the logan-walker binary")
    ap.add_argument("--skip-tree", action="store_true", help="stop after reconstruction; do not build the MSA/tree")
    ap.add_argument("--sieve", action="store_true",
                    help="OPTIONAL diversity-thinning: collapse near-duplicate reconstructed haplotypes "
                         "with the k-mer diversity sieve before the MSA/tree. HAPLOTYPE/LINEAGE-level "
                         "reducer (MinHash/LSH): it collapses sequences differing by a single SNV to ~1 "
                         "representative at any threshold, so use it only for LARGE cross-accession cohorts "
                         "of near-duplicate lineages, NOT to preserve intra-host SNVs (use --emit-vcf / "
                         "--emit-snv-table for those). Default OFF; when off, behavior is unchanged.")
    ap.add_argument("--sieve-threshold", type=float, default=0.15,
                    help="sieve LSH promotion threshold tau (default 0.15); only used with --sieve")
    ap.add_argument("--sieve-min-len", type=int, default=None,
                    help="sieve minimum sequence length (default: auto = 95%% of the query CDS length). The sieve's "
                         "own default is genome-scale (27000); gene-scale CDS need this ~= the query "
                         "length or everything is dropped. Only used with --sieve.")
    ap.add_argument("--sieve-bin", default=str(DEFAULT_SIEVE),
                    help="path to the kmer-diversity-sieve l1_pipeline_serial binary")
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

    sieve_bin = args.sieve_bin
    if args.sieve and not (Path(sieve_bin).is_file() or shutil.which(sieve_bin)):
        die(f"kmer-diversity-sieve not found at {sieve_bin}. Build it: "
            f"(cd extern/kmer-diversity-sieve && make serial), or pass --sieve-bin.")

    # length of the query CDS (sum of non-header residues) -> default sieve --min-len
    query_cds_len = sum(len(line.strip()) for line in query.open() if not line.startswith(">"))

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
    if args.no_emit_alts:
        cmd += ["--no-emit-alts"]
    if args.prune_frac is not None:
        cmd += ["--prune-frac", str(args.prune_frac)]
    if args.max_size_mb is not None:
        cmd += ["--max-size-mb", str(args.max_size_mb)]
    if args.no_prefer_majority:
        cmd += ["--no-prefer-majority"]
    if args.min_abund is not None:
        cmd += ["--min-abund", str(args.min_abund)]
    if args.min_vaf is not None:
        cmd += ["--min-vaf", str(args.min_vaf)]
    # drop a leading "--" separator that argparse.REMAINDER keeps
    extra = args.walker_args[1:] if args.walker_args[:1] == ["--"] else args.walker_args
    cmd += extra

    print(f"[disassemble] 1/2 reconstruct: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)

    n_recon = sum(1 for line in recon.open() if line.startswith(">")) if recon.is_file() else 0
    print(f"[disassemble] reconstructed {n_recon} haplotype(s) -> {recon}", flush=True)
    if n_recon == 0:
        die("no haplotypes reconstructed; nothing to align.")

    # ---- optional diversity-thinning sieve (default OFF) --------------------
    # When off, `msa_input` stays `recon` and everything downstream is byte-identical.
    # Skipped under --skip-tree: the sieve only exists to thin the MSA/tree input, so
    # there is nothing to thin for.
    msa_input = recon
    if args.sieve and not args.skip_tree:
        # Reconstructions are typically a few bases shorter than the reference (the walk
        # may not reach the exact termini), and the completeness gate is strict `< min-len`.
        # Defaulting to the exact query length would then drop every haplotype, so the
        # auto-default takes a 5% haircut. Pass --sieve-min-len to override.
        sieve_min_len = (args.sieve_min_len if args.sieve_min_len is not None
                         else int(query_cds_len * 0.95))
        sieved = out / "reconstructed.sieved.fasta"
        emit_prefix = out / "reconstructed.sieved"
        # clear any stale shards from a previous run so the concat glob can't pick them up
        for stale in out.glob("reconstructed.sieved.rank*.fasta"):
            stale.unlink()
        scmd = [sieve_bin, "--input", str(recon),
                "--min-len", str(sieve_min_len), "--max-n-frac", "0.05",
                "--threshold", str(args.sieve_threshold), "--warmup", "2000",
                "--emit-promoted", str(emit_prefix)]
        print(f"[disassemble] sieve: {' '.join(scmd)}", flush=True)
        subprocess.run(scmd, check=True)
        # the serial binary writes one shard per rank; concat PREFIX.rank*.fasta = promoted subset
        rank_files = sorted(out.glob("reconstructed.sieved.rank*.fasta"))
        with sieved.open("w") as fh:
            for rf in rank_files:
                fh.write(rf.read_text())
        n_sieved = sum(1 for line in sieved.open() if line.startswith(">")) if sieved.is_file() else 0
        print(f"[disassemble] sieve: {n_recon} -> {n_sieved} haplotypes "
              f"(tau={args.sieve_threshold}) -> {sieved}", flush=True)
        if n_sieved == 0:
            die("sieve promoted 0 haplotypes; try a lower --sieve-threshold or "
                "a lower --sieve-min-len (reconstructions may be shorter than the reference).")
        msa_input = sieved

    if args.skip_tree:
        print("[disassemble] --skip-tree: done after reconstruction.")
        return

    # ---- stage 2: codon MSA + FastTree --------------------------------------
    stem = query.stem
    aln = out / f"{stem}_codon_msa.fasta"
    tree = out / f"{stem}_fasttree.nwk"
    cmd2 = [sys.executable, str(CODON_MSA), "--ref", str(query), "--input", str(msa_input),
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
