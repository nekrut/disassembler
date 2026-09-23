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
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_WALKER = ROOT / "walker" / "target" / "release" / "logan-walker"
DEFAULT_SIEVE = ROOT / "extern" / "kmer-diversity-sieve" / "l1_pipeline_serial"
CODON_MSA = ROOT / "python" / "build_walker_codon_msa.py"


def resolve_executable(name: str):
    """Return a runnable path for `name`, or None.

    `Path(name).is_file()` accepted a file without the execute bit, so the preflight passed and
    the run died later on an uncaught PermissionError from subprocess instead of the intended
    die(). `os.access(..., X_OK)` is the same test `shutil.which` applies to PATH entries.
    """
    p = Path(name)
    if p.is_file() and os.access(p, os.X_OK):
        return str(p)
    return shutil.which(name)


def first_record_len(path: Path):
    """Residues in the FIRST record of a FASTA, plus the number of records present.

    The walker aligns only record 1 (`walker/src/io/input.rs` warns and selects it), so the
    sieve has to be sized against that same record. Summing every record instead turned the
    160-record queries/master_planetary_targets.fasta into a --min-len of 297,645 rather than
    2,573, which rejects every haplotype and then fails blaming --sieve-threshold.
    """
    n = 0
    records = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                records += 1
                continue
            if records == 1:
                n += len(line.strip())
    return n, records


def run_step(cmd, what: str) -> None:
    """Run a stage, reporting a non-zero exit in this file's `die()` style.

    `check=True` raised CalledProcessError and dumped a full traceback plus the whole argv,
    which buries the tool's own message. It matters more now that logan-walker exits non-zero
    when nothing reconstructed, since that is a routine outcome for a bad accession list, not
    an internal error.
    """
    rc = subprocess.run(cmd).returncode
    if rc != 0:
        die(f"{what} failed (exit {rc}); see its output above.")


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
    # None for the same reason as the four below: these were a second copy of the walker's
    # defaults, passed on every run whatever the walker's own defaults became.
    ap.add_argument("-w", "--workers", type=int, default=None, help="walker compute workers (walker default 8)")
    ap.add_argument("--fetch-workers", type=int, default=None, help="walker S3 fetch workers (walker default 8)")
    ap.add_argument("--threads", type=int, default=8, help="MAFFT threads (default 8)")
    ap.add_argument("--min-cov", type=float, default=0.70, help="min coverage for a sequence to enter the MSA (default 0.70)")
    ap.add_argument("--emit-vcf", action="store_true", help="also write per-accession intra-host SNP VCFs")
    ap.add_argument("--emit-snv-table", action="store_true", help="also write the long snvs.tsv cohort table")
    ap.add_argument("--no-emit-alts", action="store_true", help="do not write alt_haplotypes.fasta")
    # These four default to None ON PURPOSE. They used to carry numeric defaults copied from
    # walker/src/config.rs, and because the guards below read `is not None` they were passed on
    # every run, so changing a walker default would have left this file silently pinning the old
    # one. With None the walker's own default applies unless the user asks for something else,
    # and the guards mean what they say. Each help string names the walker default for reference
    # rather than duplicating it as behaviour.
    ap.add_argument("--prune-frac", type=float, default=None, help="unitig depth prune fraction for reference bias removal (walker default 0.01)")
    ap.add_argument("--max-size-mb", type=int, default=None, help="maximum compressed SRA size in MB (walker default 1024)")
    ap.add_argument("--no-prefer-majority", action="store_true", help="leave primary on reference allele instead of majority allele")
    ap.add_argument("--min-abund", type=float, default=None, help="minimum ka:f abundance for a variant (walker default 5.0)")
    ap.add_argument("--min-vaf", type=float, default=None, help="minimum VAF for an alternative branch (walker default 0.10)")
    ap.add_argument("--walker-bin", default=str(DEFAULT_WALKER), help="path to the logan-walker binary")
    ap.add_argument("--skip-tree", action="store_true", help="stop after reconstruction; do not build the MSA/tree")
    ap.add_argument("--sieve", action="store_true",
                    help="OPTIONAL diversity-thinning: collapse near-duplicate reconstructed haplotypes "
                         "with the k-mer diversity sieve before the MSA/tree. HAPLOTYPE/LINEAGE-level "
                         "reducer (MinHash/LSH): it collapses sequences differing by a single SNV to ~1 "
                         "representative at any threshold, so use it only for LARGE cross-accession cohorts "
                         "of near-duplicate lineages, NOT to preserve intra-host SNVs (use --emit-vcf / "
                         "--emit-snv-table for those). Default OFF; when off, behavior is unchanged. "
                         "IGNORED under --skip-tree, since the sieve exists to thin the MSA/tree input.")
    ap.add_argument("--sieve-threshold", type=float, default=0.005,
                    help="sieve LSH promotion threshold tau (default 0.005). The nucleotide-divergence "
                         "boundary is d ~= tau / (2*k) at k=15, so 0.005 puts the boundary at ~0.017%%: "
                         "gene-scale, topology-preserving on the tested cohorts. The sieve binary's own "
                         "default of 0.15 targets archive-scale (~0.5%% divergence) and collapses a "
                         "gene-scale cohort to ~1 representative. Only used with --sieve; "
                         "see docs/issue_sieve_calibration.md.")
    ap.add_argument("--sieve-sketch", type=int, default=256,
                    help="sieve MinHash sketch size m (default 256). The estimator cannot resolve "
                         "novelty finer than 1/m, so 256 sets the floor near 0.015%% divergence. The "
                         "sieve binary's own default of 64 caps resolution at ~0.05%%, coarser than "
                         "the single-SNV separations a gene-scale cohort turns on. Only used with "
                         "--sieve; see docs/issue_sieve_calibration.md.")
    ap.add_argument("--sieve-min-len", type=int, default=None,
                    help="sieve minimum sequence length (default: auto = 95%% of the query CDS length, first "
                         "record only). The sieve's "
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
    if not resolve_executable(walker):
        die(f"logan-walker not found, or not executable, at {walker}. "
            f"Build it: (cd walker && cargo build --release), or pass --walker-bin.")
    # Both MSA-stage tools, not just mafft. FastTree used to be discovered only after the walk
    # and the alignment had already run, so a long cohort completed the expensive half and then
    # died on a raw traceback, leaving a 0-byte .nwk behind.
    # Silently ignoring a flag the user asked for is the one thing worse than not having it:
    # the run exits 0 having written the UNTHINNED reconstructions.
    if args.sieve and args.skip_tree:
        print("[disassemble] warn: --sieve is ignored under --skip-tree. The sieve thins the "
              "input to the MSA/tree stage, and that stage is being skipped, so "
              "reconstructed.fasta is written unthinned. Drop --skip-tree to get "
              "reconstructed.sieved.fasta.", file=sys.stderr)

    if not args.skip_tree:
        # Spellings and search order must match build_walker_codon_msa.find_bin, or this blocks
        # a run the MSA step would have accepted. bioconda ships FastTree under either name,
        # and find_bin falls back to $MAMBA_ROOT_PREFIX/envs/codon/bin when the tool is not on
        # PATH; skipping that fallback here would fail preflight on a machine where mafft and
        # fasttree live only in the codon env.
        codon_env_bin = Path(os.environ.get("MAMBA_ROOT_PREFIX",
                                            str(Path.home() / "micromamba"))) / "envs" / "codon" / "bin"
        for names, why in ((("mafft",), "the codon MSA"), (("fasttree", "FastTree"), "the ML tree")):
            on_path = any(shutil.which(n) for n in names)
            in_codon_env = any(os.access(codon_env_bin / n, os.X_OK) for n in names)
            if not (on_path or in_codon_env):
                die(f"{names[0]} not on PATH or in {codon_env_bin} (needed for {why}). "
                    f"Activate the `disassembler` env, or --skip-tree.")

    sieve_bin = args.sieve_bin
    if args.sieve and not resolve_executable(sieve_bin):
        die(f"kmer-diversity-sieve not found, or not executable, at {sieve_bin}. Build it: "
            f"(cd extern/kmer-diversity-sieve && make serial), or pass --sieve-bin.")

    recon = out / "reconstructed.fasta"

    # ---- stage 1: reconstruct ------------------------------------------------
    cmd = [walker, "-q", str(query), "-o", str(out), "-f", str(recon)]
    if args.workers is not None:
        cmd += ["-w", str(args.workers)]
    if args.fetch_workers is not None:
        cmd += ["--fetch-workers", str(args.fetch_workers)]
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
    run_step(cmd, "logan-walker (stage 1, reconstruct)")

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
        if args.sieve_min_len is not None:
            sieve_min_len = args.sieve_min_len
        else:
            # Record 1 only: that is the record the walker reconstructed.
            query_cds_len, n_records = first_record_len(query)
            if n_records > 1:
                print(f"[disassemble] warn: {query} has {n_records} records; sizing the sieve "
                      f"on the first one ({query_cds_len} nt), which is the record the walker "
                      f"aligns. One gene per query.", file=sys.stderr)
            if query_cds_len == 0:
                die(f"{query}: first record has no sequence; cannot derive --sieve-min-len.")
            sieve_min_len = int(query_cds_len * 0.95)
        sieved = out / "reconstructed.sieved.fasta"
        emit_prefix = out / "reconstructed.sieved"
        # clear any stale shards from a previous run so the concat glob can't pick them up
        for stale in out.glob("reconstructed.sieved.rank*.fasta"):
            stale.unlink()
        # Sketch is passed explicitly at the wrapper's default of 256, NOT the sieve binary's
        # default of 64. The binary's default is calibrated for archive-scale near-duplicate
        # collapse; gene-scale cohorts need the finer sketch to resolve their single-SNV
        # separations at all. See docs/issue_sieve_calibration.md.
        scmd = [sieve_bin, "--input", str(recon),
                "--min-len", str(sieve_min_len), "--max-n-frac", "0.05",
                "--threshold", str(args.sieve_threshold), "--warmup", "2000",
                "--sketch", str(args.sieve_sketch),
                "--emit-promoted", str(emit_prefix)]
        print(f"[disassemble] sieve: {' '.join(scmd)}", flush=True)
        run_step(scmd, "kmer-diversity-sieve")
        # the serial binary writes one shard per rank; concat PREFIX.rank*.fasta = promoted subset
        rank_files = sorted(out.glob("reconstructed.sieved.rank*.fasta"))
        with sieved.open("w") as fh:
            for rf in rank_files:
                fh.write(rf.read_text())
        n_sieved = sum(1 for line in sieved.open() if line.startswith(">")) if sieved.is_file() else 0
        print(f"[disassemble] sieve: {n_recon} -> {n_sieved} haplotypes "
              f"(tau={args.sieve_threshold}, sketch={args.sieve_sketch}) -> {sieved}", flush=True)
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
    run_step(cmd2, "build_walker_codon_msa.py (stage 2, MSA + tree)")

    print(f"[disassemble] done:\n"
          f"  haplotypes : {recon}\n"
          f"  codon MSA  : {aln}\n"
          f"  phylogeny  : {tree}")


if __name__ == "__main__":
    main()
