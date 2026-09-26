"""reconstruct: accessions + query CDS -> haplotypes -> (sieve, opt) -> codon MSA + tree.

Composite: walker -> sieve (optional) -> msa. Equivalent to the old
`python/disassemble.py`. Anything after a literal `--` is passed straight through
to logan-walker.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import types
from pathlib import Path

from disassembler.cmd import msa as msa_cmd
from disassembler.cmd import sieve as sieve_cmd
from disassembler.cmd import walker as walker_cmd
from disassembler.resolve import resolve_sieve, resolve_walker


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "reconstruct",
        help="Accessions -> haplotypes -> codon MSA + tree (walker + optional sieve + msa).",
        description=(
            "End-to-end reconstruct: logan-walker, optional k-mer diversity sieve, "
            "codon MSA + FastTree GTR tree. Anything after `--` is forwarded to "
            "logan-walker verbatim."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Pass extra logan-walker flags after a literal `--`.",
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("-a", "--accession", help="A single SRA accession")
    src.add_argument("-l", "--accession-list", help="File with one accession per line")
    p.add_argument("-q", "--query", required=True, action="append",
                   help="Reference in-frame CDS FASTA. Repeatable, and each file may hold "
                        "multiple records; the walker unions the anchor pool when the "
                        "collected reference set has more than one record. In multi-anchor "
                        "mode the MSA step still needs one reference; that is the first -q "
                        "unless --msa-reference is given.")
    p.add_argument("--msa-reference", default=None,
                   help="FASTA to use as the MSA reference in multi-anchor mode "
                        "(default: first -q). Ignored under --skip-tree.")
    p.add_argument("-o", "--out-dir", required=True, help="Output directory")

    # None on purpose: numeric defaults would be a second copy of the walker's own,
    # silently pinning old values when the walker's default changed. With None the
    # flag reaches the walker only when the user asked for it.
    p.add_argument("-w", "--workers", type=int, default=None,
                   help="walker compute workers (walker default 8)")
    p.add_argument("--fetch-workers", type=int, default=None,
                   help="walker S3 fetch workers (walker default 8)")
    p.add_argument("--mafft-threads", type=int, default=8,
                   help="MAFFT threads (default 8). Distinct from walker's -t/--threads, "
                        "which sets GraphAligner's thread count and is kept at 1 for "
                        "reproducibility")
    p.add_argument("--msa-min-cov", type=float, default=0.70,
                   help="Min coverage for a sequence to enter the MSA (default 0.70). "
                        "Distinct from walker's --min-cov (default 0.85), which gates which "
                        "reconstructions are produced at all")
    p.add_argument("--emit-vcf", action="store_true",
                   help="also write per-accession intra-host SNP VCFs")
    p.add_argument("--emit-snv-table", action="store_true",
                   help="also write the long snvs.tsv cohort table")
    p.add_argument("--no-emit-alts", action="store_true",
                   help="do not write alt_haplotypes.fasta")
    p.add_argument("--prune-frac", type=float, default=None,
                   help="unitig depth prune fraction for reference bias removal "
                        "(walker default 0.01)")
    p.add_argument("--max-size-mb", type=int, default=None,
                   help="maximum compressed SRA size in MB (walker default 1024)")
    p.add_argument("--no-prefer-majority", action="store_true",
                   help="leave primary on reference allele instead of majority allele")
    p.add_argument("--min-abund", type=float, default=None,
                   help="minimum ka:f abundance for a variant (walker default 5.0)")
    p.add_argument("--min-vaf", type=float, default=None,
                   help="minimum VAF for an alternative branch (walker default 0.10)")
    p.add_argument("--walker-bin", default=None, help="path to the logan-walker binary")
    p.add_argument("--skip-tree", action="store_true",
                   help="stop after reconstruction; do not build the MSA/tree")
    p.add_argument("--sieve", action="store_true",
                   help="OPTIONAL diversity-thinning: collapse near-duplicate reconstructed "
                        "haplotypes with the k-mer diversity sieve before the MSA/tree. "
                        "HAPLOTYPE/LINEAGE-level reducer (MinHash/LSH); use only for LARGE "
                        "cross-accession cohorts of near-duplicate lineages, NOT to preserve "
                        "intra-host SNVs (use --emit-vcf / --emit-snv-table for those). "
                        "IGNORED under --skip-tree.")
    p.add_argument("--sieve-threshold", type=float, default=0.005,
                   help="sieve LSH promotion threshold tau (default 0.005, gene-scale). "
                        "Only used with --sieve.")
    p.add_argument("--sieve-sketch", type=int, default=256,
                   help="sieve MinHash sketch size m (default 256, gene-scale). Only used "
                        "with --sieve.")
    p.add_argument("--sieve-min-len", type=int, default=None,
                   help="sieve minimum sequence length (default: auto = 95%% of the query "
                        "CDS length, first record only). Only used with --sieve.")
    p.add_argument("--sieve-bin", default=None,
                   help="path to the kmer-diversity-sieve l1_pipeline_serial binary")
    p.add_argument("--unitigs-cache", default=None, metavar="DIR",
                   help="Read-through unitigs cache passed to walker. Enables offline reruns.")
    p.add_argument("--offline", action="store_true",
                   help="Refuse the S3 fallback for the walker stage. Requires --unitigs-cache.")
    p.add_argument("walker_args", nargs=argparse.REMAINDER,
                   help="extra args passed to logan-walker after `--`")
    p.set_defaults(func=run)


def _die(msg: str) -> int:
    print(f"reconstruct: {msg}", file=sys.stderr)
    return 1


def run(args) -> int:
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    queries = [Path(q) for q in args.query]
    for q in queries:
        if not q.is_file():
            return _die(f"query not found: {q}")
    # MSA/tree stage still takes a single reference. Default is the first -q, override
    # via --msa-reference; irrelevant under --skip-tree.
    msa_ref = Path(args.msa_reference) if args.msa_reference else queries[0]
    if not args.skip_tree and not msa_ref.is_file():
        return _die(f"MSA reference not found: {msa_ref}")
    if len(queries) > 1 and not args.skip_tree and args.msa_reference is None:
        print(f"[reconstruct] multi-anchor: MSA reference defaults to first -q ({msa_ref}); "
              f"pass --msa-reference to override.", flush=True)

    if not resolve_walker(args.walker_bin):
        return _die(
            f"logan-walker not found. Build it (`cd walker && cargo build --release`), "
            f"put it on PATH, or pass --walker-bin."
        )

    if args.sieve and args.skip_tree:
        print("[reconstruct] warn: --sieve is ignored under --skip-tree. The sieve thins "
              "the input to the MSA/tree stage, and that stage is being skipped, so "
              "reconstructed.fasta is written unthinned. Drop --skip-tree to get "
              "reconstructed.sieved.fasta.", file=sys.stderr)

    if not args.skip_tree:
        # Preflight MAFFT and FastTree before the expensive walk stage. Same spellings
        # and fallback as `msa.find_bin`.
        codon_env_bin = Path(os.environ.get("MAMBA_ROOT_PREFIX",
                                            str(Path.home() / "micromamba"))) / "envs" / "codon" / "bin"
        for names, why in ((("mafft",), "the codon MSA"),
                           (("fasttree", "FastTree"), "the ML tree")):
            on_path = any(shutil.which(n) for n in names)
            in_codon_env = any(os.access(codon_env_bin / n, os.X_OK) for n in names)
            if not (on_path or in_codon_env):
                return _die(
                    f"{names[0]} not on PATH or in {codon_env_bin} (needed for {why}). "
                    f"Activate the `disassembler` env, or --skip-tree."
                )

    if args.sieve and not resolve_sieve(args.sieve_bin):
        return _die(
            "kmer-diversity-sieve not found. Build it "
            "(`cd extern/kmer-diversity-sieve && make serial`), or pass --sieve-bin."
        )

    recon = out / "reconstructed.fasta"

    # ---- stage 1: walker ---------------------------------------------------
    walker_argv: list[str] = []
    for q in queries:
        walker_argv += ["-q", str(q)]
    walker_argv += ["-o", str(out), "-f", str(recon)]
    if args.workers is not None:
        walker_argv += ["-w", str(args.workers)]
    if args.fetch_workers is not None:
        walker_argv += ["--fetch-workers", str(args.fetch_workers)]
    if args.accession:
        walker_argv += ["-a", args.accession]
    else:
        walker_argv += ["-l", args.accession_list]
    if args.emit_vcf:
        walker_argv += ["--emit-vcf"]
    if args.emit_snv_table:
        walker_argv += ["--emit-snv-table"]
    if args.no_emit_alts:
        walker_argv += ["--no-emit-alts"]
    if args.prune_frac is not None:
        walker_argv += ["--prune-frac", str(args.prune_frac)]
    if args.max_size_mb is not None:
        walker_argv += ["--max-size-mb", str(args.max_size_mb)]
    if args.no_prefer_majority:
        walker_argv += ["--no-prefer-majority"]
    if args.min_abund is not None:
        walker_argv += ["--min-abund", str(args.min_abund)]
    if args.min_vaf is not None:
        walker_argv += ["--min-vaf", str(args.min_vaf)]
    if args.unitigs_cache is not None:
        walker_argv += ["--unitigs-cache", args.unitigs_cache]
    if args.offline:
        walker_argv += ["--offline"]
    # drop a leading `--` separator that argparse.REMAINDER keeps
    extra = args.walker_args[1:] if args.walker_args[:1] == ["--"] else args.walker_args
    walker_argv += list(extra)

    print(f"[reconstruct] 1/2 walker: {' '.join(walker_argv)}", flush=True)
    rc = walker_cmd.run_walker(walker_argv, walker_bin=args.walker_bin)
    if rc != 0:
        return _die(f"logan-walker (stage 1) failed (exit {rc}); see its output above.")

    n_recon = sum(1 for line in recon.open() if line.startswith(">")) if recon.is_file() else 0
    print(f"[reconstruct] reconstructed {n_recon} haplotype(s) -> {recon}", flush=True)
    if n_recon == 0:
        return _die("no haplotypes reconstructed; nothing to align.")

    # ---- optional sieve ----------------------------------------------------
    msa_input = recon
    if args.sieve and not args.skip_tree:
        sieved = out / "reconstructed.sieved.fasta"
        rc = sieve_cmd.sieve_run(
            recon, sieved,
            ref=query,
            min_len=args.sieve_min_len,
            threshold=args.sieve_threshold,
            sketch=args.sieve_sketch,
            sieve_bin=args.sieve_bin,
        )
        if rc != 0:
            return _die(f"sieve failed (exit {rc}); see its output above.")
        msa_input = sieved

    if args.skip_tree:
        print("[reconstruct] --skip-tree: done after reconstruction.")
        return 0

    # ---- stage 2: codon MSA + FastTree ------------------------------------
    stem = msa_ref.stem
    aln = out / f"{stem}_codon_msa.fasta"
    tree = out / f"{stem}_fasttree.nwk"
    print(f"[reconstruct] 2/2 codon MSA + tree: --input {msa_input} --out-aln {aln} "
          f"--out-tree {tree}", flush=True)
    msa_args = types.SimpleNamespace(
        ref=str(msa_ref),
        input=str(msa_input),
        out_aln=str(aln),
        out_tree=str(tree),
        mafft_threads=args.mafft_threads,
        msa_min_cov=args.msa_min_cov,
    )
    rc = msa_cmd.run(msa_args)
    if rc != 0:
        return _die(f"codon MSA / tree (stage 2) failed (exit {rc}); see its output above.")

    print(f"[reconstruct] done:\n"
          f"  haplotypes : {recon}\n"
          f"  codon MSA  : {aln}\n"
          f"  phylogeny  : {tree}")
    return 0
