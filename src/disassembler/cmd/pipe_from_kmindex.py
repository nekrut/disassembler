"""from-kmindex: Galaxy kmindex harvest -> stratify -> reconstruct.

Composite that chains `harvest`, `stratify` (optional), and `reconstruct` for
one target at a time. Multi-target work should compose the atomics; this
convenience wraps the single-target end-to-end path.

Accession-list selection after stratify:
  --diversity N -> <cohort-dir>/<target>.diversity_<N>.txt
  --bin ...     -> <cohort-dir>/<target>.bins.txt
  (neither)     -> <cohort-dir>/<target>.accessions.txt   (the raw harvest list)
"""

from __future__ import annotations

import argparse
import sys
import types
from pathlib import Path

from disassembler.cmd import harvest as harvest_cmd
from disassembler.cmd import pipe_reconstruct
from disassembler.cmd import stratify as stratify_cmd


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "from-kmindex",
        help="Galaxy kmindex -> stratified cohort -> reconstruct (walker + sieve + msa).",
        description=(
            "Composite: harvest a Galaxy kmindex shard collection, optionally cut a "
            "stratified cohort for one target, then run the reconstruct pipeline on it. "
            "Anything after `--` is forwarded to logan-walker."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Pass extra logan-walker flags after a literal `--`.",
    )

    # harvest inputs
    p.add_argument("--history-id", required=True,
                   help="Galaxy history ID containing the kmindex collection")
    p.add_argument("--collection-id", required=True,
                   help="Galaxy dataset collection ID for kmindex query shards")
    p.add_argument("--base-url", default="https://usegalaxy.org",
                   help="Galaxy instance URL (default: https://usegalaxy.org)")
    p.add_argument("--api-key-file", default=None,
                   help="Galaxy API key file (also honours $GALAXY_API_KEY, then /tmp/gxy.txt)")
    p.add_argument("--harvest-workers", type=int, default=6,
                   help="parallel harvest download workers (default: 6, safe for usegalaxy.org's "
                        "rate limit; raise for a self-hosted Galaxy)")
    p.add_argument("--cohort-dir", default=None,
                   help="Directory for harvested + stratified accession lists "
                        "(default: <out-dir>/cohort)")

    # stratify (optional)
    p.add_argument("--bin", action="append", dest="bins", metavar="LO:HI:N",
                   help="Custom score bins for stratify; repeatable. Nonzero scores only, "
                        "no backfill, disjoint. If omitted, reconstruct runs against the raw "
                        "harvest list <target>.accessions.txt.")

    # target selection into reconstruct
    p.add_argument("--target", required=True,
                   help="Target name (as produced by `harvest`, e.g. 'ns5' or 'HA') to "
                        "reconstruct from.")

    # reconstruct inputs
    p.add_argument("-q", "--query", required=True, action="append",
                   help="Reference in-frame CDS FASTA. Repeatable, and each file may hold "
                        "multiple records; passing more than one activates the walker's "
                        "multi-anchor mode (union anchor pool, primary is the reference "
                        "whose alignment gives the best cov/pident).")
    p.add_argument("--msa-reference", default=None,
                   help="FASTA to use as the MSA reference in multi-anchor mode "
                        "(default: first -q). Ignored under --skip-tree.")
    p.add_argument("-o", "--out-dir", required=True, help="Output directory for reconstruct")
    p.add_argument("-w", "--workers", type=int, default=None,
                   help="walker compute workers (walker default 8)")
    p.add_argument("--fetch-workers", type=int, default=None,
                   help="walker S3 fetch workers (walker default 8)")
    p.add_argument("--mafft-threads", type=int, default=8, help="MAFFT threads (default 8)")
    p.add_argument("--msa-min-cov", type=float, default=0.70,
                   help="Min coverage for a sequence to enter the MSA (default 0.70). "
                        "Distinct from walker's --min-cov (default 0.85)")
    p.add_argument("--emit-vcf", action="store_true")
    p.add_argument("--emit-snv-table", action="store_true")
    p.add_argument("--no-emit-alts", action="store_true")
    p.add_argument("--prune-frac", type=float, default=None)
    p.add_argument("--max-size-mb", type=int, default=None)
    p.add_argument("--no-prefer-majority", action="store_true")
    p.add_argument("--min-abund", type=float, default=None)
    p.add_argument("--min-vaf", type=float, default=None)
    p.add_argument("--walker-bin", default=None)
    p.add_argument("--skip-tree", action="store_true")
    p.add_argument("--sieve", action="store_true")
    p.add_argument("--sieve-threshold", type=float, default=0.005)
    p.add_argument("--sieve-sketch", type=int, default=256)
    p.add_argument("--sieve-min-len", type=int, default=None)
    p.add_argument("--sieve-bin", default=None)
    p.add_argument("--unitigs-cache", default=None, metavar="DIR",
                   help="Read-through unitigs cache passed to walker.")
    p.add_argument("--offline", action="store_true",
                   help="Refuse the S3 fallback for the walker stage. Requires --unitigs-cache.")
    p.add_argument("--skip-harvest", action="store_true",
                   help="Reuse existing <cohort-dir>/*.scores.json instead of hitting Galaxy. "
                        "Combine with --offline and a pre-warmed --unitigs-cache for a "
                        "fully offline kmindex -> tree run.")
    p.add_argument("walker_args", nargs=argparse.REMAINDER,
                   help="extra args passed to logan-walker after `--`")
    p.set_defaults(func=run)


def run(args) -> int:
    out_dir = Path(args.out_dir)
    cohort_dir = Path(args.cohort_dir) if args.cohort_dir else out_dir / "cohort"
    cohort_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1/3 harvest -------------------------------------------------------
    if args.skip_harvest:
        # Reuse existing scores; no Galaxy call.
        existing = sorted(cohort_dir.glob("*.scores.json"))
        if not existing:
            print(f"from-kmindex: --skip-harvest but no *.scores.json in {cohort_dir}",
                  file=sys.stderr)
            return 1
        print(f"[from-kmindex] 1/3 harvest skipped: reusing {len(existing)} "
              f"*.scores.json in {cohort_dir}", flush=True)
    else:
        harvest_ns = types.SimpleNamespace(
            history_id=args.history_id,
            collection_id=args.collection_id,
            base_url=args.base_url,
            api_key_file=args.api_key_file,
            out_dir=str(cohort_dir),
            download_workers=args.harvest_workers,
        )
        print(f"[from-kmindex] 1/3 harvest -> {cohort_dir}", flush=True)
        rc = harvest_cmd.run(harvest_ns)
        if rc != 0:
            return rc

    # ---- 2/3 stratify (optional) ------------------------------------------
    if args.bins:
        stratify_ns = types.SimpleNamespace(
            out_dir=str(cohort_dir),
            bins=args.bins,
        )
        print(f"[from-kmindex] 2/3 stratify -> {cohort_dir}", flush=True)
        rc = stratify_cmd.run(stratify_ns)
        if rc != 0:
            return rc
    else:
        print("[from-kmindex] 2/3 stratify skipped (no --bin)", flush=True)

    # ---- select accession list for --target -------------------------------
    if args.bins:
        acc_list = cohort_dir / f"{args.target}.bins.txt"
    else:
        acc_list = cohort_dir / f"{args.target}.accessions.txt"
    if not acc_list.is_file():
        print(f"from-kmindex: expected accession list not found: {acc_list}",
              file=sys.stderr)
        return 1
    n_acc = sum(1 for line in acc_list.open() if line.strip())
    print(f"[from-kmindex] target '{args.target}': {n_acc} accessions from {acc_list}",
          flush=True)

    # ---- 3/3 reconstruct --------------------------------------------------
    print(f"[from-kmindex] 3/3 reconstruct -> {out_dir}", flush=True)
    reconstruct_ns = types.SimpleNamespace(
        accession=None,
        accession_list=str(acc_list),
        query=args.query,
        msa_reference=args.msa_reference,
        out_dir=str(out_dir),
        workers=args.workers,
        fetch_workers=args.fetch_workers,
        mafft_threads=args.mafft_threads,
        msa_min_cov=args.msa_min_cov,
        emit_vcf=args.emit_vcf,
        emit_snv_table=args.emit_snv_table,
        no_emit_alts=args.no_emit_alts,
        prune_frac=args.prune_frac,
        max_size_mb=args.max_size_mb,
        no_prefer_majority=args.no_prefer_majority,
        min_abund=args.min_abund,
        min_vaf=args.min_vaf,
        walker_bin=args.walker_bin,
        skip_tree=args.skip_tree,
        sieve=args.sieve,
        sieve_threshold=args.sieve_threshold,
        sieve_sketch=args.sieve_sketch,
        sieve_min_len=args.sieve_min_len,
        sieve_bin=args.sieve_bin,
        unitigs_cache=args.unitigs_cache,
        offline=args.offline,
        walker_args=args.walker_args,
    )
    return pipe_reconstruct.run(reconstruct_ns)
