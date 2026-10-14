"""from-kmindex: Galaxy kmindex harvest -> stratify -> reconstruct.

Composite that chains `harvest`, `stratify` (optional), and `reconstruct` for
one target at a time. Multi-target work should compose the atomics; this
convenience wraps the single-target end-to-end path.

Accession-list selection after stratify:
  --score-bin ... -> <cohort-dir>/<target>.bins.txt
  (omitted)       -> <cohort-dir>/<target>.accessions.txt   (the raw harvest list)
"""

from __future__ import annotations

import argparse
import sys
import types
from pathlib import Path

from disassembler.cmd import harvest as harvest_cmd
from disassembler.cmd import pipe_reconstruct
from disassembler.cmd import stratify as stratify_cmd
from disassembler.validators import (
    fraction,
    positive_float,
    positive_int,
    validate_api_key_file,
    validate_unitigs_cache,
)


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "from-kmindex",
        help="Galaxy kmindex -> stratified cohort -> reconstruct (walker + sieve + msa).",
        description=(
            "End-to-end composite from a Galaxy kmindex search. Chains five stages:\n"
            "  1. harvest    : pull the kmindex shard collection from Galaxy\n"
            "  2. stratify   : optional score-bin cohort cut (opt-in with --score-bin)\n"
            "  3. walker     : reconstruct coding haplotypes from Logan unitigs\n"
            "  4. sieve      : optional diversity thinning (opt-in with --sieve)\n"
            "  5. msa        : codon MSA + FastTree GTR tree\n\n"
            "Options are grouped by pipeline stage below. The Walker group surfaces "
            "only the flags a routine user tunes per run; walker internals and a few "
            "flags with special handling live behind the `--` passthrough (see the "
            "Walker and Pass-through groups below for the exact list and rationale). "
            "`disassembler walker -h` shows the full walker CLI."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Pass extra logan-walker flags after a literal `--`.",
    )

    # ---- Stage 1: harvest (Galaxy kmindex) --------------------------------
    g_harvest = p.add_argument_group(
        "Harvest: Galaxy kmindex",
        description="Pull a kmindex shard collection from a Galaxy history. See "
                    "`disassembler harvest -h` for the underlying atomic.",
    )
    g_harvest.add_argument("--history-id", required=True,
                           help="Galaxy history ID that contains the kmindex dataset "
                                "collection to harvest.")
    g_harvest.add_argument("--collection-id", required=True,
                           help="Galaxy dataset collection ID for the kmindex query shards.")
    g_harvest.add_argument("--base-url", default="https://usegalaxy.org",
                           help="Galaxy instance URL "
                                "(default: https://usegalaxy.org).")
    g_harvest.add_argument("--api-key-file", default=None,
                           help="Path to a Galaxy API key file. If unset, honours "
                                "$GALAXY_API_KEY and then /tmp/gxy.txt.")
    g_harvest.add_argument("--harvest-workers", type=positive_int, default=6,
                           help="Parallel download workers for the harvest stage "
                                "(default: 6, safe for usegalaxy.org's rate limit; "
                                "raise for a self-hosted Galaxy).")
    g_harvest.add_argument("--cohort-dir", default=None,
                           help="Directory for harvested *.scores.json files and any "
                                "stratified accession lists "
                                "(default: <out-dir>/cohort).")
    g_harvest.add_argument("--skip-harvest", action="store_true",
                           help="Reuse existing <cohort-dir>/*.scores.json instead of "
                                "calling Galaxy. Combine with --offline and a pre-warmed "
                                "--unitigs-cache for a fully offline kmindex -> tree run.")

    # ---- Stage 2: stratify (optional) -------------------------------------
    g_stratify = p.add_argument_group(
        "Stratify (optional) + target selection",
        description="Optional score-bin cohort cut on the harvested scores, plus the "
                    "target name whose accession list is fed downstream. Without any "
                    "--score-bin, reconstruct runs against the raw harvest list "
                    "<target>.accessions.txt. See `disassembler stratify -h`.",
    )
    g_stratify.add_argument("--score-bin", action="append", dest="score_bins",
                            metavar="LO:HI:N",
                            help="Custom score bin, repeatable. Sample up to N accessions "
                                 "with score in the half-open [LO, HI). Nonzero scores "
                                 "only, no backfill, disjoint. When any --score-bin is "
                                 "given, reconstruct reads <target>.bins.txt.")
    g_stratify.add_argument("--target", required=True,
                            help="Target name to reconstruct from (as produced by "
                                 "`harvest`, e.g. 'ns5' or 'HA'). Selects which per-target "
                                 "accession list under <cohort-dir> is fed to walker.")

    # ---- Stage 3: walker (haplotype reconstruction) -----------------------
    g_walker = p.add_argument_group(
        "Walker: haplotype reconstruction",
        description=(
            "  Options for the logan-walker stage. This group surfaces only the\n"
            "  flags a routine user tunes per run. Everything else stays behind the\n"
            "  `--` passthrough (see the last group). Not surfaced here:\n"
            "\n"
            "    Walker internals (defaults are the right answer for almost every\n"
            "    run; re-declaring them here would silently pin old values when the\n"
            "    walker's default changes):\n"
            "      -k, --hops, --min-shared-kmers, --min-abund-frac,\n"
            "      --alt-max-detour, --max-alt-indel, --no-indel-alts,\n"
            "      --primary-metric\n"
            "\n"
            "    Flags with special handling in this pipeline:\n"
            "      -f/--out-fasta         hardcoded so downstream stages find the FASTA\n"
            "      -t/--threads           kept at 1 for GraphAligner reproducibility\n"
            "      --min-cov              would clash with MSA's --msa-min-cov\n"
            "      --allow-frameshift-alts   name owned by lexicmap for a different\n"
            "                                meaning\n"
            "\n"
            "  See `disassembler walker -h` for the full walker CLI.\n"
        ),
    )
    g_walker.add_argument("-q", "--query", required=True, action="append",
                          help="Reference in-frame CDS FASTA used as the walker's anchor. "
                               "Repeatable, and each file may hold multiple records; "
                               "passing more than one record activates the walker's "
                               "multi-anchor mode (union anchor pool, primary is the "
                               "reference whose alignment gives the best cov/pident).")
    g_walker.add_argument("--msa-reference", default=None,
                          help="FASTA to use as the reference for the MSA/tree stage in "
                               "multi-anchor mode (default: the first -q). Ignored under "
                               "--skip-tree.")
    g_walker.add_argument("-o", "--out-dir", required=True,
                          help="Output directory for reconstruct. Walker writes "
                               "reconstructed.fasta here; the MSA and tree land alongside.")
    g_walker.add_argument("-w", "--workers", type=positive_int, default=None,
                          help="Number of walker compute workers "
                               "(walker's own default: 8).")
    g_walker.add_argument("--fetch-workers", type=positive_int, default=None,
                          help="Parallel S3 fetch workers for the walker stage "
                               "(walker's own default: 8). Distinct from --harvest-workers, "
                               "which parallelises the Galaxy harvest above.")
    g_walker.add_argument("--emit-vcf", action="store_true",
                          help="Also write per-accession intra-host SNP VCFs from the "
                               "walker stage.")
    g_walker.add_argument("--emit-snv-table", action="store_true",
                          help="Also write the long snvs.tsv cohort table from the walker "
                               "stage.")
    g_walker.add_argument("--no-emit-alts", action="store_true",
                          help="Do not write alt_haplotypes.fasta from the walker stage.")
    g_walker.add_argument("--prune-frac", type=fraction, default=None,
                          help="Unitig depth prune fraction for reference bias removal "
                               "(walker's own default: 0.01).")
    g_walker.add_argument("--max-size-mb", type=positive_int, default=None,
                          help="Skip any accession whose compressed SRA size exceeds this "
                               "many MB (walker's own default: 1024).")
    g_walker.add_argument("--no-prefer-majority", action="store_true",
                          help="Keep the primary haplotype on the reference allele instead "
                               "of the majority allele at each variant site.")
    g_walker.add_argument("--min-abund", type=positive_float, default=None,
                          help="Minimum ka:f (k-mer abundance) for a variant to be kept "
                               "(walker's own default: 5.0).")
    g_walker.add_argument("--min-vaf", type=fraction, default=None,
                          help="Minimum variant allele fraction for an alternative branch "
                               "to be emitted (walker's own default: 0.10).")
    g_walker.add_argument("--walker-bin", default=None,
                          help="Explicit path to the logan-walker binary "
                               "(auto-discovered otherwise).")
    g_walker.add_argument("--unitigs-cache", default=None, metavar="DIR",
                          help="Read-through cache directory for downloaded unitigs. "
                               "Populate once, then rerun offline with --offline.")
    g_walker.add_argument("--offline", action="store_true",
                          help="Refuse the S3 fallback for the walker stage; every "
                               "accession must already be in --unitigs-cache. Requires "
                               "--unitigs-cache.")

    # ---- Stage 4: sieve (optional diversity thinning) ---------------------
    g_sieve = p.add_argument_group(
        "Sieve: optional diversity thinning",
        description="Optional pass between walker and MSA. Haplotype/lineage-level "
                    "MinHash/LSH reducer meant for LARGE cohorts of near-duplicate "
                    "lineages; NOT how you preserve intra-host SNVs (use --emit-vcf / "
                    "--emit-snv-table for those). All --sieve-* flags are ignored "
                    "unless --sieve is passed, and the whole stage is skipped under "
                    "--skip-tree.",
    )
    g_sieve.add_argument("--sieve", action="store_true",
                         help="Enable the diversity-thinning pass. Off by default.")
    g_sieve.add_argument("--sieve-threshold", type=fraction, default=0.005,
                         help="LSH promotion threshold tau (default 0.005, gene-scale). "
                              "Higher tau collapses more.")
    g_sieve.add_argument("--sieve-sketch", type=positive_int, default=256,
                         help="MinHash sketch size m (default 256, gene-scale). Caps the "
                              "resolution of the novelty estimator at 1/m.")
    g_sieve.add_argument("--sieve-min-len", type=positive_int, default=None,
                         help="Minimum sequence length gate for the sieve "
                              "(default: auto = 95%% of the first -q record's length).")
    g_sieve.add_argument("--sieve-bin", default=None,
                         help="Explicit path to the kmer-diversity-sieve l1_pipeline_serial "
                              "binary (auto-discovered otherwise).")

    # ---- Stage 5: MSA + tree ----------------------------------------------
    g_msa = p.add_argument_group(
        "MSA + tree",
        description="Options for the codon-MSA + FastTree GTR stage. Skipped entirely "
                    "under --skip-tree.",
    )
    g_msa.add_argument("--mafft-threads", type=positive_int, default=8,
                       help="Thread count for MAFFT in the MSA stage (default: 8).")
    g_msa.add_argument("--msa-min-cov", type=fraction, default=0.70,
                       help="Minimum coverage fraction for a reconstructed sequence to "
                            "enter the MSA (default: 0.70). Distinct from the walker's "
                            "own --min-cov (default 0.85), which gates which "
                            "reconstructions are produced at all.")
    g_msa.add_argument("--skip-tree", action="store_true",
                       help="Stop after the walker (and sieve, if enabled). Do not build "
                            "the MSA or tree.")

    # ---- Pass-through -----------------------------------------------------
    g_pass = p.add_argument_group(
        "Pass-through to logan-walker",
        description=(
            "  Everything after a literal `--` on the command line is appended\n"
            "  verbatim to the walker invocation. Reach for this to set any of the\n"
            "  walker flags NOT surfaced in the Walker group above, e.g.:\n"
            "\n"
            "    disassembler from-kmindex ... -- --primary-metric pident --hops 7\n"
            "\n"
            "  See `disassembler walker -h` for the full CLI.\n"
        ),
    )
    g_pass.add_argument("walker_args", nargs=argparse.REMAINDER,
                        help="Extra args forwarded to logan-walker (put them after `--`).")
    p.set_defaults(func=run)


def run(args) -> int:
    # Fail fast on argument bugs BEFORE spending a Galaxy round-trip on harvest.
    # These same checks fire again inside pipe_reconstruct.run at stage 3, but
    # catching them here saves the harvest.
    cache_err = validate_unitigs_cache(args.unitigs_cache, offline=args.offline)
    if cache_err:
        print(f"from-kmindex: {cache_err}", file=sys.stderr)
        return 1
    raw_extra = (
        args.walker_args[1:] if args.walker_args[:1] == ["--"] else args.walker_args
    )
    scrub_err = pipe_reconstruct._scrub_walker_passthrough(list(raw_extra))
    if scrub_err:
        print(f"from-kmindex: {scrub_err}", file=sys.stderr)
        return 1
    # API key file: only relevant when we're actually going to hit Galaxy.
    if not args.skip_harvest:
        key_err = validate_api_key_file(args.api_key_file)
        if key_err:
            print(f"from-kmindex: {key_err}", file=sys.stderr)
            return 1

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
    if args.score_bins:
        stratify_ns = types.SimpleNamespace(
            out_dir=str(cohort_dir),
            score_bins=args.score_bins,
        )
        print(f"[from-kmindex] 2/3 stratify -> {cohort_dir}", flush=True)
        rc = stratify_cmd.run(stratify_ns)
        if rc != 0:
            return rc
    else:
        print("[from-kmindex] 2/3 stratify skipped (no --score-bin)", flush=True)

    # ---- select accession list for --target -------------------------------
    if args.score_bins:
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
