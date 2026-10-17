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
from disassembler.validators import (
    fasta_has_record,
    fraction,
    positive_float,
    positive_int,
    validate_unitigs_cache,
)


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "reconstruct",
        help="Accessions -> haplotypes -> codon MSA + tree (walker + optional sieve + msa).",
        description=(
            "End-to-end reconstruct pipeline. Chains three atomic subcommands:\n"
            "  1. walker  : reconstruct coding haplotypes from Logan unitigs\n"
            "  2. sieve   : optional diversity thinning (opt-in with --sieve)\n"
            "  3. msa     : codon MSA + FastTree GTR tree\n\n"
            "Options are grouped by pipeline stage below. The Walker group surfaces "
            "only the flags a routine user tunes per run; walker internals and a few "
            "flags with special handling live behind the `--` passthrough (see the "
            "Walker and Pass-through groups below for the exact list and rationale). "
            "`disassembler walker -h` shows the full walker CLI."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Pass extra logan-walker flags after a literal `--`.",
    )

    # ---- Input source (mutually exclusive) --------------------------------
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("-a", "--accession",
                     help="Reconstruct from a single SRA accession.")
    src.add_argument("-l", "--accession-list",
                     help="Reconstruct from every accession listed (one per line) in "
                          "this file.")

    # ---- Stage 1: walker (haplotype reconstruction) -----------------------
    g_walker = p.add_argument_group(
        "Walker: haplotype reconstruction",
        description=(
            "  Options for the logan-walker stage. This group surfaces only the\n"
            "  flags a routine reconstruct user tunes per run. Everything else\n"
            "  stays behind the `--` passthrough (see the last group). Not\n"
            "  surfaced here:\n"
            "\n"
            "    Walker internals (defaults are the right answer for almost every\n"
            "    run; re-declaring them here would silently pin old values when the\n"
            "    walker's default changes):\n"
            "      -k, --hops, --min-shared-kmers, --min-abund-frac,\n"
            "      --alt-max-detour, --max-alt-indel, --no-indel-alts,\n"
            "      --primary-metric\n"
            "\n"
            "    Flags with special handling in this pipeline:\n"
            "      -f/--out-fasta         hardcoded so sieve/MSA can find the FASTA\n"
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
                               "passing more than one record (across one or more -q flags) "
                               "activates the walker's multi-anchor mode. In multi-anchor "
                               "mode the MSA step still needs one reference, defaulting to "
                               "the first -q unless --msa-reference overrides it.")
    g_walker.add_argument("--msa-reference", default=None,
                          help="FASTA to use as the reference for the MSA/tree stage in "
                               "multi-anchor mode (default: the first -q). Ignored under "
                               "--skip-tree.")
    g_walker.add_argument("-o", "--out-dir", required=True,
                          help="Output directory. Walker writes reconstructed.fasta here; "
                               "the MSA and tree land alongside.")
    # None on purpose: numeric defaults would be a second copy of the walker's own,
    # silently pinning old values when the walker's default changed. With None the
    # flag reaches the walker only when the user asked for it.
    g_walker.add_argument("-w", "--workers", type=positive_int, default=None,
                          help="Number of walker compute workers (walker's own default: 8). "
                               "Distinct from --mafft-threads, which is the MSA stage.")
    g_walker.add_argument("--fetch-workers", type=positive_int, default=None,
                          help="Parallel S3 fetch workers for the walker stage "
                               "(walker's own default: 8).")
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
                               "(walker's own default: 0.01). Passed straight to the walker.")
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
                               "(auto-discovered from PATH and the standard build "
                               "locations otherwise).")
    g_walker.add_argument("--unitigs-cache", default=None, metavar="DIR",
                          help="Read-through cache directory for downloaded unitigs. "
                               "Populate once, then rerun offline with --offline.")
    g_walker.add_argument("--offline", action="store_true",
                          help="Refuse the S3 fallback for the walker stage; every accession "
                               "must already be in --unitigs-cache. Requires --unitigs-cache.")

    # ---- Stage 2: sieve (optional diversity thinning) ---------------------
    g_sieve = p.add_argument_group(
        "Sieve: optional diversity thinning",
        description="Optional pass between walker and MSA, delegated to the "
                    "`disassembler sieve` atomic. This is a haplotype/lineage-level "
                    "reducer (MinHash/LSH) meant for LARGE cross-accession cohorts of "
                    "near-duplicate lineages; it is NOT how you preserve intra-host "
                    "SNVs (use --emit-vcf / --emit-snv-table for those). All --sieve-* "
                    "flags are ignored unless --sieve is passed, and the whole sieve "
                    "stage is skipped under --skip-tree.",
    )
    g_sieve.add_argument("--sieve", action="store_true",
                         help="Enable the diversity-thinning pass. Off by default.")
    g_sieve.add_argument("--sieve-threshold", type=fraction, default=0.005,
                         help="LSH promotion threshold tau (default 0.005, gene-scale). "
                              "A sequence is promoted only if 1 - maxJaccard >= tau against "
                              "the already-promoted set, so higher tau collapses more.")
    g_sieve.add_argument("--sieve-sketch", type=positive_int, default=256,
                         help="MinHash sketch size m (default 256, gene-scale). Caps the "
                              "resolution of the novelty estimator at 1/m.")
    g_sieve.add_argument("--sieve-min-len", type=positive_int, default=None,
                         help="Minimum sequence length gate for the sieve "
                              "(default: auto = 95%% of the first -q record's length).")
    g_sieve.add_argument("--sieve-bin", default=None,
                         help="Explicit path to the kmer-diversity-sieve l1_pipeline_serial "
                              "binary (auto-discovered otherwise).")

    # ---- Stage 3: MSA + tree ----------------------------------------------
    g_msa = p.add_argument_group(
        "MSA + tree",
        description="Options for the codon-MSA + FastTree GTR stage, delegated to the "
                    "`disassembler msa` atomic. Skipped entirely under --skip-tree.",
    )
    g_msa.add_argument("--mafft-threads", type=positive_int, default=8,
                       help="Thread count for MAFFT in the MSA stage (default: 8). "
                            "Distinct from the walker's -t/--threads, which sets "
                            "GraphAligner's thread count and is kept at 1 for "
                            "reproducibility.")
    g_msa.add_argument("--msa-min-cov", type=fraction, default=0.70,
                       help="Minimum coverage fraction for a reconstructed sequence to "
                            "enter the MSA (default: 0.70). Distinct from the walker's own "
                            "--min-cov (default 0.85), which gates which reconstructions "
                            "are produced at all.")
    g_msa.add_argument("--skip-tree", action="store_true",
                       help="Stop after the walker (and sieve, if enabled). Do not build "
                            "the MSA or tree. --sieve is also ignored under this flag.")

    # ---- Pass-through -----------------------------------------------------
    g_pass = p.add_argument_group(
        "Pass-through to logan-walker",
        description=(
            "  Everything after a literal `--` on the command line is appended\n"
            "  verbatim to the walker invocation. Reach for this to set any of the\n"
            "  walker flags NOT surfaced in the Walker group above, e.g.:\n"
            "\n"
            "    disassembler reconstruct ... -- --primary-metric pident --hops 7\n"
            "\n"
            "  See `disassembler walker -h` for the full CLI.\n"
        ),
    )
    g_pass.add_argument("walker_args", nargs=argparse.REMAINDER,
                        help="Extra args forwarded to logan-walker (put them after `--`).")
    p.set_defaults(func=run)


def _die(msg: str) -> int:
    print(f"reconstruct: {msg}", file=sys.stderr)
    return 1


# Walker flags whose value the composite computes itself. If the user re-sets
# any of these via the `--` passthrough the composite refuses, because the
# override would either silently break the stage-1 -> stage-2 handoff
# (`-f/--out-fasta` is where the sieve/MSA stage reads from) or duplicate an
# option the composite has already surfaced under its own name.
COMPOSITE_OWNED_WALKER_FLAGS = {
    "-f": "output FASTA is <out-dir>/reconstructed.fasta; change --out-dir instead",
    "--out-fasta": "output FASTA is <out-dir>/reconstructed.fasta; change --out-dir instead",
    "-o": "use the composite's own -o/--out-dir",
    "--out-dir": "use the composite's own -o/--out-dir",
    "-a": "use the composite's own -a/--accession",
    "--accession": "use the composite's own -a/--accession",
    "-l": "use the composite's own -l/--accession-list",
    "--accession-list": "use the composite's own -l/--accession-list",
    "-q": "use the composite's own -q/--query (repeatable for multi-anchor)",
    "--query": "use the composite's own -q/--query (repeatable for multi-anchor)",
}


def _scrub_walker_passthrough(extra: list[str]) -> str | None:
    """Return an error message if `extra` re-sets a composite-owned walker flag."""
    for tok in extra:
        # argparse.REMAINDER hands us the raw token stream. Match on the flag
        # itself, not its value; both `--out-fasta path` and `--out-fasta=path`
        # start with a token that begins with the flag name.
        head = tok.split("=", 1)[0]
        if head in COMPOSITE_OWNED_WALKER_FLAGS:
            return (
                f"{head} is set by the composite and cannot be overridden through "
                f"`--`. {COMPOSITE_OWNED_WALKER_FLAGS[head]}."
            )
    return None


def _warn_inert_knobs(args, prefix: str = "[reconstruct]") -> None:
    """Warn on stage-specific tuning flags whose stage will not run.

    Sieve knobs when --sieve is off; MSA/tree knobs when --skip-tree is on.
    Compared against the same defaults the parser declares; keep in sync with
    add_parser above.
    """
    if not args.sieve:
        sieve_inert = []
        if args.sieve_threshold != 0.005:
            sieve_inert.append("--sieve-threshold")
        if args.sieve_sketch != 256:
            sieve_inert.append("--sieve-sketch")
        if args.sieve_min_len is not None:
            sieve_inert.append("--sieve-min-len")
        if args.sieve_bin is not None:
            sieve_inert.append("--sieve-bin")
        if sieve_inert:
            print(f"{prefix} warn: {', '.join(sieve_inert)} set but --sieve is off, "
                  f"so these values are inert. Add --sieve to make them take effect.",
                  file=sys.stderr)

    if args.skip_tree:
        msa_inert = []
        if args.mafft_threads != 8:
            msa_inert.append("--mafft-threads")
        if args.msa_min_cov != 0.70:
            msa_inert.append("--msa-min-cov")
        if args.msa_reference is not None:
            msa_inert.append("--msa-reference")
        if msa_inert:
            print(f"{prefix} warn: {', '.join(msa_inert)} set but --skip-tree is on, "
                  f"so the MSA/tree stage will not run and these values are inert.",
                  file=sys.stderr)


def run(args) -> int:
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    queries = [Path(q) for q in args.query]
    for q in queries:
        if not q.is_file():
            return _die(f"query not found: {q}")
        if fasta_has_record(str(q)) == 0:
            return _die(
                f"query FASTA has no records (no `>` header lines): {q}. Point "
                f"-q at a CDS FASTA, not an empty or non-FASTA file."
            )
    # MSA/tree stage still takes a single reference. Default is the first -q, override
    # via --msa-reference; irrelevant under --skip-tree.
    msa_ref = Path(args.msa_reference) if args.msa_reference else queries[0]
    if not args.skip_tree and not msa_ref.is_file():
        return _die(f"MSA reference not found: {msa_ref}")
    if len(queries) > 1 and not args.skip_tree and args.msa_reference is None:
        print(f"[reconstruct] multi-anchor: MSA reference defaults to first -q ({msa_ref}); "
              f"pass --msa-reference to override.", flush=True)

    # Unitigs cache sanity: must be a directory if it exists; must exist under
    # --offline. The walker enforces most of this too, but checking here fails
    # fast before any binary is spawned or (in from-kmindex) before harvest runs.
    cache_err = validate_unitigs_cache(args.unitigs_cache, offline=args.offline)
    if cache_err:
        return _die(cache_err)

    # Scrub composite-owned flags out of the walker passthrough.
    raw_extra = args.walker_args[1:] if args.walker_args[:1] == ["--"] else args.walker_args
    scrub_err = _scrub_walker_passthrough(list(raw_extra))
    if scrub_err:
        return _die(scrub_err)

    if not resolve_walker(args.walker_bin):
        return _die(
            f"logan-walker not found. Build it (`cd walker && cargo build --release`), "
            f"put it on PATH, or pass --walker-bin."
        )

    _warn_inert_knobs(args)

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
    walker_argv += list(raw_extra)

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
        # The default --sieve-min-len is sized on the reference the MSA stage measures
        # the haplotypes against; with a single -q that is the query itself.
        rc = sieve_cmd.sieve_run(
            recon, sieved,
            ref=msa_ref,
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
