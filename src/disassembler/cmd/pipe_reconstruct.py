"""reconstruct: accessions + query CDS -> haplotypes -> (sieve, opt) -> codon MSA + tree.

Composite: walker -> sieve (optional) -> msa. Anything after a literal `--` is
passed straight through to logan-walker.
"""

from __future__ import annotations

import argparse
import time
import sys
from pathlib import Path

from disassembler.cmd import msa as msa_cmd
from disassembler.cmd import sieve as sieve_cmd
from disassembler.cmd import walker as walker_cmd
from disassembler.resolve import resolve_sieve, resolve_walker, not_found_message, EX_TOOL_NOT_FOUND, walker_version_error
from disassembler.stage import stage_args
from disassembler.validators import (

    fasta_has_record,
    fasta_record_lengths,
    fraction,
    positive_float,
    positive_int,
    validate_unitigs_cache,
)
from disassembler import summary as summary_io


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
                               "activates the walker's multi-anchor mode. Those records must "
                               "be strains of ONE gene: the walker keeps one primary per "
                               "accession, so a pool of different genes reports only the "
                               "winner. In multi-anchor mode the MSA step still needs one "
                               "reference, defaulting to the first -q unless --msa-reference "
                               "overrides it.")
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
                               "many MB (walker's own default: 4096).")
    g_walker.add_argument("--no-prefer-majority", action="store_true",
                          help="Keep the primary haplotype on the reference allele instead "
                               "of the majority allele at each variant site.")
    g_walker.add_argument("--min-abund", type=positive_float, default=None,
                          help="Minimum ka:f (k-mer abundance) for a variant to be kept "
                               "(walker's own default: 5.0). This is the allele floor only; "
                               "the prune's lower bound is the walker's --prune-min-abund, "
                               "which you can pass after a literal `--`.")
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
                         help="Minimum sequence length gate for the sieve (default: auto = "
                              "the walker's own --min-cov, 0.85 unless you pass another "
                              "after --, times the shortest -q record's length, so the "
                              "sieve does not drop reconstructions the walker published).")
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
                       help="Stop after the walker. The sieve, the MSA and the tree are "
                            "all skipped.")

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


# logan-walker's exit code for "some accessions reconstructed and some errored". Anything
# else non-zero from the walker means there is nothing worth carrying forward: 1 is nothing
# reconstructed, 3 is nothing worked at all, 127 is no binary.
EX_SOME_ERRORED = 4


def _stage_code(rc: int) -> int:
    """A failed stage's own exit code, carried out of the composite rather than flattened.

    Each stage has a vocabulary worth keeping: the walker's 3 is "no accession worked at all",
    distinct from its 1, "nothing was reconstructed"; `msa` returns 2 for a usage error, the
    same meaning argparse gives it; `resolve` returns 127 for a missing binary. Reporting all
    of them as 1 told a caller only that something went wrong.

    A code that is not positive is a stage that did not normalise a signal death (`sieve` does,
    via `report_rc`). Returning it as is would let the shell mask it to an unrelated byte, so
    it becomes 1, which is at least true.
    """
    return rc if rc > 0 else 1


def _die(msg: str, code: int = 1) -> int:
    # `code` so a missing binary can exit 127 here as it does in the stage. The
    # message already names the command, so it is not prefixed twice.
    print(msg if msg.startswith("reconstruct: ") else f"reconstruct: {msg}", file=sys.stderr)
    return code


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

# Not owned by the composite, but they make the walker exit 0 WITHOUT writing its
# `-f`, which used to leave the composite aligning whatever was already at that path. Each
# has a better form of its own, so they are refused by name rather than tolerated.
WALKER_FLAGS_WITH_NO_OUTPUT = {
    "--fetch-only": "it makes the walker exit without writing a reconstruction, and this "
                    "composite needs one. To warm a cache, run `disassembler fetch` or "
                    "`disassembler walker --fetch-only`",
    "-h": "it prints the walker's help instead of running it; "
          "run `disassembler walker -h`",
    "--help": "it prints the walker's help instead of running it; "
              "run `disassembler walker --help`",
    "-V": "it prints the walker's version instead of running it; "
          "run `disassembler walker --version`",
    "--version": "it prints the walker's version instead of running it; "
                 "run `disassembler walker --version`",
}


# The walker's own default `--min-cov` (walker/src/config.rs). The sieve's completeness
# gate is derived from it so the two stages cannot disagree about which reconstructions
# count as complete.
WALKER_DEFAULT_MIN_COV = 0.85


def _walker_min_cov(extra: list[str]) -> float:
    """The `--min-cov` the walker will actually use, given the passthrough.

    The last occurrence wins, as it does in the walker's own argument parser, and the
    walker accepts only the two-token form, never `--min-cov=X`. An unparseable or
    out-of-range value is left to the walker to reject; this returns the default so the
    sieve gate does not become the thing that fails first.
    """
    cov = WALKER_DEFAULT_MIN_COV
    for i, tok in enumerate(extra):
        if tok == "--min-cov" and i + 1 < len(extra):
            try:
                v = float(extra[i + 1])
            except (TypeError, ValueError):
                continue
            if 0.0 <= v <= 1.0:
                cov = v
    return cov


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
        if head in WALKER_FLAGS_WITH_NO_OUTPUT:
            return (
                f"{head} cannot be passed through `--`: "
                f"{WALKER_FLAGS_WITH_NO_OUTPUT[head]}."
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
    """Run the stages, then write the manifest that ties their summaries together."""
    stages: list[dict] = []
    started = time.time()
    try:
        rc = _run(args, stages)
    except BaseException as e:
        _write_manifest(args, stages, started, "failed", f"{type(e).__name__}: {e}")
        raise
    _write_manifest(args, stages, started, None if rc == 0 else "failed")
    return rc


def _write_manifest(args, stages, started, status=None, error=None) -> None:
    try:
        summary_io.manifest(
            Path(args.out_dir) / "reconstruct.summary.json",
            name="reconstruct",
            settings={k: v for k, v in sorted(vars(args).items())
                      if isinstance(v, (str, int, float, bool, type(None)))},
            source=summary_io.source("file", str(args.out_dir), len(stages)),
            stages=stages,
            started=started,
            status=status,
            error=error,
        )
    except OSError:
        # The output directory may be the thing that failed; the stage already said so.
        pass


def _run(args, stages) -> int:
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    queries = [Path(q) for q in args.query]
    for q in queries:
        if not q.is_file():
            return _die(f"query not found: {q}")
        try:
            n_records = fasta_has_record(str(q))
        except (OSError, ValueError) as e:
            return _die(f"cannot read -q: {e}")
        if n_records == 0:
            return _die(
                f"query FASTA has no records (no `>` header lines): {q}. Point "
                f"-q at a CDS FASTA, not an empty or non-FASTA file."
            )
    # MSA/tree stage still takes a single reference. Default is the first -q, override
    # via --msa-reference; irrelevant under --skip-tree.
    msa_ref = Path(args.msa_reference) if args.msa_reference else queries[0]
    if not args.skip_tree and not msa_ref.is_file():
        return _die(f"MSA reference not found: {msa_ref}")
    # Multi-anchor mode is any reference set of more than one record, from several -q or from
    # one multi-record file; the notice used to fire only for several -q.
    references = [r for q in queries for r in fasta_record_lengths(q)]
    if len(references) > 1 and not args.skip_tree and args.msa_reference is None:
        print(f"[reconstruct] multi-anchor: {len(references)} references. The MSA aligns "
              f"against the first record of the first -q ({references[0][0]} in {msa_ref}); "
              f"pass --msa-reference to choose another.", flush=True)

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

    _walker = resolve_walker(args.walker_bin)
    if not _walker:
        # The same message the `walker` stage prints, and the same exit code.
        return _die(not_found_message("walker", "reconstruct"), EX_TOOL_NOT_FOUND)
    # Refuse an old walker in the preflight, before the fetch and the walk.
    _ver_err = walker_version_error(_walker)
    if _ver_err:
        return _die(_ver_err)

    _warn_inert_knobs(args)

    if args.sieve and args.skip_tree:
        print("[reconstruct] warn: --sieve is ignored under --skip-tree. The sieve thins "
              "the input to the MSA/tree stage, and that stage is being skipped, so "
              "reconstructed.fasta is written unthinned. Drop --skip-tree to get "
              "reconstructed.sieved.fasta.", file=sys.stderr)

    if not args.skip_tree:
        # Preflight MAFFT and FastTree before the expensive walk stage, the same lookup
        # `msa.find_bin` does.
        for names, why in ((("mafft",), "the codon MSA"),
                           (("fasttree", "FastTree"), "the ML tree")):
            try:
                msa_cmd.find_bin(*names)
            except FileNotFoundError as e:
                return _die(f"{e} It is needed for {why}; or pass --skip-tree.")

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
    # The reconstruction is found by a fixed name and nothing checked that THIS run
    # wrote it, so in an --out-dir holding an earlier run the composite reported
    # "reconstructed N haplotype(s)" and built an MSA and a tree from the old cohort, with
    # exit 0, presented as this run's. Removing it first makes the check below real.
    recon.unlink(missing_ok=True)
    t0 = time.time()
    rc = summary_io.stage_entry(
        stages, "walker", t0,
        walker_cmd.run_walker(walker_argv, walker_bin=args.walker_bin),
        out / f"{queries[0].stem}.summary.json")
    if rc not in (0, EX_SOME_ERRORED):
        return _die(f"logan-walker (stage 1) failed (exit {rc}); see its output above.",
                    _stage_code(rc))
    # Exit 4 means some accessions errored while others reconstructed. Those haplotypes are
    # real, and discarding a whole cohort's alignment and tree because one accession in
    # thousands hit a network error would contradict the walker's own rule that one bad
    # accession does not stop the run. So the remaining stages run on what was produced, and
    # the walker's code is carried to this composite's own exit instead of being swallowed.
    exit_rc = 0
    if rc == EX_SOME_ERRORED:
        exit_rc = EX_SOME_ERRORED
        skipped = out / f"{queries[0].stem}.skipped.tsv"
        print(f"[reconstruct] warn: logan-walker exited {rc}: some accessions errored, others "
              f"reconstructed. Continuing with the haplotypes it did produce; this run will "
              f"still exit {rc}. The accessions that failed are listed in {skipped}.",
              file=sys.stderr)

    # The file was removed before the walker ran, so its absence here means the
    # walker exited 0 without writing one. That is not "nothing was reconstructed", it is
    # "that was not a reconstruction run", and it used to be indistinguishable from a run
    # whose output an earlier run had left behind.
    if not recon.is_file():
        return _die(f"logan-walker exited {rc} but wrote no {recon}; nothing to align. "
                    f"A walker flag passed after `--` can do this.")
    n_recon = sum(1 for line in recon.open() if line.startswith(">"))
    print(f"[reconstruct] reconstructed {n_recon} haplotype(s) -> {recon}", flush=True)
    if n_recon == 0:
        return _die("no haplotypes reconstructed; nothing to align.")

    sieved = out / "reconstructed.sieved.fasta"
    stem = msa_ref.stem
    aln = out / f"{stem}_codon_msa.fasta"
    tree = out / f"{stem}_fasttree.nwk"
    # The walker has just rewritten the reconstruction, so a sieved copy an earlier run left in
    # this --out-dir no longer matches it. With --sieve this run writes it again; kept, it sat
    # beside fresh results as if this run had made it.
    sieved.unlink(missing_ok=True)

    # ---- optional sieve ----------------------------------------------------
    msa_input = recon
    if args.sieve and not args.skip_tree:
        # The default --sieve-min-len follows the WALKER's own coverage gate, over the
        # shortest reference: each primary was reconstructed against its own anchor, which
        # need not be the MSA reference.
        #
        # Not a flat 95%: the walker publishes a primary from 85% of the reference and the
        # MSA admits a row with 70% of it in ACGT, so a flat gate here removes every
        # haplotype in between at a step documented as near-duplicate thinning, with exit 0
        # and no line in the summary. Those are the walks that stop short of a terminus and
        # the samples with a real in-frame deletion.
        #
        # No margin is added: `cov` is measured against the winning anchor, which is at
        # least as long as the shortest reference, so this cannot be stricter than the
        # walker was.
        min_len = args.sieve_min_len
        if min_len is None:
            rid, length = min(references, key=lambda r: r[1])
            frac = min(WALKER_DEFAULT_MIN_COV, _walker_min_cov(raw_extra))
            min_len = int(length * frac)
            print(f"[reconstruct] sieve --min-len {min_len}: {frac:.2f} of the shortest "
                  f"reference ({rid}, {length} nt), which is the walker's own --min-cov.",
                  flush=True)
        t0 = time.time()
        rc = summary_io.stage_entry(
            stages, "sieve", t0,
            sieve_cmd.sieve_run(
                recon, sieved,
                min_len=min_len,
                threshold=args.sieve_threshold,
                sketch=args.sieve_sketch,
                sieve_bin=args.sieve_bin,
            ),
            sieved.with_suffix(".summary.json"))
        if rc != 0:
            return _die(f"sieve failed (exit {rc}); see its output above.", _stage_code(rc))
        msa_input = sieved

    if args.skip_tree:
        # Nor does an MSA or tree an earlier run left, which would otherwise sit beside this
        # exit 0 as if this run had made it.
        for stale in (aln, tree):
            stale.unlink(missing_ok=True)
        print("[reconstruct] --skip-tree: done after reconstruction.")
        return exit_rc

    # ---- stage 2: codon MSA + FastTree ------------------------------------
    print(f"[reconstruct] 2/2 codon MSA + tree: --input {msa_input} --out-aln {aln} "
          f"--out-tree {tree}", flush=True)
    msa_args = stage_args(
        msa_cmd,
        ref=str(msa_ref),
        input=str(msa_input),
        out_aln=str(aln),
        out_tree=str(tree),
        mafft_threads=args.mafft_threads,
        msa_min_cov=args.msa_min_cov,
    )
    t0 = time.time()
    rc = summary_io.stage_entry(stages, "msa", t0, msa_cmd.run(msa_args),
                                aln.with_suffix(".summary.json"))
    if rc != 0:
        return _die(f"codon MSA / tree (stage 2) failed (exit {rc}); see its output above.",
                    _stage_code(rc))

    print(f"[reconstruct] done:\n"
          f"  haplotypes : {recon}\n"
          f"  codon MSA  : {aln}\n"
          f"  phylogeny  : {tree}")
    return exit_rc
