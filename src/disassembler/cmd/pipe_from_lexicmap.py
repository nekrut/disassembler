"""from-lexicmap: LexicMap tabular stream -> (sieve opt) -> codon MSA + tree.

Composite: lexicmap -> sieve (optional) -> msa. Lexicmap already produces a
collapsed haplotype MSA (`<prefix>.clean.msa.fasta`); this composite runs the
same reference-guided MAFFT + FastTree pipeline the walker path uses, so
downstream tools see a consistent codon MSA / Newick tree.

Sieve/MSA input defaults to `<prefix>.clean.accessions.fasta` (uncollapsed
clean cohort, one seq per accession); pass --include-partials to feed
`<prefix>.clean_expanded.accessions.fasta` instead (adds >=50% partials).
"""

from __future__ import annotations

import argparse
import sys
import types
from pathlib import Path

from disassembler.cmd import lexicmap as lexicmap_cmd
from disassembler.cmd import msa as msa_cmd
from disassembler.cmd import sieve as sieve_cmd
from disassembler.resolve import resolve_sieve


def _byte_exact_dedup(src: Path, dst: Path) -> tuple[int, int]:
    """Write to `dst` the first-encountered record per unique sequence from `src`.

    Returns (n_in, n_out). Headers are preserved verbatim from the first copy;
    sequence lines are compared case-insensitively (upper-cased) with newline
    whitespace stripped, so cosmetic line-wrapping differences between two
    otherwise identical sequences are folded together.

    Preserving the ORDER of first occurrence is what makes the subsequent
    sieve output byte-identical to sieving `src` directly: the sieve is a
    deterministic in-order filter that would have rejected every later copy
    of an exact-duplicate sequence as novelty=0 anyway, so folding them out
    beforehand is observationally equivalent while strictly reducing work.
    """
    seen: set[str] = set()
    n_in = 0
    n_out = 0
    header: str | None = None
    seq_lines: list[str] = []

    def _emit(fout) -> None:
        nonlocal n_out
        if header is None:
            return
        seq = "".join(seq_lines).replace("\n", "").upper()
        if seq and seq not in seen:
            seen.add(seq)
            fout.write(header)
            fout.writelines(seq_lines)
            n_out += 1

    with src.open() as fin, dst.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                _emit(fout)
                header = line
                seq_lines = []
                n_in += 1
            else:
                seq_lines.append(line)
        _emit(fout)
    return n_in, n_out
import os

from disassembler.validators import (
    fasta_has_record,
    fraction,
    non_negative_int,
    percentage,
    positive_int,
    validate_api_key_file,
)


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "from-lexicmap",
        help="LexicMap tabular stream -> codon MSA + tree (with optional sieve).",
        description=(
            "Composite from a LexicMap tabular stream. Chains three atomic stages:\n"
            "  1. lexicmap : tile multi-HSP hits into a clean per-accession cohort\n"
            "  2. sieve    : optional diversity thinning (opt-in with --sieve)\n"
            "  3. msa      : codon MSA + FastTree GTR tree\n\n"
            "Options are grouped by pipeline stage below."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ---- Stage 1: lexicmap (input + cohort gating) ------------------------
    g_lex = p.add_argument_group(
        "LexicMap: input & gating",
        description="Options for the lexicmap streaming stage. Exactly one of "
                    "-d/--dataset and -f/-i/--file is required. See "
                    "`disassembler lexicmap -h` for the underlying atomic.",
    )
    g_lex.add_argument("-q", "-r", "--query", "--ref", dest="ref", required=True,
                       help="Reference in-frame CDS FASTA. Only the first record is read; "
                            "its length defines the coordinate grid every HSP is projected "
                            "onto. Also used as the MSA reference downstream.")
    g_lex.add_argument("-d", "--dataset", default=None,
                       help="Galaxy dataset ID to stream from over HTTP. Mutually exclusive "
                            "with -f/-i/--file.")
    g_lex.add_argument("-f", "-i", "--file", "--input", dest="file", default=None,
                       help="Local LexicMap tabular file, or '-' to read from stdin. "
                            "Mutually exclusive with -d/--dataset.")
    g_lex.add_argument("--base-url", default="https://usegalaxy.org",
                       help="Galaxy instance base URL used with -d/--dataset "
                            "(default: https://usegalaxy.org).")
    g_lex.add_argument("--api-key", default=None,
                       help="Galaxy API key as a literal string. If unset, read from "
                            "--api-key-file, then $GALAXY_API_KEY, then /tmp/gxy.txt.")
    g_lex.add_argument("--api-key-file", default="/tmp/gxy.txt",
                       help="Path to a file whose contents are the Galaxy API key "
                            "(default: /tmp/gxy.txt). Ignored if --api-key is set.")
    g_lex.add_argument("--min-cov", type=fraction, default=0.80,
                       help="Clean full-length cohort gate: minimum fraction of the "
                            "reference an accession's merged HSPs must cover "
                            "(default: 0.80). Applied at the lexicmap stage; distinct "
                            "from --msa-min-cov below.")
    g_lex.add_argument("--min-cov-partial", type=fraction, default=0.50,
                       help="Clean partial cohort gate: minimum coverage fraction for the "
                            "expanded clean cohort (default: 0.50). Must be <= --min-cov. "
                            "Anything below this is filed as a low-coverage fragment. "
                            "Applied at the lexicmap stage.")
    g_lex.add_argument("--min-pident", type=percentage, default=60.0,
                       help="Minimum percent identity for an HSP to contribute to the "
                            "tiled reconstruction (default: 60.0). HSPs below this are "
                            "dropped before the coverage gates.")
    g_lex.add_argument("--max-internal-stops", type=non_negative_int, default=0,
                       help="Maximum premature stop codons tolerated in the CDS body "
                            "before an accession is flagged out of the clean cohort "
                            "(default: 0).")
    g_lex.add_argument("--allow-frameshift-alts", action="store_true",
                       help="Let accessions with indels/frameshifts that produce "
                            "premature stops into the clean cohort anyway. Off by default.")
    g_lex.add_argument("--require-full-length", action="store_true",
                       help="Require 100%% reference coverage (every position tiled) for "
                            "the clean full-length tier, instead of the --min-cov "
                            "threshold. Off by default.")
    g_lex.add_argument("--max-records", type=positive_int, default=None,
                       help="Stop after processing N input rows. Useful for a quick "
                            "preview on a large tabular stream.")
    g_lex.add_argument("--samples-per-haplotype", type=positive_int, default=10,
                       dest="samples_per_haplotype",
                       help="Maximum representative accessions listed per collapsed "
                            "haplotype in <prefix>.clean.haplotypes.tsv (default: 10). "
                            "Does not affect the .accessions.fasta outputs.")
    g_lex.add_argument("--no-a2m", action="store_true",
                       help="Skip writing the .a2m formatted MSA files "
                            "(<prefix>.clean.a2m and the .a2m compat alias).")
    g_lex.add_argument("--allow-empty", action="store_true",
                       help="Exit 0 even if the lexicmap clean cohort is empty. Off by "
                            "default: an empty cohort exits non-zero so a wrapping "
                            "workflow does not mistake a whole-target miss for a good run.")

    # ---- Cohort selection into sieve/MSA ----------------------------------
    g_cohort = p.add_argument_group(
        "Cohort selection",
        description="Where lexicmap writes and which of its two clean cohort FASTAs "
                    "is fed to the downstream sieve/MSA stages.",
    )
    g_cohort.add_argument("-o", "--out-dir", required=True,
                          help="Output directory. Lexicmap writes here under the prefix "
                               "`lexicmap`; the MSA and tree land alongside.")
    g_cohort.add_argument("--include-partials", action="store_true",
                          dest="include_partials",
                          help="Feed <prefix>.clean_expanded.accessions.fasta (adds >=50%% "
                               "partial hits) into sieve/MSA instead of the default "
                               "<prefix>.clean.accessions.fasta (full-length clean only).")

    # ---- Stage 2: sieve (optional diversity thinning) ---------------------
    g_sieve = p.add_argument_group(
        "Sieve: optional diversity thinning",
        description="Optional pass between lexicmap and MSA, delegated to "
                    "`disassembler sieve`. Haplotype/lineage-level MinHash/LSH reducer; "
                    "all --sieve-* flags are ignored unless --sieve is passed.",
    )
    g_sieve.add_argument("--sieve", action="store_true",
                         help="Enable the diversity-thinning pass before the MSA/tree. "
                              "Off by default.")
    g_sieve.add_argument("--sieve-threshold", type=fraction, default=0.005,
                         help="LSH promotion threshold tau (default 0.005, gene-scale). "
                              "Higher tau collapses more.")
    g_sieve.add_argument("--sieve-sketch", type=positive_int, default=256,
                         help="MinHash sketch size m (default 256, gene-scale). Caps the "
                              "resolution of the novelty estimator at 1/m.")
    g_sieve.add_argument("--sieve-min-len", type=positive_int, default=None,
                         help="Minimum sequence length gate for the sieve (default: auto "
                              "= 95%% of the reference CDS length).")
    g_sieve.add_argument("--sieve-bin", default=None,
                         help="Explicit path to the kmer-diversity-sieve l1_pipeline_serial "
                              "binary (auto-discovered otherwise).")

    # ---- Stage 3: MSA + tree ----------------------------------------------
    g_msa = p.add_argument_group(
        "MSA + tree",
        description="Options for the codon-MSA + FastTree GTR stage. See "
                    "`disassembler msa -h`.",
    )
    g_msa.add_argument("--mafft-threads", type=positive_int, default=8,
                       help="Thread count for MAFFT in the MSA stage (default: 8).")
    g_msa.add_argument("--msa-min-cov", type=fraction, default=0.70,
                       help="Minimum coverage fraction for a sequence to enter the MSA "
                            "(default: 0.70). Distinct from lexicmap's --min-cov / "
                            "--min-cov-partial above; those gate which accessions reach "
                            "the clean cohort at all.")
    p.set_defaults(func=run)


def run(args) -> int:
    if not args.dataset and not args.file:
        print("from-lexicmap: must specify either -d/--dataset <GALAXY_ID> or "
              "-f/-i/--file <PATH> (or '-' for stdin).", file=sys.stderr)
        return 2

    if args.min_cov_partial > args.min_cov:
        print(f"from-lexicmap: --min-cov-partial ({args.min_cov_partial}) must be "
              f"<= --min-cov ({args.min_cov}); otherwise the partial tier is "
              f"empty by construction.", file=sys.stderr)
        return 2

    if not Path(args.ref).is_file():
        print(f"from-lexicmap: reference FASTA not found: {args.ref}", file=sys.stderr)
        return 1
    if fasta_has_record(args.ref) == 0:
        print(f"from-lexicmap: reference FASTA has no records (no `>` header lines): "
              f"{args.ref}. Point -q at a CDS FASTA, not an empty or non-FASTA file.",
              file=sys.stderr)
        return 1

    # Galaxy API key file: only validate when streaming from a dataset AND we're
    # actually going to fall back to the file (matches get_galaxy_api_key's
    # resolution order in lexicmap.py: --api-key > env var > file).
    if args.dataset and not args.api_key and not os.environ.get("GALAXY_API_KEY"):
        key_err = validate_api_key_file(args.api_key_file)
        if key_err:
            print(f"from-lexicmap: {key_err}", file=sys.stderr)
            return 1

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / "lexicmap"

    if args.sieve and not resolve_sieve(args.sieve_bin):
        print("from-lexicmap: kmer-diversity-sieve not found. Build it "
              "(`cd extern/kmer-diversity-sieve && make serial`), or pass --sieve-bin.",
              file=sys.stderr)
        return 1

    # Warn on sieve knobs whose stage will not run. Defaults kept in sync with
    # add_parser above.
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
            print(f"[from-lexicmap] warn: {', '.join(sieve_inert)} set but --sieve is "
                  f"off, so these values are inert. Add --sieve to make them take effect.",
                  file=sys.stderr)

    # ---- 1/3 lexicmap ------------------------------------------------------
    print(f"[from-lexicmap] 1/3 lexicmap -> {prefix}.*", flush=True)
    lex_ns = types.SimpleNamespace(
        ref=args.ref,
        dataset=args.dataset,
        file=args.file,
        base_url=args.base_url,
        api_key=args.api_key,
        api_key_file=args.api_key_file,
        min_cov=args.min_cov,
        min_cov_partial=args.min_cov_partial,
        min_pident=args.min_pident,
        max_internal_stops=args.max_internal_stops,
        allow_frameshift_alts=args.allow_frameshift_alts,
        require_full_length=args.require_full_length,
        max_records=args.max_records,
        samples_per_haplotype=args.samples_per_haplotype,
        out_prefix=str(prefix),
        no_a2m=args.no_a2m,
        allow_empty=args.allow_empty,
    )
    rc = lexicmap_cmd.run(lex_ns)
    if rc != 0:
        return rc

    cohort = (
        Path(f"{prefix}.clean_expanded.accessions.fasta") if args.include_partials
        else Path(f"{prefix}.clean.accessions.fasta")
    )
    if not cohort.is_file():
        print(f"from-lexicmap: expected lexicmap cohort not found: {cohort}",
              file=sys.stderr)
        return 1
    n_cohort = sum(1 for line in cohort.open() if line.startswith(">"))
    print(f"[from-lexicmap] cohort: {n_cohort} sequences from {cohort.name}", flush=True)
    if n_cohort == 0:
        print("from-lexicmap: cohort is empty; nothing to align.", file=sys.stderr)
        return 1

    msa_input = cohort

    # ---- 2/3 sieve (optional) ---------------------------------------------
    if args.sieve:
        # Byte-exact dedup before the sieve. The sieve is a MinHash/LSH filter
        # that iterates records in order and keeps a sequence only when its
        # novelty vs. the already-promoted set is >= tau; a byte-exact later
        # copy of a sequence has novelty 0 and is always rejected. Folding
        # those duplicates out ahead of time is observationally equivalent to
        # the sieve's own fold (see tests/test_dedup_sieve_equivalence.py)
        # and strictly reduces the sieve's work.
        deduped = out_dir / f"{cohort.stem}.deduped.fasta"
        n_in, n_out = _byte_exact_dedup(cohort, deduped)
        if n_in == 0:
            print("from-lexicmap: cohort is empty after dedup; nothing to sieve.",
                  file=sys.stderr)
            return 1
        fold = n_in / n_out if n_out else 0.0
        print(f"[from-lexicmap] dedup: {n_in} -> {n_out} records "
              f"({fold:.2f}x fold) before sieve", flush=True)
        sieved = out_dir / f"{cohort.stem}.sieved.fasta"
        rc = sieve_cmd.sieve_run(
            deduped, sieved,
            ref=Path(args.ref),
            min_len=args.sieve_min_len,
            threshold=args.sieve_threshold,
            sketch=args.sieve_sketch,
            sieve_bin=args.sieve_bin,
        )
        if rc != 0:
            return rc
        msa_input = sieved
    else:
        print("[from-lexicmap] 2/3 sieve skipped (no --sieve)", flush=True)

    # ---- 3/3 codon MSA + FastTree -----------------------------------------
    stem = Path(args.ref).stem
    aln = out_dir / f"{stem}_codon_msa.fasta"
    tree = out_dir / f"{stem}_fasttree.nwk"
    print(f"[from-lexicmap] 3/3 codon MSA + tree: --input {msa_input} --out-aln {aln} "
          f"--out-tree {tree}", flush=True)
    msa_ns = types.SimpleNamespace(
        ref=args.ref,
        input=str(msa_input),
        out_aln=str(aln),
        out_tree=str(tree),
        mafft_threads=args.mafft_threads,
        msa_min_cov=args.msa_min_cov,
    )
    rc = msa_cmd.run(msa_ns)
    if rc != 0:
        return rc

    print(f"[from-lexicmap] done:\n"
          f"  cohort   : {cohort}\n"
          f"  MSA      : {aln}\n"
          f"  phylogeny: {tree}")
    return 0
