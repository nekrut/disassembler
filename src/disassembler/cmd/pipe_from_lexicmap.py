"""from-lexicmap: LexicMap tabular stream -> (sieve opt) -> codon MSA + tree.

Composite: lexicmap -> sieve (optional) -> msa. Lexicmap already produces a
collapsed haplotype MSA (`<prefix>.clean.msa.fasta`); this composite runs the
same reference-guided MAFFT + FastTree pipeline the walker path uses, so
downstream tools see a consistent codon MSA / Newick tree.

Sieve/MSA input defaults to `<prefix>.clean.accessions.fasta` (uncollapsed
clean cohort, one seq per accession); pass --use-expanded to feed
`<prefix>.clean_expanded.accessions.fasta` (adds >=50% partials).
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


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "from-lexicmap",
        help="LexicMap tabular stream -> codon MSA + tree (with optional sieve).",
        description=(
            "Composite: stream lexicmap into a clean cohort, optionally sieve, then "
            "codon MSA + FastTree GTR tree."
        ),
    )

    # lexicmap inputs (mirror `disassembler lexicmap`)
    p.add_argument("-q", "-r", "--query", "--ref", dest="ref", required=True,
                   help="Reference in-frame CDS FASTA")
    p.add_argument("-d", "--dataset", default=None,
                   help="Galaxy dataset ID to stream from over HTTP")
    p.add_argument("-f", "-i", "--file", "--input", dest="file", default=None,
                   help="Local LexicMap tabular file (or '-' for stdin)")
    p.add_argument("--base-url", default="https://usegalaxy.org")
    p.add_argument("--api-key", default=None)
    p.add_argument("--api-key-file", default="/tmp/gxy.txt")
    p.add_argument("--min-cov", type=float, default=0.80,
                   help="Clean Full-Length coverage gate (lexicmap)")
    p.add_argument("--min-cov-partial", type=float, default=0.50,
                   help="Clean Partial coverage gate (lexicmap)")
    p.add_argument("--min-pident", type=float, default=60.0)
    p.add_argument("--max-internal-stops", type=int, default=0)
    p.add_argument("--allow-frameshift-alts", action="store_true",
                   help="Allow indels/frameshifts with premature stop codons into clean cohort")
    p.add_argument("--full-length-only", action="store_true")
    p.add_argument("--max-records", type=int, default=None)
    p.add_argument("--sample-cap", type=int, default=10)
    p.add_argument("--no-a2m", action="store_true")

    # composite plumbing
    p.add_argument("-o", "--out-dir", required=True,
                   help="Output directory. Lexicmap writes here under prefix `lexicmap`; "
                        "MSA + tree land alongside.")
    p.add_argument("--use-expanded", action="store_true",
                   help="Feed <prefix>.clean_expanded.accessions.fasta into sieve/MSA "
                        "instead of <prefix>.clean.accessions.fasta.")

    # sieve
    p.add_argument("--sieve", action="store_true",
                   help="Diversity-thin the cohort before the MSA/tree.")
    p.add_argument("--sieve-threshold", type=float, default=0.005)
    p.add_argument("--sieve-sketch", type=int, default=256)
    p.add_argument("--sieve-min-len", type=int, default=None)
    p.add_argument("--sieve-bin", default=None)

    # msa
    p.add_argument("--mafft-threads", type=int, default=8, help="MAFFT threads")
    p.add_argument("--msa-min-cov", type=float, default=0.70,
                   help="min coverage for a sequence to enter the MSA (default 0.70). "
                        "Distinct from lexicmap's --min-cov / --min-cov-partial above.")
    p.set_defaults(func=run)


def run(args) -> int:
    if not args.dataset and not args.file:
        print("from-lexicmap: must specify either -d/--dataset <GALAXY_ID> or "
              "-f/-i/--file <PATH> (or '-' for stdin).", file=sys.stderr)
        return 2

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / "lexicmap"

    if args.sieve and not resolve_sieve(args.sieve_bin):
        print("from-lexicmap: kmer-diversity-sieve not found. Build it "
              "(`cd extern/kmer-diversity-sieve && make serial`), or pass --sieve-bin.",
              file=sys.stderr)
        return 1

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
        full_length_only=args.full_length_only,
        max_records=args.max_records,
        sample_cap=args.sample_cap,
        out_prefix=str(prefix),
        no_a2m=args.no_a2m,
    )
    rc = lexicmap_cmd.run(lex_ns)
    if rc != 0:
        return rc

    cohort = (
        Path(f"{prefix}.clean_expanded.accessions.fasta") if args.use_expanded
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
        sieved = out_dir / f"{cohort.stem}.sieved.fasta"
        rc = sieve_cmd.sieve_run(
            cohort, sieved,
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
