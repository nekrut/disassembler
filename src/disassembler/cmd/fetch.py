"""fetch: pre-populate the unitigs cache without running any analysis.

Thin wrapper over the walker binary in `--fetch-only` mode. Downloads every requested
accession into `<cache-dir>/<ACC>.unitigs.fa.zst` (byte-identical to the S3 payload) and
exits. Skip-existing is implicit: the walker's own cache wedge returns the on-disk file
when present, so an already-cached accession is a fast no-op.

After a `fetch`, an offline `walker --unitigs-cache DIR --offline` run hits zero network.
"""

from __future__ import annotations

import argparse

import sys

from disassembler.cmd import walker as walker_cmd
from disassembler.validators import positive_int, validate_unitigs_cache


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "fetch",
        help="Download Logan unitigs into a local cache; no analysis, no query needed.",
        description=(
            "Pre-populates a `--unitigs-cache DIR` with the compressed unitig payloads for "
            "every requested accession, then exits. Feed the same cache dir to "
            "`disassembler walker --unitigs-cache DIR --offline` (or the composites) to run "
            "as many local analyses as you want against zero network."
        ),
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("-a", "--accession", help="a single SRA accession")
    src.add_argument("-l", "--accession-list", help="file with one accession per line")
    p.add_argument("--unitigs-cache", required=True, metavar="DIR",
                   help="Directory to write <ACC>.unitigs.fa.zst files into. Created if absent.")
    p.add_argument("--fetch-workers", type=positive_int, default=None,
                   help="Parallel S3 fetch workers (walker default 8)")
    p.add_argument("--max-size-mb", type=positive_int, default=None,
                   help="Refuse accessions above this compressed size (walker default 1024)")
    p.add_argument("--allow-empty", action="store_true",
                   help="Exit 0 even if some accessions fail. Off by default.")
    p.add_argument("--walker-bin", default=None,
                   help="Path to the logan-walker binary")
    p.set_defaults(func=run)


def run(args) -> int:
    # Fetch always writes into the cache, so offline=False is the right semantics
    # for the shared validator: only reject not-a-directory, permit missing (walker
    # creates it).
    cache_err = validate_unitigs_cache(args.unitigs_cache, offline=False)
    if cache_err:
        print(f"fetch: {cache_err}", file=sys.stderr)
        return 1
    argv: list[str] = ["--fetch-only", "--unitigs-cache", args.unitigs_cache]
    if args.accession:
        argv += ["-a", args.accession]
    else:
        argv += ["-l", args.accession_list]
    if args.fetch_workers is not None:
        argv += ["--fetch-workers", str(args.fetch_workers)]
    if args.max_size_mb is not None:
        argv += ["--max-size-mb", str(args.max_size_mb)]
    if args.allow_empty:
        argv += ["--allow-empty"]
    return walker_cmd.run_walker(argv, walker_bin=args.walker_bin)
