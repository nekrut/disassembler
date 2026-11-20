"""fetch: pre-populate the unitigs cache without running any analysis.

Thin wrapper over the walker binary in `--fetch-only` mode. Downloads every requested
accession into `<cache-dir>/<ACC>.unitigs.fa.zst` (byte-identical to the S3 payload) and
exits. Skip-existing is implicit: the walker's own cache wedge returns the on-disk file
when present, so an already-cached accession is a fast no-op.

After a `fetch`, an offline `walk --unitigs-cache DIR --offline` run hits zero network.

The walker writes `fetch.summary.json` recording the size, sha256 and (for anything it
downloaded) the S3 ETag and Last-Modified of every object. That file goes into the cache
directory by default, beside the payloads it describes, so a filled cache is a frozen input
that can be cited rather than a bare directory of files.
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
            "`disassembler walk --unitigs-cache DIR --offline` (or the composites) to run "
            "as many local analyses as you want against zero network."
        ),
    )
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("-a", "--accession", help="a single SRA accession")
    src.add_argument("-l", "--accession-list", help="file with one accession per line")
    p.add_argument("--unitigs-cache", required=True, metavar="DIR",
                   help="Directory to write <ACC>.unitigs.fa.zst files into. Created if absent.")
    p.add_argument("--out-dir", default=None, metavar="DIR",
                   help="Where to write fetch.summary.json. Defaults to --unitigs-cache, so "
                        "the manifest sits beside the payloads it describes.")
    p.add_argument("--fetch-workers", type=positive_int, default=None,
                   help="Parallel S3 fetch workers (walker default 8)")
    p.add_argument("--max-size-mb", type=positive_int, default=None,
                   help="Refuse accessions above this compressed size (walker default 4096)")
    p.add_argument("--allow-empty", action="store_true",
                   help="Exit 0 even if some accessions could not be cached because Logan "
                        "has no data for them or they are over --max-size-mb. Any other "
                        "failure still exits non-zero. Off by default. Either way the cache "
                        "records which of the two it was, so the --offline run that follows "
                        "reports the same status instead of an error.")
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
    # Without an explicit -o the walker would default to ./walker_output and create a
    # directory the user never asked for, just to hold one file about the cache.
    out_dir = args.out_dir if args.out_dir is not None else args.unitigs_cache
    argv: list[str] = ["--fetch-only", "--unitigs-cache", args.unitigs_cache,
                       "-o", out_dir]
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
