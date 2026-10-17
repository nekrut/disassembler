"""Top-level `disassembler` argparse dispatcher.

Each subcommand module under `disassembler.cmd` exposes two callables:

    def add_parser(subparsers) -> None
    def run(args) -> int

`add_parser` attaches the subcommand's own argparse to the parent's `subparsers`
and calls `set_defaults(func=run)`. `run` returns the process exit code.

Atomic subcommands map 1:1 to the current scripts. Composite subcommands
(`reconstruct`, `from-kmindex`, `from-lexicmap`) call atomics in-process so
tracebacks stay coherent and stage plumbing is normal Python.
"""

from __future__ import annotations

import argparse
import sys
from typing import Sequence

from disassembler import __version__
from disassembler.cmd import (
    fetch,
    harvest,
    lexicmap,
    msa,
    pipe_from_kmindex,
    pipe_from_lexicmap,
    pipe_reconstruct,
    sieve,
    stratify,
    walker,
)

# Order controls `--help` presentation: atomics first, then composites.
_ATOMIC = [harvest, stratify, fetch, walker, lexicmap, sieve, msa]
_COMPOSITE = [pipe_reconstruct, pipe_from_kmindex, pipe_from_lexicmap]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="disassembler",
        description=(
            "Reconstruct coding haplotypes from planetary sequence indexes. "
            "Seven atomic subcommands (harvest, stratify, fetch, walker, lexicmap, sieve, msa) "
            "and three composite pipelines (reconstruct, from-kmindex, from-lexicmap)."
        ),
    )
    # The walker binary versions separately: `disassembler walker --version`.
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True, metavar="<subcommand>")
    for mod in _ATOMIC + _COMPOSITE:
        mod.add_parser(sub)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    argv = list(argv)

    # `walker` is a passthrough to a foreign binary. argparse's REMAINDER refuses to
    # collect a leading `--help`/`-h`, so a plain `disassembler walker --help` would
    # bounce off the top-level parser as "unrecognized arguments". Bypass argparse
    # entirely for that subcommand: hand every argument after `walker` to walker.run.
    if argv and argv[0] == "walker":
        from disassembler.cmd import walker as walker_mod
        ns = argparse.Namespace(passthrough=argv[1:])
        rc = walker_mod.run(ns)
        return int(rc) if rc is not None else 0

    parser = build_parser()
    args = parser.parse_args(argv)
    rc = args.func(args)
    return int(rc) if rc is not None else 0


if __name__ == "__main__":
    sys.exit(main())
