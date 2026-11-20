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
    tiler,
    msa,
    pipe_from_kmindex,
    pipe_from_lexicmap,
    pipe_reconstruct,
    sieve,
    stratify,
    walker,
)

# Order controls `--help` presentation: atomics first, then composites.
_ATOMIC = [harvest, stratify, fetch, walker, tiler, sieve, msa]
_COMPOSITE = [pipe_reconstruct, pipe_from_kmindex, pipe_from_lexicmap]

# Subcommand names from 0.1.x, mapped to what they are called now. See `main`.
_RENAMED = {"walker": "walk", "lexicmap": "tile"}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="disassembler",
        description=(
            "Reconstruct coding haplotypes from planetary sequence indexes. "
            "Seven atomic subcommands (harvest, stratify, fetch, walk, tile, sieve, msa) "
            "and three composite pipelines (reconstruct, from-kmindex, from-lexicmap)."
        ),
    )
    # The walker binary versions separately: `disassembler walk --version`.
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True, metavar="<subcommand>")
    for mod in _ATOMIC + _COMPOSITE:
        mod.add_parser(sub)
    return p


def main(argv: Sequence[str] | None = None) -> int:
    # Progress lines carry a U+2713 check mark. With PYTHONIOENCODING=ascii or a
    # legacy single-byte locale that print raised UnicodeEncodeError, in `msa` after MAFFT
    # had run and before the alignment was published, so the `finally` deleted it and the
    # run exited 1 having thrown away the work.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="backslashreplace")
        except (AttributeError, ValueError):
            pass

    if argv is None:
        argv = sys.argv[1:]
    argv = list(argv)

    # 0.1.x named two subcommands after things rather than actions, and one of those
    # things (`lexicmap`) is somebody else's program. The old names still run, with a
    # notice on stderr, and are kept out of `--help`.
    if argv and argv[0] in _RENAMED:
        new = _RENAMED[argv[0]]
        print(f"disassembler: `{argv[0]}` is now `{new}`; the old name still works but "
              f"will be removed in a later release.", file=sys.stderr)
        argv[0] = new

    # `walk` is a passthrough to a foreign binary. argparse's REMAINDER refuses to
    # collect a leading `--help`/`-h`, so a plain `disassembler walk --help` would
    # bounce off the top-level parser as "unrecognized arguments". Bypass argparse
    # entirely for that subcommand: hand every argument after `walk` to walker.run.
    if argv and argv[0] == "walk":
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
