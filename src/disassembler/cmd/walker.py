"""walker: run the Rust `logan-walker` binary.

Every argument after `disassembler walker` is forwarded to `logan-walker` verbatim,
including `--help`. The only intercepted flag is `--walker-bin <PATH>`, which the
Rust binary does not accept and which lets a caller override the resolution search.

`-q` is repeatable — pass it more than once, or hand a multi-record FASTA, to invoke
the walker's multi-anchor mode. See `logan-walker --help` for details on the
anchor pool union and per-primary winning-reference tagging.

Binary resolution follows `disassembler.resolve.resolve_walker`:
  1. --walker-bin
  2. $DISASSEMBLER_WALKER_BIN
  3. PATH
  4. walker/target/release/logan-walker under the repo root

Not re-parsed in Python. `disassemble.py` used to keep a second copy of the walker's
defaults; changing a walker default there left the wrapper silently pinning the old
one. Passthrough removes the whole class of drift.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

from disassembler.resolve import resolve_walker


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "walker",
        help="Walk Logan unitigs into coding haplotypes (Rust binary passthrough).",
        description=(
            "Forwards every argument to the `logan-walker` Rust binary. "
            "Run `disassembler walker --help` to see the walker's own CLI. "
            "The only intercepted flag is `--walker-bin <PATH>`."
        ),
        add_help=False,  # let the walker print its own --help
    )
    p.add_argument("passthrough", nargs=argparse.REMAINDER,
                   help="Arguments forwarded to logan-walker")
    p.set_defaults(func=run)


def _resolved_or_report(walker_bin: str | None) -> str | None:
    resolved = resolve_walker(walker_bin)
    if not resolved:
        print(
            "walker: logan-walker not found. Build it "
            "(`cd walker && cargo build --release`), put it on PATH, or pass "
            "--walker-bin. Resolution order: --walker-bin, "
            "$DISASSEMBLER_WALKER_BIN, PATH, walker/target/release/logan-walker.",
            file=sys.stderr,
        )
    return resolved


def run_walker(argv: list[str], walker_bin: str | None = None) -> int:
    """Run logan-walker as a subprocess. Returns its exit code.

    Composites use this so subsequent stages still get to run. The top-level
    `walker` subcommand execv's instead (see `run`) so signals reach the caller
    unmodified.
    """
    resolved = _resolved_or_report(walker_bin)
    if not resolved:
        return 127
    return subprocess.run([resolved, *argv]).returncode


def run(args) -> int:
    argv = list(args.passthrough)
    walker_bin: str | None = None
    if "--walker-bin" in argv:
        i = argv.index("--walker-bin")
        if i + 1 >= len(argv):
            print("walker: --walker-bin needs a value", file=sys.stderr)
            return 2
        walker_bin = argv[i + 1]
        del argv[i:i + 2]

    resolved = _resolved_or_report(walker_bin)
    if not resolved:
        return 127

    # execvp when possible so signals and exit codes reach the caller unmodified.
    if os.name == "posix":
        os.execv(resolved, [resolved, *argv])
        return 0  # unreachable
    return subprocess.run([resolved, *argv]).returncode
