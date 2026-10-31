"""Build the arguments a composite hands to one of its stages.

A composite calls each atomic's `run(args)` in-process. Building that `args` by hand
broke whenever the atomic gained an option: `from-kmindex` listed harvest's options
itself, harvest later grew `--shards-dir`, and every harvest through the composite died
with an AttributeError before reaching Galaxy. `stage_args` starts from the atomic's own
parser instead, so an option the composite does not set gets the atomic's default.
"""

from __future__ import annotations

import argparse


def stage_args(module, **values) -> argparse.Namespace:
    """The namespace `module.run` would get from its own parser, with `values` set.

    Raises TypeError if `values` names an option the stage does not have (a rename on
    either side fails here, not as a silently ignored setting), or leaves out one the
    stage requires.
    """
    subparsers = argparse.ArgumentParser(prog="disassembler", add_help=False).add_subparsers()
    module.add_parser(subparsers)
    (parser,) = subparsers.choices.values()

    actions = [a for a in parser._actions if a.default is not argparse.SUPPRESS]
    ns = {a.dest: a.default for a in actions}
    unknown = sorted(set(values) - set(ns))
    if unknown:
        raise TypeError(f"{parser.prog}: no option(s) {', '.join(unknown)}")
    ns.update(values)

    missing = [a.dest for a in actions if a.required and ns[a.dest] is None]
    for group in parser._mutually_exclusive_groups:
        if group.required and all(ns.get(a.dest) is None for a in group._group_actions):
            missing.append("/".join(a.dest for a in group._group_actions))
    if missing:
        raise TypeError(f"{parser.prog}: missing required {', '.join(missing)}")

    ns.update(parser._defaults)
    return argparse.Namespace(**ns)
