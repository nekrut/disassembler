#!/usr/bin/env python3
"""Strip a summary's volatile fields so two runs of the same work compare byte for byte.

`docs/summaries.md` promises that a summary is diffable once the fields that cannot be
reproduced are masked. This is that mask, used by `bench/check_golden.sh` to gate summaries
and usable by hand to diff two runs:

    python3 bench/mask_summary.py a.summary.json > a.masked
    python3 bench/mask_summary.py b.summary.json > b.masked
    diff a.masked b.masked

Masked, and why:

  run.started, run.wall_seconds, run.stage_seconds   a clock
  records[].timing                                   a clock
  tool.external[].path, tool.version, tool.disassembler
                                                     the machine, and the version under test
  inputs.source.id, inputs.source.kind               a cache or file path on this machine,
                                                     and whether the unitigs came from a
                                                     warm cache, which is a property of the
                                                     runner and not of the walk

Everything else is gated, including every count, every hash and every setting. Nothing is
reordered: the writers emit a fixed key order and that order is part of what is checked.
"""

import json
import sys

PLACEHOLDER = "<masked>"


def mask(doc):
    tool = doc.get("tool")
    if isinstance(tool, dict):
        for key in ("version", "disassembler"):
            if key in tool:
                tool[key] = PLACEHOLDER
        for ext in tool.get("external") or []:
            if isinstance(ext, dict) and "path" in ext:
                ext["path"] = PLACEHOLDER
    source = (doc.get("inputs") or {}).get("source")
    if isinstance(source, dict):
        for key in ("id", "kind"):
            if key in source:
                source[key] = PLACEHOLDER
    run = doc.get("run")
    if isinstance(run, dict):
        for key in ("started", "wall_seconds", "stage_seconds"):
            if key in run:
                run[key] = PLACEHOLDER
    for record in doc.get("records") or []:
        if isinstance(record, dict):
            for key in ("timing", "wall_seconds"):
                if key in record:
                    record[key] = PLACEHOLDER
    return doc


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: mask_summary.py <summary.json>", file=sys.stderr)
        return 2
    try:
        with open(sys.argv[1], encoding="utf-8") as fh:
            doc = json.load(fh)
    except (OSError, ValueError) as e:
        print(f"mask_summary: {sys.argv[1]}: {e}", file=sys.stderr)
        return 1
    json.dump(mask(doc), sys.stdout, indent=2, sort_keys=False)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
