"""The summary contract, shared by every Python stage.

One writer, so the stages cannot drift apart. The rules, the block order and the reasoning
are in `docs/summaries.md`; this module is that document in code, and
`walker/src/io/output.rs` is the same contract on the Rust side.

Blocks, always in this order: `schema`, `tool`, `inputs`, `settings`, `run`, `discards`,
`cohort` (optional), `records`.

Records are written one per line, compact, so a summary stays diffable and a consumer can
`grep` a single record out of a large one. `stream_records` keeps that true for a stage whose
record count is unbounded: records go to a temporary file as they are produced and are spliced
in at the end, so the summary never costs more memory than the stage already uses.
"""

from __future__ import annotations

import contextlib
import hashlib
import bisect
import json
import os
import time
from pathlib import Path

# Bumped to 2 for the 0.2.0 release. Within a version fields are only added, and 0.2.0 did
# more than add: the walker's `n_bubbles` now counts allele sets rather than the (fork, merge)
# pairs they were found between, so the same key returns a different number for the same data,
# and `no_links` is a new `status` value a reader with a closed set does not know. Both are
# listed under "Migration" in docs/summaries.md.
SCHEMA = 2

# Fields a consumer must mask before byte-comparing two summaries of the same run.
VOLATILE = ("run.started", "run.wall_seconds", "run.stage_seconds", "tool.external[].path")


def seq_sha256(seq: str | bytes) -> str:
    """Lowercase hex sha256 of a sequence, uppercased first.

    How a summary names a sequence without carrying it. Matches `seq_sha256` in the walker, so
    a hash means the same thing on both sides of the suite.
    """
    if isinstance(seq, str):
        seq = seq.encode()
    return hashlib.sha256(seq.upper()).hexdigest()


def quantiles(values) -> dict | None:
    """`p10`, `p25`, `p50`, `p75`, `p90`, `mean`, to 4 decimals, or None when empty.

    The same five quantiles under the same names everywhere. 0.1.x `lexicmap` (now `tile`)
    called the median `p50_median` and the walker reported no quantiles at all, so neither of
    a pair of summaries could be read with the other's parser.
    """
    v = sorted(float(x) for x in values)
    if not v:
        return None

    def at(q: float) -> float:
        """Linear interpolation between the two neighbouring order statistics.

        This used to be `v[round((n - 1) * q)]`, which with banker's rounding returned
        one of the two middle values of an even cohort, sometimes the lower and sometimes the
        upper. The p50s drive the tiler's divergent-regime call, so coverages of 0.55 and 0.65
        (true median 0.60) reported 0.55 and `divergent`, while 0.50/0.55/0.65/0.70 reported
        0.65 and `conserved`. This is numpy's "linear" and
        `statistics.quantiles(method="inclusive")`.
        """
        pos = (len(v) - 1) * q
        lo = int(pos)
        hi = min(lo + 1, len(v) - 1)
        return v[lo] + (v[hi] - v[lo]) * (pos - lo)

    return {
        "p10": round(at(0.10), 4),
        "p25": round(at(0.25), 4),
        "p50": round(at(0.50), 4),
        "p75": round(at(0.75), 4),
        "p90": round(at(0.90), 4),
        "mean": round(sum(v) / len(v), 4),
    }


def quantiles_from_counts(counts, scale: float = 1.0) -> dict | None:
    """`quantiles` over a multiset given as {value: how many}, each value times `scale`.

    A stage that streams millions of records cannot keep one Python float per record
    just to report five quantiles; where the values are small integers, such as how many
    reference positions an accession covered, a Counter is bounded by the range instead and
    is exact. The result has to agree with `quantiles` on the same multiset, including the
    interpolation, so the order statistics are addressed by index rather than recomputed.
    """
    total = sum(counts.values())
    if total == 0:
        return None
    values = sorted(counts)
    cumulative, running = [], 0
    for v in values:
        running += counts[v]
        cumulative.append(running)

    def nth(i: int) -> float:
        """The i-th smallest value, 0-based."""
        return values[bisect.bisect_right(cumulative, i)] * scale

    def at(q: float) -> float:
        pos = (total - 1) * q
        lo = int(pos)
        hi = min(lo + 1, total - 1)
        return nth(lo) + (nth(hi) - nth(lo)) * (pos - lo)

    mean = sum(v * n for v, n in counts.items()) * scale / total
    return {
        "p10": round(at(0.10), 4),
        "p25": round(at(0.25), 4),
        "p50": round(at(0.50), 4),
        "p75": round(at(0.75), 4),
        "p90": round(at(0.90), 4),
        "mean": round(mean, 4),
    }


def tool(name: str, version: str, external: list[dict] | None = None) -> dict:
    """The `tool` block. `external` is every binary whose version changes the result."""
    from disassembler import __version__

    return {
        "name": name,
        "version": version,
        "disassembler": __version__,
        "external": external or [],
    }


def source(kind: str, ident: str, n_requested: int) -> dict:
    """The `inputs.source` block. `kind` is one of s3, file, galaxy, stdin, cache."""
    return {"kind": kind, "id": ident, "n_requested": int(n_requested)}


def run_block(
    status: str,
    started: float,
    wall_seconds: float,
    n_records: int,
    n_ok: int,
    stage_seconds: dict | None = None,
    exit_code: int | None = None,
    error: str | None = None,
) -> dict:
    """The `run` block. `started` is a UNIX timestamp; it is written as ISO 8601 UTC."""
    if status not in ("ok", "partial", "failed"):
        raise ValueError(f"status must be ok, partial or failed, not {status!r}")
    block = {
        "status": status,
        "exit_code": (1 if status == "failed" else 0) if exit_code is None else exit_code,
        "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(started)),
        "wall_seconds": round(wall_seconds, 3),
        "stage_seconds": stage_seconds or {},
        "n_records": int(n_records),
        "n_ok": int(n_ok),
    }
    if error is not None:
        block["error"] = error
    return block


def status_of(n_records: int, n_ok: int) -> str:
    """`ok` when everything worked, `failed` when nothing did, `partial` in between."""
    if n_records > 0 and n_ok == n_records:
        return "ok"
    return "failed" if n_ok == 0 else "partial"


def _block(key: str, value) -> str:
    body = json.dumps(value, indent=2, sort_keys=False)
    return f'  "{key}": ' + body.replace("\n", "\n  ")


def _record_line(record) -> str:
    return "    " + json.dumps(record, separators=(", ", ": "), sort_keys=False)


def write(
    path,
    *,
    tool: dict,
    inputs: dict,
    settings: dict,
    run: dict,
    discards: dict,
    cohort: dict | None = None,
    records=None,
    records_file=None,
) -> Path:
    """Write one summary. `records_file` is a JSON-lines file from `stream_records`.

    A stage calls this on **every** path including failure, so the absence of a summary is
    never how a caller learns that something went wrong (`docs/summaries.md` rule 4).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write("{\n")
        f.write(f'  "schema": {SCHEMA},\n')
        for key, value in (
            ("tool", tool),
            ("inputs", inputs),
            ("settings", settings),
            ("run", run),
            ("discards", discards),
        ):
            f.write(_block(key, value) + ",\n")
        if cohort is not None:
            f.write(_block("cohort", cohort) + ",\n")
        f.write('  "records": [\n')
        if records_file is not None:
            with open(records_file, encoding="utf-8") as src:
                prev = None
                for line in src:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    if prev is not None:
                        f.write(prev + ",\n")
                    prev = "    " + line
                if prev is not None:
                    f.write(prev + "\n")
        else:
            rows = list(records or [])
            for i, record in enumerate(rows):
                f.write(_record_line(record) + ("," if i + 1 < len(rows) else "") + "\n")
        f.write("  ]\n}\n")
    return path


class Records:
    """Records collected on disk while a stage streams, for `write(records_file=...)`.

    The point of the tiler route is that it does not hold its input, so holding one record
    per accession in memory to write at the end would undo that on a cohort of millions. Each
    record is appended as one JSON line and spliced into the summary at the end.

    The path is deterministic (`<out_prefix>.records.jsonl.tmp`) so that a failing run, which
    still owes a summary, can clean up after itself without having to carry a handle around.
    """

    def __init__(self, out_prefix: str):
        self.path = f"{out_prefix}.records.jsonl.tmp"
        Path(self.path).parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(self.path, "w", encoding="utf-8")

    def add(self, record) -> None:
        self._fh.write(json.dumps(record, separators=(", ", ": "), sort_keys=False) + "\n")

    def close(self) -> None:
        if not self._fh.closed:
            self._fh.close()

    def cleanup(self) -> None:
        """Close and remove. Safe to call twice, and on a path that was never written."""
        self.close()
        with contextlib.suppress(OSError):
            os.unlink(self.path)


def discard_records(out_prefix: str) -> None:
    """Remove a leftover record file for `out_prefix`, if a run died holding one."""
    with contextlib.suppress(OSError):
        os.unlink(f"{out_prefix}.records.jsonl.tmp")


def worst(statuses) -> str:
    """The worst of several stage statuses: failed beats partial beats ok."""
    seen = list(statuses)
    if any(s == "failed" for s in seen):
        return "failed"
    if any(s == "partial" for s in seen) or not seen:
        return "partial"
    return "ok"


def manifest(path, *, name: str, settings: dict, source: dict, stages: list,
             started: float, status: str | None = None, error: str | None = None) -> Path:
    """A composite's summary: one record per stage, each naming that stage's own summary.

    A composite runs three or four stages and, before 0.2.0, wrote nothing of its own: the
    run left per-stage files with nothing tying them together, no record of the sequence or
    the arguments, and no statement of which stage failed. This is the only place a summary
    points at another summary, and it points by relative path so the set can be moved.
    """
    status = status or worst(s.get("status", "ok") for s in stages)
    return write(
        path,
        tool=tool(name, __import__("disassembler").__version__),
        inputs={"references": None, "source": source},
        settings=settings,
        run=run_block(
            status=status,
            started=started,
            wall_seconds=time.time() - started,
            n_records=len(stages),
            n_ok=sum(1 for s in stages if s.get("status") == "ok"),
            stage_seconds={s["stage"]: s.get("wall_seconds", 0.0) for s in stages},
            error=error,
        ),
        discards={"stages": {"failed": sum(1 for s in stages if s.get("status") == "failed"),
                             "skipped": sum(1 for s in stages
                                            if s.get("status") == "skipped")}},
        records=stages,
    )


def stage_entry(stages: list, name: str, started: float, rc: int,
                summary_path=None, status: str | None = None) -> int:
    """Append one stage's record to a composite's manifest and return its exit code.

    `summary_path` is recorded only when the stage actually wrote one, and as a path relative
    to the manifest where possible, so a composite's outputs can be moved as a set.
    """
    entry = {
        "stage": name,
        "status": status or ("ok" if rc == 0 else "failed"),
        "exit_code": rc,
        "wall_seconds": round(time.time() - started, 3),
    }
    if summary_path is not None:
        path = Path(summary_path)
        if path.is_file():
            entry["summary"] = path.name
    stages.append(entry)
    return rc
