"""stratify: cut a size-capped cohort from harvested kmindex scores.

`--score-bin LO:HI:N` (repeatable): sample N accessions with score in the
half-open `[LO, HI)`, except that a bin ending at 1.0 includes 1.0, the highest
containment score. Nonzero scores only. No backfill: a bin short on candidates
contributes what it has.

Score files are recognised by content, never by name. `harvest` writes
`{"target": ..., "scores": {accession: score}}`; stratify reads the files it is
given, or every file in a directory, and takes the target from inside. A file in
the shape written before 0.1.1, a bare `{accession: score}`, carries no target,
so the caller names it (`--scores NAME=FILE`, or `--target` when it is the only
one). Output files are still named after the target; a name used as a label is
fine, one used to decide what a file is was not.

The phiX174 Illumina spike-in is filtered upstream in `harvest` by target NAME,
not by score. See `harvest.safe_target_name`.

An earlier `--diversity N` mode also existed (four containment brackets,
20/30/30/20 quotas with backfill). It was removed because it capped by count
while download cost is byte-driven, so it was the wrong lever for the download-
volume problem it was meant to solve; `walker --max-size-mb` gates bytes.
Calibration-experiment workflows that used `--diversity` can express the same
shape with four `--score-bin` invocations.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import NamedTuple

from disassembler.textio import open_text

# kmindex containment scores run from 0 to 1. A bin whose upper bound is exactly this includes
# it: `0.9:1.0` used to miss every perfect match, 274,808 of the 5.2 million accessions in the
# SARS-CoV-2 spike harvest behind the 20,000-accession walker run.
TOP_SCORE = 1.0


def stratified_sample(items, k):
    """Evenly samples k items from a list sorted by score.

    A quota of zero is legitimate: `quota_canon` is `int(total_target * 0.20)`,
    which is 0 for any --diversity under 5. The old guard was `len(items) <= k`,
    false whenever the bracket was non-empty, so the run reached `len(items) / k`
    and died with ZeroDivisionError. The bracket contributing nothing is the
    right answer at that size, and `stratify_diversity`'s top-up then fills the
    cohort from the pooled remainder.
    """
    if k <= 0:
        return []
    if len(items) <= k:
        return items
    step = len(items) / k
    return [items[int(i * step)] for i in range(k)]


def parse_bin_spec(s: str):
    """Parse `LO:HI:N` into `(lo, hi, n)`. Raises ValueError with the offending spec.

    Half-open `[lo, hi)` intervals match the bracket convention (`>= 0.99`,
    `< 0.99`), except at the top (see `in_bin`). N is the count requested for the
    bin; N == 0 is legal and useful for "reserve a bin boundary but do not sample
    from it this run".
    """
    parts = s.split(":")
    if len(parts) != 3:
        raise ValueError(f"--score-bin '{s}' must be LO:HI:N (three colon-separated fields)")
    try:
        lo, hi, n = float(parts[0]), float(parts[1]), int(parts[2])
    except ValueError as e:
        raise ValueError(f"--score-bin '{s}': could not parse as LO:HI:N ({e})")
    if not (lo < hi):
        raise ValueError(f"--score-bin '{s}': LO ({lo}) must be strictly less than HI ({hi})")
    if n < 0:
        raise ValueError(f"--score-bin '{s}': N ({n}) must be non-negative")
    return (lo, hi, n)


def in_bin(score: float, lo: float, hi: float) -> bool:
    """`[lo, hi)`, or `[lo, 1.0]` when `hi` is 1.0, the highest score there is."""
    return lo <= score < hi or (hi == TOP_SCORE and score >= TOP_SCORE)


def bin_label(lo: float, hi: float) -> str:
    return f"[{lo}, {hi}]" if hi == TOP_SCORE else f"[{lo}, {hi})"


def assert_bins_disjoint(bins) -> None:
    """Raise ValueError if any two half-open bins overlap.

    Overlap is a hard refusal in the current CLI (--allow-overlap was removed with
    --diversity). If overlapping bins are intentional, invoke stratify twice on
    disjoint subsets and concatenate the outputs.
    """
    ordered = sorted(bins, key=lambda b: b[0])
    for (a_lo, a_hi, _), (b_lo, b_hi, _) in zip(ordered, ordered[1:]):
        # A bin ending at 1.0 includes 1.0, so it also overlaps one starting there.
        if b_lo < a_hi or (a_hi == TOP_SCORE and b_lo == TOP_SCORE):
            raise ValueError(
                f"--score-bin ranges overlap: {bin_label(a_lo, a_hi)} and "
                f"{bin_label(b_lo, b_hi)}. "
                f"Overlap is refused; run stratify once per disjoint bin set if "
                f"you need the same accession in more than one bin's sample."
            )


def sample_bins_from_scores(scores, bins):
    """Sample per bin. Zero-score accessions are excluded up front.

    Returns `(picks, per_bin)`; `picks` is a merged list `(acc, score)` sorted
    by descending score, `per_bin` is one dict per bin with counts and selections.

    No backfill: a bin with fewer than N candidates contributes what it has.
    """
    nonzero = {a: s for a, s in scores.items() if s > 0.0}
    picks = []
    per_bin = []
    for lo, hi, n in bins:
        in_range = sorted(
            ((a, s) for a, s in nonzero.items() if in_bin(s, lo, hi)),
            key=lambda x: -x[1],
        )
        selected = stratified_sample(in_range, n)
        per_bin.append({
            "lo": lo,
            "hi": hi,
            "requested": n,
            "available": len(in_range),
            "selected": len(selected),
            "accessions": [{"acc": a, "score": s} for a, s in selected],
        })
        picks.extend(selected)
    picks.sort(key=lambda x: -x[1])
    return picks, per_bin


class ScoreFile(NamedTuple):
    path: Path
    target: str | None  # None: the pre-0.1.1 shape, which has no target inside
    scores: dict


def _is_score_map(d) -> bool:
    return isinstance(d, dict) and all(
        isinstance(k, str) and isinstance(v, (int, float)) and not isinstance(v, bool)
        for k, v in d.items())


def read_score_file(path: Path) -> ScoreFile | None:
    """The harvest scores in `path`, or None if it holds something else. Plain or gzipped.

    `{"target": T, "scores": {accession: score}}` is what `harvest` writes; a bare, non-empty
    `{accession: score}` is the shape written before 0.1.1 (and by hand), with no target.
    Raises if the file cannot be read or is not JSON.
    """
    with open_text(path) as fh:
        start = fh.read(256)
        if not start.lstrip().startswith("{"):
            raise ValueError(f"{path}: not a JSON object")
        data = json.loads(start + fh.read())
    if (isinstance(data, dict) and isinstance(data.get("target"), str)
            and _is_score_map(data.get("scores"))):
        return ScoreFile(path, data["target"], data["scores"])
    if data and _is_score_map(data):
        return ScoreFile(path, None, data)
    return None


def find_score_files(sources, target: str | None = None) -> list[ScoreFile]:
    """The score files among `sources`, each with its target. Raises ValueError if there are
    none, or if a target cannot be told.

    A source is a file, `NAME=FILE` (naming the target of a file that has none inside), or a
    directory, where every regular, non-hidden file is read and those that are not score files
    (a harvest directory also holds accession lists and summaries) are passed over. With
    `target`, only that target is kept, and if no file names it but exactly one has no target
    inside, that one is taken to be it.
    """
    found, passed_over = [], 0
    for src in sources:
        name, path = None, Path(src)
        if "=" in str(src) and not path.exists():
            name, _, rest = str(src).partition("=")
            path = Path(rest)
        if path.is_dir():
            if name:
                raise ValueError(f"{src}: NAME= names a single score file, not a directory")
            files = (x for x in path.iterdir() if x.is_file() and not x.name.startswith("."))
            for f in sorted(files):
                try:
                    sf = read_score_file(f)
                except (OSError, ValueError):
                    sf = None
                if sf is None:
                    passed_over += 1
                else:
                    found.append(sf)
        elif path.is_file():
            sf = read_score_file(path)
            if sf is None:
                raise ValueError(f"{path}: not a harvest score file ({{\"target\": ..., "
                                 f"\"scores\": {{accession: score}}}}, or a bare "
                                 f"{{accession: score}})")
            if name:
                if sf.target not in (None, name):
                    raise ValueError(f"{src}: the file's own target is {sf.target!r}")
                sf = sf._replace(target=name)
            found.append(sf)
        else:
            raise ValueError(f"{path}: no such file or directory")

    if not found:
        where = ", ".join(str(x) for x in sources)
        other = f"; {passed_over} other file(s) there are not score files" if passed_over else ""
        raise ValueError(f"no harvest score files in {where}{other}")

    unnamed = [sf for sf in found if sf.target is None]
    if target is not None:
        named = [sf for sf in found if sf.target == target]
        if named:
            if unnamed:
                print(f"[warn] {len(unnamed)} score file(s) with no target inside were not "
                      f"used: {', '.join(str(sf.path) for sf in unnamed)}", file=sys.stderr)
            found = named
        elif len(unnamed) == 1:
            found = [unnamed[0]._replace(target=target)]
        elif unnamed:
            raise ValueError(f"no score file names target {target!r}, and "
                             f"{len(unnamed)} have no target inside; name the right one "
                             f"with --scores {target}=FILE")
        else:
            raise ValueError(f"no score file for target {target!r}; the files found are "
                             f"for {', '.join(sorted({sf.target for sf in found}))}")
    elif unnamed:
        raise ValueError(f"{len(unnamed)} score file(s) have no target inside (written before "
                         f"0.1.1, or by hand): {', '.join(str(sf.path) for sf in unnamed)}. "
                         f"Name each with --scores NAME=FILE, or re-run harvest.")

    seen = {}
    for sf in found:
        t = sf.target
        # The target names output files, so it must not reach outside the output directory.
        if not t or t in (".", "..") or "/" in t or "\\" in t:
            raise ValueError(f"{sf.path}: target {t!r} cannot be used as a file name")
        if t in seen:
            raise ValueError(f"target {t!r} is in two score files: {seen[t]} and {sf.path}")
        seen[t] = sf.path
    return found


def stratify_bins(out_dir: Path, bins, score_files: list[ScoreFile]) -> None:
    """Per-target bin stratification into `out_dir`. Emits `<target>.bins.txt`,
    `<target>.bins.json`, and a combined `bins_summary.json`.
    """
    print("\n" + "=" * 80)
    ranges = ", ".join(f"{bin_label(lo, hi)}x{n}" for lo, hi, n in bins)
    print(f" Bin stratification (nonzero scores only): {ranges}")
    print("=" * 80)

    out_dir.mkdir(parents=True, exist_ok=True)
    all_summaries = {}
    for sf in sorted(score_files, key=lambda f: f.target):
        target, scores = sf.target, sf.scores

        picks, per_bin = sample_bins_from_scores(scores, bins)

        out_txt = out_dir / f"{target}.bins.txt"
        out_json = out_dir / f"{target}.bins.json"
        out_txt.write_text(("\n".join(x[0] for x in picks) + "\n") if picks else "")
        out_json.write_text(json.dumps({
            "target": target,
            "n_available_nonzero": sum(1 for s in scores.values() if s > 0.0),
            "n_total_input": len(scores),
            "total_selected": len(picks),
            "bins": per_bin,
        }, indent=2))

        bin_str = ", ".join(
            f"{bin_label(b['lo'], b['hi'])}: {b['selected']}/{b['requested']}" for b in per_bin
        )
        print(f"  [✓] Target {target:12}: {len(picks):5d} selected | {bin_str}")
        all_summaries[target] = {
            "total_selected": len(picks),
            "bins": per_bin,
            "n_available_nonzero": sum(1 for s in scores.values() if s > 0.0),
        }

    summary_file = out_dir / "bins_summary.json"
    summary_file.write_text(json.dumps(all_summaries, indent=2))
    print(f"[✓] Bins summary saved: {summary_file}")


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "stratify",
        help="Cut a size-capped cohort from harvested kmindex scores.",
        description=(
            "Sample N accessions per --score-bin from harvest score files. Bins are half-open "
            "[LO, HI), except that a bin ending at 1.0 includes 1.0; nonzero scores only, no "
            "backfill, disjoint. Score files are recognised by content, and each names its "
            "own target."
        ),
    )
    p.add_argument("-o", "--out-dir", default="candidates",
                   help="Output directory for <target>.bins.txt, <target>.bins.json and a "
                        "combined bins_summary.json. Also where the score files are read "
                        "from when --scores is not given (pass harvest's -o).")
    p.add_argument("-s", "--scores", action="append", metavar="[NAME=]PATH",
                   help="A harvest score file, or a directory: every file in it is read, and "
                        "the score files among them are used, whatever their names. "
                        "Repeatable. NAME= gives the target of a file that has none inside "
                        "(written before 0.1.1, or by hand). Default: --out-dir.")
    p.add_argument("--target", default=None,
                   help="Stratify only this target. If no score file names it but exactly "
                        "one has no target inside, that file is taken to be this target.")
    p.add_argument("--score-bin", action="append", dest="score_bins",
                   metavar="LO:HI:N", required=True,
                   help="Sample up to N accessions whose kmindex score falls in the "
                        "half-open interval [LO, HI), or in [LO, 1.0] when HI is 1.0. "
                        "Repeatable to define multiple bins in one run (e.g. --score-bin "
                        "0.0:0.05:10 --score-bin 0.05:0.10:10). Nonzero scores only; a bin "
                        "short on candidates contributes what it has (no backfill); "
                        "overlapping bins are refused.")
    p.set_defaults(func=run)


def run(args) -> int:
    try:
        parsed_bins = [parse_bin_spec(b) for b in args.score_bins]
        assert_bins_disjoint(parsed_bins)
    except ValueError as e:
        print(f"stratify: {e}", file=sys.stderr)
        return 2

    out_dir = Path(args.out_dir)
    try:
        score_files = find_score_files(args.scores or [out_dir], args.target)
    except (OSError, ValueError) as e:
        # Finding nothing used to print "[!] No *.scores.json files found" and exit 0.
        print(f"stratify: {e}", file=sys.stderr)
        return 1
    stratify_bins(out_dir, parsed_bins, score_files)
    return 0
