"""stratify: cut a size-capped cohort from a harvested `<out-dir>/*.scores.json` tree.

`--score-bin LO:HI:N` (repeatable): sample N accessions with score in the
half-open `[LO, HI)`. Nonzero scores only. No backfill: a bin short on
candidates contributes what it has.

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
    `< 0.99`). N is the count requested for the bin; N == 0 is legal and useful
    for "reserve a bin boundary but do not sample from it this run".
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


def assert_bins_disjoint(bins) -> None:
    """Raise ValueError if any two half-open bins overlap.

    Overlap is a hard refusal in the current CLI (--allow-overlap was removed with
    --diversity). If overlapping bins are intentional, invoke stratify twice on
    disjoint subsets and concatenate the outputs.
    """
    ordered = sorted(bins, key=lambda b: b[0])
    for (a_lo, a_hi, _), (b_lo, b_hi, _) in zip(ordered, ordered[1:]):
        if b_lo < a_hi:
            raise ValueError(
                f"--score-bin ranges overlap: [{a_lo}, {a_hi}) and [{b_lo}, {b_hi}). "
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
            ((a, s) for a, s in nonzero.items() if lo <= s < hi),
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


def stratify_bins(out_dir: Path, bins) -> None:
    """Per-target bin stratification. Emits `<target>.bins.txt`, `<target>.bins.json`, and
    a combined `bins_summary.json`.
    """
    print("\n" + "=" * 80)
    ranges = ", ".join(f"[{lo}, {hi})x{n}" for lo, hi, n in bins)
    print(f" Bin stratification (nonzero scores only): {ranges}")
    print("=" * 80)

    score_files = sorted(out_dir.glob("*.scores.json"))
    if not score_files:
        print(f"[!] No *.scores.json files found in {out_dir}.")
        return

    all_summaries = {}
    for sf in score_files:
        target = sf.stem.replace(".scores", "")
        with open(sf, "r") as f:
            scores = json.load(f)

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
            f"[{b['lo']}, {b['hi']}): {b['selected']}/{b['requested']}" for b in per_bin
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
        help="Cut a size-capped cohort from a `harvest` output directory.",
        description=(
            "Sample N accessions per --score-bin from `<out-dir>/*.scores.json`. "
            "Bins are half-open [LO, HI), nonzero-scores only, no backfill, disjoint."
        ),
    )
    p.add_argument("-o", "--out-dir", default="candidates",
                   help="Directory with *.scores.json produced by `harvest` (the same "
                        "directory passed to `harvest -o`). Per-target outputs "
                        "<target>.bins.txt and <target>.bins.json are written here, "
                        "next to a combined bins_summary.json.")
    p.add_argument("--score-bin", action="append", dest="score_bins",
                   metavar="LO:HI:N", required=True,
                   help="Sample up to N accessions whose kmindex score falls in the "
                        "half-open interval [LO, HI). Repeatable to define multiple "
                        "bins in one run (e.g. --score-bin 0.0:0.05:10 "
                        "--score-bin 0.05:0.10:10). Nonzero scores only; a bin short "
                        "on candidates contributes what it has (no backfill); "
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
    stratify_bins(out_dir, parsed_bins)
    return 0
