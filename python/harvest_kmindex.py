#!/usr/bin/env python3
"""harvest_kmindex.py: Harvest Galaxy kmindex planetary shard collections into candidate cohorts.

Given a Galaxy kmindex planetary job or output dataset collection:
  1. Walks all shard datasets in parallel via Galaxy history content endpoints.
  2. Merges containment scores across all shards for each query target.
  3. Emits candidate accession lists sorted by containment score, with no threshold of its own.
  4. Optionally cuts a size-capped cohort per target, stratified across four containment brackets.
  5. Optionally chains directly into `disassemble.py` to reconstruct haplotypes and phylogenies.

The phiX174 Illumina spike-in is excluded by TARGET NAME, in `safe_target_name`, not by score.
An earlier version of this docstring credited the score brackets with that filtering, which is
probably how their lower bound came to look deliberate; it was not. See `stratify_diversity`.

Examples:
  # Harvest a Galaxy collection across 32 threads
  python/harvest_kmindex.py \
      --history-id bbd44e69cb8906b5df0f53ea81dc39ae \
      --collection-id e7f5bcc8ff5aec6f \
      -o candidates/

  # Harvest and cut a 500-accession cohort per target
  python/harvest_kmindex.py \
      --history-id bbd44e69cb8906b5df0f53ea81dc39ae \
      --collection-id e7f5bcc8ff5aec6f \
      -o candidates/ --diversity 500

  # Re-cut cohorts from a harvest already on disk, without contacting Galaxy
  python/harvest_kmindex.py --stratify-only -o candidates/ --diversity 500

  # Sample fixed counts from arbitrary score bins (e.g. 10 accessions each from three low
  # non-zero coverage bins), for walker calibration or histogram-driven experiments
  python/harvest_kmindex.py --stratify-only -o candidates/ \
      --bin 0.0:0.05:10 --bin 0.05:0.10:10 --bin 0.10:0.15:10

Size --diversity to the harvest: it caps each target, so a harvest smaller than N comes back
whole and stratification only reorders it. On a real six-target, 8,817-accession panel,
--diversity 2500 capped exactly one target and ran 77% of everything; the low hundreds are where
it starts subsetting. See the README for the measured table.

--diversity and --bin are the two cohort-cutting modes and are mutually exclusive. --diversity
imposes the fixed 20/30/30/20 quotas over the four containment brackets, with backfill;
--bin lets the caller specify arbitrary [LO, HI) ranges and per-range counts with no backfill,
so a bin whose harvest has fewer than N candidates contributes what it has. --bin only
considers nonzero-score accessions.
"""

import os
import sys
import json
import time
import argparse
import collections
import concurrent.futures
from pathlib import Path
import requests
from requests.adapters import HTTPAdapter

def get_galaxy_api_key(key_file=None):
    if key_file and Path(key_file).exists():
        return Path(key_file).read_text().strip()
    if os.environ.get("GALAXY_API_KEY"):
        return os.environ["GALAXY_API_KEY"].strip()
    default_tmp = Path("/tmp/gxy.txt")
    if default_tmp.exists():
        return default_tmp.read_text().strip()
    raise RuntimeError("Galaxy API key not found. Provide --api-key-file, set GALAXY_API_KEY, or write /tmp/gxy.txt.")

ATTEMPTS = 5


def safe_target_name(raw):
    """First token of a kmindex target name, reduced to something usable as a dict key and a
    filename. Returns None when nothing usable remains.

    The value comes from the shard JSON, so it is untrusted. `raw.split()[0]` raised
    IndexError on an empty or whitespace-only name, and since the only try/except in the
    merge loop wrapped `fut.result()`, that escaped `as_completed`, tore down the executor and
    killed a multi-hour harvest with nothing written, because merging is in memory and files
    are written only at the end.

    It is also interpolated straight into `out_dir / f"{target}.accessions.txt"`, so a name
    containing a path separator wrote outside out_dir and one containing `..` escaped it.

    The empty case is reachable without a malformed shard: the phiX174 strip below turns a
    target named exactly `phiX174` into an empty string, which then wrote dotfiles.
    """
    tokens = str(raw).split()
    if not tokens:
        return None
    name = tokens[0].replace("phiX174_", "").replace("phiX174", "")
    name = name.replace("/", "_").replace("\\", "_")
    name = name.strip(". ")
    if not name or name in (".", ".."):
        return None
    return name


def fetch_shard(item, base_url, sess):
    """Fetch one shard's JSON. Returns None when every attempt failed.

    None, not {}: an exhausted retry loop used to be indistinguishable from a shard that
    legitimately matched nothing, and nothing counted either. An expired API key therefore
    produced `Harvest complete: 0 unique SRA runs across 0 targets`, a UNION file holding a
    single newline, and exit 0, which a downstream --diversity or `disassemble.py -l` then
    ran on as an empty cohort.
    """
    # No api_key parameter: authentication rides on `sess`, whose headers are set by
    # `harvest_collection`. The argument was accepted and never read.
    name, did = item
    url = f"{base_url}/api/datasets/{did}/display"
    last = "no attempt made"
    for attempt in range(ATTEMPTS):
        try:
            r = sess.get(url, params={"raw": "true"}, timeout=30)
            if r.status_code == 200 and r.text.strip():
                return json.loads(r.text)
            last = f"HTTP {r.status_code}"
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
        # No backoff after the final attempt: it buys nothing and cost 2.5 s per shard,
        # which at 10,000 shards over 32 workers is ~40 minutes of pure sleeping on a path
        # that produces no data.
        if attempt < ATTEMPTS - 1:
            time.sleep(0.5 + attempt * 0.5)
    print(f"    [!] shard {name}: giving up after {ATTEMPTS} attempts ({last})", file=sys.stderr)
    return None

def harvest_collection(coll_id, history_id, base_url, api_key, out_dir, max_workers=32):
    """Walks dataset collection via history endpoint and fetches all shard JSONs in parallel."""
    print(f"[*] Fetching collection elements for {coll_id} in history {history_id}...")
    headers = {"x-api-key": api_key}
    url = f"{base_url}/api/histories/{history_id}/contents/dataset_collections/{coll_id}"
    r = requests.get(url, headers=headers, timeout=120)
    if r.status_code != 200:
        raise RuntimeError(f"Galaxy API error ({r.status_code}): {r.text}")
    dc = r.json()
    items = []

    def walk(elems):
        for e in elems:
            o = e.get("object", {})
            if "elements" in o and o.get("history_content_type") != "dataset":
                walk(o.get("elements", []))
            else:
                did = o.get("id")
                ident = e.get("element_identifier", did)
                if did:
                    items.append((ident, did))

    walk(dc.get("elements", []))
    print(f"[*] Total shard datasets in collection: {len(items)}")

    # Closed in the `finally` at the end of this function: the pooled connections outlive the
    # harvest otherwise, which matters at 32 workers against one host.
    sess = requests.Session()
    adapter = HTTPAdapter(pool_connections=max_workers * 2, pool_maxsize=max_workers * 2, max_retries=3)
    sess.mount("https://", adapter)
    sess.mount("http://", adapter)
    sess.headers.update(headers)

    per_target = collections.defaultdict(dict)
    n = 0
    n_failed = 0
    n_empty = 0
    n_bad_target = 0
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[*] Streaming and merging shard JSONs across {max_workers} threads...", flush=True)
    start_t = time.time()
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
            futures = {ex.submit(fetch_shard, it, base_url, sess): it for it in items}
            for fut in concurrent.futures.as_completed(futures):
                n += 1
                if n % 250 == 0 or n == len(items):
                    elapsed = time.time() - start_t
                    rate = n / elapsed if elapsed > 0 else 0
                    print(f"    Merged {n}/{len(items)} shards in {elapsed:.1f}s ({rate:.1f} shards/s)...", flush=True)
                try:
                    d = fut.result()
                except Exception as e:
                    n_failed += 1
                    print(f"    [!] shard fetch raised {type(e).__name__}: {e}", file=sys.stderr)
                    continue
                if d is None:
                    n_failed += 1
                    continue
                if not isinstance(d, dict):
                    n_failed += 1
                    continue
                if not d:
                    n_empty += 1
                    continue
                # One malformed shard must not tear down the executor and lose the whole run.
                try:
                    for shard, q in d.items():
                        if not isinstance(q, dict):
                            continue
                        for target_name, hits in q.items():
                            if not isinstance(hits, dict):
                                continue
                            clean_target = safe_target_name(target_name)
                            if clean_target is None:
                                n_bad_target += 1
                                continue
                            for acc, score in hits.items():
                                try:
                                    sc = float(score)
                                except Exception:
                                    continue
                                if sc > per_target[clean_target].get(acc, 0.0):
                                    per_target[clean_target][acc] = sc
                except Exception as e:
                    n_failed += 1
                    print(f"    [!] shard payload unusable ({type(e).__name__}: {e})", file=sys.stderr)
                    continue

    finally:
        # 32 pooled connections against one host outlive the harvest otherwise.
        sess.close()

    all_union = set()
    summary = {}
    for target, acc_dict in sorted(per_target.items()):
        accs = sorted(acc_dict, key=lambda a: -acc_dict[a])
        all_union.update(accs)
        acc_file = out_dir / f"{target}.accessions.txt"
        score_file = out_dir / f"{target}.scores.json"
        acc_file.write_text("\n".join(accs) + "\n")
        score_file.write_text(json.dumps(acc_dict, indent=2))
        summary[target] = {
            "n": len(accs),
            "max_score": max(acc_dict.values()) if acc_dict else 0.0
        }
        print(f"  [✓] Target {target:12}: {len(accs):7d} candidates (max containment: {summary[target]['max_score']:.2f})")

    union_file = out_dir / "UNION.accessions.txt"
    union_file.write_text("\n".join(sorted(all_union)) + "\n")
    summary["_union"] = {"n": len(all_union)}
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    if n_failed:
        print(f"[!] {n_failed}/{len(items)} shards failed to fetch after {ATTEMPTS} attempts "
              f"each; the cohort below is INCOMPLETE.", file=sys.stderr)
    if n_empty:
        print(f"[*] {n_empty}/{len(items)} shards were empty (no hits), which is normal.")
    if n_bad_target:
        print(f"[!] {n_bad_target} target name(s) were unusable as a key or filename and were "
              f"skipped.", file=sys.stderr)

    print(f"\n[{'!' if n_failed else '✓'}] Harvest complete: {len(all_union)} unique SRA runs "
          f"across {len(per_target)} targets"
          f"{f', {n_failed} shard(s) FAILED' if n_failed else ''}.")

    # A run that fetched nothing has produced no cohort. Saying so is the difference between
    # a downstream stage refusing to start and one silently analysing an empty accession list.
    if items and n_failed == len(items):
        raise RuntimeError(
            f"all {len(items)} shard fetches failed (last error above); "
            f"no cohort was harvested. Check the Galaxy API key and the collection id."
        )
    return per_target

def stratified_sample(items, k):
    """Evenly samples k items from a list sorted by score.

    A quota of zero is legitimate: `quota_canon` is `int(total_target * 0.20)`, which is 0 for
    any --diversity under 5, and `quota_high` is 0 for --diversity 1 or 2. The old guard was
    `len(items) <= k`, false whenever the bracket was non-empty, so the run reached
    `len(items) / k` and died with ZeroDivisionError. The bracket contributing nothing is the
    right answer at that size, and `stratify_diversity`'s top-up then fills the cohort from
    the pooled remainder.
    """
    if k <= 0:
        return []
    if len(items) <= k:
        return items
    step = len(items) / k
    return [items[int(i * step)] for i in range(k)]

# Share of each cohort drawn from each containment bracket. Preferences, not guarantees: a
# bracket holding fewer accessions than its quota contributes everything it has, and the
# shortfall is backfilled from whatever is left over.
QUOTA_CANON     = 0.20
QUOTA_HIGH      = 0.30
QUOTA_MODERATE  = 0.30
# divergent takes the remainder, so the four always sum to the requested size


def stratify_diversity(out_dir, total_target=2500):
    """Builds stratified diversity cohorts across four containment brackets.

    The brackets are canonical (>= 0.99), high homology (0.80-0.99), moderate (0.50-0.80) and
    divergent (everything below), and the last one is an `else`: it has no floor of its own and
    so inherits whatever threshold the kmindex search was run at.

    That matters. The search is launched at containment >= 0.30 specifically to reach divergent
    homologs (logan-hyphaeon, `work/launch_master_kmindex.py`), and the shard JSONs arrive
    already filtered, so nothing below it ever gets here. This function previously ended its
    chain at `elif sc >= 0.50` with no `else`, which meant an accession scoring 0.42 was filed
    nowhere at all: not rejected, not counted, just absent. Measured over six real planetary
    harvests that silently discarded 1,955 of 8,817 accessions (22%), and 189 of 230 for
    H5N1 HA (82%), whose median containment is 0.353. A second containment floor here is
    redundant with the search and was never deliberate; the `else` removes it.
    """
    print("\n" + "=" * 80)
    print(f" Generating Stratified Diversity Cohorts ({total_target:,} accessions / target)")
    print("=" * 80)

    quota_canon = int(total_target * QUOTA_CANON)
    quota_high = int(total_target * QUOTA_HIGH)
    quota_moderate = int(total_target * QUOTA_MODERATE)
    quota_div = total_target - quota_canon - quota_high - quota_moderate

    score_files = sorted(out_dir.glob("*.scores.json"))
    if not score_files:
        print(f"[!] No *.scores.json files found in {out_dir}.")
        return

    all_summaries = {}
    for sf in score_files:
        target = sf.stem.replace(".scores", "")
        with open(sf, "r") as f:
            scores = json.load(f)

        b_canon = []
        b_high = []
        b_moderate = []
        b_divergent = []

        for acc, sc in scores.items():
            if sc >= 0.99:
                b_canon.append((acc, sc))
            elif sc >= 0.80:
                b_high.append((acc, sc))
            elif sc >= 0.50:
                b_moderate.append((acc, sc))
            else:
                # No floor: the search already applied one. See the docstring.
                b_divergent.append((acc, sc))

        for b in (b_canon, b_high, b_moderate, b_divergent):
            b.sort(key=lambda x: -x[1])

        s_canon = stratified_sample(b_canon, quota_canon)
        s_high = stratified_sample(b_high, quota_high)
        s_moderate = stratified_sample(b_moderate, quota_moderate)
        s_div = stratified_sample(b_divergent, quota_div)

        selected = s_canon + s_high + s_moderate + s_div
        if len(selected) < total_target:
            used = {x[0] for x in selected}
            pool = [x for x in (b_high + b_canon + b_moderate + b_divergent)
                    if x[0] not in used]
            selected.extend(pool[:total_target - len(selected)])

        selected.sort(key=lambda x: -x[1])

        out_txt = out_dir / f"{target}.diversity_{total_target}.txt"
        out_json = out_dir / f"{target}.diversity_{total_target}.json"
        out_txt.write_text("\n".join(x[0] for x in selected) + "\n")
        out_json.write_text(json.dumps([{"acc": x[0], "score": x[1]} for x in selected], indent=2))

        smry = {
            "target": target,
            "total_selected": len(selected),
            "n_canonical_ge_99": sum(1 for x in selected if x[1] >= 0.99),
            "n_high_homology_80_99": sum(1 for x in selected if 0.80 <= x[1] < 0.99),
            # Renamed: the old key was `n_divergent_50_80`, which called this bracket
            # divergent while the genuinely divergent accessions were being dropped.
            "n_moderate_50_80": sum(1 for x in selected if 0.50 <= x[1] < 0.80),
            "n_divergent_below_50": sum(1 for x in selected if x[1] < 0.50),
            "n_available": len(scores),
            "min_score": min(x[1] for x in selected) if selected else 0.0,
            "max_score": max(x[1] for x in selected) if selected else 0.0,
        }
        all_summaries[target] = smry
        print(f"  [✓] Target {target:12}: {len(selected):5d} of {len(scores):,} "
              f"(Canon: {smry['n_canonical_ge_99']}, High: {smry['n_high_homology_80_99']}, "
              f"Mod: {smry['n_moderate_50_80']}, Div: {smry['n_divergent_below_50']}) "
              f"| score range: [{smry['min_score']:.2f}, {smry['max_score']:.2f}]")

    summary_file = out_dir / f"diversity_{total_target}_summary.json"
    summary_file.write_text(json.dumps(all_summaries, indent=2))
    print(f"[✓] Diversity summary saved: {summary_file}")


def parse_bin_spec(s):
    """Parse `LO:HI:N` into `(lo, hi, n)`. Raises ValueError with the offending spec on failure.

    Half-open `[lo, hi)` intervals match the existing bracket convention (`>= 0.99`, `< 0.99`).
    N is the count requested for the bin; N == 0 is legal and useful for "reserve a bin
    boundary but do not sample from it this run".
    """
    parts = s.split(":")
    if len(parts) != 3:
        raise ValueError(f"--bin '{s}' must be LO:HI:N (three colon-separated fields)")
    try:
        lo, hi, n = float(parts[0]), float(parts[1]), int(parts[2])
    except ValueError as e:
        raise ValueError(f"--bin '{s}': could not parse as LO:HI:N ({e})")
    if not (lo < hi):
        raise ValueError(f"--bin '{s}': LO ({lo}) must be strictly less than HI ({hi})")
    if n < 0:
        raise ValueError(f"--bin '{s}': N ({n}) must be non-negative")
    return (lo, hi, n)


def assert_bins_disjoint(bins):
    """Raise ValueError if any two half-open bins overlap. Pass --allow-overlap to skip.

    Sorted-neighbours is enough because [a,b) and [c,d) with c < d overlap iff
    c < b when sorted by lo.
    """
    ordered = sorted(bins, key=lambda b: b[0])
    for (a_lo, a_hi, _), (b_lo, b_hi, _) in zip(ordered, ordered[1:]):
        if b_lo < a_hi:
            raise ValueError(
                f"--bin ranges overlap: [{a_lo}, {a_hi}) and [{b_lo}, {b_hi}). "
                f"Pass --allow-overlap if that is intentional; the same accession will "
                f"then be able to appear in more than one bin's sample."
            )


def sample_bins_from_scores(scores, bins):
    """Sample per bin. Zero-score accessions are excluded up front.

    Returns `(picks, per_bin)` where `picks` is a merged list `(acc, score)` sorted by
    descending score, and `per_bin` is one dict per bin with counts and selections.

    No backfill: a bin with fewer than N candidates contributes what it has. This is the
    difference from `stratify_diversity`; here the caller is asking for a specific
    per-bin count, so quietly reshuffling would defeat the point.
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


def stratify_bins(out_dir, bins):
    """Per-target bin stratification. Emits `<target>.bins.txt` and `<target>.bins.json`
    per target and a combined `bins_summary.json`.
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
        all_summaries[target] = {"total_selected": len(picks), "bins": per_bin,
                                 "n_available_nonzero": sum(1 for s in scores.values() if s > 0.0)}

    summary_file = out_dir / "bins_summary.json"
    summary_file.write_text(json.dumps(all_summaries, indent=2))
    print(f"[✓] Bins summary saved: {summary_file}")


def main():
    p = argparse.ArgumentParser(description="Harvest Galaxy kmindex planetary shard collections into candidate cohorts.")
    p.add_argument("--history-id", help="Galaxy history ID containing the collection")
    p.add_argument("--collection-id", help="Galaxy dataset collection ID for kmindex query shards")
    p.add_argument("--base-url", default="https://usegalaxy.org", help="Galaxy instance URL (default: https://usegalaxy.org)")
    p.add_argument("--api-key-file", help="Path to Galaxy API key file")
    p.add_argument("-o", "--out-dir", default="candidates", help="Output directory for candidate accessions")
    p.add_argument("-w", "--workers", type=int, default=32, help="Parallel download workers (default: 32)")
    p.add_argument("--diversity", type=int, default=None,
                   help="Cut a cohort of N accessions PER TARGET, sampled across four containment "
                        "brackets (canonical >=0.99, high 0.80-0.99, moderate 0.50-0.80, divergent "
                        "below). N caps each target, so a smaller harvest is returned whole. Omit "
                        "to keep every accession. Mutually exclusive with --bin.")
    p.add_argument("--bin", action="append", dest="bins", metavar="LO:HI:N",
                   help="Sample N accessions with score in [LO, HI). Only nonzero-score "
                        "accessions are considered. Repeatable: pass one --bin per range, e.g. "
                        "--bin 0.0:0.05:10 --bin 0.05:0.10:10. Bins are strict: a bin with fewer "
                        "than N candidates contributes what it has and no other bin backfills it. "
                        "Mutually exclusive with --diversity.")
    p.add_argument("--allow-overlap", action="store_true",
                   help="Allow --bin ranges to overlap. Off by default because overlapping bins "
                        "usually mean a typo in the boundaries; on, the same accession can appear "
                        "in more than one bin's sample.")
    p.add_argument("--stratify-only", action="store_true", help="Skip Galaxy download; stratify existing *.scores.json")
    args = p.parse_args()

    if args.diversity and args.bins:
        p.error("--diversity and --bin are mutually exclusive; pick one.")

    parsed_bins = None
    if args.bins:
        try:
            parsed_bins = [parse_bin_spec(b) for b in args.bins]
            if not args.allow_overlap:
                assert_bins_disjoint(parsed_bins)
        except ValueError as e:
            p.error(str(e))

    out_dir = Path(args.out_dir)

    if not args.stratify_only:
        if not args.history_id or not args.collection_id:
            p.error("--history-id and --collection-id are required unless --stratify-only is specified.")
        api_key = get_galaxy_api_key(args.api_key_file)
        harvest_collection(args.collection_id, args.history_id, args.base_url, api_key, out_dir, max_workers=args.workers)

    if args.diversity:
        stratify_diversity(out_dir, total_target=args.diversity)
    if parsed_bins is not None:
        stratify_bins(out_dir, parsed_bins)

if __name__ == "__main__":
    main()
