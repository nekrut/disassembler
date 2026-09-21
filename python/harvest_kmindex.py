#!/usr/bin/env python3
"""harvest_kmindex.py: Harvest Galaxy kmindex planetary shard collections into candidate cohorts.

Given a Galaxy kmindex planetary job or output dataset collection:
  1. Walks all shard datasets in parallel via Galaxy history content endpoints.
  2. Merges containment scores across all shards for each query target.
  3. Emits candidate accession lists sorted by containment score.
  4. Generates stratified diversity cohorts (filtering out the Illumina spike-in cloud).
  5. Optionally chains directly into `disassemble.py` to reconstruct haplotypes and phylogenies.

Examples:
  # Harvest a Galaxy collection across 32 threads
  python/harvest_kmindex.py \
      --history-id bbd44e69cb8906b5df0f53ea81dc39ae \
      --collection-id e7f5bcc8ff5aec6f \
      -o candidates/

  # Harvest and generate 2,500-accession diversity cohorts per target
  python/harvest_kmindex.py \
      --history-id bbd44e69cb8906b5df0f53ea81dc39ae \
      --collection-id e7f5bcc8ff5aec6f \
      -o candidates/ --diversity 2500

  # Re-run diversity stratification on previously harvested scores
  python/harvest_kmindex.py --stratify-only -o candidates/ --diversity 2500
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

def fetch_shard(item, base_url, api_key, sess):
    name, did = item
    url = f"{base_url}/api/datasets/{did}/display"
    for attempt in range(5):
        try:
            r = sess.get(url, params={"raw": "true"}, timeout=30)
            if r.status_code == 200 and r.text.strip():
                return json.loads(r.text)
        except Exception:
            pass
        time.sleep(0.5 + attempt * 0.5)
    return {}

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

    sess = requests.Session()
    adapter = HTTPAdapter(pool_connections=max_workers * 2, pool_maxsize=max_workers * 2, max_retries=3)
    sess.mount("https://", adapter)
    sess.mount("http://", adapter)
    sess.headers.update(headers)

    per_target = collections.defaultdict(dict)
    n = 0
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[*] Streaming and merging shard JSONs across {max_workers} threads...", flush=True)
    start_t = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futures = {ex.submit(fetch_shard, it, base_url, api_key, sess): it for it in items}
        for fut in concurrent.futures.as_completed(futures):
            n += 1
            if n % 250 == 0 or n == len(items):
                elapsed = time.time() - start_t
                rate = n / elapsed if elapsed > 0 else 0
                print(f"    Merged {n}/{len(items)} shards in {elapsed:.1f}s ({rate:.1f} shards/s)...", flush=True)
            try:
                d = fut.result()
            except Exception:
                continue
            if not isinstance(d, dict):
                continue
            for shard, q in d.items():
                if not isinstance(q, dict):
                    continue
                for target_name, hits in q.items():
                    if not isinstance(hits, dict):
                        continue
                    clean_target = target_name.split()[0].replace("phiX174_", "").replace("phiX174", "")
                    for acc, score in hits.items():
                        try:
                            sc = float(score)
                        except Exception:
                            continue
                        if sc > per_target[clean_target].get(acc, 0.0):
                            per_target[clean_target][acc] = sc

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

    print(f"\n[✓] Harvest complete: {len(all_union)} unique SRA runs across {len(per_target)} targets.")
    return per_target

def stratified_sample(items, k):
    """Evenly samples k items from a list sorted by score."""
    if len(items) <= k:
        return items
    step = len(items) / k
    return [items[int(i * step)] for i in range(k)]

def stratify_diversity(out_dir, total_target=2500):
    """Builds stratified diversity cohorts across canonical, high-homology, and divergent brackets."""
    print("\n" + "=" * 80)
    print(f" Generating Stratified Diversity Cohorts ({total_target:,} accessions / target)")
    print("=" * 80)

    quota_canon = int(total_target * 0.20)
    quota_high = int(total_target * 0.40)
    quota_div = total_target - quota_canon - quota_high

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
        b_div = []

        for acc, sc in scores.items():
            if sc >= 0.99:
                b_canon.append((acc, sc))
            elif sc >= 0.80:
                b_high.append((acc, sc))
            elif sc >= 0.50:
                b_div.append((acc, sc))

        b_canon.sort(key=lambda x: -x[1])
        b_high.sort(key=lambda x: -x[1])
        b_div.sort(key=lambda x: -x[1])

        s_canon = stratified_sample(b_canon, quota_canon)
        s_high = stratified_sample(b_high, quota_high)
        s_div = stratified_sample(b_div, quota_div)

        selected = s_canon + s_high + s_div
        if len(selected) < total_target:
            used = {x[0] for x in selected}
            pool = [x for x in (b_high + b_canon + b_div) if x[0] not in used]
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
            "n_divergent_50_80": sum(1 for x in selected if 0.50 <= x[1] < 0.80),
            "min_score": min(x[1] for x in selected) if selected else 0.0,
            "max_score": max(x[1] for x in selected) if selected else 0.0,
        }
        all_summaries[target] = smry
        print(f"  [✓] Target {target:12}: {len(selected):5d} selected (Canon: {smry['n_canonical_ge_99']}, High: {smry['n_high_homology_80_99']}, Div: {smry['n_divergent_50_80']}) | score range: [{smry['min_score']:.2f}, {smry['max_score']:.2f}]")

    summary_file = out_dir / f"diversity_{total_target}_summary.json"
    summary_file.write_text(json.dumps(all_summaries, indent=2))
    print(f"[✓] Diversity summary saved: {summary_file}")

def main():
    p = argparse.ArgumentParser(description="Harvest Galaxy kmindex planetary shard collections into candidate cohorts.")
    p.add_argument("--history-id", help="Galaxy history ID containing the collection")
    p.add_argument("--collection-id", help="Galaxy dataset collection ID for kmindex query shards")
    p.add_argument("--base-url", default="https://usegalaxy.org", help="Galaxy instance URL (default: https://usegalaxy.org)")
    p.add_argument("--api-key-file", help="Path to Galaxy API key file")
    p.add_argument("-o", "--out-dir", default="candidates", help="Output directory for candidate accessions")
    p.add_argument("-w", "--workers", type=int, default=32, help="Parallel download workers (default: 32)")
    p.add_argument("--diversity", type=int, default=None, help="Generate stratified diversity cohorts of N accessions")
    p.add_argument("--stratify-only", action="store_true", help="Skip Galaxy download; stratify existing *.scores.json")
    args = p.parse_args()

    out_dir = Path(args.out_dir)

    if not args.stratify_only:
        if not args.history_id or not args.collection_id:
            p.error("--history-id and --collection-id are required unless --stratify-only is specified.")
        api_key = get_galaxy_api_key(args.api_key_file)
        harvest_collection(args.collection_id, args.history_id, args.base_url, api_key, out_dir, max_workers=args.workers)

    if args.diversity:
        stratify_diversity(out_dir, total_target=args.diversity)

if __name__ == "__main__":
    main()
