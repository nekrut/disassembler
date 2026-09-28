"""harvest: pull kmindex planetary shard collections from a Galaxy history.

For each shard dataset in the collection:
  1. Fetches its JSON body via the history content endpoint.
  2. Merges containment scores per target across all shards.
  3. Emits `<target>.accessions.txt` (sorted by score descending) and
     `<target>.scores.json` per target, plus `UNION.accessions.txt` and
     `summary.json`.

Stratification into a size-capped cohort is a separate subcommand (`stratify`);
this one only harvests.

The phiX174 Illumina spike-in is excluded by TARGET NAME in `safe_target_name`,
not by score.
"""

from __future__ import annotations

import argparse
import collections
import concurrent.futures
import json
import os
import sys
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter

from disassembler.validators import positive_int, validate_api_key_file

ATTEMPTS = 5


def get_galaxy_api_key(key_file: str | None = None) -> str:
    if key_file and Path(key_file).exists():
        return Path(key_file).read_text().strip()
    if os.environ.get("GALAXY_API_KEY"):
        return os.environ["GALAXY_API_KEY"].strip()
    default_tmp = Path("/tmp/gxy.txt")
    if default_tmp.exists():
        return default_tmp.read_text().strip()
    raise RuntimeError(
        "Galaxy API key not found. Provide --api-key-file, set GALAXY_API_KEY, "
        "or write /tmp/gxy.txt."
    )


def safe_target_name(raw):
    """First token of a kmindex target name, reduced to something usable as both a
    dict key and a filename. Returns None when nothing usable remains.

    The value comes from the shard JSON, so it is untrusted. `raw.split()[0]` raised
    IndexError on an empty or whitespace-only name, and the merge loop's only
    try/except wrapped `fut.result()`, so that escaped `as_completed`, tore down the
    executor, and killed a multi-hour harvest with nothing written. It is also
    interpolated straight into `out_dir / f"{target}.accessions.txt"`, so a name
    containing a path separator wrote outside `out_dir` and one containing `..`
    escaped it.

    The empty case is reachable without a malformed shard: the phiX174 strip below
    turns a target named exactly `phiX174` into an empty string, which then wrote
    dotfiles.
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


RETRY_AFTER_CAP_SECS = 120


def fetch_shard(item, base_url, sess):
    """Fetch one shard's JSON. Returns None when every attempt failed.

    None, not {}: an exhausted retry loop used to be indistinguishable from a shard
    that legitimately matched nothing, and nothing counted either. An expired API key
    therefore produced `Harvest complete: 0 unique SRA runs across 0 targets`, a
    UNION file holding a single newline, and exit 0.

    HTTP 429 is special-cased. Galaxy's dataset-display endpoint rate-limits per
    source IP and returns 429 with a `Retry-After: <seconds>` header (typically
    30-60 s on usegalaxy.org). The old policy treated 429 as a generic transient
    and slept the same 0.5/1.0/1.5/2.0 s ladder as any other retry, which was way
    too short: the retry loop just re-hit the rate limit five times in quick
    succession, all got 429, and the shard was reported failed. Now we honor
    Retry-After (capped at RETRY_AFTER_CAP_SECS so a hostile server can not
    freeze the worker) and consume an attempt, so ATTEMPTS still bounds worst-case
    per-shard runtime.
    """
    name, did = item
    url = f"{base_url}/api/datasets/{did}/display"
    last = "no attempt made"
    for attempt in range(ATTEMPTS):
        rate_limited = False
        retry_after = 0.0
        try:
            r = sess.get(url, params={"raw": "true"}, timeout=30)
            if r.status_code == 200 and r.text.strip():
                return json.loads(r.text)
            last = f"HTTP {r.status_code}"
            if r.status_code == 429:
                rate_limited = True
                try:
                    retry_after = float(r.headers.get("Retry-After", "30"))
                except ValueError:
                    retry_after = 30.0
                retry_after = min(max(retry_after, 1.0), RETRY_AFTER_CAP_SECS)
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
        # No backoff after the final attempt: it buys nothing and costs 2.5 s per
        # shard, which at 10,000 shards over 32 workers is ~40 minutes of sleeping
        # on a path that produces no data.
        if attempt < ATTEMPTS - 1:
            if rate_limited:
                print(f"    [~] shard {name}: HTTP 429, honoring Retry-After "
                      f"{retry_after:.0f}s (attempt {attempt + 1}/{ATTEMPTS})",
                      file=sys.stderr)
                time.sleep(retry_after)
            else:
                time.sleep(0.5 + attempt * 0.5)
    print(f"    [!] shard {name}: giving up after {ATTEMPTS} attempts ({last})",
          file=sys.stderr)
    return None


def harvest_collection(coll_id, history_id, base_url, api_key, out_dir, max_workers=32):
    """Walk a dataset collection via the history endpoint; fetch all shard JSONs in parallel."""
    print(f"[*] Fetching collection elements for {coll_id} in history {history_id}...")
    headers = {"x-api-key": api_key}
    url = f"{base_url}/api/histories/{history_id}/contents/dataset_collections/{coll_id}"
    # Wrap the initial GET: `requests.get` can raise ConnectionError / Timeout / etc. on a
    # bad --base-url or a dead network, which used to surface as a raw Python traceback
    # instead of the intended one-line `harvest: ...` error.
    try:
        r = requests.get(url, headers=headers, timeout=120)
    except requests.exceptions.RequestException as e:
        raise RuntimeError(
            f"could not reach Galaxy at {base_url}: {type(e).__name__}: {e}"
        ) from e
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
    adapter = HTTPAdapter(
        pool_connections=max_workers * 2,
        pool_maxsize=max_workers * 2,
        max_retries=3,
    )
    sess.mount("https://", adapter)
    sess.mount("http://", adapter)
    sess.headers.update(headers)

    per_target = collections.defaultdict(dict)
    n = 0
    n_failed = 0
    n_empty = 0
    n_bad_target = 0
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[*] Streaming and merging shard JSONs across {max_workers} threads...",
          flush=True)
    start_t = time.time()
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
            futures = {ex.submit(fetch_shard, it, base_url, sess): it for it in items}
            for fut in concurrent.futures.as_completed(futures):
                n += 1
                if n % 250 == 0 or n == len(items):
                    elapsed = time.time() - start_t
                    rate = n / elapsed if elapsed > 0 else 0
                    print(f"    Merged {n}/{len(items)} shards in {elapsed:.1f}s "
                          f"({rate:.1f} shards/s)...", flush=True)
                try:
                    d = fut.result()
                except Exception as e:
                    n_failed += 1
                    print(f"    [!] shard fetch raised {type(e).__name__}: {e}",
                          file=sys.stderr)
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
                # One malformed shard must not tear down the executor and lose the
                # whole run.
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
                    print(f"    [!] shard payload unusable "
                          f"({type(e).__name__}: {e})", file=sys.stderr)
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
            "max_score": max(acc_dict.values()) if acc_dict else 0.0,
        }
        print(f"  [✓] Target {target:12}: {len(accs):7d} candidates "
              f"(max containment: {summary[target]['max_score']:.2f})")

    union_file = out_dir / "UNION.accessions.txt"
    union_file.write_text("\n".join(sorted(all_union)) + "\n")
    summary["_union"] = {"n": len(all_union)}
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    if n_failed:
        print(f"[!] {n_failed}/{len(items)} shards failed to fetch after "
              f"{ATTEMPTS} attempts each; the cohort below is INCOMPLETE.",
              file=sys.stderr)
    if n_empty:
        print(f"[*] {n_empty}/{len(items)} shards were empty (no hits), which is normal.")
    if n_bad_target:
        print(f"[!] {n_bad_target} target name(s) were unusable as a key or filename "
              f"and were skipped.", file=sys.stderr)

    print(f"\n[{'!' if n_failed else '✓'}] Harvest complete: {len(all_union)} unique "
          f"SRA runs across {len(per_target)} targets"
          f"{f', {n_failed} shard(s) FAILED' if n_failed else ''}.")

    # A run that fetched nothing has produced no cohort. Saying so is the
    # difference between a downstream stage refusing to start and one silently
    # analysing an empty accession list.
    if items and n_failed == len(items):
        raise RuntimeError(
            f"all {len(items)} shard fetches failed (last error above); "
            f"no cohort was harvested. Check the Galaxy API key and the collection id."
        )
    return per_target


def add_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "harvest",
        help="Pull kmindex shard collections from a Galaxy history.",
        description=(
            "Streams a Galaxy kmindex planetary dataset collection and writes per-target "
            "accession lists and score files. Cohort cutting is a separate step; see "
            "`disassembler stratify`."
        ),
    )
    p.add_argument("--history-id", required=True,
                   help="Galaxy history ID containing the collection")
    p.add_argument("--collection-id", required=True,
                   help="Galaxy dataset collection ID for kmindex query shards")
    p.add_argument("--base-url", default="https://usegalaxy.org",
                   help="Galaxy instance URL (default: https://usegalaxy.org)")
    p.add_argument("--api-key-file",
                   help="Path to Galaxy API key file (also honours $GALAXY_API_KEY, "
                        "then /tmp/gxy.txt)")
    p.add_argument("-o", "--out-dir", default="candidates",
                   help="Output directory for candidate accessions")
    p.add_argument("--download-workers", type=positive_int, default=6,
                   help="Parallel download workers (default: 6). Empirically the ceiling "
                        "before usegalaxy.org's per-source rate limiter starts returning "
                        "HTTP 429; higher values still complete thanks to the Retry-After "
                        "handshake but pay 30-60 s per throttled shard. Raise it for a "
                        "self-hosted Galaxy without a rate limiter. Named --download-workers "
                        "to keep -w/--workers reserved for the walker's compute workers, "
                        "which are a different resource.")
    p.set_defaults(func=run)


def run(args) -> int:
    key_err = validate_api_key_file(args.api_key_file)
    if key_err:
        print(f"harvest: {key_err}", file=sys.stderr)
        return 1
    try:
        api_key = get_galaxy_api_key(args.api_key_file)
    except RuntimeError as e:
        print(f"harvest: {e}", file=sys.stderr)
        return 2
    out_dir = Path(args.out_dir)
    try:
        harvest_collection(
            args.collection_id,
            args.history_id,
            args.base_url,
            api_key,
            out_dir,
            max_workers=args.download_workers,
        )
    except RuntimeError as e:
        print(f"harvest: {e}", file=sys.stderr)
        return 1
    return 0
