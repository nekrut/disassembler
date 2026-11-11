# Python: kmindex harvest and stratify (cmd/harvest.py, cmd/stratify.py)

`harvest` turns kmindex shard JSON (from a Galaxy collection or a local directory) into per-target accession lists and score files; `stratify` cuts a size-capped cohort from those score files by score bins. Much of the error handling is already careful: an exhausted shard fetch is told apart from an empty one, finding no score file is an error, target names cannot escape the output directory, and the target is now stored inside each score file. The weak spot is silence. Cohorts change between identical runs because tied scores keep network arrival order, and a harvest that lost shards, read wrongly shaped JSON, or met an empty collection still exits 0 and leaves no record of the loss. `stratify` passes damaged or stale files over without saying so, and cannot read 0.1.0 harvest directories. The Galaxy client has retry, redirect and API-key handling gaps. Most fixes share one theme: order deterministically, turn every dropped input into a counted, named failure with a non-zero exit, and route the duplicated per-shard loop and the two copies of the API-key lookup through one helper each, so each fix is made once.

| ID | Severity | Where | Title | Fix effort |
|----|----------|-------|-------|------------|
| PH-01 | high | `src/disassembler/cmd/harvest.py:180` | Tied scores keep shard arrival order, so the same collection yields different lists and cohorts | trivial |
| PH-02 | medium | `src/disassembler/cmd/harvest.py:279` | A partial harvest (some shards lost) exits 0 and summary.json does not record the loss | small |
| PH-03 | medium | `src/disassembler/cmd/stratify.py:229` | 0.1.0 harvest directories are unusable by 0.1.1 stratify; `from-kmindex --skip-harvest` dead-ends or relabels the wrong target | small |
| PH-04 | medium | `src/disassembler/cmd/harvest.py:151` | Wrongly shaped shard JSON counts as a clean shard (0 targets, exit 0) | small |
| PH-05 | medium | `src/disassembler/cmd/harvest.py:179` | Re-harvesting into the same directory keeps old targets, and stratify and from-kmindex use them | small |
| PH-06 | medium | `src/disassembler/cmd/stratify.py:196` | stratify's directory scan silently drops damaged score files; a corrupt gzip escapes as a traceback | small |
| PH-07 | medium | `src/disassembler/cmd/harvest.py:75` | Targets are renamed by a phiX174 substring strip, and name collisions are merged silently | small |
| PH-08 | medium | `src/disassembler/cmd/harvest.py:385` | A Galaxy collection with zero elements harvests to an empty cohort with exit 0 | trivial |
| PH-09 | medium | `src/disassembler/cmd/harvest.py:47` | Docs steer the Galaxy API key into world-readable /tmp/gxy.txt; `--api-key` exposes it on argv | small |
| PH-10 | medium | `src/disassembler/cmd/harvest.py:324` | urllib3 retries 429/503 with Retry-After inside each attempt: 4x requests, 120 s cap bypassed | trivial |
| PH-11 | medium | `src/disassembler/cmd/stratify.py:195` | stratify parses and keeps every target's score map in memory, even with `--target` | small |
| PH-12 | low | `src/disassembler/textio.py:32` | A UTF-8 BOM makes a score file unreadable, or silently drops it from a directory | trivial |
| PH-13 | low | `src/disassembler/cmd/stratify.py:58` | stratified_sample is top-anchored: N=1 takes the bin maximum, the bin bottom is never drawn | trivial |
| PH-14 | low | `src/disassembler/cmd/harvest.py:328` | The x-api-key header is re-sent on cross-host redirects; an http:// base URL is accepted silently | small |
| PH-15 | low | `src/disassembler/cmd/harvest.py:44` | A Galaxy API key with an embedded newline is echoed in full in the error message | small |
| PH-16 | low | `src/disassembler/cmd/stratify.py:76` | stratify accepts out-of-range and non-finite bins and exits 0 having selected nothing | small |
| PH-17 | low | `src/disassembler/cmd/harvest.py:113` | Non-JSON 200 shards are retried 5 times, then the API key is blamed | small |
| PH-18 | low | `src/disassembler/cmd/harvest.py:303` | A trailing-slash `--base-url` or a non-JSON 200 listing ends in a traceback | trivial |
| PH-19 | low | `src/disassembler/cmd/harvest.py:76` | Shard names, accessions and scores are trusted: a NUL or overlong target name aborts with a traceback | small |
| PH-20 | low | `src/disassembler/cmd/stratify.py:238` | Following stratify's own advice fails: a directory plus NAME=FILE for a file inside it conflict | trivial |
| PH-21 | low | `src/disassembler/cmd/stratify.py:293` | Per-target stratify runs overwrite the "combined" bins_summary.json with the last target only | trivial |
| PH-22 | low | `src/disassembler/cmd/harvest.py:341` | Ctrl-C during a Galaxy harvest keeps fetching every queued shard before exiting | small |
| PH-23 | low | `src/disassembler/cmd/harvest.py:45` | A whitespace-only $GALAXY_API_KEY or blank /tmp/gxy.txt is sent as an empty key | small |
| PH-24 | low | `src/disassembler/cmd/harvest.py:360` | Galaxy-mode harvest does not name the shard it discards; the per-shard loop is duplicated | small |
| PH-25 | low | `src/disassembler/cmd/stratify.py:202` | `--scores` rejects pipes as "no such file"; `--shards-dir FILE` says "not found" | trivial |
| PH-26 | low | `src/disassembler/cmd/stratify.py:46` | stratified_sample's docstring, its test and cli.py describe removed code | trivial |
| PH-27 | low | `src/disassembler/cmd/harvest.py:342` | harvest keeps every parsed shard alive until the end and builds each output file as one string | small |

### PH-01 · high · Tied scores keep shard arrival order, so the same collection yields different lists and cohorts

`src/disassembler/cmd/harvest.py:180` · provenance · confirmed (reproduced) · fix effort: trivial · sources: P3-01, PC-02, SC2-04, SWD-02

**Problem.** `harvest.py:180` sorts on score alone: `accs = sorted(acc_dict, key=lambda a: -acc_dict[a])`. Python's sort is stable, so accessions with equal scores keep dict insertion order, which is shard arrival order: `as_completed(futures)` in Galaxy API mode (`harvest.py:343`), file-name order in `--shards-dir` mode (`harvest.py:233`). `harvest.py:186` writes `acc_dict` in the same order, `stratify.py:125` again sorts on `-x[1]` only, and `stratified_sample` picks by rank (`stratify.py:58`). Which tied accessions enter a cohort is therefore decided by network timing. kmindex scores tie in bulk: in the SARS-CoV-2 spike harvest 274,808 accessions score exactly 1.0 (`docs/kmindex.md:59-61`).

**Impact.** Any API-mode harvest with more than one download worker (default 6, the documented route) gives a different `<target>.accessions.txt`, `<target>.scores.json`, `<target>.bins.txt` and cohort on each run. `from-kmindex --score-bin` re-harvests every time, so a rerun that changes one walker option also changes the cohort, silently. API mode and `--shards-dir` mode disagree as well, although `docs/kmindex.md:12` says they "produce identical output". Without `--score-bin` the list order still reaches `--sieve`, which is order-greedy. Every pick is inside its bin and each run records its cohort in `<target>.bins.json`; the rating stays high (a second verifier argued medium) because a rerun of the default route gives no sign that the cohort changed, which the rubric lists as a reproducibility hole users cannot detect.

**Evidence.** Against a stand-in Galaxy server (40 shards, 10,000 accessions, 5.3% at 1.0, 0-30 ms jitter), two `harvest` runs with default workers followed by `stratify --score-bin 0.9:1.0:20`: all four output files differed, the `{accession: score}` maps were equal, and the two 20-accession cohorts shared 14 accessions (8 in an independent rerun). The same shards named `0.json…` versus `dataset_<uuid>.dat` in `--shards-dir` mode gave cohorts sharing 3 of 20. `--download-workers 1` was byte-stable.

**Proposed fix.** Break ties by accession wherever a list is ordered:
```diff
 # harvest.py, _write_outputs
-        accs = sorted(acc_dict, key=lambda a: -acc_dict[a])
+        accs = sorted(acc_dict, key=lambda a: (-acc_dict[a], a))
-        score_file.write_text(json.dumps({"target": target, "scores": acc_dict}, indent=2))
+        score_file.write_text(json.dumps(
+            {"target": target, "scores": {a: acc_dict[a] for a in accs}}, indent=2))
 # stratify.py, sample_bins_from_scores (both sorts)
-            key=lambda x: -x[1],
+            key=lambda x: (-x[1], x[0]),
-    picks.sort(key=lambda x: -x[1])
+    picks.sort(key=lambda x: (-x[1], x[0]))
```
The stratify change alone fixes the cohort; the harvest change makes the lists byte-stable and covers the no-bin path. At 1-2 M accessions the tuple key costs about 2 s and 100 MB extra; `accs = sorted(acc_dict); accs.sort(key=lambda a: -acc_dict[a])` gives the same order without the extra memory. With the fix every output was identical across runs and across modes.

**Tests and docs.** Add byte-comparison tests: `harvest_collection` with a stubbed `fetch_shard` completing in two orders, `harvest_local` under two sets of file names, `sample_bins_from_scores` on a reversed score map. Write "by score, highest first, ties by accession" in `docs/kmindex.md` (output table, stratify section) and `harvest.py:11`. Release note: 0.1.1 cohorts are not regenerated by the new order. tools-iuc PR #8800's `test-data/hiv1_pol.scores.json` holds arrival order and needs regenerating.

### PH-02 · medium · A partial harvest (some shards lost) exits 0 and summary.json does not record the loss

`src/disassembler/cmd/harvest.py:279` · reliability · confirmed (reproduced) · fix effort: small · sources: P3-04, D2-04, PA-07

**Problem.** Both modes fail only when every shard failed: `if shard_files and n_failed == len(shard_files):` (`harvest.py:279`) and `if items and n_failed == len(items):` (`harvest.py:385`). A smaller loss prints one stderr line ("the cohort below is INCOMPLETE", `harvest.py:200`) and `run()` returns 0 (`harvest.py:478`, `:502`). `summary.json` holds only per-target counts and `_union`, so nothing on disk marks the cohort as partial. `harvest_local` also keeps only `p.is_file()` entries (`harvest.py:233-234`), so a dangling symlink is not even counted.

**Impact.** Shards that exhaust their retries (429s, a key expiring mid-run, timeouts) or saved HTML error pages in `--shards-dir` give exit 0; `from-kmindex` goes on to stratify and reconstruct (`pipe_from_kmindex.py:318-320`) a cohort missing every accession of the lost shards, and the outputs cannot show it afterwards.

**Evidence.** A shards directory with one good shard, one HTML page, one empty file and one dangling symlink: `disassembler harvest --shards-dir sh -o out` printed "Found 3 shard JSON file(s)" and "2/3 shard files failed", exit 0; `from-kmindex` on the same directory ran all three stages.

**Proposed fix.** Record the loss and make it an error by default:
```python
summary["_harvest"] = {"n_shards": n_items, "n_failed": n_failed,
                       "n_empty": n_empty, "complete": n_failed == 0}
...
class PartialHarvest(RuntimeError): ...
if n_failed and not allow_partial:          # after outputs are written
    raise PartialHarvest(f"{n_failed}/{n_items} shards failed; cohort is INCOMPLETE")
```
Add `--allow-partial` to `harvest` and `from-kmindex` (forwarded via `stage_args`); `run()` maps `PartialHarvest` to a distinct documented exit status (for example 3, but note the walker already uses 3 for "an accession errored"). Count dangling links: select `p.is_file() or (p.is_symlink() and not p.exists())`. With this the example above exits non-zero and `summary.json` carries `"complete": false`.

**Tests and docs.** New test: good shard, HTML shard, dangling symlink; assert the exit status, `_harvest.n_failed == 2` and outputs present, and exit 0 with `--allow-partial`. The `fake_collection` stub in `tests/test_composites.py:149` must accept `allow_partial`. Document the exit status, the flag and `_harvest` in `docs/kmindex.md` and in the tools-iuc harvest wrapper help. This loop is shared with PH-04 and PH-24; doing PH-24's refactor first means one change.

### PH-03 · medium · 0.1.0 harvest directories are unusable by 0.1.1 stratify; `from-kmindex --skip-harvest` dead-ends or relabels the wrong target

`src/disassembler/cmd/stratify.py:229` · correctness · confirmed (reproduced) · fix effort: small · sources: PB-01

**Problem.** 0.1.0 `harvest` wrote bare `{accession: score}` maps to `<target>.scores.json`, so the file name was the only record of the target. 0.1.1 `stratify` takes the target only from inside the file. With `--target T` and no file naming T, `elif len(unnamed) == 1: found = [unnamed[0]._replace(target=target)]` (`stratify.py:229-230`) relabels the lone legacy file as T whatever its name says; with two or more it stops with "name the right one with --scores T=FILE" (`stratify.py:231-234`). `from-kmindex` forwards only `out_dir`, `target` and `score_bins` (`pipe_from_kmindex.py:324-329`) and has no `--scores` option, so that advice cannot be followed there.

**Impact.** Reusing a 0.1.0 cohort directory with `from-kmindex --skip-harvest --target hiv1_pol --score-bin …` (the documented reuse path, which `docs/kmindex.md:122-124` says works): with `hiv1_pol.scores.json` and `pfcrt.scores.json` present it exits 1 pointing at a flag the command lacks; with only `pfcrt.scores.json` present it writes `hiv1_pol.bins.txt` holding pfcrt's accessions and starts the walker on them with the hiv1_pol query (0.1.0 stopped with "expected accession list not found"). Plain `stratify -o <0.1.0 dir>` fails with "no target inside", which also breaks the tools-iuc PR #8800 wrapper's own test data.

**Evidence.** With PR #8800's bare-shape test files: two files gave `stratify: no score file names target 'hiv1_pol', and 2 have no target inside; name the right one with --scores hiv1_pol=FILE`, exit 1; one file (pfcrt) gave `hiv1_pol.bins.txt` = pfcrt's accessions. `disassembler from-kmindex -h | grep -c -- --scores` prints 0.

**Proposed fix.** In `find_score_files`' directory branch, when a file has no target inside and is named `<name>.scores.json` (the only name 0.1.0 harvest wrote), take `<name>` as its target and print `[*] <path>: taking '<name>' from the file name`. Keep `NAME=FILE` for every other bare file, and make the lone-unnamed adoption at line 229 print which file it took. Correct the advice for `from-kmindex` ("re-run harvest, or run `disassembler stratify --scores T=FILE` first"), or add a `--scores` pass-through. With the change the two-file case stratifies hiv1_pol and the one-file case ends with "no score file for target 'hiv1_pol'; the files found are for pfcrt".

**Tests and docs.** `tests/test_stratify.py:198-209` (`test_a_pre_0_1_1_file_needs_its_target_named`) pins the old behaviour; change it to expect `ns5.scores.json` to stratify as ns5 and `--target gyrA` to fail. Add a `tests/test_composites.py` case with two legacy files under `--skip-harvest`. Update `stratify.py:8-14`, the `--scores`/`--target` help and `docs/kmindex.md:66-73`, `:122-124`. Related: PH-20 (the NAME=FILE advice also fails inside `stratify`).

### PH-04 · medium · Wrongly shaped shard JSON counts as a clean shard (0 targets, exit 0)

`src/disassembler/cmd/harvest.py:151` · correctness · confirmed (reproduced) · fix effort: small · sources: P3-02

**Problem.** `_merge_payload` expects kmindex `--format json` (`{index: {query: {sample: number}}}`) but skips anything else without counting it: `if not isinstance(q, dict): continue` (`harvest.py:151-152`), `if not isinstance(hits, dict): continue` (`:154-155`), and `except Exception: continue` around `float(score)` (`:161-164`). The callers count only unparseable files, non-dict top levels and exceptions, so a well-formed object of the wrong shape is a successful shard that contributed nothing. The expected shape is not documented.

**Impact.** kmindex `--format json_vec` output (`{"R": ratio, "P": [...]}` per sample, offered by the Galaxy kmindex tool), a document one level short, or a saved Galaxy error body (`{"err_msg": …}`) gives `Harvest complete: 0 unique SRA runs across 0 targets`, exit 0, and a green Galaxy job; mixed with good shards, their accessions vanish without a warning. `from-kmindex` then fails later with the misleading "expected accession list not found".

**Evidence.** Each of `{"shard1":{"gene":{"SRR1":{"R":0.97,"P":[1,1,0]}}}}`, `{"gene":{"SRR1":0.97}}` and `{"err_msg":"History is not accessible","err_code":403002}` as the only file in `--shards-dir` gave 0 targets and exit 0.

**Proposed fix.** Validate the whole payload before merging and raise; the callers already count a raised exception as a failed, named shard:
```python
class ShardShapeError(ValueError): ...
def _merge_payload(d, per_target):
    staged = []
    for shard, q in d.items():
        if not isinstance(q, dict):
            raise ShardShapeError(f"index {shard!r} maps to a {type(q).__name__}")
        for target_name, hits in q.items():
            if not isinstance(hits, dict):
                raise ShardShapeError(f"query {target_name!r} maps to a {type(hits).__name__}")
            for acc, score in hits.items():
                if isinstance(score, bool) or not isinstance(score, (int, float)):
                    hint = (" (was kmindex run with --format json_vec?)"
                            if isinstance(score, dict) and "R" in score else "")
                    raise ShardShapeError(f"sample {acc!r} has a {type(score).__name__} "
                                          f"score, not a number{hint}")
            staged.append((target_name, hits))
    # then the existing merge over `staged`
```
With this, the three cases fail with the shard named and exit 1. Numeric strings (`"0.9"`) are now rejected; kmindex never writes them.

**Tests and docs.** Add `TestShardShape` to `tests/test_harvest.py` (json_vec with the hint, one level short, error body). State the accepted shape and `--format json` in `docs/kmindex.md` and the tools-iuc harvest help. Related: PH-02, PH-24.

### PH-05 · medium · Re-harvesting into the same directory keeps old targets, and stratify and from-kmindex use them

`src/disassembler/cmd/harvest.py:179` · correctness · confirmed (reproduced) · fix effort: small · sources: P3-12

**Problem.** `_write_outputs` writes `<target>.accessions.txt` and `<target>.scores.json` only for the current run's targets (`harvest.py:179`) and rewrites `UNION.accessions.txt` and `summary.json`, but never removes an earlier harvest's per-target files in the same `-o`. `stratify` in directory mode reads every score file it finds (`stratify.py:192-201`), and `from-kmindex` picks `cohort_dir / f"{args.target}.accessions.txt"` (or `.bins.txt`) checking only `is_file()` (`pipe_from_kmindex.py:338-345`).

**Impact.** Re-running `from-kmindex --target X -o OUT` against a new collection in which X has no hits reconstructs the old cohort of X with no warning; with `--score-bin` it cuts bins from the stale `X.scores.json`. Plain `stratify` after two harvests writes cohorts for targets that `summary.json` does not list.

**Evidence.** `harvest --shards-dir run1 -o cand` (target oldgene) then `harvest --shards-dir run2 -o cand` (newgene only): `cand/` still held `oldgene.*`, `stratify -o cand` selected both targets, and `from-kmindex --shards-dir run2 --target oldgene` printed "target 'oldgene': 2 accessions" and ran reconstruct.

**Proposed fix.**
1. `harvest.py`: at the top of `_write_outputs`, call `_drop_stale_targets(out_dir, set(per_target))`, which reads the previous `summary.json` and, for each listed target (not `_union`) absent now and with `safe_target_name(t) == t`, unlinks `t.{accessions.txt,scores.json,bins.txt,bins.json}` and says so on stderr. Only files named by the earlier summary are touched.
2. `pipe_from_kmindex.py`: after a successful harvest, if `args.target` is not a key of `cohort_dir/summary.json`, fail with `the harvest has no hits for target 'X' (targets found: …)`.
3. Optional: in `stratify` directory mode, warn about score files whose target `summary.json` does not list.

**Tests and docs.** Test two harvests with different targets into one directory (old files gone, `find_score_files` returns only the new target) and `from-kmindex` with a missing target returning 1. State in `docs/kmindex.md` (output table) that re-harvesting replaces the earlier per-target files; a user who combined two collections in one directory loses the first set, so mention it or offer `--keep-existing`. Related: PH-21.

### PH-06 · medium · stratify's directory scan silently drops damaged score files; a corrupt gzip escapes as a traceback

`src/disassembler/cmd/stratify.py:196` · correctness · confirmed (reproduced) · fix effort: small · sources: P3-08

**Problem.** In directory mode (the default `-o DIR` without `--scores`, and how the tools-iuc wrapper runs it) every read error means "not a score file": `except (OSError, ValueError): sf = None` (`stratify.py:196-197`). A truncated JSON, a zero-byte file, or a `{"target", "scores"}` file with non-numeric scores is passed over, and `passed_over` is reported only when no score file at all is found (`stratify.py:216-219`). A truncated or corrupt gzip raises `EOFError` or `zlib.error`, which neither the scan nor `run()` (`stratify.py:342`) catches, so it ends in a traceback, also for a single `--scores FILE`.

**Impact.** A directory with a good `geneA.scores.json` and a truncated `geneB.scores.json` (harvest killed during the non-atomic `write_text` at `harvest.py:186`) stratifies only geneA, with empty stderr and exit 0: the cohort silently loses a target.

**Evidence.** A directory holding a good file plus a truncated, a zero-byte and a string-scored `*.scores.json`: `disassembler stratify -o DIR --score-bin 0.9:1.0:5` wrote only geneA, exit 0, nothing on stderr. A truncated `.gz` gave `EOFError: Compressed file ended…`, a byte-flipped one `zlib.error: Error -3`, both as tracebacks.

**Proposed fix.** Split "not a score file" from "a damaged score file" in `read_score_file`:
```python
class NotJsonObject(ValueError): ...      # accession lists, binaries: passed over
class DamagedScoreFile(ValueError): ...   # never passed over
try:
    with open_text(path) as fh:
        try: start = fh.read(256)
        except UnicodeDecodeError: raise NotJsonObject(f"{path}: not text")
        if not start.strip() and path.name.endswith(".scores.json"):
            raise DamagedScoreFile(f"{path}: empty")
        if not start.lstrip().startswith("{"): raise NotJsonObject(f"{path}: not a JSON object")
        data = json.loads(start + fh.read())
except (EOFError, zlib.error, gzip.BadGzipFile, UnicodeDecodeError, json.JSONDecodeError) as e:
    raise DamagedScoreFile(f"{path}: cannot be read ({type(e).__name__}: {e})") from None
# after the two valid-shape returns:
if isinstance(data, dict) and ("scores" in data or "target" in data) and "bins" not in data:
    raise DamagedScoreFile(f"{path}: has a target/scores key but scores is not {{accession: number}}")
```
In the directory loop put `except DamagedScoreFile: raise` before `except (OSError, ValueError)`. Make `harvest` write score files atomically (temporary file plus `os.replace`). A normal 0.1.1 harvest directory, including stratify's own outputs, still stratifies.

**Tests and docs.** In `tests/test_stratify.py::TestScoreFiles` add truncated JSON, zero-byte, string-scored, truncated gzip and byte-flipped gzip next to a good file (exit 1, path in stderr, no traceback). `docs/kmindex.md:66-73`: damaged score files are an error. Related: PH-12 (a BOM is dropped the same silent way).

### PH-07 · medium · Targets are renamed by a phiX174 substring strip, and name collisions are merged silently

`src/disassembler/cmd/harvest.py:75` · correctness · confirmed (reproduced) · fix effort: small · sources: P3-06, D2-08

**Problem.** `safe_target_name` rewrites the first header token: `name = tokens[0].replace("phiX174_", "").replace("phiX174", "")` (`harvest.py:75`), then maps `/` and `\` to `_` and strips edge dots and spaces. The cleaned name is the merge key (`harvest.py:165-166`), so two records that clean to the same string are max-merged into one target without a warning. The docstrings claim the spike-in "is excluded by TARGET NAME" (`harvest.py:18-19`, `stratify.py:16-17`), but only a token exactly `phiX174` is dropped; `phiX174_x` becomes `x`, `my_phiX174_ctrl` becomes `my_ctrl`, and the RefSeq header `NC_001422.1 Escherichia phage phiX174` is kept as `NC_001422.1`. A target named `UNION` is overwritten by the union list, and one named `_union` loses its summary entry.

**Impact.** A panel with `phiX174_gp1` and `gp1` yields one `gp1` list headed by the spike-in hit, and `from-kmindex --target phiX174_gp1` fails with "expected accession list not found". `lineage/B.1` plus `lineage_B.1`, or `geneX.` plus `geneX`, merge the same way. Unusual names are needed, but the result is a wrong cohort with exit 0.

**Evidence.** One shard holding `phiX174_gp1`, `gp1`, `my_phiX174_ctrl`, `lineage/B.1`, `lineage_B.1`, `UNION`, `_union`: `harvest --shards-dir` wrote `gp1.accessions.txt` with both accessions, `my_ctrl`, one `lineage_B.1`, and a `UNION.accessions.txt` and `summary.json` that lost the `UNION` and `_union` targets; exit 0.

**Proposed fix.** (1) Drop the substring strip. (2) Exclude the spike-in explicitly: `is_spike_in(raw)` when the first token fullmatches `(?i)phix174|NC_001422(\.\d+)?`, logged as "excluded (phiX174 spike-in)". (3) Key `per_target` by the raw first token during the merge and assign file names once at write time: process keys sorted, take `safe_target_name(key)`, append `__2`, `__3` on a case-insensitive clash with an earlier name or with the reserved `union`/`_union`, and warn `target 'X' written as 'Y'`. Sorting keeps the result independent of arrival order (PH-01). Keep the `_union` key, which the PR #8800 wrapper asserts.

**Tests and docs.** Replace `tests/test_harvest.py:359-363` (which pins `phiX174_gp1` -> `gp1`) with cases for the new rule, collision suffixes, and `UNION`/`_union` targets. Fix `harvest.py:18-19`, `stratify.py:16-17` and `docs/kmindex.md:41`. Changelog: `phiX174_<gene>` targets keep their prefix. Related: PH-19 changes the same function.

### PH-08 · medium · A Galaxy collection with zero elements harvests to an empty cohort with exit 0

`src/disassembler/cmd/harvest.py:385` · correctness · confirmed (reproduced) · fix effort: trivial · sources: P3-05

**Problem.** `harvest_collection` walks `dc.get("elements", [])` (`harvest.py:317`) and never checks that it found anything; the only no-data guard, `if items and n_failed == len(items):` (`harvest.py:385`), is false for an empty list. The run writes a one-newline `UNION.accessions.txt` and `{"_union": {"n": 0}}`, prints `[✓] Harvest complete: 0 unique SRA runs across 0 targets.` and returns 0. Local mode refuses the same situation (`harvest.py:235-236`).

**Impact.** Pointing harvest at a kmindex output collection that is not yet populated (job queued or failed; Galaxy returns `elements: []`) gives a green check and exit 0 on an empty cohort; a chained script sees success, and `from-kmindex` fails later with a misleading message.

**Evidence.** A stand-in server answering `{"populated": false, "populated_state": "new", "elements": []}`: `harvest --history-id h1 --collection-id c1 --base-url http://127.0.0.1:PORT` printed "Total shard datasets in collection: 0" and "Harvest complete", exit 0.

**Proposed fix.** Right after the walk, before `out_dir.mkdir`:
```python
if not items:
    state = dc.get("populated_state") if isinstance(dc, dict) else None
    raise RuntimeError(
        f"collection {coll_id} holds no shard datasets"
        f"{f' (populated_state={state!r})' if state else ''}; has the kmindex "
        f"job finished, and is this the right collection?")
```
`run()` already maps RuntimeError to exit 1.

**Tests and docs.** Add a `TestHarvestAccounting` case with `elements: []`: `run()` returns 1 and no `UNION.accessions.txt` is written. One line in `docs/kmindex.md`: an empty or unpopulated collection is an error.

### PH-09 · medium · Docs steer the Galaxy API key into world-readable /tmp/gxy.txt; `--api-key` exposes it on argv

`src/disassembler/cmd/harvest.py:47` · security · confirmed (reproduced) · fix effort: small · sources: PX-02, SWD-04

**Problem.** Both copies of `get_galaxy_api_key` fall back to `default_tmp = Path("/tmp/gxy.txt")` and `return default_tmp.read_text().strip()` (`harvest.py:47-49`, `lexicmap.py:138-140`) with no owner or mode check and no message naming the source. The error text, the `--api-key-file` help (`harvest.py:415-417`), `validators.py:144`, `docs/kmindex.md:26` and `docs/lexicmap.md:88` all point users there. `lexicmap` and `from-lexicmap` also take `--api-key KEY` as a literal (`lexicmap.py:1212`), which lands in `/proc/PID/cmdline`, `ps` and shell history. A Galaxy API key grants full account access.

**Impact.** On a shared login or compute node, `echo KEY > /tmp/gxy.txt` under umask 022 creates a 0644 file every local user can read, and `--api-key KEY` is visible to every local user for the whole run. With no key configured, a file another user left at `/tmp/gxy.txt` is used silently and sent to `--base-url`. The PR #8800 Galaxy wrapper uses local mode and is not affected.

**Evidence.** During `disassembler lexicmap -q ref.fa -d d9 --base-url http://127.0.0.1:PORT --api-key ARGV-SECRET`, `/proc/<pid>/cmdline` (mode 0444) contained `--api-key ARGV-SECRET`. The `/tmp` fallback was confirmed by reading the code.

**Proposed fix.** One shared helper (e.g. `disassembler/galaxy.py: resolve_api_key(api_key, key_file)`) replacing both copies, which:
1. drops `/tmp/gxy.txt` (optionally reads `${XDG_CONFIG_HOME:-~/.config}/disassembler/galaxy_api_key`; for one release, read `/tmp/gxy.txt` with a deprecation warning);
2. refuses or warns when a key file has `st.st_uid != os.getuid()` or `st.st_mode & 0o077`, suggesting `chmod 600` (warn only on filesystems without POSIX modes);
3. keeps `--api-key` with a stderr warning that it is visible to other users, and accepts `--api-key-file -` for stdin;
4. logs where the key came from, never its value.
PH-15 and PH-23 belong in the same helper.

**Tests and docs.** Extend the precedence tests in `tests/test_lexicmap.py` (near line 728): `/tmp/gxy.txt` no longer read, 0644 rejected or warned, 0600 accepted, `--api-key` warns. Update `docs/kmindex.md:25-26`, `docs/lexicmap.md:88,184`, the help strings in `harvest.py`, `lexicmap.py`, `pipe_from_kmindex.py`, `pipe_from_lexicmap.py`, and `validators.py:144`.

### PH-10 · medium · urllib3 retries 429/503 with Retry-After inside each attempt: 4x requests, 120 s cap bypassed

`src/disassembler/cmd/harvest.py:324` · reliability · confirmed (reproduced) · fix effort: trivial · sources: P3-03

**Problem.** `harvest_collection` mounts `HTTPAdapter(..., max_retries=3)` (`harvest.py:324`). requests turns that into `Retry.from_int(3)`, and urllib3 retries any GET answered 413, 429 or 503 with a `Retry-After` header, sleeping the server's value (urllib3 2.x caps it at 21,600 s; older versions do not cap it). Each of `fetch_shard`'s `ATTEMPTS = 5` therefore costs up to 4 requests and 3 uncapped sleeps before `fetch_shard` sees a 429, and the docstring's promise (`harvest.py:100-102`, capped at `RETRY_AFTER_CAP_SECS` with `ATTEMPTS` bounding the per-shard time) is false.

**Impact.** Under usegalaxy.org rate limiting (429 with `Retry-After: 30-60`, which the `--download-workers` help says happens above 6 workers), a throttled shard costs 20 requests and up to about 19 minutes instead of 5 requests and 4 minutes, and the extra requests keep feeding the limiter. A server sending `Retry-After: 3600` parks a worker for hours despite the 120 s cap. Output is not wrong; time and load are.

**Evidence.** A local server answering every GET with 429 and `Retry-After: 1`, session built as in `harvest_collection` (requests 2.34.2, urllib3 2.8.0): one `fetch_shard` made 20 requests in 19.0 s; with `ATTEMPTS=1`, `RETRY_AFTER_CAP_SECS=1` and `Retry-After: 3` it made 4 requests in 9.0 s.

**Proposed fix.** Keep urllib3 to connection-level retries and leave HTTP statuses to `fetch_shard`:
```python
from urllib3.util.retry import Retry
adapter = HTTPAdapter(pool_connections=max_workers * 2, pool_maxsize=max_workers * 2,
                      max_retries=Retry(total=3, read=False, status=0,
                                        respect_retry_after_header=False))
```
Optionally widen `fetch_shard`'s special case to `r.status_code in (429, 503) and "Retry-After" in r.headers`. With the fix: 5 requests in 4 s.

**Tests and docs.** Factor the session into `_make_session(max_workers)` (current tests replace `HTTPAdapter` with a lambda and never exercise it) and add a test with a local `http.server` returning 429 plus `Retry-After: 1` asserting exactly `ATTEMPTS` requests. No doc change. Related: PH-17, PH-22.

### PH-11 · medium · stratify parses and keeps every target's score map in memory, even with `--target`

`src/disassembler/cmd/stratify.py:195` · efficiency · confirmed (reproduced) · fix effort: small · sources: PF-03, PB-03

**Problem.** `find_score_files` fully parses every file in the directory (`read_score_file`, `json.loads(start + fh.read())` at `stratify.py:164`) and appends the whole map to `found` (`stratify.py:201`); the `--target` filter runs only afterwards (`stratify.py:222-228`). All targets' dicts are alive at once. 0.1.0 loaded one score file at a time, so this is a regression. `from-kmindex --score-bin` always scans the whole cohort directory.

**Impact.** Peak memory and time scale with the sum over all targets, not the one requested: about 60-65 MB per unrequested 500,000-accession target, several GB for a large multi-record harvest (e.g. the 160-record `queries/master_planetary_targets.fasta`), enough to hit job memory limits. Output is correct.

**Evidence.** Six harvest-shaped score files of 500,000 accessions each: `stratify -o many --target g0 --score-bin 0.9:1.0:100` took 5.1 s and 457 MB max RSS; the same with only g0 in the directory, 1.5 s and 137 MB. The patched version: 1.6 s and 137 MB with identical outputs.

**Proposed fix.**
1. `ScoreFile.scores: dict | None = None` and a `load()` method that re-reads `path` when `scores` is None.
2. `read_score_file(path, only_target=None)`: after the 256-character head read, match `\s*\{\s*"target"\s*:\s*("(?:[^"\\]|\\.)*")\s*,\s*"scores"\s*:\s*\{`; if it names another target, return `ScoreFile(path, that, None)` without parsing the body (requiring `"scores"` keeps `*.bins.json` out).
3. `find_score_files` passes `target` down and keeps only `(path, target)`.
4. `stratify_bins`: `scores = sf.load()`, compute `n_available_nonzero` once (now computed twice, `stratify.py:276`, `:289`), `del scores` after each target.
Without `--target` each file is then parsed twice (about 40% slower, but 152 MB instead of 462 MB); keep the dict from discovery when only one file was found if that matters.

**Tests and docs.** Test that `read_score_file(other, only_target="a").scores is None`, that `find_score_files([dir], "a")` entries `load()` the map, and that another target's `*.bins.json` is still passed over; update `tests/test_stratify.py:189` to compare `f.load()`. No user docs change.

### PH-12 · low · A UTF-8 BOM makes a score file unreadable, or silently drops it from a directory

`src/disassembler/textio.py:32` · correctness · confirmed (reproduced) · fix effort: trivial · sources: SWP3-05

**Problem.** `text_stream` decodes with `encoding="utf-8"` (`textio.py:32`), which keeps a leading U+FEFF; `read_score_file`'s `start.lstrip().startswith("{")` (`stratify.py:162`) then fails. Named directly, a BOM-prefixed score file (PowerShell 5 `Set-Content -Encoding UTF8`, "UTF-8 with BOM" editors) gives a false "not a JSON object"; in a directory it is passed over (`stratify.py:196-197`) and its target vanishes with exit 0. A BOM on a shard makes it a counted failed shard.

**Evidence.** `printf '\xef\xbb\xbf{"target":"geneB","scores":{"SRR3":0.5}}' > bom/geneB.scores.json` next to a plain geneA file: `stratify --scores bom --score-bin 0:1.0:5 -o out` stratified only geneA, exit 0.

**Proposed fix.** `return io.TextIOWrapper(buf, encoding="utf-8-sig")`. This also covers shards and FASTA read through `textio`.

**Tests and docs.** A BOM score file found directly and in a directory; `read_text` of `b"\xef\xbb\xbf>x\nACGT\n"`, plain and gzipped, starts with `>`. Related: PH-06.

### PH-13 · low · stratified_sample is top-anchored: N=1 takes the bin maximum, the bin bottom is never drawn

`src/disassembler/cmd/stratify.py:58` · science · confirmed (reproduced) · fix effort: trivial · sources: SWP3-02

**Problem.** Bins are sorted by descending score (`stratify.py:123-126`) and `stratified_sample` returns `items[int(i * step)] for i in range(k)` (`stratify.py:58`), so index 0, the bin maximum, is always drawn and the last `n/k - 1` items never are. `docs/kmindex.md:54-55` promises a sample "spread evenly across the sorted range". The bias is systematic toward better-matching runs: for 1,000 evenly spaced scores in (0, 1], the mean pick is 1.0 at N=1 and 0.60 at N=5, against a bin mean of 0.50. Narrow bins (0.05) keep the shift small.

**Evidence.** Ten accessions scoring 0.01-0.10: `stratify --score-bin 0.0:0.11:1` picked SRR10 (0.10); `:3` picked SRR10, SRR07, SRR04.

**Proposed fix.** Take the middle of each stratum:
```diff
-    return [items[int(i * step)] for i in range(k)]
+    return [items[int((i + 0.5) * step)] for i in range(k)]
```
With it, N=1 picks SRR05 and N=3 picks SRR09, SRR05, SRR02.

**Tests and docs.** Assert `stratified_sample(list(range(10)), 3) == [1, 5, 8]` and `(…, 1) == [5]`; add "(the middle accession of each of N equal rank strata)" to `docs/kmindex.md:54-55`. Changelog: the same `--score-bin` now selects different accessions. The docstring rewrite is PH-26.

### PH-14 · low · The x-api-key header is re-sent on cross-host redirects; an http:// base URL is accepted silently

`src/disassembler/cmd/harvest.py:328` · security · confirmed (reproduced) · fix effort: small · sources: PX-03, P3-11

**Problem.** The key travels as a custom header (`harvest.py:290`, `:296`, `sess.headers.update(headers)` at `:328`). requests strips only `Authorization` when a redirect changes host, so `x-api-key` follows a 3xx anywhere; `lexicmap.py:244`, `:266` build `urllib.request.Request(url, headers={"x-api-key": key})`, which urllib also copies on redirect. Nothing checks that `--base-url` is https; `http://usegalaxy.org` answers 302 to https after the key has crossed in clear. A cross-host redirect (SSO front end, domain move, object store) hands over full account access. Stock usegalaxy.org does not redirect these calls to another host.

**Evidence.** A server on 127.0.0.1 answering 302 to `localhost:B`: both `harvest` and `lexicmap` sent `x-api-key` to B.

**Proposed fix.**
```python
def _strip_key_on_host_change(sess):          # instance patch keeps test stubs working
    orig = sess.rebuild_auth
    def rebuild_auth(prep, resp):
        orig(prep, resp)
        if sess.should_strip_auth(resp.request.url, prep.url):
            prep.headers.pop("x-api-key", None)
    sess.rebuild_auth = rebuild_auth
    return sess
# listing: requests.get(..., allow_redirects=False); on 3xx raise
#   "Galaxy at {base_url} answered with a redirect; point --base-url at the final https:// address"
# lexicmap.py:244,266: req = urllib.request.Request(url); req.add_unredirected_header("x-api-key", key)
```
Warn, or refuse without `--allow-insecure-http`, for `http://` to a non-loopback host. `add_unredirected_header` also drops the key on same-host redirects; a custom redirect handler can re-add it for same scheme and host.

**Tests and docs.** Two local servers (A redirects to B): B never sees the key; a listing redirect raises. One line in `docs/kmindex.md` and `docs/lexicmap.md`: use an https base URL.

### PH-15 · low · A Galaxy API key with an embedded newline is echoed in full in the error message

`src/disassembler/cmd/harvest.py:44` · security · confirmed (reproduced) · fix effort: small · sources: PX-05

**Problem.** Keys are only trimmed at the ends (`harvest.py:44`, `:46`, `:49`; `lexicmap.py:133-140`), and `validate_api_key_file` rejects only empty content. A key file with a comment line above the key (or a multi-line `GALAXY_API_KEY`) makes requests or urllib reject the header before any connection, quoting the whole value, and `harvest.py:299` / `lexicmap.py:271` print that exception to stderr, which may be a job, CI or scheduler log.

**Evidence.** `printf '# my usegalaxy key\nSECRET-KEY-0123\n' > k1.txt; disassembler harvest --history-id h --collection-id c --base-url http://127.0.0.1:1 --api-key-file k1.txt -o out` printed `InvalidHeader: … '# my usegalaxy key\nSECRET-KEY-0123'`.

**Proposed fix.** One check used by every key source (and by `validate_api_key_file`):
```python
def check_api_key(key, source):
    key = key.strip()
    if not key:
        raise ValueError(f"the Galaxy API key from {source} is empty")
    if any(c.isspace() or not c.isprintable() for c in key):
        raise ValueError(f"the Galaxy API key from {source} is not a single token "
                         f"({len(key.splitlines())} line(s)); the value is not shown")
    return key
```
Do not guess the key from a "first non-comment line". As defence in depth, replace the key with `<api key>` in the two `RequestException` messages. Belongs in PH-09's shared helper.

**Tests and docs.** A key file `"#c\nKEY"` and a two-line `GALAXY_API_KEY`: non-zero exit, `KEY` absent from stderr. Docs: the key file holds the key alone on one line.

### PH-16 · low · stratify accepts out-of-range and non-finite bins and exits 0 having selected nothing

`src/disassembler/cmd/stratify.py:76` · usability · confirmed (reproduced) · fix effort: small · sources: P3-13

**Problem.** `parse_bin_spec` checks only `lo < hi` (`stratify.py:76`) and `n >= 0`, and `run()` returns 0 unconditionally. A percent-style `--score-bin 90:100:50` gives an empty `<target>.bins.txt` and exit 0; `0.9:inf:5` writes the non-standard token `Infinity` into `bins.json`. Since 0.1.1 made a bin ending at 1.0 inclusive, `0.9:1.0` next to `1.0:1.01` is refused with advice that does not help, while `README.md:68` and `docs/galaxy.md:18` still show the 1.01 convention.

**Proposed fix.**
```python
if not (math.isfinite(lo) and math.isfinite(hi)):
    raise ValueError(f"--score-bin '{s}': LO and HI must be finite numbers")
if hi <= 0.0 or lo > TOP_SCORE:
    raise ValueError(f"--score-bin '{s}': kmindex scores run from 0 to 1, "
                     "so no nonzero score can fall in this bin")
```
In `assert_bins_disjoint`, for the `a_hi == 1.0 and b_lo == 1.0` clash, add: "a bin ending at 1.0 includes 1.0; end the lower bin just below 1.0 (e.g. 0.999999)". Have `run()` warn on stderr when nothing was selected (keep exit 0).

**Tests and docs.** `parse_bin_spec('90:100:5')` and `('0.9:inf:5')` raise; the overlap message mentions 0.999999. Change `README.md:68` to `0.5:1.0:100` and `docs/galaxy.md:18` to "1.0 (inclusive)".

### PH-17 · low · Non-JSON 200 shards are retried 5 times, then the API key is blamed

`src/disassembler/cmd/harvest.py:113` · reliability · confirmed (reproduced) · fix effort: small · sources: SWP3-01

**Problem.** `fetch_shard` decodes inside the retry `try` (`harvest.py:113-114`), so a `JSONDecodeError` from a body that is not JSON (for example kmindex `--format matrix` TSV) goes round the 5-attempt backoff like a network error. When all shards fail, the message says `Check the Galaxy API key and the collection id` (`harvest.py:388`), although the listing already succeeded with both. At default settings and about 2,900 shards that is about 11,500 extra GETs against the rate-limited display endpoint before a misleading error.

**Evidence.** Four TSV shards served with HTTP 200: 20 GETs, 16 sleeps, then the API-key message.

**Proposed fix.** Classify the decode failure:
```python
class ShardFormatError(ValueError): ...
if r.status_code == 200 and r.content.strip():
    text = decode_bytes(r.content)
    try:
        return json.loads(text)
    except ValueError as e:
        if not text.lstrip().startswith(("{", "[", "<")):   # truncation/proxy page: retry
            raise ShardFormatError(f"shard {name}: HTTP 200 but not JSON ({e})") from e
        raise
last = "HTTP 200 with an empty body" if r.status_code == 200 else f"HTTP {r.status_code}"
# in the loop: except ShardFormatError: raise
```
In `harvest_collection`, count these and, when all shards are format errors, say "harvest reads kmindex query output in JSON format"; reword the all-failed message to point at the per-shard errors.

**Tests and docs.** A TSV 200 gives one GET per shard and no sleeps; a `{"trunc` body still retries. `docs/kmindex.md`: harvest needs kmindex JSON output. Related: PH-04, PH-10, PH-18.

### PH-18 · low · A trailing-slash `--base-url` or a non-JSON 200 listing ends in a traceback

`src/disassembler/cmd/harvest.py:303` · reliability · confirmed (reproduced) · fix effort: trivial · sources: P3-07

**Problem.** `harvest_collection` builds `f"{base_url}/api/…"` (`harvest.py:291`, also `:105`) without stripping a trailing slash, then calls `dc = r.json()` unguarded (`harvest.py:303`); `run()` catches only RuntimeError. usegalaxy.org answers `//api/…` with `200 text/html`, so `--base-url https://usegalaxy.org/` ends in a `JSONDecodeError` traceback with no hint. A JSON list or an element with `"object": null` (`harvest.py:308`) raises `AttributeError`/`TypeError`. `lexicmap` already strips the slash (`lexicmap.py:243`).

**Evidence.** `curl -s -o /dev/null -w '%{content_type}\n' https://usegalaxy.org//api/version` prints `text/html`; the harvest above ends in `requests.exceptions.JSONDecodeError`.

**Proposed fix.**
```python
base_url = base_url.rstrip("/")
try:
    dc = r.json()
except ValueError as e:
    raise RuntimeError(f"{url} answered HTTP 200 but not with JSON (Content-Type "
                       f"{r.headers.get('Content-Type')!r}); check --base-url") from e
if not isinstance(dc, dict):
    raise RuntimeError(f"unexpected Galaxy collection answer: JSON {type(dc).__name__}")
...
o = e.get("object") or {}
```
Better: strip once in a shared argparse type for all four `--base-url` options.

**Tests and docs.** `TestCollectionAnswer`: trailing slash, HTML 200, JSON list. Related: PH-17.

### PH-19 · low · Shard names, accessions and scores are trusted: a NUL or overlong target name aborts with a traceback

`src/disassembler/cmd/harvest.py:76` · robustness · confirmed (reproduced) · fix effort: small · sources: PX-07

**Problem.** `safe_target_name` only replaces `/` and `\` (`harvest.py:76-77`). A target token with a NUL, or longer than the file-name limit minus the 15-byte suffix, makes `acc_file.write_text` (`harvest.py:184`) raise `ValueError`/`OSError`, which `run()` (catching only RuntimeError, `harvest.py:475`, `:499`) lets escape after every shard was fetched: no UNION list, no summary. `_merge_payload` also accepts `float("inf")` scores and any accession string (`bad/acc`, `SRR5\nX`, `""`), and `json.dumps` writes `Infinity`. Only malformed or hand-edited shards trigger it; a 230+ byte first header token is the plausible case.

**Evidence.** A shard with a 260-character target name: `OSError: [Errno 36] File name too long`, exit 1, empty output directory.

**Proposed fix.** Allow-list the name: `re.sub(r"[^A-Za-z0-9._+-]", "_", tokens[0])[:100].strip(".")`, None when empty. In `_merge_payload` keep a hit only when `math.isfinite(sc) and 0.0 < sc <= 1.0` and `re.fullmatch(r"[A-Za-z0-9_-]+", acc)` (the walker's own rule), and report the dropped count. Write score JSON with `allow_nan=False`. Combine with PH-07, which rewrites the same function.

**Tests and docs.** NUL and 300-character names give writable names; `inf`, `nan`, `-1`, `"bad/acc"`, `""` are dropped and counted. Document the name rule in `docs/kmindex.md`.

### PH-20 · low · Following stratify's own advice fails: a directory plus NAME=FILE for a file inside it conflict

`src/disassembler/cmd/stratify.py:238` · usability · confirmed (reproduced) · fix effort: trivial · sources: SWP3-03

**Problem.** `find_score_files` appends one entry per source walked: the directory scan adds a bare legacy file with target None, and `NAME=FILE` adds the same path again (`stratify.py:193-212`). Nothing de-duplicates before `unnamed` is computed (`stratify.py:221`), so `--scores mixed --scores geneX=mixed/geneX.scores.json` still fails with "Name each with --scores NAME=FILE" (`stratify.py:238-241`), and a directory given twice fails with `target 'geneB' is in two score files: P and P` (`stratify.py:249-250`).

**Proposed fix.** Before line 221:
```python
by_path = {}
for sf in found:
    key = sf.path.resolve()
    prev = by_path.get(key)
    if prev is None or prev.target is None:
        by_path[key] = sf
    elif sf.target is not None and sf.target != prev.target:
        raise ValueError(f"{sf.path} is named both {prev.target!r} and {sf.target!r}")
found = list(by_path.values())
```
Reword the line 238 advice: "add `--scores NAME=FILE` for each (alongside the directory)".

**Tests and docs.** Directory plus `NAME=FILE` yields both targets; a directory given twice yields one entry per file. Related: PH-03.

### PH-21 · low · Per-target stratify runs overwrite the "combined" bins_summary.json with the last target only

`src/disassembler/cmd/stratify.py:293` · usability · confirmed (reproduced) · fix effort: trivial · sources: SWP3-07

**Problem.** `stratify_bins` writes `bins_summary.json` from this run's targets only (`stratify.py:293`), while earlier runs' `<target>.bins.*` stay. `from-kmindex` always stratifies one `--target` into the cohort directory, so after adding geneA then geneB the "combined" summary (`stratify.py:257`, `:310`, `docs/kmindex.md:78`) lists only geneB. No code reads the file; a reader would conclude geneA was never stratified. Same in 0.1.0.

**Proposed fix.**
```python
combined = {}
if summary_file.is_file():
    try:
        combined = json.loads(summary_file.read_text())
    except (OSError, ValueError):
        combined = {}
    if not isinstance(combined, dict):
        combined = {}
combined.update(all_summaries)          # a re-stratified target replaces its own entry
summary_file.write_text(json.dumps(dict(sorted(combined.items())), indent=2))
```
(Alternative: rebuild from every `*.bins.json` in `out_dir`.)

**Tests and docs.** `--target geneA` then `--target geneB` into one `-o`: summary keys are both. Reword `docs/kmindex.md:78`, `stratify.py:257`, `:310`. Related: PH-05.

### PH-22 · low · Ctrl-C during a Galaxy harvest keeps fetching every queued shard before exiting

`src/disassembler/cmd/harvest.py:341` · usability · confirmed (reproduced) · fix effort: small · sources: P3-10

**Problem.** All shards are submitted up front inside `with ThreadPoolExecutor(...) as ex:` (`harvest.py:341-342`); on KeyboardInterrupt `__exit__` calls `shutdown(wait=True)` without `cancel_futures`, so the interrupt takes effect only after the whole collection is fetched, then ends in a traceback with nothing written. With 60 shards at 1 s each and 6 workers, SIGINT at 2.5 s exited 8.2 s later after all 60 requests.

**Proposed fix.**
```python
ex = concurrent.futures.ThreadPoolExecutor(max_workers=max_workers)
try:
    ...  # submit and as_completed loop
except BaseException:
    ex.shutdown(wait=False, cancel_futures=True); raise
else:
    ex.shutdown()
finally:
    sess.close()
```
In `run()`: `except KeyboardInterrupt: print("harvest: interrupted; nothing was written"); return 130`. For a prompt exit during a 429 back-off, pass a `threading.Event` and use `stop.wait(delay)` instead of `time.sleep`. With the first two steps: exit 130 within 1.1 s, 18 of 60 requests.

**Tests and docs.** Patched `fetch_shard` counting calls plus a merge that raises KeyboardInterrupt: call count well below the item count, `run()` returns 130.

### PH-23 · low · A whitespace-only $GALAXY_API_KEY or blank /tmp/gxy.txt is sent as an empty key

`src/disassembler/cmd/harvest.py:45` · usability · confirmed (reproduced) · fix effort: small · sources: SWP3-04

**Problem.** `if os.environ.get("GALAXY_API_KEY"): return os.environ["GALAXY_API_KEY"].strip()` (`harvest.py:45-46`) tests the raw value and returns the stripped one, so `' '` or `$'\r'` yields `''`; the `/tmp/gxy.txt` branch (`:48-49`) has no check at all; `lexicmap.py:133-140` repeats this, including `--api-key " "`. Galaxy treats the empty key as anonymous and the user sees only 403 "History is not accessible by the current user". A truly empty variable falls through correctly.

**Evidence.** `GALAXY_API_KEY='   ' disassembler harvest --history-id H --collection-id C --base-url http://127.0.0.1:PORT -o out`: the server logged `x-api-key=''`, the run ended with the 403 message.

**Proposed fix.** Check the stripped value in each branch and raise `RuntimeError("$GALAXY_API_KEY is set but empty (or whitespace-only); set it to the Galaxy API key, or unset it.")` (and the same for the file); fold into PH-09's shared helper with PH-15's `check_api_key`.

**Tests and docs.** `mock.patch.dict(os.environ, {"GALAXY_API_KEY": " \n"})` makes `get_galaxy_api_key()` raise; same for a blank default file and `lexicmap`'s `api_key=" "`.

### PH-24 · low · Galaxy-mode harvest does not name the shard it discards; the per-shard loop is duplicated

`src/disassembler/cmd/harvest.py:360` · reuse · confirmed (reproduced) · fix effort: small · sources: PS-03

**Problem.** `harvest_local` and `harvest_collection` carry separate copies of the per-shard loop, and they have drifted: local mode names a non-object shard (`harvest.py:261-263`) and a failed merge (`:272`), while collection mode does `if not isinstance(d, dict): n_failed += 1; continue` silently (`:360-362`) and prints `shard payload unusable` without a name (`:372`), although the name is in `futures[fut]`. Users cannot tell which of thousands of shards to re-fetch.

**Proposed fix.**
```python
def _fold_shard(label, d, per_target) -> tuple[str, int]:
    """Return ('failed'|'empty'|'ok', n_bad_target); report failures naming `label`."""
    if not isinstance(d, dict):
        print(f"    [!] shard {label}: JSON {type(d).__name__}, not a kmindex result object",
              file=sys.stderr)
        return "failed", 0
    if not d:
        return "empty", 0
    try:
        return "ok", _merge_payload(d, per_target)
    except Exception as e:
        print(f"    [!] shard {label} payload unusable ({type(e).__name__}: {e})",
              file=sys.stderr)
        return "failed", 0
```
Use `label = futures[fut][0]` in collection mode, name the shard in "shard fetch raised", and move the progress line into a `_progress()` helper. Do this before PH-02 and PH-04 so their changes are made once.

**Tests and docs.** A list payload and a payload whose merge raises: stderr names the element identifier; local and collection modes give the same counts and messages.

### PH-25 · low · `--scores` rejects pipes as "no such file"; `--shards-dir FILE` says "not found"

`src/disassembler/cmd/stratify.py:202` · usability · confirmed (reproduced) · fix effort: trivial · sources: SWP3-06

**Problem.** `find_score_files` accepts only `is_dir()` or `is_file()` (`stratify.py:189`, `:202`); anything else, including `<(zcat x)` and `/dev/stdin` on a pipe, gets `no such file or directory` (`:214`), although `textio` is written to read pipes. `harvest --shards-dir merged.json` reports `shards directory not found` for an existing file (`harvest.py:231-232`).

**Proposed fix.** `stratify.py:202`: `elif path.exists():` (the read is single-pass). In `harvest_local`:
```python
if not shards_dir.exists():
    raise RuntimeError(f"shards directory not found: {shards_dir}")
if not shards_dir.is_dir():
    raise RuntimeError(f"--shards-dir {shards_dir} is not a directory; "
                       f"pass the directory that holds the shard files")
```

**Tests and docs.** `find_score_files([fifo])` fed by a writer thread returns one entry; `harvest_local(<file>)` raises "not a directory".

### PH-26 · low · stratified_sample's docstring, its test and cli.py describe removed code

`src/disassembler/cmd/stratify.py:46` · docs · confirmed (by reading) · fix effort: trivial · sources: PS-14

**Problem.** The `stratified_sample` docstring justifies the `k <= 0` guard with `quota_canon`, `--diversity` and `stratify_diversity`'s top-up (`stratify.py:46-51`), all removed (the module docstring says so at `stratify.py:19-24`) and contradicting the "no backfill" contract; `tests/test_stratify.py:21` repeats it. `cli.py:11` says subcommands "map 1:1 to the current scripts", but there are no scripts. `grep -rn "quota_canon\|stratify_diversity" src tests` finds only the docstring. The guard is still needed because `LO:HI:0` is accepted.

**Proposed fix.**
```python
def stratified_sample(items, k):
    """k items spread evenly over a score-sorted list. k == 0 (a `LO:HI:0` bin, which
    parse_bin_spec accepts) returns nothing rather than dividing by zero; a list no
    longer than k is returned whole. There is no backfill."""
```
`tests/test_stratify.py:21`: "A zero quota is legitimate: --score-bin LO:HI:0 is accepted." `cli.py:11`: "Atomic subcommands are the individual pipeline stages (harvest, stratify, fetch, walker, lexicmap, sieve, msa)." If PH-13 lands, describe the middle-of-stratum rule here too.

**Tests and docs.** The existing zero-quota test pins the guard; no user docs change.

### PH-27 · low · harvest keeps every parsed shard alive until the end and builds each output file as one string

`src/disassembler/cmd/harvest.py:342` · efficiency · confirmed (reproduced) · fix effort: small · sources: PF-06

**Problem.** `harvest_collection` keeps a dict of all futures, `futures = {ex.submit(fetch_shard, it, base_url, sess): it for it in items}` (`harvest.py:342`), for the whole executor block. Each finished `Future` holds its parsed shard payload, so every shard dict stays resident after it has been merged into `per_target`. `_write_outputs` then builds each file in memory before writing it: `acc_file.write_text("\n".join(accs) + "\n")` (`harvest.py:184`), the score JSON via `json.dumps(..., indent=2)` (`:186`) and the union list (`:195`), and fills `all_union` (`:181`) even for a single target, where it duplicates that target's key set. `harvest_local` has no futures but shares the write step.

**Impact.** At the scale `docs/kmindex.md:5` and `:128-129` describe (about 2,869 shards, about 5.1 million accessions for one target) peak RSS is about 1.42 GB, of which about 0.15 GB is retained futures and about 0.37 GB the whole-string write step. It grows linearly and is not a blow-up, but it is avoidable cost on the modest Galaxy job slots this step is meant for. Output is correct.

**Evidence.** A stand-in Galaxy server serving 2,869 shards with 5,163,936 accessions for one target: the current code peaked at 1,423.6 MB max RSS (53.6 s), the same run with the futures released and the writes streamed at 1,274.7 MB (54.4 s) with an identical `UNION.accessions.txt`. `_write_outputs` alone on 5,163,936 accessions: +662 MB over its starting RSS now, +296 MB when writing through file handles; all four output files byte-identical.

**Proposed fix.**
1. Release each future as it is consumed. `as_completed` keeps its own set and drops each future as it yields it (CPython 3.13), so the caller's dict is the only lasting reference. Pop it there, which also keeps the shard identity PH-24 needs for its messages:
```python
for fut in concurrent.futures.as_completed(futures):
    item = futures.pop(fut)          # drop our reference; `item` names the shard
```
2. In `_write_outputs`, write through handles: `fh.writelines(a + "\n" for a in accs)`, `json.dump({...}, fh, indent=2)`; build `all_union` only when `len(per_target) > 1` and otherwise reuse the single target's keys sorted. Keep writing a lone `"\n"` for an empty union so output stays byte-identical.

**Tests and docs.** Existing harvest tests compare output files and should pass unchanged. Add a test that `_write_outputs` on zero, one and two targets produces files byte-identical to the current implementation (the empty case must still write a single `"\n"`). No docs change. Related: PH-22 and PH-24 rework the same executor loop, and PH-01 and PH-06 change the same writes (tie order, atomic replace); make these changes together.
