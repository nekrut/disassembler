# Galaxy: tools-iuc PR #8800 wrappers, docs/galaxy.md, and the Galaxy delivery plan sections 2 to 4

This section covers the six tool wrappers in galaxyproject/tools-iuc PR #8800 (paths written `tools-iuc#8800/<file>` are relative to `tools/disassembler/` at PR head b9a16ad), the repository's `docs/galaxy.md` and `docs/packaging.md` release step, and sections 2 to 4 of the Galaxy delivery plan (`GALAXY_PLAN.md`, 2026-11-04). The wrappers fit the 0.1.0 packages they pin: commands render as intended, outputs are collected correctly, and most help text is accurate for that release. They are fragile in three ways. First, they encode 0.1.0 behaviour that 0.1.1 changed (test data, the primary-metric default, help text), so the bump that `docs/packaging.md` describes as a token change is in fact a revision of the wrappers. Second, they make decisions the user cannot see or control: the memory fit silently lowers the size cap, and the sieve applies walker-tuned gates to LexicMap rows. Third, one staging choice lets the walker overwrite its own input. The plan's UDT and workflow sketches repeat several of these problems and misstate what Galaxy and the Tool Shed provide. The fixes share one theme: expose as parameters what the CLI composites derive internally, diff wrapper defaults and help against the CLI at every release, and move the memory model into the walker (`--memory-mb`) so that no wrapper keeps its own copy.

| ID | Severity | Where | Title | Fix effort |
|---|---|---|---|---|
| GX-01 | high | `tools-iuc#8800/disassembler_walker.xml:55` | Walker wrapper overwrites the Reference CDS dataset when it is named `reconstructed` | trivial |
| GX-02 | medium | `tools-iuc#8800/disassembler_sieve.xml:21` | Galaxy LexicMap route drops partly covered rows: sieve and msa run with walker-tuned gates | small |
| GX-03 | medium | `tools-iuc#8800/disassembler_sieve.xml:36` | Sieve wrapper's default 95% length gate drops walker primaries admitted at 85-95% coverage | small |
| GX-04 | medium | `tools-iuc#8800/disassembler_walker.xml:81` | Walker wrapper forces `--primary-metric cov` and cannot select 0.1.1's `matched` | trivial |
| GX-05 | medium | `tools-iuc#8800/disassembler_walker.xml:34` | Walker wrapper lowers `--max-size-mb` to fit job memory, so the cohort depends on the destination | small |
| GX-06 | medium | `walker/src/stages/stream.rs:110` | Fetch-then-offline runs report no-data and over-cap accessions as errors | small |
| GX-07 | medium | `src/disassembler/cmd/harvest.py:162` | `harvest` silently drops kmindex `json_vec` scores and reports an empty harvest as success | small |
| GX-08 | medium | `tools-iuc#8800/disassembler_walker.xml:26` | Walker wrapper's memory fit counts 43 MB per graph; chosen settings exceed the allocation | trivial |
| GX-09 | medium | `tools-iuc#8800/disassembler_walker.xml:100` | `DISASSEMBLER_WALKER_MAX_SIZE_MB` cannot be overridden from job_conf | trivial |
| GX-10 | medium | `docs/packaging.md:66` | The 0.1.1 bump is not a token change: three PR tests fail, and galaxy.md already documents 0.1.1 | small |
| GX-11 | medium | `tools-iuc#8800/disassembler_stratify.xml:38` | Stratify wrapper names targets after datasets, which 0.1.1 no longer honours | small |
| GX-12 | medium | `tools-iuc#8800/macros.xml:3` | `@WALKER_VERSION@` token defeats `planemo autoupdate` | small |
| GX-13 | medium | `GALAXY_PLAN.md:125` | Plan's interim sieve UDT (bare `l1_pipeline_serial`) returns 0 or 1 records with exit 0 | trivial |
| GX-14 | medium | `GALAXY_PLAN.md:122` | Plan's walker UDT cannot reuse the IUC memory block; ported, it cuts the cap to 497 MB | small |
| GX-15 | medium | `GALAXY_PLAN.md:138` | Plan's workflows reference UDTs by id, which Galaxy cannot resolve; import caveat reversed | small |
| GX-16 | low | `tools-iuc#8800/disassembler_walker.xml:91` | Walker wrapper loses `walker_summary.json` in failed jobs | trivial |
| GX-17 | low | `docs/galaxy.md:80` | galaxy.md offers `--allow-empty` as an alternative to filtering failed jobs | trivial |
| GX-18 | low | `src/disassembler/cmd/lexicmap.py:955` | With `--require-full-length`, log and wrapper help still say the partial tier is 50-80% | small |
| GX-19 | low | `tools-iuc#8800/disassembler_msa.xml:30` | msa wrapper stale at 0.1.1: `--msa-min-cov` mislabelled, no `--min-identity` | small |
| GX-20 | low | `tools-iuc#8800/disassembler_lexicmap.xml:41` | `--allow-frameshift-alts` help describes 0.1.0 semantics | trivial |
| GX-21 | low | `tools-iuc#8800/disassembler_msa.xml:90` | msa/lexicmap help "gap-free, only adds the tree" becomes false at 0.1.1 | trivial |
| GX-22 | low | `tools-iuc#8800/disassembler_stratify.xml:51` | Stratify N help says "evenly across the score range"; sampling is by rank | trivial |
| GX-23 | low | `docs/galaxy.md:86` | MAFFT: galaxy.md calls it an IUC tool; the plan says no 7.526 container exists | trivial |
| GX-24 | low | `docs/galaxy.md:79` | No native Galaxy step splits an accession list; galaxy.md misnames tool owners | trivial |
| GX-25 | low | `GALAXY_PLAN.md:195` | Plan caveats omit UDT destination limits (`docker_net: none`; usegalaxy.org 8 cores/29 GB, CVMFS only) | trivial |
| GX-26 | low | `GALAXY_PLAN.md:123` | Plan's fetch UDT runs `disassembler` in the logan-walker container, which has none | trivial |
| GX-27 | low | `GALAXY_PLAN.md:144` | WF1 sketch omits the reference CDS and target selection | trivial |
| GX-28 | low | `GALAXY_PLAN.md:178` | WF3-frozen quota estimate (median x 50) understates skewed cohorts | trivial |

### GX-01 · high · Walker wrapper overwrites the Reference CDS dataset when it is named `reconstructed`

`tools-iuc#8800/disassembler_walker.xml:55` · galaxy · confirmed (reproduced) · fix effort: trivial · sources: G1-01

**Problem.** The wrapper links the reference into the job directory under its sanitised dataset name, `ln -s '$query' '${reference_name}.fasta'` (:55; the name is built at :51 with a trailing `.fasta|.fas|.fa|.fna` stripped), and passes `--out-fasta reconstructed.fasta` (:66) into the same directory. For a dataset named `reconstructed`, `reconstructed.fasta`, `.fa`, `.fas` or `.fna`, both paths are the same symlink. The walker reads the query at start-up and later writes the output with `std::fs::File::create(path)` (`walker/src/io/output.rs:72`), which follows the symlink and truncates the Galaxy dataset file. That write happens at `walker/src/main.rs:558`, before `check_run` (:597), so the input is destroyed even when the job then fails.

**Impact.** `disassembler reconstruct` and `from-kmindex` name their output `reconstructed.fasta` (`src/disassembler/cmd/pipe_reconstruct.py:350`), so a user who uploads a CLI result and uses it as a reference or multi-anchor panel triggers this. On local or cluster runners with writable inputs (the default), the dataset is overwritten silently; every history and library copy shares the file, its metadata goes stale, and reruns use the wrong reference. With read-only container mounts the job fails only at the final write, after all downloads and alignment. `fasta.gz` inputs are safe because gunzip writes a copy.

**Evidence.** The wrapper's command was rendered for a dataset named `reconstructed.fasta` and run with walker 0.1.1 offline and a stand-in GraphAligner. A successful job replaced the 128-byte reference with the walker's output (`>SRRMATCH0_primary ...`); a job that found nothing (exit 1) left it at 0 bytes; under a read-only bind mount the job failed with `Read-only file system (os error 30)`. A hard-linked copy changed with it every time. The CLI has the same hole: `disassembler reconstruct -q out/reconstructed.fasta -o out` truncated the query to 0 bytes.

**Proposed fix.** Stage the reference in its own directory. The stem, and every name the walker derives from it, stays the same:
```diff
+        mkdir ref &&
         #if $query.is_of_type('fasta.gz'):
-            gunzip -c '$query' > '${reference_name}.fasta' &&
+            gunzip -c '$query' > 'ref/${reference_name}.fasta' &&
         #else:
-            ln -s '$query' '${reference_name}.fasta' &&
+            ln -s '$query' 'ref/${reference_name}.fasta' &&
         #end if
 ...
-            --query '${reference_name}.fasta'
+            --query 'ref/${reference_name}.fasta'
```
An equivalent fix writes `--out-fasta walker_out/reconstructed.fasta` and points `from_work_dir` there. With either, all three cases keep the input intact, and the control run's outputs are byte-identical. As defence in depth for 0.2, the walker should refuse, before any download, an `--out-fasta` that is the same file (same device and inode) as any `--query` (about 25 lines); this also closes the CLI variant.

**Tests and docs.** PR: a walker test whose Reference CDS is `test-data/reconstructed.fasta` (a copy of `hiv1_pol.fasta`); under containerised `planemo test` the unfixed wrapper fails. Walker: a `tests/cli.rs` case with `--out-fasta` symlinked to the query, asserting a usage error and an intact query. `docs/galaxy.md:43`: say the reference is linked inside `ref/` under its dataset name; `docs/walker.md`: note the new usage error.

### GX-02 · medium · Galaxy LexicMap route drops partly covered rows: sieve and msa run with walker-tuned gates

`tools-iuc#8800/disassembler_sieve.xml:21` · galaxy · confirmed (reproduced) · fix effort: small · sources: G1-07, IC-06, G2-02

**Problem.** The sieve wrapper runs `disassembler sieve --input input.fa --output promoted.fasta --ref ref.fa|--min-len N --threshold T --sketch S` (:21-30) and exposes no gap limit, so the completeness gate uses `DEFAULT_MAX_GAP_FRAC = 0.05` (`src/disassembler/cmd/sieve.py:30`), which counts `-` as a bad base. LexicMap clean rows sit on the reference grid with `-` wherever no hit reached and are admitted at 80% coverage, so every clean row at 80-95% coverage is dropped. Likewise `msa` keeps its 0.70 default and drops the 50-70% partials of the expanded cohort. The `from-lexicmap` composite derives both gates from the cohort gate (`src/disassembler/cmd/pipe_from_lexicmap.py:339`, `:362`), but `docs/galaxy.md:59-63` draws `lexicmap -> sieve (optional) -> msa` without them, and the lexicmap wrapper calls its clean FASTA "ready for **disassembler sieve**" (`disassembler_lexicmap.xml:256`). At the 0.1.0 pin the gate is hard-coded, so the fix needs the 0.1.1 bump.

**Impact.** Green jobs with thinned cohorts. `docs/sieve.md:61-63` puts the loss at 22% of a default C. auris FKS1 cohort and a third of the distinct SARS-CoV-2 spike haplotypes.

**Evidence.** A 20-row cohort (600 nt, every second row ending in 60 `-`) run with the wrapper's command line, `disassembler sieve --input cohort.fa --output p1.fa --ref ref.fa --threshold 0.005 --sketch 256`, printed `20 -> 10 haplotypes` and kept no gapped row; adding `--max-gap-frac 0.25` kept all 20. On 10 full rows plus 10 at 60% coverage, `from-lexicmap --include-partials` kept 20/20, while `lexicmap` followed by `msa --input clean_expanded` kept 10/20.

**Proposed fix.** At the 0.1.1 bump, in `disassembler_sieve.xml` (or a select "input comes from: walker / LexicMap clean / LexicMap expanded / custom"):
```diff
+            --max-gap-frac $max_gap_frac
             --threshold $threshold
 ...
+        <param argument="--max-gap-frac" type="float" value="0.05" min="0" max="1"
+               label="Maximum share of gaps or other non-ACGT characters"
+               help="0.05 suits walker haplotypes. LexicMap rows carry '-' where no hit reached:
+                     use 1 - coverage gate + 0.05 (0.25 for a default clean cohort)."/>
```
In the CLI, move the composite's two formulas into a shared helper and add `--cohort-gate F` to `sieve` and `msa` (used only when the explicit flag is absent):
```python
def sieve_gap_limit(gate): return round(min(1.0, 1.0 - gate + GATE_MARGIN), 6)
def msa_cov_gate(gate):    return round(min(0.70, gate - GATE_MARGIN), 6)
```
Have `lexicmap` record its gates in `lexicmap.summary.json`, and expose `--cohort-gate` in `disassembler_msa.xml`.

**Tests and docs.** Wrapper test: 8 rows of 600 nt, 4 with 15% `-`, asserting 4 promoted at the default and 8 at 0.25. CLI: `msa --cohort-gate 0.5` keeps 20/20 on the fixture above; `from-lexicmap --include-partials` equals `lexicmap` plus `msa --cohort-gate 0.5`. Docs: `docs/galaxy.md:57-64` (name the derived parameters; model the optional sieve as a boolean input with `when:`), `docs/sieve.md:74`, both wrappers' help.

### GX-03 · medium · Sieve wrapper's default 95% length gate drops walker primaries admitted at 85-95% coverage

`tools-iuc#8800/disassembler_sieve.xml:36` · galaxy · confirmed (reproduced) · fix effort: small · sources: SWG-01

**Problem.** The sieve tool's default length gate is "95% of the reference CDS length" (:36, rendered as `--ref ref.fa`; `auto_min_len` returns `int(length * 0.95)`, `src/disassembler/cmd/sieve.py:52`). The walker tool admits primaries down to `--min-cov` 0.85 (`disassembler_walker.xml:119`), where cov is haplotype length over reference length (`walker/src/stages/stitch.rs:198`). Every primary with 0.85 <= cov < 0.95 is removed by the completeness gate before thinning, although `msa` (0.70) would keep it.

**Impact.** On the suite's walker -> sieve -> msa route, accessions whose walk stops short of a terminus or carries a net deletion vanish from the alignment and tree; the job is green and one stdout line records it. How often real cohorts have such primaries was not measured.

**Evidence.** A 900-nt reference with primaries of 900, 900, 810 and 783 nt, run as `disassembler sieve --input input.fa --output promoted.fasta --ref ref.fa --threshold 0.005 --sketch 256`: `completeness gate drops 2 of 4: 2 shorter than 855 nt`, `4 -> 2 haplotypes`, exit 0. With the fix: gate `<765 nt`, `4 -> 4`.

**Proposed fix.**
```diff
+DEFAULT_MIN_LEN_FRAC = 0.85   # walker --min-cov default; cov = len(hap)/len(ref)
-def auto_min_len(ref: Path) -> int | None:
+def auto_min_len(ref: Path, frac: float = DEFAULT_MIN_LEN_FRAC) -> int | None:
 ...
-    return int(length * 0.95)
+    return int(length * frac)
```
Add `--min-len-frac` to `disassembler sieve`. In the wrapper, replace the 95% option with a float "Minimum length as a fraction of the shortest reference record" (default 0.85; "set it to the walker's Minimum reference coverage"). Because the sieve keeps the first sufficiently novel record, have `byte_exact_dedup` emit records longest first, so a partial copy does not represent a lineage whose complete copy comes later.

**Tests and docs.** Update `tests/test_sieve_integration.py:327-334` (expectations 760/950) and `tests/test_composites.py:277`; add a 0.87-length row that survives the default. Docs: `docs/sieve.md:40`, `:44`, the `auto_min_len` docstring, the wrapper help.

### GX-04 · medium · Walker wrapper forces `--primary-metric cov` and cannot select 0.1.1's `matched`

`tools-iuc#8800/disassembler_walker.xml:81` · galaxy · confirmed (reproduced) · fix effort: trivial · sources: IC-01, G1-03, RB-04

**Problem.** The wrapper always renders `--primary-metric $anchoring.primary_metric` (:81) from a select of `cov` (selected, :128), `pident` and `cov_x_pident`; there is no `matched`. That matches the pinned 0.1.0, whose default was `cov` and which rejects `matched`. Walker 0.1.1 made `PrimaryMetric::Matched` the default (`walker/src/config.rs:456`) because `cov` is an uncapped ratio that lets a reference shorter than the sample's gene win (`walker/src/config.rs:7-17`, `docs/walker.md:76-79`). A token-only bump, as `docs/packaging.md:66-67` prescribes, keeps every Galaxy multi-anchor job on `cov`, and users cannot choose `matched`.

**Impact.** Latent until the bump, which is why this is medium. After it, for multi-record references (e.g. `panels/sarscov2_spike_11strains.fa`), primaries, `anchor=` tags and VCF/`snvs.tsv` coordinates follow the shortest covering reference (`config.rs` cites 62% + 21% of a 20,000-accession run landing on the two shortest spike references under `cov`), and Galaxy disagrees with the CLI at the same version. Single-record references, including all PR tests, are unaffected. Traces: a `[warn] ... 1.031x the query` line and `settings.primary_metric` in the summary.

**Evidence.** One accession whose graph is a 300-nt gene; panel `exact` (the gene) plus `short` (the gene minus 9 nt in frame); walker 0.1.1 offline. CLI default: `anchor=exact cov=1.000`, `"primary_metric": "matched"`. Wrapper flags: `anchor=short cov=1.031`. A 0.1.0 build exits 1 on `--primary-metric matched`.

**Proposed fix.** In the same commit as `@WALKER_VERSION@` 0.1.1, not earlier:
```xml
-  <option value="cov" selected="true">Reference coverage (cov)</option>
+  <option value="matched" selected="true">Fraction of the reference's bases matched (matched)</option>
+  <option value="cov">Reference coverage (cov)</option>
```
Help: "`cov` was the default before 0.1.1 and favours shorter references; use it only to reproduce old runs." Add a release step that diffs wrapper defaults against `logan-walker --help` (only `primary_metric` drifted from 0.1.0 to 0.1.1). Workflows saved against the 0.1.0 tool keep `cov` after an upgrade.

**Tests and docs.** `test-data/hiv1_pol_panel.fasta` = `hiv1_pol` plus an in-frame 9-nt-deletion copy, with SRR8209082; assert `^>SRR8209082_primary .* anchor=hiv1_pol$` and `"primary_metric": "matched"`. Add `<assert_command><has_text text="--primary-metric matched"/></assert_command>` to walker test 1. Docs: `docs/packaging.md:66-67`; `docs/galaxy.md:8-9` ("so the two stay in sync" does not hold for defaults).

### GX-05 · medium · Walker wrapper lowers `--max-size-mb` to fit job memory, so the cohort depends on the destination

`tools-iuc#8800/disassembler_walker.xml:34` · provenance · confirmed (reproduced) · fix effort: small · sources: G1-05

**Problem.** When one fetch and one compute worker do not fit 80% of `GALAXY_MEMORY_MB`, the wrapper rewrites the cap, `max_size_mb=\$(( (budget - 43) / 6 ));` (:34). The cap is not a tool parameter, the help says only that threads and memory follow the job's resources, and the notice is a stderr line (:42). The summary does record the effective cap (`settings.max_size_mb`) and each skip as `"status": "too_large"`, but `too_large` is not an error, so the job is green and the FASTA simply has fewer records.

**Impact.** Same inputs and parameters, different cohort per destination. With the default 1024 MB cap, any allocation below about 7.7 GB (`(6 x 1024 + 43) / 0.8` MB) lowers it: to 265 MB at 2 GB and 538 MB at 4 GB. In `bench/pfcrt_50.txt`, 17 of 50 accessions exceed 265 MB and 5 exceed 538 MB, so MSA and tree silently lose up to a third of the cohort.

**Evidence.** A small accession plus an 80 MB one, run through the wrapper's shell block and walker 0.1.1 offline: `GALAXY_MEMORY_MB=600` printed `trimmed to 1 compute and 1 fetch worker(s) with a 72 MB payload cap`, exit 0, the large accession `too_large`; `GALAXY_MEMORY_MB=16384` processed it.

**Proposed fix.**
1. `<param name="max_size_mb" argument="--max-size-mb" type="integer" value="1024" min="64" label="Skip accessions whose Logan payload exceeds (MB)"/>` in an advanced section; `max_size_mb=$max_size_mb` at :23.
2. A boolean `fit_cap` (default false, "Lower the cap to fit the job's memory instead of failing"), used in the :32-35 branch:
```
#if not $fit_cap:
  echo "a \${max_size_mb} MB cap needs at least \$(( (6*max_size_mb+43)*5/4 )) MB of job memory" >&2; exit 1;
#end if
max_size_mb=\$(( (budget - 43) / 6 ));
```
3. Document the cap, the memory rule and the `too_large` status in the help. For 0.2, a `skipped.tsv` output (backlog item 1(c)). Apply together with GX-08 (memory model) and GX-09 (admin ceiling).

**Tests and docs.** A tool test with `max_size_mb=64` asserting `"max_size_mb": 64` in the summary; a shell unit test of the block at 600 MB (exit 1, "needs at least 7733 MB") and with `fit_cap`. `docs/galaxy.md:40-41`: the wrapper also cuts `-w` and the cap, not only `--fetch-workers`.

### GX-06 · medium · Fetch-then-offline runs report no-data and over-cap accessions as errors

`walker/src/stages/stream.rs:110` · galaxy · confirmed (reproduced) · fix effort: small · sources: G2-04

**Problem.** `--fetch-only` classifies an accession with no Logan data, or over `--max-size-mb`, as a negative (`walker/src/pipeline.rs:226-238`) but leaves nothing in the cache to say so. A later `--offline` run sees a missing file and returns `Error::Fetch { ... "cache miss under --offline; expected {}" }` (`stream.rs:110-117`), which `check_run` (`walker/src/main.rs:401`) counts as an error that `--allow-empty` does not cover. The documented two-step flow (`docs/walker.md:161-168`) and the plan's WF3-frozen variant therefore do not reproduce a live run: `no_logan_data` and `too_large` become `error`.

**Impact.** With one accession per job, a live `--allow-empty` job exits 0 where the frozen job exits 3 and fails. With mixed lists, `summary.json` says `error` where the truth is `no_logan_data`. `fetch` without `--allow-empty` exits 1 on such input, and the documented example omits that flag.

**Evidence.** Walker 0.1.1, a cache holding one fabricated accession, plus `SRR0000001` (404 on S3), stand-in GraphAligner: live with `--allow-empty` exits 0 with `[no_path, no_logan_data]`; `--offline --allow-empty` exits 1 with `[no_path, error]`; `SRR0000001` alone offline exits 3, "cache miss under --offline".

**Proposed fix.** Record negatives in the cache and honour them offline:
```rust
// stream.rs
pub fn absent_path(dir, acc) -> PathBuf   // dir/"{acc}.unitigs.absent"
pub fn record_absent(dir, acc, e)         // writes "no_logan_data" or "too_large\t{detail}"
// Fetcher::fetch, offline cache-miss branch, before Error::Fetch:
//   marker "no_logan_data" -> Err(NoLoganData); "too_large" -> Err(TooLarge { detail })
// write_cache_atomic: remove a stale .absent after a successful fetch
// pipeline::run_fetch_only: on a negative, record_absent(cache, acc, &e)
```
With this patch both offline runs match the live statuses and exit codes. Store the remote size in the `too_large` marker and re-check it against the current cap. For Galaxy, have the fetch step also emit `fetch_status.tsv` (accession, status, size, ETag).

**Tests and docs.** `stream.rs` unit tests for both markers; a `tests/cli.rs` case: `--fetch-only --allow-empty`, then `--offline --allow-empty` on a no-data accession, exits 0. Docs: `docs/walker.md:153`, `:161-168` (add `--allow-empty`, describe the markers), the `main.rs` help text, the plan's WF3-frozen paragraph.

### GX-07 · medium · `harvest` silently drops kmindex `json_vec` scores and reports an empty harvest as success

`src/disassembler/cmd/harvest.py:162` · correctness · confirmed (reproduced) · fix effort: small · sources: G1-11

**Problem.** `_merge_payload` converts each hit with `sc = float(score)` (:162) and on any exception does `continue` (:163-164), with no counter or message. kmindex `--format json_vec` (an option of the Galaxy `kmindex_query` 0.7.1+galaxy0 tool, output datatype `json` like the default) writes each hit as `{"R": ratio, "P": [positions]}`. Every value fails `float()`, and harvest exits 0 with zero targets. Separately, the PR's harvest help (`disassembler_harvest.xml:27`, `:115-117`) describes only a per-shard collection, while `kmindex_query` defaults to `merge=true` and emits one merged JSON; that file parses correctly, so this half is a doc gap.

**Impact.** A user who chose json_vec gets `Harvest complete: 0 unique SRA runs across 0 targets`, empty collections and a union file with one blank line; downstream steps run on nothing and the green job hides why.

**Evidence.** A shard in the json_vec layout: `disassembler harvest --shards-dir vec/shards --out-dir cand` exits 0 with `summary.json` `{"_union": {"n": 0}}`; the same hits as plain json give 2 targets and 5 runs.

**Proposed fix.**
```python
if isinstance(score, dict) and "R" in score:
    score = score["R"]
try:
    sc = float(score)
except (TypeError, ValueError):
    n_unparsed += 1
    continue
```
Propagate `n_unparsed`, warn in `_write_outputs`, and in `harvest_local`/`harvest_collection` raise when `n_unparsed > 0` and no numeric score was merged ("hit values not numeric; expected kmindex --format json or json_vec"). Do not fail on genuinely empty shards (`{"idx": {"q": {}}}`). With the unwrap, the json_vec shard yields the same lists as json and `tests/test_harvest.py` passes.

**Tests and docs.** A json_vec shard gives the same lists as its json twin; a shard with only non-numeric values raises. Wrapper: label "kmindex query results (json or json_vec)", help naming the merged dataset, a test with a merged input. Docs: `docs/galaxy.md:17`, the harvest input description in `docs/kmindex.md`.

### GX-08 · medium · Walker wrapper's memory fit counts 43 MB per graph; chosen settings exceed the allocation

`tools-iuc#8800/disassembler_walker.xml:26` · galaxy · confirmed (reproduced) · fix effort: trivial · sources: G1-04, IC-05, RB-05

**Problem.** The wrapper fits `(4 + fetch + compute) x cap + 43 MB x compute` into 80% of `GALAXY_MEMORY_MB`: `fit=\$(( (budget - 43 * workers) / max_size_mb - 4 - workers ))` (:26), then `/ (max_size_mb + 43)` (:30) and `(budget - 43) / 6` (:34). The walker's own bound counts `GRAPH_PEAK_PER_COMPRESSED_BYTE = 1.7` times the cap per parsed graph (`walker/src/config.rs:445`, used at :424-434). The 43 MB was one 33 MB accession's array sizes; a 348 MB accession's graph peaked at 498 MB. `stages/parse.rs` has changed only in comments since v0.1.0, so this also holds for the pinned 0.1.0.

**Impact.** 4 slots/16 GB: the wrapper keeps 4 compute, 4 fetch, 1024 MB, and the walker reports an 18.8 GB worst case (118% of the allocation); 8/32 GB: 105%; 1/4 GB and 2/8 GB: 96-100%, so the reserve for GraphAligner is gone. The worst case needs several near-cap accessions in flight; the scheduler then kills the job (a visible failure). An actual OOM was not reproduced.

**Evidence.** The PR's shell block was run for each slots/memory pair and its choice passed to walker 0.1.1, which prints its bound: `slots=4 mem=16384 -> W=4 F=4 cap=1024`, bound 18.8 GB; `slots=8 mem=32768 -> W=8 F=8 cap=1024`, bound 33.6 GB.

**Proposed fix.** Count 2.7 x cap per compute worker (one buffer plus 1.7 for the graph), in tenths of a cap to stay in shell integers, and re-grow fetch workers after trimming compute:
```sh
units=\$(( budget * 10 / max_size_mb ));
fit=\$(( (units - 40 - 27 * workers) / 10 ));                  # was :26
...
fit=\$(( (units - 50) / 27 ));                                  # was :30
if [ "\$fit" -lt "\$workers" ]; then workers=\$fit; fi;
if [ "\$workers" -ge 1 ]; then fit=\$(( (units - 40 - 27 * workers) / 10 ));
  if [ "\$fit" -gt 1 ]; then fetch_workers=\$(( fit < workers ? fit : workers )); fi; fi;
...
max_size_mb=\$(( budget * 10 / 77 ));                           # was :34
```
Every tested pair then stays within 80% (4/16 GB gives 2/2/1024 at 71%; 8/32 GB gives 7/2/1024 at 78%). Long term, a walker `--memory-mb N` option should derive workers and cap from `Config::max_resident_bytes`, so no wrapper (or UDT, see GX-14) copies the model.

**Tests and docs.** A shell unit test with `GALAXY_SLOTS=4 GALAXY_MEMORY_MB=16384` asserting the "trimmed to" line. Rewrite the comment at :9-20 to cite `Config::max_resident_bytes` and `GRAPH_PEAK_PER_COMPRESSED_BYTE`; the :38 floor ("roughly 500 MB") becomes about 620 MB. `docs/walker.md:17` "about 43 MB" describes array size, not peak.

### GX-09 · medium · `DISASSEMBLER_WALKER_MAX_SIZE_MB` cannot be overridden from job_conf

`tools-iuc#8800/disassembler_walker.xml:100` · galaxy · confirmed (by reading) · fix effort: trivial · sources: G1-06

**Problem.** The cap is set with `<environment_variable name="DISASSEMBLER_WALKER_MAX_SIZE_MB">1024</environment_variable>` (:100), with the comment "admins can override per destination in job_conf.yml" (:97-99). Galaxy writes destination `env` entries first and the tool's variables after them (release_25.0 `lib/galaxy/jobs/runners/__init__.py:500-501`: `envs = destination.get("env", []); envs.extend(job_wrapper.environment_variables)`; the same order in release_24.1 and dev), so the tool's 1024 always wins. The error message at :38 gives the same ineffective advice.

**Impact.** An admin who sets the variable to 4096 on a large-memory destination still gets `--max-size-mb 1024`; accessions over 1 GB (2 of 50 in `bench/pfcrt_50.txt`) stay `too_large`. Only the undocumented `docker_env_<NAME>` and `singularity_env_<NAME>` destination parameters override it, and nothing does on conda destinations.

**Evidence.** Read in Galaxy release_24.1, release_25.0 and dev. Replaying the emitted order (destination value 4096, then the tool's value) leaves `max_size_mb=1024`.

**Proposed fix.** Remove the `<environment_variables>` block. Minimal version: `max_size_mb=\${DISASSEMBLER_WALKER_MAX_SIZE_MB:-1024}` at :23 and the same default at :41. Together with GX-05's parameter, treat the variable as an optional admin ceiling:
```sh
max_size_mb=$max_size_mb &&
if [ -n "\${DISASSEMBLER_WALKER_MAX_SIZE_MB:-}" ] && [ "\$DISASSEMBLER_WALKER_MAX_SIZE_MB" -lt "\$max_size_mb" ]; then
    echo "server ceiling: lowering --max-size-mb to \$DISASSEMBLER_WALKER_MAX_SIZE_MB" >&2;
    max_size_mb=\$DISASSEMBLER_WALKER_MAX_SIZE_MB; fi &&
```
Fix the comment at :97-99 and the message at :38.

**Tests and docs.** A shell unit test of the command block with the variable preset, expecting it to be honoured. Admin notes in `docs/galaxy.md`.

### GX-10 · medium · The 0.1.1 bump is not a token change: three PR tests fail, and galaxy.md already documents 0.1.1

`docs/packaging.md:66` · galaxy · confirmed (reproduced) · fix effort: small · sources: D2-03, G1-02, G1-08, BL1-12, IC-09, RB-06

**Problem.** The wrappers pin `@TOOL_VERSION@` and `@WALKER_VERSION@` to 0.1.0 (`tools-iuc#8800/macros.xml:2-3`), and `docs/packaging.md:66-67` says the Galaxy release step is only "Update the version tokens". At 0.1.1 that bump fails three PR tests:
- harvest test 1: harvest now writes `{"target": ..., "scores": ...}` (`src/disassembler/cmd/harvest.py:186`), while the goldens are bare maps (`disassembler_harvest.xml:58-59`);
- stratify test 1: the same bare maps are refused (`2 score file(s) have no target inside`, exit 1; see GX-11);
- walker test 3: a 404 is now `no_logan_data` (`walker/src/pipeline.rs:260`), so `<has_text text="HTTP 404"/>` (`disassembler_walker.xml:255`) fails.

Help and comments go stale too: the `--allow-empty` help (:122) implies success whenever nothing reconstructs, but at 0.1.1 a job whose accessions all errored exits 3; :165 says "fails at the fetch (HTTP 404)"; :67-69 name `/dev/shm/logan-walker` as the default. Meanwhile `docs/galaxy.md:46-50` already describes 0.1.1 ("a job in which every accession errored fails with exit code 3 regardless"), which is false for the pinned 0.1.0: there `--allow-empty` exits 0 when every accession errored (v0.1.0 `walker/src/main.rs:454`).

**Impact.** No wrong output reaches users from the test failures: CI goes red at the bump and blocks delivery. The forced `cov` (GX-04) and the memory model (GX-08) would survive a token bump. Until then the repository docs misdescribe the wrapper as submitted. The PR is unmerged, so no users are affected yet.

**Evidence.** With 0.1.1: harvest over the PR's three test shards differs from `test-data/hiv1_pol.scores.json` by the `"target"`/`"scores"` wrapper; stratify on the PR test data exits 1 as quoted; `logan-walker --accession SRR0000001 ... --allow-empty` exits 0 with `"status": "no_logan_data"` and no "HTTP 404". A v0.1.0 build with every accession errored and `--allow-empty` exits 0 with a 0-byte FASTA; 0.1.1 exits 3.

**Proposed fix.** In the PR, at the bump:
1. Regenerate `test-data/hiv1_pol.scores.json` and `pfcrt.scores.json` with 0.1.1 harvest; stratify test 1 then yields the expected `bins_summary`.
2. `disassembler_walker.xml:255`: `<has_text text="no_logan_data"/>`; :165: "has no Logan unitigs (status no_logan_data)".
3. :122 help: "Succeed when nothing was found (no Logan data, over the size cap, no anchors, no path). Accessions that error still fail the job." Update :67-69.
4. Apply GX-04, GX-08 and GX-11 in the same commit; re-run walker tests 1-2 under GraphAligner 1.0.17b.

In the repository: `docs/packaging.md:66-67` becomes "Update the tokens, run `planemo test` against the new packages, and update parameters, defaults, help and test data to match; diff wrapper defaults against `logan-walker --help`." `docs/galaxy.md:3-6`: note that the PR head still pins 0.1.0. `docs/galaxy.md:46-50`: prefix "From logan-walker 0.1.1" and say that 0.1.0 exits 0, so summaries must be checked for `"status": "error"`.

**Tests and docs.** The edited harvest test 1, stratify test 1 and walker test 3 pin the bump; add `<assert_command><has_text text="--primary-metric matched"/></assert_command>` to walker test 1.

### GX-11 · medium · Stratify wrapper names targets after datasets, which 0.1.1 no longer honours

`tools-iuc#8800/disassembler_stratify.xml:38` · galaxy · confirmed (reproduced) · fix effort: small · sources: IC-02, G1-12, PB-04

**Problem.** The wrapper stages each dataset as `candidates/<sanitised, de-duplicated dataset name>.scores.json` (`ln -s '$scores_file' 'candidates/${target}.scores.json'`, :38) and says "Each dataset's name is used as the target name" (:47). That held at 0.1.0, which took the target from the file stem. 0.1.1 recognises score files by content and takes the target from inside (`src/disassembler/cmd/stratify.py:165-169`, `docs/kmindex.md:67-73`); bare 0.1.0 files found by a directory scan are refused (`stratify.py:238-241`), and the staged name is ignored.

**Impact.** After the bump: (a) the PR's test data and every 0.1.0 score dataset in users' histories fail with exit 1, and the wrapper has no way to pass `NAME=`; (b) a renamed 0.1.1 dataset yields outputs named after the inner target, contrary to the help; (c) two harvests of one gene fail with "target ... is in two score files". (a) and (c) are loud; (b) misleads.

**Evidence.** With the 0.1.1 CLI and the wrapper's staging: the PR test data gives `stratify: 2 score file(s) have no target inside (written before 0.1.1, or by hand)`, exit 1; a 0.1.1 file staged as `my_gene.scores.json` writes `hiv1_pol.bins.txt`; the same file staged twice exits 1.

**Proposed fix.** CLI, `find_score_files`: let `NAME=` name only a file that has no target inside (what its help already says):
```diff
-                if sf.target not in (None, name):
-                    raise ValueError(f"{src}: the file's own target is {sf.target!r}")
-                sf = sf._replace(target=name)
+                if sf.target is None:
+                    sf = sf._replace(target=name)
+                elif sf.target != name:
+                    print(f"[warn] {src}: the file names its own target {sf.target!r}; "
+                          f"using that, not {name!r}", file=sys.stderr)
```
Wrapper: drop the symlink staging and the uniquifying loop; pass `--scores '${name}=${scores_file}'` per dataset with `--out-dir candidates` (the sanitiser already removes `=`). With both, the PR test data gives exactly the expected lists, and a renamed 0.1.1 file warns and proceeds. Falling back to file stems in directory scans is not recommended: it contradicts the documented design.

**Tests and docs.** `tests/test_stratify.py`: `--scores other=<0.1.1 file>` exits 0, warns, and names outputs after the inner target. Add a 0.1.1-shape wrapper test. Update the `--scores` help in `stratify.py`, `docs/kmindex.md:71-72`, and the wrapper help at :47 and :147: "the target is read from each file; older files take the dataset name; two datasets of one target are refused".

### GX-12 · medium · `@WALKER_VERSION@` token defeats `planemo autoupdate`

`tools-iuc#8800/macros.xml:3` · galaxy · confirmed (by reading) · fix effort: small · sources: SWG-02

**Problem.** The walker's tool version is `@WALKER_VERSION@+galaxy@VERSION_SUFFIX@` (`disassembler_walker.xml:1`), and its requirement `version="@WALKER_VERSION@"` (`macros.xml:13`) sits in the macros file every tool imports. planemo autoupdate (`planemo/autoupdate.py`, master) skips any tool whose version does not contain `@TOOL_VERSION@` (:265-270) and treats only `version="@TOOL_VERSION@"` as the main requirement (:93); for other requirements it compares the version text with the newest conda version (:232-239). The text `@WALKER_VERSION@` never equals a version, so the first autoupdate of any disassembler tool rewrites the shared requirement to a literal `0.1.1`.

**Impact.** After merge, an autoupdate bumps `@TOOL_VERSION@` to 0.1.1 and installs logan-walker 0.1.1 under a walker tool that still says `0.1.0+galaxy0`, a provenance hole users cannot see. Later walker releases never move the walker tool's version. That the tools-iuc bot runs autoupdate on every tool directory was not checked.

**Evidence.** A re-implementation of `autoupdate.py:93`, `:232-239` and `:265-270`, run against the PR directory with current Bioconda versions: the walker gives "Tool version does not contain @TOOL_VERSION@" (skipped); the logan-walker requirement text `@WALKER_VERSION@` against conda `0.1.1` is rewritten to the literal. planemo itself was not run.

**Proposed fix.** One `@TOOL_VERSION@` per package, in a file only that package's tools import:
- `macros.xml`: `@PROFILE@`, help, citations, `reference_input`;
- `disassembler_macros.xml`: `@TOOL_VERSION@` (disassembler), `@VERSION_SUFFIX@`, `requirements`, `msa_requirements`, `version_command`;
- `walker_macros.xml`: `@TOOL_VERSION@` (logan-walker), `@VERSION_SUFFIX@`, `<requirement type="package" version="@TOOL_VERSION@">logan-walker</requirement>`; the walker tool uses `version="@TOOL_VERSION@+galaxy@VERSION_SUFFIX@"`.

No tool may import both package files. Alternative: move the walker to `tools/logan_walker/` with its own macros and `.shed.yml`. Bump both packages to 0.1.1 together (GX-10).

**Tests and docs.** Run `planemo autoupdate --dry-run tools/disassembler` and `planemo lint` in the PR's CI; every tool should be reported as handled.

### GX-13 · medium · Plan's interim sieve UDT (bare `l1_pipeline_serial`) returns 0 or 1 records with exit 0

`GALAXY_PLAN.md:125` · galaxy · confirmed (reproduced) · fix effort: trivial · sources: G2-06

**Problem.** The plan's sieve UDT wraps `l1_pipeline_serial` in `biocontainers/kmer-diversity-sieve:0.1.0` with an embedded dedup pre-pass, and later "wraps the binary directly" (:98-99, :192-193). `disassembler sieve` does more than dedup: it derives `--min-len` from the reference (`src/disassembler/cmd/sieve.py:33-52`), passes gene-scale `--threshold 0.005 --sketch 256 --warmup 2000 --max-n-frac 0.05` (:208-217), concatenates the `*.rank*.fasta` shards (:224-226), and exits 1 when nothing is promoted (:237-242). The binary's defaults are archive-scale (`long min_len=27000`, `extern/kmer-diversity-sieve/scripts/l1_pipeline.c:335`; `#define DEF_THRESHOLD 0.15`, :79), and it exits 0 on an empty result.

**Impact.** Any gene-scale cohort through a UDT built as the plan describes ends as a green job with 0 or 1 sequences; codon alignment or the tree then fails, or runs on one sequence.

**Evidence.** On the PR's 33-record `hiv1_pol_reconstructed.fasta`: `disassembler sieve ... --ref hiv1_pol.fasta` keeps 28 records; `l1_pipeline_serial --input IN --emit-promoted b` rejects all 33 at the `<27000 nt` gate, exit 0; with `--min-len 2602`, 1 record; with all five gene-scale options and the shards concatenated, the output is byte-identical to `disassembler sieve`.

**Proposed fix.** Plan change only. The interim UDT is `disassembler-sieve`, wrapping `disassembler sieve -i IN --ref REF -o sieved.fasta [--min-len N] [--max-gap-frac F] [--threshold T] [--sketch S]` in `biocontainers/disassembler:<ver>--pyhdfd78af_0`, which has Python and the sieve as a run dependency (so no dedup in awk). A future bare-binary UDT must set all five options, derive `--min-len` from the reference, concatenate the shards, and fail on an empty result (`[ -s out.fasta ] || exit 1`).

**Tests and docs.** A UDT test on `hiv1_pol_reconstructed.fasta` plus `hiv1_pol.fasta` asserting 28 records; an impossible `--min-len` expecting a non-zero exit. Plan lines 98-99, 125, 192-193.

### GX-14 · medium · Plan's walker UDT cannot reuse the IUC memory block; ported, it cuts the cap to 497 MB

`GALAXY_PLAN.md:122` · galaxy · confirmed (reproduced) · fix effort: small · sources: G2-01

**Problem.** The plan says the walker UDT "carries the GALAXY_SLOTS/MEMORY memory-fitting shell from the IUC XML". A UDT has no `<environment_variables>`, so `DISASSEMBLER_WALKER_MAX_SIZE_MB` is undefined, and the XML's Cheetah escapes (`\$workers`) survive into the shell: the block pasted verbatim builds, then fails at run time ("division by 0", or `--workers '$workers'` and a usage error). Ported correctly, the block still lowers the cap (GX-05) and still uses the 43 MB model (GX-08). With no resource requirement, usegalaxy.org runs a UDT at 1 core and 3788 MB, where the cap drops to 497 MB; at its UDT ceiling (8 cores, 29696 MB) the block keeps 8/8/1024, whose walker bound is 33.6 GB.

**Impact.** Pasted versions fail loudly. The ported version succeeds while skipping 11 of the 150 bench accessions as `too_large` instead of 5 (pfcrt_50: 8 instead of 2), visible only on stderr and in the summary.

**Evidence.** Galaxy's UDT interpolation was simulated (cwl-utils logic) and the result run in bash: the ported block prints `trimmed to 1 compute and 1 fetch worker(s) with a 497 MB payload cap`; `logan-walker --max-size-mb 497` skips DRR344521 (553 MB) as `too_large`.

**Proposed fix.** Rewrite the plan row: declare `requirements: [{type: resource, cores_min: 4, ram_min: 24576}]`, make the cap an integer input (default 1024), fit workers with the walker's own model, and fail rather than lower the cap:
```sh
cap=$(inputs.max_size_mb) && workers=$GALAXY_SLOTS && [ -n "$workers" ] || workers=1; fetch_workers=$workers;
if [ -n "$GALAXY_MEMORY_MB" ]; then budget=\$(( GALAXY_MEMORY_MB * 4 / 5 ));
  need=\$(( cap * (4 + fetch_workers + workers) + cap * workers * 17 / 10 ));
  while [ "$need" -gt "$budget" ] && [ \$(( fetch_workers + workers )) -gt 2 ]; do
    if [ "$fetch_workers" -ge "$workers" ] && [ "$fetch_workers" -gt 1 ] || [ "$workers" -eq 1 ]; then fetch_workers=\$(( fetch_workers - 1 )); else workers=\$(( workers - 1 )); fi;
    need=\$(( cap * (4 + fetch_workers + workers) + cap * workers * 17 / 10 )); done;
  if [ "$need" -gt "$budget" ]; then echo "... needs at least \$(( need * 5 / 4 + 1 )) MB ..." 1>&2; exit 1; fi; fi;
```
Results: 1 core/3788 MB exits 1 ("needs at least 9856 MB"); 4 cores/24576 MB gives 4/4 (18.8 GiB); 8 cores/29696 MB gives 5/5 (22.5 GiB). Without `ram_min` every default job would exit 1, so the requirement is part of the fix. The durable fix is the walker `--memory-mb` option (GX-08).

**Tests and docs.** A UDT test with a small `GALAXY_MEMORY_MB` asserting the exit-1 message; for `--memory-mb`, a `config.rs` test that the chosen settings satisfy `max_resident_bytes <= budget`. Plan line 122.

### GX-15 · medium · Plan's workflows reference UDTs by id, which Galaxy cannot resolve; import caveat reversed

`GALAXY_PLAN.md:138` · galaxy · confirmed (by reading) · fix effort: small · sources: G2-05

**Problem.** The plan authors workflows with "UDT ids for disassembler steps" (:138-139) and says they "import cleanly but refuse at invoke on accounts that haven't registered them" (:187-188). Galaxy resolves a user-defined tool only by `tool_uuid` or by an embedded `GalaxyUserTool` definition: `ToolModule.__init__` looks up unprivileged tools by `tool_uuid` and otherwise falls back to `toolbox.get_tool(tool_id, ...)`, which holds no UDTs (release_26.2 `lib/galaxy/workflow/modules.py:2713-2718`). A bare UDT id resolves to nothing, even for an account that registered the tool. With an embedded definition the behaviour is the reverse of the caveat: import creates a private copy (`resolve_tool_uuid`, :2744-2757) and is refused for users without the Custom Tool Execution role (`lib/galaxy/managers/tools.py:97-101`, `:207-208`).

**Impact.** Workflows written per section 3 cannot be tested by planemo in iwc-lab CI, and Dockstore users get missing-tool steps. The cost is rework of the delivery format, not wrong output.

**Evidence.** Read in the Galaxy 26.1 and 26.2 sources, gxformat2 (round-trips `run: {class: GalaxyUserTool}`), and planemo (`ensure_user_tool_execute` for workflows with embedded user tools). Not run: no Galaxy server was available.

**Proposed fix.** Plan text:
1. Section 3: keep `udt/*.gxtool.yml` as the reviewed source; generate each workflow with the definitions inlined (`run: {class: GalaxyUserTool, ...}`) by a small script; commit the generated `.gxwf.yml` with its tests, `.dockstore.yml`, README and CHANGELOG.
2. Section 4 caveat: "Workflows embed the UDT definitions. Importing creates private copies for users with the Custom Tool Execution role on a server with `enable_beta_tool_formats`; other users are refused at import."
3. Confirm that the planemo used by planemo-ci-action supports embedded user tools.

**Tests and docs.** A planemo `workflow_test` of WF3 on Galaxy 26.1 or later with beta tool formats enabled. Plan lines 136-139 and 186-191.

### GX-16 · low · Walker wrapper loses `walker_summary.json` in failed jobs

`tools-iuc#8800/disassembler_walker.xml:91` · galaxy · confirmed (reproduced) · fix effort: trivial · sources: SWG-03

**Problem.** The command ends `$output_options.emit_snv_table &&` (:91), then `mv walker_out/*_walker_summary.json walker_summary.json` (:94). The walker writes its summary before `check_run`, but a non-zero exit skips the `mv`, so the Summary dataset is empty in exactly the failed jobs whose per-accession status and error text would explain the failure. An offline cache miss reproduced it: exit 3, no `walker_summary.json`.

**Proposed fix.**
```diff
-            $output_options.emit_snv_table &&
-        mv walker_out/*_walker_summary.json walker_summary.json
+            $output_options.emit_snv_table ;
+        rc=\$? ;
+        find walker_out -type f -name '*_walker_summary.json' -exec mv {} walker_summary.json ';' ;
+        exit \$rc
```

**Tests and docs.** Add `<assert_stderr>` to the `expect_failure` test; Summary help: "written for failed runs too".

### GX-17 · low · galaxy.md offers `--allow-empty` as an alternative to filtering failed jobs

`docs/galaxy.md:80` · docs · confirmed (reproduced) · fix effort: trivial · sources: G2-08

**Problem.** "With one accession per job, either set `--allow-empty` or filter failed datasets before collapsing." Since 0.1.1, `check_run` returns `AllErrored` before it looks at `allow_empty` (`walker/src/main.rs:408-409`), so an accession that errors (fetch failure, GraphAligner crash) fails its job regardless, and one red element pauses the collapse for the whole cohort. Lines 46-50 were updated for 0.1.1; line 80 was not. Reproduced with an unreachable proxy and `--allow-empty`: exit 3.

**Proposed fix.**
```diff
-or MSA. With one accession per job, either set `--allow-empty` or filter failed datasets before
-collapsing.
+or MSA. With one accession per job, filter failed datasets before collapsing: an accession that
+errored still fails its job (exit 3). Setting `--allow-empty` as well keeps accessions where
+nothing was found as empty, succeeded datasets.
```

**Tests and docs.** Doc only; see GX-10 for the wrapper's `--allow-empty` help.

### GX-18 · low · With `--require-full-length`, log and wrapper help still say the partial tier is 50-80%

`src/disassembler/cmd/lexicmap.py:955` · docs · confirmed (reproduced) · fix effort: small · sources: SWG-04

**Problem.** The full-length gate is `full_cov = 1.0 if full_length_only else min_coverage` (:747), so under `--require-full-length` partials span 50% to under 100%. The log still prints `full_length_cov={min_coverage*100:.1f}%` (:549) and `Clean Partials ({min_coverage_partial:.0%}-{min_coverage:.0%} cov)` (:955); docstrings (:22, :531) and the PR wrapper ("the 50-80 percent partial tier", `disassembler_lexicmap.xml:47`, `:264-265`) say the same. `summary.json` records no gates. Data are correct; the cohort can be mislabelled in a write-up.

**Proposed fix.** Compute `full_cov` once before the banner; print `full_length_cov={full_cov*100:.1f}%` plus " (--require-full-length)" and `Clean Partials ({min_coverage_partial:.0%}-<{full_cov:.0%} cov)`; add `summary_data["gates"]` (full_length_cov, partial_cov, require_full_length, min_pident, max_internal_stops, allow_frameshift_alts). Wrapper: "from the partial gate up to the full-length gate: 50-80% by default, 50-100% with *Require 100% coverage*".

**Tests and docs.** A test with `full_length_only=True` and one 83% row asserting the new wording and `gates`; `docs/lexicmap.md:55`.

### GX-19 · low · msa wrapper stale at 0.1.1: `--msa-min-cov` mislabelled, no `--min-identity`

`tools-iuc#8800/disassembler_msa.xml:30` · wrapper-drift · confirmed (reproduced) · fix effort: small · sources: G1-14

**Problem.** At 0.1.1 the gate counts A/C/G/T only (`if real_bases(seq) < ref_len * args.msa_min_cov:`, `src/disassembler/cmd/msa.py:242`) and `--min-identity` (default 0.60, :166) fails runs with low median identity. The wrapper, correct for 0.1.0, labels the gate "Minimum length relative to the reference" and lacks `--min-identity`, which the CLI's error message tells users to adjust. Reproduced: a reference-length row with 50% `-` is dropped at 0.70.

**Proposed fix.** At the bump:
```xml
<param argument="--msa-min-cov" ... label="Minimum coverage of the reference"
       help="A/C/G/T bases as a fraction of the reference length; '-' and N do not count."/>
<param argument="--min-identity" type="float" value="0.60" min="0" max="1"
       label="Minimum median identity to the reference" help="0 disables the check."/>
```
plus `--min-identity $min_identity` in the command. For 0.2 (backlog item d): `msa --prep-only` and an optional `--out-tree`, so a workflow can use stock MAFFT and FastTree tools.

**Tests and docs.** A wrapper test with `min_identity=0` on divergent input; `docs/codon_msa.md` option table.

### GX-20 · low · `--allow-frameshift-alts` help describes 0.1.0 semantics

`tools-iuc#8800/disassembler_lexicmap.xml:41` · docs · confirmed (reproduced) · fix effort: trivial · sources: G1-10

**Problem.** The help says "Let accessions whose indels or frameshifts produce premature stops into the clean cohort anyway", true at 0.1.0. At 0.1.1 the flag only waives the new frameshift tier (`frame_ok = not frameshifts or allow_frameshifts`, `src/disassembler/cmd/lexicmap.py:749`); stops still flag, and the tier descriptions (:48, :219-220, :258-259) never mention frameshifts. With the PR test data and the flag, SRR003 stays `FLAGGED_PREMATURE_STOPS` at 0.1.1 but becomes clean at 0.1.0.

**Proposed fix.** At the bump:
```
help="Let accessions whose tiled sequence carries a frameshift (reading frame shifted for 30 or
      more reference bases) into the clean tiers instead of flagging them. Premature stops still
      flag; use 'Maximum premature stop codons' to tolerate those."
```
List "flagged (premature stops, frameshift, or low coverage)" at :48, :219-220, :258-259.

**Tests and docs.** A test with `allow_frameshift_alts=true` asserting `flagged_table` still contains SRR003.

### GX-21 · low · msa/lexicmap help "gap-free, only adds the tree" becomes false at 0.1.1

`tools-iuc#8800/disassembler_msa.xml:90` · docs · confirmed (reproduced) · fix effort: trivial · sources: G1-13, BL2-13

**Problem.** The msa help (:90-93) and lexicmap help (`disassembler_lexicmap.xml:246-247`) say clean rows under *Require 100% coverage* are gap-free and `msa` "only adds the tree". At 0.1.1 trusted in-frame deletions count as covered (`src/disassembler/cmd/lexicmap.py:420`), so clean rows can carry `---` and MAFFT may move them. Reproduced: a deletion in a GAT repeat placed at column 154 by lexicmap was moved to column 145 by `msa`.

**Proposed fix.** In both files:
```
With *Require 100% coverage ...* on, the only gaps in clean rows are in-frame deletions;
rows without a deletion come through MAFFT unchanged, and a deletion in a repeat may be
placed at a different column.
```

**Tests and docs.** Optional wrapper test on a deletion fixture. The repository's own docs are already correct.

### GX-22 · low · Stratify N help says "evenly across the score range"; sampling is by rank

`tools-iuc#8800/disassembler_stratify.xml:51` · docs · confirmed (reproduced) · fix effort: trivial · sources: SWG-05

**Problem.** `stratified_sample` takes every `len(items)/k`-th of the score-sorted candidates (`src/disassembler/cmd/stratify.py:57`), so picks follow candidate density; the wrapper's long help (:137-138) says so, the parameter help contradicts it. Reproduced: one bin `0.0:1.01:4` over 4 strong and 96 weak hits returned 1 strong and 3 weak.

**Proposed fix.** Help only:
```
help="Taken at even steps through the bin's candidates ranked by score (every len/N-th, starting
      with the best), so densely populated score ranges contribute more; use narrower bins to
      control the score distribution. A bin with fewer candidates contributes what it has."
```

**Tests and docs.** A unit test pinning indices 0, 25, 50, 75 on the skewed list; reword `docs/kmindex.md:55`.

### GX-23 · low · MAFFT: galaxy.md calls it an IUC tool; the plan says no 7.526 container exists

`docs/galaxy.md:86` · docs-galaxy · confirmed (reproduced) · fix effort: trivial · sources: BL2-06, PK2-12

**Problem.** galaxy.md says "MAFFT and FastTree already have well-maintained IUC tools". There is no IUC MAFFT; the maintained one is `rnateam/mafft` `rbc_mafft_add` 7.526+galaxy2 (MAFFT 7.526). The plan (`GALAXY_PLAN.md:126`, `:194`) says "no 7.526 biocontainer exists", but `quay.io/biocontainers/mulled-v2-027f30e88c6da4ab4b0170ff590dbab158b07243:43ff19c348e21674cd51004463f1ea0b2e5f2336-0` (mafft 7.526 + fasta3 36.3.8) does. IUC `fasttree` is 2.1.10 and defaults to protein: on a codon MSA it builds a JTT tree with only a stderr warning. `rbc_mafft_add`'s default command gave output byte-identical to `msa`'s MAFFT step on a 22-row fixture.

**Proposed fix.** galaxy.md:85-89: name `toolshed.g2.bx.psu.edu/repos/rnateam/mafft/rbc_mafft_add/7.526+galaxy2` and IUC `fasttree` with "set Nucleotide and GTR+CAT". Plan: drop the mafft UDT row and its fallbacks; note that `msa` passes `--auto` and `rbc_mafft_add` does not.

**Tests and docs.** Doc only; a future split workflow should assert FastTree runs with `-nt`.

### GX-24 · low · No native Galaxy step splits an accession list; galaxy.md misnames tool owners

`docs/galaxy.md:79` · docs · confirmed (by reading) · fix effort: trivial · sources: G2-07

**Problem.** The plan relies on "native collection ops for the accession scatter" (`GALAXY_PLAN.md:137-138`), but Galaxy 26.2's collection operations cannot turn lines of a text dataset into elements; a Tool Shed tool is needed. galaxy.md:71 calls it "Galaxy's **Split file** tool" (it is `bgruening/split_file_to_collection`), and :79 says "the IUC **Collapse Collection** tool" (it is `nml/collapse_collections`). The MAFFT owner is GX-23.

**Proposed fix.** Name both tools with owners. In the plan: "split step (bgruening/split_file_to_collection, chunk mode with *Maintain record order*), map-over, Filter failed datasets, nml Collapse Collection". The split tool defaults to round-robin, so name chunk mode to keep the collapsed order reproducible.

**Tests and docs.** Doc only.

### GX-25 · low · Plan caveats omit UDT destination limits

`GALAXY_PLAN.md:195` · galaxy · confirmed (reproduced) · fix effort: trivial · sources: G2-10

**Problem.** Galaxy's documented UDT destination sets `docker_net: "none"` (release_26.2 `doc/source/admin/user_defined_tools.md:75`, `:111`); without network the walker exits 3 with "Connection refused" fetch errors (reproduced with an unreachable proxy). The plan also omits usegalaxy.org's UDT destination limits (usegalaxy-playbook `tpv/environments.yaml.j2:413-428`): at most 8 cores and 29 GB, CVMFS-published BioContainers only. So a "mulled-build" MAFFT image works only once on CVMFS, and the walker's memory fit must target 29 GB.

**Proposed fix.** Add both limits to section 4, citing the files and the date; make a one-accession walker UDT run on usegalaxy.org the first step. Optionally, have the walker's exit-3 message suggest `--unitigs-cache`/`--offline` when every fetch failed to connect.

**Tests and docs.** Plan sections 2 and 4; `docs/galaxy.md:36-37` should name the `docker_net` setting.

### GX-26 · low · Plan's fetch UDT runs `disassembler` in the logan-walker container, which has none

`GALAXY_PLAN.md:123` · galaxy · confirmed (by reading) · fix effort: trivial · sources: PK2-05

**Problem.** The row reads `disassembler-fetch | disassembler fetch | same walker container`. The logan-walker Bioconda recipe's only run requirement is `graphaligner ==1.0.17b`, so the image has no `disassembler`; the UDT, and WF3-frozen's first step (:173), fail with `command not found`. `disassembler fetch` is a thin wrapper around `logan-walker --fetch-only` (`src/disassembler/cmd/fetch.py:58`).

**Proposed fix.** Wrap `logan-walker --fetch-only --unitigs-cache cache -l accessions.txt [--max-size-mb N] [--allow-empty] [--fetch-workers N]` in the walker image, or keep `disassembler fetch` in `biocontainers/disassembler:<ver>`, which includes logan-walker.

**Tests and docs.** A UDT test with a pre-populated cache (no network needed). Plan lines 123 and 173.

### GX-27 · low · WF1 sketch omits the reference CDS and target selection

`GALAXY_PLAN.md:144` · galaxy · confirmed (reproduced) · fix effort: trivial · sources: G2-03

**Problem.** WF1's only input is "kmindex shard JSONs" (:144), but the walker requires `--query` (`walker/src/main.rs:454`), and the composite it mirrors is single-target (`src/disassembler/cmd/pipe_from_kmindex.py:3-5`). harvest and stratify emit one file per target (a two-target shard gave one `<target>.bins.txt` each), so a per-target collection appears with no step to pick one. Any author notices the missing reference at once.

**Proposed fix.**
```text
inputs: kmindex shard JSONs (list), target name (text), reference CDS (fasta), score bins (text)
  -> disassembler-harvest -> __EXTRACT_DATASET__ (identifier = sanitised target name)
  -> disassembler-stratify (--target) -> split -> walker (mapped) -> filter failed -> collapse
  -> [optional] sieve -> codon-align -> fasttree
```
Better: make WF3 a subworkflow and have WF1 add only the head.

**Tests and docs.** A workflow test with a two-target shard; plan lines 141-150.

### GX-28 · low · WF3-frozen quota estimate (median x 50) understates skewed cohorts

`GALAXY_PLAN.md:178` · galaxy · confirmed (reproduced) · fix effort: trivial · sources: G2-09

**Problem.** "Median ~71 MB ... (~3.5 GB per 50-accession cohort)" multiplies the median by 50, but Logan sizes are heavily skewed (`docs/kmindex.md:82`). Measured with HEAD requests: `bench/pfcrt_50.txt` totals 29.6 GiB, 9.1 GiB under the default 1024 MB cap (2.6 times the estimate); `bench/hiv1_50.txt` 7.0 / 3.1 GiB; `bench/sarscov2_50.txt` 1.15 GiB. usegalaxy.org has no zstd datatype, so fetched files land as `binary`. The backlog (`BACKLOG_0_2.md:234`) endorses the figure; `docs/galaxy.md:27` says "median about 50 MB each".

**Proposed fix.** Replace with dated per-organism totals and the note that the cap bounds one object, not the cohort.

**Tests and docs.** Plan :178-180, backlog :234-235, `docs/galaxy.md:27`.
