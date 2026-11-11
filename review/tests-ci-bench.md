# Tests, CI and bench: Python tests, the CI workflow, the golden gate and the simulation (tests/, .github/workflows/ci.yml, bench/)

The every-push CI is sound in shape: clippy with `-D warnings`, offline walker tests, a three-version Python matrix that discovers tests instead of listing them, and a real-tool MSA job that treats a skip as a failure. The weak side is what the gates leave out. The only real-data gate (`bench/check_golden.sh`) runs only on manual dispatch, has not run in CI since 2026-09-12, and even then checks only the primary under non-default flags. Several modules (resolver, validators, `fetch`, the walker passthrough) and many cohort-shaping rules (sampling, merging, flag forwarding, sieve defaults) have no test, and mutation runs show that each can be broken with the suite still green. The k-mer simulation that justifies the default `-k 21` overstates its own results and cannot be rerun. Most fixes share one theme: make an existing expectation executable, by adding a pinning test, a `--locked`, a required check or a seed, rather than changing runtime code.

| ID | Severity | Where | Title | Fix effort |
|---|---|---|---|---|
| TC-01 | medium | `.github/workflows/ci.yml:101` | Golden gate CI job runs only on manual dispatch; no CI run for either release | small |
| TC-02 | medium | `bench/check_golden.sh:59` | No end-to-end regression check of the default walk (pruning, majority preference, alts, VCF) | medium |
| TC-03 | medium | `src/disassembler/resolve.py:47` | Resolver, validators, `fetch` and walker passthrough are untested; the gap hides a silent binary swap and three wrapper defects | small |
| TC-04 | medium | `bench/simulated_kmer_sensitivity/README.md:25` | The k=21 simulation misreports its own results and measures anchoring only | small |
| TC-05 | medium | `src/disassembler/cmd/pipe_reconstruct.py:357` | Forwarding of 12 walker options from the composites is untested | small |
| TC-06 | medium | `tests/test_sieve.sh:38` | Sieve gene-scale defaults (tau 0.005, sketch 256) are pinned by no test | small |
| TC-07 | medium | `tests/test_stratify.py:28` | Even-spread sampling, max-merge, zero-drop, 429 handling and nested collections are unpinned | small |
| TC-08 | low | `.github/workflows/ci.yml:5` | `main` is unprotected: CI is advisory and fork PRs were merged untested | trivial |
| TC-09 | low | `.github/workflows/ci.yml:34` | CI builds without `--locked` and checks no release identifiers; Bioconda builds with `--locked` | small |
| TC-10 | low | `bench/check_golden.sh:32` | Golden gate green-lights a stale prebuilt walker without rebuilding or naming it | trivial |
| TC-11 | low | `tests/test_msa_real_tools.py:75` | Real-tool MSA test passes with `--keeplength` removed | trivial |
| TC-12 | low | `tests/test_panels.py:40` | Reference-quality test skips `queries/` and FASTA extensions other than `.fa`/`.fasta` | small |
| TC-13 | low | `tests/test_sieve.sh:51` | Sieve tests skip green when the submodule is absent or the build fails | trivial |
| TC-14 | low | `tests/test_sieve_integration.py:35` | Sieve test bootstrap trusts any existing binary and calls `make` unguarded | small |
| TC-15 | low | `tests/test_dedup_sieve_equivalence.py:62` | No test pins `sieve_run`'s binary flags; the equivalence test hard-codes its own copy | small |
| TC-16 | low | `bench/simulated_kmer_sensitivity/harness.py:143` | Conserved-motif model protects 30-nt windows, shorter than a 31-mer | small |
| TC-17 | low | `bench/simulated_kmer_sensitivity/harness.py:488` | Simulation harness cannot regenerate its committed results | small |
| TC-18 | low | `bench/check_golden.sh:85` | Golden gate labels a walk that found nothing as "NOT a walk regression" | trivial |
| TC-19 | low | `bench/check_golden.sh:29` | Golden gate ignores a failed `mktemp` and writes its work files into `/` | trivial |
| TC-20 | low | `.github/workflows/ci.yml:23` | CI lints with floating `stable` Rust under `-D warnings`; cache key ignores the compiler | trivial |
| TC-21 | low | `.github/workflows/ci.yml:3` | CI hygiene: no dependency audit, no job timeouts or permissions, no tag-time version check | small |
| TC-22 | low | `.github/workflows/ci.yml:16` | CI runs only on Linux and has no coverage signal, while macOS packages are planned | small |
| TC-23 | low | `.github/workflows/ci.yml:62` | No Python lint in CI; the existing debt is 13 auto-fixable lines | trivial |
| TC-24 | low | `.github/workflows/ci.yml:33` | Rust sources are not rustfmt-clean and CI has no format check | small |
| TC-25 | low | `bench/golden/README.md:48` | `bench/golden/README.md` still documents the deleted archive, verifier and harness | small |
| TC-26 | low | `tests/test_composites.py:59` | Test scaffolding is duplicated across files; one fixture reader assumes single-line FASTA | small |
| TC-27 | info | `bench/golden/GRAPHALIGNER_VERSION:1` | Two GraphAligner 1.0.17b builds are in use and no record tells them apart | small |

### TC-01 · medium · Golden gate CI job runs only on manual dispatch; no CI run for either release

`.github/workflows/ci.yml:101` · ci · confirmed (reproduced) · fix effort: small · sources: T1-01, X1-01

**Problem.** The only CI job that exercises real reconstructions (`bench/check_golden.sh` plus the 17 `#[ignore]` data-backed walker tests) is guarded by `if: github.event_name == 'workflow_dispatch'` (`.github/workflows/ci.yml:101`), so every push and pull request reports it as skipped. It has been dispatched twice, both on 2026-09-12 (257401e and ff9e539), 63 commits before v0.1.1. Since then the fixtures were recaptured, GraphAligner was re-pinned to 1.0.17b (00e3b97), the data-backed step with its SKIP check was added (e36970e), and `walker/src` changed by +2004/−513 lines in 17 files. The job as it now stands has never executed. Release gating rests on step 2 of `docs/packaging.md:60-62`, a local run that leaves no record.

**Impact.** A walker change that shifts the primary on ERR2028820 or SRR12363134, or breaks a data-backed test, is merged and tagged with every check green. A broken job definition (micromamba solve, the `--fetch-only` step, the SKIP grep) is found only when someone first dispatches it. The comment at `ci.yml:119-122` already records three data-backed tests that stayed broken unnoticed until 0.1.1. A maintainer may have run the gate locally; what is missing is any CI record of it, not proof of bad output.

**Evidence.** `gh run list --workflow ci.yml --event workflow_dispatch` lists two runs (#4 and #6, 2026-09-12); `ci.yml` is the only workflow. `git rev-list --count ff9e539..25c27f2` prints 63; `git diff --stat f655edd HEAD -- walker/src` shows 17 files, +2004/−513.

**Proposed fix.** Move the job to its own workflow, `.github/workflows/golden.yml`, drop the `if:` and trigger it without a human. A separate file is needed because `paths` filters apply to a whole workflow:

```yaml
on:
  workflow_dispatch:
  schedule: [{cron: '23 4 * * 1'}]
  push:
    branches: [main]
    tags: ['v*', 'logan-walker-v*']   # no paths filter here: it would drop tag pushes
  pull_request:
    branches: [main]
    paths: ['walker/**', 'bench/golden/**', 'bench/check_golden.sh',
            'environment.yml', '.github/workflows/golden.yml']
```

Add `actions/cache` keyed on the accession list for `$RUNNER_TEMP/logan` (Logan objects are immutable) so S3 hiccups do not turn the gate red. Dispatch it once on 25c27f2 now. A tag-push run is a post-hoc check, so packaging step 2 should also require the URL of a green run for the release commit.

**Tests and docs.** Validate with one dispatch and one PR touching `walker/`. Update the `ci.yml:98-100` comment, `bench/README.md:19-20` and `docs/packaging.md:60-62`. GitHub disables schedules after 60 days without repository activity; the first run may fail on the untested 1.0.17b environment, which is the point. TC-02 covers what the gate checks; TC-10, TC-18 and TC-19 cover defects in the script itself.

### TC-02 · medium · No end-to-end regression check of the default walk (pruning, majority preference, alts, VCF)

`bench/check_golden.sh:59` · test-coverage · confirmed (by reading) · fix effort: medium · sources: X1-09

**Problem.** The only real-data gate runs the walker with `--prune-frac 0.0 --no-prefer-majority --no-emit-alts` (`bench/check_golden.sh:58-59`) and compares only the bare primary sequence (`:61`, `:91`). So the code that shapes every default primary (the prune/realign acceptance in `walker/src/pipeline.rs`, `stages::variants::prefer_majority`) and the alt-haplotype, VCF and SNV-table outputs never run end to end in any test. The data-backed `#[ignore]` tests call stage functions directly; the only one that touches pruning, `emits_gfa_for_cross_check` (`walker/src/stages/extract.rs:851`), asserts nothing and CI skips it. The offline `walker/tests/cli.rs` tests check exit codes and scratch handling, never output content. Each rule does have synthetic unit tests.

**Impact.** A change that alters default primaries through pruning or the majority swap, or that drops or duplicates alt haplotypes or VCF records, passes every gate. Every user runs the defaults. The script header says SRR8209082's primary "is stabilised by the default pruning, which this old-behaviour gate turns off on purpose" (`:24-25`), so the accession where the defaults matter most is the one left ungated.

**Evidence.** Read from the gate's flags and comparison, the list of `#[ignore]` tests and `walker/tests/cli.rs`. Not run: it needs GraphAligner 1.0.17b and S3 access.

**Proposed fix.**
1. Capture a default-flag fixture set under GraphAligner 1.0.17b (`-t 1`) with `--emit-vcf`: per accession, the primary FASTA, `alt_haplotypes.fasta`, the VCF body without date and command header lines, and the SNV table, under `bench/golden/default/`.
2. In `check_golden.sh` add a `DEFAULT_GATES` loop that runs the walker without the old-behaviour flags and `cmp`s each normalised file; keep the existing gate unchanged. Gate SRR8209082 only if two consecutive default runs reproduce byte for byte; otherwise report it as INFO.
3. Offline complement for every push: a `walker/tests/cli.rs` test with a fabricated unitig fixture and a stand-in GraphAligner that asserts primary, alt and VCF content under defaults.

**Tests and docs.** The fixtures and loop above. Update `bench/golden/README.md` ("What is asserted"), the script header and `docs/walker.md:284-289` to describe both gate sets and how to regenerate them. A default gate on a non-deterministic subgraph would flap, hence the reproducibility check first.

### TC-03 · medium · Resolver, validators, `fetch` and walker passthrough are untested; the gap hides a silent binary swap and three wrapper defects

`src/disassembler/resolve.py:47` · test-coverage · confirmed (reproduced) · fix effort: small · sources: T1-05, P4-14, BL2-07

**Problem.** No file in `tests/` imports `disassembler.resolve`, `disassembler.validators` or `disassembler.cmd.fetch`, nothing calls the `walker` passthrough, and nothing sets `DISASSEMBLER_WALKER_BIN`/`DISASSEMBLER_SIEVE_BIN` (backlog section 5 still holds). Behind the gap sit four defects:
- `resolve_walker` chains `_executable(explicit) or _executable(env) or _from_path(...) or <in-tree>` (`resolve.py:47-52`; `resolve_sieve` at `:57-63`), so a named path that is missing or not executable silently falls through to another binary.
- `positive_float` checks only `if v <= 0.0:` (`validators.py:67`), so `nan` and `inf` pass; the walker writes them into the summary as bare `NaN`/`inf` (`walker/src/io/output.rs:223-229`).
- `run_walker` returns `subprocess.run(...).returncode` (`cmd/walker.py:71`) and `sieve_run` returns `rc` (`cmd/sieve.py:222`), so a child killed by signal N exits 256−N (247 for SIGKILL).
- `run` intercepts only the exact token `--walker-bin` (`cmd/walker.py:77`), so `--walker-bin=PATH` reaches logan-walker as an unknown flag.

**Impact.** A user who pins a walker with `--walker-bin` or `DISASSEMBLER_WALKER_BIN` for reproducibility can run a different build with exit 0, likely where several conda builds coexist. `--min-cov nan` silently disables the coverage gate; `--min-abund inf` produces a summary that `json.loads` rejects. OOM kills are reported as 247 instead of 137.

**Evidence.** With no `logan-walker` on `PATH`, `disassembler walker --walker-bin /nonexistent/logan-walker --version` prints `logan-walker 0.1.1`, rc 0 (the in-tree build). `python -c "from disassembler.validators import positive_float as f; print(f('nan'), f('inf'))"` prints `nan inf`. A `fetch` with a `--walker-bin` script that runs `kill -9 $$` exits 247.

**Proposed fix.**

```diff
 def resolve_walker(explicit=None):
+    named = explicit or os.environ.get("DISASSEMBLER_WALKER_BIN")
+    if named:                      # a path the caller named wins or fails
+        return _executable(named)
     return (
-        _executable(explicit)
-        or _executable(os.environ.get("DISASSEMBLER_WALKER_BIN"))
-        or _from_path("logan-walker")
+        _from_path("logan-walker")
         or _executable(_repo_root() / "walker" / "target" / "release" / "logan-walker"))
```

Same change in `resolve_sieve`; `_resolved_or_report` should say "--walker-bin X is not an executable file". In `positive_float` use `if not math.isfinite(v) or v <= 0.0:`. Map negative return codes with `128 - rc if rc < 0 else rc` in `run_walker` and `sieve_run`. Accept and strip `--walker-bin=PATH` in `run`. In the walker, reject non-finite values and out-of-range `--min-cov`/`--min-vaf`/`--min-abund` in `parse_args`, and write non-finite floats as `null`.

**Tests and docs.** `tests/test_resolve.py` (explicit wins; env beats `PATH`; `PATH` beats in-tree; a missing or non-executable named path gives `None`; same for the sieve). `tests/test_validators.py` (edges of all validators, including `nan`/`inf`). `tests/test_fetch_walker.py` (argv forwarding, exit-code passthrough, self-killing stubs expecting 137/132, `--walker-bin=PATH`, missing value gives 2, no binary gives 127). A LexicMap fixture with minus-strand rows. Update the `resolve.py` docstring, `docs/walker.md:44` and `docs/sieve.md:42`; release notes must say a stale `DISASSEMBLER_WALKER_BIN` now exits 127. TC-22 proposes the coverage report that would keep such gaps visible.

### TC-04 · medium · The k=21 simulation misreports its own results and measures anchoring only

`bench/simulated_kmer_sensitivity/README.md:25` · science · confirmed (reproduced) · fix effort: small · sources: SC1-10, D2-10, X1-06

**Problem.** `docs/walker.md:344-345` cites this study as the basis of the default `-k 21`, but its text does not match its committed data and overstates what it measured. `README.md:25` (repeated at `:145`) says k=21, m=3 "Extends 100% gene recovery down to **80% identity**", while its own Table 2 (`:114`) and the CSV give 80.0 at 80% for bla_kpc and 44–80% across the six gene × fragmentation cells. Five min_shared=1 cells at 80% differ from the CSV (Table 1: k15 83 vs 85, k21 34 vs 39, k31 1 vs 5; Table 2: k21 92 vs 80, k31 4 vs 20). "Gene recovery" is `gene_found = 1 if n_ta > 0 else 0` (`harness.py:411`), that is one anchored unitig; the "5-hop BFS" (`:417-423`) is `for delta in range(-5, 6)` over list indices, not a graph walk, and no GraphAligner or `--min-cov 0.85` step is simulated, which the README does not say. Smaller errors: `:18` cites a retired `subgraph_fraction_guard`, and the figure `theoretical_vs_empirical_decay.png` plots per-k-mer `(1-d)^k` against per-unitig retention. `docs/walker.md:286-288` says SRR8209082's subgraph "is about 10%" (true at k=31; 42.4% at the default k=21) and that the walker warns above 1% (the guard is 3 unitigs per base, `walker/src/config.rs:469`).

**Impact.** A user choosing `-k` expects 100% detection at 80% identity; the data say 44–80%, and for anchoring only. k=19 is equal or better than k=21 in 5 of 6 cells at 80% with the same zero false positives, so the study does not single out 21.

**Evidence.** A pivot of `simulation_results.csv` (bla_kpc, moderate, uniform) against README Tables 1–2 shows the five differing cells and no others. Running the shipped walker on SRR8209082 with `queries/query_hiv1_pol.fasta`: k=31 gives 8,466 subgraph nodes (10.38%), k=21 gives 34,582 (42.41%).

**Proposed fix.** Generate README Tables 1–2 from the CSV with a script. Reword `:25`/`:145` to "100% anchoring down to 85%; 44–80% at 80%"; rename the metric "anchoring rate (≥1 unitig with ≥M shared k-mers)"; state that neither graph BFS on the target nor alignment is simulated; fix `:18` (retired guard) and `:63` ("Near Limit", although 27 of 81 rows exceed the guard); plot `1-(1-I^k)^(L-k+1)`. In `docs/walker.md:286-289` write "about 42% at the default `-k 21` (10% at `-k 31`)" and describe the 3-unitigs-per-base guard; at `:344` say the default rests on anchoring recall only. Later, add an end-to-end arm with GraphAligner at `--min-cov 0.85`.

**Tests and docs.** A CI script check that the README tables equal the CSV pivot. Rerunning needs TC-17; TC-16 covers the motif model behind the 60%-identity claim. Rerunning on the real background may change numbers and the `-k 21` default may need re-justifying.

### TC-05 · medium · Forwarding of 12 walker options from the composites is untested

`src/disassembler/cmd/pipe_reconstruct.py:357` · test-coverage · confirmed (reproduced) · fix effort: small · sources: T1-04

**Problem.** `reconstruct` turns 12 optional settings into logan-walker flags by hand (`walker_argv += ["-w", str(args.workers)]` and so on, `pipe_reconstruct.py:357-384`), and `from-kmindex` copies the same 12 values into the reconstruct namespace (`pipe_from_kmindex.py:358-378`). No test passes any of them: the stub walker records its argv, but tests assert only on `-l`. The docstring at `tests/test_composites.py:4-6` claims an option a composite fails to hand a stage "fails here", which holds only for misspelled `stage_args` keys. No CI job runs the composites against the real walker.

**Impact.** Nothing is broken today (all 12 flags exist in `walker/src/main.rs`). A renamed flag would make `reconstruct --prune-frac 0.02` fail at stage 1. A key dropped from `from-kmindex` makes `stage_args` fill in `None`, so `--min-vaf` or `--offline` is silently ignored: different output, or S3 traffic despite `--offline`, with CI green.

**Evidence.** In a copy of the repository, renaming all 12 flag strings in `pipe_reconstruct.py`, or deleting all 12 `key=args.key` lines in `pipe_from_kmindex.py`, still gives `Ran 168 tests ... OK`.

**Proposed fix.** Add `TestWalkerOptionForwarding` to `tests/test_composites.py`:

```python
FORWARDED = [(["-w", "3"], ["-w", "3"]),
             (["--min-abund", "2"], ["--min-abund", "2.0"]), ...]   # all 12
def _check(self, argv):  # mkdir the cache dir (--offline requires it); chdir, restore via addCleanup
    rc, log = _run(argv + [t for c, _ in FORWARDED for t in c] + ["--skip-tree"])
    got = self.walker_argv().split()
    for _, want in FORWARDED:
        i = got.index(want[0]); self.assertEqual(got[i:i + len(want)], want)
# test_reconstruct_forwards_every_walker_option, test_from_kmindex_forwards_every_walker_option
# test_the_walker_parses_every_forwarded_flag (no build needed):
#   assertTrue(re.search(rf'"{re.escape(flag)}"\s*(=>|\|)', main_rs_text), flag)
```

With these tests the renaming mutant fails 2 tests and dropping only `prune_frac` fails 1. Optionally add a CI step that runs `reconstruct --skip-tree` with all 12 options against the built walker (`--unitigs-cache DIR --offline`, stand-in GraphAligner) and checks the `settings` block of `*_walker_summary.json`.

**Tests and docs.** The three tests above. Narrow the docstring at `tests/test_composites.py:4-6` or keep it once the tests make it true. The source-grep test must change if the walker moves to clap.

### TC-06 · medium · Sieve gene-scale defaults (tau 0.005, sketch 256) are pinned by no test

`tests/test_sieve.sh:38` · test-coverage · confirmed (reproduced) · fix effort: small · sources: SWT-02

**Problem.** `docs/sieve.md` justifies the calibrated defaults threshold 0.005 and sketch 256, but no test exercises or pins them. `tests/test_sieve.sh` runs the binary at `THRESHOLD=0.15`, `WARMUP=10` and no `--sketch` (`:38-39`), and its header (`:13-23`) implies the pipeline collapses the fixture to one representative, which is true only at those archive-scale settings. The seam test asserts only `1 <= n < 33` (`tests/test_sieve_integration.py:145`); the sketch test passes an explicit 128; the equivalence tests hard-code `0.005`/`256` instead of reading the wrapper's defaults. The defaults are also repeated in five places (`cmd/sieve.py`, and the three composites).

**Impact.** If the default moves consistently, or an upstream sieve bump shifts the calibration, the fixture's promoted count moves from the documented 28 (`docs/sieve.md:117`) to 16 (sketch 64) or 2 (tau 0.03) with CI green. On real cohorts, lineages vanish before the MSA and tree. Drift in only one of the five places makes the entry points disagree, also untested.

**Evidence.** On `tests/fixture_reconstructed.fasta` with `--min-len 2000`, `l1_pipeline_serial` promotes 1 (0.15/10/no sketch), 28 (0.005/2000/256), 16 (0.005/2000/64) and 2 (0.03/2000/256). With the sketch default set to 64 in all five places, `python -m unittest discover -s tests -t tests` gives `Ran 168 ... OK` and `bash tests/test_sieve.sh` passes.

**Proposed fix.** Add `tests/test_sieve_defaults.py` (about 45 lines): (1) the `sieve`, `reconstruct`, `from-kmindex` and `from-lexicmap` subparsers' `get_default(...)` for threshold and sketch equal `inspect.signature(sieve.sieve_run)` defaults; (2) `sieve_run(FIXTURE, out, min_len=2000, sieve_bin=...)` with no threshold or sketch promotes 28. The sketch-64 mutant then fails with `16 != 28`, and a one-place tau drift in `from-kmindex` fails test 1. Moving the defaults to named constants with `default=None` in the composites would make test 1 structural.

**Tests and docs.** The new file; relabel `test_sieve.sh` as a smoke test of the binary's own defaults and fix its header. A legitimate submodule bump will fail test 2 by design; re-derive 28 and update `docs/sieve.md:117` together. TC-15 covers the warmup flag.

### TC-07 · medium · Even-spread sampling, max-merge, zero-drop, 429 handling and nested collections are unpinned

`tests/test_stratify.py:28` · test-coverage · confirmed (reproduced) · fix effort: small · sources: P3-09

**Problem.** The logic in `src/disassembler/cmd/harvest.py` and `src/disassembler/cmd/stratify.py` that decides which accessions reach the walker is not pinned. `test_normal_quotas_are_unchanged` checks only the length of `stratified_sample` on eight equal scores (`tests/test_stratify.py:28-31`), so the even spread at `stratify.py:58` could become top-N. The merge test feeds 0.4 then 0.7 in file order, which last-wins also satisfies (`tests/test_harvest.py:216-219`). `FakeResponse` has no `headers`, so the 429/Retry-After branch (`harvest.py:116-122`) never runs. Every fake collection is flat (recursion at `harvest.py:309-310` untested). The only order test uses one bin, so the cross-bin `picks.sort` (`stratify.py:137`) can be removed. The session-close test passes without the `finally` it claims to test, because the error it sees is raised at `harvest.py:386`, after the `finally`.

**Impact.** A refactor to top-N sampling, last-wins merging, keeping zero scores, dropping the Retry-After cap or skipping `list:list` collections changes the cohort silently with CI green. No current defect.

**Evidence.** Nine single-line mutants of `harvest.py`/`stratify.py`, each run through the full `unittest discover`: eight pass (`OK (skipped=3)`); only the stratify zero-score mutant fails.

**Proposed fix.** Tests only:
- `stratified_sample(list(range(10)), 5) == [0, 2, 4, 6, 8]` and `(range(9), 3) == [0, 3, 6]`.
- Two bins `[(0, 0.05, 2), (0.05, 0.10, 2)]` give `bins.txt` reading `hi2 hi1 lo2 lo1`.
- `harvest_local` with A=0.9 in `s1` and A=0.4 in `s2` keeps 0.9; a shard with `A: 0.0` leaves A out of `per_target` and `UNION.accessions.txt`.
- `FakeResponse(headers=...)`: `[429 Retry-After: 7, 200 "{}"]` sleeps `[7.0]`; `99999`, `soon` and a missing header sleep `[120, 30, 30]`.
- A `list:list` payload whose two inner datasets both reach `fetch_shard`.
- Session test: patch `hk.concurrent.futures.as_completed` to raise (restore in `finally`) and assert the session closed.

**Tests and docs.** Extend `tests/test_harvest.py` and `tests/test_stratify.py` (about 100 lines); all nine mutants then fail. Confirm with the maintainer that harvest should drop zero scores before pinning it.

### TC-08 · low · `main` is unprotected: CI is advisory and fork PRs were merged untested

`.github/workflows/ci.yml:5` · ci · confirmed (reproduced) · fix effort: trivial · sources: T1-07

**Problem.** CI triggers on push and PR to `main` (`ci.yml:5-8`), but nothing makes its result binding: the GitHub API reports `main` as `"protected": false` and there are no rulesets. The five fork PRs (#6, #8, #10, #11, #12) have zero check runs on their head commits and were merged; #12 broke clippy and the sieve integration test on `main`. Runs #20–#26 (2026-09-17 to 09-25) never started because of a billing failure, and six pushes landed with no test run, unnoticed. Both release commits did have green CI.

**Evidence.** `gh api repos/nekrut/disassembler/branches/main` shows `protected: false`; `gh api repos/nekrut/disassembler/rulesets` returns `[]`.

**Proposed fix.** Settings only: a ruleset on `main` requiring a PR and the checks `clippy + test (offline)`, `python tests (3.11)`, `python tests (3.12)`, `python tests (3.13)` and `msa against real mafft and fasttree`; approve pending fork runs before review; fix the spending limit so an outage shows as "expected check missing". Optionally add a CI badge to `README.md`.

**Tests and docs.** State the rule in `docs/packaging.md` or a CONTRIBUTING note. Keep required names in sync with job names.

### TC-09 · low · CI builds without `--locked` and checks no release identifiers; Bioconda builds with `--locked`

`.github/workflows/ci.yml:34` · ci · confirmed (reproduced) · fix effort: small · sources: SWK-03, PK1-04

**Problem.** Every cargo call in CI (`ci.yml:34`, `:36`, `:114`, `:137`), in `docs/packaging.md:60-61` and in `bench/check_golden.sh:33` runs without `--locked`, so cargo silently rewrites a stale `Cargo.lock`. The Bioconda recipe builds with `cargo install --no-track --locked`. A `Cargo.toml` change without the regenerated lock passes CI and fails only in the recipe build, after the tag is public. Nothing checks that tag, `__version__`, `Cargo.toml`, `Cargo.lock` and the sieve submodule pin agree; the pin regression (`39d90c0` instead of `70022da`, commit b1b6acf) shipped in v0.1.1 with CI green. Versions agree today.

**Evidence.** In a copy of `walker/`, bump `Cargo.toml` to 0.1.2: `cargo metadata --locked` exits 101 ("cannot update the lock file … because --locked was passed"); without `--locked` it exits 0 and rewrites `Cargo.lock`.

**Proposed fix.**

```yaml
on:
  push:
    branches: [main]
    tags: ['v*', 'logan-walker-v*']
...
      - run: cargo clippy --locked --all-targets --release -- -D warnings
      - run: cargo test --locked --release
  release-consistency:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v5
      - run: bash .github/release_consistency.sh
```

Also add `--locked` at `ci.yml:114` and `:137`, `bench/check_golden.sh:33` and packaging step 2. `release_consistency.sh` checks `Cargo.toml` = `Cargo.lock`, that the submodule pin matches the commit named in `docs/packaging.md:42` (a one-line `extern/SIEVE_COMMIT` file would be sturdier than grepping docs), and on a tag ref the version check from TC-21 (the two tag series are independent). Re-pin the submodule to `70022da` first, or the job is red from the start.

**Tests and docs.** A PR editing `Cargo.toml` without the lock must go red. Add "`cargo check --locked` must pass before tagging" to packaging step 1.

### TC-10 · low · Golden gate green-lights a stale prebuilt walker without rebuilding or naming it

`bench/check_golden.sh:32` · test-defect · confirmed (reproduced) · fix effort: trivial · sources: X1-03

**Problem.** `if [ ! -x "$WALKER" ]; then … cargo build --release …; fi` (`bench/check_golden.sh:32-34`) builds only when no binary exists, and the script prints the GraphAligner version but not the walker's. A developer who edits `walker/src` and runs the gate locally, as its header advises, gates the previous build; even a tree that does not compile passes with `golden: concordant`. The CI job builds first, so only local runs are affected.

**Proposed fix.**

```diff
-if [ ! -x "$WALKER" ]; then
-  echo "building walker..."; ( cd "$ROOT/walker" && cargo build --release ) || exit 1
+if [ -z "${LOGAN_WALKER_BIN:-}" ]; then
+  echo "building walker..."; ( cd "$ROOT/walker" && cargo build --release --locked ) || exit 1
 fi
+[ -x "$WALKER" ] || { echo "FAIL: no walker binary at $WALKER"; exit 1; }
+echo "walker: $("$WALKER" --version 2>&1 | head -n 1) at $WALKER"
```

**Tests and docs.** Document `LOGAN_WALKER_BIN` (currently undocumented) in the header (`:7-10`) and `bench/golden/README.md`.

### TC-11 · low · Real-tool MSA test passes with `--keeplength` removed

`tests/test_msa_real_tools.py:75` · test-defect · confirmed (reproduced) · fix effort: trivial · sources: T1-03, P2-10

**Problem.** `test_every_row_is_reference_length` calls itself "the only check that would notice MAFFT quietly ignoring" `--keeplength` (`:76-77`), but no fixture row has an insertion relative to the 2739-nt reference, so rows are reference-width with or without the flag (`src/disassembler/cmd/msa.py:282`). The other MSA tests stub MAFFT. A regression is not silent: the length guard at `msa.py:298-304` stops the run on the first cohort with an insertion. Deleting the `"--keeplength",` line still gives `Ran 16 tests ... OK` for `tests.test_msa_real_tools tests.test_msa`.

**Proposed fix.** In `setUpClass`, append a row with an in-frame insertion and run `msa` on that copy:

```python
src = read_fasta(INPUT); cls.ins_src = next(iter(src)); seq = src[cls.ins_src]
cls.input = cls.work / "input_with_insertion.fasta"
cls.input.write_text(INPUT.read_text().rstrip("\n")
                     + f"\n>INS_test\n{seq[:1200]}GGTACCGGT{seq[1200:]}\n")
# then pass --input cls.input, and assert rows["INS_test"] == rows[cls.ins_src]
```

**Tests and docs.** Add a stub test in `tests/test_msa.py`: an over-long MAFFT row gives exit 1, the "length != ref length" message, and an untouched earlier alignment. Confirm the exact-equality assertion on MAFFT 7.526, or relax it to width plus "no `GGTACCGGT`".

### TC-12 · low · Reference-quality test skips `queries/` and FASTA extensions other than `.fa`/`.fasta`

`tests/test_panels.py:40` · test-coverage · confirmed (reproduced) · fix effort: small · sources: X1-08, DA-06

**Problem.** `panels()` returns only `PANEL_DIR.glob("*.fa")` and `glob("*.fasta")` with `PANEL_DIR = REPO_ROOT / "panels"` (`:23`, `:40`). Nothing under `queries/` is checked, including `master_planetary_targets.fasta`, which ships seven records with N/IUPAC codes; a future `panels/*.fna`, `*.fas` or `*.fa.gz` would be skipped silently. A `panels/bad.fna` containing `ATGNNNTAA` produces no failure.

**Proposed fix.** Enumerate every walker-readable FASTA and fail on unknown files:

```python
FASTA_EXT = (".fa", ".fasta", ".fna", ".fas")
def panels():
    ps = sorted(p for p in PANEL_DIR.iterdir() if p.name != "README.md")
    bad = [p.name for p in ps if not p.name.removesuffix(".gz").endswith(FASTA_EXT)]
    assert not bad, f"unrecognised file(s) in panels/: {bad}"
    return ps
```

Open `.gz` with `gzip.open(..., "rt")`. Add `tests/test_query_fastas.py` over `queries/**`: ACGT only, length % 3 == 0, no internal stop, unique ids (not ATG/terminal stop: 60 NS5 mat_peptide records lack them).

**Tests and docs.** The new test fails on exactly the seven records until they are fixed; land both together. Update `panels/README.md:5-6`.

### TC-13 · low · Sieve tests skip green when the submodule is absent or the build fails

`tests/test_sieve.sh:51` · ci · confirmed (reproduced) · fix effort: trivial · sources: T1-08

**Problem.** `skip() { echo "SKIP: $*"; exit 0; }` (`tests/test_sieve.sh:46`) runs when the submodule Makefile is missing (`:51`) or `cc` is absent (`:52`). The Python modules raise `SkipTest` on "no-submodule" and on "no-compiler", which `ensure_sieve` returns after any failed `make`, not only a missing compiler (`tests/test_sieve_integration.py:40-41`, `:54-55`). CI accepts these skips, unlike the real-tools job. Today CI checks out with `submodules: recursive`, so the gap is latent; a non-recursive clone or release tarball reports green with no sieve code run.

**Proposed fix.** Opt-in strictness, set in CI:

```bash
# ci.yml python-tests: env: { DISASSEMBLER_REQUIRE_SIEVE: "1" }
skip() { [[ -n "${DISASSEMBLER_REQUIRE_SIEVE:-}" ]] && fail "sieve required, but: $*"; echo "SKIP: $*"; exit 0; }
```

In Python, route skips through one helper that raises `RuntimeError` when the variable is set. Report a failed `make` with `cc` present as a failure with make's stderr.

**Tests and docs.** Update the comments at `ci.yml:58-61` and `test_sieve.sh:48-50`. Share the helper with TC-14.

### TC-14 · low · Sieve test bootstrap trusts any existing binary and calls `make` unguarded

`tests/test_sieve_integration.py:35` · test-defect · confirmed (reproduced) · fix effort: small · sources: T1-09, BL2-10

**Problem.** Both Python sieve modules use `if SIEVE_BIN.exists(): return "ok"` and otherwise call `subprocess.run(["make", "serial"], ...)` unwrapped (`tests/test_sieve_integration.py:35-40`, `tests/test_dedup_sieve_equivalence.py:30-34`). A binary built with `-march=native` on another CPU is used as is and never rebuilt (make sees it as current): 6 failures and 4 errors that report `SIGILL`, not the cause. On a host without `make`, `setUpClass` raises `FileNotFoundError` (3 errors) although the docstrings promise a clean skip. The in-place build can also race between two runs on one tree.

**Proposed fix.** One `ensure_sieve()` in a shared test module: no Makefile gives "no-submodule"; an existing binary that survives a probe run on a 2-record FASTA gives "ok"; otherwise `make -B serial BIN_SERIAL=.l1_pipeline_serial.<pid>` inside `try/except OSError`, then `os.replace` into place; a failed build or missing `make` skips; a binary still dying on a signal gives "unusable". In CI these outcomes fail rather than skip, through the TC-13 helper. In `test_sieve.sh`, probe after the `-x` check and rebuild once if the exit status is above 128.

**Tests and docs.** Fix the "Skips cleanly" docstrings in both modules and backlog item 5.

### TC-15 · low · No test pins `sieve_run`'s binary flags; the equivalence test hard-codes its own copy

`tests/test_dedup_sieve_equivalence.py:62` · test-defect · confirmed (reproduced) · fix effort: small · sources: PS-10

**Problem.** `_run_sieve` builds the command itself ("the same flags the composite uses", `:63`) with a literal `"--warmup", "2000"` (`:73`), while production builds it separately (`src/disassembler/cmd/sieve.py:208-217`, warmup at `:214`); `tests/test_sieve_integration.py:243-246` holds a third copy. No test asserts the warmup, so the documented below-warmup premise (`sieve.py:65-69`, `docs/sieve.md:54`) is tied to nothing. Changing `sieve.py:214` to `"--warmup", "10"` leaves all 18 sieve tests green.

**Proposed fix.** In `sieve.py` add `WARMUP = 2000` and `sieve_command(binary, input, emit_prefix, *, min_len, max_gap_frac, threshold, sketch) -> list[str]`, used by `sieve_run`. The equivalence test calls `sieve.sieve_command(...)` and asserts `len(records) < sieve.WARMUP`; one `sieve_run` test asserts `"--threshold 0.005 --warmup 2000 --sketch 256"` in the logged command. With `WARMUP = 10` three tests then fail.

**Tests and docs.** The two assertions. Related: TC-06 (threshold and sketch defaults).

### TC-16 · low · Conserved-motif model protects 30-nt windows, shorter than a 31-mer

`bench/simulated_kmer_sensitivity/harness.py:143` · science · confirmed (reproduced) · fix effort: small · sources: X1-07

**Problem.** `mutate_clustered_motifs` hard-codes three protected windows `(int(0.20 * L), int(0.20 * L) + 30)` (`:143-145`), 10% of an 882-nt gene, not the "~20%" its docstring states (`:132`). A 30-nt window holds ten 21-mers and no complete 31-mer, so the headline "k = 21 maintains 96–100% … down to 60% identity (whereas k = 31 drops to 48%)" (`README.md:26`, `:146`) follows mostly from that width. With the same code, bla_kpc, m=3, 200 replicates, k=31 at 60% identity gives 48.5% at width 30, 61.0% at 31, 86.0% at 45 and 97.0% at 60. The k=21-over-k=31 conclusion still holds on the uniform and codon models.

**Proposed fix.** Add a `motif_width` parameter (`--motif-widths 12,30,45,60`) as a grid axis, or place blocks at each gene's real motif coordinates. Fix the docstring ("~10%, three 30-nt blocks") and reword README §1.4, §5 and §6.1 to state the block width and show the sweep.

**Tests and docs.** A unit test that the function protects exactly 3×width positions. Rerun with TC-17 in place. Related: TC-04.

### TC-17 · low · Simulation harness cannot regenerate its committed results

`bench/simulated_kmer_sensitivity/harness.py:488` · reproducibility · confirmed (reproduced) · fix effort: small · sources: X1-05, BL2-11, X1-02

**Problem.** `REPO_ROOT = HERE.parent.parent.parent` (`:28`) is above the checkout and `CDS_DIR`/`ECOLI_PATH` (`:29-30`) point into a `work/` tree, so the script dies as committed. With the gene path fixed, a missing genome is silently replaced by a random sequence (`:488-491`), `--replicates` defaults to 30 while the committed CSV used 25 (`:471`), and `--out-dir` defaults to the script directory, overwriting the tracked CSV/JSON (`:473`). The background lengths come from `np.random.lognormal` (`:296`) before the only `np.random.seed` (`:335`), so even with the real genome three builds give 19/295, 23/243 and 22/265 k=15 false-positive anchors/nodes.

**Evidence.** `python3 bench/simulated_kmer_sensitivity/harness.py --replicates 1 --out-dir /tmp/x` raises `FileNotFoundError: …/work/cds/individual/bla_vim.fasta`.

**Proposed fix.**

```diff
-REPO_ROOT = HERE.parent.parent.parent
-CDS_DIR = REPO_ROOT / "work" / "cds" / "individual"
+REPO_ROOT = HERE.parents[1]
+CDS_DIR = REPO_ROOT / "queries" / "markers"
-    lengths = np.random.lognormal(mean=5.0, sigma=0.4, size=n_unitigs).astype(int)
+    lengths = np.random.default_rng(seed).lognormal(mean=5.0, sigma=0.4, size=n_unitigs).astype(int)
```

Add `seed=42` to `build_background_graph` and remove the dead global seed; default `--replicates` to 25; require `--background-fasta FASTA` or `--synthetic-background`; refuse to overwrite without `--force`; write a `simulation_run.json` sidecar (replicates, background path or "synthetic", its sha256, seed).

**Tests and docs.** CI smoke run with `--replicates 1 --bg-unitigs 300 --synthetic-background`. Regenerate README figures (`:18`, `:63`, `:77`); document fetching NC_000913.3.

### TC-18 · low · Golden gate labels a walk that found nothing as "NOT a walk regression"

`bench/check_golden.sh:85` · usability · confirmed (reproduced) · fix effort: trivial · sources: SWT-04

**Problem.** Every non-zero walker exit prints `ERROR $acc (walker exited non-zero: a run failure, NOT a walk regression)` (`:85`, same claim at `:83-84`). But the walker exits 3 only when every accession errored and 1 when the walk ran and ended `no_anchors`/`no_path` (`walker/src/main.rs:401-414`); with one accession per gate, exit 1 is a walk regression. The gate still fails, but sends the developer to S3 or GraphAligner, and `:110` then counts it as a regression.

**Proposed fix.**

```bash
  run_primary "$acc" "$query"; rc=$?
  if [ "$rc" -ne 0 ]; then
    status="$(sed -n 's/.*"status": *"\([a-z_]*\)".*/\1/p' "$WORK/$acc"/*_walker_summary.json 2>/dev/null | head -n 1)"
    case "$rc:$status" in
      1:no_anchors|1:no_path) echo "FAIL $acc  (walk reconstructed nothing: status=$status)" ;;
      *) echo "ERROR $acc  (walker exit $rc, status=${status:-unknown}: a run failure, NOT a walk regression)" ;;
    esac
    show_log "$acc"; fails=$((fails+1)); continue
  fi
```

`no_logan_data` (404) stays ERROR. Fix the `:83-84` comment; make `:110` print `golden: $fails gate(s) failed`.

**Tests and docs.** A shell test with an offline cache sharing no k-mers with the queries and a stub GraphAligner.

### TC-19 · low · Golden gate ignores a failed `mktemp` and writes its work files into `/`

`bench/check_golden.sh:29` · reliability · confirmed (reproduced) · fix effort: trivial · sources: X1-12

**Problem.** The script runs under `set -uo pipefail` without `-e` and never checks `WORK="$(mktemp -d)"` (`:29`). With `TMPDIR` pointing at a missing directory, `WORK` is empty, so logs, outputs and comparisons go to `/ERR2028820.log`, `/ERR2028820/` and so on (`:57-61`, `:90`), and the EXIT trap removes nothing. As root the gate reports `golden: concordant` and leaves the files in `/`; as non-root every accession is reported as a walker failure.

**Proposed fix.**

```diff
-WORK="$(mktemp -d)"
+WORK="$(mktemp -d)" && [ -d "$WORK" ] || { echo "FAIL: cannot create a temp directory (TMPDIR=${TMPDIR:-unset})"; exit 1; }
```

Adding `-e` is wrong: the script relies on non-zero statuses from `run_primary` and `grep`.

### TC-20 · low · CI lints with floating `stable` Rust under `-D warnings`; cache key ignores the compiler

`.github/workflows/ci.yml:23` · ci · confirmed (by reading) · fix effort: trivial · sources: X1-11

**Problem.** Both Rust jobs install whatever `stable` is that day (`ci.yml:23`, `:106`) and the gate is `clippy -- -D warnings` (`:34`); there is no `rust-toolchain.toml` and no `rust-version`. Run 39180037252 (2026-10-19) picked up rustc 1.99.1 and failed on new lints (`manual_inspect`, `cloned_ref_to_slice_refs`) with no code change. The cache key `${{ runner.os }}-cargo-${{ hashFiles('walker/Cargo.lock') }}` (`:31`) has no compiler component, so after a toolchain bump every run rebuilds and the cache is never re-saved until `Cargo.lock` changes.

**Proposed fix.** Add `walker/rust-toolchain.toml`:

```toml
[toolchain]
channel = "1.99.1"
components = ["clippy", "rustfmt"]
profile = "minimal"
```

Run `rustup toolchain install` in `walker/` in both jobs; key the cache on `rustc -V` plus the lock hash with a `restore-keys` prefix, or use `Swatinem/rust-cache@v2`; set `rust-version` in `walker/Cargo.toml`. Optionally a weekly non-blocking job on floating `stable`.

**Tests and docs.** Note the pinned toolchain in `walker/README.md` and the bump step in `docs/packaging.md`.

### TC-21 · low · CI hygiene: no dependency audit, no job timeouts or permissions, no tag-time version check

`.github/workflows/ci.yml:3` · ci · confirmed (by reading) · fix effort: small · sources: T1-14

**Problem.** There is no `.github/dependabot.yml` and no `cargo audit`/`cargo deny` step for the 61 locked crates (rustls, ring, webpki-roots, zstd-sys); no job sets `timeout-minutes` (a hang holds a runner for the 6-hour default, and the golden job downloads Logan data); there is no top-level `permissions:`; `on:` (`:3-8`) has no tag trigger, so nothing checks a tag against its version; `actions/setup-python@v5` (`:53`) already draws a Node 20 deprecation annotation.

**Proposed fix.** `dependabot.yml` for `cargo` (`/walker`), `pip` and `github-actions`, weekly; a weekly `cargo audit` job; `timeout-minutes: 20` (90 for golden); `permissions: { contents: read }`; `actions/setup-python@v6`. On tag pushes check only the matching series, since the two tag series are independent (`docs/packaging.md:63`):

```bash
case "$GITHUB_REF_NAME" in
  logan-walker-v*) want=${GITHUB_REF_NAME#logan-walker-v}; have=$(grep -m1 '^version' walker/Cargo.toml | cut -d'"' -f2);;
  v*) want=${GITHUB_REF_NAME#v}; have=$(python -c 'import sys; sys.path.insert(0,"src"); import disassembler; print(disassembler.__version__)');;
esac; test "$want" = "$have"
```

**Tests and docs.** Mention the tag check and audit job in the `docs/packaging.md` checklist. Related: TC-09, TC-20.

### TC-22 · low · CI runs only on Linux and has no coverage signal, while macOS packages are planned

`.github/workflows/ci.yml:16` · ci · confirmed (by reading) · fix effort: small · sources: T1-13

**Problem.** All four jobs use `runs-on: ubuntu-latest` (`ci.yml:16`, `:40`, `:70`, `:102`) and none measures coverage. Bioconda publishes linux-64 only today, but backlog item 0.2 drops the osx skips, so macOS packages would ship code no upstream job has compiled there. The walker has macOS-relevant paths (`/dev/shm` fallback at `walker/src/scratch.rs:106-114`, `#[cfg(unix)]` libc calls) and the sieve Makefile has an arm64 branch. No Linux-only API was found, so a break is plausible, not observed.

**Proposed fix.** Matrix `build-test` over `[ubuntu-latest, macos-13]` (x86_64, matching the only osx GraphAligner build), tests only, skipping `empty_subgraph_yields_nothing` until the known `/dev/null` test bug is fixed; add one `macos-latest` entry to `python-tests`. In one Python entry run `coverage run -m unittest discover -s tests -t tests && coverage report`, report-only until a baseline exists, then `--fail-under` at or below it.

**Tests and docs.** State in `docs/packaging.md` which platforms CI covers. Related: TC-03.

### TC-23 · low · No Python lint in CI; the existing debt is 13 auto-fixable lines

`.github/workflows/ci.yml:62` · ci · confirmed (reproduced) · fix effort: trivial · sources: BL2-12

**Problem.** The Python job runs only the unit tests and `test_sieve.sh` (`ci.yml:62-66`); `pyproject.toml` has no linter config, so an undefined name in an untested error branch surfaces only as a user's `NameError`. `ruff check --isolated --no-cache --select F,E9 src tests bench` reports 13 issues (7 unused imports, 6 f-strings without placeholders), all auto-fixable, none a bug. In the walker, `#![allow(dead_code, unused_variables)]` (`walker/src/main.rs:24`) hides 9 dead items from clippy, including `Anchors::fraction_of_graph`, which a comment at `anchor.rs:311` says is "still reported".

**Proposed fix.** Run `ruff check --fix --select F,E9 src tests bench`; add `[tool.ruff.lint] select = ["F", "E9", "B006", "B904", "PLE"]` to `pyproject.toml`; add a CI step `pipx run ruff==0.16.9 check src tests bench` (pinned). Replace the crate-wide allow with item-level `#[allow(dead_code)]` or delete the items.

### TC-24 · low · Rust sources are not rustfmt-clean and CI has no format check

`.github/workflows/ci.yml:33` · ci · confirmed (reproduced) · fix effort: small · sources: K3

**Problem.** `cd walker && cargo fmt --check | grep -c '^Diff in'` prints 332, across all 17 source files and `tests/cli.rs`; there is no `rustfmt.toml`, and CI runs clippy but no format check. Any contributor whose editor formats on save produces hundreds of unrelated changed lines, which hides real changes in review and pollutes `git blame`.

**Proposed fix.** Choose the style once (rustfmt defaults, or a `walker/rustfmt.toml` such as `max_width = 120`, `use_small_heuristics = "Max"`), run `cargo fmt` in one formatting-only commit right after a release, list it in `.git-blame-ignore-revs`, and add `cargo fmt --check` before clippy in `build-test`.

### TC-25 · low · `bench/golden/README.md` still documents the deleted archive, verifier and harness

`bench/golden/README.md:48` · docs · confirmed (by reading) · fix effort: small · sources: X1-10, D2-12

**Problem.** Commit f5bcb4c deleted `bench/archive/` (`capture_golden.py`, `verify_golden.py`, `run_benchmark.sh`), but `bench/golden/README.md` still points there (`:47-48`, `:54`), and `:58-96` (asserted fields, hash-seed notes, "Regression coverage (verified)", `run_benchmark.sh`) describe tooling that no longer exists. They read as current guarantees, while the gate compares only the primary. The two `*.golden.json` files (about 3 MB) are read by nothing.

**Proposed fix.** At `:47-48` say the tools were removed in f5bcb4c (`git show f5bcb4c^:bench/archive/verify_golden.py`); keep one paragraph stating the gate byte-compares the primary only and the JSONs are a record of the 2026-09-10 capture; collapse `:58-96` into a short history note, keeping the SRR8209082 caveat (`:77-81`). Optionally delete the JSONs. A CI grep for `\.\./archive/|run_benchmark` would catch such pointers.

### TC-26 · low · Test scaffolding is duplicated across files; one fixture reader assumes single-line FASTA

`tests/test_composites.py:59` · reuse · confirmed (reproduced) · fix effort: small · sources: PS-12

**Problem.** The tests carry five FASTA readers, three writers, two byte-identical MAFFT/FastTree stub builders (`tests/test_msa.py:38-46` and `tests/test_composites.py:59-67`), two fake walkers, two sieve bootstraps (TC-14) and three manual `PATH` swaps. `TestSieveRun.setUpClass` appends every non-header line as a record (`tests/test_sieve_integration.py:199-204`), so rewrapping `tests/fixture_reconstructed.fasta` at 60 columns makes five of its tests fail with misleading messages ("completeness gate drops 70 of 70").

**Proposed fix.** Add `tests/_helpers.py` (importable under `discover -t tests`) with a multi-line `read_fasta`, `write_fasta`, `count_records`, the shared `ensure_sieve` from TC-14, and a `ToolStubs` mixin that writes stubs and patches `PATH` with `mock.patch.dict`. First switch `TestSieveRun.setUpClass` to the shared reader; then replace the copies file by file. Test counts must stay the same.

### TC-27 · info · Two GraphAligner 1.0.17b builds are in use and no record tells them apart

`bench/golden/GRAPHALIGNER_VERSION:1` · provenance · plausible · fix effort: small · sources: PK2-13, PK1-06

**Problem.** Bioconda has three linux-64 builds of `graphaligner 1.0.17b`. Environments with Python ≥ 3.11 (`environment.yml`, the disassembler image) resolve `hdcf5f25_1`; environments without Python, such as `biocontainers/logan-walker:0.1.1--hfa8f182_0` and probably the CI golden job (`ci.yml:112`), take `h21ec9f0_2`, which pulls protobuf 3.14.0 and Python 3.9. All builds print `GraphAligner bioconda 1.0.17b-`, which is all that `GRAPHALIGNER_VERSION`, `check_golden.sh` and the walker summary record. No output difference has been shown.

**Proposed fix.** Run `bench/check_golden.sh` once in each environment and record the result in section 10 of `docs/issue_graphaligner_end_clipping.md`. Record the build string (`… (conda build hdcf5f25_1)`) in `GRAPHALIGNER_VERSION`; have `check_golden.sh` print `conda list '^graphaligner$'` when `CONDA_PREFIX` is set; optionally add the binary's sha256 to the summary JSON. Do not pin `*_1` in the recipe: osx-64 has only `_0`/`_2`.
