# Walker test suite: unit tests in walker/src and the CLI tests in walker/tests/cli.rs

The walker has a large unit suite (about 220 tests) that covers parsing, anchoring, extraction, stitching and the writers piece by piece, and the parsers and graph code are well pinned on small fixtures. The suite is weak in two places. First, it is not hermetic: `impl Drop for Subgraph` unlinks any path a test fixture hands it, so as root `cargo test` deletes `/dev/null`. Other tests also fail or leak depending on `LOGAN_KEEP_SCRATCH`, proxy variables or the shared `/dev/shm`. Second, no automated test goes through a successful reconstruction. The CLI stand-in aligners never write a GAF, the golden gate runs only on manual dispatch and with non-default flags, and several tests check a copy of the logic (a sort, a predicate, a dead parser) rather than the code `run` executes. Two default-path defects already shipped through these gaps. Most fixes share one theme: tests must call the production function and own their files, and one GAF-writing stand-in aligner in `tests/cli.rs` would close most of the end-to-end gap.

| ID | Severity | Where | Title | Fix effort |
|---|---|---|---|---|
| TW-01 | high | `walker/src/stages/phase.rs:444` | `Subgraph`'s Drop unlinks paths that test fixtures never created; as root `cargo test` deletes /dev/null | small |
| TW-02 | medium | `walker/src/stages/variants.rs:735` | No test covers consensus correction, alt recovery or band widening with realistic alignment windows | medium |
| TW-03 | medium | `walker/tests/cli.rs:37` | No automated test reaches a successful reconstruction or `--fetch-only` | small |
| TW-04 | medium | `walker/src/main.rs:695` | No test checks that each walker flag sets the field it names, or that unknown flags are refused | small |
| TW-05 | medium | `walker/src/pipeline.rs:961` | The input-order test sorts its own vector; removing `run`'s sort passes every test | small |
| TW-06 | low | `walker/src/stages/phase.rs:450` | Phase and extract tests do not pin strands, ordering, the fork filter, the prune-floor tie or the edge guards | small |
| TW-07 | low | `walker/src/stages/align.rs:326` | `parse_gaf` is dead code; seven parser tests pin it instead of the live `parse_gaf_multi` | small |
| TW-08 | low | `walker/tests/cli.rs:67` | `LOGAN_KEEP_SCRATCH` in the environment fails four walker tests | trivial |
| TW-09 | low | `walker/src/stages/stream.rs:413` | `size_cap_tests` fail when a proxy variable is set without `NO_PROXY` covering 127.0.0.1 | trivial |
| TW-10 | low | `walker/src/stages/extract.rs:742` | Walker tests leak empty temp directories and use the real `/dev/shm` scratch base | small |
| TW-11 | low | `walker/src/io/input.rs:601` | Two tests in the data-backed CI step pass when their inputs are broken or absent | trivial |
| TW-12 | low | `walker/src/stages/extract.rs:566` | The recall-study subgraph test passes on a corrupt fixture and reports it as "no data" | trivial |
| TW-13 | low | `walker/src/stages/extract.rs:884` | `emits_gfa_for_cross_check` deletes the GFA it reports having written | trivial |

### TW-01 · high · `Subgraph`'s Drop unlinks paths that test fixtures never created; as root `cargo test` deletes /dev/null

`walker/src/stages/phase.rs:444` · test-defect · confirmed (reproduced) · fix effort: small · sources: K1, T1-12

**Problem.** `impl Drop for Subgraph` (`walker/src/types.rs:205-212`) runs `std::fs::remove_file(&self.gfa_path)` whenever `LOGAN_KEEP_SCRATCH` is unset, whoever built the value. Unit tests build `Subgraph` by hand with paths they do not own. The test `empty_subgraph_yields_nothing` uses `gfa_path: PathBuf::from("/dev/null")` (`walker/src/stages/phase.rs:444`). `stitch.rs:253` uses `/tmp/logan-stitch-test-nonexistent.gfa`, under the comment "Points at nothing; Drop's remove_file is allowed to fail", and `variants.rs:721` uses `/tmp/logan-variants-test-nonexistent.gfa`. `align.rs:641` already works around the same hazard with `std::mem::forget(sub)`.

**Impact.** As a normal user the unlink of `/dev/null` fails silently, so GitHub Actions never shows it. As root (Docker, CI or Bioconda build containers, `sudo cargo test`) it succeeds. All 10 `tests/cli.rs` tests then fail with `NotFound`, because they open `/dev/null` for child stdio. `cargo test --test cli` alone passes, so the failure looks flaky. The next shell redirect `> /dev/null` creates a regular file that stores everything later written to it, until someone runs `mknod -m 666 /dev/null c 1 3`. Any file that happens to exist at the two fixed `/tmp` names is deleted on every run.

**Evidence.** As root in a container, `cd walker && cargo test --release` gave `222 passed` for the unit tests, then `0 passed; 10 failed` for `tests/cli.rs`. Afterwards `ls -la /dev/null` showed a regular file, and bisecting by test filter pointed to `stages::phase`. In a second run, `echo precious > /tmp/logan-stitch-test-nonexistent.gfa` (and the same for the variants name) followed by the stitch and variants tests left neither file in place.

**Proposed fix.** Make ownership explicit, so Drop removes only a file that `extract::write_gfa` created (it is the only production constructor, `extract.rs:302`):

```rust
// types.rs
pub struct Subgraph { /* ... */ pub n_edges: usize,
    /// True only when extract::write_gfa created gfa_path; Drop unlinks only then.
    pub owns_gfa: bool }
impl Drop for Subgraph { fn drop(&mut self) {
    if self.owns_gfa && std::env::var_os("LOGAN_KEEP_SCRATCH").is_none() {
        let _ = std::fs::remove_file(&self.gfa_path); } } }
// extract.rs write_gfa, align.rs:675 and :793 (tests that write a real temp GFA): owns_gfa: true
// phase.rs:442, stitch.rs:248, variants.rs:717, align.rs:627: owns_gfa: false
// align.rs:641: delete the std::mem::forget(sub) workaround
```

This was applied and tested: the full unit suite passes, and the planted `/tmp` files survive. A stronger variant puts the path in its own type, `pub gfa: Option<ScratchFile>`, where only `ScratchFile` has a Drop. Alternatively, keep the field private behind a constructor so a future production constructor cannot forget `owns_gfa: true`.

**Tests and docs.** Add `dropping_a_hand_built_subgraph_leaves_its_path_alone`: write a temp file, drop a `Subgraph` with `owns_gfa: false` and assert the file exists, then drop one with `true` and assert it is gone. Add a final CI step `test -c /dev/null`, or run the Rust job once as root in a container. Update the doc comment on the Drop impl (`types.rs:198-204`) and remove the stale comment at `stitch.rs:252`.

### TW-02 · medium · No test covers consensus correction, alt recovery or band widening with realistic alignment windows

`walker/src/stages/variants.rs:735` · test-coverage · confirmed (reproduced) · fix effort: medium · sources: R5-07

**Problem.** Every `variants.rs` test uses the fixture `walked()`, which has `path_start: 0, path_end: 7` with `path_len: 7` (`variants.rs:734-736`). The aligned window therefore always equals the whole path. In real GAFs, `path_start > 0` and `path_end < path_len`. The CLI stand-in aligners never write a GAF (see TW-03). The golden gate runs with `--prune-frac 0.0 --no-prefer-majority --no-emit-alts` (`bench/check_golden.sh:59`) and only on manual dispatch (`.github/workflows/ci.yml:101`). In `io/output.rs`, the random band test makes at most 3 indels of at most 9 nt (`output.rs:914-916`). The test `a_wildly_divergent_pair_falls_back_to_the_exact_dp` is accepted by the first band (66), so it never exercises the fallback.

**Impact.** Regressions in the default-on paths (majority correction, alt-window shifting, band-edge detection) pass `cargo test` and CI. Two such defects shipped in 0.1.0 and 0.1.1 this way: majority correction is skipped when the alignment window lies inside the path, and the banded aligner can return a non-optimal alignment.

**Evidence.** A new test with a 40/60 diamond and window [1,5) failed on the current code (`majority not swapped in`), while the existing suite passed. Disabling the band-edge refusal in `align_banded` (`if false && ...`) left all 16 `io::output` tests green. Instrumenting the existing output tests counted `widened=0 fallbacks=0`. For an ins40/del40 pair, the banded and exact alignments disagree at 948 of 2400 positions.

**Proposed fix.**
1. Parameterise `walked(start, end)`. Add window-inside-path cases ([1,5), [0,5)) for a substitution and an in-frame indel, for both `prefer_majority` and `find_alts`.
2. In `tests/cli.rs`, add a fixture with the 40/60 diamond graph and a stand-in aligner that writes a GAF line with `path_start=40`. Assert that the primary carries the 60x base, the alt has `vaf<=0.5`, the VCF row has the expected POS/REF/ALT, and `summary.json` records `prefer_majority`.
3. In `output.rs`, add an ins40/del40 test compared against `align_full`, and a direct test that `align_banded` returns `None` when the optimum touches the band edge. Rewrite the divergent-pair test to assert `align_banded(.., first_band).is_none()`.
4. Optionally add a golden capture under default flags.

Steps 1-2 fail until the majority-correction defect is fixed, and step 3 fails until the banded-alignment defect is fixed, so land each test with its fix.

**Tests and docs.** Document in `bench/README` or `check_golden.sh` which flags the golden gate covers. Correct the "never a wrong answer" doc comment at `output.rs:714-719` once the band fix lands.

### TW-03 · medium · No automated test reaches a successful reconstruction or `--fetch-only`

`walker/tests/cli.rs:37` · test-coverage · confirmed (reproduced) · fix effort: small · sources: T1-02, R6-08

**Problem.** The CLI stand-in aligners never produce an alignment: `for (name, body) in [("ga_fast", "exit 0"), ("ga_slow", "sleep 5")]` (`walker/tests/cli.rs:37`). Even the `matching` fixtures stop at `no_path`, and every exit-0 assertion is paired with `--allow-empty`. No test calls `pipeline::run`, `run_one` (only caller: `pipeline.rs:481`) or `run_fetch_only` (`pipeline.rs:186`), and none passes `--fetch-only`. The only end-to-end success check is the manual golden gate, which runs with non-default flags (see TW-02).

**Impact.** A regression in the post-alignment half of `run_one` (stitch, second-pass prune, majority correction, alts), in the FASTA, summary or bubbles writers, or in `--fetch-only` counts and exit codes passes `cargo test` and push CI. Users find it after release.

**Evidence.** On an offline 3-accession `matching` fixture, the suite's `ga_fast` gave `rc=1 ... 0/3 accessions reconstructed`. A stand-in that writes one GAF row gave `rc=0 ... 3/3 accessions reconstructed`, with `r.fasta` records equal to the query. `--fetch-only` reported `3 fetched, 0 not available, 0 failed` with rc 0.

**Proposed fix.** In `Fixture::build`, add a third stand-in `ga_ok` that reads the record name from `-f` and writes this line to the `-a` path:

```
name\t63\t0\t63\t+\t>0\t63\t0\t63\t63\t63\t60
```

Then add two tests. `a_matching_accession_is_reconstructed` runs with `-w 3` and default flags, without `--allow-empty`. It asserts exit 0, "3/3 accessions reconstructed" on stdout, `r.fasta` holding exactly `SRRFAKE{0,1,2}_primary` equal to QUERY, and three `"status": "success"` entries in the summary. `fetch_only_counts_hits_and_misses` asserts that a hit gives 0 and "2 fetched, 0 not available", and that an offline miss gives a non-zero exit. With both added, the CLI suite passed 12 of 12. Follow-ups: a two-branch bubble fixture (abundances 60/40) checking the alt record, `bubbles.json` and `snvs.tsv`, and a mixed success/error run. Also move the survivor predicate (`pipeline.rs:744-747`, retyped in the test at `:1076-1080`) into one function that both production and the test call (see TW-06).

**Tests and docs.** Update the header comment at `tests/cli.rs:1-4`. Optionally note in `docs/walker.md:284` that the golden gate covers only legacy flags. The stand-in hard-codes 63 nt and node `>0`, so it must change when QUERY or the fixture changes.

### TW-04 · medium · No test checks that each walker flag sets the field it names, or that unknown flags are refused

`walker/src/main.rs:695` · test-coverage · confirmed (reproduced) · fix effort: small · sources: SWT-01

**Problem.** `parse_args` (`main.rs:159-311`) is a hand-written `match` over about 40 flag spellings. The only unit test that calls it, `worker_and_thread_counts_must_be_at_least_one` (`main.rs:695`), checks `-w`, `--fetch-workers` and `-t`. The CLI tests pass only a handful of path flags. Nothing checks `--min-vaf`, `--min-cov`, `--[no-]prefer-majority`, `--emit-vcf`, `--max-size-mb` and the rest. Nothing checks the arm `other => return Err(... "unknown flag {other}")` (`main.rs:307`). The parser is correct today.

**Impact.** A refactor of `parse_args`, which the 0.2 flags will require, can swap or drop a mapping. Users would then get a different reconstruction, or no VCF, with exit 0 and green CI.

**Evidence.** A mutant parser swapped `--min-vaf`/`--min-cov`, inverted `--prefer-majority`, ignored `--emit-vcf`, scaled `--max-size-mb` by 1024 and skipped unknown flags. It passed all unit and CLI tests. Given a misspelled `--emit-vfc`, the mutant exited 0 with `"min_cov":0.3,"prefer_majority":true,"min_vaf":0.5` in the summary. The current binary exits 1 with `unknown flag --emit-vfc`.

**Proposed fix.** Add two unit tests in `main.rs` `mod tests`:
- `each_flag_sets_exactly_the_field_it_names` renders every flag-settable field through `snap(&Args)`. For each of about 40 cases (both spellings, each boolean in both directions, `--max-size-mb 1` → `1048576`), it asserts that exactly the expected `(field, value)` pairs differ from `parse_args(&[])`. It also asserts that out-of-range values such as `--prune-frac 1.0` are refused.
- `every_flag_in_usage_is_accepted_and_unknown_ones_are_refused` scans `USAGE` for `--long` tokens, asserts that none gives "unknown flag", and asserts that `--emit-vfc` does.

Both fail on the mutant and pass on the current sources.

**Tests and docs.** The snapshot list must grow with each new flag. No doc change.

### TW-05 · medium · The input-order test sorts its own vector; removing `run`'s sort passes every test

`walker/src/pipeline.rs:961` · test-defect · confirmed (reproduced) · fix effort: small · sources: SWT-03

**Problem.** The doc comment of `sorting_on_index_restores_input_order` says that "`run` sorts on `index` before returning" (`pipeline.rs:939-941`). The test body scrambles a `Vec` and calls `results.sort_by_key(|r| r.index);` itself (`pipeline.rs:961`). It never calls `run`, so it tests `slice::sort_by_key`. The sort that matters is at `pipeline.rs:349`. No CLI test checks output order with more than one accession.

**Impact.** If the sort at `:349` is dropped or moved, `summary.json`, the FASTA, the VCFs and `snvs.tsv` follow completion order, which changes from run to run. `--sieve` would then keep whichever duplicate finished first. Reproducibility is lost silently.

**Evidence.** With `pipeline.rs:349` deleted, all unit and CLI tests passed. On an 8-accession offline fixture (`-w 4 --fetch-workers 4`), the summary order was `SRRFAKE0 4 5 6 7 1 3 2`. The current binary gives `0..7`.

**Proposed fix.** Add a CLI test that goes through `run`. Use `Fixture::matching("order", 4)` and a stand-in aligner with `case "$*" in *SRRFAKE0*) sleep 2;; esac`. The scratch GFA path contains the accession (`extract.rs:295`), so SRRFAKE0 finishes last. Run it with `-w 4 --fetch-workers 4 --allow-empty`. Assert from the progress lines that SRRFAKE0 was reported last, and assert that the `"acc"` entries in `q_walker_summary.json` are in 0..3 order. The test fails without the sort and passed 4 of 4 runs with it (about 2 s each). Then either delete the unit test, or move the sort into `fn restore_input_order(&mut Vec<AccessionResult>)`, which `run` calls and the unit test exercises. The input-order assertion proposed in TW-03 is not enough on its own, because completion order can already match input order.

**Tests and docs.** Remove or rewrite the unit test and its doc comment.

### TW-06 · low · Phase and extract tests do not pin strands, ordering, the fork filter, the prune-floor tie or the edge guards

`walker/src/stages/phase.rs:450` · test-coverage · confirmed (reproduced) · fix effort: small · sources: R3-08

**Problem.** No phase test reads `Branch::strand`. No fixture exercises `c != fork as u32` (`phase.rs:117`). `output_is_deterministic` (`phase.rs:450`) uses a one-bubble fixture, so it cannot detect ordering changes. `prune_nodes` (`>= floor`, `extract.rs:178`) has no tie case, although `pipeline.rs:744-748` copies the predicate. The out-of-range edge guards in `expand` (`extract.rs:95-100`) and `write_gfa` (`extract.rs:332-333`) have no malformed-input test. A regression in any of these reaches `bubbles.json` or the prune decision, or turns a clean `Err` into an index panic. A combined mutant that breaks all five passed every unit and CLI test.

**Proposed fix.** Add tests:
- a two-fork fixture asserted against a literal `(fork, strand, merge, strand, nodes)` list;
- a reverse-strand diamond;
- a tandem `0+ -> {1+,2+} -> 0+` fixture with no bubble;
- a unitig exactly at `floor`, through one shared `extract::survives(graph, id, floor)` that `pipeline.rs` also calls;
- an edge to `L:+:99:+` in a 2-unitig graph, so that `expand` returns `Err(Error::Parse)` and `write_gfa` skips the L line.

### TW-07 · low · `parse_gaf` is dead code; seven parser tests pin it instead of the live `parse_gaf_multi`

`walker/src/stages/align.rs:326` · test-defect · confirmed (reproduced) · fix effort: small · sources: R4-10, RS-02

**Problem.** `align::run` returns `parse_gaf_multi(&bytes, &expects)` (`align.rs:246`). `pub fn parse_gaf` (`align.rs:326`) is called only by seven tests. The two functions differ: `parse_gaf` matches the whole name and falls back to every row, while `parse_gaf_multi` matches the first word and drops other rows (`align.rs:282-285`). `unknown_query_name_still_parses` therefore asserts behaviour the walker does not have. The crate-wide `#![allow(dead_code, unused_variables)]` (`main.rs:24`) hides the warning, which appears as `function 'parse_gaf' is never used` once the allow is narrowed.

**Proposed fix.** Delete `parse_gaf` (`align.rs:322-390`) and the "same logic as [`parse_gaf`]" sentence (`align.rs:253`). Point the tests at the live parser:

```rust
fn one(gaf: &[u8], name: &str) -> Result<Option<GafPath>> {
    Ok(parse_gaf_multi(gaf, &[name])?.remove(name))
}
```

Rewrite the unknown-name test as `unknown_query_name_is_dropped`, which asserts `None`. Add stand-in tests for `run` with exit 1, `kill -ABRT $$`, and exit 0 with no GAF. Replace the crate-wide allow with item-level allows. Doing that surfaces about ten other dead items, which must be resolved in the same change under `clippy -D warnings`.

### TW-08 · low · `LOGAN_KEEP_SCRATCH` in the environment fails four walker tests

`walker/tests/cli.rs:67` · test-defect · confirmed (reproduced) · fix effort: trivial · sources: SWT-06

**Problem.** The code under test reads `LOGAN_KEEP_SCRATCH` from the environment at each use (`scratch.rs:39`, `align.rs:241`, `types.rs:208`, `pipeline.rs:855`), and any value, even `0`, turns cleanup off. The CLI fixture's `Command::new(BIN)` (`tests/cli.rs:67`) does not remove the variable. A developer who follows `docs/walker.md:268` and exports it sees four cleanup tests fail on a correct build. `LOGAN_KEEP_SCRATCH=0 cargo test` reproduces it: 2 unit and 2 CLI failures.

**Proposed fix.**

```rust
// scratch.rs: one reader, used at all four sites
pub fn keep_scratch() -> bool { !cfg!(test) && std::env::var_os("LOGAN_KEEP_SCRATCH").is_some() }
// tests/cli.rs cmd_ga
let mut c = Command::new(BIN); c.env_remove("LOGAN_KEEP_SCRATCH");
```

Note in `docs/walker.md:268` that `cargo test` ignores the variable. Optionally add one CLI test that sets it and asserts that the directory is kept.

### TW-09 · low · `size_cap_tests` fail when a proxy variable is set without `NO_PROXY` covering 127.0.0.1

`walker/src/stages/stream.rs:413` · test-defect · confirmed (reproduced) · fix effort: trivial · sources: R1-14

**Problem.** The two loopback-server tests build their client with `Fetcher::new` (`stream.rs:413-415`). That function leaves ureq's default `Proxy::try_from_env()` in place, so with `HTTPS_PROXY`, `HTTP_PROXY` or `ALL_PROXY` set and no `NO_PROXY` entry for loopback, the request goes to the proxy. `env -u NO_PROXY -u no_proxy HTTPS_PROXY=http://127.0.0.1:9 cargo test size_cap_tests` fails both tests with "refused at the cap", which reads like a size-cap regression. Production behaviour (honouring the proxy) is correct.

**Proposed fix.** Use a test-only agent with no proxy:

```rust
fn fetcher(cap: u64) -> Fetcher {
    let mut f = Fetcher::new(&Config { max_compressed_bytes: cap, ..Config::default() });
    f.agent = ureq::Agent::new_with_config(
        ureq::Agent::config_builder().proxy(None).http_status_as_error(false).build());
    f
}
```

Optionally document in `docs/walker.md` that the walker honours `HTTPS_PROXY` and `NO_PROXY`.

### TW-10 · low · Walker tests leak empty temp directories and use the real `/dev/shm` scratch base

`walker/src/stages/extract.rs:742` · test-defect · confirmed (reproduced) · fix effort: small · sources: T1-11, RS-09

**Problem.** `tmp_cfg` (`extract.rs:742-746`) and phase's `cfg()` (`phase.rs:221-225`) point `scratch_dir` at per-PID directories under `$TMPDIR` that are never removed. Every `cargo test` leaves seven empty directories behind. `cmd_ga` (`tests/cli.rs:66-77`) passes no `--scratch-dir`, so the CLI tests use, and sweep, the machine-wide `/dev/shm`. The SIGPIPE test leaves a `logan-walker.<pid>.*` directory there, because `scratch.rs:72` handles only SIGINT, SIGTERM and SIGHUP. The product has the same gap: `logan-walker ... | head` leaves its scratch directory in `/dev/shm` until the next run sweeps it.

**Proposed fix.** Add a `#[cfg(test)] struct TestDir(PathBuf)` whose Drop calls `remove_dir_all`, or add `tempfile` as a dev-dependency, and hold one in each test that sets `scratch_dir`. In `cmd_ga`, add `.arg("--scratch-dir").arg(self.dir.join("scratch"))` before `.args(extra)`. For the product, as a separate change, keep SIGPIPE ignored and treat `ErrorKind::BrokenPipe` on stdout as an orderly stop so that `Scratch::drop` runs. Registering SIGPIPE with signal-hook would make `println!` panic instead. Then assert in the SIGPIPE test that the scratch base is empty.

### TW-11 · low · Two tests in the data-backed CI step pass when their inputs are broken or absent

`walker/src/io/input.rs:601` · test-defect · confirmed (reproduced) · fix effort: trivial · sources: SWT-05

**Problem.** The manual golden job fails on any `SKIP:` line (`ci.yml:119-122`, `:139`). Two tests in that run get past this guard:
- `surveys_real_query_fastas` logs bad FASTAs (`Err(err) => eprintln!(...)`, `input.rs:601`) and asserts only `n > 0`. An empty or malformed FASTA under `queries/` therefore still passes.
- `phase_diagnostic::shows_why_the_count_drops` returns without printing anything when `LOGAN_TEST_DIR` is unset or its file is missing (`phase.rs:555`, `:560`), and it asserts nothing.

`recall_study_kmers` (`anchor.rs:640`) needs only repository files but is `#[ignore]`d, so it never runs on push.

**Proposed fix.** In the survey, collect the errors and `assert!(errors.is_empty(), "unreadable query FASTAs: {errors:#?}")`, then drop `#[ignore]`. Add `--skip shows_why_the_count_drops` next to `--skip cross_check` in `ci.yml`, or print `SKIP:` on both early exits. Un-ignore `recall_study_kmers` and make a missing query file panic. Update the `ci.yml` comment to list every skipped diagnostic.

### TW-12 · low · The recall-study subgraph test passes on a corrupt fixture and reports it as "no data"

`walker/src/stages/extract.rs:566` · test-defect · confirmed (reproduced) · fix effort: trivial · sources: SWR3-02

**Problem.** `recall_study_subgraphs::load` turns read, decode and parse failures into `None` (`extract.rs:566`, `:569`, `:570`). The caller then prints `SKIP: no data for {acc}` and continues (`:591-593`). A truncated or non-zstd `ERR2028820.unitigs.fa.zst` therefore gives `ok` locally with a false reason. The sibling tests (for example `phase.rs:496-502`) panic instead. In CI the `SKIP:` guard still fails the job, but its message points at missing data.

**Proposed fix.**

```rust
let compressed = match std::fs::read(&p) {
    Ok(b) => b,
    Err(e) if e.kind() == std::io::ErrorKind::NotFound => return None,
    Err(e) => panic!("{p}: {e}"),
};
let rd = stream::decoder(&f).unwrap_or_else(|e| panic!("{p}: {e}"));
Some(parse::build(&f.accession, rd).unwrap_or_else(|e| panic!("{p}: {e}")))
```

Add a non-ignored `#[should_panic]` test on 24 non-zstd bytes.

### TW-13 · low · `emits_gfa_for_cross_check` deletes the GFA it reports having written

`walker/src/stages/extract.rs:884` · test-defect · confirmed (reproduced) · fix effort: trivial · sources: R3-07

**Problem.** The ignored helper exists to emit a GFA and an anchor list for diffing against another implementation (`extract.rs:835-841`). It holds `let sub = write_gfa(...)` (`extract.rs:884`) and prints `wrote <gfa_path>`. When `sub` drops, `Drop for Subgraph` unlinks the file, so `LOGAN_GFA_OUT` ends up holding only `anchors.txt`. The file name also embeds the PID and a hard-coded `ERR2028820`. This is distinct from TW-01: here the walker did create the file, so the ownership fix does not help.

**Proposed fix.** Copy the GFA to a stable name before `sub` drops:

```rust
let stable = PathBuf::from(&outdir).join("subgraph.gfa");
std::fs::copy(&sub.gfa_path, &stable).unwrap();
assert!(stable.exists());
```

Name `DIR/subgraph.gfa` and `DIR/anchors.txt` in the doc comment.
