# Walker: pipeline, command line, memory model and exit status (pipeline.rs, main.rs, config.rs, tests/cli.rs)

The per-accession pipeline in `logan-walker` is sound in its normal path: fetch, parse, align, stitch and prune produce the expected primaries, and the "never a regression" rules for the prune pass hold for the outcomes the code anticipated. The weak points sit at the edges of a run. The prune pass can still lose or fail to correct a haplotype when the second alignment errors or coverage is uneven; the exit status, stderr and summary do not tell a scheduler that most of a shard failed; all results are held in memory and written only after the last accession; and the printed memory "worst case" is neither a bound nor checked against the host. The command line has several sharp edges (a summary path built from raw FASTA ids, a fixed 900 s fetch deadline, `--max-size-mb 0`, multi-anchor runs over unrelated genes). Most fixes share one theme: make the run's state explicit and durable as it happens (stream results to disk, count errors and skips, enforce one memory budget) instead of reconstructing it at the end.

| ID | Severity | Where | Title | Fix effort |
|----|----------|-------|-------|------------|
| WP-01 | medium | `walker/src/pipeline.rs:729` | Global prune floor removes thin gene stretches; prune is refused and reference-pull errors stay | small |
| WP-02 | medium | `walker/src/pipeline.rs:667` | Multi-anchor mode keeps one primary even when the `-q` references are different genes | small |
| WP-03 | medium | `walker/src/main.rs:405` | Exit 0 whenever one accession succeeds, however many errored | small |
| WP-04 | medium | `walker/src/pipeline.rs:792` | An error in the prune pass throws away the good first-pass haplotype | small |
| WP-05 | medium | `walker/src/main.rs:578` | Summary file name embeds the raw first record id; `/` nests it, `..` escapes `--out-dir` | trivial |
| WP-06 | medium | `walker/src/main.rs:557` | Outputs are written only after the last accession; a bad `-o`, a kill or a stdout error loses the run | medium |
| WP-07 | medium | `walker/src/stages/stream.rs:49` | Fixed 900 s whole-transfer fetch deadline, no flag, retried from byte 0 | small |
| WP-08 | medium | `walker/src/config.rs:445` | Printed "worst case" memory bound is an unenforced estimate that real runs exceed | medium |
| WP-09 | medium | `walker/src/pipeline.rs:335` | Every per-accession result is kept in memory until the run ends | medium |
| WP-10 | medium | `walker/src/main.rs:497` | Memory bound is only printed, never checked against available RAM, and wrong for `--fetch-only` and short lists | medium |
| WP-11 | low | `walker/src/main.rs:532` | A GraphAligner known to clip the last base (1.0.19+) only draws a warning | small |
| WP-12 | low | `walker/src/main.rs:562` | Earlier runs' `snvs.tsv`, `alt_haplotypes.fasta` and VCFs survive in a reused `--out-dir` | small |
| WP-13 | low | `walker/src/pipeline.rs:416` | Compute-stage failures reach stdout only; stderr says "see the log above" | trivial |
| WP-14 | low | `walker/src/error.rs:82` | A run where every accession was skipped as `too_large` is reported as "nothing was found" | small |
| WP-15 | low | `walker/src/pipeline.rs:640` | Single-anchor `no_path` drops the rejected alignment; log and summary cannot say why | small |
| WP-16 | low | `walker/src/pipeline.rs:655` | Summary candidates below `--min-cov` report null `cov`/`length`/`nodes` | small |
| WP-17 | low | `walker/src/main.rs:211` | `--max-size-mb 0` is accepted as a 0-byte cap that refuses every accession | trivial |
| WP-18 | low | `walker/src/pipeline.rs:438` | A panic in a fetch worker aborts the whole run with exit 101 | small |
| WP-19 | low | `walker/src/config.rs:471` | Default of 8 compute workers oversubscribes hosts with fewer cores | small |
| WP-20 | low | `walker/src/pipeline.rs:281` | Hard-wired 4-slot fetch queue adds four payloads to peak memory | trivial |
| WP-21 | low | `walker/src/pipeline.rs:580` | Help text, README, comments and docs disagree with the code in about fifteen places | small |
| WP-22 | low | `walker/src/pipeline.rs:273` | `pipeline::run` doc says downstream stages get only the first reference | trivial |
| WP-23 | low | `walker/src/main.rs:24` | Crate-wide `allow(dead_code)` hides nine dead or test-only items | small |
| WP-24 | low | `walker/src/pipeline.rs:743` | Prune survivor predicate is written three times, kept in step by a comment | trivial |
| WP-25 | info | `walker/src/pipeline.rs:32` | Graph-peak ratio comments give 1.5x and 1.6x where the constant is 1.7 | trivial |
| WP-26 | low | `walker/src/main.rs:331` | A closed stdout (SIGPIPE) skips scratch cleanup and leaves GraphAligner children running | small |
| WP-27 | low | `walker/src/main.rs:197` | `--min-cov`, `--min-vaf`, `--min-abund`, `--min-shared-kmers` accept NaN, inf and out-of-range values | trivial |
| WP-28 | low | `walker/src/main.rs:179` | The argument parser keeps only the last repeated `-a`/`-l` and takes the next option as a value | small |

### WP-01 · medium · Global prune floor removes thin gene stretches; prune is refused and reference-pull errors stay

`walker/src/pipeline.rs:729` · science · confirmed (reproduced) · fix effort: small · sources: SC1-06

**Problem.** The prune floor is one number for the whole gene (`let floor = stages::extract::depth_floor(cfg, depth);`, `walker/src/pipeline.rs:729-730`), and `prune_nodes` drops every subgraph unitig below it (`walker/src/stages/extract.rs:178`), including route unitigs with no deeper alternative. A gene stretch genuinely covered below 1% of the depth disappears, the realignment cannot span the reference, the prune is refused (`walker/src/pipeline.rs:800`) and the reference-pulled first-pass route stands. `prefer_majority` skips such sites (`walker/src/stages/variants.rs:188`). Separately, `path_depth` is an unweighted median over steps (`ab[ab.len() / 2]`, `walker/src/stages/extract.rs:133`), so short error unitigs can drag the floor to `--min-abund`.

**Impact.** A single-strain sample at 10,000x with reference-matching error branches at 30x at 3 variant sites and its last 200 nt at 60x (amplicon dropout, 3' fall-off) gets a primary carrying 3 reference alleles it does not have, `provisional: false`, exit 0; the only signal is `prune.outcome = refused_*`. Uneven coverage is routine in ARTIC amplicon data.

**Evidence.** Simulated 900-nt gene, shipped binary, stand-in max-identity aligner: even coverage gives a primary matching the strain at 900/900; with the 60x tail the log says `floor 100.0x cost the alignment; keeping unpruned` and the primary matches at 897/900.

**Proposed fix.** In `walker/src/stages/extract.rs` and `walker/src/pipeline.rs`:
1. `path_depth`: k-mer-weighted median, weighting each step by `seq_len - 30`.
2. New `protected_route_nodes(graph, steps, floor) -> HashSet<UnitigId>`: pass-1 route unitigs below the floor with no sibling at or above the floor leaving the predecessor's oriented end or entering the successor (a thin stretch of the gene, not an error branch).
3. `prune_nodes(graph, nodes, seqs, floor, &protected)` keeps protected nodes; the survivor count at `walker/src/pipeline.rs:743-748` must use the same predicate (see WP-24).
4. When the outcome is `refused_*`, set `provisional = true` or add a `bias_uncorrected` flag.

Do not relax `prefer_majority`'s precondition instead (its doc comment records the all-forks variant as unstable). A patched build gave `applied` and 900/900 in all simulated cases and passed all 221 unit tests.

**Tests and docs.** Unit tests for `protected_route_nodes` (thin tail protected, error branch with a deep sibling not) and the weighted `path_depth`; a CLI case asserting `applied`. Recheck the figures at `walker/src/pipeline.rs:699-719` with `bench/check_golden.sh`. `docs/walker.md:129-131,204`: a `refused_*` primary is not bias-corrected.

### WP-02 · medium · Multi-anchor mode keeps one primary even when the `-q` references are different genes

`walker/src/pipeline.rs:667` · correctness · confirmed (reproduced) · fix effort: small · sources: SWR6-01

**Problem.** Multi-anchor mode treats every record across all `-q` files as a strain of one gene: `let winner = stages::stitch::pick_primary(candidates, cfg.primary_metric);` (`walker/src/pipeline.rs:667`) ranks all references together and keeps one. Nothing checks that they are related. The `-q` help (`walker/src/main.rs:55-58`) only says two or more references switch on multi-anchor mode, and `README.md:74-75` lists the 160-record, 25-gene `queries/master_planetary_targets.fasta` among the references to start from. The "strains of the same gene" restriction appears only in `docs/walker.md:53`. Version 0.1.0 printed "One gene per query" for a multi-record `-q`; 0.1.1 removed the warning without a replacement.

**Impact.** `-q geneA.fa -q geneB.fa`, or the master panel, on an accession carrying several genes: one `<ACC>_primary ... anchor=<winner>`, exit 0, no warning. The other genes appear only as unselected `candidates` in the summary JSON. A resistance-gene panel (bla_KPC plus bla_NDM) yields one gene per sample.

**Evidence.** Two unrelated CDSs (90 and 120 nt) in one cached accession, stand-in GraphAligner returning full-length rows, `-q geneA.fa -q geneB.fa`: `DONE`, one record `anchor=geneA`, geneB `selected:false` with cov 1, pident 100. A k-mer component check (at least 1% shared canonical 21-mers) splits the master panel into 49 groups.

**Proposed fix.** After the multi-anchor notice in `walker/src/main.rs`, group references by shared canonical 21-mers (union-find, link when shared at least 1% of the smaller set) and warn when there is more than one group:
```
[warn] multi-anchor: the 2 references fall into 2 groups that share under 1% of their 21-mers ([geneA] [geneB]); one primary is kept per accession
```
Warn, do not refuse: divergent but legitimate panels also split. Later, an opt-in `--per-group-primary` runs `pick_primary` per group (`<ACC>_primary_<group>`). A patched build warned for the two-gene case, stayed silent for the spike strain panel, and passed all 221 unit tests.

**Tests and docs.** Unit test for the grouping (two random sequences give 2 groups; a sequence and a 1%-mutated copy give 1). CLI test asserting the warning on stderr. `-q` help, `README.md:74-75`, `docs/walker.md:51-53` and the Python `-q` help: references must be strains of one gene; describe the master panel as a catalogue, not a `-q` panel.

### WP-03 · medium · Exit 0 whenever one accession succeeds, however many errored

`walker/src/main.rs:405` · reliability · confirmed (reproduced) · fix effort: small · sources: RA-03, R6-01

**Problem.** `check_run` returns success as soon as one accession was reconstructed: `if total == 0 || success > 0 { return Ok(()); }` (`walker/src/main.rs:405`). The closing line prints only `"{}/{} accessions reconstructed"` (`walker/src/main.rs:583-586`), and the summary JSON carries `total` and `success` but no error count. The rule is inconsistent with its neighbours: one `no_anchors` plus one error exits non-zero even with `--allow-empty`, and `--fetch-only` exits 1 on any error. It is documented (`docs/walker.md:304`) and pinned by the test `one_success_is_enough_whatever_else_happened` (`walker/src/main.rs:639`).

**Impact.** A shard that lost 9,999 of 10,000 accessions to an S3/DNS outage, damaged cache files or GraphAligner OOM kills exits 0. Snakemake, Nextflow and Galaxy (`detect_errors="exit_code"` in the tools-iuc PR #8800 wrapper) record a good shard and never retry; the cohort shrinks without notice. The project's own analysis lost 8 of 187 accessions to transient DNS failures (`docs/issue_intra_accession_variance_lost.md:216`).

**Evidence.** One good cached accession plus nine truncated `.zst` files, `--offline`, stand-in GraphAligner: nine `ERROR: ... incomplete frame` lines, `1/10 accessions reconstructed`, exit 0, summary `{"total":10,"success":1}`.

**Proposed fix.** Non-breaking for 0.1.x:
```rust
// main.rs: --max-errors <N>  -> args.max_errors: Option<usize>
// error.rs: Error::TooManyErrored { total, success, errors, max }  => exit_code 4
fn check_run(results, allow_empty, max_errors) -> Result<()> {
    if total == 0 { return Ok(()) }
    if errors == total { return Err(AllErrored { .. }) }            // 3
    if let Some(m) = max_errors { if errors > m { return Err(TooManyErrored { .. }) } } // 4
    if success > 0 { return Ok(()) }
    /* existing nothing-found logic */
}
```
Always print `"{s}/{t} accessions reconstructed, {n} found nothing or skipped, {e} errored"` and add `errors` and `negatives` to the summary top level. For 0.2, consider failing with exit 4 on any error by default and mapping 4 to a warning in the Galaxy `<stdio>` block. A patched build gave exit 4 with `--max-errors 0` and exit 0 with `--max-errors 9`.

**Tests and docs.** Unit test `max_errors_fails_a_run_that_lost_too_many_accessions`; CLI test with one good and one corrupt cached accession and `--max-errors 0` expecting exit 4 and outputs present. Exit table in `docs/walker.md:304-312`, USAGE, `error.rs` exit-code doc.

### WP-04 · medium · An error in the prune pass throws away the good first-pass haplotype

`walker/src/pipeline.rs:792` · reliability · confirmed (reproduced) · fix effort: small · sources: R3-02, RE-07

**Problem.** The default second pass (`--prune-frac 0.01`) drops the first-pass subgraph (`drop(subgraph);`, `walker/src/pipeline.rs:787`) and then runs `write_gfa(...)?`, `stages::align::run(cfg, &pruned, ...)?` (`:792`) and `stitch(...)?`. With `?`, any error in this pass fails the whole accession although a valid haplotype is already held in `haplotype`/`final_path`. That breaks the contract stated at `:721-725` ("Pruning can improve an accession and cannot break one") and in `docs/issue_intra_accession_variance_lost.md` §5.1. Two side defects: the comment at `:780-781` ("Both subgraphs live at the same scratch path") is false, since `scratch_stem` (`walker/src/stages/extract.rs:256-261`) names every call uniquely; and the refusal branch (`:840-845`) re-runs the BFS, decode and GFA write, all with `?`, to rebuild a file nothing reads.

**Impact.** GraphAligner aborting or being OOM-killed on the pruned graph, or tmpfs scratch filling during the second write, turns a reconstructed accession into `ERROR`: missing from the FASTA, `prune: null` in the summary, exit 3 for a one-accession run. The trigger is uncommon but sits on the default path for every accession with sub-floor unitigs.

**Evidence.** Two-unitig fixture (abundance 100 and 5), a stand-in GraphAligner that exits 134 only when its `-a` path ends in `.pass2.gaf`, `--prune-frac 0.5`: `ERROR: ... GraphAligner failed: exited with exit status: 134`, exit 3, empty FASTA. With `--prune-frac 0`: exit 0 and a 63-nt primary.

**Proposed fix.** In the `PruneOutcome::Prune` arm:
1. Keep `subgraph`; delete only its GFA file (unless `LOGAN_KEEP_SCRATCH` is set) and prune from clones of `nodes`/`seqs` (or make `prune_nodes` take slices).
2. Run write_gfa + align + stitch in a closure returning `Result<(Subgraph, Option<(GafPath, Haplotype)>)>`; on `Err(e)` log `prune: realignment failed (e); keeping the unpruned result` and record a new `PruneResult::RefusedError` (`"refused_error"`).
3. Delete the rebuild at `:840-845`, since pass 1 is still alive.
4. Add `RefusedError` to `types.rs` `as_str` and `main.rs` `decision_lines`; fix the comment at `:780-781`.

A patched build (about 60 changed lines) printed `prune: realignment failed (...)`, wrote the 63-nt primary with `outcome: refused_error`, exited 0, left scratch empty, and passed all 10 CLI tests.

**Tests and docs.** CLI test with the pass-2-failing stand-in asserting exit 0, a primary and `"outcome": "refused_error"`. `docs/walker.md:129-131`: add `refused_error`.

### WP-05 · medium · Summary file name embeds the raw first record id; `/` nests it, `..` escapes `--out-dir`

`walker/src/main.rs:578` · security · confirmed (reproduced) · fix effort: trivial · sources: RX-01, R5-03, R6-03, R2-02, G1-09, IC-03, R1-09, RA-04, RC-05, D1-07

**Problem.** The summary path is `args.out_dir.join(format!("{}_walker_summary.json", queries[0].name))` (`walker/src/main.rs:578`). For a multi-record `-q` file, `Query::name` is `format!("{stem}__{record_id}")` (`walker/src/io/input.rs:295`), where `record_id` is the first header word, unvalidated. `write_summary` then calls `create_dir_all(parent)` (`walker/src/io/output.rs:387`), so header text chooses the directory. Accessions are validated for path use (`walker/src/io/input.rs:96-104`); record ids are not. `docs/walker.md:108` documents `<out-dir>/<query stem>_walker_summary.json`, which is wrong for every multi-record reference.

**Impact.** Standard strain ids (`>A/California/07/2009(H1N1)`, `>hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124`) nest the summary: the walker exits 0, `out/*_walker_summary.json` matches nothing, and the tools-iuc PR #8800 wrapper's final `mv` fails the Galaxy job after the whole run. A first id such as `x/../../../NAME` creates or overwrites `NAME_walker_summary.json` outside `--out-dir`, at any path the job user can write (limited to that suffix; no code execution). A 300-character id fails the run at the end with ENAMETOOLONG. Single-record references are unaffected; FASTA, VCF and TSV outputs are correct.

**Evidence.** One cached accession, stand-in GraphAligner, `-q panel.fasta` with first header `>A/Puerto_Rico/8/1934`: exit 0, file at `out/panel__A/Puerto_Rico/8/1934_walker_summary.json`, `ls out/*_walker_summary.json` fails. With first id `x/../../../ESCAPED`, run from `job/work` with `-o out`, the file lands at `job/ESCAPED_walker_summary.json`. An id of 40 `../` plus an absolute path overwrote an existing summary there.

**Proposed fix.** Name the summary after the first `-q` file's stem in every mode; the record ids are already in the summary's `references` array:
```diff
-    let summary_path = args.out_dir.join(format!("{}_walker_summary.json", queries[0].name));
+    let stem = args.query[0].file_stem().map(|s| s.to_string_lossy().into_owned())
+        .ok_or_else(|| error::Error::Usage(format!("{}: cannot derive a name", args.query[0].display())))?;
+    let summary_path = args.out_dir.join(format!("{stem}_walker_summary.json"));
```
Also drop `create_dir_all(parent)` from `write_summary` (`main` already creates `--out-dir`). Keep `Query::name` as it is; it labels anchor messages and the fallback VCF CHROM. Alternative that keeps today's names for safe ids: map every character outside `[A-Za-z0-9._-]` in the name to `_`. Release note: multi-record runs change name from `panel__<id>_...` to `panel_...`; the wrapper can then `mv 'walker_out/${reference_name}_walker_summary.json'` instead of globbing.

**Tests and docs.** CLI test `the_summary_name_never_contains_a_record_id` over ids `refA`, a GISAID id, `../../../ESCAPED` and `/../../../ESCAPED`, asserting exactly one `<out>/panel_walker_summary.json`. `docs/walker.md:108` ("first `-q` file's name minus its last extension; `q.fasta.gz` gives `q.fasta`"), `walker/src/types.rs:32-33`, `walker/src/io/input.rs:216-218`.

### WP-06 · medium · Outputs are written only after the last accession; a bad `-o`, a kill or a stdout error loses the run

`walker/src/main.rs:557` · reliability · confirmed (reproduced) · fix effort: medium · sources: R6-04

**Problem.** `main` runs the whole cohort (`let summary = pipeline::run(&args.cfg, accessions, &queries)?;`, `walker/src/main.rs:554`) and only then touches the output location (`std::fs::create_dir_all(&args.out_dir)...`, `:557`). Results are held in memory until then (WP-09). So (a) an unusable `-o` or `-f` is reported only after every accession was fetched and aligned; (b) any termination before the end (scheduler walltime SIGTERM, Ctrl-C, OOM kill) leaves no FASTA, no summary and nothing to resume from; (c) progress uses `println!` (`walker/src/pipeline.rs:394-416`, `:218`), which panics on a write error, so a failing stdout aborts with exit 101 and no outputs. The comment at `walker/src/pipeline.rs:324-326` acknowledges (c).

**Impact.** Hours of fetching and alignment lost: (b) is likely for long HPC shards; (a) occurs whenever `-o` names a regular file or `-f` (default `reconstructed_haplotypes.fasta` in the CWD, not in `-o`) has an unwritable parent.

**Evidence.** Four cached accessions, stand-in GraphAligner sleeping 2 s, `-w 1`: `-o afile` (a regular file) prints four `DONE` lines, then after 8 s `error: cannot create .../afile: File exists`, exit 1. SIGTERM after `[3/4]`: exit 143, no output directory. `> /dev/full`: `failed printing to stdout: No space left on device`, exit 101, no outputs.

**Proposed fix.**
1. (trivial) Validate outputs before the run:
```rust
std::fs::create_dir_all(&args.out_dir).map_err(error::file_err("create", &args.out_dir))?;
std::fs::OpenOptions::new().create(true).append(true).open(&args.out_fasta)
    .map_err(error::file_err("open", &args.out_fasta))?;
let summary = pipeline::run(&args.cfg, accessions, &queries)?;
```
2. (trivial) Make progress non-fatal: `let _ = writeln!(std::io::stdout().lock(), ...)` in `report` and `run_fetch_only`.
3. (medium) Stream results in input order to the FASTA, alt FASTA and a JSON-lines journal in `--out-dir`, assemble `summary.json`, `bubbles.json`, VCFs and the SNV table from the journal at the end, write final files under temp names and rename; add `--resume` to skip accessions already journalled. This is the same change as WP-09.

With step 1 applied, a bad `-o` failed in 0.02 s before any work and a normal run wrote all outputs.

**Tests and docs.** CLI test: `-o` pointing at a regular file and a stand-in aligner that writes a marker must fail before the aligner runs. Unit test that `report` survives a closed stdout. For step 3, kill after N results and check the journal holds N records; keep output byte-identical (golden gate). `docs/walker.md` outputs section: when files are written, how to resume.

### WP-07 · medium · Fixed 900 s whole-transfer fetch deadline, no flag, retried from byte 0

`walker/src/stages/stream.rs:49` · reliability · confirmed (reproduced) · fix effort: small · sources: R6-06

**Problem.** `Fetcher::new` sets `.timeout_global(Some(Duration::from_secs(cfg.fetch_timeout_secs)))` (`walker/src/stages/stream.rs:49`) with `fetch_timeout_secs: 900` (`walker/src/config.rs:474`). In ureq 3.4 this deadline covers the whole request including the body, so a download still arriving steadily is aborted at 900 s as `Transient("body read failed: timeout: global")` (`stream.rs:216`). `fetch()` retries with the same plain GET (`stream.rs:124-125`, no `Range`), so each attempt restarts from byte 0 and hits the same wall. `parse_args` has no flag for the timeout or the attempt count. The justification (`config.rs:347-348`) uses a single-stream rate, while the default runs 8 streams sharing the link.

**Impact.** A 1 GiB file (the default cap) needs about 1.19 MB/s per stream; with 8 fetch workers that is about 76 Mbit/s aggregate. On a 50 Mbit/s link every accession above about 700 MB fails; at 20 Mbit/s, above about 270 MB. Accessions over 350 MB were about 7.5% of the panel. Each costs 3 x 900 s of a fetch worker and up to about 2 GB of wasted transfer, and is reported as an error. Workaround: lower `--fetch-workers`.

**Evidence.** Unit test against a local HTTP server trickling a 40 kB body over about 4 s with `fetch_timeout_secs: 2`: `Transient(body read failed: timeout: global)` after 2.0 s; the same body completes with a longer deadline. Retry from zero confirmed by reading `try_once` (`self.agent.get(url).call()` with no headers).

**Proposed fix.**
1. Resume on retry: keep the bytes a failed attempt received; send `Range: bytes=<n>-` with `If-Range: <ETag of the first response>`; on `206` with a matching `Content-Range` append, on `200` restart. Read through `body_mut().as_reader()` into the kept buffer so a timeout keeps progress, and count accumulated bytes against `--max-size-mb`.
2. Add `--fetch-timeout <SECS>` and `--fetch-attempts <N>` (at least 1) to `parse_args` and pass them through `pipe_reconstruct.py`/`pipe_from_kmindex.py`.
3. Scale the default, e.g. `max(900, max_compressed_bytes / 256 KiB/s)`, and correct the `config.rs:347` comment.

ureq 3 has no idle (stall) timeout, so resume plus deadline is the practical design.

**Tests and docs.** Unit test: a server that trickles and honours `Range`, `fetch_timeout_secs: 2`, `fetch_attempts: 3`, a 4 s body: content byte-identical and the second request carries `Range`. Parse tests for the new flags. `docs/walker.md` option table, `docs/codon_msa.md:90`, `--help`.

### WP-08 · medium · Printed "worst case" memory bound is an unenforced estimate that real runs exceed

`walker/src/config.rs:445` · memory · confirmed (reproduced) · fix effort: medium · sources: BL1-03, D1-05, RA-06

**Problem.** The startup line prints `... GB worst case` (`walker/src/main.rs:497-511`) from compressed buffers x cap plus `compute_workers x cap x GRAPH_PEAK_PER_COMPRESSED_BYTE` (`= 1.7`, `walker/src/config.rs:445`, measured one accession per run). `docs/walker.md:252` also calls it a worst case. Nothing enforces it: `parse::build` grows its arrays without a budget, and extract/phase structures, the GFA on tmpfs scratch (default `/dev/shm`), glibc arena retention across accessions and GraphAligner are outside the model. Several comments still say the GFA is "~38 KB" (`walker/src/types.rs:193`, `walker/src/config.rs:364`, `walker/src/stages/extract.rs:19,278`, `walker/src/stages/align.rs:11`).

**Impact.** Anyone who sizes a job from the line with the cap near accession size (the tools-iuc PR #8800 wrapper does) gets an OOM kill, not a `too_large` row. Typical runs with a 1 GB cap and small accessions stay inside it.

**Evidence.** `-w 1 --fetch-workers 1`, offline cache, stand-in GraphAligner: 8 copies of ERR2028820 with `--max-size-mb 33` printed 0.2 GB (254 MiB) and peaked at 303 MiB; 6 copies of a 73 MiB HIV accession with `queries/query_hiv1_pol.fasta` printed 0.6 GB and peaked at 838 MiB; a format-valid dense file (3 M four-base unitigs, 2.4 MiB compressed, `--max-size-mb 4`) printed 0.0 GB (about 31 MiB) and peaked at 138 MiB.

**Proposed fix.**
1. `parse::build_with_budget(acc, reader, budget)`: count `4*(offsets+targets+abundance+seq_offset capacities) + seq_2bit capacity` and return `Error::TooLarge` before a doubling would exceed `1.7 x max_compressed_bytes`; map it to `Outcome::TooLarge` in the compute worker (`walker/src/pipeline.rs:482`) and reword the line ("over --max-size-mb or graph budget").
2. In `main()` on linux-gnu: `libc::mallopt(libc::M_MMAP_THRESHOLD, 1 << 20)` (`libc` is already a dependency).
3. A node or byte cap on the subgraph with its own outcome, next to the existing 3-unitigs-per-base guard.
4. Startup line and docs: say GraphAligner and tmpfs scratch are excluded, or add about 1.2x compressed per worker when scratch is tmpfs; fix the "~38 KB" comments.

Steps 1-2 in a patched build: dense file became `too_large ... over the 6 MB graph budget` at 15.9 MiB peak; 8x ERR2028820 peaked at 243 MiB (inside 254); 221 unit tests passed. The subgraph case needs step 3. Longer term, one byte budget (WP-10 step 3).

**Tests and docs.** Unit test: `build_with_budget` on the dense fixture returns `TooLarge`; pipeline test that it becomes `Outcome::TooLarge`. `docs/walker.md` Memory section: what is and is not bounded.

### WP-09 · medium · Every per-accession result is kept in memory until the run ends

`walker/src/pipeline.rs:335` · memory · confirmed (reproduced) · fix effort: medium · sources: RF-02

**Problem.** The collector gathers every `AccessionResult` into one `Vec` (`result_rx.into_iter().inspect(..).collect::<Vec<_>>()`, `walker/src/pipeline.rs:335-341`). Each result holds the primary twice (`haplotype` and `haplotypes`, cloned at `:894,907`), all alts and all bubbles with their branches. Every writer runs after `pipeline::run` returns (`walker/src/main.rs:554-581`), and `snv_sites_all` (`walker/src/io/output.rs:444`) materialises every accession's VCF rows at once. The diagram says "results channel (unbounded; results are small)" (`walker/src/pipeline.rs:40`), and the printed bound omits this term.

**Impact.** Memory grows linearly with finished accessions on top of the printed "worst case": golden summaries range from 53 KB (ERR2028820, 202 bubbles) to 2.9 MB (SRR8209082, 11,825 bubbles), so a 20,000-accession run of deep HIV samples can retain several to about 20 GB. Nothing is on disk until the end (WP-06).

**Evidence.** A fixture with 1,176 bubbles and 32 alts per accession, symlinked N times, `-w 4 --fetch-workers 4`, stand-in GraphAligner: peak RSS 6.7 MiB (N=1), 100 MiB (250), 327 MiB (1,000), 628 MiB (2,000); about 0.30 MiB per finished accession.

**Proposed fix.**
1. Replace the collect with an in-order streaming sink: `next: usize` plus a `BTreeMap<usize, AccessionResult>` reorder buffer; drain while the map holds `next`, passing each result to a `Sink` holding the open writers (primary FASTA, alt FASTA, `bubbles.json` element, `snvs.tsv` rows, per-accession VCF) and appending the summary JSON fragment to a temp file. Keep only a small stub (index, accession, outcome, error, prune/provisional fields, bubble count) for `decision_lines`, `bubble_events`, `check_run` and `success`.
2. `write_summary` writes the header, `total`, `success`, splices the fragment file, then `bubble_events`; output stays byte-identical because order is still input order.
3. Compute SNV sites per accession in the sink.
4. (trivial, independent) Drop the duplicate `AccessionResult::haplotype`; derive it from `haplotypes.first()`.
Bound the reorder window by blocking compute workers when the map exceeds, for example, 4x workers.

**Tests and docs.** Test over N cached accessions with shuffled completion order that outputs are byte-identical to today's; unit test that the reorder buffer emits in index order; `bench/check_golden.sh`. Fix the diagram at `walker/src/pipeline.rs:40` and the `docs/walker.md` Memory section.

### WP-10 · medium · Memory bound is only printed, never checked against available RAM, and wrong for `--fetch-only` and short lists

`walker/src/main.rs:497` · memory · confirmed (reproduced) · fix effort: medium · sources: R6-07, BL1-08

**Problem.** `Config::max_resident_bytes()` and its parts only format the startup line (`walker/src/main.rs:497-511`); nothing compares it with `MemAvailable` or the cgroup limit, and nothing limits bytes in flight. At the defaults the walker prints `33.6 GB worst case` and proceeds on any machine; the remedy in `docs/walker.md:249-257` is manual. The line is also printed before the `--fetch-only` branch (`walker/src/main.rs:515`) and built from worker counts only (`peak_compressed_buffers()` = `channel_capacity + fetch_workers + compute_workers`, `walker/src/config.rs:408`), so it charges `--fetch-only` for a queue and 8 parsed graphs that `run_fetch_only` (`walker/src/pipeline.rs:186`) never creates, and charges a one-accession run for 20 buffers and 8 graphs.

**Impact.** A cohort of accessions near 350 MB at default flags needs about 11.8 GB, enough to OOM a 16 GB node, and with WP-06 nothing survives the kill. Conversely, `--fetch-only --fetch-workers 1` prints 26.6 GB where 1 GB is the bound, and `-a ACC` prints 33.6 GB where 2.7 GB is; schedulers sizing from the line over-request by up to about 25x. The PR #8800 Galaxy wrapper re-implements this arithmetic in shell with the outdated 43 MB-per-graph model.

**Evidence.** On a 15 GB host with 12 GB available, a default run printed `33.6 GB worst case` and proceeded; `grep -rn "MemAvailable\|cgroup\|memory.max" walker/src` finds nothing. `--fetch-only --fetch-workers 1` on one accession printed `13 compressed buffers (4 queued + 1 fetching + 8 parsing) ... 26.6 GB worst case`.

**Proposed fix.**
1. (trivial) Make the bound mode- and count-aware and print it after the branch is chosen:
```rust
pub fn peak_compressed_buffers_for(&self, n: usize) -> usize {
    let h = if self.fetch_only { self.fetch_workers } else { self.peak_compressed_buffers() };
    h.min(n)
}
pub fn peak_graphs_for(&self, n: usize) -> usize {
    if self.fetch_only { 0 } else { self.compute_workers.min(n) }
}
```
plus `max_resident_bytes_for(n)`; compute `n = accessions.len()` before it is moved.
2. (small) Read `/sys/fs/cgroup/memory.max`, else `MemAvailable` from `/proc/meminfo`, and print a `[warn]` naming settings that would fit when the bound exceeds it.
3. (medium) `--max-memory-mb` (default 0.8 x that limit) backed by a byte semaphore: a fetch worker acquires `len x (1 + GRAPH_PEAK_PER_COMPRESSED_BYTE)` once the size is known and before reading the body; `run_one` releases the compressed share after parsing and the rest when the graph drops; an oversize item waits to run alone. Acquire before reading the body and release on every error path to avoid deadlock with the bounded channel.

**Tests and docs.** `memory_bound_tests::bound_respects_accession_count_and_fetch_only`; semaphore tests with fake sizes (never over budget, oversize item runs alone). `docs/walker.md:249-256` and the `pipeline.rs` module doc. Tell the Galaxy wrapper authors to pass `--max-memory-mb` instead of their shell model.

### WP-11 · low · A GraphAligner known to clip the last base (1.0.19+) only draws a warning

`walker/src/main.rs:532` · provenance · confirmed (reproduced) · fix effort: small · sources: R4-07

**Problem.** The only version check is `if !version.contains("1.0.17b") { eprintln!("[warn] ... 1.0.19 and later clip the query's last base ...") }` (`walker/src/main.rs:532-537`). The project's own diagnosis (`docs/issue_graphaligner_end_clipping.md`, lines 66-69) shows 1.0.19/1.0.20 drop the last reference base of every alignment, yet such a run proceeds, passes the coverage gate and exits 0. The same text is printed for 1.0.18, which does not clip. Bioconda and `environment.yml:19` pin 1.0.17b, so only `--graphaligner` pointing at another install triggers it; then every primary is one base short (`cov=0.998`).

**Evidence.** A stand-in reporting `GraphAligner bioconda 1.0.20-` and emitting a clipped GAF: the warning, `DONE`, `>SRRFAKE0_primary cov=0.998`, exit 0.

**Proposed fix.** In `walker/src/stages/align.rs` add `version_class(banner) -> GaVersion { Validated, KnownClipping(v), Unvalidated(v) }` (extract the numeric token; 1.0.17b / at least 1.0.19 / other). In `main.rs`, refuse `KnownClipping` with a usage error unless `--allow-unvalidated-graphaligner`; give a version-specific warning for `Unvalidated`. Optionally add the reference span to the FASTA header when the GAF row does not span the reference.

**Tests and docs.** Unit test over 17b/18/19/20 banners; CLI test with a 1.0.20 stand-in expecting exit 1, and exit 0 with the override. `docs/walker.md` options table and line 40.

### WP-12 · low · Earlier runs' `snvs.tsv`, `alt_haplotypes.fasta` and VCFs survive in a reused `--out-dir`

`walker/src/main.rs:562` · provenance · confirmed (reproduced) · fix effort: small · sources: IC-11, SWD-06

**Problem.** Optional outputs are written only inside `if args.cfg.emit_vcf || args.cfg.emit_snv_table { ... }` (`walker/src/main.rs:562`) and `if args.cfg.emit_alts { ... }` (`:575`); nothing removes them when the flag is off, and `write_vcfs` removes a stale VCF only for an accession of the current run with no sites (`walker/src/io/output.rs:470-480`). A second command-line run into the same `-o` with another list and without `--emit-snv-table` leaves run 1's `snvs.tsv`, `alt_haplotypes.fasta` and VCFs beside run 2's FASTA and summary, describing accessions run 2 never processed. Galaxy jobs use fresh directories and are unaffected.

**Proposed fix.** After `create_dir_all(out_dir)`:
```rust
if !cfg.emit_snv_table { remove_if_exists(out_dir.join("snvs.tsv"))?; }
if !cfg.emit_alts      { remove_if_exists(out_dir.join("alt_haplotypes.fasta"))?; }
io::output::remove_stale_vcfs(&summary, out_dir, cfg.emit_vcf)?;
```
`remove_stale_vcfs` deletes `*.vcf` whose header contains `##source=logan-walker` and that this run will not rewrite; ignore NotFound. Optionally list written files in the summary (`"outputs": [...]`).

**Tests and docs.** CLI test: two runs into one `-o`, the second without the flags; stale files gone, a non-walker `.vcf` kept. `docs/walker.md:106-110`.

### WP-13 · low · Compute-stage failures reach stdout only; stderr says "see the log above"

`walker/src/pipeline.rs:416` · usability · confirmed (reproduced) · fix effort: trivial · sources: R6-13, RC-10

**Problem.** `report()` prints failures with `println!` (`walker/src/pipeline.rs:416`). The fetch pool also echoes them to stderr (`:446-449`), but `compute_worker` does not, so parse errors, GraphAligner failures and caught panics appear only on stdout. The final stderr error, "all N accessions failed with an error (see the log above)" (`walker/src/error.rs:94-98`), points at nothing when only stderr is read (scheduler `.err` files, Galaxy's failure view). The panic message also lacks the accession.

**Proposed fix.** In `compute_worker`, after the `catch_unwind` result:
```rust
.unwrap_or_else(|_| AccessionResult::failed(index, panic_acc.clone(),
    format!("[{panic_acc}] internal error: a stage panicked")));
if result.outcome == Outcome::Failed {
    eprintln!("{}", result.error.as_deref().unwrap_or("unknown error"));
}
```
A patched build showed `[SRR900000] GraphAligner failed: ... boom` on stderr with stdout discarded.

**Tests and docs.** CLI test with a stand-in GraphAligner exiting 7: exit 3, stderr contains "GraphAligner failed". One sentence in `docs/walker.md` on stdout versus stderr.

### WP-14 · low · A run where every accession was skipped as `too_large` is reported as "nothing was found"

`walker/src/error.rs:82` · usability · confirmed (reproduced) · fix effort: small · sources: BL1-06

**Problem.** `too_large` is classed with the "nothing found" outcomes (`walker/src/types.rs:429`), so a run where `--max-size-mb` skipped every accession ends with "0/N accessions reconstructed, and none errored: nothing was found" (`walker/src/error.rs:82-87`), although no data was read; with `--allow-empty` it exits 0 with an empty FASTA. There is no tabular per-accession status a workflow could branch on, and the compressed size exists only inside the `error` detail string. The per-accession `too_large` lines above the error do say what happened.

**Proposed fix.** Give `NothingReconstructed` a `skipped` count and print "0/2 reconstructed: 2 skipped by --max-size-mb, 0 found nothing, 0 errored" (exit codes unchanged). Write `<out-dir>/status.tsv` before `check_run` on every run with columns `accession, status, class (reconstructed|not_found|skipped|error), compressed_bytes, n_unitigs, n_anchors, subgraph_nodes, provisional, detail`, carrying `compressed_bytes: Option<u64>` on the result. Add per-status counts to the summary (shared with WP-03).

**Tests and docs.** Message test for `[TooLarge, TooLarge]`; `status.tsv` writer test; CLI test with `--max-size-mb 1` and accessions above it. `docs/walker.md` outputs and exit-status sections; Galaxy wrapper outputs.

### WP-15 · low · Single-anchor `no_path` drops the rejected alignment; log and summary cannot say why

`walker/src/pipeline.rs:640` · usability · confirmed (reproduced) · fix effort: small · sources: SWR6-02

**Problem.** Candidate reports are built only when `queries.len() > 1` (`walker/src/pipeline.rs:639-640`), "so single-anchor `summary.json` stays byte-identical". That reason no longer holds: 0.1.1 already added `settings`, `provisional` and `subgraph_fraction` to single-anchor summaries, and `bench/check_golden.sh` compares only FASTA sequences. A single-anchor accession whose best alignment covered 79% of the reference is printed as a bare `no_path` (`:406`), identical to one with no alignment, so a near miss cannot be told from an absent gene without rerunning with a lower `--min-cov`.

**Proposed fix.** Build `candidate_reports` unconditionally; add `ref_span: Option<f64>` (`path.query_span() / ref_len` from the GAF) to `CandidateReport` and write it; in `report()` print `no_path (best alignment spans 0.79 of the reference, under --min-cov)` or `no_path (no alignment)`. Related: WP-16.

**Tests and docs.** CLI test with a stand-in returning a 50-of-63-base row asserting the reason line and `candidates[0].ref_span` about 0.79. Comments at `pipeline.rs:635-638`, `types.rs:476-479`, `output.rs:284-285`; `docs/walker.md:59-62`.

### WP-16 · low · Summary candidates below `--min-cov` report null `cov`/`length`/`nodes`

`walker/src/pipeline.rs:655` · provenance · confirmed (reproduced) · fix effort: small · sources: RC-08, R2-08, R4-08

**Problem.** In multi-anchor mode, candidate `cov`, `length` and `nodes` come from the stitched haplotype (`cov: c.haplotype.as_ref().map(|h| h.cov)`, `walker/src/pipeline.rs:655-657`), and `stitch` discards it below the floor (`if cov < cfg.min_cov { return Ok(None); }`, `walker/src/stages/stitch.rs:219`). So exactly the references that lost on coverage show `"cov": null, "pident": 100, "length": null`, contrary to `walker/src/types.rs:376-379` ("Populated even for candidates that did not clear `min_cov`") and `docs/walker.md:61`. Selection and haplotypes are unaffected.

**Proposed fix.** Separate measurement from gate:
```rust
pub fn stitch(..) -> Result<Option<Haplotype>> {
    Ok(stitch_unfiltered(..)?.filter(|h| h.cov >= cfg.min_cov))
}
// candidate loop
let measured = path.as_ref().map(|p| stitch_unfiltered(cfg, &graph, p, &subgraph, q)).transpose()?.flatten();
let measure = measured.as_ref().map(|h| (h.cov, h.seq.len(), h.nodes));
let hap = measured.filter(|h| h.cov >= cfg.min_cov);
```
Optionally add `"passed_min_cov": bool`. A patched build reported `cov 0.588, length 290, nodes 7` for a sub-floor reference.

**Tests and docs.** `stitch_unfiltered` returns `Some` where the gated form returns `None`; CLI test with one reference under the floor. `types.rs:390-393`.

### WP-17 · low · `--max-size-mb 0` is accepted as a 0-byte cap that refuses every accession

`walker/src/main.rs:211` · usability · confirmed (reproduced) · fix effort: trivial · sources: SWR6-03

**Problem.** `let mb: u64 = value()?.parse()...` (`walker/src/main.rs:211`) has no lower bound, and `mb.saturating_mul(1024 * 1024)` makes 0 a 0-byte cap, although the comment beside it says an absurd value should mean "no cap". Every accession, cached ones included, becomes `too_large`; with `--allow-empty` the run exits 0 with an empty FASTA. The Python front ends already reject 0 (`validators.py:45`), and every other count flag uses `at_least_one`.

**Proposed fix.**
```diff
-   let mb: u64 = value()?.parse().map_err(|e| error::Error::Usage(format!("--max-size-mb: {e}")))?;
+   let mb = at_least_one("--max-size-mb", &value()?)? as u64;
```

**Tests and docs.** Add `--max-size-mb` to the flag list in `worker_and_thread_counts_must_be_at_least_one`; `docs/walker.md:197`: at least 1.

### WP-18 · low · A panic in a fetch worker aborts the whole run with exit 101

`walker/src/pipeline.rs:438` · robustness · confirmed (reproduced) · fix effort: small · sources: R6-12

**Problem.** `compute_worker` wraps each accession in `catch_unwind`, but `fetch_worker` calls `match fetcher.fetch(index, &accession) {` (`walker/src/pipeline.rs:438`) unguarded. A panic there makes `std::thread::scope` (`:306`) re-panic in `main`: exit 101 (not a documented exit code), no outputs, and the fallback at `:351-367` ("no result was produced") never runs. `run_fetch_only` (`:203-213`) has the same shape. No reachable panic was found in the shipped fetch path, so this is defence in depth; an injected panic lost all finished results.

**Proposed fix.**
```rust
let fetched = match std::panic::catch_unwind(AssertUnwindSafe(|| fetcher.fetch(index, &accession))) {
    Ok(r) => r,
    Err(p) => {
        let msg = format!("internal error: fetch panicked: {}", panic_message(&*p));
        eprintln!("[{accession}] {msg}");
        let _ = result_tx.send(AccessionResult::failed(index, accession, msg));
        continue;
    }
};
fn panic_message(p: &(dyn Any + Send)) -> &str {
    p.downcast_ref::<&str>().copied()
     .or_else(|| p.downcast_ref::<String>().map(String::as_str))
     .unwrap_or("non-string panic payload") }
```
Use `panic_message` on the compute path too (`:484-490` drops the payload), and wrap `run_fetch_only` the same way.

**Tests and docs.** Unit test with a panicking stub fetcher (needs a small `trait Fetch` seam): one `Failed` result, others complete.

### WP-19 · low · Default of 8 compute workers oversubscribes hosts with fewer cores

`walker/src/config.rs:471` · efficiency · confirmed (reproduced) · fix effort: small · sources: RF-08

**Problem.** `compute_workers: 8,` (`walker/src/config.rs:471`) applies whatever the core count. Parse and align are CPU-bound, so on a 4-core host 8 workers add no throughput but keep more payloads and graphs resident: on 24 cached 44 MiB accessions, `-w 4` peaked at 981 MiB against 1,177 MiB at the default, for the same CPU time, and the printed bound fell from 33.6 to 22.8 GB.

**Proposed fix.** Keep `Config::default()` at 8 (three tests assume it) and cap in `parse_args` only when `-w` was not given:
```rust
if !workers_set {
    if let Ok(n) = std::thread::available_parallelism() {
        args.cfg.compute_workers = args.cfg.compute_workers.min(n.get());
    }
}
```
`available_parallelism` respects cgroup CPU quotas.

**Tests and docs.** `parse_args` test: no `-w` gives `min(8, cores)`, `-w 8` gives 8. USAGE `[default: min(8, cores)]`, `docs/walker.md:199`, `config.rs:305-312`, `pipeline.rs:26`.

### WP-20 · low · Hard-wired 4-slot fetch queue adds four payloads to peak memory

`walker/src/pipeline.rs:281` · memory · confirmed (reproduced) · fix effort: trivial · sources: RF-03

**Problem.** `sync_channel::<FetchedUnitigs>(cfg.channel_capacity)` (`walker/src/pipeline.rs:281`) with `channel_capacity: 4` (`walker/src/config.rs:472`) and no flag. When compute is the bottleneck every fetch worker already waits in `send()` holding a finished payload, so the four slots only let four more downloads start. They add `4 x --max-size-mb` to the bound (20% at defaults) and the PR #8800 wrapper's `/ 6` cap formula. On 12 cached 44 MiB accessions, `-w 1 --fetch-workers 1` peaked at 356 MiB with the queue and 180 MiB with depth 0, wall time not worse; S3 throughput was not measured.

**Proposed fix.** `channel_capacity: 0` by default plus `"--queue-depth" => cfg.channel_capacity = parse_usize(..)?`. `sync_channel(0)` is a rendezvous; dropping `fetch_rx` (`:328`) still unblocks senders. Benchmark against S3 at 8/8 before changing the default.

**Tests and docs.** `peak_buffers_counts_all_three_holders` expects `(0, 8, 8)` and 16; CLI test that depth 0 with one worker each completes. `pipeline.rs` diagram, `config.rs` comments, `docs/walker.md:232-234`.

### WP-21 · low · Help text, README, comments and docs disagree with the code in about fifteen places

`walker/src/pipeline.rs:580` · docs · confirmed (reproduced) · fix effort: small · sources: R6-14, RB-10, D1-10

**Problem.** Two items affect runtime: the neighbourhood warning hard-codes `"WARNING: the 5-hop neighbourhood ..."` whatever `--hops` is (`walker/src/pipeline.rs:580`), and `LOGAN_KEEP_SCRATCH=0` keeps scratch because only presence is tested (`var_os(..).is_some()`, `walker/src/scratch.rs:39`; also `pipeline.rs:855`, `align.rs:241`, `types.rs:208`). Text errors: `walker/README.md:19` names two runtime crates (Cargo.toml has five) and both source maps omit `scratch.rs`; `-q` help says the primary is "best cov/pident" though `--primary-metric` defaults to `matched` (`main.rs:57-58`); `pipeline.rs:717` says a 180x isolate gets a 1.8x floor (it gets 5x via `--min-abund`, which `--help` does not mention); `.cargo/config.toml:4` calls the `-t 2` smoke alias "same configuration every time" though `-t 2` is non-deterministic (`config.rs:331-334`); `docs/walker.md:288` (guard is 3 unitigs per base, not 1%) and `:341` (3 attempts, not 3 retries); stale comments at `main.rs:220-222`, `main.rs:146`, `align.rs:124`, `input.rs:4,206`, `extract.rs:253`, `pipeline.rs:164`; displaced doc blocks (`main.rs:364-370`, `input.rs:11-19`); no defaults for `-o`/`-f` in USAGE.

**Proposed fix.** Pass `cfg.hops` into the warning; add `scratch::keep_scratch()` (false for unset, empty, `0`, `false`, `no`, `off`) and call it at all four sites; correct each text item as listed. Note the `LOGAN_KEEP_SCRATCH=0` change in the changelog.

**Tests and docs.** Unit test for `keep_scratch`; CLI test with `--hops 2` grepping stderr for "2-hop".

### WP-22 · low · `pipeline::run` doc says downstream stages get only the first reference

`walker/src/pipeline.rs:273` · docs · confirmed (by reading) · fix effort: trivial · sources: SWR6-04

**Problem.** "The downstream stages are handed the first [`Query`] for now" (`walker/src/pipeline.rs:273-274`) and "until Phase 4 fans them out" (`:297`) are false: `run_one` aligns against every reference (`:613`), picks the primary (`:667`) and uses the winner for prune, majority correction and alternatives. A maintainer is told coordinates are computed against `queries[0]`, the assumption behind an earlier bug recorded at `io/output.rs:449-451`.

**Proposed fix.**
```diff
-/// anchor index. The downstream stages are handed the first [`Query`] for now — later phases
-/// of the multi-anchor feature (align/stitch/pick_primary) fan out over the whole pool.
+/// anchor index. `run_one` aligns each subgraph against every reference, picks the primary
+/// with `--primary-metric`, then prunes, corrects and finds alternatives against the winner.
```
and reword `:297` accordingly.

**Tests and docs.** None required.

### WP-23 · low · Crate-wide `allow(dead_code)` hides nine dead or test-only items

`walker/src/main.rs:24` · simplification · confirmed (reproduced) · fix effort: small · sources: RS-08

**Problem.** `#![allow(dead_code, unused_variables)]` (`walker/src/main.rs:24`) disables the dead-code lint for the whole binary, so CI's `clippy -D warnings` cannot catch it. Removing it gives nine warnings: dead even in tests are `Error::NotImplemented`, `QueryIndex::len/is_empty`, `stream::recommended_chunk`, `View.local`; test-only are `align::parse_gaf`, `anchor::index_query`, `UnitigGraph::seq_len/base_at`; never read are `Anchors.shared_kmers`, `Anchors.fraction_of_graph`, `Subgraph.n_edges`. Comments at `anchor.rs:311` and `types.rs:150` describe uses that do not exist.

**Proposed fix.** Delete line 24; delete the dead items; mark test-only ones `#[cfg(test)]` (move `parse_gaf`'s tests to `parse_gaf_multi`); emit or remove the three fields (interim `#[cfg_attr(not(test), expect(dead_code, reason = "..."))]`). A patched build was clippy-clean with 221 tests passing.

### WP-24 · low · Prune survivor predicate is written three times, kept in step by a comment

`walker/src/pipeline.rs:743` · reuse · confirmed (by reading) · fix effort: trivial · sources: RS-05

**Problem.** `graph.abundance[id as usize] >= floor` appears at `walker/src/pipeline.rs:747` (survivor count, under `// Predicate must match \`prune_nodes\` exactly.`), `walker/src/stages/extract.rs:178` (the prune) and the test at `pipeline.rs:1077`. They agree today; if `prune_nodes` changes (WP-01 changes it) and the count does not, the outcome is classified on the wrong count and `unitigs_kept` stops matching the GFA.

**Proposed fix.**
```rust
pub fn survives_prune(graph: &UnitigGraph, id: UnitigId, floor: f32) -> bool {
    graph.abundance[id as usize] >= floor
}
```
Use it at all three sites, drop the comment, and optionally log `surviving` instead of `before - dropped` (`pipeline.rs:786`).

**Tests and docs.** One test that a unitig exactly at the floor survives.

### WP-25 · info · Graph-peak ratio comments give 1.5x and 1.6x where the constant is 1.7

`walker/src/pipeline.rs:32` · docs · confirmed (by reading) · fix effort: trivial · sources: BL1-11

**Problem.** `walker/src/pipeline.rs:32` says `UnitigGraph (CSR, ~1.5x compressed)`, an unexplained figure, and `:77` and `:505` say "up to ~1.6x ... (`Config::GRAPH_PEAK_PER_COMPRESSED_BYTE`)", which reads as if the constant were 1.6. It is 1.7 (`walker/src/config.rs:445`); 1.62x is the measured maximum. `docs/walker.md:237-251` already separates the two.

**Proposed fix.** `:32` "<=1.7x compressed"; `:77` and `:505` "measured up to 1.62x, counted at 1.7x (`Config::GRAPH_PEAK_PER_COMPRESSED_BYTE`)". Do not rewrite 1.62 as 1.7; the distinction matters.

### WP-26 · low · A closed stdout (SIGPIPE) skips scratch cleanup and leaves GraphAligner children running

`walker/src/main.rs:331` · reliability · confirmed (reproduced) · fix effort: small · sources: R6-09

**Problem.** `main` restores the default SIGPIPE action (`libc::signal(libc::SIGPIPE, libc::SIG_DFL);`, `walker/src/main.rs:331`) so that `| head` ends the run, but the cleanup thread listens only for `Signals::new([SIGINT, SIGTERM, SIGHUP])` (`walker/src/scratch.rs:73`). A SIGPIPE therefore kills the process without removing `logan-walker.<pid>.<suffix>`. GraphAligner, started with `cmd.output()` (`walker/src/stages/align.rs:201`), is not tracked, so it outlives the walker on SIGPIPE and on SIGTERM alike. And `create_dir_all(&cfg.scratch_dir)` runs for every accession (`walker/src/stages/extract.rs:293`), so a worker still running can recreate the directory after the signal thread has removed it. The CLI test `a_closed_stdout_ends_the_run_instead_of_hanging` (`walker/tests/cli.rs:107`) checks only the exit signal.

**Impact.** `logan-walker -l list.txt -q ref.fa | head -1` exits 141 as intended, but the GFAs and query files of the accessions in flight stay in `/dev/shm`, held in RAM until a later run over the same base sweeps them, and the aligners keep running and write GAFs into the directory. Only interactive use is affected; Galaxy and schedulers redirect stdout to a file.

**Evidence.** Offline cache, 6 accessions, `-w 4`, a stand-in aligner that sleeps 2 s, output piped into `head -1`: exit 141, and the scratch directory held 9 files, including `*.sub.gaf` written after the walker died. Stock binary: 5 of 5 runs left the directory and 4 of 5 left a live aligner. A SIGTERM run removed the directory, but the stand-in aligners kept running.

**Proposed fix.**
1. `scratch.rs`: add `SIGPIPE` to the `Signals::new` set, so this ending goes through the same remove-then-re-raise path.
2. Replace progress `println!` with a crate-level macro that, on `ErrorKind::BrokenPipe`, waits for the signal thread instead of panicking (with step 1 alone the collector panicked with "failed printing to stdout: Broken pipe"). This overlaps WP-06 step 2.
3. `align.rs`, Linux: `pre_exec` with `prctl(PR_SET_PDEATHSIG, SIGTERM)` plus a `getppid()` check, so the kernel stops the aligner whenever the walker dies, SIGKILL and the OOM killer included. For macOS, a PID registry that the signal thread sends SIGTERM to (on its own it races with spawns in flight).
4. `extract.rs:293`: stop recreating the scratch directory outside tests (test helpers create their own); make `scratch::remove` retry up to 20 times at 5 ms intervals.

A patched build left no directory and no aligner in 16 of 16 SIGPIPE runs, kept exit 141, killed all 4 aligners on SIGTERM, and passed the unit and CLI tests. `pre_exec` forces fork instead of posix_spawn (negligible here); PDEATHSIG fires when the forking thread exits, and compute threads outlive their children.

**Tests and docs.** Give `a_closed_stdout_ends_the_run_instead_of_hanging` its own `--scratch-dir` and a slow stand-in aligner, and assert the directory is empty and no child survives; add a matching SIGTERM test. `docs/walker.md:266` (add SIGPIPE, say aligners are stopped) and `:322`.

### WP-27 · low · `--min-cov`, `--min-vaf`, `--min-abund`, `--min-shared-kmers` accept NaN, inf and out-of-range values

`walker/src/main.rs:197` · input-validation · confirmed (reproduced) · fix effort: trivial · sources: R6-10, RS-13

**Problem.** `--prune-frac` and `--min-abund-frac` are range-checked (`walker/src/main.rs:283-302`); these four options are only parsed, e.g. `args.cfg.min_cov = value()?.parse()...` (`walker/src/main.rs:197`). Rust's `f64::from_str` accepts `nan`, `inf` and `-inf`. With `--min-cov nan` the gate `if cov < cfg.min_cov` (`walker/src/stages/stitch.rs:219`) is always false, so any alignment passes. The summary writes `format!("{:?}", c.min_cov)` (`walker/src/io/output.rs:223`), giving a bare `NaN` or `inf`, which is invalid JSON. `--min-shared-kmers 0` is recorded as 0 but runs as 1 (`cfg.min_shared_kmers.max(1)`, `walker/src/stages/anchor.rs:263`), and the help's "validated 1..=5" (`walker/src/main.rs:75`) is not enforced.

**Impact.** `--min-cov 85` (a percentage) turns every accession into `no_path`: exit 1 with "nothing was found", or exit 0 with an empty FASTA under `--allow-empty`. `--min-cov nan` reports DONE for any alignment. `--min-vaf inf --min-abund -inf` writes a summary that `json.load` rejects. All of these need a user error. The Python wrappers check their own `--min-cov` but do not pass these walker options through.

**Evidence.** Offline fixture, `-a SRR900000` (reconstructs at cov 1.00), stock binary: `--min-cov 85 --allow-empty` gave `no_path`, exit 0; `--min-cov nan` gave DONE and `"min_cov": NaN` in the summary; `--min-vaf inf --min-abund -inf` gave a `JSONDecodeError` at line 18; `--min-shared-kmers 0` gave DONE with `"min_shared_kmers": 0`.

**Proposed fix.** One helper for all four flags, also used by `--min-abund-frac` and `--prune-frac` so the messages match:
```rust
fn number_in(flag: &str, raw: &str, lo: f64, hi: f64) -> error::Result<f64> {
    let v: f64 = raw.parse().map_err(|_| Usage(format!("{flag} needs a number, got {raw:?}")))?;
    if !v.is_finite() || v < lo || v > hi {
        return Err(Usage(format!("{flag} must be in [{lo}, {hi}], got {raw}")));
    }
    Ok(v)
}
"--min-cov" => args.cfg.min_cov = number_in(flag, &value()?, 0.0, 1.0)?,
"--min-vaf" => args.cfg.min_vaf = number_in(flag, &value()?, 0.0, 1.0)?,
"--min-abund" => args.cfg.min_abund = number_in(flag, &value()?, 0.0, f32::MAX as f64)? as f32,
"--min-shared-kmers" => args.cfg.min_shared_kmers = at_least_one(flag, &value()?)? as u32,
```
A patched build rejected all four cases with a usage error (exit 1) and still parsed valid values. Command lines with `--min-shared-kmers 0` now fail fast; `--min-cov 0` is still accepted. Related: WP-17.

**Tests and docs.** Unit test `thresholds_are_range_checked` in `main.rs`. Help at `main.rs:75`: "1 or more; tested 1..=5", or enforce 1..=5, and `docs/walker.md:195` to match; state the ranges of `--min-cov` and `--min-vaf` at `docs/walker.md:202,212`.

### WP-28 · low · The argument parser keeps only the last repeated `-a`/`-l` and takes the next option as a value

`walker/src/main.rs:179` · usability · confirmed (reproduced) · fix effort: small · sources: R6-11

**Problem.** `parse_args` overwrites single-value options when they are repeated: `"-a" | "--accession" => args.accession = Some(value()?)` (`walker/src/main.rs:179`), and likewise `-l`, `-o` and `-f` (`:180-183`). `value()` returns `argv.get(i + 1).cloned()` (`:174-175`), whatever that token is, another option included. `-q` is documented as repeatable, so repeating `-a` is a natural mistake; the code guards against `-a` together with `-l` (`:472-473`, "`-a` used to win silently") but not against `-a` twice. `std::env::args()` (`:338`) panics on a non-UTF-8 argument.

**Impact.** `-a SRR900000 -a SRR900001` processes only SRR900001 and exits 0; `-l l1.txt -l l2.txt` processes only `l2.txt` and exits 0; `-o --emit-vcf` writes the outputs into a directory named `--emit-vcf`, drops the VCF request and exits 0. `-a --offline` takes `--offline` as the accession and fails with `no_logan_data`, exit 1. A `-q` path with a non-UTF-8 byte panics with exit 101. The worst outcome is a silently partial run that exits 0, and it needs a malformed command line.

**Evidence.** Offline fixture, stock binary: each case above reproduced as stated (`[1/1] SRR900001: DONE ... 1/1 accessions reconstructed`; `[1/1] SRR900002: DONE` for lists of 2 and 1 accessions; `panicked at std/src/env.rs` for `-q $'q\xff.fa'`).

**Proposed fix.**
```rust
let value = || match argv.get(i + 1) {
    None => Err(Usage(format!("{flag} needs a value"))),
    Some(v) if v.starts_with('-') && v.parse::<f64>().is_err() =>
        Err(Usage(format!("{flag} needs a value, got the option {v}"))),
    Some(v) => Ok(v.clone()),
};
// -a, -l, -o, -f, --scratch-dir, --unitigs-cache, --graphaligner: usage error when given twice,
// with the hint "put several accessions in one --accession-list file".
let argv: Vec<String> = match std::env::args_os().skip(1).map(|a| a.into_string()).collect() {
    Ok(v) => v,
    Err(bad) => { eprintln!("error: argument {bad:?} is not valid UTF-8"); return ExitCode::FAILURE }
};
```
Numeric values such as `-1` still parse, so WP-27's range checks can name them. Optionally accept `--opt=value` for a value that really starts with `-`. A patched build failed every case above with a usage error and exit 1. A path that begins with `-` (`-o -out`) is now rejected; `./-out` works.

**Tests and docs.** Unit tests `repeated_single_value_options_are_rejected` and `an_option_is_not_taken_as_a_value` in `main.rs`. USAGE: mark `-a` as "once; use -l for several".
