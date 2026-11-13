# Walker: GraphAligner invocation and stitching (align.rs, stitch.rs)

The GraphAligner wrapper and the stitcher are carefully written against defects the project has already met: the query is never passed under the user's file name, a gzipped query is decompressed, GAF rows are kept per query name, and scratch files are keyed on PID and task so no worker can read another's GAF. What remains fragile falls into two groups. First, the single-record "byte-compat" path in `align::run` re-opens and copies the user's query file instead of using the sequence `read_query` already parsed, so the aligner and the anchoring stage see different inputs (WA-01, WA-07, WA-08, WA-11). Second, the walker trusts GraphAligner's GAF fields, and the primary's labels, more than they deserve: ambiguous reference bases count as matches, `--min-cov` measures walked length rather than reference span, the path is not checked against the graph, and a mosaic of two strains is labelled `vaf=1.00` (WA-02 to WA-05, WA-12). The common theme of the fixes is to use the data the walker already holds (the parsed query, the graph, the reference span, route depth) instead of re-reading or trusting external bytes, and to report in the output what the gates actually measured. Most defects need a pipe, an unusual query file, a mixed sample, indels near the coverage threshold or a misbehaving aligner. GraphAligner was not installed on the review machine: every reproduction used the shipped release binary with a small stand-in aligner script, and GraphAligner's own behaviour was taken from its v1.0.17b source.

| ID | Severity | Where | Title | Fix effort |
|----|----------|-------|-------|------------|
| WA-01 | medium | `walker/src/stages/align.rs:169` | A piped `-q` is re-read for every alignment, so every accession becomes `no_path` | small |
| WA-02 | medium | `walker/src/stages/stitch.rs:139` | N and IUPAC reference bases count as matches: N-laden panel records win and report inflated pident | small |
| WA-03 | medium | `walker/src/stages/stitch.rs:165` | Multi-anchor primary choice ignores route depth; on a tie, panel order can pick a 5% strain | small |
| WA-04 | medium | `walker/src/stages/stitch.rs:64` | "A connected path makes chimeras impossible" is false for mixtures; `_primary` has no mixture flag | small |
| WA-05 | medium | `walker/src/stages/stitch.rs:198` | `--min-cov` gates on walked length, not reference span: indels flip accept and reject | small |
| WA-06 | medium | `walker/src/stages/align.rs:201` | No deadline on GraphAligner: one stuck alignment stalls the run and loses all finished results | small |
| WA-07 | medium | `walker/src/stages/align.rs:169` | Single-record query header reaches GraphAligner verbatim: a tab gives `no_path`, a non-UTF-8 byte fails every accession | trivial |
| WA-08 | medium | `walker/src/stages/align.rs:169` | Single-record query: whitespace inside sequence lines reaches GraphAligner, giving `no_path` or SIGABRT | small |
| WA-09 | low | `walker/src/pipeline.rs:800` | Prune pass accepts a realignment that only trims 1-3 end bases | small |
| WA-10 | low | `walker/src/stages/align.rs:231` | Exit 0 without a GAF is reported as the biological negative `no_path` | trivial |
| WA-11 | low | `walker/src/io/input.rs:180` | Query bytes outside GraphAligner's alphabet pass validation, then abort every alignment | trivial |
| WA-12 | low | `walker/src/stages/stitch.rs:92` | GAF paths are not checked against the graph; an inconsistent row is published as a primary | small |
| WA-13 | low | `walker/src/stages/align.rs:211` | `LOGAN_KEEP_SCRATCH` ignored on failure paths, `0`/empty enable it, stderr tail reversed | small |
| WA-14 | low | `walker/src/stages/align.rs:190` | Multi-anchor mode re-logs the full reference-name list for every accession | trivial |

### WA-01 · medium · A piped `-q` is re-read for every alignment, so every accession becomes `no_path`

`walker/src/stages/align.rs:169` · correctness · confirmed (reproduced) · fix effort: small · sources: RE-03, D1-01, R4-02

**Problem.** `read_query` parses `-q` once at startup, but when the query file holds one record, `align::run` opens the user's path again on every call: `let bytes = crate::io::input::read_input(&queries[0].path, "query")?;` (`walker/src/stages/align.rs:169`). This runs once per accession and once more per prune realignment (`walker/src/pipeline.rs:792`). A pipe (`<(zcat ref.fa.gz)`, `/dev/stdin`) is empty after the first read, so GraphAligner receives a 0-byte query FASTA and finds nothing. A multi-record query through a pipe fails differently: `walker/src/main.rs:460-462` reads the path a second time in `read_query_set` and stops with "no FASTA records found". `README.md:96-99` explicitly allows "a pipe such as `<(zcat hits.tsv.gz)`" for every input file.

**Impact.** Anchoring uses the in-memory sequence and succeeds, so the run looks normal; every accession is then `no_path`, the run ends "nothing was found" (exit 1), and with `--allow-empty` it exits 0 with empty output that reads as a biological negative. Editing or deleting the query file during a long run also changes the reference for the remaining accessions. Galaxy passes regular files and is not affected; shell users following the README are.

**Evidence.** Release binary, offline cache with 3 accessions, a stand-in aligner logging the size of its `-f` file: `-q q.fasta` reconstructs 3/3 (query FASTA 70 bytes each time); `-q <(cat q.fasta)` and `cat q.fasta | logan-walker -q /dev/stdin …` give `no_path` for all 3, query FASTA 0 bytes, exit 1; `-q <(cat two_records.fa)` gives `usage error: /dev/fd/63: no FASTA records found`.

**Proposed fix.** Never reopen the query path after startup.
1. In `align::run`, delete the byte-copy branch and always write the normalised records the walker already holds, `>{record_id}\n{seq}\n`, which is the existing merged-pool branch. This also fixes WA-07 and WA-08, which have the same line as their cause.
```rust
// align.rs, replacing the `if queries.len() == 1 && queries[0].n_records == 1 { ... } else { ... }`
let mut body = Vec::with_capacity(queries.iter().map(|q| q.seq.len() + q.record_id.len() + 3).sum());
for q in queries {
    body.push(b'>'); body.extend_from_slice(q.record_id.as_bytes()); body.push(b'\n');
    body.extend_from_slice(&q.seq); body.push(b'\n');
}
std::fs::write(&query_fa, &body).map_err(file_err("write", &query_fa))?;
```
2. In `main.rs:460-462`, keep the decompressed bytes from the first read and parse the record set from them: split `read_query_set` into a reader and `read_query_set_from(Vec<(PathBuf, Vec<u8>)>)`, and call the latter when `first.n_records > 1`.

If the maintainers want the GraphAligner input to stay byte-identical for the golden fixtures, the narrower alternative is to store the decompressed bytes on `Query` (`raw: Arc<[u8]>`) at startup and write those; that fixes the re-read but leaves WA-07 and WA-08 open. Either way, run `bench/check_golden.sh` once under GraphAligner 1.0.17b before release.

**Tests and docs.** Add a `walker/tests/cli.rs` test that feeds the fixture query on stdin (`-q /dev/stdin`) and asserts exit 0 and a non-empty haplotype FASTA, plus the same with a 2-record query. Update `graphaligner_gets_a_fasta_named_byte_copy_of_the_query` (`walker/src/stages/align.rs:649`) and the "byte-for-byte copy" doc comments at `walker/src/stages/align.rs:133-136` and `walker/src/types.rs:39-44`. `README.md:96-99` then becomes true as written.

### WA-02 · medium · N and IUPAC reference bases count as matches: N-laden panel records win and report inflated pident

`walker/src/stages/stitch.rs:139` · science · confirmed (reproduced) · fix effort: small · sources: SWR4-02

**Problem.** The walker uses GAF column 10 (`num_matches`) as given: `pick_primary` scores `p.num_matches as f64 / c.ref_len as f64` (`walker/src/stages/stitch.rs:139`) and `stitch` reports `pident = num_matches / aln_len * 100` (`walker/src/stages/stitch.rs:214`). GraphAligner 1.0.17b counts that column with `characterMatch`, whose `ambiguousMatch` returns true for `N` against any base and for R/Y/K/M/S/W/B/D/H/V against their member bases (`GraphAlignerCommon.h:240-283` in its source). Every N in a reference is therefore a guaranteed match. `tests/test_panels.py:12-14` states that Ns "count against identity", and `docs/walker.md:71` says a reference loses a point for every base the walk mismatches; both are false for N.

**Impact.** In multi-anchor mode a reference with a run of Ns beats the true nearest reference wherever the sample differs inside the N window, and then drives pruning, majority correction, alternatives and VCF coordinates. In any mode, pident is inflated. The haplotype sequence itself is unaffected. `queries/master_planetary_targets.fasta` contains records with 659 N (zika_ns5_var3) and 862 N/Y (ebola_gp_var2); the shipped `panels/sarscov2_spike_11strains.fa` is currently clean, but `tests/test_panels.py` records that two of its records carried 268 and 477 Ns in 0.1.0.

**Evidence.** Panel [refA, refN], refN = refA with 120 N at positions 90-209, sample differing from refA at 6 sites in that window, stand-in aligner counting matches as GraphAligner does: refA 294 matches (pident 98), refN 300 (pident 100), output `>SRRFAKE0_primary cov=1.000 pident=100.0 … anchor=refN`.

**Proposed fix.** Correct the count once, right after `parse_gaf_multi` in `align::run`, so the first pass, the prune pass, `matched`, `pident` and the candidate report all see it:
```rust
let mut map = parse_gaf_multi(&bytes, &expects)?;
for q in queries {
    if let Some(p) = map.get_mut(q.record_id.as_str()) {
        let amb = q.seq[p.query_start..p.query_end].iter()
            .filter(|b| !matches!(b, b'A' | b'C' | b'G' | b'T')).count();
        p.num_matches = p.num_matches.saturating_sub(amb);
    }
}
```
The correction is exact for N aligned as a match and conservative for N placed in an insertion or IUPAC codes that mismatched. With it, the run above gives `anchor=refA pident=98.0` and refN 180 matches. Optionally warn about, or reject, panel records with more than a set fraction of non-ACGT bases.

**Tests and docs.** Unit test for the discount and a `pick_primary` test with the refA/refN numbers above. Update the `GafPath::num_matches` doc (`walker/src/types.rs:244`) and `docs/walker.md:71` to say ambiguous reference bases never count as matched; extend `tests/test_panels.py` to `queries/*.fasta`. Golden pident values change only if a fixture query contains IUPAC codes.

### WA-03 · medium · Multi-anchor primary choice ignores route depth; on a tie, panel order can pick a 5% strain

`walker/src/stages/stitch.rs:165` · science · confirmed (reproduced) · fix effort: small · sources: SC1-07

**Problem.** `pick_primary` ranks candidates only by how closely each walk resembles its own reference (`PrimaryMetric::Matched => cmpf(a.4, b.4).then_with(|| cmpf(ap, bp))`, `walker/src/stages/stitch.rs:156`) and breaks ties by pool order (`primary.then(a.0.cmp(&b.0))`, `walker/src/stages/stitch.rs:165`). The route's depth, `Haplotype::mean_ka` (`walker/src/stages/stitch.rs:233`), never enters the ranking. When a sample holds two strains that each equal a panel record, both candidates score 1.0 matched and pident 100, and the record listed first wins. Nothing repairs it later: `prefer_majority` skips a fork whose taken branch is below `min_vaf` (`walker/src/stages/variants.rs:188`), so a 5% route is never moved to the 95% one, and the candidate report (`CandidateReport`, `walker/src/types.rs:381`) carries no depth, so a downstream tool cannot re-rank as `docs/walker.md:61-62` suggests.

**Impact.** Multi-anchor runs on mixed samples (co-infection, wastewater): the primary, its `anchor=` tag and every VCF/SNV coordinate follow the minority strain, the majority strain leaves no trace in alternatives or the VCF, and the result depends silently on panel order.

**Evidence.** Synthetic graph with strain A at 9,500x and strain B at 500x (5 SNVs apart), stand-in aligner returning the maximum-identity path: panel order (refA, refB) gives `anchor=refA ka=9772.7`; order (refB, refA) gives `anchor=refB ka=5681.8`; both candidates report 900 matches and pident 100; 0 alternatives, 0 VCF rows.

**Proposed fix.** Rank by the metric, treat candidates within a small band of the best as tied, and choose among them by `mean_ka`; keep the current order when depths are equal:
```rust
let best = key(&with_haps[0]);                              // key = the metric's primary value
let mut band: Vec<_> = with_haps.into_iter().filter(|c| best - key(c) <= cfg.primary_tie_eps).collect();
band.sort_by(|a, b| cmpf(a.3.mean_ka, b.3.mean_ka));        // cmpf sorts descending; stable, so equal depth keeps metric order
```
Add `--primary-tie-eps` (default 0.002; `0` restores 0.1.1 behaviour) and add `mean_ka` and the route's minimum ka to `CandidateReport` and the summary JSON. With this, both panel orders select refA. Measure the band width on the 11-strain panel before release and document it as a 0.2 behaviour change.

**Tests and docs.** Unit test `matched_tie_goes_to_the_deeper_route_not_panel_order` (900/900 matches, ka 500 listed first versus 9,500: the deeper wins; a 98%-identity deep candidate still loses to an exact one). The existing `primary_metric_still_breaks_ties_on_pool_order` still holds. Update `pick_primary`'s doc comment (`walker/src/stages/stitch.rs:114-129`) and the "Multi-anchor mode" section of `docs/walker.md`. Related: WA-04.

### WA-04 · medium · "A connected path makes chimeras impossible" is false for mixtures; `_primary` has no mixture flag

`walker/src/stages/stitch.rs:64` · science · confirmed (reproduced) · fix effort: small · sources: SC1-09

**Problem.** `walker/src/stages/stitch.rs:64` says "Requiring a connected path makes them [chimeras] impossible" (also `:14`). That holds only for joins no k-mer supports. In a de Bruijn graph with k=31, allele sites more than 30 bases apart are joined through shared unitigs, so every combination of alleles is a connected path that no read links; `walker/src/stages/variants.rs:32-34` already says this for alternatives. The primary is such a combination, and the header hard-codes `vaf: 1.0` (`walker/src/stages/stitch.rs:234`), written by `walker/src/io/output.rs:84`, so a mosaic of two strains looks clonal in the FASTA that feeds `disassembler msa` and the tree.

**Impact.** Accessions holding two strains that differ at two or more sites more than 30 nt apart, with depths that do not separate cleanly, produce a `_primary` that may match neither strain; downstream it reads as excess substitutions and recombination. Alternatives and the VCF still carry the variation, but the main FASTA gives no warning.

**Evidence.** Synthetic graph, strain A at 60% over the 5' half and 40% over the 3' half, strain B the reverse, 6 sites: header `>SRRMIX_primary cov=1.000 pident=100.0 nodes=13 ka=723.1 vaf=1.00`; the primary differs from each strain at 3 of 900 sites (alleles BBBAAA, carried by no molecule); same with `--no-prefer-majority`.

**Proposed fix.**
1. Reword `stitch.rs:14` and `:59-64`, `docs/walker.md:3` and `README.md:7`: a connected path rules out joins without k-mer support; it does not phase sites more than k-1 apart; `_primary` is a majority-consensus path, not a phased haplotype.
2. Count on-path forks where both branches clear `min_abund`/`min_vaf` (`variants::view` has the numbers) and record `two_allele_sites` and `min_site_vaf` in `Haplotype` and the summary JSON. Add `min_vaf=` and `mixed=1` (when `two_allele_sites >= 2`) at the end of the primary header, keeping `vaf=` for one release.
3. Optionally give `disassembler msa`/`sieve` an `--exclude-mixed` option that reads the tag.

**Tests and docs.** Unit test in `variants.rs`: two 60/40 forks 200 nt apart give `two_allele_sites == 2` and `mixed=1`; a clonal diamond gives 0. CLI test parsing the header; update the Outputs table in `docs/walker.md` and the `bench/` golden headers. Related: WA-03.

### WA-05 · medium · `--min-cov` gates on walked length, not reference span: indels flip accept and reject

`walker/src/stages/stitch.rs:198` · correctness · confirmed (reproduced) · fix effort: small · sources: RA-02, R4-05

**Problem.** `stitch` computes `let cov = target.len() as f64 / query.seq.len() as f64;` (`walker/src/stages/stitch.rs:198`), haplotype length over reference length, and rejects with `if cov < cfg.min_cov { return Ok(None); }` (`:219`). `--min-cov` is documented as the "minimum fraction of the reference a primary must cover" (`docs/walker.md:202`, `walker/src/main.rs:71`, `walker/src/types.rs:267`). The reference span (GAF columns 3-4, `GafPath::query_span()`) is parsed and used by the prune rule but not by this gate, so deletions count against coverage and insertions for it. The published `cov` value itself is documented as walked length (`docs/walker.md:112`); the gate is what disagrees.

**Impact.** With the default 0.85 on a 300-nt gene: a sample carrying the whole gene with a 60-nt in-frame deletion gets cov 0.80 and is reported `no_path` ("gene absent"), biasing presence calls against deletion alleles; an alignment that reaches only reference base 240 but walks a 20-nt insertion gets cov 0.867 and is published.

**Evidence.** Stand-in aligner returning canned GAF rows: `SRRDEL1: no_path` (full span, deletion); `>SRRINS1_primary cov=0.867` with reference span 0-240 of 300. With the gate on reference span, `SRRDEL1` is reconstructed and `SRRINS1` is `no_path`; all existing tests still pass, so none pins the gate.

**Proposed fix.**
```rust
-    if cov < cfg.min_cov {
+    let ref_cov = path.query_span() as f64 / query.seq.len() as f64;
+    if ref_cov < cfg.min_cov {
```
Keep `cov` (walked/reference) in the output for compatibility, and add `ref_span=<qs>-<qe>` to the primary header (`walker/src/io/output.rs`) and `ref_cov`/`ref_span` to each record and candidate in summary.json. Replace the `COV_EXCESS` heuristic (`:204`) with an exact check that GAF column 2 (query length) equals `query.seq.len()`.

**Tests and docs.** Unit tests: span 0-300 over a 240-nt target is accepted at 0.85; span 0-240 over a 280-nt target is rejected. CLI test with the deletion/insertion fixture. Say "reference span" at `docs/walker.md:202`/`:299` and `walker/src/types.rs:267`; add `ref_span=` to the header table. Some reruns change status; regenerate golden fixtures if needed.

### WA-06 · medium · No deadline on GraphAligner: one stuck alignment stalls the run and loses all finished results

`walker/src/stages/align.rs:201` · reliability · confirmed (reproduced) · fix effort: small · sources: R4-03, RE-09

**Problem.** `align::run` waits with `let output = cmd.output().map_err(...)` (`walker/src/stages/align.rs:201`) and no deadline; the only timeouts in `walker/src/` are the HTTP ones in `stages/stream.rs`. The subgraph-size guard only marks a result "provisional". All outputs are written after `pipeline::run` returns (`walker/src/main.rs:554-581`), and the signal handler (`walker/src/scratch.rs:75-83`) only removes scratch and exits.

**Impact.** One GraphAligner call that never returns (a pathological neighbourhood, a stalled filesystem, a wrapper waiting on a lock) holds its worker forever; a scheduler wall-time kill then loses every accession already finished in the shard, with no FASTA, summary or VCF. Whether real GraphAligner 1.0.17b can run unbounded on a Logan subgraph was not observed.

**Evidence.** Six cached accessions, `--offline -w 2`, a stand-in aligner that sleeps 300 s for one accession: five accessions finish in under a second, no output directory exists, SIGTERM exits with nothing written, and the `sleep` child survives with parent PID 1.

**Proposed fix.** Replace `cmd.output()` with `spawn()`, drain stderr on a helper thread (so a chatty child cannot block on a full pipe), and poll `try_wait()` every 50 ms against a new `Config::graphaligner_timeout_secs` (`--graphaligner-timeout <SECS>`, default 3600, `0` = none). On expiry:
```rust
child.kill().ok(); child.wait().ok();
let _ = std::fs::remove_file(gaf_out); let _ = std::fs::remove_file(&query_fa);
return Err(Error::GraphAligner { accession: subgraph.accession.clone(),
    detail: format!("timed out after {secs} s (--graphaligner-timeout)") });
```
Only that accession becomes `error`; with `--graphaligner-timeout 2` the run above exits 0 with 5/6 reconstructed and all outputs written. As a separate, larger follow-up, write per-accession results incrementally or flush on signal so an interrupted shard keeps finished work.

**Tests and docs.** Unit test with a fake aligner sleeping 10 s and a 1 s timeout: `Err` containing "timed out" within 3 s and no leftover GAF. CLI test using the existing `ga_slow` stand-in with `--graphaligner-timeout 1` expecting outputs present. Add the flag to the options table in `docs/walker.md` and the Galaxy wrapper, and list a timeout as `error` in the outcomes table.

### WA-07 · medium · Single-record query header reaches GraphAligner verbatim: a tab gives `no_path`, a non-UTF-8 byte fails every accession

`walker/src/stages/align.rs:169` · correctness · confirmed (reproduced) · fix effort: trivial · sources: R4-01, RE-04

**Problem.** For a one-record `-q` file, `align::run` writes the user's file byte for byte (`walker/src/stages/align.rs:169-170`), header included. GraphAligner 1.0.17b takes the read name as the whole header line and writes it unescaped as GAF column 1 (`fastqloader.h:58`, `GraphAlignerGAFAlignment.h:193` in its source). `parse_gaf_multi` then requires UTF-8 (`std::str::from_utf8(gaf_bytes)`, `walker/src/stages/align.rs:268`) and splits on tabs (`:278`). A tab in the header shifts every column, the path column becomes `+`, and the row is dropped. `read_query` accepted the same header (lossy UTF-8, first word). GraphAligner 1.0.19 and later truncate the name at the first blank, so this is specific to the pinned version.

**Impact.** `>gene<TAB>description`: every accession `no_path`, exit 1 (or empty exit 0 with `--allow-empty`). A Latin-1 header such as `>gene isolate S\xe3o`: every accession fails with "GAF is not valid UTF-8", exit 3. In multi-anchor mode with one record per file (`-q a.fa -q b.fa`), pass 1 works (the merged FASTA carries only `>record_id`), but the prune realignment of the winner (`walker/src/pipeline.rs:792`) goes through the byte-copy path and is always refused.

**Evidence.** Stand-in aligner reproducing 1.0.17b's naming: space in header gives DONE, exit 0; tab gives `no_path`, exit 1; Latin-1 gives `GAF is not valid UTF-8`, exit 3; multi-anchor with a tab gives `prune: 1 refused (realignment gave no haplotype)`.

**Proposed fix.** Write `>{record_id}` instead of the user's header. The fix proposed in WA-01 (always write the normalised record from memory) does this and also covers WA-08. If the maintainers prefer to keep the sequence bytes unchanged, replace only the header line:
```rust
let hdr = bytes.iter().position(|&b| b == b'>').unwrap_or(0);
let body = bytes[hdr..].iter().position(|&b| b == b'\n').map_or(bytes.len(), |i| hdr + i + 1);
let mut fa = Vec::with_capacity(bytes.len() + queries[0].record_id.len() + 2);
fa.push(b'>'); fa.extend_from_slice(queries[0].record_id.as_bytes()); fa.push(b'\n');
fa.extend_from_slice(&bytes[body..]);
```
Parsing the GAF as bytes is then optional hardening.

**Tests and docs.** Update `graphaligner_gets_a_fasta_named_byte_copy_of_the_query` (`walker/src/stages/align.rs:649`); add a test with a tab and a Latin-1 header using a 1.0.17b-like stand-in. Fix the "byte-for-byte copy" wording at `walker/src/stages/align.rs:133-136`, `walker/src/types.rs:40-41` and `docs/walker.md:90-95`. Confirm with `bench/check_golden.sh` under 1.0.17b.

### WA-08 · medium · Single-record query: whitespace inside sequence lines reaches GraphAligner, giving `no_path` or SIGABRT

`walker/src/stages/align.rs:169` · correctness · confirmed (reproduced) · fix effort: small · sources: SWR4-01

**Problem.** `read_query` removes all whitespace inside sequence lines (`walker/src/io/input.rs:143`), and anchoring uses that sequence. The single-record branch of `align::run` copies the raw file instead (`walker/src/stages/align.rs:169-170`). GraphAligner 1.0.17b strips only a trailing `\r` and trailing spaces from each line (`fastqloader.h:66-68`), so internal spaces and tabs stay in the read. Its minimizer seeding restarts at every non-ACGT byte, and in the DP a blank reaches `ambiguousMatch`'s `default:` branch, which calls `std::abort()` (`GraphAlignerCommon.h:284-286`).

**Impact.** Codon-spaced sequence (`ATG GCT AAA …`) gives no 19-base seed and every accession is `no_path`, with nothing pointing at the query. GenBank-style 10-base blocks or a stray tab give SIGABRT on every accession, exit 3. The same file works inside a multi-anchor pool, but the prune pass realigns the winner alone through the byte path and fails or aborts.

**Evidence.** Stand-in aligner following these three 1.0.17b rules: plain file DONE; codon-spaced `no_path`, exit 1; 10-base blocks and trailing tab `GraphAligner failed: exited with signal: 6 (SIGABRT)`, exit 3; the same files in a pool DONE. With the WA-01 fix all finish DONE and the prune pass succeeds.

**Proposed fix.** The fix in WA-01: drop the byte-copy branch and always write `>{record_id}\n{seq}\n` from the parsed query. Any variant that keeps the sequence lines byte for byte (the header-only fix in WA-07, or storing raw bytes) leaves this defect.

**Tests and docs.** Rewrite `graphaligner_gets_a_fasta_named_byte_copy_of_the_query` (`walker/src/stages/align.rs:649`) to assert `>pol\nACGTACGTAAAA\n` for both plain and gzipped sources; add a CLI test with a codon-spaced query and a stand-in that fails on internal blanks; remove the "byte-compat path" comments at `walker/src/stages/align.rs:133-136` and `:163-168`. Run `bench/check_golden.sh` once under 1.0.17b; uppercasing and dropping the header description are expected to be neutral.

### WA-09 · low · Prune pass accepts a realignment that only trims 1-3 end bases

`walker/src/pipeline.rs:800` · correctness · confirmed (reproduced) · fix effort: small · sources: RA-05

**Problem.** The pruned realignment replaces the first pass when `pruned.query_span() + PRUNE_SPAN_TOLERANCE >= unpruned.query_span()` (`walker/src/pipeline.rs:159`, tolerance 3 at `:174`). The guard at `:800` does not check whether pruning changed anything else, so a realignment that is the first-pass haplotype minus up to 3 end bases is accepted and logged as "prune: 1 applied": pure loss, for example a truncated stop codon with exit 0. The maintainers' own count in `pipeline.rs:165-167` (33 of 511 refusals in a 20,000-accession run lost 3 bases or fewer) suggests the case is real under 1.0.17b.

**Evidence.** Stand-in aligner replaying a first pass `[0,1275)` and a second pass `[0,1274)` on the same route: `prune: 1 applied`, primary length 1274 ending `CTA`; with `--prune-frac 0`, length 1275 ending `TAA`.

**Proposed fix.** Refuse a shorter pruned result that is only a trimmed copy of the first-pass haplotype:
```rust
fn trimmed_copy(pruned: &[u8], unpruned: &[u8]) -> bool {
    pruned.len() < unpruned.len() && unpruned.windows(pruned.len().max(1)).any(|w| w == pruned)
}
// pipeline.rs:800
Some((p, h)) if pruned_covers_as_much(path, &p)
    && !trimmed_copy(&h.seq, haplotype.as_ref().map(|u| u.seq.as_slice()).unwrap_or(&[])) => {
```
Give this case its own `PruneResult` and log message (the current refusal text says "past the 3-base tolerance"). The related one-base end clipping with GraphAligner 1.0.19 and later is a recorded design decision (`docs/issue_graphaligner_end_clipping.md`) and is only warned about at `walker/src/main.rs:532-537`; failing unless an override flag is given is optional.

**Tests and docs.** Unit test: a pruned haplotype equal to `unpruned[..n-1]` is refused; one with an interior change is accepted. Say "only if it changes the haplotype beyond trimming" at `docs/walker.md:204`.

### WA-10 · low · Exit 0 without a GAF is reported as the biological negative `no_path`

`walker/src/stages/align.rs:231` · correctness · confirmed (reproduced) · fix effort: trivial · sources: R4-04

**Problem.** A missing GAF after a clean exit maps to `Ok(HashMap::new())`, i.e. `no_path` (`walker/src/stages/align.rs:231-234`). GraphAligner 1.0.17b always creates the GAF before aligning and aborts if it cannot (`Aligner.cpp:198-203, 868` in its source), so this state means the configured program is not GraphAligner or a wrapper swallowed a failure. The startup version probe accepts any program that prints a line and exits 0. `logan-walker … --graphaligner /usr/bin/true --allow-empty` reports `no_path` and exits 0.

**Proposed fix.**
```rust
Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
    if !crate::scratch::keep() { let _ = std::fs::remove_file(&query_fa); }
    return Err(Error::GraphAligner { accession: subgraph.accession.clone(),
        detail: format!("{} exited 0 but wrote no GAF", bin.display()) });
}
```
Do not require "GraphAligner" in the `--version` line: source builds print `Version Branch …`. (`scratch::keep()` is introduced in WA-13.)

**Tests and docs.** Change the `ga_fast` stand-in in `walker/tests/cli.rs:37` from `"exit 0"` to `": > \"$6\""` (create the `-a` file). Add a CLI test where a stand-in exits 0 without writing the GAF (expect error, exit 3, even with `--allow-empty`) and keep one where an empty GAF gives `no_path`. In the `docs/walker.md` status table, define `no_path` as "GraphAligner ran and found no alignment".

### WA-11 · low · Query bytes outside GraphAligner's alphabet pass validation, then abort every alignment

`walker/src/io/input.rs:180` · usability · confirmed (reproduced) · fix effort: trivial · sources: SWR4-03

**Problem.** `read_query` and `read_query_set` accept any non-whitespace byte; non-ACGT bytes only produce a warning about k-mer anchoring (`walker/src/io/input.rs:180`, `:281`). GraphAligner 1.0.17b aborts (`GraphAlignerCommon.h:284-286`) on any read byte outside ACGTU, IUPAC codes and `-`. A trailing stop `*`, a `.`, `X` or GenBank position digits therefore make every accession fail with SIGABRT, exit 3, after all downloads, with no hint at the query; in a pool, one bad record fails the clean records too.

**Proposed fix.** Add `check_alphabet(path, record_id, seq)` in `walker/src/io/input.rs`, called in both readers after the empty-sequence check, that rejects the first byte not in `b"ACGTUNRYKMSWBDHV-"` as `Error::Usage`, naming the record, byte and 1-based position (for example `record "gene" has '*' at position 301`). Related: WA-08.

**Tests and docs.** Unit tests: `ACGT*` and `AC.GT` rejected, `ACGTN` and `ACGTR` accepted. Document the accepted alphabet in `docs/walker.md` and the `-q` help.

### WA-12 · low · GAF paths are not checked against the graph; an inconsistent row is published as a primary

`walker/src/stages/stitch.rs:92` · robustness · confirmed (reproduced) · fix effort: small · sources: R4-06

**Problem.** `stitch_steps` appends `&s[k_overlap..]` for each step (`walker/src/stages/stitch.rs:92`) without checking that consecutive steps are a graph edge; `parse_path` skips unknown bytes and ids above `u32::MAX` (`walker/src/stages/align.rs:402-414`); `stitch` clamps coordinates (`stitch.rs:189`) and never compares GAF `path_len` with the stitched length. A GAF that disagrees with the GFA (a different aligner build, a wrapper, a truncated file) becomes a chimera labelled `cov=1.000 pident=100.0`. With a stand-in emitting `>0>1>2` where the true walk is `>0<1>2`, the walker published a primary with 89 differences and pident 100.

**Proposed fix.** Add `stitch::check_gaf_walk(graph, subgraph, path) -> Result<(), String>` that requires each consecutive pair to be an edge (in either orientation), the stitched length to equal `path.path_len`, and `path_start <= path_end <= len`. Call it before `stitch` at the two GAF ingestion points in `walker/src/pipeline.rs` (`:623`, `:794`), not on spliced paths, and map a failure to `Error::GraphAligner`. Make `parse_path` return an error on any byte other than `>`/`<`/digits and on oversize ids. Do not validate inside `stitch_steps`: toy test graphs there are not true k-1 walks.

**Tests and docs.** Add unit tests for a non-edge step, a `path_len` mismatch and an out-of-range end, and a CLI test where the stand-in emits `>0>1>2` and expects exit 3. List "GAF inconsistent with the graph" under the error outcomes in `docs/walker.md`, and run the real-data golden gate once.

### WA-13 · low · `LOGAN_KEEP_SCRATCH` ignored on failure paths, `0`/empty enable it, stderr tail reversed

`walker/src/stages/align.rs:211` · usability · confirmed (reproduced) · fix effort: small · sources: R4-09

**Problem.** The debug switch is honoured only on success (`walker/src/stages/align.rs:241`); when GraphAligner fails, `:211-212` (and `:232`, `:236-237`) delete the GAF and query FASTA unconditionally, so the inputs needed to rerun the failing alignment are gone. All sites test presence only, so `LOGAN_KEEP_SCRATCH=0` or an empty value keeps every file in `/dev/shm` for the whole run. `:214` prints the last three stderr lines newest first ("Aborted | Assertion failed: …").

**Proposed fix.** One helper in `walker/src/scratch.rs`, used at every site, wrapping the failure-path `remove_file` calls in `if !crate::scratch::keep()`:
```rust
pub fn keep() -> bool {
    matches!(std::env::var_os("LOGAN_KEEP_SCRATCH"), Some(v) if !v.is_empty() && v != "0")
}
// align.rs:214
let mut tail: Vec<&str> = tail.lines().rev().take(3).collect();
tail.reverse();
```

**Tests and docs.** Test `keep()` through a pure `keep_from(Option<OsString>)`; add a CLI test that a failing stand-in with `LOGAN_KEEP_SCRATCH=1` leaves the `.queries.fa`. Update `docs/walker.md:268` and `scratch.rs:13` ("any value other than empty or 0").

### WA-14 · low · Multi-anchor mode re-logs the full reference-name list for every accession

`walker/src/stages/align.rs:190` · efficiency · confirmed (reproduced) · fix effort: trivial · sources: R4-11

**Problem.** `walker/src/stages/align.rs:188-193` prints `[ACC] multi-anchor: aligning against N references [names…]` for every accession, although the list is constant for the run and `walker/src/main.rs:544` already announces multi-anchor mode. With a 160-record panel each line is 3,241 bytes, about 65 MB of identical stderr over 20,000 accessions, burying per-accession warnings.

**Proposed fix.** Drop the per-accession line (or keep `[ACC] multi-anchor: N references` without the list) and print the names once at startup inside main.rs's `if queries.len() > 1` block:
```rust
let names: Vec<&str> = queries.iter().map(|q| q.record_id.as_str()).collect();
eprintln!("[logan-walker] multi-anchor references: {}", names.join(", "));
```

**Tests and docs.** Update the comment at `align.rs:186-187`; a CLI test with a 2-record `-q` can assert the names appear once.
