# Walker: graph parsing, anchoring, subgraph extraction and phasing (parse.rs, types.rs, anchor.rs, extract.rs, phase.rs)

The core of this area is sound. The CSR parser rejects gapped ids, non-ACGT bases and u32 overflow without wrapping or panicking. Strand bits are kept in the graph, so each bubble's mirror traversal is deduplicated correctly. No defect was found in the anchoring or BFS logic itself; the entries there concern speed, memory and stale figures. The weak points are at the edges. The parser trusts its input too much: it reads only one spelling of the abundance tag, accepts NaN and inf, accepts files with no links or cut at a zstd frame boundary, and puts no limit on decoded size. Data-dependent conditions (an empty Logan assembly, more than 4.29 Gbp of sequence) are reported as `error` rather than as a negative status. `bubbles.json` counts flank pairs instead of variant sites, keeps error-level branches on deep samples and admits foldbacks. Many comments and docs still quote k=31 figures, an old guard and "~38 KB" GFAs. Most fixes follow one pattern: validate at parse time, classify data-driven failures by cause, and gate `bubbles.json` with the same depth-relative floor the rest of the walker uses. Only three entries can change the primary haplotype, VCF or SNV table: WG-02 and WG-07 (Logan format variants) and WG-05 (non-default flags).

| ID | Severity | Where | Title | Fix effort |
|---|---|---|---|---|
| WG-01 | medium | `walker/src/stages/parse.rs:206` | Empty Logan assemblies (55-byte objects) are reported as `error`, defeating `--allow-empty` | small |
| WG-02 | medium | `walker/src/stages/parse.rs:144` | Parser ignores Logan's documented `km:f` abundance alias; depth 0 silently disables majority and variant logic | trivial |
| WG-03 | medium | `walker/src/stages/phase.rs:184` | `bubbles.json` is gated only by an absolute 5x floor on the unpruned subgraph; deep samples get noise bubbles | small |
| WG-04 | medium | `walker/src/stages/phase.rs:129` | `bubbles.json` repeats one allele set once per flank pair, inflating bubble counts | small |
| WG-05 | medium | `walker/src/stages/extract.rs:159` | Prune floor is `max(prune_frac x depth, --min-abund)`: the allele gate moves the primary, undocumented | small |
| WG-06 | medium | `walker/src/stages/parse.rs:92` | Decoded size is unbounded: a crafted small cache file defeats `--max-size-mb` and the printed memory bound | small |
| WG-07 | medium | `walker/src/stages/parse.rs:154` | Parser silently accepts unitig files with no `L:` links; the walk runs on isolated unitigs | trivial |
| WG-08 | low | `walker/src/stages/parse.rs:137` | Undocumented 4.29 Gbp parser ceiling: larger accessions fail as `error`, not `too_large` | trivial |
| WG-09 | low | `walker/src/stages/extract.rs:132` | Non-finite `ka:f` values accepted; NaN can panic the median sort, inf gives invalid JSON | trivial |
| WG-10 | low | `walker/src/stages/parse.rs:209` | A cached unitig file cut exactly at a zstd frame boundary parses as a smaller, valid graph | trivial |
| WG-11 | low | `walker/src/stages/phase.rs:117` | Foldbacks and self-loops pass as two-allele bubbles | trivial |
| WG-12 | low | `walker/src/stages/extract.rs:232` | Neighbourhood sequence copy and tmpfs GFA are unbounded and outside the memory model; "~38 KB" is wrong | small |
| WG-13 | low | `walker/src/stages/parse.rs:80` | glibc arena retention after the first accession inflates per-worker memory beyond the 1.7x calibration | trivial |
| WG-14 | low | `walker/src/stages/extract.rs:6` | Subgraph-size figures, the guard description and the "5-hop" warning are stale against the defaults | small |
| WG-15 | low | `walker/src/stages/parse.rs:124` | Parse error names a nonexistent fallback; BCALM2 misattribution; stale anchor.rs input list | trivial |
| WG-16 | low | `walker/src/types.rs:178` | Comments name `subgraph_fraction` as the guard, call multi-anchor unfinished, and give `pident` the wrong unit | trivial |
| WG-17 | low | `walker/src/types.rs:462` | `AccessionResult` stores the primary haplotype three times; writers read different copies | small |
| WG-18 | low | `walker/src/stages/parse.rs:182` | `parse::build` packs bases one at a time; a lookup table with 4-base packing halves parse time | small |
| WG-19 | info | `walker/src/stages/anchor.rs:284` | Anchor scan pays an unpredictable branch per base; a bitmap prefilter is about 2.5x faster | small |

### WG-01 · medium · Empty Logan assemblies (55-byte objects) are reported as `error`, defeating `--allow-empty`

`walker/src/stages/parse.rs:206` · correctness · confirmed (reproduced) · fix effort: small · sources: SWR1-01

**Problem.** For accessions whose Logan assembly is empty, S3 returns HTTP 200 with a valid 55-byte zstd object (an empty frame plus the seekable-format footer). The fetcher rejects only 0-byte payloads (`walker/src/stages/stream.rs:96`, `:219`), the decoder yields nothing, and `parse::build` fails with `return Err(bad("no unitig records found".into()))` (`walker/src/stages/parse.rs:206-207`). `compute_worker` maps every `run_one` error to `AccessionResult::failed` (`walker/src/pipeline.rs:482`), so a property of the data becomes `status: error`.

**Impact.** About 1 in 100 to 1 in 400 accessions: a listing of `u/ERR1000*` showed 10 of 1,000 objects at exactly 55 bytes, and 14 of about 5,000 across eight prefixes. Every cohort summary counts them as errors. A shard where everything else found nothing exits 1 under `--allow-empty` instead of 0, which fails a Galaxy job, and a shard made only of such accessions exits 3.

**Evidence.** Put `https://logan-pub.s3.amazonaws.com/u/ERR10000001/ERR10000001.unitigs.fa.zst` (55 bytes) in a cache directory and run `logan-walker -a ERR10000001 -q q.fa --unitigs-cache cache --offline --allow-empty --graphaligner ./ga.sh`. Output: `ERROR: [ERR10000001] parse failed: no unitig records found`, `--allow-empty does not apply to errors`, exit 3.

**Proposed fix.** Distinguish "decoded to zero bytes" from "non-empty stream without records":
```rust
// error.rs
EmptyAssembly { accession: String },  // "[acc] Logan's unitig file for this accession is empty"
// parse.rs build(): count bytes read; after the loop
if bytes_seen == 0 { return Err(Error::EmptyAssembly { accession: accession.clone() }); }
if expected_id == 0 { return Err(bad("no unitig records found".into())); }   // unchanged
// pipeline.rs compute_worker
.unwrap_or_else(|e| match e {
    Error::EmptyAssembly { .. } => AccessionResult::empty(index, accession, Outcome::NoLoganData,
        Some("Logan's unitig file for this accession is empty".into())),
    e => AccessionResult::failed(index, accession, e.to_string()),
})
```
Make the `NoLoganData` console line (`walker/src/pipeline.rs:407`) print `r.error` when set. With this patch the run above reports `no_logan_data` and exits 0. A separate `empty_assembly` status is the alternative if downstream tools must tell 404 from empty. WG-08 changes the same line in `pipeline.rs`; implement both as one error classifier.

**Tests and docs.** A parse unit test with the 55-byte object's bytes asserting `Error::EmptyAssembly`; keep `build_str("")` an error. A `tests/cli.rs` case: `--allow-empty`, exit 0, status `no_logan_data`. `docs/walker.md:300`: "Logan has no unitigs for the accession (HTTP 404, or an empty assembly)"; same for the doc comments at `error.rs:24-26` and `types.rs:415`.

### WG-02 · medium · Parser ignores Logan's documented `km:f` abundance alias; depth 0 silently disables majority and variant logic

`walker/src/stages/parse.rs:144` · correctness · confirmed (reproduced) · fix effort: trivial · sources: R2-03, SC1-05

**Problem.** Abundance is read only from `if f.starts_with(b"ka:f:") {` (`walker/src/stages/parse.rs:144`); a record without that tag keeps `let mut ab = 0.0f32;` (`:142`) with no warning. Logan's format page (IndexThePlanet/Logan `Unitigs.md`, line 17) says that "in some large accessions, the `ka:f:[xx]` field is replaced by `km:f:[xx]` ... consider that they are the same information". For such an accession every unitig has depth 0, so the prune, majority correction, bubble detection and variant calling are all disabled, and the run still reports `status: success`.

**Impact.** The primary follows the GraphAligner route even where that route carries the minority allele; `alt_haplotypes.fasta` and the VCF are empty; the summary shows `"mean_ka": 0.0`. The only trace is a stderr line `prune: depth 0x ... skipping the prune`. Frequency is unknown: about 240 real objects probed contained none, and the default `--max-size-mb 1024` filters out many large objects.

**Evidence.** One diamond (major allele C at ka 60, minor G at ka 40, G is the reference) written once with `ka:f` and once with `km:f`, run `--offline` with a stand-in GraphAligner routed through the minor allele. `ka:f`: depth 100x, site corrected to C, `bubbles=1`, one alt. `km:f`: `depth 0x`, `bubbles=0`, `ka=0.0`, empty alts, primary keeps G. With the patch below the `km:f` output equals the `ka:f` output.

**Proposed fix.**
```diff
-                if f.starts_with(b"ka:f:") {
+                // Logan Unitigs.md: km:f is a known alias of ka:f in some large accessions.
+                if f.starts_with(b"ka:f:") || f.starts_with(b"km:f:") {
                     ab = std::str::from_utf8(&f[5..])
                         .ok()
                         .and_then(|s| s.parse::<f32>().ok())
+                        .filter(|v| v.is_finite() && *v >= 0.0)
```
The `filter` line also fixes WG-09. Second step: count records without an abundance tag and, if any, print `[warn] [ACC] N of M unitigs have no ka:f/km:f abundance; depth-based stages are disabled` and add `"abundance_missing": N` to the summary. Do not make it a hard error; many unit fixtures use tagless headers.

**Tests and docs.** Unit test `km_f_is_read_as_abundance` (`km:f:60.0` gives 60.0; `nan`, `inf`, `-7.5`, `1e40` are errors). Optionally a CLI case with a `km:f` fixture asserting `bubbles=1`. Mention the alias at `docs/walker.md:15` and `:333`.

### WG-03 · medium · `bubbles.json` is gated only by an absolute 5x floor on the unpruned subgraph; deep samples get noise bubbles

`walker/src/stages/phase.rs:184` · science · confirmed (reproduced) · fix effort: small · sources: R3-01

**Problem.** `phase` keeps an allele when `ab >= cfg.min_abund && vaf >= cfg.min_vaf` (`walker/src/stages/phase.rs:184`), where `min_abund` is an absolute 5x and the VAF is computed only among the bubble's own branches. It is called once, before alignment and before the prune (`walker/src/pipeline.rs:600`), so neither depth-relative floor (`--prune-frac`, `--min-abund-frac`) reaches it. `config.rs:148-151` itself says 5x "admits branches carried by a handful of erroneous reads" on a 15,000x sample. The module doc (`phase.rs:5-7`) still says such a fork "represents two populations genuinely present in one sample", and the tools-iuc PR #8800 wrapper help presents the output as the fork structure along the walked path.

**Impact.** On SRR8209082 (HIV-1 pol, about 15,000x; also the Galaxy wrapper's test accession) `bubbles.json` holds 3,632 bubbles. In every one, every allele is at most 66.3x (median strongest allele 11.9x), under 0.5% of gene depth; the real 3,600x second population is not among them. A user reading the file as quasispecies evidence reads noise. The primary, VCF and SNV table are unaffected, and `docs/walker.md:144-146` already steers users to the VCF.

**Evidence.** `logan-walker -a SRR8209082 -q queries/query_hiv1_pol.fasta --unitigs-cache cache --offline --graphaligner <stand-in replaying the route>` prints `bubbles=3632`; adding `--min-abund-frac 0.01` gives a byte-identical `bubbles.json`.

**Proposed fix.** Add `phase_with_floor(cfg, graph, sub, floor)` gating on `ab >= floor.max(cfg.min_abund)`, keep `phase` as a wrapper, and move the call after the prune block in `pipeline.rs` with:
```rust
let bubble_floor = final_path.as_ref().map_or(cfg.min_abund, |p| {
    let depth = stages::extract::path_depth(&graph, &p.steps);
    (cfg.prune_frac.max(cfg.min_abund_frac) as f32 * depth).max(cfg.min_abund)
});
```
With this patch the default run gives `bubbles=0`, the golden flags (`--prune-frac 0 --no-prefer-majority --no-emit-alts`) give a byte-identical `bubbles.json`, and the 9 phase tests pass. A `no_path` run still reports all bubbles because there is no depth estimate. Optionally add `on_route: bool` to `Bubble`. See WG-04 and WG-11 for the other two `bubbles.json` defects.

**Tests and docs.** Unit test: a diamond of 6x and 8x branches beside a 15,000x route is dropped at `prune_frac` 0.01 and kept when both fractions are 0. Update `phase.rs:5-7`, the walker.md outputs table and `--prune-frac`/`--min-abund-frac` help, the "14 bubbles" in `docs/issue_intra_accession_variance_lost.md:87`, the wrapper help, and the changelog.

### WG-04 · medium · `bubbles.json` repeats one allele set once per flank pair, inflating bubble counts

`walker/src/stages/phase.rs:129` · correctness · confirmed (reproduced) · fix effort: small · sources: R3-03

**Problem.** `for &merge in &reconverged {` (`walker/src/stages/phase.rs:129`) emits one `Bubble` per (fork, merge) pair, and the only dedup, `if (f, merge) > (merge ^ 1, f ^ 1) { continue; }` (`:149`), removes only the mirror traversal. The same alleles reachable between several forks or merges are emitted once per flank pair with identical branches and VAFs. `n_bubbles`, the summary `bubbles` field, `bubble_events` and the DONE line therefore count flank pairs, and the inflation depends on local topology.

**Impact.** SRR8209082 with defaults: 3,632 bubbles but 317 distinct allele sets (11.5x), one set repeated 16 times. Two SNPs exactly k bases apart form a bow-tie `0 -> {1,2} -> {3,4} -> 5` and are reported as 4 bubbles. Comparing bubble counts across samples as a diversity signal, the use `phase.rs:19-20` states, is misleading. VCF, SNV table and alts are unaffected.

**Evidence.** A unit fixture for the bow-tie returns 4 bubbles (`0+->3+` and `0+->4+` with identical branches, `1+->5+` and `2+->5+` likewise). A count of distinct oriented allele sets over SRR8209082's subgraph gives 317.

**Proposed fix.** Key records on the mirror-canonical allele set:
```rust
let mut a = in_bubble.clone(); a.sort_unstable();
let mut m: Vec<u32> = a.iter().map(|x| x ^ 1).collect(); m.sort_unstable();
let key = a.min(m);
```
Collect into a `BTreeMap<Vec<u32>, Bubble>`, add `flanks: Vec<(UnitigId, Strand, UnitigId, Strand)>` to `Bubble` (written as `"flanks"`), and keep the smallest pair as `fork`/`merge`. A keep-first version gave bow-tie 4 -> 2, SRR8209082 3,632 -> 317, ERR2028820 1 -> 1.

**Tests and docs.** Rewrite `a_fork_with_two_merges_reports_both` to expect one record with two flanks; add the bow-tie test expecting 2. `docs/walker.md:107`: one record per allele set. State the changed meaning of `n_bubbles` and the schema change in the release notes.

### WG-05 · medium · Prune floor is `max(prune_frac x depth, --min-abund)`: the allele gate moves the primary, undocumented

`walker/src/stages/extract.rs:159` · usability · confirmed (reproduced) · fix effort: small · sources: R3-04, IC-07, RC-07

**Problem.** `depth_floor` returns `(cfg.prune_frac as f32 * depth).max(cfg.min_abund)` (`walker/src/stages/extract.rs:159`). Every user-facing description says otherwise: `--prune-frac` "drop unitigs below this fraction of the gene's depth" (`main.rs:77`, `docs/walker.md:204`, wrapper help) and `--min-abund` "Minimum ka:f for a variant" (`main.rs:67`, `docs/walker.md:211`). `pipeline.rs:717-718` says "a 180x isolate gets a 1.8x floor and nothing is removed" and `config.rs:124` says "0 nodes lost"; the floor is 5.0, and ERR2028820 keeps 479 of 747 subgraph nodes. The doc comment at `extract.rs:156` ("can only tighten the existing floor") is also wrong.

**Impact.** Raising `--min-abund` to be stricter about variants raises the prune floor, removes unitigs on the first route, and can replace the primary haplotype. Any `--prune-frac` below 5/depth is the same setting, so a sweep from 0.0001 to 0.02 is a no-op on genes under 500x. Default runs are as intended.

**Evidence.** Offline fixture at 60x with a SNP bubble (ref 18x, alt 42x), `--no-prefer-majority`: floor 5.0, nothing removed, primary carries the reference allele. Adding `--min-abund 20`: floor 20.0, primary carries the other allele. An 8x fixture gives floor 5.0 for `--prune-frac` 0.0001, 0.01 and 0.5.

**Proposed fix.** Decouple, keeping today's defaults:
```rust
// config.rs
pub prune_min_abund: f32,                       // default 5.0
// extract.rs depth_floor
(cfg.prune_frac as f32 * depth).max(cfg.prune_min_abund)   // else cfg.prune_min_abund
// main.rs: "--prune-min-abund" => cfg.prune_min_abund
// output.rs settings: add ("prune_min_abund", ...)
```
With this patch all 222 unit tests pass and `--min-abund 20` leaves the primary unchanged. If decoupling is rejected, fix the docs at least.

**Tests and docs.** New test `min_abund_does_not_move_the_prune_floor`; update `the_floor_is_off_when_prune_frac_is_zero`, `a_gene_below_min_abund_would_empty_the_subgraph` and the settings test. Docs: USAGE, `docs/walker.md:204` plus a new row, `pipeline.rs:717-719`, `config.rs:124`, `extract.rs:156`, composite and wrapper help.

### WG-06 · medium · Decoded size is unbounded: a crafted small cache file defeats `--max-size-mb` and the printed memory bound

`walker/src/stages/parse.rs:92` · robustness/security · confirmed (reproduced) · fix effort: small · sources: RX-02

**Problem.** `--max-size-mb` is checked only against the compressed size (`walker/src/stages/stream.rs:82`, `:189`, `:208`). `stream::decoder` (`stream.rs:298-309`) returns an unbounded zstd reader, and `parse::build` reads it with `rd.read_until(b'\n', &mut line)` (`walker/src/stages/parse.rs:92-94`) with no line-length limit. The u32 base-count guard runs only at the next header (`parse.rs:137-140`). Header bytes are copied whole into error messages (`parse.rs:115-119`), which go to stdout and the summary JSON. The startup "memory bound ... worst case" (`main.rs:497`) assumes real Logan decode ratios (3.5-5.5x).

**Impact.** Real S3 objects do not trigger it. A crafted or corrupt cache entry does, for example a Galaxy dataset fed to `--unitigs-cache` or a shared cache directory: one file under the 1 GiB cap can expand to terabytes, per worker, and OOM the job or node. A 1 GiB sequence line yields `no_anchors` with exit 0; a 256 MiB header yields a 268 MB stdout line and a 268 MB summary JSON.

**Evidence.** A valid header followed by 1 GiB of `A` (32,813 bytes compressed), run with `-w 1 --max-size-mb 1 --offline`: startup prints `0.0 GB worst case`, peak RSS 1,292 MB, status `no_anchors`. A 256 MiB header with no newline (8,212 bytes): peak 524 MB, stdout 268,435,590 bytes.

**Proposed fix.**
1. In `stream::decoder`, wrap the reader in a budgeted `Read` (budget `max(64 MiB, 12 x compressed_len)`) returning `InvalidData("decode budget exceeded")`, mapped to `Error::TooLarge`.
2. In `parse::build`, replace `read_until` with a `fill_buf`-based `read_line_capped(rd, &mut line, 64 MiB)`, or pack sequence straight from `fill_buf()`.
3. Pass untrusted bytes in errors through `clip()` (80 bytes plus "... (N bytes)").
4. Optional: enforce a graph-byte budget (`GRAPH_PEAK_PER_COMPRESSED_BYTE x compressed_len`) in the loop, so record-count bombs are caught too.

With steps 1-3 applied both fixtures fail in 0.3 s at 76-77 MB peak. Real ratios measured are 4.5x and 5.5x, so 12x leaves about 2x margin; report the ratio in the `too_large` detail.

**Tests and docs.** Parse tests: a 100 MiB line without newline is rejected with a short message; a long header is clipped. A stream test: a frame expanding past 12x gives `Outcome::TooLarge`. State the 12x assumption in the walker.md memory section and `pipeline.rs:60-90`.

### WG-07 · medium · Parser silently accepts unitig files with no `L:` links; the walk runs on isolated unitigs

`walker/src/stages/parse.rs:154` · data/correctness · plausible · fix effort: trivial · sources: R2-05

**Problem.** `parse::build` adds edges only for header fields that `starts_with(b"L:")` (`walker/src/stages/parse.rs:154`), and its only post-loop check is `expected_id == 0` (`:206`). A file in which no record carries a link parses into an edgeless graph without a warning. The BFS then adds no neighbours, the GFA has `S` lines and no `L` lines, and the accession ends as `no_path` or a one-unitig fragment, which reads as a biological negative. `docs/walker.md:335` states that every edge appears in both endpoints' headers; the code relies on it without checking.

**Impact.** Such objects exist on Logan: the first 256 KiB of ERR3988483 (2.53 GB) has 0 of 15,864 headers with links (and headers without the accession prefix, `>0 ka:f:2.2`), ERR3988478 (9.56 GB) 0 of 21,370; control SRR8209082 has 9,752 of 10,070. Both are over the default cap (`too_large`) and, with the cap raised, over the 4.29 Gbp guard (WG-08), so no silent false negative was demonstrated. It becomes one for any link-less object within the limits, or once backlog item 1(b) removes the size cap.

**Evidence.** A three-record fixture without `L:` tokens runs to `no_path`, exit 0, GFA with 2 `S` and 0 `L` lines, no warning; the same records with links give 3 `S` and 4 `L` lines.

**Proposed fix.** After the `expected_id == 0` check:
```rust
if targets.is_empty() && expected_id >= MIN_RECORDS_FOR_LINK_CHECK /* 1000 */ {
    return Err(bad(format!("no L: link tokens in any of {expected_id} records: \
        this object has no adjacency and cannot be walked")));
}
```
The threshold keeps tiny link-less fixtures (`parse.rs:323`) valid. A dedicated `Outcome::NoLinks` (status `no_links`) is the better long-term choice so cohort reports list these separately.

**Tests and docs.** Parse test: 1,000 link-less records give the new error; 2 still parse. Document the link-less variant and prefix-less headers in `docs/walker.md:335`; note it under backlog 1(b).

### WG-08 · low · Undocumented 4.29 Gbp parser ceiling: larger accessions fail as `error`, not `too_large`

`walker/src/stages/parse.rs:137` · reliability · confirmed (reproduced) · fix effort: trivial · sources: R2-06, BL1-07

**Problem.** `seq_offset` is `Vec<u32>` (`walker/src/types.rs:126`), guarded by `u32::try_from(total_bases).map_err(|_| bad("total sequence exceeds u32 addressing".into()))?` (`walker/src/stages/parse.rs:137-140`, again at `:213-215`). The guard is correct, but `bad` builds `Error::Parse` and `pipeline.rs:482` makes it `error`. Raising `--max-size-mb` past roughly 1,500-2,400 MB (the retry backlog item 1(c) proposes) therefore downloads and parses up to 4.29 Gbp, then reports `error`; a single-accession job exits 3 even with `--allow-empty`. Neither `--help` (`main.rs:76`) nor `docs/walker.md:197` mentions the ceiling. A synthetic 2^32 + 2^20-base fixture reproduces this offline.

**Proposed fix.** Return `Error::TooLarge` and classify compute errors like fetch errors:
```rust
let too_big = || Error::TooLarge { accession: accession.clone(),
    detail: format!("more than {} bases of unitig sequence, the parser's limit; \
                     --max-size-mb cannot raise it", u32::MAX as u64 + 1) };
seq_offset.push(u32::try_from(total_bases).map_err(|_| too_big())?);   // both sites
// pipeline.rs:482
.unwrap_or_else(|e| classify_failure(index, accession, e.to_string(), &e))   // renamed `unfetched`
```
Combine with WG-01, which edits the same line. Warn at startup when `--max-size-mb` exceeds about 1,500.

**Tests and docs.** Add `build_with_limit(acc, reader, max_bases)` so a unit test can hit the limit with a 10-base cap; extend `a_fetch_failure_is_classified_by_its_cause` (`pipeline.rs:968`). Update `main.rs:76` and `docs/walker.md:197`, `:301` (`too_large` is no longer always "skipped before downloading").

### WG-09 · low · Non-finite `ka:f` values accepted; NaN can panic the median sort, inf gives invalid JSON

`walker/src/stages/extract.rs:132` · reliability · confirmed (reproduced) · fix effort: trivial · sources: R3-10

**Problem.** The parser accepts anything `f32::from_str` accepts (`walker/src/stages/parse.rs:147`), so `ka:f:NaN`, `ka:f:inf` and `ka:f:-3` load. Both median helpers sort with `ab.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));` (`walker/src/stages/extract.rs:132`, `walker/src/stages/variants.rs:546`), which is not a total order once NaN is present; since Rust 1.81 the sort may panic, and `catch_unwind` (`pipeline.rs:480`) turns the accession into "internal error: a stage panicked". One `inf` on the route makes `"mean_ka": inf` in summary.json (`walker/src/io/output.rs:271`), which is not JSON. Needs corrupted or hand-made input; Logan has not been seen to emit such values. A standalone program with this comparator panics on 21-200 values containing 3 or 9 NaNs.

**Proposed fix.** Reject at the source with `.filter(|v| v.is_finite() && *v >= 0.0)` after the parse (the same line as WG-02's fix), and make both sorts total with `ab.sort_by(f32::total_cmp);`.

**Tests and docs.** Next to the `ka:f:notanumber` test (`parse.rs:348`) assert `NaN`, `inf`, `-1` are errors; a `path_depth` test over 64 steps with several NaNs does not panic.

### WG-10 · low · A cached unitig file cut exactly at a zstd frame boundary parses as a smaller, valid graph

`walker/src/stages/parse.rs:209` · robustness · confirmed (reproduced) · fix effort: trivial · sources: R2-04

**Problem.** Logan objects contain several zstd frames. A file ending exactly after one frame decodes without error, and `parse::build` sets `n_unitigs = expected_id` (`walker/src/stages/parse.rs:209`) without checking that every `L:` target is below it. The run reports `no_anchors` or walks a partial graph, exit 0 under `--allow-empty`, no warning; `write_gfa` drops dangling links silently (`extract.rs:332`). The decoder doc (`stream.rs:290-291`) claims truncation always surfaces as an error. The walker cannot write such a file itself (atomic rename, ureq rejects short bodies), so it needs a hand-filled cache and a cut that lands on a frame boundary. A synthetic two-frame file cut after frame 1 reproduced it (`no_anchors`, 1,000 of 2,000 unitigs).

**Proposed fix.** After `parse.rs:209`:
```rust
if let Some(&p) = targets.iter().find(|&&p| (p & ID_MASK) as usize >= n_unitigs) {
    return Err(bad(format!(
        "link to unitig {} but the file holds only {n_unitigs} records; \
         the unitig file is truncated or incomplete", p & ID_MASK)));
}
```

**Tests and docs.** `file_cut_at_a_frame_boundary_is_an_error` with two `zstd::encode_all` frames. Correct `stream.rs:290-291` and `:68-69`.

### WG-11 · low · Foldbacks and self-loops pass as two-allele bubbles

`walker/src/stages/phase.rs:117` · correctness · confirmed (reproduced) · fix effort: trivial · sources: R3-06

**Problem.** The reconvergence filter excludes only the fork's own oriented node, `.filter(|&(c, count)| count >= 2 && c != fork as u32)` (`walker/src/stages/phase.rs:117`), and `in_bubble` (`:159-163`) is never checked for distinct unitigs. A unitig end linking to both strands of one neighbour (an inverted repeat) yields a bubble whose two alleles are the same unitig at VAF 0.5 each; a self-looping unitig can be both allele and merge. The data-backed test at `phase.rs:534` asserts the invariant the code does not enforce. ERR2028820's full graph has 61 foldback ends and 1,746 self-edge ends, none in the golden subgraphs. Only `bubbles.json` and bubble counts are affected (see WG-03, WG-04).

**Proposed fix.** After the `in_bubble.len() < 2` check:
```rust
let mut u: Vec<u32> = in_bubble.iter().map(|&b| b >> 1).collect();
u.sort_unstable(); u.dedup();
if u.len() < in_bubble.len() || u.contains(&(f >> 1)) || u.contains(&(merge >> 1)) {
    continue;
}
```
Do not reject `merge == f ^ 1` outright: a hairpin with two distinct arms is a real bubble.

**Tests and docs.** Foldback and self-loop fixtures expecting no bubble; a two-arm hairpin expecting one. Change the assertion at `phase.rs:518/534` to compare oriented nodes and add a distinct-alleles check; extend the `Bubble` doc (`types.rs:349-358`).

### WG-12 · low · Neighbourhood sequence copy and tmpfs GFA are unbounded and outside the memory model; "~38 KB" is wrong

`walker/src/stages/extract.rs:232` · memory · confirmed (reproduced) · fix effort: small · sources: R3-05

**Problem.** `collect_seqs` decodes every neighbourhood unitig to one byte per base (`walker/src/stages/extract.rs:232-240`) and keeps the copy in `Subgraph.seqs`; `write_gfa` writes all of it to `cfg.scratch_dir` (tmpfs by default) with no cap. `Config::max_resident_bytes` (`config.rs:432-434`) counts neither, and `exceeds_guard` only warns. Five comments say the GFA "is ~38 KB" (`extract.rs:19`, `:278`, `types.rs:193`, `config.rs:364`, `align.rs:11`). On SRR8209082 at defaults the neighbourhood is 42% of the graph and the GFA is 4.86 MB, 2.6x the compressed input. Larger targeted libraries with `-w 8` on a 64 MB `/dev/shm` (Docker's default) can hit ENOSPC; not measured at scale.

**Proposed fix.** 1. Correct the five comments to the measured range (about 80 KB on ERR2028820 to 4.9 MB on SRR8209082), scaling with the neighbourhood. 2. Add `--max-subgraph-bytes` (default well above 4.9 MB, e.g. 256 MB): sum S-line bytes from `seq_offset` before writing in `write_gfa_file`, fail the accession with a clear error when over, and add `compute_workers x cap` to `max_resident_bytes` and the startup line. 3. Optional: drop `Subgraph.seqs` and decode on demand in `stitch.rs:80`, `variants.rs:595/642`, `pipeline.rs:783`.

**Tests and docs.** `write_gfa` with a cap below the subgraph size errors and leaves no file; a `max_resident_bytes` test. Update the walker.md memory section and `--scratch-dir` help (`main.rs:116`).

### WG-13 · low · glibc arena retention after the first accession inflates per-worker memory beyond the 1.7x calibration

`walker/src/stages/parse.rs:80` · memory · confirmed (reproduced) · fix effort: trivial · sources: RF-05

**Problem.** The five CSR arrays start empty and grow by doubling (`walker/src/stages/parse.rs:80-84`). After the first accession frees its large arrays, glibc raises its dynamic mmap threshold (up to 32 MiB), so later arrays grow inside the heap arenas and leave freed chunks resident. `GRAPH_PEAK_PER_COMPRESSED_BYTE = 1.7` (`config.rs:445`) was measured "one accession per run" (`config.rs:439-440`). On 24 cached 44 MiB accessions with `-w 4 --fetch-workers 1`, peak RSS was 708 MiB by default and 564 MiB with `MALLOC_MMAP_THRESHOLD_=131072`, about 25% waste; the printed bound (0.7 GB) was only marginally exceeded.

**Proposed fix.** In `main()` before any thread starts:
```rust
#[cfg(all(target_os = "linux", target_env = "gnu"))]
unsafe { libc::mallopt(libc::M_MMAP_THRESHOLD, 1 << 20); }
```
(`libc` is already a dependency; this gave 538 MiB.) Optionally reserve array capacity from the compressed length. Re-measure the constant on multi-accession runs.

**Tests and docs.** An `#[ignore]` or bench script asserting peak RSS stays under the printed bound over N accessions. Update `config.rs:438-444`, `docs/walker.md:238-250`, and `pipeline.rs:78` ("~1.6x" against the 1.7x constant).

### WG-14 · low · Subgraph-size figures, the guard description and the "5-hop" warning are stale against the defaults

`walker/src/stages/extract.rs:6` · docs · confirmed (reproduced) · fix effort: small · sources: R3-09, BL1-09, R2-09, SWR2-03

**Problem.** `extract.rs` quotes three different sets of neighbourhood sizes: lines 6-7 (356 and 26,323 nodes, "about a third"), lines 49-52 and the table at 199-200 (422 and 8,466, "about a tenth"). The second set is what `-k 31` gives; none matches the default `k: 21` (`config.rs:451`), and `anchor.rs:7-17` has the same k=31 figures. `docs/walker.md:288` says the walker "warns whenever a subgraph exceeds 1% of the graph"; the guard is 3 unitigs per reference base (`extract.rs:630`, `docs/walker.md:121-128`). The warning at `walker/src/pipeline.rs:580` hard-codes "the 5-hop neighbourhood" whatever `--hops` is, as does the doc at `types.rs:187`. `pipeline.rs:780-781` and `:840-841` describe behaviour the code no longer has. Readers sizing scratch space or the stream-decode work (backlog (b)) expect under 2% of the graph; SRR8209082 keeps 42%.

**Evidence.** Offline runs with a stand-in aligner: ERR2028820/pfcrt k=21 gives 220 anchors and 747 nodes, k=31 gives 77 and 422; SRR8209082/HIV pol k=21 gives 12,876 and 34,582 (42.4%), k=31 gives 92 and 8,466. `--hops 2` prints `WARNING: the 5-hop neighbourhood ...`.

**Proposed fix.** Replace the extract.rs and anchor.rs figures with one table naming k (k=31: 77/422 and 92/8,466; k=21: 220/747 of 1,277,362 and 12,876/34,582 of 81,534). `docs/walker.md:286-289`: "about 42% at the default k=21" and "flagged `provisional` above 3 unitigs per reference base". `bench/check_golden.sh:19`: "~35k-node (k=21) ... ~42%". Correct the comments at `pipeline.rs:780-781` and `:840-841`. The "~38 KB" comments are covered by WG-12. In the warning:
```diff
-            "[{}] WARNING: the 5-hop neighbourhood holds ...
+            "[{}] WARNING: the {}-hop neighbourhood holds ...
             graph.accession,
+            cfg.hops,
```

**Tests and docs.** A CLI or pipeline test with `--hops 2` asserting "2-hop". Drop "still accurate" and "~99.9%" from the backlog and plan.

### WG-15 · low · Parse error names a nonexistent fallback; BCALM2 misattribution; stale anchor.rs input list

`walker/src/stages/parse.rs:124` · docs · confirmed (reproduced) · fix effort: trivial · sources: R2-09

**Problem.** A file with a gap in unitig ids fails with "... CSR construction relies on record order; this file needs the fallback path." (`walker/src/stages/parse.rs:123-124`); no fallback exists. `types.rs:12-13` says "21 bits covers the range ... fall back to a map if it ever fails", while ids are packed in 30 bits (`parse.rs:36`) and a gap is a hard error. `types.rs:11`, `types.rs:110`, `config.rs:61` and `docs/walker.md:326` attribute the unitigs to BCALM2; Logan's `Unitigs.md` says they were built with a modified Cuttlefish2 and carry BCALM2-style links. `anchor.rs:19` lists "decompressed unitig bytes" as an input that `anchor::find` does not take. The k=31 figures in extract.rs and anchor.rs are covered in WG-14.

**Proposed fix.** Message: "unitig ids are not dense and sequential: expected {expected_id}, found {id}. The walker builds its graph in record order and cannot read this file." `types.rs:11-13`: ids are packed in 30 bits (`parse::ID_BITS`) and `parse::build` rejects gaps. Replace "BCALM2" with "Logan's modified Cuttlefish2 (BCALM2-style `L:` links)" at the four sites. Delete "decompressed unitig bytes," from `anchor.rs:19`.

**Tests and docs.** Optionally assert in the gap-id test that the message does not mention "fallback".

### WG-16 · low · Comments name `subgraph_fraction` as the guard, call multi-anchor unfinished, and give `pident` the wrong unit

`walker/src/types.rs:178` · docs · confirmed (by reading) · fix effort: trivial · sources: RS-11, SWR3-01, SWR2-03

**Problem.** `types.rs:178-179` says `extract::subgraph_fraction` "is what `pipeline::run_one` actually checks", and `extract.rs:188` calls it "the guard signal"; the gate is `stages::extract::exceeds_guard` (`pipeline.rs:577`), and `config.rs:298` and `anchor.rs:312` already say `subgraph_fraction` is not a gate. `pipeline.rs:273` ("handed the first [`Query`] for now"), `pipeline.rs:297` and `input.rs:559` ("Phase 4") and `input.rs:4` ("single record") describe multi-anchor as unfinished, though `pipeline.rs:613-634` aligns every query and runs `pick_primary`. `types.rs:269-270` documents `Haplotype::pident` as `num_matches / aln_len`, a 0-1 ratio, while `stitch.rs:214` stores a percentage and `types.rs:387` says so; a threshold written from the field doc (`pident >= 0.95`) is off by 100x. `types.rs:174` and `anchor.rs:311` call `fraction_of_graph`/`shared_kmers` "diagnostics", but nothing outside anchor.rs reads them. User output and `docs/walker.md` are correct.

**Proposed fix.**
```diff
-    /// Identity of the GraphAligner alignment that produced this path: `num_matches /
-    /// aln_len` from the GAF.
+    /// Percent identity (0-100) of the GraphAligner alignment that produced this path:
+    /// `num_matches / aln_len * 100` from the GAF.
```
Point `types.rs:176-179` and `extract.rs:186-190` at `exceeds_guard` (unitigs per reference base against `Config::subgraph_nodes_per_base_guard`). Rewrite the multi-anchor comments to say align and stitch run over every query and `pick_primary` chooses the winner. Either emit `fraction_of_graph` in summary.json or delete it, `shared_kmers` and the `anchor.rs:311` sentence (and the fixture at `extract.rs:442-443`).

**Tests and docs.** Optional stitch test asserting `pident == 100.0` when `num_matches == aln_len`.

### WG-17 · low · `AccessionResult` stores the primary haplotype three times; writers read different copies

`walker/src/types.rs:462` · simplification · confirmed (by reading) · fix effort: small · sources: RS-10

**Problem.** `haplotype: Option<Haplotype>` (`walker/src/types.rs:462`), `haplotypes: Vec<Haplotype>` (built as `haplotype.clone().into_iter().collect()` at `pipeline.rs:894`, so never more than one entry) and `anchor_used` hold the same data, and `pipeline.rs:907` rebuilds `haplotype` from `haplotypes.first()`. `write_fasta` and `write_result` read `haplotypes` (`output.rs:74`, `:257-262`), `snv_sites` and `report` read `haplotype`, `reference_of` reads `anchor_used`. Production fills all three coherently, but a constructor or fixture that fills one copy produces a FASTA record with no VCF or SNV rows; the test fixtures already set both, and `variant_writers_use_the_winning_reference` (`output.rs:1041`) builds a state production never builds. The `haps=` branch at `pipeline.rs:389` is dead.

**Proposed fix.** Delete `haplotypes` and `anchor_used`; iterate `r.haplotype.iter()` in the writers, keeping the JSON `"haplotypes": [...]` array shape; `reference_of` reads `r.haplotype.as_ref().and_then(|h| h.anchor_used.as_deref())`; tag `haplotype` in place in multi-anchor mode; build the `NoAnchors` result at `pipeline.rs:549` as `AccessionResult { n_unitigs, timing, ..AccessionResult::empty(index, accession, Outcome::NoAnchors, None) }`.

**Tests and docs.** The existing FASTA and summary.json byte-identity tests pin the output; rewrite `variant_writers_use_the_winning_reference`; run one multi-anchor CLI case.

### WG-18 · low · `parse::build` packs bases one at a time; a lookup table with 4-base packing halves parse time

`walker/src/stages/parse.rs:182` · efficiency · confirmed (reproduced) · fix effort: small · sources: RF-04

**Problem.** For every base the sequence loop (`walker/src/stages/parse.rs:182-202`) runs a five-way `match`, a conditional push, a recomputed `seq_2bit.len() - 1` and a bounds-checked OR. Parse is the largest compute stage once the payload is local. On an ERR2028820-sized synthetic (1.28M unitigs, 136 MB decoded) `build` took 3.28 s and a table-driven version 1.34 s with byte-identical arrays (zstd decode alone: 0.46 s). For S3 runs the gain is CPU time, not wall time.

**Proposed fix.** `static LUT: [u8; 256]` (A/C/G/T to 0-3, else 4). Finish the partial byte while `total_bases & 3 != 0`; `reserve(rest.len() / 4 + 1)`; for each `chunks_exact(4)` look up four codes, report the first byte with the existing "non-ACGT base" error if `(a|b|c|d) > 3`, else `push(a | b<<2 | c<<4 | d<<6)`; pack the 0-3 remaining bases into one byte; update `total_bases`.

**Tests and docs.** Build the same records with line wraps of 1 to 9 bases and compare `seq_2bit`/`seq_offset` with a reference packer; add non-ACGT cases inside a chunk and in the remainder.

### WG-19 · info · Anchor scan pays an unpredictable branch per base; a bitmap prefilter is about 2.5x faster

`walker/src/stages/anchor.rs:284` · efficiency · confirmed (reproduced) · fix effort: small · sources: RF-06

**Problem.** `anchor::find` probes the query table for every base: `if j + 1 >= k && index.contains(fwd.min(rev))` (`walker/src/stages/anchor.rs:284`). For a single gene the table runs at about 31% load (PfCRT: 1,255 k-mers in 4,096 slots), so the "empty slot or keep probing" branch is unpredictable. Measured 6-14 ns per base against the "a few nanoseconds" in the module doc (`anchor.rs:48-51`); about 0.4-0.6 s of CPU per 68 Mbp accession. Results are correct.

**Proposed fix.** Build a 1-bit-per-position filter once in `index_query_set`, sized at 16 bits per table slot and indexed by high bits of the multiply hash (independent of the table's `>> 32 & mask` bits), and test it first:
```rust
let b = (canon.wrapping_mul(MULT) >> index.fshift) as usize;
if (index.filter[b >> 6] >> (b & 63)) & 1 != 0 && index.contains(canon) { count += 1; }
```
Roll the first k-1 bases in a separate loop to drop `j + 1 >= k` from the hot loop. On a synthetic accession this gave 2.5-2.7x on PfCRT with identical anchors; panels at about 5% load gain less. The filter must be built from every stored canonical k-mer, or anchors are lost silently.

**Tests and docs.** A property test that the new path returns the same `ids` and `shared_kmers` as the current loop on random graphs and query pools. Update `anchor.rs:48-51`.
