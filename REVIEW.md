# Review of disassembler 0.1.1

**Reviewed:** `nekrut/disassembler` at commit `25c27f2` (tags `v0.1.1` and `logan-walker-v0.1.1`), with the maintainers' Galaxy delivery plan (2026-11-04), their 0.2 backlog (2026-11-10) and tools-iuc PR #8800 at head `b9a16ad`.
**When:** 2026-11-10 and 2026-11-11. Nothing in the repository was changed.
**In this file:** the defects that matter most, a suggested order of work, the verdict on the backlog and the Galaxy plan, all 314 findings with links to the detailed sections in [`review/`](review/), the proposed patches, and how the review was done and what it could not check.

## What matters most
1. **Default reconstructions are not majority-corrected.** With default flags the walker discards its majority correction unless the GraphAligner row spans the whole stitched path, which the real rows recorded in the repository do not; where the aligner walked the rarer, reference-like strain, `_primary` publishes that strain and no output says so. This is in 0.1.0 and 0.1.1. WV-01 (critical); small fix, but it switches on swap code that has never run on real rows, so land it with a default-flags golden run (TC-02) and a panel re-run.
2. **Minority alleles in mixed and deep samples are misreported.** A walked allele at 1-10% of a fork is neither corrected nor reported, so the primary can carry a contaminant or minority strain at `pident=100` (WV-02). Logan's `ka:f` stops at 25,000, so on deep viral or amplicon samples minor-allele VAFs are inflated and pass `--min-vaf` as `PASS`; the shipped Zika example has `RO=25000` in all 11 VCF rows (WV-03). Both small.
3. **The published 0.1.1 container runs walker 0.1.0.** `biocontainers/disassembler:0.1.1--pyhdfd78af_0` bundles logan-walker 0.1.0, so container users get the old `cov` primary metric and no prune-span guard, and nothing in their outputs shows it. The recipe floor admits 0.1.0, neither the recipe test nor the Python wrapper checks the walker version, and Bioconda's autobump bot published version-only recipes (PK-01, PK-04). Small: recipe build 1, autobump off, a run-time version floor.
4. **On the LexicMap route, reference-matching contigs hide variants, and diversity measures coverage.** Where an accession's HSPs overlap, the HSP closer to the reference wins every base, so SNPs, stops and frameshifts on bubble contigs vanish and the accession is tiered as clean wild type (PL-01, high; a medium-sized fix, and the two candidate rules need a test on real HSPs). Haplotype counts, entropy and wild-type fraction treat untiled positions as differences, so they track assembly fragmentation (PL-02, small).
5. **The sieve silently drops lineages, and can hang.** The wrapper leaves the binary's admission gate armed, so after record 2,000 a divergent lineage is rejected with exit 0 and which lineages survive depends on file order (PM-01, small: pass `--warmup 2147483647` and check the `REJECT` count). The vendored binary loops forever once more than 65,536 divergent records are promoted (PM-02, small upstream patch; refuse such inputs until it ships). The sieve is opt-in.
6. **Cohorts change between identical runs.** `harvest` and `stratify` sort on score alone, so tied accessions keep network arrival order; two default harvests of the same collection gave 20-accession cohorts sharing 14 accessions. PH-01; trivial (break ties by accession).
7. **Shared scratch and `nohup` lose work.** The 0.1.1 scratch sweep judges liveness by PID, so a run in another container or on another node that shares the `--scratch-dir` deletes a live run's files, which is then recorded as `no_path` with exit 0 (WF-01, high; a medium-sized fix: lock-file liveness). The new signal handler overrides `nohup`, so a run dies with all its results when the SSH session drops; this is a 0.1.1 regression (WF-02, trivial).
8. **Partial input is reported as success.** `lexicmap -d` takes a Galaxy download cut short as the whole table, with status ok and exit 0 (PL-03). The walker exits 0 when one accession of many succeeds (WP-03), and `harvest` exits 0 when some shards were lost (PH-02). Small each.
9. **Two spike panel records carry the wrong lineage.** `Omicron_EG51` is a 2020 B.1.2 spike one base from Wuhan, so D614G samples are labelled EG.5.1 and real EG.5.1 samples are not (DD-01); `Omicron_BA5` is a BA.2 reporter clone (DD-02). Small: replace the records and add a lineage-marker test.
10. **tools-iuc PR #8800 can be neither merged as it is nor bumped by changing tokens.** A Reference CDS dataset named `reconstructed` is staged under the same name as the walker's output, so the output overwrites the user's dataset through the symlink (GX-01, trivial). At 0.1.1 three PR tests fail and the wrapper cannot select `matched` (GX-10, GX-04). Its memory model counts 43 MB per graph and silently lowers the size cap, so the cohort depends on the Galaxy destination (GX-08, GX-05). Trivial to small each.
11. **Three of the Galaxy plan's UDT rows cannot be built as written.** The bare-binary sieve returns 0 or 1 records with exit 0 (GX-13), the walker UDT cannot reuse the IUC memory block and, ported, cuts the cap to 497 MB on usegalaxy.org's default job (GX-14), and workflows that reference UDTs by id do not resolve (GX-15). Trivial to small, all in the plan.
12. **The test gates would not catch these, and running them as root deletes `/dev/null`.** The golden gate runs only on manual dispatch and with non-default flags, and no test reaches a realistic alignment window or a successful reconstruction (TC-01, TC-02, TW-02, TW-03). `cargo test` as root unlinks `/dev/null` (TW-01, high). Small to medium.

## Suggested order of work
This follows the release plan in [review/backlog-verification.md](review/backlog-verification.md#5-recommended-release-plan). Low and info findings not named here can go into whichever release touches their file.

**Hotfix now: Bioconda 0.1.1 build 1, no source change.**
- disassembler recipe: walker floor `>=0.1.1,<0.2` (PK-01), `mafft >=7.526` (PK-02, PK-11), run dependencies trimmed (PK-06), an `msa` smoke test; mark the `_0` image do-not-use.
- logan-walker recipe: build 1, keep the osx skip with the GraphAligner reason (PK-05). All three recipes: `extra: autobump: enable: false` (PK-04). Sieve recipe stays at `70022da` (PK-08). Docs: `docs/packaging.md` and `environment.yml` (PK-02, PK-09, PK-05) and the release checklist (PK-04).

**0.1.2: disassembler and logan-walker, fixes and additive flags, no new walker output files.**
- Release gates first: TC-01, TC-02, TW-02, TW-01, TW-03.
- Walker results: WV-01 (with WV-07, WV-11, WV-12), WV-02, WV-03, WV-04 to WV-10; WA-01 to WA-05; WG-02, WG-03, WG-04, WG-05, WG-07; WP-01, WP-02, WP-04.
- Walker robustness: WF-01, WF-02, WF-04, WF-03, WF-06, WG-01, WG-08, WP-03, PK-03. Memory, backlog (b) stages 0-1: WP-20, WG-13, WP-10 and WP-08 (first steps), WP-25, WF-08, WF-16. Cheap provenance: the cache identity check from WF-05.
- Python results: PK-01 (run-time walker floor); PL-01, PL-02, PL-03; PH-01; PM-01; PM-02 (interim limit); PC-01 to PC-08; PH-02, PH-04 to PH-08; GX-07; PL-04 to PL-09; PM-03, PM-04, PM-06, PM-07, PM-08, PM-17.
- Data: DD-01, DD-02, DD-05. Backlog (d), `msa --no-tree`, `--prep-only`, `--aligned`: PM-18, PC-22. Hygiene: PK-08 (submodule forward to `70022da`), TC-23, WP-23, TC-14. Recipes `logan-walker >=0.1.2,<0.2`. tools-iuc #8800 revised against 0.1.2: GX-01, GX-04, GX-08, GX-05, GX-09, GX-11, GX-12, GX-19, GX-21.

**0.2.0.**
- Backlog (a), provenance in full: WF-05, PK-12, PK-07, WV-20, PM-16, TC-27. Backlog (c), per-accession status table: WP-14 and the summary counts of WP-03.
- Backlog (b) stages 2-3: PK-13, WP-09, the later steps of WP-10. Recipe floor `logan-walker >=0.2,<0.3`; a sieve floor once upstream ships the PM-02 fix.
- Galaxy: status output and retry branch, `--max-memory-mb` in the wrapper, corrected UDT rows (GX-13, GX-14, GX-15, GX-26). Publish workflows only against a release that has the critical and high fixes on their path.

## The maintainers' backlog and Galaxy plan
- **"The 0.1.1 package itself is correct" is false.** WV-01 (critical) and 14 high findings stand at 0.1.1, including a 0.1.1 regression (WF-02). The review recommends shipping 0.1.2 with these fixes first rather than waiting for the backlog items, and moving backlog (d) and (b) stages 0-1 into it.
- **Item 0 (Bioconda) found the recipe problem but not the image problem.** The 0.1.1 image has mafft 7.526, which is fine, and logan-walker 0.1.0, which is not (PK-01). "A stock install has a broken `msa`" holds only with bioconda ranked above conda-forge under flexible priority; the documented order gets 7.526 (PK-11). Autobump can be switched off per recipe, so opening the PR first is not "the only defence" (PK-04). Dropping the osx skip cannot make macOS installable, because of GraphAligner 1.0.17b (PK-05).
- **(a) provenance: the gap is real, the scope is not ~10 lines.** Cache hits, `--offline` and Galaxy-staged inputs have no HTTP response, so capturing headers records nothing on the reproducible route. Record a content hash of the parsed bytes (Logan's ETag is the plain MD5) and a cache sidecar: 150 to 250 lines. The Logan manifests are 5.1 and 11.8 GB SRA tables that pin no unitig bytes (WF-05, PK-12).
- **(b) memory: three supporting claims are wrong, and the cheap saving is elsewhere.** At 0.1.1 defaults SRR8209082 keeps 34,582 of 81,534 unitigs (42%), not "~26k of ~1.28M" (WG-14); 1.7x is not a bound, the steady state is 3.17x (WP-08); one decode pass cannot reduce to the subgraph (PK-13). Of the 33.6 GB default bound, 20.0 GB is compressed buffers; stages 0-1 cut a measured peak from 530.8 to 50.7 MiB with no output change.
- **(c) status output: real, but not trivial, and not the largest silent drop.** The compressed size is not a field, `--fetch-only` writes nothing machine-readable, fetch-then-offline turns negatives into errors (WF-03), and the PR wrapper has no cap parameter (GX-05). A run that exits 0 when one accession succeeds hides more (WP-03).
- **(d) `msa` prep-only holds as written.** Add `--no-tree`, `--prep-only` and an `--aligned` check mode in 0.1.2 (PM-18). The Tool Shed already has MAFFT 7.526 with `--addfragments --keeplength` (`rnateam/mafft`), so no MAFFT container work is needed (GX-23).
- **(e) sieve: the ancestry test was read backwards, and the real fix is local.** `70022da` is upstream `main`; the submodule moved back by accident, so do not move the recipe to `39d90c0` (PK-08). Disarm the admission gate in the wrapper now (PM-01); the upstream request that matters is the hang above 65,536 records (PM-02). The interim UDT with an embedded dedup pre-pass cannot run: the sieve container has no Python (GX-13).
- **Galaxy plan: the direction is sound, the sketches do not run as drawn.** UDTs plus Format-2 workflows in iwc-lab, with PR #8800 as the Tool Shed route, is feasible. Three UDT rows cannot be built as written and the mafft row is unnecessary (GX-13, GX-14, GX-26, GX-23); the workflows need embedded tool definitions, a reference input and a real split step (GX-15, GX-27, GX-24); walker steps must use `logan-walker:0.1.1--hfa8f182_0`. First probe usegalaxy.org for Logan egress, element identifiers and a 4-core, 24 GiB job.
- **PR #8800** cannot be merged as it is (GX-01) or token-bumped (GX-10); three of the four open review comments are right on the facts. The bump checklist, the answers to the comments and the recommended route are in [review/galaxy-delivery.md](review/galaxy-delivery.md). Claim-by-claim verdicts, replacement text ready to paste, and designs for (a) to (d) are in [review/backlog-verification.md](review/backlog-verification.md).

## All findings
Severity follows a fixed rubric. **critical**: silently wrong scientific output in default usage, data destruction, or a remotely triggerable security hole. **high**: wrong output, crash or hang in plausible usage, a broken published install, or a reproducibility hole users cannot detect. **medium**: wrong output or crash on edge-case input, misleading documentation, resource blow-ups, significant test gaps. **low**: minor robustness, minor doc errors, cleanup. **info**: observations with no defect behind them.

| Area | critical | high | medium | low | info | total |
|---|---:|---:|---:|---:|---:|---:|
| [WV · Walker: majority correction, variants, VCF and output writers](review/walker-variants-output.md) | 1 | 2 | 7 | 14 | - | 24 |
| [WA · Walker: GraphAligner invocation and stitching](review/walker-align-stitch.md) | - | - | 8 | 6 | - | 14 |
| [WG · Walker: graph parsing, anchoring, subgraph extraction and phasing](review/walker-graph.md) | - | - | 7 | 11 | 1 | 19 |
| [WP · Walker: pipeline, command line, memory model and exit status](review/walker-pipeline-cli.md) | - | - | 10 | 17 | 1 | 28 |
| [WF · Walker: fetch, cache, inputs, scratch directory and signals](review/walker-fetch-scratch.md) | - | 2 | 9 | 11 | - | 22 |
| [TW · Walker test suite](review/tests-walker.md) | - | 1 | 4 | 8 | - | 13 |
| [PL · Python: LexicMap route](review/python-lexicmap.md) | - | 3 | 9 | 17 | - | 29 |
| [PM · Python: codon MSA and diversity sieve](review/python-msa-sieve.md) | - | 2 | 8 | 17 | - | 27 |
| [PH · Python: kmindex harvest and stratify](review/python-harvest-stratify.md) | - | 1 | 10 | 16 | - | 27 |
| [PC · Python: composites and plumbing](review/python-composites.md) | - | - | 9 | 13 | - | 22 |
| [TC · Tests, CI and bench](review/tests-ci-bench.md) | - | - | 7 | 19 | 1 | 27 |
| [DD · Documentation and reference data](review/docs-data.md) | - | 1 | 4 | 14 | - | 19 |
| [PK · Packaging, release and Bioconda](review/packaging-release.md) | - | 1 | 4 | 9 | 1 | 15 |
| [GX · Galaxy: tools-iuc PR #8800 wrappers, docs/galaxy.md and the plan](review/galaxy-wrappers.md) | - | 1 | 14 | 13 | - | 28 |
| **All areas** | **1** | **14** | **110** | **185** | **4** | **314** |

In the Where column, paths are relative to the repository root, except: `tools-iuc#8800/` is `tools/disassembler/` in PR #8800 at `b9a16ad`; `bioconda:` is the bioconda-recipes master branch; `GALAXY_PLAN.md`, `plan:` and `backlog:` are the maintainers' planning documents. Five pairs describe the same defect from two areas and should be fixed once: DD-03 and WV-13 (VCF docs call ALT the minor allele), GX-06 and WF-03 (fetch-then-offline negatives become errors), PK-03 and WF-10 (compiled-in TLS roots only), WP-07 and WF-09 (fixed 900 s fetch deadline), WP-12 and WV-17 (stale variant outputs in a reused `--out-dir`).

### [WV · Walker: majority correction, variants, VCF and output writers](review/walker-variants-output.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| WV-01 | critical | Default-on majority correction is discarded unless the alignment spans the whole stitched path | `walker/src/stages/variants.rs:259` | small |
| WV-02 | high | A walked allele at 1-10% of a fork is neither corrected nor reported | `walker/src/stages/variants.rs:188` | small |
| WV-03 | high | Logan `ka:f` is clipped at 25,000; VAFs and gates use it unflagged | `walker/src/stages/variants.rs:385` | small |
| WV-04 | medium | Fork VAF counts exits that never rejoin: deflated AF, lost alleles, blocked majority swap | `walker/src/stages/variants.rs:385` | small |
| WV-05 | medium | A refused prune reopens the variant gate: sub-floor branches published as PASS alleles | `walker/src/stages/variants.rs:553` | small |
| WV-06 | medium | Hard-coded cap of 32 alternatives silently truncates alt FASTA, VCF and snvs.tsv | `walker/src/stages/variants.rs:375` | small |
| WV-07 | medium | ">12 corrections, treat as provisional" is stderr-only; summary says `provisional: false` | `walker/src/stages/variants.rs:291` | small |
| WV-08 | medium | VCF and snvs.tsv drop every allele carrying two or more substitutions | `walker/src/io/output.rs:612` | small |
| WV-09 | medium | Stale 20 kb guard silently drops all SNV/VCF output for longer genes | `walker/src/io/output.rs:581` | small |
| WV-10 | medium | Banded alignment accepts a clipped in-band optimum; SNVs get wrong POS/REF | `walker/src/io/output.rs:726` | trivial |
| WV-11 | low | `--allow-frameshift-alts` also lets the majority swap reframe `_primary` | `walker/src/stages/variants.rs:230` | trivial |
| WV-12 | low | `prefer_majority` tries only the most abundant alternative at a fork | `walker/src/stages/variants.rs:192` | small |
| WV-13 | low | VCF docs call ALT/AF the minor allele; ALT can be the majority allele | `docs/vcf_example.md:28` | trivial |
| WV-14 | low | VCF text not validated: raw id as contig, IUPAC REF, per-sample values in INFO | `walker/src/io/output.rs:615` | small |
| WV-15 | low | Non-finite floats written as bare `inf`/`NaN`, making summary.json invalid | `walker/src/io/output.rs:223` | small |
| WV-16 | low | `--emit-vcf` deletes or overwrites any `<out-dir>/<ACC>.vcf`, walker-written or not | `walker/src/io/output.rs:475` | small |
| WV-17 | low | Reruns into one `--out-dir` keep stale variant outputs beside fresh primaries | `walker/src/main.rs:562` | small |
| WV-18 | low | Outputs written in place with the summary last; a late write error loses the summary | `walker/src/main.rs:558` | small |
| WV-19 | low | An over-long accession makes `--emit-vcf` abort before the summary is written | `walker/src/io/output.rs:476` | trivial |
| WV-20 | low | summary.json identifies references by id and length only; omits emit flags and variant decisions | `walker/src/io/output.rs:207` | small |
| WV-21 | low | Code docs say alt haplotypes are substitutions only and opt-in; indels are on by default | `walker/src/stages/variants.rs:41` | trivial |
| WV-22 | low | `find_alts` ignores the shared `view()`; oriented adjacency is built three times | `walker/src/stages/variants.rs:320` | small |
| WV-23 | low | variants.rs re-implements extract's median depth and relative floor | `walker/src/stages/variants.rs:538` | trivial |
| WV-24 | low | output.rs repeats writer setup and JSON arms; pident formula duplicated; header doc incomplete | `walker/src/io/output.rs:67` | small |

### [WA · Walker: GraphAligner invocation and stitching](review/walker-align-stitch.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| WA-01 | medium | A piped `-q` is re-read for every alignment, so every accession becomes `no_path` | `walker/src/stages/align.rs:169` | small |
| WA-02 | medium | N and IUPAC reference bases count as matches: N-laden panel records win and report inflated pident | `walker/src/stages/stitch.rs:139` | small |
| WA-03 | medium | Multi-anchor primary choice ignores route depth; on a tie, panel order can pick a 5% strain | `walker/src/stages/stitch.rs:165` | small |
| WA-04 | medium | "A connected path makes chimeras impossible" is false for mixtures; `_primary` has no mixture flag | `walker/src/stages/stitch.rs:64` | small |
| WA-05 | medium | `--min-cov` gates on walked length, not reference span: indels flip accept and reject | `walker/src/stages/stitch.rs:198` | small |
| WA-06 | medium | No deadline on GraphAligner: one stuck alignment stalls the run and loses all finished results | `walker/src/stages/align.rs:201` | small |
| WA-07 | medium | Single-record query header reaches GraphAligner verbatim: a tab gives `no_path`, a non-UTF-8 byte fails every accession | `walker/src/stages/align.rs:169` | trivial |
| WA-08 | medium | Single-record query: whitespace inside sequence lines reaches GraphAligner, giving `no_path` or SIGABRT | `walker/src/stages/align.rs:169` | small |
| WA-09 | low | Prune pass accepts a realignment that only trims 1-3 end bases | `walker/src/pipeline.rs:800` | small |
| WA-10 | low | Exit 0 without a GAF is reported as the biological negative `no_path` | `walker/src/stages/align.rs:231` | trivial |
| WA-11 | low | Query bytes outside GraphAligner's alphabet pass validation, then abort every alignment | `walker/src/io/input.rs:180` | trivial |
| WA-12 | low | GAF paths are not checked against the graph; an inconsistent row is published as a primary | `walker/src/stages/stitch.rs:92` | small |
| WA-13 | low | `LOGAN_KEEP_SCRATCH` ignored on failure paths, `0`/empty enable it, stderr tail reversed | `walker/src/stages/align.rs:211` | small |
| WA-14 | low | Multi-anchor mode re-logs the full reference-name list for every accession | `walker/src/stages/align.rs:190` | trivial |

### [WG · Walker: graph parsing, anchoring, subgraph extraction and phasing](review/walker-graph.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| WG-01 | medium | Empty Logan assemblies (55-byte objects) are reported as `error`, defeating `--allow-empty` | `walker/src/stages/parse.rs:206` | small |
| WG-02 | medium | Parser ignores Logan's documented `km:f` abundance alias; depth 0 silently disables majority and variant logic | `walker/src/stages/parse.rs:144` | trivial |
| WG-03 | medium | `bubbles.json` is gated only by an absolute 5x floor on the unpruned subgraph; deep samples get noise bubbles | `walker/src/stages/phase.rs:184` | small |
| WG-04 | medium | `bubbles.json` repeats one allele set once per flank pair, inflating bubble counts | `walker/src/stages/phase.rs:129` | small |
| WG-05 | medium | Prune floor is `max(prune_frac x depth, --min-abund)`: the allele gate moves the primary, undocumented | `walker/src/stages/extract.rs:159` | small |
| WG-06 | medium | Decoded size is unbounded: a crafted small cache file defeats `--max-size-mb` and the printed memory bound | `walker/src/stages/parse.rs:92` | small |
| WG-07 | medium | Parser silently accepts unitig files with no `L:` links; the walk runs on isolated unitigs | `walker/src/stages/parse.rs:154` | trivial |
| WG-08 | low | Undocumented 4.29 Gbp parser ceiling: larger accessions fail as `error`, not `too_large` | `walker/src/stages/parse.rs:137` | trivial |
| WG-09 | low | Non-finite `ka:f` values accepted; NaN can panic the median sort, inf gives invalid JSON | `walker/src/stages/extract.rs:132` | trivial |
| WG-10 | low | A cached unitig file cut exactly at a zstd frame boundary parses as a smaller, valid graph | `walker/src/stages/parse.rs:209` | trivial |
| WG-11 | low | Foldbacks and self-loops pass as two-allele bubbles | `walker/src/stages/phase.rs:117` | trivial |
| WG-12 | low | Neighbourhood sequence copy and tmpfs GFA are unbounded and outside the memory model; "~38 KB" is wrong | `walker/src/stages/extract.rs:232` | small |
| WG-13 | low | glibc arena retention after the first accession inflates per-worker memory beyond the 1.7x calibration | `walker/src/stages/parse.rs:80` | trivial |
| WG-14 | low | Subgraph-size figures, the guard description and the "5-hop" warning are stale against the defaults | `walker/src/stages/extract.rs:6` | small |
| WG-15 | low | Parse error names a nonexistent fallback; BCALM2 misattribution; stale anchor.rs input list | `walker/src/stages/parse.rs:124` | trivial |
| WG-16 | low | Comments name `subgraph_fraction` as the guard, call multi-anchor unfinished, and give `pident` the wrong unit | `walker/src/types.rs:178` | trivial |
| WG-17 | low | `AccessionResult` stores the primary haplotype three times; writers read different copies | `walker/src/types.rs:462` | small |
| WG-18 | low | `parse::build` packs bases one at a time; a lookup table with 4-base packing halves parse time | `walker/src/stages/parse.rs:182` | small |
| WG-19 | info | Anchor scan pays an unpredictable branch per base; a bitmap prefilter is about 2.5x faster | `walker/src/stages/anchor.rs:284` | small |

### [WP · Walker: pipeline, command line, memory model and exit status](review/walker-pipeline-cli.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| WP-01 | medium | Global prune floor removes thin gene stretches; prune is refused and reference-pull errors stay | `walker/src/pipeline.rs:729` | small |
| WP-02 | medium | Multi-anchor mode keeps one primary even when the `-q` references are different genes | `walker/src/pipeline.rs:667` | small |
| WP-03 | medium | Exit 0 whenever one accession succeeds, however many errored | `walker/src/main.rs:405` | small |
| WP-04 | medium | An error in the prune pass throws away the good first-pass haplotype | `walker/src/pipeline.rs:792` | small |
| WP-05 | medium | Summary file name embeds the raw first record id; `/` nests it, `..` escapes `--out-dir` | `walker/src/main.rs:578` | trivial |
| WP-06 | medium | Outputs are written only after the last accession; a bad `-o`, a kill or a stdout error loses the run | `walker/src/main.rs:557` | medium |
| WP-07 | medium | Fixed 900 s whole-transfer fetch deadline, no flag, retried from byte 0 | `walker/src/stages/stream.rs:49` | small |
| WP-08 | medium | Printed "worst case" memory bound is an unenforced estimate that real runs exceed | `walker/src/config.rs:445` | medium |
| WP-09 | medium | Every per-accession result is kept in memory until the run ends | `walker/src/pipeline.rs:335` | medium |
| WP-10 | medium | Memory bound is only printed, never checked against available RAM, and wrong for `--fetch-only` and short lists | `walker/src/main.rs:497` | medium |
| WP-11 | low | A GraphAligner known to clip the last base (1.0.19+) only draws a warning | `walker/src/main.rs:532` | small |
| WP-12 | low | Earlier runs' `snvs.tsv`, `alt_haplotypes.fasta` and VCFs survive in a reused `--out-dir` | `walker/src/main.rs:562` | small |
| WP-13 | low | Compute-stage failures reach stdout only; stderr says "see the log above" | `walker/src/pipeline.rs:416` | trivial |
| WP-14 | low | A run where every accession was skipped as `too_large` is reported as "nothing was found" | `walker/src/error.rs:82` | small |
| WP-15 | low | Single-anchor `no_path` drops the rejected alignment; log and summary cannot say why | `walker/src/pipeline.rs:640` | small |
| WP-16 | low | Summary candidates below `--min-cov` report null `cov`/`length`/`nodes` | `walker/src/pipeline.rs:655` | small |
| WP-17 | low | `--max-size-mb 0` is accepted as a 0-byte cap that refuses every accession | `walker/src/main.rs:211` | trivial |
| WP-18 | low | A panic in a fetch worker aborts the whole run with exit 101 | `walker/src/pipeline.rs:438` | small |
| WP-19 | low | Default of 8 compute workers oversubscribes hosts with fewer cores | `walker/src/config.rs:471` | small |
| WP-20 | low | Hard-wired 4-slot fetch queue adds four payloads to peak memory | `walker/src/pipeline.rs:281` | trivial |
| WP-21 | low | Help text, README, comments and docs disagree with the code in about fifteen places | `walker/src/pipeline.rs:580` | small |
| WP-22 | low | `pipeline::run` doc says downstream stages get only the first reference | `walker/src/pipeline.rs:273` | trivial |
| WP-23 | low | Crate-wide `allow(dead_code)` hides nine dead or test-only items | `walker/src/main.rs:24` | small |
| WP-24 | low | Prune survivor predicate is written three times, kept in step by a comment | `walker/src/pipeline.rs:743` | trivial |
| WP-25 | info | Graph-peak ratio comments give 1.5x and 1.6x where the constant is 1.7 | `walker/src/pipeline.rs:32` | trivial |
| WP-26 | low | A closed stdout (SIGPIPE) skips scratch cleanup and leaves GraphAligner children running | `walker/src/main.rs:331` | small |
| WP-27 | low | `--min-cov`, `--min-vaf`, `--min-abund`, `--min-shared-kmers` accept NaN, inf and out-of-range values | `walker/src/main.rs:197` | trivial |
| WP-28 | low | The argument parser keeps only the last repeated `-a`/`-l` and takes the next option as a value | `walker/src/main.rs:179` | small |

### [WF · Walker: fetch, cache, inputs, scratch directory and signals](review/walker-fetch-scratch.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| WF-01 | high | Scratch sweep judges liveness by PID and deletes live runs' directories across PID namespaces or hosts | `walker/src/scratch.rs:151` | medium |
| WF-02 | high | Walker overrides an inherited SIG_IGN, so a `nohup`'d run dies on hangup (0.1.1 regression) | `walker/src/scratch.rs:73` | trivial |
| WF-03 | medium | Fetch then `--offline` turns `no_logan_data` / `too_large` accessions into errors (exit 1 or 3) | `walker/src/stages/stream.rs:110` | small |
| WF-04 | medium | Scratch sweep recursively deletes any user directory named `logan-walker.<integer>[.*]` | `walker/src/scratch.rs:139` | trivial |
| WF-05 | medium | summary.json records no input identity; cache files are used without any identity check | `walker/src/io/output.rs:245` | medium |
| WF-06 | medium | `--fetch-only` reports a failed cache write as "cached" and exits 0 | `walker/src/stages/stream.rs:128` | small |
| WF-07 | medium | 1.5 s fixed retry budget turns brief outages into per-accession errors; 408/429 never retried | `walker/src/stages/stream.rs:161` | small |
| WF-08 | medium | Compressed unitig files held whole in RAM in 4+F+W places; cache hits copied needlessly | `walker/src/stages/stream.rs:92` | small |
| WF-09 | medium | Fixed 900 s whole-request fetch timeout discards partial downloads and restarts from byte 0 | `walker/src/stages/stream.rs:49` | small |
| WF-10 | medium | Fetch trusts only compiled-in Mozilla roots; `SSL_CERT_FILE` ignored behind TLS-inspecting networks | `walker/src/stages/stream.rs:48` | small |
| WF-11 | medium | Query FASTA read twice: piped multi-record panel refused, duplicate warnings, two parsers | `walker/src/main.rs:460` | small |
| WF-12 | low | Lowercase SRA accession is fetched verbatim, 404s, and is recorded as `no_logan_data` | `walker/src/io/input.rs:107` | trivial |
| WF-13 | low | Unitigs cache accepts and trusts non-zstd entries; a bad entry fails forever, even online | `walker/src/stages/stream.rs:96` | small |
| WF-14 | low | `LOGAN_KEEP_SCRATCH` directory is deleted by the next run's sweep | `walker/src/scratch.rs:151` | trivial |
| WF-15 | low | Signal handler removes scratch and exits without stopping running GraphAligner children | `walker/src/scratch.rs:80` | small |
| WF-16 | low | `--fetch-only` reads every cached file whole, validates nothing, counts hits as fetched | `walker/src/stages/stream.rs:92` | small |
| WF-17 | low | Scratch directory left behind when the run dies by SIGPIPE, SIGQUIT, SIGUSR1/2 or SIGXCPU | `walker/src/scratch.rs:73` | small |
| WF-18 | low | Runtime failures labelled "usage error"; non-zstd payloads "parse failed"; non-UTF-8 GAF fails every accession | `walker/src/error.rs:75` | small |
| WF-19 | low | Killed cache writer leaves orphan hidden `.tmp` files; entry not fsynced before rename | `walker/src/stages/stream.rs:262` | small |
| WF-20 | low | As PID 1 of a container, SIGTERM/SIGINT ends the walker with SIGSEGV (exit 139) | `walker/src/scratch.rs:80` | trivial |
| WF-21 | low | Unitig fetcher follows redirects to any origin, including plain http | `walker/src/stages/stream.rs:48` | trivial |
| WF-22 | low | Inputs read with no size limit (gzip inflated whole); accession length unbounded | `walker/src/io/input.rs:43` | small |

### [TW · Walker test suite](review/tests-walker.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| TW-01 | high | `Subgraph`'s Drop unlinks paths that test fixtures never created; as root `cargo test` deletes /dev/null | `walker/src/stages/phase.rs:444` | small |
| TW-02 | medium | No test covers consensus correction, alt recovery or band widening with realistic alignment windows | `walker/src/stages/variants.rs:735` | medium |
| TW-03 | medium | No automated test reaches a successful reconstruction or `--fetch-only` | `walker/tests/cli.rs:37` | small |
| TW-04 | medium | No test checks that each walker flag sets the field it names, or that unknown flags are refused | `walker/src/main.rs:695` | small |
| TW-05 | medium | The input-order test sorts its own vector; removing `run`'s sort passes every test | `walker/src/pipeline.rs:961` | small |
| TW-06 | low | Phase and extract tests do not pin strands, ordering, the fork filter, the prune-floor tie or the edge guards | `walker/src/stages/phase.rs:450` | small |
| TW-07 | low | `parse_gaf` is dead code; seven parser tests pin it instead of the live `parse_gaf_multi` | `walker/src/stages/align.rs:326` | small |
| TW-08 | low | `LOGAN_KEEP_SCRATCH` in the environment fails four walker tests | `walker/tests/cli.rs:67` | trivial |
| TW-09 | low | `size_cap_tests` fail when a proxy variable is set without `NO_PROXY` covering 127.0.0.1 | `walker/src/stages/stream.rs:413` | trivial |
| TW-10 | low | Walker tests leak empty temp directories and use the real `/dev/shm` scratch base | `walker/src/stages/extract.rs:742` | small |
| TW-11 | low | Two tests in the data-backed CI step pass when their inputs are broken or absent | `walker/src/io/input.rs:601` | trivial |
| TW-12 | low | The recall-study subgraph test passes on a corrupt fixture and reports it as "no data" | `walker/src/stages/extract.rs:566` | trivial |
| TW-13 | low | `emits_gfa_for_cross_check` deletes the GFA it reports having written | `walker/src/stages/extract.rs:884` | trivial |

### [PL · Python: LexicMap route](review/python-lexicmap.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| PL-01 | high | Overlapping reference-matching HSPs silently override variants, stops and frameshifts | `src/disassembler/cmd/lexicmap.py:363` | medium |
| PL-02 | high | Diversity summary counts coverage differences as distinct haplotypes | `src/disassembler/cmd/lexicmap.py:761` | small |
| PL-03 | high | `-d` reads a Galaxy download cut short of Content-Length as the whole table | `src/disassembler/cmd/lexicmap.py:811` | small |
| PL-04 | medium | Coding status from length alone; genetic code 11 only; reference never checked for stops | `src/disassembler/cmd/lexicmap.py:537` | small |
| PL-05 | medium | Last codon is never stop-scanned, even when the reference has no terminal stop | `src/disassembler/cmd/lexicmap.py:164` | trivial |
| PL-06 | medium | Indel inside an HSP overlap is marked an edge artefact, so a frameshift can pass as clean | `src/disassembler/cmd/lexicmap.py:410` | trivial |
| PL-07 | medium | Untrusted (edge) deletion's `-` displaces an overlapping HSP's real base | `src/disassembler/cmd/lexicmap.py:389` | small |
| PL-08 | medium | Whitespace inside reference sequence lines becomes reference positions | `src/disassembler/cmd/lexicmap.py:221` | trivial |
| PL-09 | medium | Frameshift tier applied to references classified as non-coding | `src/disassembler/cmd/lexicmap.py:749` | trivial |
| PL-10 | medium | Raw FASTA record id used as a file name: `/` ids crash, `../` or absolute ids escape the output directory | `src/disassembler/cmd/lexicmap.py:578` | small |
| PL-11 | medium | `clean.a2m` is a nucleotide alignment, but docs offer it to protein language models | `src/disassembler/cmd/lexicmap.py:56` | trivial |
| PL-12 | medium | Ungrouped input aborts mid-stream leaving misclassified files; the advised `sort` breaks the header | `src/disassembler/cmd/lexicmap.py:710` | small |
| PL-13 | low | `--max-records` cuts through an accession's HSP block and mis-tiers it | `src/disassembler/cmd/lexicmap.py:890` | trivial |
| PL-14 | low | Rows with too few columns are dropped silently | `src/disassembler/cmd/lexicmap.py:821` | small |
| PL-15 | low | Quantiles pick a rounded index; even-cohort medians are wrong and the regime call flips | `src/disassembler/cmd/lexicmap.py:178` | trivial |
| PL-16 | low | Empty or untiled cohort is called "divergent" with invented 0.0 medians | `src/disassembler/cmd/lexicmap.py:941` | small |
| PL-17 | low | Frameshift span counts deleted reference bases | `src/disassembler/cmd/lexicmap.py:452` | trivial |
| PL-18 | low | `--min-cov-partial 0` files zero-coverage accessions as clean all-gap rows | `src/disassembler/cmd/lexicmap.py:785` | trivial |
| PL-19 | low | Under `--require-full-length`, `--min-cov` still drives a check, labels and the rescue count | `src/disassembler/cmd/lexicmap.py:1284` | small |
| PL-20 | low | A row's `qseq` is never checked against the reference slice it is projected onto | `src/disassembler/cmd/lexicmap.py:868` | small |
| PL-21 | low | Missing `qseq`/`sseq` and truncated or corrupt gzip end in tracebacks | `src/disassembler/cmd/lexicmap.py:636` | small |
| PL-22 | low | Outputs written in place; failed or `--no-a2m` reruns leave new and stale files mixed | `src/disassembler/cmd/lexicmap.py:646` | small |
| PL-23 | low | Galaxy `deferred` and `paused` datasets exit 75 "try again later" | `src/disassembler/cmd/lexicmap.py:252` | trivial |
| PL-24 | low | Memory grows about 180 B per accession, contrary to the documented model | `src/disassembler/cmd/lexicmap.py:693` | small |
| PL-25 | low | Tiling builds a tuple per reference position; slice assignment is about 3.5x faster | `src/disassembler/cmd/lexicmap.py:343` | small |
| PL-26 | low | `stream_lexicmap_to_msa(file_path="-")` closes the process's `sys.stdin` | `src/disassembler/cmd/lexicmap.py:924` | trivial |
| PL-27 | low | Galaxy key lookup and HTTP access duplicated; tool-not-found messages and exit codes differ | `src/disassembler/cmd/lexicmap.py:128` | small |
| PL-28 | low | Test-only wrapper, a counter never read, triplicated flagged-row writes | `src/disassembler/cmd/lexicmap.py:426` | small |
| PL-29 | low | Haplotype table never reports the protein change of an in-frame indel | `src/disassembler/cmd/lexicmap.py:1068` | small |

### [PM · Python: codon MSA and diversity sieve](review/python-msa-sieve.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| PM-01 | high | Sieve wrapper leaves the binary's admission gate armed; divergent records after record 2,000 are dropped | `src/disassembler/cmd/sieve.py:214` | small |
| PM-02 | high | Sieve binary hangs forever once its fixed 2^20-slot LSH table fills (more than 65,536 divergent records) | `extern/kmer-diversity-sieve/scripts/l1_pipeline.c:155` | small |
| PM-03 | medium | msa silently hides frameshifts and trimmed insertions; docs promise codon-wise gaps it does not make | `src/disassembler/cmd/msa.py:282` | small |
| PM-04 | medium | Sieve judges novelty only within a month parsed from `\|YYYY-MM` headers; undocumented | `extern/kmer-diversity-sieve/scripts/l1_pipeline.c:217` | small |
| PM-05 | medium | sieve's shard glob deletes look-alike user files and misses shards for names with glob characters | `src/disassembler/cmd/sieve.py:182` | small |
| PM-06 | medium | msa temp names built from user basenames collide; MAFFT fails or a 2-row alignment replaces `--out-aln` | `src/disassembler/cmd/msa.py:229` | small |
| PM-07 | medium | sieve passes FASTA through unnormalised: a BOM drops record 1, whitespace counts as bad bases | `src/disassembler/cmd/sieve.py:88` | trivial |
| PM-08 | medium | Sieve sketches only the first 65,536 k-mers of a record, silently merging long loci | `extern/kmer-diversity-sieve/scripts/l1_pipeline.c:264` | trivial |
| PM-09 | medium | msa renames temp files over special-file outputs: `/dev/null` replaced as root, traceback otherwise | `src/disassembler/cmd/msa.py:326` | small |
| PM-10 | medium | msa aligns every exact-duplicate row with MAFFT; aligning distinct rows once gives identical output | `src/disassembler/cmd/msa.py:251` | small |
| PM-11 | low | Codon MSA is written in lower case; docs claim rows come through "byte for byte" | `src/disassembler/cmd/msa.py:279` | trivial |
| PM-12 | low | msa keeps a minus-strand row it cannot orient; it is aligned as noise with only a warning | `src/disassembler/cmd/msa.py:248` | small |
| PM-13 | low | Sieve output re-selected by 255-character header: same-header records over-emitted, long headers truncated | `extern/kmer-diversity-sieve/scripts/l1_pipeline.c:226` | small |
| PM-14 | low | A sieve run that promotes nothing truncates the previous output and blames `--threshold` | `src/disassembler/cmd/sieve.py:224` | small |
| PM-15 | low | Sieve dedup records no member ids and counts empty records as exact duplicates | `src/disassembler/cmd/sieve.py:199` | small |
| PM-16 | low | Python stages record no tool versions; msa discards the MAFFT and FastTree version output | `src/disassembler/cmd/msa.py:289` | small |
| PM-17 | low | msa runs whatever MAFFT is on PATH; no guard against the known-bad 7.525 | `src/disassembler/cmd/msa.py:186` | small |
| PM-18 | low | msa cannot write an alignment without FastTree; composite `--skip-tree` help misstates the sieve | `src/disassembler/cmd/msa.py:160` | medium |
| PM-19 | low | Sieve binary killed by a signal is reported as "exited -4" with exit status 252 | `src/disassembler/cmd/sieve.py:219` | trivial |
| PM-20 | low | sieve reports an unwritable output directory as "cannot read the input FASTA" | `src/disassembler/cmd/sieve.py:194` | trivial |
| PM-21 | low | msa and other subcommands crash on a non-UTF-8 stdout because of the U+2713 marker | `src/disassembler/cmd/msa.py:297` | trivial |
| PM-22 | low | msa `sanitize_id` is quadratic when many ids share their first 40 characters | `src/disassembler/cmd/msa.py:57` | trivial |
| PM-23 | low | msa holds cohort and alignment in memory and checks rows with per-character Python loops | `src/disassembler/cmd/msa.py:197` | small |
| PM-24 | low | Sieve dedup holds every distinct sequence in memory; the wrapper re-reads files it could stream | `src/disassembler/cmd/sieve.py:71` | small |
| PM-25 | low | msa coverage gate counts raw bases before alignment; flanked partial-gene rows enter MSA and tree | `src/disassembler/cmd/msa.py:242` | small |
| PM-26 | low | msa truncates record ids to 40 characters and renames clean ids order-dependently, with no name map | `src/disassembler/cmd/msa.py:54` | small |
| PM-27 | low | msa and sieve scratch files use predictable PID names and follow planted symlinks | `src/disassembler/cmd/msa.py:229` | small |

### [PH · Python: kmindex harvest and stratify](review/python-harvest-stratify.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| PH-01 | high | Tied scores keep shard arrival order, so the same collection yields different lists and cohorts | `src/disassembler/cmd/harvest.py:180` | trivial |
| PH-02 | medium | A partial harvest (some shards lost) exits 0 and summary.json does not record the loss | `src/disassembler/cmd/harvest.py:279` | small |
| PH-03 | medium | 0.1.0 harvest directories are unusable by 0.1.1 stratify; `from-kmindex --skip-harvest` dead-ends or relabels the wrong target | `src/disassembler/cmd/stratify.py:229` | small |
| PH-04 | medium | Wrongly shaped shard JSON counts as a clean shard (0 targets, exit 0) | `src/disassembler/cmd/harvest.py:151` | small |
| PH-05 | medium | Re-harvesting into the same directory keeps old targets, and stratify and from-kmindex use them | `src/disassembler/cmd/harvest.py:179` | small |
| PH-06 | medium | stratify's directory scan silently drops damaged score files; a corrupt gzip escapes as a traceback | `src/disassembler/cmd/stratify.py:196` | small |
| PH-07 | medium | Targets are renamed by a phiX174 substring strip, and name collisions are merged silently | `src/disassembler/cmd/harvest.py:75` | small |
| PH-08 | medium | A Galaxy collection with zero elements harvests to an empty cohort with exit 0 | `src/disassembler/cmd/harvest.py:385` | trivial |
| PH-09 | medium | Docs steer the Galaxy API key into world-readable /tmp/gxy.txt; `--api-key` exposes it on argv | `src/disassembler/cmd/harvest.py:47` | small |
| PH-10 | medium | urllib3 retries 429/503 with Retry-After inside each attempt: 4x requests, 120 s cap bypassed | `src/disassembler/cmd/harvest.py:324` | trivial |
| PH-11 | medium | stratify parses and keeps every target's score map in memory, even with `--target` | `src/disassembler/cmd/stratify.py:195` | small |
| PH-12 | low | A UTF-8 BOM makes a score file unreadable, or silently drops it from a directory | `src/disassembler/textio.py:32` | trivial |
| PH-13 | low | stratified_sample is top-anchored: N=1 takes the bin maximum, the bin bottom is never drawn | `src/disassembler/cmd/stratify.py:58` | trivial |
| PH-14 | low | The x-api-key header is re-sent on cross-host redirects; an http:// base URL is accepted silently | `src/disassembler/cmd/harvest.py:328` | small |
| PH-15 | low | A Galaxy API key with an embedded newline is echoed in full in the error message | `src/disassembler/cmd/harvest.py:44` | small |
| PH-16 | low | stratify accepts out-of-range and non-finite bins and exits 0 having selected nothing | `src/disassembler/cmd/stratify.py:76` | small |
| PH-17 | low | Non-JSON 200 shards are retried 5 times, then the API key is blamed | `src/disassembler/cmd/harvest.py:113` | small |
| PH-18 | low | A trailing-slash `--base-url` or a non-JSON 200 listing ends in a traceback | `src/disassembler/cmd/harvest.py:303` | trivial |
| PH-19 | low | Shard names, accessions and scores are trusted: a NUL or overlong target name aborts with a traceback | `src/disassembler/cmd/harvest.py:76` | small |
| PH-20 | low | Following stratify's own advice fails: a directory plus NAME=FILE for a file inside it conflict | `src/disassembler/cmd/stratify.py:238` | trivial |
| PH-21 | low | Per-target stratify runs overwrite the "combined" bins_summary.json with the last target only | `src/disassembler/cmd/stratify.py:293` | trivial |
| PH-22 | low | Ctrl-C during a Galaxy harvest keeps fetching every queued shard before exiting | `src/disassembler/cmd/harvest.py:341` | small |
| PH-23 | low | A whitespace-only $GALAXY_API_KEY or blank /tmp/gxy.txt is sent as an empty key | `src/disassembler/cmd/harvest.py:45` | small |
| PH-24 | low | Galaxy-mode harvest does not name the shard it discards; the per-shard loop is duplicated | `src/disassembler/cmd/harvest.py:360` | small |
| PH-25 | low | `--scores` rejects pipes as "no such file"; `--shards-dir FILE` says "not found" | `src/disassembler/cmd/stratify.py:202` | trivial |
| PH-26 | low | stratified_sample's docstring, its test and cli.py describe removed code | `src/disassembler/cmd/stratify.py:46` | trivial |
| PH-27 | low | harvest keeps every parsed shard alive until the end and builds each output file as one string | `src/disassembler/cmd/harvest.py:342` | small |

### [PC · Python: composites and plumbing](review/python-composites.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| PC-01 | medium | from-lexicmap re-aligns reference-grid rows with MAFFT, moving bases off their LexicMap coordinates | `src/disassembler/cmd/pipe_from_lexicmap.py:372` | small |
| PC-02 | medium | `reconstruct --sieve` drops walker primaries at 85-95% of reference length | `src/disassembler/cmd/pipe_reconstruct.py:414` | small |
| PC-03 | medium | `from-lexicmap --sieve` drops clean rows whose in-frame deletion exceeds the 5% margin | `src/disassembler/cmd/pipe_from_lexicmap.py:339` | small |
| PC-04 | medium | Composites align a stale handoff FASTA from an earlier run when this run did not write it | `src/disassembler/cmd/pipe_reconstruct.py:392` | small |
| PC-05 | medium | from-kmindex reconstructs an earlier harvest's cohort when the new harvest lacks `--target` | `src/disassembler/cmd/pipe_from_kmindex.py:342` | small |
| PC-06 | medium | Unusable `--walker-bin`/`--sieve-bin` or `DISASSEMBLER_*_BIN` silently falls back to another binary | `src/disassembler/resolve.py:47` | small |
| PC-07 | medium | Relative `--walker-bin`/`--sieve-bin` (`./name`) is re-resolved through PATH by subprocess | `src/disassembler/resolve.py:37` | trivial |
| PC-08 | medium | Composites cannot set msa's `--min-identity`; `from-lexicmap --min-pident` below 60 fails at stage 3 | `src/disassembler/cmd/pipe_from_lexicmap.py:363` | small |
| PC-09 | medium | Child-tool wrappers do not forward SIGTERM/SIGHUP, orphaning walker/MAFFT and misreporting signal exits | `src/disassembler/cmd/walker.py:71` | small |
| PC-10 | low | reconstruct/from-kmindex turn walker exit 3 into 1 and have no working `--allow-empty` | `src/disassembler/cmd/pipe_reconstruct.py:390` | small |
| PC-11 | low | Stale-output cleanup keyed on this run's reference stem misses an earlier run's MSA and tree | `src/disassembler/cmd/pipe_reconstruct.py:432` | small |
| PC-12 | low | `-d` and `-f` are documented as mutually exclusive but both are accepted; Galaxy wins | `src/disassembler/cmd/pipe_from_lexicmap.py:206` | trivial |
| PC-13 | low | from-kmindex and from-lexicmap reject bad inputs and missing tools only after stage 1 | `src/disassembler/cmd/pipe_from_kmindex.py:258` | small |
| PC-14 | low | A stray positional token sends it and all later options to the walker passthrough | `src/disassembler/cmd/pipe_reconstruct.py:198` | small |
| PC-15 | low | Passthrough and validator ranges disagree with the walker (`-t`, `--prune-frac`, `--min-abund`, repeated flags) | `src/disassembler/cmd/pipe_reconstruct.py:213` | small |
| PC-16 | low | Empty `-a ""` crashes reconstruct and fetch with a TypeError traceback | `src/disassembler/cmd/pipe_reconstruct.py:361` | trivial |
| PC-17 | low | textio's gzip reader never closes the file it opened; truncated gzip ends in a traceback | `src/disassembler/textio.py:31` | small |
| PC-18 | low | Inputs are read whole with no decompressed-size cap; small gzip inputs can cost GBs of RAM | `src/disassembler/textio.py:48` | small |
| PC-19 | low | Sieve/MSA defaults copied as literals into composites and inert-knob checks | `src/disassembler/cmd/pipe_reconstruct.py:251` | small |
| PC-20 | low | Composites re-declare about 40 stage options by hand; help text has already drifted | `src/disassembler/cmd/pipe_from_kmindex.py:235` | medium |
| PC-21 | low | `cli.main` has no top-level error handling: OS errors, truncated gzip, a closed stdout and Ctrl-C end in tracebacks | `src/disassembler/cli.py:75` | small |
| PC-22 | low | `--skip-tree` help says the sieve still runs; `--sieve --skip-tree` still requires a sieve binary | `src/disassembler/cmd/pipe_from_kmindex.py:236` | trivial |

### [TC · Tests, CI and bench](review/tests-ci-bench.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| TC-01 | medium | Golden gate CI job runs only on manual dispatch; no CI run for either release | `.github/workflows/ci.yml:101` | small |
| TC-02 | medium | No end-to-end regression check of the default walk (pruning, majority preference, alts, VCF) | `bench/check_golden.sh:59` | medium |
| TC-03 | medium | Resolver, validators, `fetch` and walker passthrough are untested; the gap hides a silent binary swap and three wrapper defects | `src/disassembler/resolve.py:47` | small |
| TC-04 | medium | The k=21 simulation misreports its own results and measures anchoring only | `bench/simulated_kmer_sensitivity/README.md:25` | small |
| TC-05 | medium | Forwarding of 12 walker options from the composites is untested | `src/disassembler/cmd/pipe_reconstruct.py:357` | small |
| TC-06 | medium | Sieve gene-scale defaults (tau 0.005, sketch 256) are pinned by no test | `tests/test_sieve.sh:38` | small |
| TC-07 | medium | Even-spread sampling, max-merge, zero-drop, 429 handling and nested collections are unpinned | `tests/test_stratify.py:28` | small |
| TC-08 | low | `main` is unprotected: CI is advisory and fork PRs were merged untested | `.github/workflows/ci.yml:5` | trivial |
| TC-09 | low | CI builds without `--locked` and checks no release identifiers; Bioconda builds with `--locked` | `.github/workflows/ci.yml:34` | small |
| TC-10 | low | Golden gate green-lights a stale prebuilt walker without rebuilding or naming it | `bench/check_golden.sh:32` | trivial |
| TC-11 | low | Real-tool MSA test passes with `--keeplength` removed | `tests/test_msa_real_tools.py:75` | trivial |
| TC-12 | low | Reference-quality test skips `queries/` and FASTA extensions other than `.fa`/`.fasta` | `tests/test_panels.py:40` | small |
| TC-13 | low | Sieve tests skip green when the submodule is absent or the build fails | `tests/test_sieve.sh:51` | trivial |
| TC-14 | low | Sieve test bootstrap trusts any existing binary and calls `make` unguarded | `tests/test_sieve_integration.py:35` | small |
| TC-15 | low | No test pins `sieve_run`'s binary flags; the equivalence test hard-codes its own copy | `tests/test_dedup_sieve_equivalence.py:62` | small |
| TC-16 | low | Conserved-motif model protects 30-nt windows, shorter than a 31-mer | `bench/simulated_kmer_sensitivity/harness.py:143` | small |
| TC-17 | low | Simulation harness cannot regenerate its committed results | `bench/simulated_kmer_sensitivity/harness.py:488` | small |
| TC-18 | low | Golden gate labels a walk that found nothing as "NOT a walk regression" | `bench/check_golden.sh:85` | trivial |
| TC-19 | low | Golden gate ignores a failed `mktemp` and writes its work files into `/` | `bench/check_golden.sh:29` | trivial |
| TC-20 | low | CI lints with floating `stable` Rust under `-D warnings`; cache key ignores the compiler | `.github/workflows/ci.yml:23` | trivial |
| TC-21 | low | CI hygiene: no dependency audit, no job timeouts or permissions, no tag-time version check | `.github/workflows/ci.yml:3` | small |
| TC-22 | low | CI runs only on Linux and has no coverage signal, while macOS packages are planned | `.github/workflows/ci.yml:16` | small |
| TC-23 | low | No Python lint in CI; the existing debt is 13 auto-fixable lines | `.github/workflows/ci.yml:62` | trivial |
| TC-24 | low | Rust sources are not rustfmt-clean and CI has no format check | `.github/workflows/ci.yml:33` | small |
| TC-25 | low | `bench/golden/README.md` still documents the deleted archive, verifier and harness | `bench/golden/README.md:48` | small |
| TC-26 | low | Test scaffolding is duplicated across files; one fixture reader assumes single-line FASTA | `tests/test_composites.py:59` | small |
| TC-27 | info | Two GraphAligner 1.0.17b builds are in use and no record tells them apart | `bench/golden/GRAPHALIGNER_VERSION:1` | small |

### [DD · Documentation and reference data](review/docs-data.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| DD-01 | high | Spike panel record `Omicron_EG51` is a 2020 B.1.2 spike (Wuhan + D614G), not EG.5.1 | `panels/sarscov2_spike_11strains.fa:586` | small |
| DD-02 | medium | Spike panel record `Omicron_BA5` is a BA.2 reporter-clone spike, not BA.5 | `panels/sarscov2_spike_11strains.fa:456` | small |
| DD-03 | medium | VCF docs call ALT the minor allele; the walker also emits the sample's majority allele as ALT | `docs/vcf_example.md:28` | trivial |
| DD-04 | medium | The documented Zika worked example is flagged `provisional` by the 0.1.1 tangle guard, and no doc says so | `walker/src/pipeline.rs:577` | medium |
| DD-05 | medium | Master target file: `erg11` and `avr3a` var1..9 are unrelated genes from the reference's own contig | `queries/master_planetary_targets.fasta:61` | small |
| DD-06 | low | kmindex docs omit the required `--format json`; other formats harvest to an empty result with exit 0 | `docs/kmindex.md:36` | small |
| DD-07 | low | Quick start uses repository-relative `queries/` paths that a Bioconda install does not provide | `README.md:60` | trivial |
| DD-08 | low | README install line omits the linux-64-only and Python 3.11/3.12 limits of the Bioconda stack | `README.md:30` | trivial |
| DD-09 | low | Per-sample values written to VCF INFO, so the documented `bcftools merge` mixes INFO across samples | `walker/src/io/output.rs:493` | small |
| DD-10 | low | Seven master-file records contain N or IUPAC codes, two of them 32-42% N; no test reads `queries/` | `queries/master_planetary_targets.fasta:181` | small |
| DD-11 | low | README promises pipe inputs; `reconstruct`, `from-lexicmap` and `lexicmap` reject a piped reference | `README.md:98` | small |
| DD-12 | low | `vcf_example.md` SQL reads an undefined table `snv`; output block and sample-group text do not match the data | `docs/vcf_example.md:83` | trivial |
| DD-13 | low | `--min-shared-kmers 0` runs as 1 but the summary records 0; docs imply an enforced 1-5 range | `walker/src/main.rs:208` | trivial |
| DD-14 | low | Sieve docs: wrong O(N^2) motive, misstated effective threshold, 9-vs-10 count, unreproducible tables | `docs/sieve.md:4` | small |
| DD-15 | low | `lexicmap.md` says directory-mode outputs use the reference file stem; the code uses the record id | `docs/lexicmap.md:122` | trivial |
| DD-16 | low | `walker.md` Reproducibility describes the removed 1%-of-graph warning and a wrong 10% share | `docs/walker.md:288` | trivial |
| DD-17 | low | `contigs_vs_unitigs.md`: stale 2,708 nt primaries and an understated 15-33x size ratio | `docs/contigs_vs_unitigs.md:87` | trivial |
| DD-18 | low | LexicMap links point to a nonexistent `IndexThePlanet/LexicMap` repository | `docs/lexicmap.md:3` | trivial |
| DD-19 | low | Smaller doc drift: retry count, walker crate list, Galaxy tag naming | `docs/walker.md:341` | trivial |

### [PK · Packaging, release and Bioconda](review/packaging-release.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| PK-01 | high | Biocontainer `disassembler:0.1.1--pyhdfd78af_0` ships logan-walker 0.1.0; nothing checks the walker version | `bioconda:recipes/disassembler/meta.yaml:34` | small |
| PK-02 | medium | Docs hide that the published 0.1.1 can resolve MAFFT 7.525; `environment.yml` gives a bioconda-first install command | `docs/packaging.md:40` | trivial |
| PK-03 | medium | logan-walker trusts only compiled-in roots, so every fetch fails behind TLS-inspecting proxies | `walker/src/stages/stream.rs:48` | small |
| PK-04 | medium | Release procedure lets Bioconda's autobump bot publish version-only recipes | `docs/packaging.md:63` | small |
| PK-05 | medium | packaging.md blames the macOS gap on a fixed CI bug; GraphAligner 1.0.17b makes macOS uninstallable | `docs/packaging.md:27` | trivial |
| PK-06 | low | Published 0.1.1 recipe still pulls zstandard, numpy and pandas | `bioconda:recipes/disassembler/meta.yaml:30` | trivial |
| PK-07 | low | Versions are bumped only at release, so source builds report the previous release | `docs/packaging.md:58` | small |
| PK-08 | low | Sieve submodule pin was reverted by a docs commit; backlog (e) would move the recipe backwards | `extern/kmer-diversity-sieve` (gitlink) | trivial |
| PK-09 | low | Docs attribute the MAFFT 7.525 bug to `--add`; `msa` uses `--addfragments`, which 7.526 fixed | `docs/packaging.md:39` | trivial |
| PK-10 | low | `environment.yml` still installs zstandard, which nothing imports | `environment.yml:31` | trivial |
| PK-11 | low | Backlog item 0's "stock install has a broken `msa`" holds only with bioconda ahead of conda-forge | `backlog:61` | trivial |
| PK-12 | low | Backlog 1(a) understates provenance work: cached and offline runs have no headers, manifests are not free | `walker/src/stages/stream.rs:102` | small |
| PK-13 | low | Plan 1(b): O(subgraph) memory needs two decodes and a spooled adjacency, not one pass | `plan:1(b)` / `walker/src/stages/extract.rs:93` | large |
| PK-14 | low | Backlog calls Logan objects unindexed zstd; they are seekable with a seek table and record index | `backlog:1(b)` (line 148) | trivial |
| PK-15 | info | Backlog's ~34 B/unitig CSR figure is right for one accession but varies with graph shape | `walker/src/types.rs:97` | trivial |

### [GX · Galaxy: tools-iuc PR #8800 wrappers, docs/galaxy.md and the plan](review/galaxy-wrappers.md)
| ID | Severity | Title | Where | Fix effort |
|---|---|---|---|---|
| GX-01 | high | Walker wrapper overwrites the Reference CDS dataset when it is named `reconstructed` | `tools-iuc#8800/disassembler_walker.xml:55` | trivial |
| GX-02 | medium | Galaxy LexicMap route drops partly covered rows: sieve and msa run with walker-tuned gates | `tools-iuc#8800/disassembler_sieve.xml:21` | small |
| GX-03 | medium | Sieve wrapper's default 95% length gate drops walker primaries admitted at 85-95% coverage | `tools-iuc#8800/disassembler_sieve.xml:36` | small |
| GX-04 | medium | Walker wrapper forces `--primary-metric cov` and cannot select 0.1.1's `matched` | `tools-iuc#8800/disassembler_walker.xml:81` | trivial |
| GX-05 | medium | Walker wrapper lowers `--max-size-mb` to fit job memory, so the cohort depends on the destination | `tools-iuc#8800/disassembler_walker.xml:34` | small |
| GX-06 | medium | Fetch-then-offline runs report no-data and over-cap accessions as errors | `walker/src/stages/stream.rs:110` | small |
| GX-07 | medium | `harvest` silently drops kmindex `json_vec` scores and reports an empty harvest as success | `src/disassembler/cmd/harvest.py:162` | small |
| GX-08 | medium | Walker wrapper's memory fit counts 43 MB per graph; chosen settings exceed the allocation | `tools-iuc#8800/disassembler_walker.xml:26` | trivial |
| GX-09 | medium | `DISASSEMBLER_WALKER_MAX_SIZE_MB` cannot be overridden from job_conf | `tools-iuc#8800/disassembler_walker.xml:100` | trivial |
| GX-10 | medium | The 0.1.1 bump is not a token change: three PR tests fail, and galaxy.md already documents 0.1.1 | `docs/packaging.md:66` | small |
| GX-11 | medium | Stratify wrapper names targets after datasets, which 0.1.1 no longer honours | `tools-iuc#8800/disassembler_stratify.xml:38` | small |
| GX-12 | medium | `@WALKER_VERSION@` token defeats `planemo autoupdate` | `tools-iuc#8800/macros.xml:3` | small |
| GX-13 | medium | Plan's interim sieve UDT (bare `l1_pipeline_serial`) returns 0 or 1 records with exit 0 | `GALAXY_PLAN.md:125` | trivial |
| GX-14 | medium | Plan's walker UDT cannot reuse the IUC memory block; ported, it cuts the cap to 497 MB | `GALAXY_PLAN.md:122` | small |
| GX-15 | medium | Plan's workflows reference UDTs by id, which Galaxy cannot resolve; import caveat reversed | `GALAXY_PLAN.md:138` | small |
| GX-16 | low | Walker wrapper loses `walker_summary.json` in failed jobs | `tools-iuc#8800/disassembler_walker.xml:91` | trivial |
| GX-17 | low | galaxy.md offers `--allow-empty` as an alternative to filtering failed jobs | `docs/galaxy.md:80` | trivial |
| GX-18 | low | With `--require-full-length`, log and wrapper help still say the partial tier is 50-80% | `src/disassembler/cmd/lexicmap.py:955` | small |
| GX-19 | low | msa wrapper stale at 0.1.1: `--msa-min-cov` mislabelled, no `--min-identity` | `tools-iuc#8800/disassembler_msa.xml:30` | small |
| GX-20 | low | `--allow-frameshift-alts` help describes 0.1.0 semantics | `tools-iuc#8800/disassembler_lexicmap.xml:41` | trivial |
| GX-21 | low | msa/lexicmap help "gap-free, only adds the tree" becomes false at 0.1.1 | `tools-iuc#8800/disassembler_msa.xml:90` | trivial |
| GX-22 | low | Stratify N help says "evenly across the score range"; sampling is by rank | `tools-iuc#8800/disassembler_stratify.xml:51` | trivial |
| GX-23 | low | MAFFT: galaxy.md calls it an IUC tool; the plan says no 7.526 container exists | `docs/galaxy.md:86` | trivial |
| GX-24 | low | No native Galaxy step splits an accession list; galaxy.md misnames tool owners | `docs/galaxy.md:79` | trivial |
| GX-25 | low | Plan caveats omit UDT destination limits (`docker_net: none`; usegalaxy.org 8 cores/29 GB, CVMFS only) | `GALAXY_PLAN.md:195` | trivial |
| GX-26 | low | Plan's fetch UDT runs `disassembler` in the logan-walker container, which has none | `GALAXY_PLAN.md:123` | trivial |
| GX-27 | low | WF1 sketch omits the reference CDS and target selection | `GALAXY_PLAN.md:144` | trivial |
| GX-28 | low | WF3-frozen quota estimate (median x 50) understates skewed cohorts | `GALAXY_PLAN.md:178` | trivial |

## Proposed patches
[`review/patches/`](review/patches/) holds 70 patches for this repository, named by the finding ids they address; none is on any branch. Each applies on its own to `25c27f2` with `git apply` (the PM-02 patch is rooted at the sieve repository and applies inside `extern/kmer-diversity-sieve`); patches that touch the same lines conflict with each other and need a hand merge. The [README](review/patches/README.md) lists, for each patch, the files it touches, what it was tested with, and which pairs conflict. Four more patches, for the Galaxy wrappers of tools-iuc PR #8800 (GX-01, GX-02, GX-04, GX-10), are in [`review/patches/tools-iuc-pr8800/`](review/patches/tools-iuc-pr8800/README.md) and apply to the PR head `b9a16ad`. All of them are proposals to review, not changes that were made.

## How the review was done and its limits
The walker was built, clippy and both test suites were run, and the external hosts were probed. Then 38 reviewers, each with one assignment, read the code and the external inputs: line-by-line readers for each group of walker and Python source files; cross-cutting reviewers for contracts between stages and tools, wrappers and resource guards, behaviour removed since 0.1.0, reuse, efficiency, root causes and security; and reviewers for tests and CI, documentation, reference data, scientific method, packaging, the Galaxy wrappers and plan, and the backlog. Fourteen more reviewers each received the issues already reported for one area and looked only for new ones. Every candidate went to an independent verifier whose job was to refute it, reproduce it where possible, re-rate it against the rubric and work out a fix; clusters rated high or critical got a second, sceptical verifier. One editor per area merged duplicates, re-checked every cited line against `25c27f2` and wrote the section.

Reviewers reported 492 candidates: 161 were folded into others as duplicates, 1 was refuted ([review/appendix-refuted.md](review/appendix-refuted.md)), 3 remain plausible but unconfirmed (PL-12, WG-07, TC-27), and 1 was not verified on its own because it reports the same defect as WV-01. Of the 314 findings, 284 were reproduced by running code and 27 were confirmed by reading. Baseline on the review machine (Ubuntu 24.04, Rust 1.97, Python 3.13): `cargo build --release` succeeds and `cargo clippy --all-targets --release -- -D warnings` is clean; `cargo test --release` passes 222 unit tests with 17 data-backed tests ignored, and the 10 CLI tests pass only while `/dev/null` is intact (TW-01); the Python suite passes 168 tests with 4 skipped; `cargo fmt --check` reports 332 diff hunks (TC-24).

Limits:
- GraphAligner was not available. The golden concordance gate (`bench/check_golden.sh`) and the walker's data-backed ignored tests were not run, and neither were the PR's walker tests 1 and 2. The alignment stage was driven by stand-in aligner scripts and by GAF rows recorded in the repository; findings about behaviour on real alignments rest on those.
- MAFFT 7.505 and FastTree 2.1.11 from the distribution were used instead of the pinned 7.526 and 2.2.0, and the MAFFT 7.525 crash was not reproduced directly. conda was not available: what a solver resolves was reasoned from channel repodata and checked against the contents of the published containers, not by running a solve.
- GitHub pages and the API of other repositories could not be read. The state of the Bioconda pull requests the backlog cites (#71066, #71067, #71068) was read only from a summary of their public pages, not from the API, so details such as who may edit them were not checked. Of PR #8800, the wrapper files and the second reviewer's open comments were read; the first reviewer's resolved threads and nine collapsed conversations were not, so the review cannot say whether they already raised any defect listed here.
- Live Logan S3 was reachable: object headers were read and a few small unitig files were fetched. No large-cohort run was made.
- The patches were checked to apply to `25c27f2`. Rust patches were not rebuilt when they were collected; the test results quoted for them are those reported by the verifier who wrote each one.
- Some fixes need a maintainer decision that the review could not settle: the minority-branch rule in WV-02, the overlap-resolution rule in PL-01 (to be tested on real LexicMap HSPs), and whether the sieve should keep month buckets (PM-04, PM-13).
