# The 0.2 backlog and the Galaxy plan, checked against 0.1.1

This section checks the maintainers' own planning documents against the code at `25c27f2` (tags `v0.1.1`, `logan-walker-v0.1.1`) and against the live services they cite. The documents are "disassembler 0.2 backlog" (2026-11-10) and section 1 of "Galaxy delivery plan for disassembler" (2026-11-04). Each claim was checked by reading the code, by running the release binary (offline, with real Logan objects and a stand-in GraphAligner), or by querying api.anaconda.org, quay.io, the Bioconda recipes and `logan-pub.s3.amazonaws.com`. Finding ids (WV-01, PK-01, ...) refer to entries in the other sections of this review; GX ids are in the Galaxy section.

## 1. Summary

- **Item 0 is right that the Bioconda 0.1.1 recipes shipped without the maintainer's changes, but the image check it asked for missed the real problem.** `biocontainers/disassembler:0.1.1--pyhdfd78af_0` contains mafft 7.526, which is fine, and **logan-walker 0.1.0**, which is not (PK-01, high). The backlog's headline, "a stock install has a broken `msa`", holds only when bioconda is ranked before conda-forge with flexible priority (PK-11). The bot can be switched off per recipe (PK-04), so "the only defence" is wrong. Dropping the osx skip cannot make macOS installable, because of GraphAligner (PK-05).
- **"The 0.1.1 package itself is correct" is false.** The default-on majority correction is inert in practice (WV-01, critical). The review has 14 high findings, including a 0.1.1 regression (WF-02) and live scratch directories deleted across PID namespaces (WF-01). These defects should set the next release; the backlog items should not.
- **(a) provenance: the gap is real, but the scope is wrong.** The "~10 lines" estimate is really 150 to 250 lines. Cache hits, `--offline` and Galaxy-staged inputs have no HTTP response, so capturing headers alone would record nothing on the reproducible route. Logan's ETag is the plain MD5 of the object, so a content hash both works everywhere and can be checked against S3. The manifests do not give provenance "for free": they are 5.1 and 11.8 GB SRA tables that pin no unitig bytes (WF-05, PK-12).
- **(b) memory: the backlog corrected the plan's 700 B/unitig, but three of its own claims are wrong.** The subgraph figures it calls "still accurate" are wrong: at the defaults SRR8209082 keeps 34,582 of 81,534 unitigs (42%). The 1.7x constant is not a bound (3.17x steady state, WP-08). One decode pass cannot reduce memory to the subgraph (PK-13). And the argument has moved: of the 33.6 GB default bound, 20.0 GB is compressed buffers, not graphs.
- **(b): the cheapest large saving does not touch the parser.** Stages 0 and 1 below (queue depth 0, allocator threshold, an accession-count-aware bound, compressed bytes read from disk instead of held in RAM) cut a measured peak from 530.8 to 50.7 MiB. They change no output and fit a patch release.
- **(c) status output: the gap is real, but it is not "trivial", and it is not the largest silent drop.** The size is not a field anywhere. `--fetch-only` writes nothing machine-readable. Fetch-then-offline turns negatives into errors (WF-03). The PR #8800 wrapper has no cap parameter (GX-05). A raised cap stops at 2^32 bases (WG-08). Separately, a run exits 0 as soon as one accession succeeds, however many errored (WP-03).
- **(d) msa prep-only holds exactly as written.** A tested design adds `--no-tree`, `--prep-only` and an `--aligned` check mode. The plan's "no MAFFT 7.526 container" problem is moot: the ToolShed tool `rnateam/mafft` already runs 7.526 with `--addfragments --keeplength`.
- **(e) sieve: the ancestry claim reads the test backwards.** The recipe's `70022da` is upstream `main`; the submodule moved back by accident (PK-08). The more important error is the diagnosis. The admission gate also drops unique divergent records after record 2,000, and the wrapper can disarm it today (PM-01, high), so this is not blocked on upstream. The upstream request that matters is the hang above 65,536 records (PM-02, high).
- **Sections 4 and 5:** most of the corrections to the plan are right. "Subgraph figures unchanged" is wrong, and the PR #8800 0.1.1 bump is not a token change (GX-10). The wrapper's 43 MB memory model should be fixed now (GX-08), not after (b).
- **Priority change:** ship the recipe build-1 hotfix with the walker floor and autobump off now. Then ship a 0.1.2 that carries WV-01 and the other critical, high and silent-result fixes, (b) stages 0 and 1, the cheap cache identity check, and (d). (a), (c) and (b) stages 2 and 3 make 0.2.0.

## 2. Claim-by-claim verdicts

### Framing and "Why 0.2 and not 0.1.2"

| claim | verdict | evidence | related |
|---|---|---|---|
| "Every item below was checked in code or against a live service" | partly | Several checks were misread: submodule ancestry, subgraph figures, the image's walker version, autobump | PK-08, PK-01, PK-04 |
| The plan was "audited against the repo on 2026-11-04 at 0.1.0" | partly | Plan 1(a) describes 0.1.1's `summary.json` (versions, references, 19 settings), but v0.1.0's `write_summary` wrote none of them (`git show v0.1.0:walker/src/io/output.rs`, lines 298-328) | — |
| Recipe has `run_exports: pin_subpackage(name, max_pin="x.x")`; the #71068 branch walker floor is `>=0.1.1,<0.2` | true | bioconda master recipe; the `as042:disassembler-0.1.1` branch recipe | — |
| (a) and (c) add `summary.json` fields or a file; (d) adds a CLI mode; all additive | true | No file in `src/`, `tests/` or `bench/` reads the walker summary | — |
| (b) "may retire the `too_large` class outright" | false | Only stage 3 below, plus a subgraph cap, makes `too_large` a transfer guard. Stages 0 and 1 change no output | PK-13, WG-08 |
| "If only (d) lands, 0.1.2 is defensible; if (a), (b) or (c) land, go to 0.2.0" | partly | The rule is sound for schema changes. But WV-01 and the high findings cannot wait for 0.2, and (b) stages 0-1 and (d) are output-neutral or additive | WV-01 |

### Item 0: verified bullets and todo points

| claim | verdict | evidence | related |
|---|---|---|---|
| Autobump PRs #71066/#71067 opened minutes after the tags and were merged 2026-11-05, ahead of maintainer PR #71068, which is open | true | Tags 2026-11-02 03:55 UTC; packages built 04:23 and 04:24 UTC; uploaded 2026-11-05 04:14 UTC (api.anaconda.org). PR states were read from the PR pages through a summarising fetch, not the API | PK-04 |
| The published package carries none of the recipe hygiene | true | The master recipe is the 0.1.0 recipe with a new version and sha256 | PK-04, PK-06 |
| disassembler 0.1.1 build 0 depends on the ten listed packages | true | api.anaconda.org and repodata list the same ten | PK-06 |
| `conda create --dry-run -c bioconda -c conda-forge disassembler=0.1.1` resolves mafft 7.525 | true | Reasoned from repodata; conda was not run. Bioconda first + flexible gives 7.525. Bioconda first + strict does not solve: bioconda's biopython stops at 1.70 (Python ≤3.6) | PK-11 |
| 7.525 "segfaults on the `--add` reference-anchoring path `msa` uses" | partly | `msa.py:284` passes `--addfragments`. The only 7.525→7.526 change is upstream b4cc76f0, "Fixed a bug in --addfragments" (`core/pairlocalalign.c`). Reading the source is consistent with a NULL dereference; the crash was not reproduced here | PK-09, PM-17 |
| Bioconda's newest mafft is 7.525; 7.526 is conda-forge only | true | bioconda: 7.525 (2024-03); conda-forge: 7.526 on linux-64, osx-64, osx-arm64, linux-aarch64 | — |
| "A stock install of the published package has a broken `msa`" | partly | The documented order `-c conda-forge -c bioconda` (`README.md:33`, `docs/packaging.md:20`) gets 7.526 under either priority, and both disassembler images contain mafft 7.526 h4bc722e_0. 7.525 results only from bioconda first + flexible, which `environment.yml:3` (`conda install -c bioconda disassembler`) does lead to | PK-11, PK-02 |
| logan-walker 0.1.0 and 0.1.1 published, build 0, linux-64 only; osx skip still on master | true | api.anaconda.org; master recipe | — |
| "so macOS still cannot install disassembler at all" | partly | True, but the stated cause (the skip) is incomplete. GraphAligner 1.0.17b's two osx-64 builds need protobuf 3.14.0 (Python ≤3.9); disassembler needs Python ≥3.11; there is no osx-arm64 GraphAligner. Removing the skip does not help | PK-05 |
| kmer-diversity-sieve 0.1.0 build 0 only (linux-64, linux-aarch64) | true | api.anaconda.org | — |
| Todo 1: bump build 0→1 for disassembler and logan-walker; sieve is already build 1 | true | The branch still has `number: 0`, `number: 0`, and sieve `number: 1` | PK-01 |
| Todo 2: "Everything else in it is still correct" (mafft ≥7.526, deps trimmed, floor ≥0.1.1, osx skips dropped, osx-arm64 on the sieve) | partly | The description is accurate. Dropping logan-walker's osx skip publishes an osx-64 walker with which disassembler cannot be installed. The branch also lacks the autobump opt-out | PK-05, PK-04 |
| Todo 3: add a recipe test on the MAFFT path; "this is what let the unpinned mafft ship unnoticed" | partly | The tests are as listed. Bioconda's build solve ranks conda-forge first and got 7.526 in both images, so an `msa` test would have passed too: only the pin protects users. The test that did let a defect ship is `disassembler walker --version \| grep -q "logan-walker"`, which accepts 0.1.0. The model test, `tests/test_msa_real_tools.py`, still passes with `--keeplength` removed | PK-01, PK-11, TC-11 |
| Todo 4: mafft inside `disassembler:0.1.1--pyhdfd78af_0` is unknown | true | Now answered. conda-meta: mafft 7.526 h4bc722e_0 (conda-forge) and **logan-walker 0.1.0** hfa8f182_0. The binary prints `logan-walker 0.1.0` and its help shows `[default: cov]`. The image was built 2026-11-02; walker 0.1.1 reached the channel 2026-11-05 | PK-01 |
| Todo 5: bot branches live in `bioconda/bioconda-recipes` | true | Heads `bioconda:bump/disassembler`, `bioconda:bump/logan_walker` | — |
| Todo 5: `maintainerCanModify` is false; a recipe maintainer has `pull` only | not verifiable from here | Needs the GitHub API | — |
| Todo 5: "The only defence is to open the maintainer PR before or immediately with the tag push" | **false** | `extra: autobump: enable: false` is honoured by bioconda-utils (`autobump.py:399-402`, `ExcludeDisabled`). None of the three recipes sets it | PK-04 |

### (a) Logan provenance

| claim | verdict | evidence | related |
|---|---|---|---|
| `summary.json` records walker version, GraphAligner version, references and exactly 19 settings (`write_run_info`) | true | `walker/src/io/output.rs:199`, `:216` (`let settings: [(&str, String); 19]`) | — |
| Per accession: acc, status, n_unitigs, n_anchors, subgraph_nodes, subgraph_fraction, provisional, haplotypes, bubbles, candidates | partly | `output.rs:245-365`. The key is `reconstructed`, not `haplotypes`; `candidates` appears only in multi-anchor mode; the list omits `prune`, `error` and `timing` | — |
| "Nothing about what was actually downloaded" | true | No size, URL, ETag, source or hash. `pipeline.rs:536` drops `fetched`, the only holder of the size, before any result exists. A live-S3 run and an `--offline` run of SRR8209082 give summaries identical apart from timings | WF-05 |
| `u/{acc}/` objects are unversioned and mutable | true | Bucket versioning is off (`GET /?versioning` empty, `VersionId null`). `c/` was rewritten in place in January 2025 (old data kept under `c1.0/`), and `u/` gained objects on 2025-12-14 | WF-05 |
| A rerun after a Logan regeneration changes results with no detectable difference in any output | true | As above; the cache is also trusted unchecked: ERR2028820's bytes filed as SRR8209082 are analysed as SRR8209082 with exit 0 (`stream.rs:67-68`, `parse.rs:111`) | WF-05 |
| `try_once` reads `Content-Length` and discards the other headers | true | `walker/src/stages/stream.rs:189` is the only header read in the file | — |
| ETag and Last-Modified are served (SRR8209082: `d2dd288b…`, 2024-03-31) | true | HEAD: `"d2dd288b238c3c88b46acdcd78860dce"`, `Sun, 31 Mar 2024 17:38:36 GMT`. The ETag is the plain MD5 of the content: 50 of 50 sampled objects (up to 1.2 GB) have plain 32-hex ETags, and `md5sum` matched for each object downloaded | PK-12 |
| Capture ETag + Last-Modified + fetch timestamp, "~10 lines" | **false** | The cache-hit branch (`stream.rs:74-109`) has no response, so `--unitigs-cache`, `--offline`, `disassembler fetch` and the plan's WF3-frozen would record nothing. `FetchedUnitigs` is dropped before results exist. Five `AccessionResult` constructors need the field. There is no timestamp formatter and no JSON reader in the crate. Realistic size: 150-250 lines with tests | PK-12, WF-05 |
| Last-Modified as provenance | partly | It is a bulk-copy date: 47 of 50 objects are dated 30/31 Mar 2024, the rest 2025-12-14 | PK-12 |
| Logan release-versions the project (`stats/` manifests v1.1, v1.2; `c1.0/`) | true | Both manifests exist; `c1.0/` exists | — |
| Recording the manifest "gives release-level provenance" | **false** | The manifests are 5,103,319,821 and 11,789,623,482 bytes of SRA run metadata (35 columns) with no object hash or Logan version. `c1.0/` holds contigs only; `u/` has no versioned prefix. A release name pins no unitig bytes | PK-12 |

### (b) Stream-decode

| claim | verdict | evidence | related |
|---|---|---|---|
| The whole compressed object is held in RAM (`read_to_vec`, cap + 1) | true | `stream.rs:205-209`; cache hits also use `std::fs::read` (`stream.rs:92`), so the cache path copies too | WF-08 |
| `parse.rs` builds a full CSR over every unitig before anchors are consulted | true | `pipeline.rs:534` (`parse::build`), `:542` (`anchor::find`) | — |
| `extract.rs` keeps "~300 to 26k nodes out of ~1.28M"; "still accurate" | **false** | At 0.1.1 defaults (`-k 21`, 3 shared k-mers, 5 hops), ERR2028820/pfcrt keeps 747 of 1,277,362 (0.06%) and SRR8209082/hiv1_pol keeps 34,582 of 81,534 (42.4%). The 26k figure (`extract.rs:6-8`) is SRR8209082 under the earlier BLAST anchoring; its whole graph is 81,534 unitigs. The old figures reproduce only at `-k 31` | WG-14 |
| The plan's ~700 B/unitig is the old Python walker's figure | true | `parse.rs:3-5` | — |
| "The Rust CSR is ~34 B/unitig" | partly | True for ERR2028820.  Arrays: 33.54 B/unitig on ERR2028820 and 38.17 on SRR8209082 (formula 12 + 4 × links/unitig + bases/4). By peak RSS: 35.7 B/unitig | PK-15 |
| The plan overstates by ~20x; already wrong at 0.1.0 | true | 18-21x. `parse.rs` changed only in comments since v0.1.0 | — |
| `GRAPH_PEAK_PER_COMPRESSED_BYTE = 1.7`, "measured 1.25x and 1.62x, plus the parser's doubling slack" | partly | `config.rs:445` is right. The measured ratios already include growth slack, and 1.7 is 1.62 rounded up. It was measured one accession per run; from the second accession on the steady state is 3.17x (102.7 MiB against 43.5 MiB) because glibc's mmap threshold rises | WP-08, WG-13 |
| `max_resident_bytes()` prints ~20 GB in buffers plus graphs, ~34 GB at the defaults | true | `config.rs:432`, `main.rs:497-511`. Printed: 20 buffers (4 queued + 8 fetching + 8 parsing) × 1024 MB = 20.0 GB, plus 8 graphs = 13.6 GB, total 33.6 GB | — |
| A capped 1 GB accession "now implies ~1.7 GB of graph peak, not 43 MB"; "the wrapper's memory ask got larger" | partly | That describes the model only. The walker's real memory did not change at 0.1.1; PR #8800's 43 MB was already too low at 0.1.0. Measured runs exceed the walker's own formula by 7 to 19% | GX-08, WP-08 |
| Payoffs: `too_large` largely disappears; the wrapper stops needing a worst-case estimate | **false** | Only a design that keeps nothing proportional to the graph in RAM (stage 3 below) delivers this, and then the subgraph needs a cap of its own (42% of the graph on SRR8209082). Until then a raised cap stops at 2^32 bases | PK-13, WG-08 |
| Ruled out: range-fetch, because "the object is a flat unindexed `.zst`" | partly | The conclusion holds: there is no content index, and anchors are known only after every sequence is scanned. The reason is wrong: objects are multi-frame (about 128 MiB decompressed per frame) with a seek-table frame and a frame-to-first-record index | PK-14 |
| Ruled out: shrinking the cap (1 GB compressed is several GB decompressed) | true | Decompression ratios measured: 4.28 (ERR2028820), 5.31 (SRR8209082) | — |
| `pipeline.rs` doc says ~1.6x where the constant is 1.7 | true | Incomplete: ~1.5x at `pipeline.rs:32`; ~1.6x at `:77`, `:505` and `docs/walker.md:237`; 1.7 at `docs/walker.md:249`; 1.24x-1.62x at `types.rs:94` | WP-25 |

Not in the backlog: the four-slot queue holds payloads for nothing (WP-20). Decoded size is unbounded, so a small crafted cache file defeats the cap (WG-06). Every result is retained until the run ends (WP-09). The bound is printed but never enforced, and it ignores `--fetch-only` and short lists (WP-10).

### (c) Per-accession status output

| claim | verdict | evidence | related |
|---|---|---|---|
| No `skipped.tsv`; outputs are the FASTA, `bubbles.json`, the summary, optional VCFs, `snvs.tsv`, `alt_haplotypes.fasta` | true | `main.rs:556-581` | — |
| `too_large` and `no_logan_data` are separate statuses at 0.1.1 (both `error` at 0.1.0) | true | Six statuses: `output.rs:21-30`, `types.rs:411`; v0.1.0 had four | — |
| `--allow-empty` no longer hides an all-errored run (exit 3, `check_run`) | true | `main.rs:401-415`; reproduced with and without the flag | — |
| `--allow-empty` plus every accession skipped still exits 0 with an empty FASTA | true | Two cached accessions, `--max-size-mb 1 --allow-empty`: both `too_large`, exit 0, 0-byte FASTA. Without the flag: exit 1, "none errored: nothing was found", although nothing was examined | WP-14 |
| A `skipped.tsv` (accession, outcome, size) is "trivial" | partly | The size is not a field, only text inside `error` (`stream.rs:191-194`), and it is unknown without a Content-Length. `--fetch-only` writes no per-accession file. Fetch-then-`--offline` reports `no_logan_data` and `too_large` as `error` | WF-03, WP-14 |
| It becomes "a natural retry branch, re-invoking with a raised `--max-size-mb`" | partly | PR #8800 has no cap parameter (only an admin variable that job_conf cannot override), so a workflow cannot raise it. At one slot the bound charges a one-accession retry 7.7 × cap, where 2.7 × cap is true. Above 2^32 unitig bases (2.1 to 2.4 GB compressed) parsing fails as `error`, exit 3: a synthetic 4.3 Gbp file failed after 12.9 s and 1,031 MiB | GX-05, GX-09, WP-10, WG-08 |

Not in the backlog: a run exits 0 as soon as one accession succeeds, whatever the error count, and the summary has no error count (WP-03). Empty 55-byte Logan objects are reported as `error` (WG-01).

### (d) `msa` prep-only

| claim | verdict | evidence | related |
|---|---|---|---|
| `--out-tree` is `required=True`; FastTree always runs | true | `msa.py:160`; FastTree is looked up before any work (`:187`) and run unconditionally (`:315-326`) | PM-18 |
| `--skip-tree` exists only in `reconstruct` and `from-kmindex` | true | `pipe_reconstruct.py:181`, `pipe_from_kmindex.py:235`; `from-lexicmap --skip-tree` exits 2 | — |
| There it means "stop after the walker", not "align without a tree" | true | `pipe_reconstruct.py:429-435`. Its help text is also wrong: it says the sieve still runs, but the sieve is skipped (`:408`) | PC-22 |
| `docs/galaxy.md` Open Questions still says "needs a preparation-only mode in the CLI first" | true | `docs/galaxy.md:85-89`. Line 88 is itself wrong: `--skip-tree` drops the MSA too, so the split does not "turn `--skip-tree` into simply not connecting FastTree" | — |
| A prep-only mode lets a workflow attach "the real IUC `fasttree`" | true | `iuc/fasttree` 2.1.10+galaxy3. Its model defaults to Protein, so a workflow must set Nucleotide and GTR. Its version is not the 2.2.0 the CLI pins | — |

### (e) Sieve upstream

| claim | verdict | evidence | related |
|---|---|---|---|
| The dedup pre-pass (`byte_exact_dedup`) is still in `cmd/sieve.py`, and its docstring explains why | true | `sieve.py:55-69` | — |
| Duplicates corrupt the warmup statistics, causing 26-273 spurious `REJECT (error)` drops, so upstream must make admission dedup-aware | partly | The duplicate mechanism is real, but the diagnosis is incomplete. With `--warmup 2000` (`sieve.py:214`) the gate also rejects unique records: a 0.67%-divergent accession after 2,100 unique haplotypes is dropped (kept when placed first), and a 5-member clade 2% away loses all 5. A rejected record is never added, so its clade follows it. The wrapper can disarm the gate today with `--warmup 2147483647` | PM-01 |
| A real git tag is still missing (`git ls-remote --tags` empty) | true | Checked 2026-11-10 | — |
| The recipe defines 0.1.0 as `70022da`; the submodule is at `39d90c0` | true | `git ls-tree HEAD extern/kmer-diversity-sieve` → `39d90c0` | PK-08 |
| "`70022da` is no longer an ancestor of upstream main"; "upstream rewrote that history"; "the recipe pins a commit off the upstream branch" | **false** | Upstream `main` is `70022da`; its parent is `39d90c0`. `git merge-base --is-ancestor 70022da origin/main` exits 0, meaning it is an ancestor (it is the tip); the backlog's test had the arguments reversed. The submodule went backwards: 00e3b97 advanced it to `70022da`, then the docs commit b1b6acf set it back to `39d90c0` | PK-08 |
| The two commits differ only in `README.md` (56 deletions) | partly | README only, but it is 56 insertions going from `39d90c0` to `70022da` | PK-08 |
| Options: recipe → `39d90c0` / submodule back / press for a tag | partly | Option 1 moves the package backwards. Option 2 is right if "back" means forward to `70022da`. Option 3 stands | PK-08 |
| Interim: the sieve UDT keeps the dedup pre-pass in its `shell_command` | **false** | The sieve container has no Python (busybox only). The bare binary with archive-scale defaults (`min_len` 27,000, threshold 0.15) keeps 0 of 33 records of the HIV-1 pol test cohort, with exit 0 | GX-13 |

Not in the backlog: the binary hangs forever once more than 65,536 divergent records are promoted, because its 2^20-slot table never grows (PM-02). It sketches only the first 65,536 k-mers of a record (PM-08). It compares records only within a month parsed from `|YYYY-MM` headers (PM-04).

### Section 4: corrections owed to the plan

| claim | verdict | evidence | related |
|---|---|---|---|
| 1(b): drop "~700 B/unitig"; replace with the CSR cost and the 1.7x model, "which makes the same argument more strongly" | partly | The first half is right. The second misleads: with real numbers the dominant term is the compressed buffers (20.0 of 33.6 GB), not the CSR, and removing them needs no parser change | WF-08 |
| 1(f): the fixes were never in #70419 (merged); they are in #71068 (open) | true | #70419 (merged 2026-10-26) added the 0.1.0 recipes. The master recipe still has bare `mafft`, `zstandard`, `numpy`, `pandas` and `logan-walker >=0.1.0` | PK-04 |
| Section 2: the 0.1.0 container pins are stale; the 0.1.1 tags exist | partly | Both tags exist (pushed 2026-11-05), but `disassembler:0.1.1--pyhdfd78af_0` carries walker 0.1.0, so pinning it is wrong | PK-01 |
| Section 2, mafft row: "run inside the disassembler container's env" is not safe until 0.4 is checked | partly | The solve is unpinned, but both disassembler images have mafft 7.526. The row is moot anyway: ToolShed `rnateam/mafft` ships `rbc_mafft_add` 7.526+galaxy2 on a mulled mafft 7.526 image | GX-23 |
| 1(c): 0.1.1 already split the statuses and made an all-errored run exit 3 | true | `output.rs:21-30`, `main.rs:401-415` | — |
| Unchanged and still correct: subgraph scale figures | **false** | See (b) row 3 | WG-14 |
| Unchanged and still correct: ~71 MB median compressed accession | partly | `stream.rs:7` says ~71 MB, `docs/galaxy.md:27` about 50 MB (gene-hit cohorts), `docs/kmindex.md:82` a 0.6 MB median for 45 random accessions. Median × 50 understates skewed cohorts | GX-28 |
| Unchanged and still correct: walker jobs need outbound network; `--offline`, `--unitigs-cache`, `--fetch-only` are the escape hatches | true | `main.rs:255-267`. Behind TLS-inspecting proxies every fetch fails (compiled-in roots only) | PK-03 |
| Unchanged and still correct: live-S3 tests stay behind `#[ignore]` | true | `.github/workflows/ci.yml:33-37`, `stream.rs:336-384` | TC-01 |

### Section 5: carried forward

| claim | verdict | evidence | related |
|---|---|---|---|
| Sieve tests ERROR rather than skip on a binary built for another CPU; cheap fix: skip on a signal death | partly | Reproduced with a SIGILL stand-in: `test_sieve_integration` gives 6 FAIL and 2 ERROR of 15, `test_dedup_sieve_equivalence` 2 ERROR of 3. A skip in CI tests nothing; rebuild first, and skip only outside CI | TC-14 |
| No Python linter in CI; "would need its own pass" | partly | True that there is none. The debt is 13 one-line fixes (unused imports, f-strings without placeholders). Clippy is also blinded by the crate-wide `#![allow(dead_code, unused_variables)]` (`walker/src/main.rs:24`) | TC-23, WP-23 |
| Tests missing: `fetch`, walker passthrough, `resolve`, `validators`, minus-strand HSPs | true | No test imports them. The gap hides a silent binary swap (PC-06), PATH re-resolution of `./name` (PC-07) and signal exits reported as 256-N (PC-09) | TC-03 |
| Prune tolerance not measured at walker20k scale | not verifiable from here | No GraphAligner and no cohort. `PRUNE_SPAN_TOLERANCE = 3` (`pipeline.rs:174`); boundary 3 is kept and 4 refused (unit tests) | WA-09 |
| GraphAligner end clipping not reported upstream | not verifiable from here | Consistent with `docs/issue_graphaligner_end_clipping.md:4`, `:274` | WP-11 |
| `harness.py` hardcodes `work/` paths | true | Worse than stated: `REPO_ROOT` resolves to the directory above the checkout (`harness.py:28-30`); a missing genome silently becomes a random 4.6 Mb sequence (`:488-491`) | TC-17 |
| PR #8800 still ships 0.1.0 wrappers; the 0.1.1 bump list is still due | true | `macros.xml` tokens are 0.1.0. The bump is not a token change: three PR tests fail at 0.1.1, `--primary-metric matched` cannot be selected, and a Reference CDS dataset named `reconstructed` is overwritten | GX-10, GX-04, GX-01 |
| Its memory block uses the flat 43 MB model, which (b) "supersedes twice over" | partly | It does use 43 MB (`disassembler_walker.xml:26-34`). The fix should not wait for (b): count 2.7 × cap per compute worker now. Its settings exceed the allocation by up to 118% | GX-08 |
| The "msa is a no-op on lexicmap output" text was scoped locally; check whether it went up | partly | It went up (PR head `disassembler_msa.xml:90-93`, `disassembler_lexicmap.xml:246`), and it is wrong again at 0.1.1: a trusted in-frame deletion counts as covered, so the rows are not gap-free | GX-21 |

### "Verified already fixed at 0.1.1"

| claim | verdict | evidence | related |
|---|---|---|---|
| lexicmap frameshift tiering | true | `lexicmap.py:778-804`; `TestCoverageTakesPrecedence`. Neighbouring defects remain: an indel inside an HSP overlap passes as an edge artefact; non-coding references get the frameshift tier | PL-06, PL-09 |
| Multi-anchor primary ranking defaults to `matched` | true | In the walker: `config.rs:456`, `stitch.rs:156`. But the published 0.1.1 image runs walker 0.1.0 (`cov`); PR #8800 forces `cov`; N/IUPAC reference bases count as matches; ties fall to panel order | PK-01, GX-04, WA-02, WA-03 |
| Multi-anchor VCF coordinates | true | `output.rs:452` (`reference_of`); unit-tested only. A separate banded-alignment defect can still give wrong POS/REF | WV-10 |
| `from-kmindex` harvest crash | true | `tests/test_composites.py:144` | — |
| Walker hang | true | `pipeline.rs:328`, `main.rs:316-328`. Residue: death by SIGPIPE leaves scratch behind. A different hang remains: GraphAligner has no deadline | WF-17, WA-06 |
| Walker exit status (`no_logan_data`, `too_large`, exit 3 on all-errored) | partly | The walker holds. The composites turn 3 into 1 and have no working `--allow-empty`; empty Logan objects are `error`; exit is 0 with one success however many errored | PC-10, WG-01, WP-03 |
| `msa` strand orientation and identity check | partly | A minus-strand row at ≥30% divergence often shares no 15-mer, stays unoriented and is kept with a warning. No composite can set `--min-identity` | PM-12, PC-08 |
| Scratch-directory handling (`scratch.rs`) | **false** | A walker in another PID namespace deleted a live run's directory (2 of 4 accessions lost, exit 0). `nohup` no longer protects the walker (0.1.1 regression). User directories named `logan-walker.<n>` are deleted | WF-01, WF-02, WF-04 |
| `provisional` guard replacement | true | `config.rs:299`, `extract.rs:620-631`. Its warning says "5-hop" whatever `--hops` is, and the documented Zika example is now flagged provisional without the docs saying so | WG-14, DD-04 |
| 3-base prune tolerance | true | `pipeline.rs:158-174`: boundary 3 kept, 4 refused | WA-09 |
| Three defective spike panel records | true | For those three. Two other records in the same panel are wrong: `Omicron_EG51` is a 2020 B.1.2 spike, `Omicron_BA5` a BA.2 reporter clone. `queries/` has 7 N/IUPAC records and unrelated "variants" | DD-01, DD-02, DD-05, DD-10 |
| "The 0.1.1 package itself is correct; only its recipe is not" | **false** | The default majority correction is discarded unless the GAF row spans the whole path (WV-01, critical), and 14 high findings stand at 0.1.1 | WV-01, WV-02, WV-03, WF-01, WF-02 |

### Galaxy plan section 1

| claim | verdict | evidence | related |
|---|---|---|---|
| Preamble: changes "ordered by effort"; (a) and (b) are the blockers | false | (a) is medium (150-250 lines); (b) is large and staged; (c) small; (d) medium. Silently missing accessions are fixed by (c) plus a wrapper cap parameter, not by (b) | GX-05 |
| (a) `summary.json` has versions, references, 19 settings, "nothing about what was actually downloaded" | true | At 0.1.1; see the (a) table; v0.1.0 recorded none of these | WF-05 |
| (a) every `u/` object serves ETag and Last-Modified; objects are unversioned and mutable | true | 50 of 50 sampled; versioning off | — |
| (a) capturing the two headers + timestamp "is a small change to `stream.rs`/`output.rs`" | **false** | Cache, offline and fetch-collection routes have no response; 150-250 lines | PK-12 |
| (a) recording the manifest "gives release-level provenance for free" | **false** | 5.1 and 11.8 GB SRA tables, multipart ETags, no object hashes | PK-12 |
| (b) the walker reads the whole compressed object, then parses every unitig into a full CSR | true | `stream.rs:205-209`, `pipeline.rs:534` | WF-08 |
| (b) "~700 B/unitig RAM" | **false** | 33.5-38.2 B/unitig of arrays; 700 B was the retired Python parser (`parse.rs:3-5`) | PK-15 |
| (b) "discarding ~99.9%", "~300–26k nodes out of ~1.28M" | partly | 99.94% on ERR2028820, but 57.6% on SRR8209082 (34,582 of 81,534 kept) | WG-14 |
| (b) decoding incrementally and reducing to the subgraph "during the pass" gives O(subgraph) | **false** | The stream is already decoded incrementally (`stream.rs:298-309`, `parse.rs:77`). One pass cannot know which unitigs are within 5 hops; see 4(b) | PK-13 |
| (b) `too_large` largely disappears; no worst-case estimate needed | **false** | Only at stage 3 with a subgraph cap | PK-13 |
| (b) range-fetch impossible (flat unindexed `.zst`); shrinking the cap does not help | partly | The conclusion holds; the objects carry a seek table | PK-14 |
| (b) "the single change that fixes the two most-quoted Galaxy objections at their root" | partly | Memory headroom: mostly stages 0-1. Silently omitted accessions: (c) and GX-05 | — |
| (c) `--allow-empty` + all-skipped gives a green job with an empty output; `skipped.tsv` is "trivial" | partly | The gap is real; "trivial" is not. See the (c) table | WP-14, WF-03 |
| (d) `--skip-tree` "exists only inside the composites"; atomic `msa` requires `--out-tree` | partly | Two of three composites, and there it skips the sieve and MSA too | PM-18 |
| (d) `docs/galaxy.md` already lists the split | true | `docs/galaxy.md:85-89` | — |
| (e) dedup protects the warmup statistics; ask upstream for dedup-aware statistics | partly | The gate also drops unique records; disarm it downstream | PM-01 |
| (e) the recipe labels untagged `70022da` as 0.1.0 | true | sieve recipe on master | PK-08 |
| (e) interim: embed the dedup pre-pass in the sieve UDT's `shell_command` | false | No Python in that image; archive-scale defaults return 0-1 records | GX-13 |
| (f) header "(bioconda #70419, merged)" | **false** | #70419 added the 0.1.0 recipes; the hygiene is in open #71068 | PK-04 |
| (f) `mafft >=7.526` floor: 7.525 segfaults on the `--add` path | partly | The option is `--addfragments` | PK-09 |
| (f) 7.525 is bioconda's newest; "a bioconda-priority solve resolves the documented-broken version" | partly | Flexible priority only; strict does not solve; the documented order gets 7.526 | PK-11 |
| (f) drop `zstandard`, `numpy`, `pandas` | true | Nothing in `src/` imports them | PK-06 |
| (f) "Recipe test doesn't exercise the MAFFT path — a bad solve passes CI" | partly | True that it doesn't. Bioconda CI solves conda-forge-first, so the test would pass; the missing check was the walker version | PK-01, PK-11 |

## 3. Corrections to make

Replacement wording, ready to paste. "B" items are in the backlog, "P" items in the plan.

**B1. Front matter `description`.** Replace "because the published disassembler 0.1.1 resolves the segfaulting mafft 7.525 and has a broken `msa`" with:
> because the published 0.1.1 recipe has no walker floor (the only 0.1.1 biocontainer runs logan-walker 0.1.0) and an unpinned `mafft` that resolves 7.525 when bioconda is ranked above conda-forge.

**B2. "Why 0.2 and not 0.1.2", bullet 2 and closing line.** Replace the (b) bullet with:
> Item (b) is staged. Stages 0 and 1 (queue depth, allocator threshold, an accession-count-aware bound, compressed bytes read from disk instead of RAM) change no output and go into 0.1.2. Stages 2 and 3 rewrite parse/extract. Only stage 3, with a subgraph-size cap, turns `too_large` into a transfer and disk guard.

Replace "If only (d) lands and nothing else, 0.1.2 is defensible. If (a), (b) or (c) land, go to 0.2.0." with:
> 0.1.2 ships first and does not wait for any item here. It carries the critical and high fixes from the 2026-11 review, (b) stages 0-1 and (d). (a), (c) and (b) stages 2-3 make 0.2.0, with the walker floor moved to `>=0.2,<0.3`.

**B3. Item 0, second verified bullet.** Replace from "`conda create --dry-run`" to "has a broken `msa`." with:
> `conda create --dry-run -c bioconda -c conda-forge disassembler=0.1.1` (bioconda first, flexible priority) resolves mafft 7.525. 7.525 has a bug in `--addfragments`, the option `disassembler msa` passes, fixed in 7.526 by upstream commit b4cc76f0. With strict priority that channel order does not solve at all. The documented order (`-c conda-forge -c bioconda`) gets 7.526 under either priority, and both published disassembler images contain conda-forge mafft 7.526. Users who followed `environment.yml:3` (`conda install -c bioconda disassembler`) can have 7.525.

**B4. Item 0, third verified bullet.** After "cannot install disassembler at all", add:
> Removing the skip does not fix this. GraphAligner 1.0.17b's two osx-64 builds depend on protobuf 3.14.0, which exists only for Python ≤3.9, and there is no osx-arm64 GraphAligner. disassembler stays linux-64 only until a protobuf-free osx-64 GraphAligner build exists.

**B5. Todo 2.** Replace with:
> Keep in #71068: `mafft >=7.526`, run deps trimmed to biopython + requests, walker floor `>=0.1.1,<0.2`, osx-arm64 on the sieve. Change: keep `skip: True  # [osx]` on logan-walker, with the GraphAligner reason (dropping it publishes an osx-64 walker with which disassembler cannot be installed), and add `extra: autobump: enable: false` to all three recipes.

**B6. Todo 3.** Replace "none of which touch MAFFT, so a bad solve passes CI. This is what let the unpinned mafft ship unnoticed." with:
> none of which touch MAFFT. Add a `disassembler msa` smoke test on a three-record fixture (needs `test: source_files:`). It will not catch a channel-order problem: Bioconda's build solve ranks conda-forge first and got 7.526, so only the `mafft >=7.526` pin protects users. What did ship unnoticed is logan-walker 0.1.0 inside the 0.1.1 image. `disassembler walker --version | grep -q logan-walker` accepts any version, so the floor belongs in the recipe and in a run-time check (`MIN_WALKER_VERSION` in `src/disassembler/resolve.py`).

**B7. Todo 4.** Replace with:
> Done 2026-11-10. `biocontainers/disassembler:0.1.1--pyhdfd78af_0` contains mafft 7.526 (conda-forge h4bc722e_0) and **logan-walker 0.1.0**: it was built 2026-11-02, before walker 0.1.1 reached the channel on 2026-11-05. `walker`, `fetch`, `reconstruct` and `from-kmindex` run from that image behave like 0.1.0 (`--primary-metric cov`, no prune-span guard). After the build-1 recipe merges, read the `_1` image's `conda-meta`, move every pin to `0.1.1--pyhdfd78af_1`, and list `_0` as do-not-use in `docs/packaging.md`. `biocontainers/logan-walker:0.1.1--hfa8f182_0` contains walker 0.1.1 and is fine.

**B8. Todo 5, last sentence.** Replace with:
> Switch the bot off per recipe with `extra: autobump: enable: false` (bioconda-utils' `ExcludeDisabled` honours it). With the bot off, a forgotten manual PR means no release reaches Bioconda. The checklist must therefore open the recipe PR right after tagging, then verify the published package: a dry-run solve, plus the image's `conda-meta`.

**B9. (a) heading and second paragraph.** Heading: "(a) Record what each accession's analysis read. About 150 to 250 lines with tests. OPEN." Replace "The fetcher already reads `Content-Length` … into the per-accession `summary.json` entry." with:
> Record, per accession, the source (S3 or cache), URL, compressed size, ETag, Last-Modified, fetch time, and the MD5 and SHA-256 of the bytes parsed. The hashes are the part that works everywhere: a `--unitigs-cache`/`--offline` run, `disassembler fetch` and a Galaxy-staged collection have no HTTP response, so headers alone would record nothing on the reproducible route. On `u/` objects the ETag is the plain MD5 of the content (50 of 50 sampled), so a recorded MD5 can be checked against S3 later with a HEAD. Last-Modified is a bulk-copy date (31 Mar 2024 for almost every object), not a build date. Downloads write a sidecar next to the cache entry so that cache hits can report the original headers. Also refuse a cache file whose records name another accession.

**B10. (a), "Also worth taking for free".** Replace the paragraph with:
> The `stats/logan_accessions_v1.1_SRA2023.csv.zst` (5.1 GB) and `v1.2_SRA2025.csv.zst` (11.8 GB) manifests are SRA metadata tables with no object hash. `u/` has no versioned prefix (`c1.0/` is contigs only), and it gained objects on 2025-12-14. Naming a manifest pins none of the bytes the walker reads. At most, record one `ListObjectsV2` of `stats/` (key, ETag, size, date) at start-up as a label of which releases existed; the per-object hash is the pin.

**B11. (b), point 3.** Replace with:
> 3. `stages/extract.rs` then BFS-extracts the anchor neighbourhood and discards the rest. At 0.1.1 defaults (k=21, 3 shared k-mers, 5 hops) it keeps 747 of 1,277,362 unitigs on ERR2028820/pfcrt (0.06%) but 34,582 of 81,534 on SRR8209082/hiv1_pol (42.4%). The old "26k" figure was SRR8209082 under BLAST anchoring, and its whole graph is 81,534 unitigs, not 1.28M.

**B12. (b), "The Rust CSR is ~34 B/unitig".** Replace the sentence and the next one with:
> The Rust CSR's arrays cost 12 + 4 × (links per unitig) + (bases per unitig)/4 bytes per unitig: 33.5 B on ERR2028820 and 38.2 B on the link-dense SRR8209082 (35.7 B by peak RSS on ERR2028820). The plan overstates this 18 to 21 times.

**B13. (b), the constant.** Replace "(measured 1.25x on ERR2028820, 1.62x on SRR15477033, plus the parser's doubling slack)" with:
> (1.62x on SRR15477033 rounded up; the measured ratios already include Vec growth). It was measured one accession per run. From the second accession on, glibc's raised mmap threshold makes the steady-state graph 3.17x its compressed size (ERR2028820: 102.7 MiB against 43.5 MiB), and real runs exceed the printed bound by 7 to 19%. Pinning `M_MMAP_THRESHOLD` restores 1.34x. Nothing after parse is counted (subgraph sequences, `phase`, the GFA on tmpfs scratch, GraphAligner, retained results), so the printed figure is an estimate, not a bound.

**B14. (b), "So a capped 1 GB accession … got larger, not smaller."** Replace with:
> The walker's own estimate for one in-flight 1 GB accession is 2.7 GB (one buffer plus the graph). The walker's real memory did not change at 0.1.1 (`parse.rs` changed only in comments), so PR #8800's 43 MB model was already too low at 0.1.0. Of the 33.6 GB printed at the defaults, 20.0 GB is compressed buffers held in 4 + F + W places and 13.6 GB is graphs. The larger term is the buffers, and removing them needs no parser change.

**B15. (b), "Payoffs, unchanged".** Replace with:
> Payoffs, by stage. Stages 0 and 1 cut the default bound from 33.6 to 13.6 GB and measured peaks by up to 10x (530.8 → 50.7 MiB over 20 cached accessions at `-w 1 --fetch-workers 8`), with no output change. Only stage 3 (adjacency spooled to disk) makes memory independent of accession size, and then the subgraph needs its own cap (`--max-subgraph-nodes`), because on amplicon-like data the subgraph is most of the graph. Until stage 2 removes the `u32` sequence offsets, a raised cap also stops at 2^32 unitig bases (2.1 to 2.4 GB compressed).

**B16. (b), first "ruled out" item.** Replace "(the object is a flat unindexed `.zst`)" with:
> (there is no content index, and the anchors are known only after every unitig's sequence has been scanned). The objects are not flat: they are multi-frame zstd, about 128 MiB decompressed per frame with frames ending on record boundaries, followed by a seek-table frame and a frame-to-first-record index. That allows frame-parallel decoding, knowing the decompressed size before decoding, and detecting truncation. It is not a published Logan contract, so any use needs a full-stream fallback.

**B17. (b), trivial adjacent fix.** Replace with:
> The constant 1.7 is quoted as ~1.5x (`pipeline.rs:32`), ~1.6x (`pipeline.rs:77`, `:505`, `docs/walker.md:237`), 1.7 (`docs/walker.md:249`) and 1.24x to 1.62x (`types.rs:94`). Make all five agree, and say what the figure leaves out.

**B18. (c) heading and last paragraph.** Heading: "(c) Emit a per-accession status table. Small (about 150 lines with tests). OPEN." Replace "A `skipped.tsv` (accession, outcome, size) …" with:
> Write `<out-dir>/status.tsv` on every run: one row per input accession in input order (accession, status, class, compressed size, source and hash from (a), counts, detail). Write it before the exit-status check, so that it also exists on exit 1 and 3, and add a top-level `counts` object to `summary.json`. Also needed:
> - the compressed size as a field (today it exists only inside the `error` text);
> - the same table from `--fetch-only` (`--status-tsv FILE`);
> - negative markers in the cache, so that a fetch-then-`--offline` run keeps `no_logan_data`/`too_large` instead of `error`;
> - `too_large`, not `error`, for accessions over the parser's 2^32-base limit;
> - in the Galaxy wrapper, a `--max-size-mb` parameter, without which a workflow cannot retry with a raised cap.
>
> The larger silent drop is that a run exits 0 as soon as one accession succeeds, however many errored (`check_run`, `main.rs:405`). Print per-status counts in the closing line and add `--max-errors N` (exit 4).

**B19. (d), after the second paragraph.** Add:
> Also give `msa` an `--aligned` check mode, so a split workflow keeps the 0.1.1 length and identity gates. The Galaxy MAFFT tool is `rnateam/mafft` (`rbc_mafft_add` 7.526+galaxy2, `--addfragments --keeplength`), not an IUC tool. IUC `fasttree` defaults to the Protein model. `docs/galaxy.md:88` is wrong as well: `--skip-tree` drops the MSA, not just the tree.

**B20. (e), whole section.** Retitle "3. Sieve: fix the wrapper now; upstream asks". Replace the first two paragraphs and the bullet list with:
> The admission gate is not an upstream blocker. `sieve_run` passes `--warmup 2000` (`sieve.py:214`), so after 2,000 records the binary rejects any record whose count of new k-mers is an outlier. A rejected record is never added, so its clade follows it. This happens with unique records too: a 0.67%-divergent accession after 2,100 unique haplotypes is dropped, and is kept when placed first. Byte-exact dedup does not cure it. Fix it here: pass `--warmup 2147483647` so the gate never arms, and fail when the binary reports `REJECT (error)` above 0. Keep the dedup pre-pass as a speed-up only.
>
> Wanted from upstream (`stevenweaver/kmer-diversity-sieve`):
> - grow the promoted-set LSH table: it has a fixed 2^20 slots, and more than 65,536 divergent records hang the binary forever (`l1_pipeline.c:155`); a patch exists;
> - lift or document the 65,536 k-mer sketch cap (`:264`) and the month bucketing from `|YYYY-MM` headers;
> - a real tag. Checked 2026-11-10: `git ls-remote --tags` is still empty.

Replace the "New, and worth resolving in 0.2" paragraph with:
> The recipe builds `70022da`, which is upstream `main`; `39d90c0` is its parent, and the two differ only by 56 README lines. The disassembler submodule went backwards by accident: 00e3b97 advanced it to `70022da`, and the docs commit b1b6acf set it back to `39d90c0`. Move the submodule forward to `70022da` in its own commit, leave the recipe alone, and keep pressing for a tag.

Replace the last paragraph with:
> The interim sieve UDT wraps `disassembler sieve` in the disassembler container, which has Python and the sieve. The sieve container has no Python, and the bare binary with its archive-scale defaults keeps 0 of 33 records of the HIV-1 pol test cohort, with exit 0. A later bare-binary tool must set all five gene-scale options, derive `--min-len` from the reference, concatenate the rank shards, pass `--warmup 2147483647`, and fail on an empty result.

**B21. Section 4.** In the 1(b) bullet, replace "which makes the same argument more strongly" with:
> At the defaults, 20.0 of the 33.6 GB bound is compressed buffers and 13.6 GB is graphs. That moves the argument: the cheaper, larger saving is keeping compressed bytes out of RAM (stages 0-1), not a single-pass parser.

To the UDT-table bullet, add: "but `disassembler:0.1.1--pyhdfd78af_0` bundles logan-walker 0.1.0; pin `0.1.1--pyhdfd78af_1` once rebuilt." Replace the mafft-row bullet with:
> Settled: both disassembler images contain mafft 7.526. The row is moot, because ToolShed `rnateam/mafft` ships `rbc_mafft_add` 7.526+galaxy2 on `mulled-v2-027f30e88c6da4ab4b0170ff590dbab158b07243:43ff19c348e21674cd51004463f1ea0b2e5f2336-0`.

In "Unchanged and still correct", delete "the subgraph scale figures", and replace "the ~71 MB median compressed accession" with:
> the per-accession size: the repository gives ~71 MB (`stream.rs:7`), about 50 MB (`docs/galaxy.md:27`, gene-hit cohorts) and a 0.6 MB median for random accessions (`docs/kmindex.md:82`). Size frozen-input quotas from the cohort's own sizes.

**B22. Section 5.**
- Sieve-tests bullet: replace "Cheap fix: …" with: "Fix: in `setUpClass`, run the binary once. On a signal death, rebuild with `make serial` and retry. Skip only outside CI and only if the rebuild fails: a skip in CI tests nothing."
- Linter bullet: replace "Would need its own pass over the existing code." with: "Small: pyflakes rules find 13 one-line issues. Also remove the crate-wide `#![allow(dead_code, unused_variables)]` (`walker/src/main.rs:24`), which hides nine dead or test-only items from clippy."
- `harness.py` bullet: add "`REPO_ROOT` resolves to the directory above the checkout, and a missing genome silently becomes a random sequence."
- PR #8800 bullet: replace "Its whole … list is still due, and its walker memory block … supersedes twice over." with:
> The 0.1.1 bump is not a token change. Three PR tests fail at 0.1.1: the harvest and stratify test data are 0.1.0-shaped, and walker test 3 asserts `HTTP 404`, which is now `no_logan_data`. The walker tool cannot select `--primary-metric matched`, and a Reference CDS dataset named `reconstructed` is overwritten by the output. The memory block counts 43 MB per graph where the walker's model is 1.7 × cap plus one buffer: fix it now (2.7 × cap per compute worker), and later pass a walker `--max-memory-mb` instead of copying the model.
- Replace "check whether that ever went up" with: "It went up (PR head `disassembler_msa.xml:90-93`, `disassembler_lexicmap.xml:246`) and is wrong again at 0.1.1: under `--require-full-length` a trusted in-frame deletion counts as covered, so the rows are not gap-free."

**B23. "Verified already fixed at 0.1.1".** Replace the last sentence with:
> These fixes hold, with the exceptions below, but 0.1.1 is not correct. The default majority correction never applies unless the alignment spans the whole stitched path (WV-01, critical). The 2026-11 review lists 14 high findings, among them a 0.1.1 regression (`nohup` no longer protects the walker, WF-02). Scratch handling is not fixed across PID namespaces (WF-01). The composites still turn exit 3 into 1 (PC-10). The spike panel has two more mislabelled records (DD-01, DD-02).

**P1. Section 1 preamble.** Replace "Ordered by effort. (a) and (b) are the ones we consider blockers…" with:
> Effort: (a) medium (150-250 lines), (b) large and staged, (c) small, (d) medium. For Galaxy, silently omitted accessions are fixed by (c) plus a wrapper cap parameter, and memory headroom by (b) stages 0-1. (a) is required for any reproducibility claim.

**P2. 1(a).** Replace the second paragraph's last sentence and the "Additionally" paragraph with the texts of B9 and B10.

**P3. 1(b).** Replace "(~700 B/unitig RAM)" with "(33.5 to 38.2 B/unitig of arrays)". Replace "discarding ~99.9% … ~1.28M unitigs)" with "discarding most of it (at 0.1.1 defaults: 747 of 1,277,362 unitigs kept on ERR2028820, but 34,582 of 81,534 on SRR8209082)". Replace the paragraph "If the zstd stream were decoded incrementally …" and its two payoff bullets with:
> The stream is already decoded incrementally; what is buffered is the compressed object. One pass cannot reduce to the subgraph: anchors are known only once every sequence is scanned, and the 5-hop BFS needs the edges of unitigs anywhere in the stream (expanding as it goes, a single pass reaches 354 of 747 nodes on ERR2028820). O(subgraph) memory needs two decodes and the adjacency spooled to disk; adjacency kept in RAM is still O(all unitigs), at 12 to 20 B/unitig. Plan it in four stages (below); stages 0-1 remove the larger memory term without touching the parser.

Replace "What does NOT work (checked)…" with the text of B16. Replace the Rationale with: "Rationale: stages 0-1 fix memory headroom cheaply; silently omitted accessions are fixed by (c) and a wrapper cap parameter."

**P4. 1(c).** Replace "— trivial" with "— small", and append the B18 list.

**P5. 1(d).** Replace "`--skip-tree` exists only inside the composites" with:
> `--skip-tree` exists only in `reconstruct` and `from-kmindex`, and there it stops after the walker: the sieve, the MSA and the tree are all skipped.

**P6. 1(e).** Replace the body with the B20 texts.

**P7. 1(f).** Retitle "(f) Recipe hygiene (bioconda-recipes #71068, open)". Replace "7.525 segfaults on the `--add` path `msa` uses" with "7.525 has a bug in `--addfragments`, which `msa` uses (fixed in 7.526, upstream b4cc76f0)". Replace "A bioconda-priority solve resolves the documented-broken version." with "A bioconda-first solve with flexible priority resolves it; the documented conda-forge-first order gets 7.526." Replace the test bullet with: "Recipe tests check neither MAFFT nor the walker version, and the walker check is the one that mattered: the 0.1.1 image ships walker 0.1.0." Add bullets: "`logan-walker >=0.1.1,<0.2`", "`extra: autobump: enable: false` on all three recipes", "keep logan-walker's osx skip (GraphAligner 1.0.17b)".

**P8. Section 2 rows the backlog depends on.**
- Container column: replace `0.1.0--pyhdfd78af_0` with `0.1.1--pyhdfd78af_1` (after the rebuild; never `_0`).
- `disassembler-fetch` row: replace "same walker container" with "disassembler container (`fetch` is the Python wrapper; the walker image has no `disassembler`), or wrap `logan-walker --fetch-only --unitigs-cache` directly in the walker container" (GX-26).
- `disassembler-codon-align` row: replace "needs change (d); until then wrap `msa` whole and discard the tree" with "after (d): `msa --prep-only` → `rnateam/mafft` `rbc_mafft_add` → `msa --aligned --no-tree` → IUC `fasttree` (Nucleotide, GTR). Never discard the tree with `--out-tree /dev/null`: run as root, `msa` replaces the device node" (PM-09).
- sieve row: the B20 UDT text (GX-13).
- mafft row: the B21 mafft text (GX-23).

## 4. Designs

### 4(a) Provenance: record what each accession's analysis read

**What to record.** A per-accession `input` object in `summary.json`, written before `timing` (nothing existing moves):

```json
"input": {
  "source": "s3",
  "url": "https://logan-pub.s3.amazonaws.com/u/ERR2028820/ERR2028820.unitigs.fa.zst",
  "cache_path": null,
  "compressed_bytes": 34016719,
  "etag": "89d0255e8d9b262f4f0ec26bcf53ae20",
  "last_modified": "Sun, 31 Mar 2024 15:36:42 GMT",
  "fetched_at": "2026-11-10T22:36:41Z",
  "md5": "89d0255e8d9b262f4f0ec26bcf53ae20",
  "sha256": "a03b08561b2f7233f3ac48c08825199d37ed801cfa8264151e04b6f7c92f187f"
}
```

Rules:
- `source` is `s3` or `cache`. `md5` and `sha256` are always of the bytes this run parsed, so they are filled on every route: S3, cache, `--offline`, a Galaxy-staged collection.
- `etag`, `last_modified` and `fetched_at` come from the response, or from the sidecar on a cache hit. They are `null` for a cache filled by 0.1.x, by `curl` or by Galaxy.
- `"input": null` when no bytes and no response exist (`no_logan_data`, transport error).
- A `too_large` refusal made by Content-Length has a response. It records `compressed_bytes`, `etag`, `last_modified` and `fetched_at`, with null hashes. That also gives 4(c) its size column.
- ETag quotes are stripped. A multipart ETag (`…-N`) is recorded but never compared with the MD5.

Run level, a new top-level object next to `settings` (so "exactly 19 settings" stays true):

```json
"input_source": {"base_url": "https://logan-pub.s3.amazonaws.com/u/", "unitigs_cache": "/data/cache", "offline": true, "content_hash": true}
```

**Cache sidecar.** `write_cache_atomic` also writes `<ACC>.unitigs.fa.zst.meta`, using the same tempfile-and-rename. The crate writes JSON by hand and has no JSON reader, so the sidecar is `key<TAB>value` lines:

```text
schema	1
accession	ERR2028820
url	https://logan-pub.s3.amazonaws.com/u/ERR2028820/ERR2028820.unitigs.fa.zst
compressed_bytes	34016719
etag	89d0255e8d9b262f4f0ec26bcf53ae20
last_modified	Sun, 31 Mar 2024 15:36:42 GMT
fetched_at	2026-11-10T22:36:41Z
md5	89d0255e8d9b262f4f0ec26bcf53ae20
sha256	a03b08561b2f7233f3ac48c08825199d37ed801cfa8264151e04b6f7c92f187f
logan_walker	0.2.0
```

On a cache hit:
- Copy `etag`, `last_modified` and `fetched_at` from the sidecar.
- Fail the accession ("cache entry does not match its sidecar") if the size or a hash differs.
- Without a sidecar, those three fields are null, and the hashes are still there.

`--revalidate-cache` (online only) sends a HEAD and compares ETag and length with the sidecar, or with the computed MD5. On a difference it warns and refetches.

**Two checks that cost nothing** (both fit 0.1.2):
- `parse::build` keeps the first header's name prefix, which it discards today (`parse.rs:111`). `run_one` then fails the accession when the prefix is present and is not the accession: today ERR2028820's bytes filed as SRR8209082 are analysed as SRR8209082 with exit 0 and no warning; with the check they are refused. Files with bare numeric ids (DRR344513: `>0 ka:f:9500.0`) are not checked.
- Add the object MD5 to `bench/golden/*.golden.json` and print it in `bench/check_golden.sh`, so the gate can tell "Logan changed the object" from "the walker regressed". Values: ERR2028820 `89d0255e8d9b262f4f0ec26bcf53ae20`, SRR8209082 `d2dd288b238c3c88b46acdcd78860dce`, SRR12363134 ETag `6818293beadb1d53dd73fba79b266183` (not downloaded).

**Release-level label (optional).** Do not resolve cohorts against the 5 to 12 GB manifests. If a label is wanted, record one `ListObjectsV2` of `stats/` at start-up (about 2 KB) as `logan_stats`: key, ETag, size and date of each manifest. It shows which releases existed when the run was made, which decides `no_logan_data`. It does not pin bytes; the hash does.

**Code to touch.**

| where | change |
|---|---|
| `walker/src/types.rs:62` `FetchedUnitigs` | add `input: InputProvenance` (`source, url, cache_path, compressed_bytes, etag, last_modified, fetched_at, md5, sha256`) |
| `types.rs:448` `AccessionResult`, `:546` `AccessionResult::empty` | `input: Option<InputProvenance>`; `empty` takes it |
| `stages/stream.rs:172` `try_once`, `:189` | read `etag` and `last-modified` beside `content_length()`; `FetchAttempt::TooLarge` and `Error::TooLarge` carry `compressed_bytes: Option<u64>` and the headers |
| `stream.rs:74` `fetch` (cache branch `:77-109`), `:249` `write_cache_atomic` | hash on both paths, while reading or copying; write and read the sidecar; `--revalidate-cache` |
| `pipeline.rs:536` | clone `fetched.input` before `drop(fetched)`; set it on both returns (`:549`, `:914`) |
| `pipeline.rs:457` `compute_worker`, `:252` `unfetched`, `:186` `run_fetch_only` | keep a copy outside `catch_unwind` so a parse failure or panic still reports what it was fed; pass the refusal's headers; return per-accession records |
| `stages/parse.rs:111` | keep the first header prefix in `UnitigGraph` |
| `io/output.rs:245` `write_result`, `:199` `write_run_info` | print `input` and `input_source`; add a small RFC 3339 formatter (about 15 lines) |
| `main.rs` | `--revalidate-cache`, `--no-content-hash`, usage text |
| `walker/Cargo.toml` | move `sha2 = "0.11.0"` from `[dev-dependencies]` to `[dependencies]`; add `md-5` from the same family |

**Staged order.**
1. (0.1.2) Prefix check and golden MD5s.
2. (0.2.0) `InputProvenance`, hashing, headers, `input` and `input_source`.
3. (0.2.0) Sidecar write and read, then `--revalidate-cache`.
4. Optional: `logan_stats`.

**Tests.**
- `stream.rs` unit tests, using the existing local-server helper with headers: `etag_and_last_modified_are_captured`, `a_too_large_refusal_records_size_and_etag`, `a_download_writes_a_sidecar`, `a_cache_hit_reports_its_sidecar`, `a_cache_hit_without_a_sidecar_still_hashes`, `a_cache_entry_that_disagrees_with_its_sidecar_is_an_error`.
- `output.rs`: `summary_json_records_the_input`.
- `tests/cli.rs`: `an_offline_run_records_cache_source_and_hash`, `a_cache_file_for_another_accession_is_refused`.
- Ignored network test: `etag_is_the_md5_of_a_real_object`.

**Docs.**
- `docs/walker.md`: Outputs (the `input` object), Offline cache (sidecar, revalidation), and a new first Reproducibility bullet: Logan inputs are mutable, and the hash identifies them.
- `docs/galaxy.md`: walker notes.
- `bench/README.md`: object MD5 in the fixtures.

**Compatibility.** Additive. No file in `src/`, `tests/` or `bench/` reads the walker summary. PR #8800's tests assert the `success` property and the text `HTTP 404`; the latter already fails at 0.1.1 (GX-10). `fetched_at` makes two summaries differ, but `timing` already does. The Galaxy wrapper needs no change. For the frozen-input workflow, the hash is what makes a collection-staged run auditable, because Galaxy-staged files have no sidecar.

**Effort.** Medium: 150 to 250 lines with tests. Hashing costs 0.064 s (MD5) and 0.028 s (SHA-256) on the 34 MB ERR2028820 object, against 1.27 s of parse; `--no-content-hash` opts out.

**Risks.**
- Multipart ETags are not MD5s; objects written 2025-12-14 also expose `x-amz-checksum-crc32`.
- Last-Modified is weak and should be documented as such.
- Comparing MD5 with the ETag at download time is optional: Content-Length and zstd frame checksums already catch a cut transfer. The sidecar hash is what catches a hand-filled cache file cut at a frame boundary (WG-10).

The same release should finish the provenance picture: the walker commit in `--version` and the summary (PK-07), reference hashes rather than id and length only (WV-20), tool versions from `msa`, `sieve` and the composites (PM-16), and which GraphAligner build ran (TC-27).

### 4(b) Memory: what streaming can and cannot do, and the staged design

**What a single decode pass can and cannot do.**

The walker already decodes the stream incrementally (`stream.rs:298-309`, `parse.rs:77`). What it buffers is the compressed object, held in up to 4 + F + W places (WF-08).

In one decode pass it can:
- check every record (ACGT, dense ids);
- count `n_unitigs`;
- count anchor k-mers on the text as it streams;
- record every unitig's id and links, at 12.1 B/unitig (ERR2028820) or 19.9 B/unitig (SRR8209082) as a u32 CSR.

It cannot know which sequences to keep. The anchors are complete only at the end of the stream, and the 5-hop BFS reaches unitigs anywhere in the file, including before the anchor that reaches them. A pass that expands as it goes reaches 354 of the 747 nodes on ERR2028820 and holds the sequence of only 289. On SRR8209082 it reaches 23,709 and holds 18,424 of 34,582.

So a single pass must keep either every sequence (today's design: 13.45 / 10.30 B/unitig of 2-bit sequence plus 4 B offsets and 4 B abundance) or the bytes for a second decode. Those bytes belong on disk (the cache file or a spool), not in RAM.

Feeding the HTTP body straight into the decoder (design D) removes the compressed buffer but keeps the full CSR. It also ties a compute worker to a socket, and the bytes still have to land somewhere for the cache, for retries and for a second pass. It is the download half of stage 1, not a design by itself.

**Measured bytes per unitig.** CSR array sizes from a streaming count; `types.rs:97-127`.

| | ERR2028820 | SRR8209082 |
|---|---|---|
| compressed / decompressed bytes (ratio) | 34,016,719 / 145,673,706 (4.28) | 1,861,058 / 9,873,585 (5.31) |
| unitigs / links / bases | 1,277,362 / 2,582,206 / 68,745,081 | 81,534 / 323,387 / 3,359,745 |
| compressed B/unitig; decompressed B/unitig | 26.6; 114.0 | 22.8; 121.1 |
| adjacency (`offsets` + `targets`) | 4.00 + 8.09 = **12.09** | 4.00 + 15.87 = **19.87** |
| sequence (`seq_2bit`) | **13.45** | **10.30** |
| `seq_offset` | **4.00** | **4.00** |
| abundance | **4.00** | **4.00** |
| total | **33.54** B/unitig = 42.8 MB = 1.26x compressed | **38.17** B/unitig = 3.1 MB = 1.67x compressed |
| walker peak RSS, one accession, no anchors | 81.6 MiB (graph 43.5 MiB = 35.7 B/unitig = 1.34x) | 12.9 MiB (graph 5.4 MiB = 3.05x) |

**Measured peak memory of today's walker.** Offline runs over 20 hard links to ERR2028820 (32.44 MiB), with a query that finds no anchors. Units are MiB.

| run | 0.1.1 peak | walker's own model at the true size |
|---|---|---|
| 1 accession, `-w 1 --fetch-workers 1` | 81.6 | 93.3 |
| ≥8 accessions, `-w 1 --fetch-workers 1` | **303.0** | 255.5 (exceeded by 19%) |
| 20 accessions, `-w 1 --fetch-workers 8` | **519.8** | 482.6 (exceeded by 8%) |
| 20 accessions, `-w 8 --fetch-workers 8` | 967.0 | 1,095.7 |
| as row 2 with `MALLOC_MMAP_THRESHOLD_=131072` | 243.9 | 243.8 |
| 8 × SYNHIV40 (73.0 MiB each), `-w 1 --fetch-workers 1` | **607.8** | 568.0 (exceeded by 7%) |

SYNHIV40 is SRR8209082 replicated 40 times with shifted ids (3,261,360 unitigs). With the hiv1_pol query its subgraph is 1,383,280 nodes. Its single-accession peak is 472.2 MiB, against 200.0 MiB with no anchors and a 202.9 MiB model: the stages after parse are not counted.

**Designs prototyped outside the walker.** Each uses the walker's header parsing, k-mer index, BFS and GFA layout. Every design wrote a GFA byte-identical to the 0.1.1 walker's on ERR2028820/pfcrt (81,340 bytes) and SRR8209082/hiv1_pol (4,852,423 bytes), and all designs agree on SYNHIV40 (216,839,199 bytes). Times are parse + anchor + extract.

| | walker 0.1.1 | A: compressed in RAM, two passes | B: compressed on disk, two passes | C: frontier passes from disk | E: B + adjacency spooled to disk |
|---|---|---|---|---|---|
| ERR2028820/pfcrt (747 nodes) | 81.5 MiB; 1.99 s | 54.4; 1.67 s | 22.1; 1.55 s | **7.6**; 2.43 s, 6 decodes | **9.2**; **1.47 s** |
| SRR8209082/hiv1_pol (34,582 nodes, 42%) | 15.5 | 13.8 | 12.3 | 16.2, 5 decodes | 12.3 |
| SYNHIV40/hiv1_pol (1.38M nodes, 42%) | 472.2 (includes `phase`); 4.45 s | 281.3; 4.35 s | 208.7; 4.38 s | 349.6; 7.60 s | 208.4; 5.29 s |
| SYNHIV40, no anchors | 200.0 | | 69.1 | 7.1 | 9.0 |

One extra light decode costs 0.25 s per 146 MB of text (about 0.55 GB/s). The two-pass designs are not slower than today, because pass 1 no longer packs 2-bit sequence or parses abundances for every unitig. C is not recommended: its decode count grows with `--hops` (up to hops + 1), and it has the worst memory on subgraph-heavy accessions.

**The stages.** F is fetch workers, W compute workers, cap is `--max-size-mb`; defaults are F = W = 8 and cap = 1024 MB.

| stage | change | code | measured saving | bound afterwards | release |
|---|---|---|---|---|---|
| 0 | `channel_capacity` 4 → 0 (a fetch worker blocked in `send` already is a queue slot). `mallopt(M_MMAP_THRESHOLD, 1 << 20)` at the top of `main` on glibc. Count-aware bound: `min(n, F + W)` buffers, `min(n, W)` graphs, and in `--fetch-only` F buffers and no graphs, printed after the mode is known and naming what it leaves out. Budgeted parse: return `TooLarge` before a doubling exceeds 1.7 × cap. Make the five 1.5x/1.6x/1.7 quotes agree | `config.rs:472`, `:407`, `:424`, `:432`; `main.rs:335`, `:497-515`; `parse.rs:71`; `pipeline.rs:482`; `pipeline.rs:32`, `:77`, `:505`, `types.rs:94`, `docs/walker.md:237` | 20 × ERR2028820 at `-w 1 --fetch-workers 1`: 304.0 → 174.3 MiB (queue, wall time unchanged at 34.1 s) → 114.7 MiB (allocator). Budgeted parse: a dense 2.4 MiB file became `too_large` at 15.9 MiB peak; 8 × ERR2028820 peaked at 243 MiB, inside the printed 254 | (F + W + 1.7 W) × cap: 29.6 GB at defaults, computed from the walker's formula; 2.7 × cap for a one-accession run (today 33.6 GB at defaults, 7.7 × cap at one slot) | 0.1.2 |
| 1 | Never hold compressed bytes in RAM. `FetchedUnitigs.compressed: Vec<u8>` becomes `enum Payload { File(std::fs::File), Memory(Vec<u8>) }` (the second for tests). A cache hit opens the file instead of `std::fs::read`. A download is copied from the body reader, limited to cap + 1, into the cache tempfile or into a spool file, hashing as it goes (4(a)). Open once and keep the handle. The spool defaults to the cache directory, else the system temp directory, never `/dev/shm`. `--fetch-only` stops reading cached files whole | `types.rs:67`; `stream.rs:92`, `:205-209`, `:249`, `:298`; `pipeline.rs:186` | Cache hits on a build with stage 0's allocator setting: `-w 1 --fetch-workers 8` **530.8 → 50.7 MiB**; `-w 8 --fetch-workers 8` **954.8 → 348.0 MiB** (exactly eight 43.5 MiB graphs). Summaries identical apart from timings. A second measurement with a 46.4 MB file: 355.6 → 92-97 MiB at one slot | 1.7 × W × cap: **13.6 GB** at defaults, 1.7 GB at one slot. Disk: up to (F + W) × cap of spool when no cache is set | 0.1.2 |
| 2 | Two passes, subgraph-local data downstream (design B). Split `parse::build` into `scan(reader, sink, on_record)` (ids, links into an adjacency sink, ACGT check, k-mer count on the text) and `collect(reader, nodes)`, a second decode that stops at the last wanted id and returns `SubgraphData { nodes, seqs, abundance, edges, n_unitigs }`. `anchor::find` becomes a per-record `shared_kmers(seq)`. Replace `&UnitigGraph` by `&SubgraphData` in the 19 places it is taken. Keep the unpruned nodes and sequences for the refused-prune branch instead of re-expanding. Keep the old path behind `--legacy-parse` for one release | `parse.rs:71`; `anchor.rs:255`; `extract.rs:53`, `:127`, `:146`, `:169`, `:203`, `:220`, `:280`; `phase.rs:65`; `variants.rs:83`, `:154`, `:304`, `:538`; `pipeline.rs:729-845` (re-expansion at `:843-844`) | Design B (stages 1 + 2) against today: 81.5 → 22.1 MiB on ERR2028820. Against design A (two passes, compressed still in RAM): 54.4 → 22.1 MiB on ERR2028820, 281.3 → 208.7 MiB on SYNHIV40 | About 0.9 × W × cap (adjacency is 0.45 to 0.87 × compressed). The 2^32-base ceiling goes, leaving 2^32 edges and 2^30 ids | 0.2.0 |
| 3 | Adjacency to the spool (design E). The sink writes `offsets` (u64) and `targets` to two spool files, and `expand` reads only the rows it visits (`read_exact_at`). Add `--max-subgraph-nodes` with a status of its own. Then raise or drop the default cap, and state the disk need: compressed size plus 0.6 to 1.05 times it for the adjacency, per in-flight accession | as stage 2, plus the sink and `main.rs` | ERR2028820 9.2 MiB; SYNHIV40 with no anchors 9.0 MiB. But SYNHIV40 with anchors (1.38M-node subgraph) still peaks at 208.4 MiB: where the subgraph is most of the graph, O(subgraph) is O(graph) | W × subgraph: about 150 B per node in RAM, plus 157 B per node of GFA on scratch, plus GraphAligner. No size-cap term | 0.2.0 or later |

**Largest saving for the least work: stage 1, with stage 0.** Stage 1 removes the 20.0 GB buffer term of the 33.6 GB default bound. The core is about 40 lines plus the spool, needs no parser change, measured a 10x lower peak on cache-heavy runs, and leaves outputs byte-identical. Stage 2 saves another ~0.8 × W × cap for a large rewrite. Only stage 3 delivers the plan's two payoffs.

**Separate memory items, same release as stages 0 and 1:**
- Results are retained until the end of the run, about 0.30 MiB per finished deep accession: stream them through an in-order sink (WP-09).
- Decoded size is unbounded: the stage 0 budgeted parse covers the graph, and the seek table (PK-14) can give the decompressed size up front (WG-06).
- The bound is never compared with `MemAvailable` or the cgroup limit: warn first, then add `--max-memory-mb` backed by a byte semaphore (WP-10).

**Tests.**
- Stage 0: `memory_bound_tests::bound_respects_accession_count_and_fetch_only`; `build_with_budget` on a dense fixture returns `TooLarge` and maps to `Outcome::TooLarge`.
- Stage 1: change `a_cache_hit_returns_the_bytes_on_disk` and `a_cached_file_over_the_size_cap_is_too_large` to decode and compare; a cache hit returns `Payload::File`; a spool file is removed on every error path; summaries over N cached accessions are identical to today's apart from timings.
- Stage 2: data-backed `#[ignore]` tests comparing GFA bytes between `--legacy-parse` and the new path on ERR2028820, SRR8209082 and SRR12363134 (the last not yet tried), plus `bench/check_golden.sh`; a CI job that compares the two paths on the CLI fixtures.
- Stage 3: the subgraph cap yields its own status; a full spool disk fails the accession with a file error and removes the partial file.

**Docs.**
- `docs/walker.md` Memory section: what is bounded and what is not; the spool and its disk need; the 2^32 ceiling until stage 2.
- The `pipeline.rs` module doc and `config.rs:395-445`.
- `extract.rs:6-8`: subgraph figures at the defaults (WG-14).
- The "~38 KB GFA" comments (WG-12).
- `docs/galaxy.md`: memory notes.

**Compatibility.**
- Stages 0 and 1 change no output. The startup line changes (and becomes correct for `--fetch-only` and short lists).
- Stage 2 folds `anchor` time into `parse`. Keep both `timing` keys, and document that `anchor` is 0.
- Stage 3 adds a status value. Consumers that enumerate statuses, including the 4(c) table's `class` column, must accept it.
- The Galaxy wrapper should fix its 43 MB model now (GX-08: 2.7 × cap per compute worker). After stage 1 it over-requests. The durable fix is a walker `--max-memory-mb` that the wrapper passes, instead of a second copy of the model in shell (GX-08, GX-14). After stage 3 the job needs scratch disk scaled by the cap, not RAM.

**Effort.** Stage 0: trivial to small. Stage 1: small. Stage 2: large. Stage 3: medium.

**Risks.**
- **Golden identity.** The GFA decides GraphAligner's output; it was reproduced byte for byte on two real accessions only, and GraphAligner itself could not be run here. Gate stage 2 on GFA identity for all three golden accessions.
- **Validation moves.** Today a malformed `ka:f:` anywhere fails the accession. In two passes, abundance is parsed only for subgraph nodes unless pass 1 keeps parsing it. Decide, and document the choice.
- **One open handle for both passes**, so that a concurrent rename of the cache entry cannot change the bytes between passes.
- **tmpfs.** A spool on `/dev/shm` is RAM under another name. Galaxy already points `--scratch-dir` at the job directory; the CLI default does not.
- **Throughput.** With queue depth 0 a fetch worker cannot bank a download while compute is busy. This made no difference offline. After stage 1 the queue holds file handles, so its depth can go back up.
- **Page cache.** In stage 3 the adjacency lives in page cache, which cgroups count but can reclaim. A destination that kills on total memory, including dirty pages, needs the spool flushed.

### 4(c) Per-accession status table

**What to write.** `<out-dir>/status.tsv` on every run. One header line, then one row per input accession in input order. It is written before `check_run` (`main.rs:401`), so it exists on exit 1 and 3 as well. In the example, the `too_large` row and the counts of the `success` row come from runs in this review; `…` marks values that depend on the alignment, and `HTTP 404` is new detail text (today a `no_logan_data` result carries none).

```text
accession	status	class	compressed_bytes	source	md5	n_unitigs	n_anchors	subgraph_nodes	provisional	primary_length	cov	pident	prune	detail
ERR2028820	success	reconstructed	34016719	cache	89d0255e8d9b262f4f0ec26bcf53ae20	1277362	220	747	false	…	…	…		
SRR8209082	too_large	skipped	1861058	s3					false					compressed size 1861058 exceeds max_compressed_bytes 1048576
SRR0000001	no_logan_data	not_found		s3					false					HTTP 404
```

- `class` is `reconstructed`, `not_found` (`no_anchors`, `no_path`, `no_logan_data`), `skipped` (`too_large`, and the stage 3 subgraph cap) or `error`.
- `source` and `md5` are empty until 4(a) lands.
- In `detail`, tabs and newlines are replaced by spaces.
- New columns are only ever appended.

**Same release, same data:**
- `summary.json` gets a top-level `"counts": {"success": 0, "no_anchors": 0, "no_path": 0, "no_logan_data": 0, "too_large": 2, "error": 0}`. A wrapper test or workflow conditional can assert on it (closes the summary half of WP-03).
- The closing line counts every class: `[logan-walker] 0 success, 2 too_large (skipped by --max-size-mb, not examined), 0 errored; see status.tsv`.
- Without `--allow-empty`, an all-skipped run says "every accession was skipped by --max-size-mb", not "nothing was found" (`error.rs:82-87`, WP-14).
- Exit codes do not change. A failed Galaxy step cannot feed a retry branch, so all-skipped stays exit 0 under `--allow-empty`. Add `--max-errors N` (exit 4) for runs that must not lose accessions silently (WP-03).

**`--fetch-only`.** It writes the same table when `--status-tsv FILE` is given (`disassembler fetch --status-tsv`), with status values `cached`, `no_logan_data`, `too_large` and `error`.

**Negatives survive fetch-then-offline (WF-03, prototyped).** When the fetch stage has a cache directory and gets a 404 or a size refusal, it writes `<ACC>.unitigs.absent` in the sidecar format of 4(a):

```text
schema	1
accession	SRR8209082
status	too_large
compressed_bytes	1861058
max_compressed_bytes	1048576
checked_at	2026-11-10T22:40:02Z
```

Under `--offline`, a miss with a marker reports `no_logan_data`, or `too_large` while the recorded size still exceeds the current cap. Otherwise it stays `error`. Online runs ignore markers and refresh them, and a successful fetch removes the marker. In a patched build, the offline run reported `no_logan_data` and `too_large` and exited 0, where 0.1.1 exits 3.

**Parser ceiling.** Until 4(b) stage 2, the `u32` overflow at `parse.rs:137-140` and `:211-214` returns a `TooLarge`-class error with the detail "over the parser's 4.29 Gbp limit", not `error`/exit 3 after a full download (WG-08).

**Code to touch.**

| where | change |
|---|---|
| `types.rs:411` `Outcome` | `class()` next to `is_negative()`; `compressed_bytes: Option<u64>` on `AccessionResult` (or via `input`) |
| `stream.rs:189-194`, `:110-117` | `TooLarge` carries the size (Content-Length, or cap + 1 when the body limit found it); `.absent` write and offline lookup |
| `io/output.rs` | `write_status_tsv`; `counts` in `write_summary` (`:384`) |
| `main.rs:556-581`, `:401` | write `status.tsv` before `check_run`; closing line; `--status-tsv`, `--max-errors` |
| `error.rs:82-87` | all-skipped wording |
| `pipeline.rs:186` | `run_fetch_only` returns per-accession records |
| `parse.rs:137-140`, `:211-214` | `TooLarge` class |
| `src/disassembler/cmd/fetch.py` | pass `--status-tsv` through |

**Tests.** `status_tsv_has_one_row_per_accession_in_input_order`, `status_tsv_detail_has_no_tabs_or_newlines`, `status_tsv_is_written_when_the_run_exits_3`, `all_too_large_says_skipped_not_nothing_found`, `fetch_only_writes_status_tsv`, `an_offline_miss_with_a_negative_sidecar_keeps_its_status`, `max_errors_fails_a_run_that_lost_too_many_accessions`, and `over_4gbp_is_too_large_not_error` (the synthetic file, `#[ignore]`: 13 s, 1 GB).

**Docs.**
- `docs/walker.md`: Outputs, Exit status, and Offline cache (`.absent` files; deleting one forces a refetch). The `too_large` row's "skipped before downloading" is untrue without a Content-Length or for a cache file.
- `docs/galaxy.md`: the new output, the retry branch, and what the wrapper must expose.

**Compatibility.** One new file and one new top-level key in `summary.json`; no existing key changes. For the Galaxy wrapper:
- add a `tabular` output, collected on failed jobs too (GX-16 shows the summary is lost there today);
- replace walker test 3's `<has_text text="HTTP 404"/>`, which already fails at 0.1.1, with a `status.tsv` line match or `counts.no_logan_data`.

A retry branch (filter `status == too_large`, cut column 1, a second walker step on a larger allocation) needs three more things:
- a `--max-size-mb` parameter in the wrapper, with the admin variable as a ceiling (GX-05, GX-09);
- the count-aware bound from 4(b) stage 0;
- the parser ceiling classified as above.

**Effort.** Small: about 150 lines with tests. The negative-marker half is already prototyped.

**Risks.** The column set becomes an interface, so only append to it. The all-skipped exit 0 under `--allow-empty` is deliberate and must be documented as such. Stale `.absent` markers are refreshed only by online runs.

### 4(d) `msa` without a tree, and a split-workflow mode

**CLI of the atomic `msa`** (a tested prototype exists):

| mode | flags | tools needed | writes |
|---|---|---|---|
| align + tree (today) | `--out-aln A --out-tree T` | MAFFT, FastTree | A, T |
| align only | `--out-aln A --no-tree` | MAFFT | A |
| prep only | `--prep-only --out-prep-ref R --out-prep-seqs S` | none | R: one record `>Reference`, upper case. S: the rows that passed `--msa-min-cov`, reverse-complemented where needed, ids sanitised. Prints `mafft --auto --keeplength --thread N --addfragments S R > aligned.fasta` |
| check an external alignment | `--aligned ALN --out-aln A` plus `--out-tree T` or `--no-tree` | FastTree only with `--out-tree` | A (ALN after the checks), T |
| optional, any mode | `--out-id-map TSV` | none | `sanitised_id  original_id  strand  acgt_bases`, so tree tips map back to records (ids are cut to 40 characters today, PM-26) |

Rules:
- Exactly one of `--out-tree` / `--no-tree` unless `--prep-only` (exit 2 otherwise), so a forgotten `--out-tree` cannot silently drop the tree.
- Tools are looked up only for the stages that run.
- `--aligned` requires the first record to be `Reference`, equal to `--ref` ungapped. Every row must have the reference's length ("was MAFFT run with --keeplength?"). The 0.1.1 identity check (`identity_ok`, `msa.py:109`) runs on it, and the codon QC of PM-03 belongs here too.

Without `--aligned`, a split workflow loses the length and identity gates that 0.1.1 added.

**Code to touch.**

| where | change |
|---|---|
| `src/disassembler/cmd/msa.py:138-172` `add_parser` | `--out-tree` optional; add `--no-tree`, `--prep-only`, `--out-prep-ref`, `--out-prep-seqs`, `--aligned`, `--out-id-map` |
| `msa.py:175-335` `run` | split into `prepare_rows` (`:222-276`, using `orientation` `:91`, `sanitize_id` `:50`, `real_bases` `:64`), `run_mafft` (`:278-293`), `check_alignment` (`:295-306`) and `build_tree` (`:312-328`); call `find_bin` (`:185-190`) per stage |
| `pipe_reconstruct.py:181`, `pipe_from_kmindex.py:235` | `add_argument("--skip-msa", "--skip-tree", dest="skip_tree")` with correct help ("Stop after the walker; the sieve, MSA and tree are skipped", PC-22); add `--no-tree` (MSA yes, FastTree no: skip the FastTree preflight and remove an earlier `<stem>_fasttree.nwk`) |
| `pipe_from_lexicmap.py` | `--no-tree` only (it has no walker stage to stop after) |
| all three composites | forward `--min-identity` (PC-08) |
| `src/disassembler/stage.py` `stage_args` | mode-dependent requirements move into `run()`, because `stage_args` reads argparse's `required` |

**Galaxy pipeline after (d):**

```text
disassembler msa --prep-only        (disassembler container; needs no MAFFT, so the
                                     recipe's unpinned mafft stops mattering here)
  -> rbc_mafft_add 7.526+galaxy2    (ToolShed rnateam/mafft; "Fragments"; keep length on)
  -> disassembler msa --aligned … --no-tree   (length and identity gates)
  -> fasttree 2.1.10+galaxy3        (ToolShed iuc/fasttree; set Nucleotide and GTR)
```

The ToolShed MAFFT tool passes neither `--auto` nor `--thread`. On the 34-row HIV-1 pol fixture with MAFFT 7.505, three outputs were byte-identical: `msa`'s own, `mafft --auto --keeplength --thread 2 --addfragments`, and the wrapper's command line. This was not checked with 7.526. IUC FastTree is 2.1.10, not the pinned 2.2.0, so a workflow tree need not equal the CLI's.

**Staged order.**
1. The help-text fix (trivial).
2. `--no-tree`.
3. `--prep-only` and `--aligned`.
4. The composite flags.
5. `--out-id-map`.

All of it is additive and fits 0.1.2.

**Tests.** Five prototype tests pass with real MAFFT and FastTree:
- `--no-tree` with no FastTree on PATH;
- `--prep-only` with neither tool;
- prep, then the Galaxy wrapper's command line, then `--aligned`, gives the same rows as the fused `msa`;
- `--aligned` refuses an alignment made without `--keeplength`, and refuses a noise alignment;
- `reconstruct --no-tree` writes the MSA and no tree, and `--skip-tree` still writes neither.

They go into `tests/test_msa.py` (stubs), `tests/test_msa_real_tools.py` (split equals fused) and `tests/test_composites.py`. The existing suite passes against the prototype except `TestStageArgs.test_a_missing_required_option_is_refused` (`tests/test_composites.py:95-99`), which pins the old required set and must change.

**Docs.**
- `docs/codon_msa.md`: a modes table and "Running it".
- `docs/galaxy.md`: the tools table, workflow sketches, Open Questions (delete the item), the MAFFT owner (`rnateam`), and line 88.
- `README.md`, `msa -h` and the three composites' help.
- Plan section 2: the codon-align and mafft rows.
- PR #8800's msa help.

**Compatibility.** `msa -a A -t T` behaves exactly as today. `msa -a A` alone still exits 2, now with a message naming `--no-tree`. `summary.json` is unaffected. The existing Galaxy `disassembler_msa` tool keeps working, and new tools are optional.

**Effort.** Medium.

**Risks.**
- The ToolShed MAFFT omits `--auto` and `--thread`. That made no difference on the tested cohort, but it is unverified with 7.526.
- FastTree versions differ between Galaxy and the CLI.
- The id round trip needs `--out-id-map`.
- The `--aligned` checks are the only guard in a split workflow and must not be weakened.

Never use `--out-tree /dev/null` as the interim "discard the tree": run as root, `msa` replaces the device node (PM-09).

## 5. Recommended release plan

### (i) Bioconda build-number hotfix (0.1.1, build 1; no source change)

1. On the #71068 branch, change the disassembler recipe:
   - `number: 1`;
   - `logan-walker >=0.1.1,<0.2` (PK-01);
   - `mafft >=7.526` (PK-02, PK-11);
   - run deps = biopython, requests and the binaries, dropping `zstandard`, `numpy`, `pandas` (PK-06);
   - a `disassembler msa` smoke test on a three-record fixture, with `test: source_files:`.

   Closes backlog 0.1, 0.2 and 0.3.
2. logan-walker recipe: `number: 1`. Keep `skip: True  # [osx]` and give the GraphAligner protobuf/Python reason (PK-05). Re-enable osx only on kmer-diversity-sieve.
3. `extra: autobump: enable: false` on all three recipes (PK-04; replaces backlog 0.5's "only defence").
4. Sieve recipe: leave the source at `70022da` (PK-08; drops backlog (e) option 1).
5. After the merge:
   - check that `conda create --dry-run -c conda-forge -c bioconda disassembler=0.1.1` shows logan-walker 0.1.1 and mafft 7.526;
   - read `conda-meta` in `biocontainers/disassembler:0.1.1--pyhdfd78af_1`;
   - move all pins to `_1` and list `_0` as do-not-use.

   Closes backlog 0.4, which is answered.
6. Repository docs, as soon as possible: `docs/packaging.md:39-41` and `environment.yml:3` (PK-02, PK-09); `docs/packaging.md:27-30` (PK-05); the "Making a release" checklist (PK-04).

### (ii) 0.1.2: disassembler 0.1.2 and logan-walker 0.1.2 (fixes and additive flags; no new walker output files)

Each group's fixes ship together; the order is by risk.

1. **Release gates first.**
   - Run the golden gate on PRs and tags (TC-01) and add a default-flags golden run (TC-02).
   - Add a test with realistic alignment windows (TW-02); this is what would have caught WV-01.
   - Stop tests unlinking `/dev/null` (TW-01, high).
   - Add a CLI test that reaches a successful reconstruction (TW-03).
2. **Walker results.**
   - WV-01 (critical: majority correction inert);
   - WV-02 (1-10% alleles neither corrected nor reported);
   - WV-03 (`ka:f` clipped at 25,000);
   - WV-04, WV-05, WV-06, WV-07, WV-08, WV-09, WV-10;
   - WA-01, WA-02, WA-03, WA-04, WA-05;
   - WG-02 (`km:f` alias), WG-03, WG-04, WG-05, WG-07;
   - WP-01, WP-02 (warn on unrelated `-q` references), WP-04.
3. **Walker robustness.**
   - WF-01 (high: live scratch deleted across PID namespaces), WF-02 (high, 0.1.1 regression), WF-04 (data loss);
   - WF-03 (negative markers, part of backlog (c)), WF-06, WG-01, WG-08 (parser ceiling as `too_large`);
   - WP-03 (error count in the closing line and `--max-errors`);
   - PK-03 (`SSL_CERT_FILE`).
4. **Memory, backlog (b) stages 0-1.** WP-20, WG-13, WP-10 step 1, WP-08 step 1, WP-25, WF-08, WF-16. Closes the backlog's "trivial adjacent fix".
5. **Cheap provenance.** The cache-entry identity check and the golden object MD5s (part of WF-05, backlog (a)).
6. **Python results.**
   - Run-time walker floor `MIN_WALKER_VERSION = (0, 1, 2)` (PK-01);
   - PL-01, PL-02, PL-03 (high);
   - PH-01 (high: tie order changes cohorts);
   - PM-01 (high: disarm the sieve gate; replaces backlog (e)'s upstream ask);
   - PM-02 interim: refuse more than 65,536 admitted records until upstream grows the table;
   - PC-01 to PC-08;
   - PH-02, PH-04, PH-05, PH-06, PH-07, PH-08, and GX-07 (`json_vec` harvest is silently empty);
   - PL-04 to PL-09;
   - PM-03 (codon QC warnings), PM-04, PM-06, PM-07, PM-08, PM-17.
7. **Data.** DD-01 (high), DD-02, DD-05.
8. **Backlog (d).** `msa --no-tree`, `--prep-only`, `--aligned` and the composite flags (PM-18, PC-22).
9. **Repository hygiene.**
   - Move the sieve submodule forward to `70022da` (PK-08);
   - Python lint in CI (TC-23) and remove the crate-wide `allow` (WP-23);
   - sieve tests rebuild on a signal death (TC-14).
10. **Recipes for 0.1.2.** `logan-walker >=0.1.2,<0.2`, a new sha256, build 0, and a published-package check (PK-04 procedure).
11. **tools-iuc #8800**, revised against 0.1.2 rather than token-bumped (GX-10): GX-01 (high: Reference CDS dataset overwritten), GX-04, GX-08 (2.7 × cap memory model), GX-05 and GX-09 (cap parameter, admin ceiling), GX-11, GX-12, GX-19, GX-21.

### (iii) 0.2.0

1. **Backlog (a) in full:** `input`, `input_source`, sidecars and `--revalidate-cache` (WF-05, PK-12). In the same release: walker commit id (PK-07), reference hashes (WV-20), Python tool versions (PM-16), GraphAligner build (TC-27).
2. **Backlog (c) in full:** `status.tsv`, `counts`, `--status-tsv` in fetch, all-skipped wording (WP-14, WP-03 summary half).
3. **Backlog (b) stage 2** (two passes, subgraph-local; `--legacy-parse` for one release; gated on GFA identity), and stage 3 if ready (PK-13). Also the streaming result sink (WP-09) and `--max-memory-mb` (WP-10 steps 2-3).
4. **Recipe floors:** `logan-walker >=0.2,<0.3`. A sieve version floor once upstream tags a release with the table growth (PM-02).
5. **Galaxy:** a `status.tsv` output and a retry branch; the wrapper passes `--max-memory-mb`; the plan's UDT rows are corrected (GX-13, GX-14, GX-15, GX-26).

### Backlog items to drop or reword

- **Drop:** 0.4, which is answered (replace it with PK-01's do-not-use note). In (e): "upstream rewrote history", option 1 (move the recipe), and the dedup-in-UDT interim. In section 4: "subgraph scale figures unchanged". In section 5: "check whether it went up". In the already-fixed list: "The 0.1.1 package itself is correct".
- **Reword:**
  - 0 headline, 0.2, 0.3, 0.5 (B1, B3-B8);
  - (a) scope and effort (B9, B10);
  - (b) rationale, figures, payoffs and staging (B11-B17);
  - (c) to a full `status.tsv`, plus WP-03 (B18);
  - (e) as a downstream fix plus a new upstream list (B20);
  - section 5's sieve-test fix, linter estimate and PR #8800 bump (B22).
- **Move:** (d) and (b) stages 0-1 from 0.2.0 to 0.1.2.
