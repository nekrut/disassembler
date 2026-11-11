# Galaxy delivery: assessment of the plan and recommended route

This section assesses how disassembler should reach Galaxy users. It covers the "Galaxy delivery plan" (2026-11-04; sections 2 to 4: user-defined tools, Format-2 workflows, caveats), tools-iuc PR #8800 at head b9a16ad, the 0.2 backlog where it touches Galaxy, and `docs/galaxy.md`, all at disassembler 0.1.1 (25c27f2). Findings are cited by id; the detail behind GX-NN is in the Galaxy wrappers section. External facts (quay.io tags, the Tool Shed, usegalaxy.org) were last checked on 2026-11-10 and 2026-11-11. "UDT" means a Galaxy user-defined tool (`class: GalaxyUserTool`).

## 1. Summary

- **The direction is sound.** UDTs plus Format-2 workflows in iwc-lab as an interim home, with PR #8800 kept as the route to the Tool Shed, is feasible. iwc-lab exists for exactly this: it accepts `.gxwf.yml`, has CI, and installs nothing on usegalaxy.* servers. Every container the route needs is on quay.io and in the CVMFS Singularity cache, which is the only image source usegalaxy.org's UDT destination accepts.
- **It is a lab route, not public availability.** The plan says this correctly. UDTs need `enable_beta_tool_formats` and a role of type Custom Tool Execution, and each user who imports a workflow gets a private copy of its tools.
- **Three UDT rows cannot be built as written, and one is unnecessary.** The walker cannot reuse the IUC memory shell; if ported, it cuts the size cap to 497 MB on usegalaxy.org's default 1-core, 3.7 GB UDT job (GX-14). The bare-binary sieve returns 0 or 1 records with exit 0 (GX-13). The fetch UDT calls `disassembler` in an image that does not contain it, and is not needed anyway (GX-26). No mafft UDT is needed: a MAFFT 7.526 image and Tool Shed tool already exist (GX-23).
- **The workflow sketches do not run as drawn:**
  - Steps referenced by UDT id cannot be resolved; definitions must be embedded with `run:` (GX-15).
  - WF1 has no reference input and no target selection (GX-27).
  - No built-in Galaxy step turns an accession list into a collection (GX-24).
  - WF2 assembled from the atomic commands keeps 11 rows where the composite keeps 21 (GX-02).
  - WF3-frozen reports skipped accessions as errors (GX-06).
- **Containers.** The plan's 0.1.0 tags are stale. The only 0.1.1 `disassembler` image ships logan-walker 0.1.0 (PK-01), so walker steps must use `logan-walker:0.1.1--hfa8f182_0`. MAFFT 7.526 is inside both `disassembler` images.
- **PR #8800 is open and still pins 0.1.0.** It cannot be merged as it is: GX-01 overwrites a reference dataset. It cannot be bumped by changing tokens either: at 0.1.1 three tests fail, and GX-04, GX-08, GX-02 and GX-20 would stay wrong with every test green (GX-10). Three of the four open review comments are right on the facts. The fourth asks to understand the size cap; the answer is a job parameter now and walker changes later.
- **What must happen first.** Probe usegalaxy.org for three things: egress to the Logan bucket, element identifiers inside a job, and a 4-core / 24 GiB resource request. Rebuild the Bioconda recipes (PK-01, PK-04). Bump the PR to 0.1.1 with the checklist in section 2.
- **CLI changes come before the workflows:**
  - `msa --no-tree` (backlog (d), PM-18);
  - walker `--out-summary` (WP-05) and `--memory-mb` (GX-08);
  - a per-accession status table that also covers `--fetch-only` (backlog (c), GX-06);
  - content-hash provenance (backlog (a), PK-12);
  - sieve fixes PM-01 and GX-03, and `--cohort-gate` (GX-02).
- **Publish workflows only against a release that also fixes the critical and high findings on their path:** WV-01 (critical), WV-02, WV-03, PL-01, PL-02, PM-01, PM-02, PH-01. A workflow runs at default settings, so it would hand these defects to every user.

## 2. tools-iuc PR #8800

### State

- **PR status.** "Disassembler" by as042 is open and not merged. It is labelled ready-for-review, has reviewers bgruening and d-callan, and adds 1,339 lines across 19 files.
- **Commits.** There are nine commits. The last two, d7c9d02 ("ftype, assertions, optional outputs, and other misc changes") and b9a16ad ("Filename fix and memory fix"), answer bgruening's reviews of 27 October and 1 November. CI ("Galaxy Tool Linting and Tests") passed on b9a16ad with 13 jobs.
- **Pins.** `macros.xml:2-3` pins `@TOOL_VERSION@` and `@WALKER_VERSION@` to 0.1.0. The msa tool also pins `mafft 7.526` and `fasttree 2.2.0`, so the Galaxy tool does not inherit the unpinned MAFFT of the published recipe.
- **Open review.** d-callan's review of 4 November is open: a general comment and four inline comments, with no replies and no suggestion blocks.
- **Repository docs.** `docs/galaxy.md:46-50` already describes 0.1.1 behaviour, which the PR does not install (GX-10).

What the capture of the discussion has, and what it lacks. The open comments by the second reviewer (d-callan) are captured: short quotes, the rest paraphrased, and the walker comment is cut off mid-sentence. The first reviewer's (bgruening's) 12 resolved inline threads could not be retrieved; only their anchors are known (macros 1, harvest 3, lexicmap 1, stratify 4, walker 3). Nine collapsed conversations that the page did not load could not be retrieved either. So this review cannot say whether a resolved thread already raised any defect listed below.

### The open comments and the answers

**General comment** ("I worry we're mixing concerns here a bit"; "I might recommend UDTs and IWC-lab as an interim home so you aren't blocked"). It proposes a sync call and an interim home. It does not say what should happen to the PR. The delivery plan adopts the interim home. What remains to settle with the reviewer is whether the PR should continue as six tools, or be narrowed (for example, msa kept fused or removed until the CLI can split it).

**harvest, on the Input help** ("kmindex tool can w v0.7.1 return a single merged file rather than only the per shard collection"). The comment is correct. `kmindex_query 0.7.1+galaxy0` (Tool Shed revision 8:e7bd625dfaa6, 2026-10-10) has a `merge` option that is on by default and writes one `output_merged_json`. The per-sub-index collection is now the non-default mode. The merged document has one top-level key per sub-index, which is the shape `_merge_payload` already iterates (`harvest.py:150-166`). A merged file built from the PR's three test shards gives six outputs byte-identical to the three-shard run.

What harvest should accept:
- the merged JSON (one dataset) or the per-sub-index collection, both typed `json`, through one `type="data" multiple="true"` input;
- not `json_vec`: it is also typed `json`, and today harvest skips every value silently and exits 0 with no targets (GX-07; fix: unwrap `"R"` and fail when nothing numeric was merged);
- not `matrix`.

Harvest is not made redundant by the merge. It still takes the maximum per accession across sub-indexes, drops zeros, ranks, and splits per target. Wrapper changes:
- description and label;
- the Input help (`disassembler_harvest.xml:2`, `:27`, `:115-117`);
- a third test with a merged file;
- regenerated goldens (GX-10, PH-01).

`docs/galaxy.md:17` changes in the same way.

**msa, on "alignment keeps exactly the reference's length"** ("This is a tool using mafft and fasttree? Is it actually a workflow in disguise?"). Largely yes. Both external steps exist on the Tool Shed:
- `rnateam/mafft` "MAFFT add" `rbc_mafft_add 7.526+galaxy2` runs `--addfragments --keeplength`, with neither `--auto` nor `--thread`. Its output was byte-identical to `msa`'s MAFFT step on the fixtures tried (MAFFT 7.505).
- `iuc/fasttree` is `2.1.10+galaxy3`, not the 2.2.0 pinned here. Its model defaults to Protein, so a workflow must set Nucleotide and GTR (GX-23).

The parts specific to disassembler come before and after those two tools:
- before: the A/C/G/T coverage filter, strand orientation by shared 15-mers, header sanitising and de-duplication, and the `Reference` row;
- after: the row-length and identity checks.

The CLI cannot run only those parts today. `--out-tree` is required (`msa.py:160`), and FastTree is looked up before any work (PM-18; backlog (d)). For tools-iuc the split needs two modes: `msa --prep-only` and a check mode `msa --aligned ALN --no-tree`; a prototype gave the same rows as the fused run. For a UDT one `--no-tree` flag is enough. Until then the tool stays fused, and its help should say so.

**sieve, on "This is a lineage-level reducer"** ("Maybe what we want is to wrap this tool then?"). Wrapping `kmer-diversity-sieve` directly is possible (Bioconda 0.1.0 for linux-64 and linux-aarch64; `biocontainers/kmer-diversity-sieve:0.1.0--hab16a5f_0`), but premature:
- The binary's defaults are for whole SARS-CoV-2 genomes (`min_len=27000`, threshold 0.15, sketch 64), and it exits 0 on an empty result. On the PR's 33-record fixture the binary at defaults keeps 0 records, and with only `--min-len 2602` it keeps 1. With all five gene-scale options and its rank shards concatenated, it matches `disassembler sieve` byte for byte at 28 records (GX-13).
- Upstream has no tag; Bioconda's 0.1.0 is commit 70022da.
- The duplicate pre-pass in `cmd/sieve.py` works around an upstream statistics defect. Even with it, `disassembler sieve` leaves the binary's admission gate armed, and divergent records after record 2,000 are dropped with exit 0 (PM-01). A direct wrapper must pass `--warmup 2147483647`.
- The binary hangs once more than 65,536 divergent records fill its table (PM-02). It also groups records by a month parsed from `|YYYY-MM` headers (PM-04).

Recommendation: keep `disassembler_sieve`, fix GX-02, GX-03 and PM-01, and state in the help which `l1_pipeline` options it sets. Add a `kmer_diversity_sieve` tool when upstream tags a release with duplicate-aware admission statistics and the PM-02 fix.

**walker, on the per-accession size cap** ("I would like to understand this better, see what if anything else we can come up w…").

What the cap is: `--max-size-mb` is compared with the object's Content-Length before download. It protects memory. Walker 0.1.1 holds compressed objects whole, up to `4 + F + W` of them at once (F fetch workers, W compute workers), and each parsed graph peaks at about 1.7 times the compressed size. Its own bound is `cap × (4 + F + 2.7 W)`.

What is wrong in the wrapper:
- The memory block counts 43 MB per graph, `cap × (4 + F + W) + 43 × W`, and so under-asks by 1,698 MB per compute worker at the default cap. On 4 slots / 16 GB it keeps 4 + 4 workers, whose bound is 18.8 GB (GX-08).
- Below about 7.7 GB of job memory it lowers the cap silently, to 265 MB at 2 GB and 538 MB at 4 GB, so the cohort depends on the destination (GX-05).
- The admin variable that is meant to change the cap cannot be overridden from `job_conf` (GX-09).
- A skipped accession shows only as a missing record.

What it costs. On the repository's 150 benchmark accessions (median 41.1 MB) the counts over each cap are:

| cap | accessions over the cap (of 150) | pfcrt_50 |
|---|---|---|
| 1024 MB | 5 | 2 |
| 538 MB | 8 | 5 |
| 265 MB | 21 | 17 |

Raising the cap cannot rescue objects over 2^32 bases, about 2.1 to 2.4 GB compressed (WG-08). pfcrt_50's largest object is 10.5 GB.

What can be done, cheapest first:
1. **Wrapper, now.** Make the cap an advanced integer parameter (default 1024), so it is recorded with the job, and never lower it silently. If one worker of each kind does not fit, fail and name the memory needed: 9,856 MB at 1024 MB.
2. **Wrapper, at 0.1.1.** Say in the help that skipped accessions appear as `too_large` in the summary.
3. **Walker, small.** Set the fetch queue depth to 0 (`config.rs:472`, WP-20). This removes four buffers: 304 to 174 MiB measured, wall time unchanged.
4. **Walker, small.** Add `--memory-mb N`, so the walker fits workers to a budget itself and no wrapper copies the model (GX-08).
5. **Walker, backlog (c).** A status table, so a workflow can retry `too_large` accessions on a larger destination.
6. **Walker, backlog (b).** Stream-decode. It turns the cap into a transfer limit only at its last stage (two decodes plus adjacency spooled to disk; PK-13).

### Wrapper defects that matter

| id | sev. | wrapper | defect | present at the 0.1.0 pin? |
|---|---|---|---|---|
| GX-01 | high | walker | a reference dataset named `reconstructed(.fasta)` is truncated through the symlink and overwritten | yes |
| GX-08 | medium | walker | memory fit counts 43 MB per graph; chosen settings exceed the allocation | yes |
| GX-05 | medium | walker | cap lowered silently to fit memory; cohort depends on the destination | yes |
| GX-09 | medium | walker | `DISASSEMBLER_WALKER_MAX_SIZE_MB` cannot be overridden from job_conf | yes |
| GX-12 | medium | macros | `@WALKER_VERSION@` defeats `planemo autoupdate`; first autoupdate installs walker 0.1.1 under a 0.1.0 tool version | yes |
| GX-03 | medium | sieve | default 95% length gate drops walker primaries at 85-95% coverage | yes |
| GX-02 | medium | sieve, msa | LexicMap rows judged by walker-tuned gates; 10 of 20 dropped | yes (gate hard-coded at 0.1.0) |
| GX-10 | medium | all | bump fails harvest test 1, stratify test 1, walker test 3 | at the bump |
| GX-04 | medium | walker | forces `--primary-metric cov`; cannot select 0.1.1's `matched` | at the bump |
| GX-11 | medium | stratify | target taken from the dataset name; 0.1.1 reads it from the file | at the bump |
| GX-16, WP-05 | low, medium | walker | summary `mv` skipped on failure; a `/` in the first record id nests the summary path and the glob fails | yes |
| GX-19, GX-20, GX-21, GX-22, GX-18 | low | msa, lexicmap, stratify | labels and help describe 0.1.0 behaviour or are wrong | at the bump |
| GX-07, PH-01 | medium, high | harvest | `json_vec` gives an empty green harvest; tied scores keep file-name order, so test goldens and cohorts depend on staging names | yes |

No shell-injection path was found: every user string that reaches the shell is regex-validated or reduced to `[\w.-]`.

### Checklist: the 0.1.1 bump (one commit; none of it is a token change)

| # | file | change | finding |
|---|---|---|---|
| 1 | macros | one `@TOOL_VERSION@` per package, in `disassembler_macros.xml` and `walker_macros.xml`; both set to 0.1.1 | GX-12, GX-10 |
| 2 | walker | `mkdir ref` and link or decompress the reference to `ref/<name>.fasta` | GX-01 |
| 3 | walker | `matched` option, selected; `<assert_command>` for `--primary-metric matched` | GX-04 |
| 4 | walker | memory block for `cap × (4 + F + 2.7 W)`; cap as an advanced integer parameter; fail instead of lowering, unless an explicit "fit cap" boolean is set | GX-08, GX-05 |
| 5 | walker | delete `<environment_variables>`; treat `${DISASSEMBLER_WALKER_MAX_SIZE_MB:-}` as an optional admin ceiling | GX-09 |
| 6 | walker | keep the exit status, then collect the summary with `find ... -name '*_walker_summary.json'` even on failure | GX-16, WP-05 |
| 7 | walker | test 3 asserts `no_logan_data` (not `HTTP 404`); test 4 asserts exit 1 and "none errored"; fix the `--allow-empty` help and the `/dev/shm` comment; add a test whose reference is named `reconstructed.fasta` | GX-10, GX-01 |
| 8 | harvest | regenerate both `*.scores.json`; describe the merged input; add a merged-file test | GX-10, GX-07 |
| 9 | stratify | drop the staging; pass `--scores` per dataset; regenerate test data; fix the help (target read from the file; 1.01 no longer needed; sampling by rank) | GX-11, GX-22 |
| 10 | sieve | add `--max-gap-frac` (default 0.05, help: `1 − gate + 0.05` after lexicmap); warn that the 95% gate drops 85-95% primaries until `--min-len-frac` exists | GX-02, GX-03 |
| 11 | msa | add `--min-identity`; relabel `--msa-min-cov` (A/C/G/T over reference length); correct the LexicMap note | GX-19, GX-21 |
| 12 | lexicmap | correct the `--allow-frameshift-alts` text and add the frameshift tier; partial-tier wording under *Require 100%*; accept `tabular,tabular.gz` | GX-20, GX-18 |
| 13 | all | stderr assertions on the bare `expect_failure` tests (walker 4, lexicmap 3, sieve 3); optional `<stdio>` describing exit 3 and mapping 137 to `fatal_oom` | |

Walker tests 1 and 2 must be rerun under GraphAligner 1.0.17b; they could not be run in this review.

### Checklist: 0.2

| CLI change (backlog item, finding) | wrapper change |
|---|---|
| (a) per-object provenance in the summary (WF-05, PK-12) | help names the recorded hashes; a test asserts the `input` object |
| (b) stage 0 and `--memory-mb` (GX-08, WP-10, WP-20) | replace the memory block with `--memory-mb $GALAXY_MEMORY_MB`; keep `--workers` from `GALAXY_SLOTS` |
| (c) status table, written on exit 1 and 3 too (WP-14) | new `tabular` output; test on `SRR0000001` |
| `--out-summary FILE` (WP-05) | drop the rename and the glob |
| (d) `msa --prep-only`, `--aligned`, `--no-tree` (PM-18) | split into prep and check tools; workflows use MAFFT add 7.526 and FastTree (Nucleotide, GTR) |
| `--cohort-gate` on sieve and msa; lexicmap records its gates (GX-02, GX-18) | expose `--cohort-gate`; help |
| `sieve --min-len-frac` default 0.85; admission gate off (GX-03, PM-01) | replace the 95% option |
| harvest accepts `json_vec`, fails on non-numeric input (GX-07, PH-04) | label "json or json_vec" |
| optional cached-unitigs collection input (`--unitigs-cache --offline`) with fetch-time negatives honoured (GX-06) | hidden-by-default input; offline tests become possible |
| (e) tagged upstream sieve with duplicate-aware statistics | add `kmer_diversity_sieve`, retire `disassembler_sieve` |

## 3. The plan, row by row

### What a UDT can and cannot do

These facts come from Galaxy source (25.0 to 26.2) and usegalaxy.org's configuration. They drive the tables below.

- **Format.** `class: GalaxyUserTool` with `id`, `name`, `version`, `container` and `shell_command`. Since Galaxy 25.0; usegalaxy.org runs 26.2 rc1.
- **Outputs.** Several data outputs and collections (`discover_datasets` by pattern) are allowed. There is no `<filter>`, so every output is fixed.
- **Not available.** `<environment_variables>`, reference data, data tables, metadata files.
- **Shell.** `$(...)` and `${...}` in `shell_command` are Galaxy expressions. Shell substitutions must be escaped (`\$(`, `\${`) or written as unbraced variables.
- **Tests.** A UDT's `tests:` block is stored but never run. The only executable test is a workflow test.
- **Resources on usegalaxy.org.** A UDT that requests nothing runs with `GALAXY_SLOTS=1` and `GALAXY_MEMORY_MB=3788`. `requirements: [{type: resource, cores_min, ram_min}]` is honoured up to 8 cores and 29 GB.
- **Images.** Only BioContainers already on CVMFS can be used; nothing is pulled.
- **Network.** Galaxy's reference configuration for UDT destinations sets `docker_net: "none"`. usegalaxy.org's UDT pool is described as egress-only (open outbound); this was not confirmed by a run.

### Containers

| image | on quay / CVMFS | contents and consequence |
|---|---|---|
| `disassembler:0.1.0--pyhdfd78af_0` | yes / yes | disassembler 0.1.0, logan-walker 0.1.0, mafft 7.526 (conda-forge), fasttree 2.2.0, GraphAligner 1.0.17b (hdcf5f25_1), sieve 0.1.0. Stale |
| `disassembler:0.1.1--pyhdfd78af_0` | yes / yes | disassembler 0.1.1 with **logan-walker 0.1.0** (PK-01), mafft 7.526 (h4bc722e_0, conda-forge), fasttree 2.2.0, GraphAligner 1.0.17b (hdcf5f25_1), sieve 0.1.0, Python 3.12, plus unused numpy, pandas, zstandard (PK-06). Built 2026-11-02, pushed 2026-11-05. No `_1` build exists yet. Fine for harvest, stratify, lexicmap, sieve and msa; wrong for anything that runs the walker |
| `logan-walker:0.1.1--hfa8f182_0` | yes / yes | logan-walker 0.1.1, GraphAligner 1.0.17b (h21ec9f0_2, with Python 3.9 and protobuf 3.14); **no `disassembler`**. The walker image to use |
| `logan-walker:0.1.0--hfa8f182_0` | yes / yes | stale |
| `kmer-diversity-sieve:0.1.0--hab16a5f_0` | yes / yes | `l1_pipeline_serial` and busybox; **no Python** (6.8 MB) |
| `mafft` | newest `7.525--h031d066_1` | no 7.526 tag |
| `mulled-v2-027f30e88c6da4ab4b0170ff590dbab158b07243:43ff19c348e21674cd51004463f1ea0b2e5f2336-0` | yes / yes | mafft 7.526 plus fasta3 36.3.8; the image behind `rbc_mafft_add 7.526+galaxy2` |
| `fasttree:2.2.0--h7b50bb2_1` | yes / yes | FastTree 2.2.0 |

The two images carry different GraphAligner 1.0.17b builds, and nothing records which one a result came from (TC-27).

### UDT rows

| UDT (plan) | buildable as described? | what breaks | evidence | fix |
|---|---|---|---|---|
| `disassembler-harvest` | yes | stale tag; input described as a list collection only; `UNION.accessions.txt` matches the per-target pattern; `json_vec` and empty collections exit 0 with nothing (GX-07, PH-08); ties keep staging-name order (PH-01) | merged kmindex file gives byte-identical outputs | `shards`: data, `multiple`, json; move the union out before discovery; tag 0.1.1 |
| `disassembler-stratify` | yes, at 0.1.1 | at the plan's 0.1.0 tag the target comes from the file name | 0.1.1 reads the target from the file (`stratify.py:165-170`) | `--scores PATH`; bins as a repeat; one target per job |
| `disassembler-lexicmap` | yes | `--out-dir` names every output after the reference record id, so `from_work_dir` cannot be fixed (and a `/` in the id crashes, PL-10); twelve fixed outputs; gates not recorded (GX-18) | `gene600.clean.msa.fasta` from a test reference | `-o out/lexicmap`; declare all twelve outputs |
| `disassembler-walker` | **no** | the IUC shell has 8 braced and 4 arithmetic expansions that become Galaxy expressions, uses a variable from `<environment_variables>`, and renames the reference with Cheetah. Pasted, it fails at run time; ported, it cuts the cap to 497 MB at 1 core / 3788 MB (11 of 150 bench accessions skipped instead of 5). At the 8-core / 29,696 MB ceiling it keeps 8 + 8 workers, a 33.6 GB bound (116%). Without the block, the walker's own defaults (8 + 8) bound 33.6 GB on a 3.7 GB job. No `<filter>` for optional outputs. No network where `docker_net: none` | GX-14, GX-25 | `logan-walker:0.1.1` image; `resource: cores_min 4, ram_min 24576`; integer `max_size_mb` input; fit with the 0.1.1 model and fail rather than lower; `ref/` directory (GX-01); `find` for the summary (GX-16); `--scratch-dir` in the job directory; `-t 1`; `primary_metric` matched; always pass `--emit-alts --emit-vcf --emit-snv-table`. Durable: walker `--memory-mb` |
| `disassembler-fetch` | **no**, and not needed | runs `disassembler` in the walker image, which has none (GX-26); usegalaxy.org 26.2 has 809 datatypes and none for zstd, so payloads become `binary`; offline runs then turn `no_logan_data` and `too_large` into `error` (GX-06); a failed cache write is reported as cached (WF-06) | recipe run requirements; datatype list; reproduced statuses | drop. Upload by URL (`https://logan-pub.s3.amazonaws.com/u/<ACC>/<ACC>.unitigs.fa.zst`) with the rule-based uploader or `location:` in a job file. Galaxy records the URL, and since `u/` ETags are the objects' MD5 the upload can carry a checksum |
| `disassembler-codon-align` | yes (interim: `msa` whole) | `--mafft-threads` defaults to 8 on a 1-core job; FastTree always runs; `--out-tree /dev/null` would replace the device node when run as root (PM-09) | `msa.py:160-162` | `--mafft-threads $GALAXY_SLOTS` (unbraced); write the tree to a file and keep it as an output; after (d), `--no-tree` |
| `kmer-diversity-sieve` | **no** | binary defaults keep 0 of 33 records; `--min-len` alone keeps 1; exit 0 both times; no Python in the image for the duplicate pre-pass; shards to concatenate; PM-01 and PM-02 apply | GX-13 | `disassembler-sieve` in the `disassembler:0.1.1` image, with `--max-gap-frac`, threshold and sketch; after PM-01, admission off |
| `mafft` | not needed | "no 7.526 biocontainer exists" holds only for `biocontainers/mafft`. The mulled 7.526 image and `rbc_mafft_add 7.526+galaxy2` exist. A locally mulled image cannot run on usegalaxy.org's UDT destination | GX-23 | drop the row; MAFFT 7.526 is already inside the `disassembler` images |
| fasttree (Tool Shed) | yes, with caveats | `iuc/fasttree` is `2.1.10+galaxy3` (FastTree 2.1.10, Protein by default); the CLI pins 2.2.0, so the workflow's tree is not the composite's | Tool Shed metadata | keep the tree inside codon-align until (d); afterwards set Nucleotide and GTR and state the version in the README |

### Workflows

| workflow | what works | what breaks | evidence |
|---|---|---|---|
| WF1 `from-kmindex` | harvest → stratify chains at 0.1.1; the IUC `kmindex_query` tool exists | no reference input. harvest emits one list per target, so mapping the walker over it runs every target's cohort against one reference: the PR's two-target shards would stream `pfcrt` (29.6 GiB for `pfcrt_50`) to find no anchors. "Map-over cohort" is really a map-over of targets. No built-in step splits a list into elements. The kmindex tool's default output is one merged file, not a collection. Steps are referenced by UDT id. Tied cohorts depend on staging names | GX-27, GX-24, GX-15, PH-01; `pipe_from_kmindex.py:3-5` (single target) |
| WF2 `from-lexicmap` | lexicmap → codon-align needs no network and is the easiest to test | the composite derives the sieve's `--max-gap-frac` as `1 − gate + 0.05` and msa's `--msa-min-cov` as `min(0.70, gate − 0.05)` (`pipe_from_lexicmap.py:339`, `:362`); the atomics use 0.05 and 0.70. Expanded cohort: 11 rows against 21 for the composite, with or without sieve. The plan has no sieve, while `docs/galaxy.md:59` and the CLI have an optional one. Gates are not recorded. MAFFT can move in-frame deletions in repeats (GX-21) | GX-02, GX-18 |
| WF3 `reconstruct` | walker map-over, then collapse, then codon-align is the right shape | needs a split tool (GX-24). `--allow-empty` does not cover errors at 0.1.1, so `__FILTER_FAILED_DATASETS__` is also needed (GX-17). An optional sieve needs a boolean input, `when:` and `pick_value`. Default UDT job too small (GX-14). No network on a by-the-book UDT destination (GX-25). With N accessions per job a job is green if one accession succeeded, however many errored (WP-03). The sieve drops 85-95% primaries (GX-03) and late divergent lineages (PM-01) | reproduced exit 3 under `--allow-empty` with no network |
| WF3-frozen | an offline walker on staged unitig files is testable without S3 in the job (`SRR8209082` is 1.8 MB) | the fetch UDT cannot run (GX-26). Offline, `no_logan_data` and `too_large` become `error`: a single no-data accession exits 3 where the live run exits 0, and a mixed list exits 1 (GX-06). `--fetch-only` without `--allow-empty` exits 1, so Galaxy discards what was fetched. Quota is understated (GX-28): `pfcrt_50` totals 29.59 GiB (9.09 GiB under the 1024 MB cap), `hiv1_50` 7.04 / 3.13 GiB, `sarscov2_50` 1.15 GiB, against the plan's 3.5 GB. Offline runs record no input identity (PK-12) | reproduced statuses; HEAD sizes 2026-11-10 |

## 4. Contradictions

| topic | what the sources say | what is right |
|---|---|---|
| delivery route | plan: UDTs and iwc-lab now, tools-iuc later. `docs/galaxy.md`: six tools-iuc tools, UDTs not mentioned | not a conflict: the plan front-runs the PR. `docs/galaxy.md` should gain an "Interim delivery" section |
| fetch | plan: a UDT and a frozen variant. `docs/galaxy.md:26-28`: not a tool, because of quota | neither a fetch tool nor nothing: upload by URL for the frozen entry. The quota point in `docs/galaxy.md` stands |
| scatter | plan: "native collection ops". `docs/galaxy.md:71`, `:79`: "Galaxy's Split file", "IUC Collapse Collection" | a tool is needed. Name the real tools: `bgruening/split_file_to_collection` (chunk mode, keep order) and `nml/collapse_collections` (GX-24) |
| MAFFT | plan: no 7.526 container; run inside the disassembler image or mulled-build. Backlog §4: the in-image fallback is unsafe until backlog item 0.4 (MAFFT inside the image) is checked. `docs/galaxy.md:86`: an "IUC" MAFFT tool | checked: both `disassembler` images hold MAFFT 7.526. The Tool Shed tool is `rnateam/mafft`. A 7.526 image exists (GX-23) |
| `msa` and `--skip-tree` | `docs/galaxy.md:88`: the split turns `--skip-tree` into "not connecting FastTree" | wrong: composite `--skip-tree` also skips the MSA (and the sieve, despite its help; PC-22, PM-18). Backlog (d) matches the code |
| sieve | plan and backlog (e): bare `l1_pipeline_serial` UDT with an embedded duplicate pre-pass. `docs/galaxy.md`: `disassembler_sieve` wraps `sieve` | `docs/galaxy.md` (GX-13). The pre-pass is not the whole story: PM-01 |
| sieve on the LexicMap route | plan WF2: none. `docs/galaxy.md:59`: optional. CLI: optional with derived gates | keep it optional, with the derived gates (GX-02) |
| harvest input | plan: "shard JSONs (list collection)". `docs/galaxy.md:17`: "the shard collection" | merged file or collection; the kmindex tool defaults to merged |
| object size | plan: median ~71 MB, so 3.5 GB per 50. `docs/galaxy.md:27`: median about 50 MB. `docs/kmindex.md:82`: median 0.6 MB. Backlog §4: 71 MB "still correct" | the medians describe different populations (`kmindex.md` samples random accessions; the others are gene-hit cohorts), so they do not contradict each other. What is wrong is median × N for quota: give totals (GX-28) |
| `--allow-empty` | `docs/galaxy.md:46-50` describes 0.1.1 (all-errored exits 3); `:80` offers `--allow-empty` as an alternative to filtering failed jobs | `:46-50` is right for 0.1.1 but not for the 0.1.0 the PR installs (GX-10). `:80` is wrong at 0.1.1 (GX-17) |
| wrapper memory | `docs/galaxy.md:40-41`: trims `--fetch-workers` so buffers fit. Wrapper: also trims compute workers and the cap, at 43 MB per graph | the 0.1.1 model is `cap × (4 + F + 2.7 W)` (GX-08). Lowering the cap changes the cohort (GX-05) |
| "thin wrapper, so the two stay in sync" | `docs/galaxy.md:8-9`; `docs/packaging.md:66-67`: a bump is "update the version tokens" | defaults and test data drift (GX-04, GX-10) |
| importing UDT workflows | plan §4: workflows "import cleanly but refuse at invoke" | bare UDT ids do not resolve at all. Embedded definitions are refused at import for users without the role (GX-15) |
| plan 1(b) memory | plan: ~700 B per unitig; reduce to the subgraph "during the pass". Backlog: ~34 B, and the 1.7× model "makes the argument more strongly"; Logan objects are "flat unindexed .zst" | CSR is 33.5 to 38.2 B per unitig (PK-15). At defaults 20.0 of the 33.6 GB bound is compressed buffers, which the cheap stages remove. One pass cannot give O(subgraph) (PK-13). The objects are seekable multi-frame zstd (PK-14) |
| plan 1(a) provenance | "~10 lines"; the manifest gives "release-level provenance for free" | cache and offline runs have no response headers. The manifests are 5.1 and 11.8 GB SRA-metadata tables that pin no `u/` bytes. Record a content hash of the parsed bytes plus a cache sidecar (PK-12, WF-05) |
| plan 1(c) status file | "trivial"; retry by raising `--max-size-mb` | size is not a field; `--fetch-only` writes nothing; the wrapper has no cap parameter; objects over 2^32 bases cannot be rescued (WG-08) |
| plan 1(f), backlog item 0 | plan: hygiene in "#70419, merged". Backlog: fixes in #71068 (open); "a stock install has a broken msa"; check the MAFFT in the image; "the only defence" against autobump is to open the PR first | the backlog is right about #71068. "Broken msa" holds only for bioconda-first channels with flexible priority (PK-11). The image has MAFFT 7.526 but walker 0.1.0 (PK-01). `extra: autobump: enable: false` is a defence (PK-04) |
| backlog (e) submodule | "70022da is no longer an ancestor of upstream main" | wrong: upstream main is 70022da and 39d90c0 is its parent. Moving the recipe to 39d90c0 goes backwards (PK-08) |
| LexicMap "msa is a no-op" note | backlog §5: scoped locally, perhaps never pushed | it is in the PR head, scoped to *Require 100%*, and still overstated (GX-21) |

## 5. Recommended delivery

### Tools

All tools are kept as `udt/*.gxtool.yml` in the repository and embedded in the workflows. Use exact image tags; move to `_1` builds once they are published and on CVMFS.

| UDT | container | inputs | outputs | depends on |
|---|---|---|---|---|
| `disassembler-harvest` | `disassembler:0.1.1--pyhdfd78af_0` | `shards`: data, multiple, json (merged file or collection) | `accessions` (list, `(?P<designation>.+)\.accessions\.txt`), `scores` (list), `union`, `summary` | GX-07, PH-01, PH-08 before publishing |
| `disassembler-stratify` | same | `scores`: one json; `bins`: repeat of lo, hi, n | `cohort` (txt), `bins` (json) | |
| `disassembler-split-accessions` | same, or `bgruening/split_file_to_collection` 0.5.2 | `accessions`: txt; `chunk`: integer, default 1 | `chunks` (list; the element name is the accession when chunk is 1) | |
| `disassembler-walker` | `logan-walker:0.1.1--hfa8f182_0`; `resource: cores_min 4, ram_min 24576` | `accessions`: txt; `query`: fasta; optional `cache`: one `.unitigs.fa.zst` (frozen entry); `max_size_mb` (1024); the IUC tuning options with `primary_metric` matched | `reconstructed`, `alt_haplotypes`, `vcfs` (list), `snv_table`, `bubbles`, `summary`; `status` after (c) | `--out-summary`, `--memory-mb`, (c), (a), GX-06 |
| `disassembler-lexicmap` | `disassembler` image | `table`: tabular; `reference`: fasta; gates | the twelve `-o out/lexicmap` files | gates in summary (GX-18) |
| `disassembler-sieve` | `disassembler` image | `input`, `reference`: fasta; `max_gap_frac`, `threshold`, `sketch` (then `min_len_frac`, `cohort_gate`) | `sieved` | PM-01, GX-03, GX-02 |
| `disassembler-codon-align` | `disassembler` image | `haplotypes`, `reference`: fasta; `msa_min_cov`, `min_identity` (then `cohort_gate`) | `alignment`, and `tree` (FastTree 2.2.0) until (d) | (d) |

Stock steps:
- built-ins: `__EXTRACT_DATASET__`, `__FILTER_FAILED_DATASETS__`, `__FLATTEN__`;
- `iuc/pick_value`;
- `nml/collapse_collections` (or `cat1`);
- after (d): `iuc/fasttree` 2.1.10+galaxy3 (Nucleotide, GTR);
- optional: `iuc/bcftools_merge` 1.24+galaxy0, since the walker's VCFs carry a contig line and one sample column.

### WF-A `disassembler-reconstruct` (the plan's WF3; also the subworkflow of WF-B)

Inputs: accession list (txt), reference CDS (fasta), chunk size (integer, default 1), run sieve (boolean, default false).

1. `disassembler-split-accessions` → a list of chunks.
2. `disassembler-walker`, mapped over the chunks, with `--allow-empty`.
3. `__FILTER_FAILED_DATASETS__` on `reconstructed`. An errored accession still fails its job (exit 3). With chunks larger than 1, errors inside a green job are visible only in the summary until the status table exists (WP-03).
4. Collapse to one FASTA, in collection order.
5. `disassembler-sieve` with `when:` bound to the boolean.
6. `pick_value`: the first non-null of steps 5 and 4.
7. `disassembler-codon-align` → alignment and tree.

Outputs: reconstructed FASTA, alignment, tree, per-chunk summaries, SNV tables; the status table after (c).

**Frozen entry** (a second test and a documented way to run the same workflow, not a separate workflow). The user uploads a list collection of unitig files by URL. The walker maps over that collection with `cache` connected, `-a` set to the element identifier, and `--offline`. So no accession reaches the walker without its file, and GX-06 cannot arise. Galaxy keeps the source URL (and a checksum, if one is given at upload), which the walker does not record for cached input (PK-12).

### WF-B `disassembler-from-kmindex` (WF1)

Inputs: the kmindex query result (json, one or many), target name (text), the reference CDS for that target (fasta), score bins (text), plus WF-A's options.

1. `disassembler-harvest`.
2. `__EXTRACT_DATASET__` on `scores` by element identifier. The identifier is harvest's sanitised target name, which can differ from the query name (PH-07).
3. `disassembler-stratify` → cohort list.
4. WF-A as a subworkflow.

### WF-C `disassembler-from-lexicmap` (WF2)

Inputs: LexicMap table, reference CDS, cohort (clean or expanded), run sieve (boolean).

1. `disassembler-lexicmap`.
2. `disassembler-sieve` (conditional) with `--cohort-gate`, then `pick_value`.
3. `disassembler-codon-align` with `--cohort-gate`.

Until `--cohort-gate` exists: expose `max_gap_frac` and `msa_min_cov` as workflow inputs, write the two formulas in the README, and do not offer the expanded cohort.

### Packaging for iwc-lab

- Generate each `.gxwf.yml` from `udt/*.gxtool.yml` with the definitions inlined as `run: {class: GalaxyUserTool, ...}`, using a short script, and commit the generated files. gxformat2 supports this since 0.25.0, and Galaxy 26.1 and 26.2 pin 0.28.0.
- Each workflow directory needs `<name>.gxwf-tests.yml`, `.dockstore.yml`, `README.md`, `CHANGELOG.md`, `license`, `release` and `creator`. Test files over 100 KB go to Zenodo.
- CI runs `planemo-ci-action` against a planemo-managed Galaxy `release_26.1`. Planemo master grants the test user the Custom Tool Execution role and enables beta tool formats. Whether the planemo that the action installs already does this is unverified.
- Walker tests use the frozen entry with `location: https://logan-pub.s3.amazonaws.com/u/SRR8209082/SRR8209082.unitigs.fa.zst`, so no job needs S3 access.

### What has to change in the CLI first

| change | backlog item | findings | needed by |
|---|---|---|---|
| `msa --no-tree`; `--out-tree` optional; FastTree looked up only when used | (d) | PM-18 | codon-align UDT; tools-iuc split also needs `--prep-only` and `--aligned` |
| `logan-walker --out-summary FILE` | new | WP-05, GX-16 | walker UDT and XML |
| `logan-walker --memory-mb N`: fit workers to the budget, never lower the cap, exit 2 with the amount needed; bound counts accessions and mode | (b) stage 0 | GX-08, GX-14, WP-10, WP-08 | walker UDT on small jobs |
| fetch queue depth 0; on glibc, `mallopt(M_MMAP_THRESHOLD, 1 << 20)` (fixes the memory allocator's mmap threshold so freed graph memory is returned) | (b) stage 0 | WP-20, WF-08 | lowers `ram_min` |
| per-accession `status.tsv` (accession, status, class, compressed size, ETag/MD5, detail), written on exit 1 and 3 too, also by `--fetch-only`; `counts` in the summary | (c) | WP-14, WP-03 | retry branch; chunked jobs |
| `--fetch-only` records negatives; `--offline` honours them | (c) | GX-06, WF-03 | frozen entry via CLI; tools-iuc cache input |
| `input` object per accession (URL, size, ETag, Last-Modified, fetch time, MD5 and SHA-256 of parsed bytes) and a cache sidecar | (a) | WF-05, PK-12 | publishing the live-S3 entries |
| sieve: admission gate off (`--warmup 2147483647`), REJECT count checked; `--min-len-frac` (0.85) | (e) is upstream; these are local | PM-01, GX-03 | any workflow with the sieve |
| `--cohort-gate` on sieve and msa; lexicmap writes its gates | new | GX-02, GX-18 | WF-C |
| harvest: tie-break by accession; `json_vec`; fail on empty collection or non-numeric input | new | PH-01, GX-07, PH-04, PH-08 | WF-B |
| walker correctness on the default path | | WV-01 (critical), WV-02, WV-03; PL-01, PL-02 for WF-C; PM-02 for the sieve | publishing any workflow |

Release with the recipe procedure fixed first (PK-04: autobump off, floors set, image contents checked after publishing).

### What stays in tools-iuc

PR #8800, bumped to 0.1.1 with the checklist in section 2, is the route to public availability. It keeps six tools: `disassembler_harvest`, `_stratify`, `_walker`, `_lexicmap`, `_sieve` and `_msa`.
- `disassembler_sieve` stays until upstream tags a duplicate-aware release.
- `disassembler_msa` stays fused until `--prep-only` and `--aligned` exist; then it splits into a prep tool and a check tool.
- The walker gains the status output and an optional cached-unitigs input at 0.2.

The UDTs use the same command lines. Once the Tool Shed tools are installed on usegalaxy.*, each workflow swaps its `run:` blocks for Tool Shed ids and can move to IWC proper as `.ga`.

### Order of work

| step | work | depends on | unblocks |
|---|---|---|---|
| 1 | Bioconda: disassembler and logan-walker 0.1.1 build 1 (walker floor `>=0.1.1`, `mafft >=7.526`, trimmed dependencies); `extra: autobump: enable: false`; check the new image's contents | none | a correct `disassembler` image (PK-01) |
| 2 | probe usegalaxy.org with three throwaway UDTs: egress to `logan-pub.s3.amazonaws.com`, element identifiers of a list input, a `resource` request of 4 cores / 24576 MiB | none | decides whether the live entry runs there, and how the frozen entry passes accessions |
| 3 | PR #8800: 0.1.1 bump per the checklist; answer the four open comments; agree the PR's scope with the reviewer | none | public route at 0.1.1 |
| 4 | CLI and walker changes in the table above, small ones first: (d), `--out-summary`, PM-01, GX-03, GX-02, harvest items; then `--memory-mb` with stage 0, (c) with GX-06, (a); WV-01, WV-02, WV-03 | none | steps 6 to 8; tools-iuc 0.2 |
| 5 | release 0.2 with the fixed procedure; new images on quay and CVMFS | 1, 4 | final pins |
| 6 | `udt/*.gxtool.yml`, the inlining script, WF-A with a frozen-entry test on `SRR8209082`; can start on 0.1.1 images | 2 | WF-B, WF-C |
| 7 | WF-C (no network), then WF-B | 6; `--cohort-gate` and harvest items from 4 | |
| 8 | iwc-lab PR (planemo version confirmed) | 5, 6, 7 | lab use |
| 9 | tools-iuc 0.2: status output, msa split, cached input; later `kmer_diversity_sieve`; then workflows on Tool Shed ids and into IWC | 3, 5; upstream sieve tag | public workflows |
