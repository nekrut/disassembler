# Plan: port disassembler to Galaxy

Status (2026-10-17): the six Galaxy tool wrappers are drafted and pass their tests locally (section 2). The Bioconda packages (section 1) and the `.ga` workflows (section 3) are not built yet. Section 5 tracks what is left.

Goal: expose the disassembler suite as a set of Galaxy tools and workflows so a user can go from a kmindex planetary query (or a LexicMap tabular search) to a codon MSA and tree without leaving Galaxy. This document is the port's shopping list: the Bioconda packages that have to exist, the Galaxy tools that wrap them, and the knobs each tool exposes.

The port keeps the CLI's interface. Every Galaxy tool is a thin wrapper around the same subcommand a terminal user runs today, so the two surfaces stay in sync as the CLI evolves. Where Galaxy exposed a problem in the CLI, it was fixed in the CLI itself (section 1.7), so terminal users get the fix too. Flags that only make sense outside Galaxy (local caches, S3 fallbacks, API keys, binary path overrides) are not exposed by the wrappers; they stay on the CLI for workstation use.

---

## 1. Bioconda packages

The current tree ships three independently buildable artifacts: a Python package, a Rust binary, and a C binary. They sit in one repo for development reasons but have no shared build system, which makes the one-package-per-artifact split obvious. The master package is `disassembler`; it depends on the other two plus the existing Bioconda tools.

Bioconda builds every recipe from a downloadable tagged tarball, so each source repository must be public and tagged before its recipe can be submitted.

### 1.1 `logan-walker` (new)

The Rust path walker, built from `walker/`. One binary, `logan-walker`.

- Source: `walker/` in `nekrut/disassembler` at a tagged version, with its own tag series (for example `logan-walker-v0.1.0`) so its tags do not collide with the Python package's.
- Build: `cargo install --locked --path walker`. Two crate dependencies (`ureq`, `zstd`), but both compile C underneath (`ring` for TLS, `zstd-sys` for a bundled libzstd), so the recipe needs `{{ compiler('c') }}` and `{{ stdlib('c') }}` alongside `{{ compiler('rust') }}`, plus `cargo-bundle-licenses` to ship the crates' licenses.
- Shipped artifact: the `logan-walker` executable on PATH.
- Runtime dependency: GraphAligner (shelled out via PATH). Expressed as a Bioconda `run` dependency so a user installing `logan-walker` picks up `graphaligner` automatically.
- Version pin: `graphaligner=1.0.20` for now, with a known regression noted in `docs/issue_graphaligner_end_clipping.md`. Switch to `1.0.17b` if the end-clipping fix matters to a given analysis.
- Why its own package: the walker is useful without the rest of disassembler (anyone with a Logan accession list and a reference CDS can use it), and keeping it independent lets Galaxy tool wrappers depend on only the walker when that is all they need.

### 1.2 `kmer-diversity-sieve` (new)

The MinHash/LSH diversity sieve from `extern/kmer-diversity-sieve`, currently a git submodule.

- Source: the upstream repo the submodule tracks. It has no tagged release yet, which a recipe needs.
- Build: `make serial` produces `l1_pipeline_serial`. Single C source file, no system deps beyond a C compiler. The Makefile hard-sets `CFLAGS := -O3 -march=native`, which tunes the binary to the build machine, so it can crash with an illegal instruction elsewhere. The recipe overrides it by passing `CFLAGS` on the `make` command line, unless upstream relaxes it first.
- Shipped artifact: the `l1_pipeline_serial` executable on PATH.
- Why its own package: it was built for a general problem (archive-scale deduplication of near-identical consensus genomes) and has utility outside disassembler. If upstream already has a Bioconda recipe we use theirs; otherwise we contribute one.
- Open question: whether to ship `l1_pipeline_threads` and the MPI build too. The disassembler uses only `_serial`; a general-purpose Bioconda package should probably ship all three but default to serial.

### 1.3 `disassembler` (new, master)

The Python package providing the `disassembler` CLI and every subcommand.

- Source: `src/disassembler/` in `nekrut/disassembler`.
- Build: standard `pyproject.toml` wheel, `noarch: python`. The version is read from `disassembler.__version__`.
- Python runtime deps (already in `pyproject.toml`): `biopython`, `requests`, `zstandard`, `numpy`, `pandas`.
- External-tool `run` dependencies on Bioconda:
  - `logan-walker` (new, section 1.1)
  - `kmer-diversity-sieve` (new, section 1.2)
  - `graphaligner` (existing; transitive via logan-walker, listed here for clarity)
  - `mafft` (existing)
  - `fasttree` (existing)
- Not shipped: Rust and C toolchains. End users never build anything.
- Why "master": installing `disassembler` from Bioconda pulls in everything the CLI expects to find on PATH. A workstation user runs `conda install -c bioconda disassembler` and the full suite works. A Galaxy wrapper that needs only one piece (the walker) depends on the narrower package directly.

### 1.4 Existing Bioconda packages we lean on

Not new, but worth listing so the dependency graph is explicit:

| package | used by |
|---|---|
| `graphaligner` | logan-walker (shells out) |
| `mafft` | disassembler msa |
| `fasttree` | disassembler msa |
| `zstd` | nothing at build or run time: the walker's `zstd` crate compiles its own bundled libzstd. The CLI is only handy for inspecting cached `.zst` files by hand |

No existing Bioconda package covers kmindex or LexicMap search themselves; those live in Galaxy and are run as Galaxy tools, not pulled into a disassembler install.

### 1.5 Dependency graph

```
                         disassembler (python wheel)
                        /       |         |         \
             logan-walker    kmer-diversity-sieve   mafft    fasttree
                 |
            graphaligner
```

Everything the CLI calls out to is reachable from one `conda install disassembler`. The Galaxy tool wrappers in section 2 depend on narrower slices of this graph.

### 1.6 Versioning policy

- Each package versions independently and tags its own releases.
- `disassembler` pins the minor versions of its dependencies, so a `disassembler 0.1.x` run reproduces regardless of whether logan-walker 0.1.0 or 0.1.1 is the latest.
- Golden-fixture concordance (`bench/check_golden.sh`) is the gate for a `disassembler` release.
- `disassembler --version` prints `disassembler.__version__`, which is also the installed package version. `logan-walker --version` prints the version in `walker/Cargo.toml`. The Galaxy wrappers use these as their version commands.
- The first public release is `0.1.0` for both packages (`disassembler.__version__` and `walker/Cargo.toml`). The drafted wrappers require exactly these versions through the tokens in their `macros.xml`, so a different release number means changing those tokens too.

### 1.7 CLI changes the port needed

Made on 2026-10-17. The first tagged releases must include them, because the wrappers depend on them.

- `disassembler --version`. Galaxy runs each tool's version command inside the job script, which runs under `set -e`, so a version command that fails aborts the job before the tool runs.
- File formats are recognised by content, never by file name. Galaxy stages every dataset as `dataset_<uuid>.dat`:
  - `logan-walker` always hands GraphAligner a temporary `.fa` copy of the reference, because GraphAligner picks its parser from the file name and refuses `.dat`. The copy is byte for byte for a single-record reference, so the golden gate is unchanged.
  - `lexicmap` detects gzip from its magic bytes instead of a `.gz` suffix.
  - `harvest --shards-dir` reads every non-hidden file as a shard instead of only `*.json`. A file that is not valid JSON counts as a failed shard and is named in the log.
- `reconstruct --sieve` and `from-kmindex --sieve` no longer crash after the walker stage. Found while testing: the sieve call still used the `query` variable that multi-anchor support had replaced.

---

## 2. Galaxy tools

Six tools cover both arms of the workflow. One existing subcommand (`fetch`) does not become a Galaxy tool; its purpose (local filesystem caching) is subsumed by Galaxy's data staging.

Convention used below: "exposed" means the flag appears in the tool XML as a user-visible param. "Not exposed" means the wrapper sets it from Galaxy context, hardcodes it, or omits it entirely.

The wrappers are drafted in a fork of `galaxyproject/tools-iuc` (`as042/tools-iuc`, branch `disassembler`, directory `tools/disassembler/`), which is where they will be submitted. They share a `macros.xml` and publish as a Tool Shed suite (`suite_disassembler`). Fifteen planemo tests pass on Galaxy 25.1 against the working tree. They have not run against real Bioconda packages, which do not exist yet. The walker tests stream about 2 MB from the public Logan bucket, so they need network access. Every Python tool's version command is `disassembler --version`; the walker's is `logan-walker --version`.

### 2.1 `disassembler_harvest`

Merge kmindex planetary shard JSONs into per-target scored accession lists.

- Depends on: `disassembler`.
- Driven by: `disassembler harvest --shards-dir` (the local-files mode added for the port).

Inputs:
- Shard JSON collection (`data_collection_type="list"`, format `json`). This is the direct output of the IUC `kmindex_query` tool.

Exposed params: none beyond the input collection. Harvest's output shape is deterministic once the inputs are set.

Not exposed:
- `--history-id`, `--collection-id`, `--base-url`, `--api-key-file`: Galaxy API mode does not apply inside Galaxy.
- `--shards-dir`: set by the wrapper, which stages each collection element as `shards/<index>.json`.
- `--download-workers`: ignored in local mode.

Outputs:
- A `list` collection of per-target `<target>.accessions.txt`.
- A `list` collection of per-target `<target>.scores.json`.
- `UNION.accessions.txt` as a single `txt` dataset. The wrapper moves it out of the output directory first, so the per-target collection does not pick it up.
- `summary.json` as a single `json` dataset.

### 2.2 `disassembler_stratify`

Cut a size-capped cohort from a harvest by score bin.

- Depends on: `disassembler`.
- Driven by: `disassembler stratify`.

Inputs:
- Per-target scores collection (from harvest). Element identifiers are the target names. `stratify` reads only the scores, so the accession-list collection is not an input.

Exposed params:
- `--score-bin LO:HI:N`: as a repeat of (LO, HI, N) in the tool XML, so the workflow can request multiple bins. Bins are half-open, so a bin meant to include a perfect containment score of 1.0 needs HI above 1. The wrapper defaults HI to 1.01.

Not exposed:
- `-o/--out-dir`: Galaxy-managed.

Outputs:
- Per-target subsampled accession list collection (`<target>.bins.txt`).
- Per-target bins-JSON collection (`<target>.bins.json`).
- `bins_summary.json`.

### 2.3 `disassembler_walker`

Stream unitigs, build the local subgraph, walk haplotypes, emit variants.

- Depends on: `logan-walker` (which transitively depends on `graphaligner`).
- Driven by: `logan-walker` directly. The wrapper does not need the Python package.

Design decision: the tool takes **an accession list dataset** (or a single accession typed in), and the scatter strategy is chosen by how the list is split before it reaches the tool (section 3.4). An earlier draft of this plan exposed one accession per job as a scalar input, but Galaxy cannot map a text parameter over a collection. A `txt` dataset maps naturally: split the list into one-line files and every accession gets its own job, retry and provenance; split it into N-line files for scatter-per-shard; or pass it whole for one batch job. One tool covers all three strategies.

Inputs:
- Accessions: an accession list (`txt`, one per line, for example a harvest or stratify output), or a single accession (text, restricted to letters, digits, `_` and `-`).
- Reference CDS FASTA (`-q`). More than one record switches on multi-anchor mode. The wrapper links it under the dataset's name, because the walker names its summary file and multi-anchor tags after the file stem, which would otherwise be `dataset_<uuid>`.

Exposed params:
- `--max-size-mb` (default 1024): per-accession size cap, important user policy knob.
- `--min-cov` (default 0.85): reconstruction quality gate.
- `--prune-frac` (default 0.01): reference-bias correction strength.
- `--prefer-majority` / `--no-prefer-majority` (default on): majority-allele preference at biallelic sites.
- `--allow-empty` (default off): exposed so a workflow can decide whether a job that reconstructs nothing is a soft or hard failure.
- Anchoring (collapsed section): `-k`, `--hops`, `--min-shared-kmers`, and `--primary-metric` (`cov`, `pident` or `cov_x_pident`), which picks the primary haplotype in multi-anchor mode.
- Intra-host variants (collapsed section): `--min-abund`, `--min-vaf`, `--min-abund-frac`, `--alt-max-detour`, `--no-indel-alts`, `--allow-frameshift-alts`, `--max-alt-indel`. Most users leave them alone.
- Optional outputs: `--emit-alts` / `--no-emit-alts` (default on), `--emit-vcf` and `--emit-snv-table` (default off). Each controls whether the corresponding output is created.

Not exposed:
- `--unitigs-cache`, `--offline`, `--fetch-only`: Galaxy handles data staging. The walker streams S3 inside its job.
- `--walker-bin`: binary is on PATH via the `logan-walker` package.
- `--scratch-dir`: set by the wrapper to a directory inside the job. The walker's default, `/dev/shm/logan-walker`, is shared by every job on the node and is only 64 MB inside a default Docker container.
- `-o`, `-f`: output paths Galaxy-managed.
- `-w, --workers`: set from `\${GALAXY_SLOTS:-1}`.
- `--fetch-workers`: starts at the slot count and is trimmed so the walker's worst-case resident memory, (4 queued + fetch workers + compute workers) x `--max-size-mb`, fits in `GALAXY_MEMORY_MB` when the job runner sets it. Never below 1.
- `-t, --threads`: left at 1. GraphAligner's output is not reproducible with more than one thread, so it does not follow the slot count.

Outputs:
- `reconstructed.fasta` (always; empty when nothing reconstructed and `--allow-empty` is on).
- `alt_haplotypes.fasta` (filtered on `--emit-alts`).
- A `list` collection of `<acc>.vcf` (filtered on `--emit-vcf`). The walker writes a VCF only for an accession with at least one recovered SNV, so the collection can be empty.
- `snvs.tsv` (filtered on `--emit-snv-table`).
- `bubbles.json`, and `summary.json` (the walker's `<stem>_walker_summary.json`, renamed).

The walker streams from `logan-pub.s3.amazonaws.com` inside the job, so the Galaxy server needs outbound network access.

### 2.4 `disassembler_lexicmap`

Stream a LexicMap tabular dataset, do multi-HSP coordinate tiling and codon QC, emit the cohort ledger and the collapsed clean cohort.

- Depends on: `disassembler`.
- Driven by: `disassembler lexicmap -f`.

Inputs:
- LexicMap tabular dataset (`tabular`), from LexicMap Search run with *Output all columns* (`--all`) on, so `qseq` and `sseq` are present. Gzipped input is detected by content.
- Reference CDS FASTA (`-q`).

Exposed params:
- `--min-cov` (default 0.80): clean full-length cohort gate.
- `--min-cov-partial` (default 0.50): clean partial cohort gate.
- `--min-pident` (default 60.0): HSP identity floor.
- `--max-internal-stops` (default 0): codon QC strictness.
- `--require-full-length` (off): require 100% coverage for the full-length cohort instead of `--min-cov`.
- `--allow-frameshift-alts` (off): whether frameshift candidates reach the clean cohort.
- `--samples-per-haplotype` (default 10): representative accession cap in the haplotypes TSV.
- `--no-a2m` (off): disable the `.a2m` output.
- `--allow-empty` (off): succeed even when the clean cohort is empty.

Not exposed:
- `-d/--dataset`, `--base-url`, `--api-key`, `--api-key-file`: Galaxy API mode does not apply.
- `-o/--out-prefix`: Galaxy-managed.
- `--max-records`: a testing/preview limit, not a user policy.
- `-f/--file` with stdin: the wrapper always points at a staged tabular dataset.

Outputs:
- `clean.msa.fasta`, `clean.haplotypes.tsv`, `clean.accessions.fasta`, `clean_expanded.accessions.fasta`, `flagged.accessions.fasta`, `flagged.tsv`, `cohort_ledger.tsv`, `summary.json`, `clean.a2m`.
- Not exposed: the compatibility aliases `msa.fasta`, `haplotypes.tsv` and `.a2m`, which duplicate the `clean.*` files.

### 2.5 `disassembler_sieve`

k-mer diversity sieve: thin a cohort of near-duplicate reconstructions before the O(N^2) MSA/tree stage.

- Depends on: `disassembler`, which pulls in `kmer-diversity-sieve` (section 1.2). The gene-scale defaults live in the `disassembler sieve` subcommand, so the wrapper inherits them instead of restating them.
- Driven by: `disassembler sieve`.

Inputs:
- FASTA of candidate sequences to thin.
- The length gate, chosen by a conditional: either a reference CDS FASTA (`--ref`), from which the subcommand derives `--min-len` as 95% of the first record, or an explicit `--min-len`.

Exposed params:
- `--threshold` (tau, default 0.005): the gene-scale default, not the sieve binary's own default.
- `--sketch` (M, default 256): resolution floor.
- `--min-len` (in the explicit branch of the conditional): integer length gate. The wrapper never computes it itself.

The composites spell these three knobs `--sieve-threshold`, `--sieve-sketch` and `--sieve-min-len`. The wrapper drives the standalone subcommand, so it uses the short names.

Not exposed:
- `--sieve-bin`: binary is on PATH via the `kmer-diversity-sieve` package.

Outputs:
- Thinned FASTA.

Decision: `disassembler_sieve` ships on Galaxy as part of the disassembler suite, alongside the other five tools, rather than as a separate general-purpose IUC sieve tool. The sieve binary keeps its own Bioconda package (section 1.2), so a standalone IUC tool with the binary's archive-scale defaults can still be added later without affecting this one. Section 4 records the decision.

### 2.6 `disassembler_msa`

Reference-guided in-frame codon MSA (MAFFT `--keeplength`) plus FastTree GTR tree.

- Depends on: `disassembler`, which pulls in `mafft` and `fasttree`.
- Driven by: `disassembler msa`.

Inputs:
- Reconstructed haplotypes FASTA.
- Reference CDS FASTA.

Exposed params:
- `--msa-min-cov` (default 0.70): MSA entry coverage gate.

Not exposed:
- `--mafft-threads`: set from `\${GALAXY_SLOTS:-1}`.
- `-a/--out-aln`, `-t/--out-tree`: Galaxy-managed.

`msa` always writes a tree (`--out-tree` is required). `--skip-tree` exists only on the `from-kmindex` and `reconstruct` composites and has no tool-level equivalent. Under the decomposition recommended below it becomes a workflow choice: connect the FastTree box or not.

Outputs:
- Codon MSA FASTA.
- FastTree Newick.

Drafted as a single tool. The decomposition into "reference-guided codon MSA prep" + IUC `mafft` + IUC `fasttree` that section 4 recommends needs a prep-only mode in the CLI first: `disassembler msa` filters and renames the sequences, runs MAFFT, checks every aligned length against the reference and runs FastTree in one call, so today there is no disassembler-specific step to wrap on its own. Section 4 carries the full note.

### 2.7 Not ported: `disassembler_fetch`

`fetch` populates a local filesystem cache for `--unitigs-cache` / `--offline`. Galaxy's data model subsumes this: tool outputs are managed datasets in the history and are already cached across downstream jobs. Shipping fetch as a tool would require either (a) producing a collection of 51.6 MB median unitig files that charges against user storage quota, (b) a data manager writing into a shared admin-managed cache, which only makes sense on a dedicated Logan Galaxy instance.

Decision: drop `fetch` from the Galaxy port. The CLI keeps it for workstation use.

### 2.8 Not ported: the composites

`reconstruct`, `from-kmindex`, and `from-lexicmap` are the three CLI composites. On Galaxy they become `.ga` workflow files (section 3), not tools.

---

## 3. Workflows

Three `.ga` files, one per CLI composite.

### 3.1 `from-kmindex.ga`

Mirror of `disassembler from-kmindex`. Full planetary-scale triage.

```
kmindex_query (IUC)
  -> disassembler_harvest
  -> disassembler_stratify  (optional, controlled by a workflow conditional)
  -> disassembler_walker    (over the stratified accession list, split per section 3.4)
  -> disassembler_sieve     (optional)
  -> disassembler_msa
```

Workflow-level choices that were CLI flags:
- `--sieve` on the CLI (off by default) becomes "connect the sieve box or wire walker to msa directly" in the workflow editor.
- `--skip-tree` becomes "stop at msa output" versus "connect FastTree".
- `--skip-harvest` does not translate; replaced by workflow entry-point choice (start from an accession list dataset instead of a kmindex collection).

### 3.2 `from-lexicmap.ga`

Mirror of `disassembler from-lexicmap`.

```
lexicmap_search (IUC or custom)
  -> disassembler_lexicmap
  -> disassembler_sieve     (optional)
  -> disassembler_msa
```

Note that `disassembler from-lexicmap` currently does a byte-exact dedup pass before the sieve. On Galaxy this is a one-line Python step; either fold it into `disassembler_lexicmap`'s output path or ship a small `disassembler_dedup` tool. Minor decision.

### 3.3 `reconstruct.ga`

Entry point for users who already have an accession list (from an SRA query, a published panel, a previous harvest cached elsewhere).

```
accession list (user upload or history item)
  -> disassembler_walker  (scattered)
  -> disassembler_sieve   (optional)
  -> disassembler_msa
```

### 3.4 Scatter policy

The walker takes an accession list (section 2.3), so the scatter strategy is chosen by how the list is split before it reaches the walker. Galaxy's existing **Split file** tool does the splitting; no disassembler splitter tool is needed. Three strategies are selectable at the workflow level:

1. **Scatter-per-accession** (default): split into one-line files. Natural Galaxy pattern, per-accession retry and provenance, breaks at >~10k accessions due to scheduler overhead.
2. **Scatter-per-shard**: split into files of N lines; each job runs the walker's own internal list loop. Clean up to >100k accessions.
3. **Single-job batch**: pass the list unsplit, so one walker call processes the whole list internally. Loses per-accession isolation, scales to any size memory and time allow.

For SARS-CoV-2-scale cohorts (10^5-10^6 accessions) only option 2 is realistic. Workflow-level routing selects the strategy.

Downstream of the walker, the per-job reconstruction FASTAs are concatenated with the IUC **Collapse Collection** tool before the sieve or msa. A job that reconstructs nothing fails by default, so with one job per accession either turn on `--allow-empty` or filter failed datasets before collapsing.

---

## 4. Open questions and deferred decisions

Resolving these is part of the plan, not a prerequisite to it. Each question lists its status, any outside owner the answer waits on, and the current recommendation when there is one.

- **Sieve tool boundary.** Section 2.5. Ship `disassembler_sieve` as a thin wrapper with gene-scale defaults baked in, or contribute a general-purpose `kmer_diversity_sieve_l1` IUC tool with upstream defaults and set gene-scale values at the workflow level. *Status: resolved.* `disassembler_sieve` ships in the disassembler suite for now, which keeps the gene-scale defaults in one place and is less for a workflow user to understand. `kmer-diversity-sieve` remains its own Bioconda package (section 1.2), so a general-purpose IUC sieve tool can still be contributed later without conflicting with this one.
- **MSA tool boundary.** Section 2.6. Ship `disassembler_msa` as one tool, or decompose into a small reference-guided codon prep tool + IUC MAFFT + IUC FastTree. One tool is easier to use; decomposed is more inspectable. *Status: open; current recommendation is to decompose.* Rationale: avoid re-wrapping MAFFT and FastTree, which are already well-maintained in the IUC. The disassembler-specific value is the reference-guided in-frame codon preservation layer, and that is what gets its own tool. The decomposition also makes `--skip-tree` disappear naturally (just do not connect the FastTree box in the workflow). It waits on a prep-only mode in the CLI (section 2.6), so the drafted wrapper is a single tool until that exists.
- **Metadata enrichment.** Not currently in the CLI. For SARS-CoV-2 at scale the harvest's 5M candidate list needs an SRA metadata join (date, platform, country) before stratify, so that a stratified cohort is informative, not random. This likely lives in the sibling `sra-metadata-enrichment` project and connects into the workflow, rather than in disassembler. *Status: open, out of scope for this port.*
- **VCF merging downstream of walker.** The CLI's `snvs.tsv` long table is the scalable output; a merged multi-sample VCF is a convenience. On Galaxy the natural shape is `bcftools merge` (IUC) wired to the walker's per-accession VCF outputs. Not a disassembler tool; worth calling out in the workflow docs. *Status: open, workflow-doc concern not a tool concern.*
- **GraphAligner version.** `environment.yml` pins 1.0.20 but 1.0.19/1.0.20 regressed query end-clipping (`docs/issue_graphaligner_end_clipping.md`). The `logan-walker` Bioconda package should default to a known-good version; 1.0.17b for correctness-sensitive work, 1.0.20 for parity with the committed golden fixtures (which were captured under the affected build). The right default depends on whether a Galaxy run is expected to byte-reproduce the committed panels, which it probably is not. *Status: open.* Affects correctness of every Galaxy output, so the resolution has to land before the first published Galaxy tool version; probably forces a goldens rebuild.
- **Logan Galaxy instance vs usegalaxy.org.** The port is drafted assuming usegalaxy.org, which charges user storage and has rate limits. A dedicated Logan Galaxy instance changes two things: the `fetch` tool becomes attractive again as a data manager, and the walker's `--max-size-mb` policy could default to something larger. *Status: open; deferable.* The port works on usegalaxy.org as drafted; a Logan-instance variant can land later as an addendum.
- **Who maintains the Bioconda recipes.** The disassembler side is natural for the disassembler maintainers. `kmer-diversity-sieve` is upstream of us; preferred path is to contribute the recipe upstream and track their releases, not vendor it. *Status: open.*

---

## 5. Order of work

A reasonable sequence; none of the steps block the next past where noted. The questions in section 4 appear here as explicit resolve-this steps so they do not get skipped.

1. **Settle who authors the `kmer-diversity-sieve` Bioconda recipe** (section 4, who maintains the Bioconda recipes). The preferred path is a recipe contributed upstream that tracks upstream releases. **Blocks step 2.**
2. Contribute (or receive) a Bioconda recipe for `kmer-diversity-sieve`. Direction set by step 1. Needs a tagged upstream release first (section 1.2). **Blocks step 5.**
3. Decide the GraphAligner version pin for `logan-walker` (section 4). If 1.0.17b wins, rebuild `bench/golden/` before the port ships so the concordance gate stays live. **Blocks step 4.**
4. `nekrut/disassembler` is public (done 2026-10-18), so Bioconda can download it. Tag `logan-walker 0.1.0` from `walker/`, including the section 1.7 changes. Contribute Bioconda recipe depending on step 3's version pin. **Blocks steps 5 and 6.**
5. Tag `disassembler 0.1.0`, including the section 1.7 changes. Contribute Bioconda recipe depending on steps 2 and 4.
6. **Drafted.** The six tool XMLs are in `as042/tools-iuc`, branch `disassembler`, `tools/disassembler/` (section 2). Remaining: change the requirement versions if the tags from steps 4 and 5 differ from `0.1.0`, add authors to the suite's citation, and re-run `planemo test` with normal dependency resolution once the packages from steps 2, 4 and 5 exist. Decomposing `msa` (section 4) is a later follow-up.
7. Author the three `.ga` workflow files, wiring the disassembler tools together with IUC `mafft`, IUC `fasttree`, and (per section 4) `bcftools merge` where appropriate.
8. Submit tools and workflows to the IUC (`tools-iuc`) and IWC (`iwc`) repositories respectively, following their house styles.

The two open questions that do not appear as blocking steps are deferable: metadata enrichment (section 4) is a separate project connected into the workflow later, and Logan-vs-usegalaxy.org (section 4) is addressed as an addendum once a dedicated Logan instance exists.
