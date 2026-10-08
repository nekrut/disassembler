# Plan: port disassembler to Galaxy

Status: planning document, nothing in this file is built yet.

Goal: expose the disassembler suite as a set of Galaxy tools and workflows so a user can go from a kmindex planetary query (or a LexicMap tabular search) to a codon MSA and tree without leaving Galaxy. This document is the port's shopping list: the Bioconda packages that have to exist, the Galaxy tools that wrap them, and the knobs each tool exposes.

The port deliberately keeps the CLI unchanged. Every Galaxy tool is a thin wrapper around the same subcommand a terminal user runs today, so the two surfaces stay in sync as the CLI evolves. Flags that only make sense outside Galaxy (local caches, S3 fallbacks, API keys, binary path overrides) are not exposed by the wrappers; they stay on the CLI for workstation use.

---

## 1. Bioconda packages

The current tree ships three independently buildable artifacts: a Python package, a Rust binary, and a C binary. They sit in one repo for development reasons but have no shared build system, which makes the one-package-per-artifact split obvious. The master package is `disassembler`; it depends on the other two plus the existing Bioconda tools.

### 1.1 `logan-walker` (new)

The Rust path walker, built from `walker/`. One binary, `logan-walker`.

- Source: `walker/` in `as042/disassembler` at a tagged version.
- Build: `cargo build --release`. Pure Rust, two crate dependencies (`ureq`, `zstd`), no system libs beyond what Rust pulls in.
- Shipped artifact: the `logan-walker` executable on PATH.
- Runtime dependency: GraphAligner (shelled out via PATH). Expressed as a Bioconda `run` dependency so a user installing `logan-walker` picks up `graphaligner` automatically.
- Version pin: `graphaligner=1.0.20` for now, with a known regression noted in `docs/issue_graphaligner_end_clipping.md`. Switch to `1.0.17b` if the end-clipping fix matters to a given analysis.
- Why its own package: the walker is useful without the rest of disassembler (anyone with a Logan accession list and a reference CDS can use it), and keeping it independent lets Galaxy tool wrappers depend on only the walker when that is all they need.

### 1.2 `kmer-diversity-sieve` (new)

The MinHash/LSH diversity sieve from `extern/kmer-diversity-sieve`, currently a git submodule.

- Source: the upstream repo the submodule tracks.
- Build: `make serial` produces `l1_pipeline_serial`. Single C source file, no system deps beyond a C compiler.
- Shipped artifact: the `l1_pipeline_serial` executable on PATH.
- Why its own package: it was built for a general problem (archive-scale deduplication of near-identical consensus genomes) and has utility outside disassembler. If upstream already has a Bioconda recipe we use theirs; otherwise we contribute one.
- Open question: whether to ship `l1_pipeline_threads` and the MPI build too. The disassembler uses only `_serial`; a general-purpose Bioconda package should probably ship all three but default to serial.

### 1.3 `disassembler` (new, master)

The Python package providing the `disassembler` CLI and every subcommand.

- Source: `src/disassembler/` in `as042/disassembler`.
- Build: standard `pyproject.toml` wheel.
- Python runtime deps (already in `pyproject.toml`): `biopython`, `requests`, `zstandard`, `numpy`, `pandas`.
- External-tool `run` dependencies on Bioconda:
  - `logan-walker` (new, section 1.1)
  - `kmer-diversity-sieve` (new, section 1.2)
  - `graphaligner` (existing; transitive via logan-walker, listed here for clarity)
  - `mafft` (existing)
  - `fasttree` (existing)
- Not shipped: Rust and C toolchains. End users never build anything.
- Why "master": installing `disassembler` from Bioconda pulls in everything the CLI expects to find on PATH. A workstation user runs `conda install -c bioconda disassembler` and the full suite works. A Galaxy wrapper that needs only one piece (walker, sieve) depends on the narrower package directly.

### 1.4 Existing Bioconda packages we lean on

Not new, but worth listing so the dependency graph is explicit:

| package | used by |
|---|---|
| `graphaligner` | logan-walker (shells out) |
| `mafft` | disassembler msa |
| `fasttree` | disassembler msa |
| `zstd` | logan-walker (link-time via the `zstd` crate; and runtime for CLI inspection of cached files) |

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
- `disassembler` pins the minor versions of its dependencies, so a `disassembler 0.3.x` run reproduces regardless of whether logan-walker 0.0.2 or 0.0.3 is the latest.
- Golden-fixture concordance (`bench/check_golden.sh`) is the gate for a `disassembler` release.

---

## 2. Galaxy tools

Six tools cover both arms of the workflow. One existing subcommand (`fetch`) does not become a Galaxy tool; its purpose (local filesystem caching) is subsumed by Galaxy's data staging.

Convention used below: "exposed" means the flag appears in the tool XML as a user-visible param. "Not exposed" means the wrapper sets it from Galaxy context, hardcodes it, or omits it entirely.

### 2.1 `disassembler_harvest`

Merge kmindex planetary shard JSONs into per-target scored accession lists.

- Depends on: `disassembler`.
- Driven by: `disassembler harvest --shards-dir` (the local-files mode added for the port).

Inputs:
- Shard JSON collection (`data_collection_type="list"`, format `json`). This is the direct output of the IUC `kmindex_query` tool.

Exposed params: none beyond the input collection. Harvest's output shape is deterministic once the inputs are set.

Not exposed:
- `--history-id`, `--collection-id`, `--base-url`, `--api-key-file`: Galaxy API mode does not apply inside Galaxy.
- `--shards-dir`: set by the wrapper after staging the collection into a directory.
- `--download-workers`: ignored in local mode.

Outputs:
- A `list` collection of per-target `<target>.accessions.txt`.
- A `list` collection of per-target `<target>.scores.json`.
- `UNION.accessions.txt` as a single `txt` dataset.
- `summary.json` as a single `json` dataset.

### 2.2 `disassembler_stratify`

Cut a size-capped cohort from a harvest by score bin.

- Depends on: `disassembler`.
- Driven by: `disassembler stratify`.

Inputs:
- Per-target accession lists collection (from harvest).
- Per-target scores collection (from harvest).

Exposed params:
- `--score-bin LO:HI:N`: as a repeat element in the tool XML, so the workflow can request multiple bins.

Not exposed:
- `-o/--out-dir`: Galaxy-managed.

Outputs:
- Per-target subsampled accession list collection.
- Per-target bins-JSON collection.
- `bins_summary.json`.

### 2.3 `disassembler_walker`

Stream unitigs, build the local subgraph, walk haplotypes, emit variants.

- Depends on: `logan-walker` (which transitively depends on `graphaligner`).
- Driven by: `disassembler walker`, or `logan-walker` directly (the subcommand is a passthrough).

Design decision carried forward from `plan_rust_path_walker.md` thinking: expose the tool as **one-accession-per-job scattered over a list collection**, not as a batch tool reading an accession list file. Rationale: per-accession job isolation, retries, resource allocation, and provenance all become Galaxy-native. The internal-list mode of the CLI stays available for terminal users and for the scatter-per-shard variant discussed in section 3.3.

Inputs:
- One SRA accession (string or single-element collection, scatterable).
- Reference CDS FASTA (`-q`).

Exposed params:
- `--max-size-mb` (default 1024): per-accession size cap, important user policy knob.
- `--min-cov` (default 0.85): reconstruction quality gate.
- `--prune-frac` (default 0.01): reference-bias correction strength.
- `--prefer-majority` / `--no-prefer-majority` (default on): majority-allele preference at biallelic sites.
- Intra-host variant params: `--min-abund`, `--min-vaf`, `--min-abund-frac`, `--alt-max-detour`, `--no-indel-alts`, `--allow-frameshift-alts`, `--max-alt-indel`. Group in an "advanced" conditional block in the tool XML; most users leave them alone.
- Anchoring: `-k`, `--min-shared-kmers`, `--hops`. Advanced block, same reasoning.
- Output policy: `--emit-alts` / `--no-emit-alts`, `--emit-vcf`, `--emit-snv-table`. Expose as booleans; each controls whether the corresponding output dataset is emitted.
- `--allow-empty` (default off): exposed so a workflow can decide whether an all-fail shard is a soft or hard failure.

Not exposed:
- `--unitigs-cache`, `--offline`: Galaxy handles data staging. The walker streams S3 inside its job.
- `--walker-bin`: binary is on PATH via the `logan-walker` package.
- `--scratch-dir`: Galaxy job sandbox handles this.
- `-o`, `-f`: output paths Galaxy-managed.
- `-w, --workers`, `--fetch-workers`, `-t, --threads`: set from `\${GALAXY_SLOTS:-1}`. The walker's memory bound depends on `--max-size-mb * --fetch-workers * --workers`, so the wrapper either caps these at job-resource limits or does the arithmetic to pick safe values from GALAXY_SLOTS and GALAXY_MEMORY_MB.

Outputs:
- `reconstructed.fasta` (always).
- `alt_haplotypes.fasta` (filtered on `--emit-alts`).
- `<acc>.vcf` (filtered on `--emit-vcf`).
- `snvs.tsv` (filtered on `--emit-snv-table`).
- `bubbles.json`, `summary.json`.

### 2.4 `disassembler_lexicmap`

Stream a LexicMap tabular dataset, do multi-HSP coordinate tiling and codon QC, emit the cohort ledger and the collapsed clean cohort.

- Depends on: `disassembler`.
- Driven by: `disassembler lexicmap -f`.

Inputs:
- LexicMap tabular dataset (`tabular`).
- Reference CDS FASTA (`-q`).

Exposed params:
- `--min-cov` (default 0.80): clean full-length cohort gate.
- `--min-cov-partial` (default 0.50): clean partial cohort gate.
- `--min-pident` (default 60.0): HSP identity floor.
- `--max-internal-stops` (default 0): codon QC strictness.
- `--allow-frameshift-alts` (off): whether frameshift candidates reach the clean cohort.
- `--samples-per-haplotype` (default 10): representative accession cap in the haplotypes TSV.
- `--no-a2m` (off): disable the `.a2m` output.

Not exposed:
- `-d/--dataset`, `--base-url`, `--api-key`, `--api-key-file`: Galaxy API mode does not apply.
- `-o/--out-prefix`: Galaxy-managed.
- `--max-records`: a testing/preview limit, not a user policy.
- `-f/--file` with stdin: the wrapper always points at a staged tabular dataset.

Outputs:
- `clean.msa.fasta`, `clean.haplotypes.tsv`, `clean.accessions.fasta`, `clean_expanded.accessions.fasta`, `flagged.accessions.fasta`, `flagged.tsv`, `cohort_ledger.tsv`, `summary.json`, `clean.a2m`.

### 2.5 `disassembler_sieve`

k-mer diversity sieve: thin a cohort of near-duplicate reconstructions before the O(N^2) MSA/tree stage.

- Depends on: `kmer-diversity-sieve` (and `disassembler` for the wrapper's gene-scale defaults, or inline the defaults).
- Driven by: `disassembler sieve`.

Inputs:
- FASTA of candidate sequences to thin.

Exposed params:
- `--sieve-threshold` (tau, default 0.005): the gene-scale default, not the sieve binary's own default.
- `--sieve-sketch` (M, default 256): resolution floor.
- `--sieve-min-len` (default: 95% of query length, which the wrapper cannot compute without the query, so expose as an explicit integer).

Not exposed:
- `--sieve-bin`: binary is on PATH via the `kmer-diversity-sieve` package.

Outputs:
- Thinned FASTA.

Open question: whether to ship `disassembler_sieve` at all, or to encourage a separate IUC-grade `kmer_diversity_sieve_l1` tool with the sieve binary's own defaults and set the gene-scale defaults at the workflow level. The sieve itself is Stephen Weaver's project; the Galaxy-tool-shape question has not been raised with him yet. Section 4 carries the full note.

### 2.6 `disassembler_msa`

Reference-guided in-frame codon MSA (MAFFT `--keeplength`) plus FastTree GTR tree.

- Depends on: `disassembler`, `mafft`, `fasttree`.
- Driven by: `disassembler msa`.

Inputs:
- Reconstructed haplotypes FASTA.
- Reference CDS FASTA.

Exposed params:
- `--msa-min-cov` (default 0.70): MSA entry coverage gate.
- `--mafft-threads`: set from `\${GALAXY_SLOTS:-1}`.

Not exposed:
- Output paths: Galaxy-managed.
- `--skip-tree`: let the workflow decide whether to connect FastTree, instead of a tool-level flag.

Outputs:
- Codon MSA FASTA.
- FastTree Newick.

Open question: whether to ship `disassembler_msa` as one tool, or decompose into "reference-guided codon MSA prep" (small disassembler-specific step) + IUC `mafft` + IUC `fasttree`. Current recommendation is to decompose, so the Galaxy port does not re-wrap community tools that already exist in the IUC. Section 4 carries the full note.

### 2.7 Not ported: `disassembler_fetch`

`fetch` populates a local filesystem cache for `--unitigs-cache` / `--offline`. Galaxy's data model subsumes this: tool outputs are managed datasets in the history and are already cached across downstream jobs. Shipping fetch as a tool would require either (a) producing a collection of 51.6 MB median unitig files that charges against user storage quota, (b) a data manager writing into a shared admin-managed cache, which only makes sense on a dedicated Logan Galaxy instance.

Decision: drop `fetch` from the Galaxy port. The CLI keeps it for workstation use.

### 2.8 Not ported: the composites

`reconstruct`, `from-kmindex`, and `from-lexicmap` are the three CLI composites. On Galaxy they become `.ga` workflow files (section 3), not tools.

---

## 3. Workflows

Three `.ga` files cover the current CLI composites, with one extra that falls out for free.

### 3.1 `from-kmindex.ga`

Mirror of `disassembler from-kmindex`. Full planetary-scale triage.

```
kmindex_query (IUC)
  -> disassembler_harvest
  -> disassembler_stratify  (optional, controlled by a workflow conditional)
  -> disassembler_walker    (scattered over the stratified accession list)
  -> disassembler_sieve     (optional)
  -> disassembler_msa
```

Workflow-level choices that were CLI flags:
- `--sieve` / `--no-sieve` on the CLI becomes "connect the sieve box or wire walker to msa directly" in the workflow editor.
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

The walker tool is defined as one-accession-per-job (section 2.3). Three scatter strategies become selectable at the workflow level, matching the discussion in `plan_galaxy_port.md`'s predecessor conversation:

1. **Scatter-per-accession** (default): natural Galaxy pattern, per-accession retry and provenance, breaks at >~10k accessions due to scheduler overhead.
2. **Scatter-per-shard**: pre-split the accession list into shards of N, scatter N jobs, each shard runs the walker with its own internal list loop. Clean up to >100k accessions. Requires a tiny `disassembler_shard_list` splitter tool, or using Galaxy's built-in collection ops.
3. **Single-job batch**: one walker call processes the whole list internally. Loses per-accession isolation, scales to any size memory and time allow.

For SARS-CoV-2-scale cohorts (10^5-10^6 accessions) only option 2 is realistic. Workflow-level routing selects the strategy.

---

## 4. Open questions and deferred decisions

Resolving these is part of the plan, not a prerequisite to it. Each question lists its status, any outside owner the answer waits on, and the current recommendation when there is one.

- **Sieve tool boundary.** Section 2.5. Ship `disassembler_sieve` as a thin wrapper with gene-scale defaults baked in, or contribute a general-purpose `kmer_diversity_sieve_l1` IUC tool with upstream defaults and set gene-scale values at the workflow level. The second is better Galaxy citizenship; the first is less for a workflow user to understand. *Status: open.* The sieve project is Stephen Weaver's; we have been working closely with him. He has confirmed that `kmer-diversity-sieve` should have its own Bioconda package (section 1.2), but the equivalent question about a dedicated IUC Galaxy tool has not been raised with him yet. **Blocked on asking him.** No recommendation until that conversation happens.
- **MSA tool boundary.** Section 2.6. Ship `disassembler_msa` as one tool, or decompose into a small reference-guided codon prep tool + IUC MAFFT + IUC FastTree. One tool is easier to use; decomposed is more inspectable. *Status: open; current recommendation is to decompose.* Rationale: avoid re-wrapping MAFFT and FastTree, which are already well-maintained in the IUC. The disassembler-specific value is the reference-guided in-frame codon preservation layer, and that is what gets its own tool. The decomposition also makes `--skip-tree` disappear naturally (just do not connect the FastTree box in the workflow).
- **Metadata enrichment.** Not currently in the CLI. For SARS-CoV-2 at scale the harvest's 5M candidate list needs an SRA metadata join (date, platform, country) before stratify, so that a stratified cohort is informative, not random. This likely lives in the sibling `sra-metadata-enrichment` project and connects into the workflow, rather than in disassembler. *Status: open, out of scope for this port.*
- **VCF merging downstream of walker.** The CLI's `snvs.tsv` long table is the scalable output; a merged multi-sample VCF is a convenience. On Galaxy the natural shape is `bcftools merge` (IUC) wired to the walker's per-accession VCF outputs. Not a disassembler tool; worth calling out in the workflow docs. *Status: open, workflow-doc concern not a tool concern.*
- **GraphAligner version.** `environment.yml` pins 1.0.20 but 1.0.19/1.0.20 regressed query end-clipping (`docs/issue_graphaligner_end_clipping.md`). The `logan-walker` Bioconda package should default to a known-good version; 1.0.17b for correctness-sensitive work, 1.0.20 for parity with the committed golden fixtures (which were captured under the affected build). The right default depends on whether a Galaxy run is expected to byte-reproduce the committed panels, which it probably is not. *Status: open.* Affects correctness of every Galaxy output, so the resolution has to land before the first published Galaxy tool version; probably forces a goldens rebuild.
- **Logan Galaxy instance vs usegalaxy.org.** The port is drafted assuming usegalaxy.org, which charges user storage and has rate limits. A dedicated Logan Galaxy instance changes two things: the `fetch` tool becomes attractive again as a data manager, and the walker's `--max-size-mb` policy could default to something larger. *Status: open; deferable.* The port works on usegalaxy.org as drafted; a Logan-instance variant can land later as an addendum.
- **Who maintains the Bioconda recipes.** The disassembler side is natural for the disassembler maintainers. `kmer-diversity-sieve` is upstream of us (Stephen Weaver's project); preferred path is to contribute the recipe upstream and track their releases, not vendor it. *Status: open, pending the same conversation as the sieve tool boundary.*

---

## 5. Order of work

A reasonable sequence; none of the steps block the next past where noted. The questions in section 4 appear here as explicit resolve-this steps so they do not get skipped.

1. **Ask Stephen about sieve Galaxy-tool shape and recipe ownership.** Resolves two open questions at once (section 4, sieve tool boundary; section 4, who maintains the `kmer-diversity-sieve` recipe). Determines whether a `disassembler_sieve` wrapper is written at step 4 or dropped in favour of a general IUC tool, and who authors the Bioconda recipe. **Blocks steps 2 and the sieve part of step 4.**
2. Contribute (or receive) a Bioconda recipe for `kmer-diversity-sieve`. Direction set by step 1. **Blocks step 4.**
3. Decide the GraphAligner version pin for `logan-walker` (section 4). If 1.0.17b wins, rebuild `bench/golden/` before the port ships so the concordance gate stays live. **Blocks step 5.**
4. Tag `logan-walker 0.0.2` from `walker/` with the current tree. Contribute Bioconda recipe depending on step 3's version pin. **Blocks step 6.**
5. Tag `disassembler 0.3.0`. Contribute Bioconda recipe depending on steps 2 and 4.
6. Author the Galaxy tool XMLs in a new `galaxy/` directory in this repo. Current count is six, pending step 1's answer on sieve. Order of least-to-most-tricky: `sieve` (if we ship it), `msa` prep (decomposed per section 4's recommendation), `lexicmap`, `harvest`, `stratify`, `walker`.
7. Author the three `.ga` workflow files, wiring the disassembler tools together with IUC `mafft`, IUC `fasttree`, and (per section 4) `bcftools merge` where appropriate.
8. Submit tools and workflows to the IUC (`tools-iuc`) and IWC (`iwc`) repositories respectively, following their house styles.

The two open questions that do not appear as blocking steps are deferable: metadata enrichment (section 4) is a separate project connected into the workflow later, and Logan-vs-usegalaxy.org (section 4) is addressed as an addendum once a dedicated Logan instance exists.
