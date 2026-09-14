<p align="center">
  <img src="assets/logo.jpg" alt="logan-disassembler" width="640">
</p>

# disassembler

Planetary-scale sequence extraction, de novo graph path-walking, and streaming multi-HSP coordinate assembly from the [Logan](https://github.com/IndexThePlanet/Logan) SRA petabase.

`disassembler` provides two primary, complementary engines:

1. **`logan-walker`** (Native Rust graph path-walker) — Streams unitig cDBGs from S3, anchors a reference coding sequence by canonical $k$-mer containment, extracts the local subgraph, walks paths with GraphAligner, and reconstructs full-length coding haplotypes, intra-host alternative alleles, and biallelic SNVs (`--emit-vcf`, `--emit-snv-table`) without collapsing intra-host variation or suffering reference-bias loss.
2. **`LexicMapStreamer`** (`python/lexicmap_streamer.py`) — Streams planetary tabular search outputs from [LexicMap](https://github.com/IndexThePlanet/LexicMap) over the Galaxy HTTP API, local TSV, or STDIN. Resolves Logan unitig fragmentation by collinear multi-HSP coordinate tiling on the reference CDS, enforces reading-frame codon QC (screening premature stop codons like Sanger $am3$ amber mutants), partitions accessions into a 4-tier Clean and Flagged cohort ledger, and collapses full-length clean sequences into unique haplotypes online with mutation calling and Shannon entropy.

This repository is a decoupled extraction of the path-walking and sequence assembly components from
[`nekrut/logan-hypyaeon`](https://github.com/nekrut/logan-hypyaeon). It carries only the
pieces needed to go from an SRA accession list or LexicMap search plus a query gene to reconstructed haplotypes and phylogenies.
The planetary-scan orchestration, HyphAeon neural selection suite, and Roland PLM branch are intentionally left behind.

## Layout

```
walker/     logan-walker: native Rust local graph path-walker (the core walker engine)
python/     lexicmap_streamer.py (streaming multi-HSP coordinate tiler, codon QC & haplotype collapser)
            disassemble.py (end-to-end driver: accessions -> walker -> codon MSA -> FastTree)
            build_walker_codon_msa.py (reference-guided in-frame codon MSA + FastTree ML tree)
            harvest_kmindex.py (parallel Galaxy kmindex planetary shard harvester & cohort triage)
queries/    reference coding sequences (pfcrt, hiv1_pol, 25 planetary markers, master panel)
bench/      golden fixtures, check_golden.sh concordance gate, k-mer sensitivity simulations
docs/       codon MSA/phylogeny step, VCF example, contigs-vs-unitigs, variance-loss investigation
```

## What logan-walker does

For each SRA accession, given a reference CDS:

1. **Stream** the compressed unitig file from `s3://logan-pub/u/<ACC>/<ACC>.unitigs.fa.zst`
   directly into memory, rejecting oversized accessions on `Content-Length` first.
2. **Parse** the BCALM2 FASTA in a single streaming pass into a compact CSR graph
   (adjacency + `ka:f` abundance + 2-bit sequences, ~43 MB for a 1.28M-unitig accession).
3. **Anchor** the query by canonical k-mer containment (k=21), replacing BLAST.
4. **Extract** the 5-hop neighbourhood around the anchors as a local GFA.
5. **Align** the reference through the local graph with `GraphAligner -x dbg`.
6. **Stitch** the walked path into a primary haplotype, optionally **prune** reference-bias
   error branches and re-align, **prefer** the majority allele at biallelic sites, and
   **emit** intra-host alternative haplotypes to a separate file.
7. Write `reconstructed.fasta`, `alt_haplotypes.fasta`, `bubbles.json`, and a summary JSON;
   optionally the intra-host SNP calls (`--emit-vcf`, `--emit-snv-table`).

Then, optionally, an eighth downstream step turns the reconstructed haplotypes into a
phylogeny:

8. **Codon-align and build a tree.** Reference-guided in-frame codon MSA (MAFFT `--keeplength`,
   so reference coordinates and the reading frame are preserved) and a FastTree GTR
   maximum-likelihood tree. Run it as one command with the walker via `python/disassemble.py`
   (see below), or on its own with `python/build_walker_codon_msa.py`. See
   [`docs/codon_msa.md`](docs/codon_msa.md).

### Why unitigs, not Logan contigs

Logan also publishes assembled contigs, which are far smaller. On internally diverse samples,
though, the assembler shatters the gene: measured on nine deep Zika NS5 accessions, no single
contig held more than 8% of the gene and it took ~50 fragments to tile it, while the unitig walk
reconstructed the full 2,708 nt from every one. A clonal control assembled into one clean
contig, so the difference is intra-host diversity, not a Logan quirk. Full numbers in
[`docs/contigs_vs_unitigs.md`](docs/contigs_vs_unitigs.md).

## Build and run

```bash
# GraphAligner must be on PATH (bioconda). k-mer anchoring replaces BLAST; MAFFT is only for
# the downstream codon MSA in python/.
cd walker
cargo build --release

./target/release/logan-walker \
  -a SRR8209082 \
  -q ../queries/query_hiv1_pol.fasta \
  -o out -f out/reconstructed.fasta
```

Reproduce the old one-haplotype-per-accession behaviour (matches the golden fixtures) with
`--prune-frac 0.0 --no-prefer-majority --no-emit-alts`.

### End to end: accessions to a phylogeny

`python/disassemble.py` chains the walker and the codon-MSA/FastTree step into one command, so
you go straight from accessions to a tree. It needs `GraphAligner`, `MAFFT`, `FastTree`, and
biopython — all in the `disassembler` conda env (`environment.yml`).

```bash
conda env create -f environment.yml && conda activate disassembler
python/disassemble.py -a SRR6251967 -q queries/markers/zika_ns5.fasta -o out/zika --emit-snv-table
#   out/zika/reconstructed.fasta       reconstructed coding haplotypes
#   out/zika/zika_ns5_codon_msa.fasta  in-frame codon alignment (MAFFT --keeplength)
#   out/zika/zika_ns5_fasttree.nwk     GTR maximum-likelihood tree
```

Use `-l accessions.txt` for a cohort, `--skip-tree` to stop after reconstruction, and pass any
extra walker flags after a literal `--`. How the codon alignment stays in frame, and what the
outputs contain, is in [`docs/codon_msa.md`](docs/codon_msa.md).

### Planetary triage: Galaxy kmindex harvest to candidate cohorts

When screening planetary datasets indexed in Logan (e.g. via Galaxy `kmindex query` across thousands of shard databases), `python/harvest_kmindex.py` streams and merges the shard collections in parallel, extracts candidate SRA accessions by containment score, and creates stratified diversity cohorts:

```bash
# Stream and merge ~2,869 Galaxy shard JSONs across 32 threads
python/harvest_kmindex.py \
    --history-id <GALAXY_HISTORY_ID> \
    --collection-id <GALAXY_COLLECTION_ID> \
    -o candidates/ --diversity 2500

# Reconstruct a target gene cohort directly with disassemble.py
python/disassemble.py \
    -l candidates/query_gene.diversity_2500.txt \
    -q queries/query_gene.fasta \
    -o out/query_gene \
    --emit-vcf --emit-snv-table --no-emit-alts -w 16 --fetch-workers 16
```

## LexicMapStreamer: Streaming Multi-HSP Assembly & Cohort Partitioning

`python/lexicmap_streamer.py` provides streaming multi-HSP coordinate tiling, reading-frame codon QC, 4-tier cohort stratification, and online haplotype collapsing for planetary search results from [LexicMap](https://github.com/IndexThePlanet/LexicMap).

### The Multi-HSP Unitig-Split Challenge

In Logan unitig collections, sequence diversity, coverage fluctuations, and de Bruijn graph branchings frequently break a coding gene across multiple collinear unitigs. In standard LexicMap tabular outputs, each unitig is reported as an individual high-scoring pair (HSP) with its own `hsp`, `qcovHSP`, `alenHSP`, `qstart`, `qend`, and `pident`:

- A gene split across two unitigs (e.g. 55% coverage on unitig 1 and 45% on unitig 2) fails naive single-row filters (like `qcovHSP >= 80%`).
- Filtering by row coverage discards complete biological sequences or produces fragmented assemblies.

### How LexicMapStreamer Works

1. **Collinear Multi-HSP Coordinate Tiling:**
   Groups search rows on the fly by SRA accession (`sgenome`). For each accession, projects all qualifying HSPs onto the reference coordinate system (`qstart` to `qend`). Overlapping positions between HSPs are resolved by selecting the base from the HSP with the higher percent identity (`pident`), seamlessly stitching split unitigs into full-length sequences without introducing spurious indels.
2. **Reading-Frame Codon QC & Stop-Codon Screening:**
   Translates CDS in-frame using Genetic Code Table 11 (bacterial/standard). Scans the coding body (excluding the natural terminal stop codon) for internal premature stop codons (e.g. amber mutants like Sanger $\Phi$X174 $am3$ dual-coding $gpE\text{ }W7^*$, sequencing errors, or out-of-frame indels).
3. **Four-Tier Cohort Partitioning:**
   - **Tier 1: Clean Full-Length** (`status=CLEAN_FULL_LENGTH`): Coverage $\ge 80\%$ (configurable via `--min-coverage`) and 0 premature stop codons.
   - **Tier 2: Clean Partial** (`status=CLEAN_PARTIAL`): Coverage between $50\%\text{--}80\%$ and 0 premature stop codons. Rescues partial hits without contaminating full-length alignments.
   - **Tier 3: Flagged Premature Stops** (`flag=PREMATURE_STOPS`): Premature stop codons detected; outputs exact 1-based codon coordinates in an audit ledger.
   - **Tier 4: Low-Coverage Fragments** (`flag=LOW_COVERAGE`): Coverage $<50\%$.
4. **Online Haplotype Collapser:**
   Collapses clean full-length sequences in $O(1)$ streaming memory, calculates exact frequencies, calls nucleotide (`198G>A`) and amino acid (`V66V`, `I61V`) mutations against the reference, classifies mutation types (`WT`, `Synonymous`, `Missense`, `Nonsense`, `Complex`), computes Shannon entropy, and records representative sample accessions.
5. **Regime Diagnostics & Cross-Tool Triage:**
   Computes coverage and identity quantiles ($p_{10}$, median, $p_{90}$) and stop-codon rates. If a target is in a **divergent regime** (median identity $<85\%$ or coverage $<0.60$), reference-guided projection fragments; `LexicMapStreamer` issues an automated alert recommending de novo graph path-walking via `logan-walker`.

### Usage Examples

```bash
# 1. Stream directly from Galaxy over HTTP (using Galaxy dataset ID):
python/lexicmap_streamer.py \
    -q queries/markers/zika_ns5.fasta \
    -d f9cad7b01a47213500416a19986dca5c \
    -o out/zika_lexicmap

# 2. Process a local LexicMap tabular file (or gzipped TSV):
python/lexicmap_streamer.py \
    -q queries/query_hiv1_pol.fasta \
    -f /data/hiv1_pol_lexicmap.tsv.gz \
    -o out/hiv1_pol

# 3. Stream from standard input in a shell pipeline:
curl -s "https://usegalaxy.org/api/datasets/<ID>/display?key=<API_KEY>" | \
    python/lexicmap_streamer.py -q queries/master/pfcrt.fasta -f - -o out/pfcrt
```

### Outputs Generated

| Output File | Description |
|---|---|
| `<prefix>.clean.msa.fasta` | Collapsed clean haplotypes with rank, count, and frequency |
| `<prefix>.clean.haplotypes.tsv` | Tabular haplotype catalog with amino acid/nucleotide mutations, mutation class, and representative accessions |
| `<prefix>.clean.accessions.fasta` | Uncollapsed full-length clean sequences for every qualifying accession |
| `<prefix>.clean_expanded.accessions.fasta` | Expanded clean cohort combining full-length and clean partials ($\ge 50\%$) |
| `<prefix>.flagged.accessions.fasta` | All flagged sequences (premature stops or low coverage) |
| `<prefix>.flagged.tsv` | Audit ledger of flagged accessions with stop codon positions |
| `<prefix>.cohort_ledger.tsv` | Complete audit ledger recording coverage, identity, status, tier, and stops for every accession |
| `<prefix>.summary.json` | Machine-readable summary: quantiles, cohort counts, multi-HSP rescue statistics, and Shannon entropy |
| `<prefix>.clean.a2m` | A2M formatted alignment for protein language models and evolutionary couplings |

### LexicMapStreamer CLI Reference

| Flag | Default | Meaning |
|---|---|---|
| `-q, -r, --query, --ref <FASTA>` | — | Reference in-frame CDS FASTA (required) |
| `-d, --dataset <ID>` | — | Galaxy dataset ID to stream over HTTP |
| `-f, -i, --file, --input <FILE>` | — | Path to local tabular file, `.gz` file, or `-` for stdin |
| `-o, --out-prefix <PATH>` | `msa/<ref>` | Output directory or file prefix |
| `--base-url <URL>` | `https://usegalaxy.org` | Galaxy server base URL |
| `--api-key <KEY>` | — | Galaxy API key string |
| `--api-key-file <PATH>` | `/tmp/gxy.txt` | Path to Galaxy API key file |
| `--min-coverage <F>` | `0.80` | Minimum merged query coverage for clean full-length cohort |
| `--min-coverage-partial <F>` | `0.50` | Minimum coverage for clean partial cohort |
| `--min-pident <F>` | `60.0` | Minimum percent identity for candidate HSPs |
| `--max-internal-stops <N>` | `0` | Maximum premature internal stops allowed in CDS |
| `--allow-frameshifts` | off | Allow indels/premature stops into clean cohort |
| `--sample-cap <N>` | `10` | Maximum sample accessions recorded per haplotype in TSV |
| `--max-records <N>` | — | Stop after processing $N$ lines (for testing/previews) |
| `--no-a2m` | off | Disable `.a2m` output |

## logan-walker Command-line reference

`logan-walker [OPTIONS]`, one of `-a`/`-l` and `-q` required. Boolean options that are **on by
default** are turned off with their `--no-…` form; boolean options that are **off by default**
are turned on with their bare form (and back off with `--no-…`, e.g. inside a wrapper). Run
`logan-walker --help` for the same list.

### Input and output

| flag | default | meaning |
|---|---|---|
| `-a, --accession <ACC>` | — | a single SRA accession |
| `-l, --accession-list <FILE>` | — | file with one accession per line (use instead of `-a`) |
| `-q, --query <FASTA>` | — | reference in-frame CDS FASTA (required) |
| `-o, --out-dir <DIR>` | `walker_output` | output directory (also holds `bubbles.json`, summary, VCFs, `snvs.tsv`, GFA scratch) |
| `-f, --out-fasta <FILE>` | `reconstructed_haplotypes.fasta` | reconstructed primary haplotypes |

### Anchoring and fetch

| flag | default | meaning |
|---|---|---|
| `-k, --k <N>` | `21` | anchor k-mer size, 1–32 (the Logan graph k stays 31; this is the query-anchoring k) |
| `--min-shared-kmers <N>` | `3` | shared canonical k-mers a unitig needs to be an anchor (validated 1–5) |
| `--hops <N>` | `5` | hop radius of the local-subgraph BFS around the anchors |
| `--max-size-mb <N>` | `1024` | reject an accession whose compressed size exceeds this, on `Content-Length`, before download |
| `--fetch-workers <N>` | `8` | parallel S3 fetch streams (bandwidth saturates ~8) |
| `--scratch-dir <DIR>` | `/dev/shm/logan-walker` | where per-accession GFA/GAF go; should be tmpfs |

### Walk and consensus

| flag | default | meaning |
|---|---|---|
| `--min-cov <F>` | `0.85` | minimum fraction of the reference a haplotype must cover to be kept |
| `--prune-frac <F>` | `0.01` | drop unitigs below this fraction of the gene's depth and re-align, removing the aligner's reference bias; `0.0` disables |
| `--prefer-majority` / `--no-prefer-majority` | on | put the majority allele on the primary at two-allele sites (off = leave the reference-matching allele) |
| `-t, --threads <N>` | `1` | threads per GraphAligner call; above 1 is not reproducible |
| `-w, --workers <N>` | `8` | compute workers (size to cores) |

### Intra-host variants and their outputs

| flag | default | meaning |
|---|---|---|
| `--min-abund <F>` | `5.0` | minimum `ka:f` for a branch to count as a real allele |
| `--min-vaf <F>` | `0.10` | minimum variant allele frequency for an alternative branch |
| `--min-abund-frac <F>` | `0.0` (off) | raise the variant floor to this fraction of gene depth (median `ka:f` on the walked path) |
| `--alt-max-detour <N>` | `100` | unitigs an alternative allele may traverse before it must rejoin the route |
| `--no-indel-alts` | (indels on) | substitutions only; refuse alternatives that change length |
| `--allow-frameshift-alts` | off | also accept indel alleles whose length change is not a multiple of 3 (these reframe the CDS) |
| `--max-alt-indel <N>` | `30` | largest length change still treated as an indel allele rather than a rerouting |
| `--emit-alts` / `--no-emit-alts` | on | write intra-host alternative haplotypes to `alt_haplotypes.fasta` |
| `--emit-vcf` | off | also write `<out-dir>/<acc>.vcf`: intra-host SNPs with `ka:f` frequency, in query coordinates |
| `--emit-snv-table` | off | also write `<out-dir>/snvs.tsv`: the long cohort table, one row per (accession, site, allele) |

### Misc

| flag | meaning |
|---|---|
| `-h, --help` | print help |
| `-V, --version` | print version |

## Provenance and status

The walker began as `logan-hypyaeon` PR #4 (`Haplotype recovery`), which layered intra-accession
variant recovery on top of the PR #3 Rust port. Intra-host alleles are recovered by
`stages::variants`; the parallel upstream `stages::detour` engine has been retired in favour of
it (`docs/issue_intra_accession_variance_lost.md` section 17 records why). Concordance: the Rust
primary is byte-identical to `bench/golden/ERR2028820.pfcrt.primary.fasta` under old-behaviour
flags, checked by `bench/check_golden.sh`. Offline unit tests pass; network- and
data-dependent tests are `#[ignore]` by default.

See `docs/issue_intra_accession_variance_lost.md` for the central methodological finding: the
reference-matching-minority-allele bias, its measurement, and the fix.

## Credits

The Rust walker and intra-host haplotype-recovery method originate from `logan-hypyaeon` PR #4, authored by Andrew ([@as042](https://github.com/as042)), building on the PR #3 Rust port. See `NOTICE` and `docs/` for full provenance.

## Intra-host SNP VCF

`--emit-vcf` writes `<out-dir>/<accession>.vcf` alongside the primary: the on-path
biallelic sites (intra-host SNPs) with their k-mer-coverage frequency, in query
coordinates (`CHROM` = gene, `POS` = 1-based position in the CDS). The primary FASTA is
unchanged. Because every accession is called against the same gene, the VCFs share a
coordinate frame. Frequencies are `ka:f` ratios, a proxy for allele frequency, and are
per-site, not phased.

For a cohort, `--emit-snv-table` writes `<out-dir>/snvs.tsv`, one long ("tidy") row per
`(accession, site, allele)`. This is the output that scales: a gene can be matched by millions
of runs, and a merged multi-sample VCF would be a million-column near-empty matrix, whereas the
long table is sparse-native, concatenates across cluster shards, and loads straight into DuckDB
for per-site summaries and cross-sample queries. Merging the per-accession VCFs into a dense
matrix is a small-cohort convenience only.

A worked example on the Zika NS5 cohort — an example VCF, the long table with DuckDB queries,
and (for the small cohort) the merged matrix and a heatmap — is in
[`docs/vcf_example.md`](docs/vcf_example.md).
