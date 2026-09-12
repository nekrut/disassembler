# logan-walker

The native Rust engine of [disassembler](../README.md). Streams Logan cDBG unitigs from S3,
anchors a reference coding sequence against them, walks the local subgraph, and reconstructs the
full-length coding haplotypes present in the sample.

## Build

```bash
cargo build --release            # target/release/logan-walker
cargo test --release             # offline tests (network/data tests are #[ignore])
```

`GraphAligner` must be on `PATH` at run time (bioconda; see [`../environment.yml`](../environment.yml)).
Rust itself installs via [rustup](https://rustup.rs). Only three runtime crates: `ureq`, `zstd`,
and `sha2` (tests only).

## Run

```bash
logan-walker -a SRR8209082 -q ../queries/query_hiv1_pol.fasta -o out -f out/reconstructed.fasta
```

Key flags (`logan-walker --help` for the full list):

| flag | meaning |
|---|---|
| `-a` / `-l` | one accession, or a file listing accessions |
| `-q` | reference CDS FASTA |
| `-o` / `-f` | output directory / reconstructed-haplotypes FASTA |
| `-k` | anchor k-mer size (default 21) |
| `--prune-frac` | drop unitigs below this fraction of gene depth, then re-align (default 0.01) |
| `--prefer-majority` | put the majority allele on the primary at biallelic sites (default on) |
| `--emit-alts` | write intra-host alternative haplotypes to `alt_haplotypes.fasta` (default on) |
| `--emit-vcf` | write per-accession `<acc>.vcf` of intra-host SNPs (default off) |
| `--emit-snv-table` | write one long `snvs.tsv` across the run — the cohort output that scales (default off) |
| `-w` / `--fetch-workers` | compute workers / S3 fetch workers |

Reproduce the old one-haplotype-per-accession behaviour with
`--prune-frac 0.0 --no-prefer-majority --no-emit-alts`.

## Outputs

`reconstructed.fasta` (primaries), `alt_haplotypes.fasta`, `bubbles.json`, a summary JSON, and
optionally `<acc>.vcf` / `snvs.tsv`.

## Source map

```
src/
  main.rs        CLI, config, IN -> pipeline -> OUT
  pipeline.rs    two-pool threading and stage wiring   <- read this for the flow
  types.rs       what moves between stages             <- read this for the shapes
  config.rs      tunables
  io/            accession list + query in; fasta, vcf, snv table, bubbles, summary out
  stages/
    stream.rs    S3 fetch + streaming zstd decode
    parse.rs     BCALM2 records -> CSR graph + 2-bit sequences
    anchor.rs    canonical k-mer anchoring (replaces BLAST)
    extract.rs   5-hop BFS, sequence collection, GFA emit, prune
    align.rs     GraphAligner wrapper
    stitch.rs    path -> haplotype
    variants.rs  prefer-majority + intra-host alt recovery (SNPs -> alts, vcf, snv table)
    phase.rs     bubble detection (bubbles.json)
```

The walker reconstructs haplotypes; to codon-align them and build a tree, feed
`reconstructed.fasta` to `../python/build_walker_codon_msa.py`, or run the whole thing in one
command with `../python/disassemble.py` (see [`../docs/codon_msa.md`](../docs/codon_msa.md)).

See [`../docs/`](../docs/) for the design plan, as-built notes, the variance-loss
investigation, the codon-MSA step, and the VCF / contig experiments.
