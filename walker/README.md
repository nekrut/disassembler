# logan-walker

The Rust walker of [disassembler](../README.md). For each SRA accession it streams the Logan
unitig graph from S3, anchors a reference coding sequence by shared k-mers, walks the local
subgraph with GraphAligner, and reconstructs the full-length coding haplotypes present in the
sample. The user documentation, including the full option reference, is
[`../docs/walker.md`](../docs/walker.md).

## Build

```bash
cargo build --release            # target/release/logan-walker
cargo test --release             # offline tests; network and data tests are #[ignore]
cargo clippy --all-targets --release -- -D warnings
```

GraphAligner 1.0.17b must be on `PATH` at run time (see [`../environment.yml`](../environment.yml)
and [`../docs/issue_graphaligner_end_clipping.md`](../docs/issue_graphaligner_end_clipping.md)).
Rust installs via [rustup](https://rustup.rs). Two runtime crates, `ureq` and `zstd`, plus
`sha2` for tests. The packaged build is on Bioconda as `logan-walker`.

## Run

```bash
logan-walker -a SRR8209082 -q ../queries/query_hiv1_pol.fasta -o out -f out/reconstructed.fasta
```

`disassembler walker` passes its arguments straight to this binary, and `disassembler
reconstruct` runs it as its first stage.

## Source map

```
src/
  main.rs        CLI, config, IN -> pipeline -> OUT
  pipeline.rs    two-pool threading and stage wiring   <- read this for the flow
  types.rs       what moves between stages             <- read this for the shapes
  config.rs      tunables
  error.rs       error type
  io/            accession list and query in; FASTA, VCF, SNV table, bubbles, summary out
  stages/
    stream.rs    S3 fetch and streaming zstd decode
    parse.rs     BCALM2 records -> CSR graph and 2-bit sequences
    anchor.rs    canonical k-mer anchoring
    extract.rs   hop-limited BFS, sequence collection, GFA emit, prune
    align.rs     GraphAligner wrapper
    stitch.rs    path -> haplotype, primary selection in multi-anchor mode
    variants.rs  prefer-majority and intra-host alternative recovery (alts, VCF, SNV table)
    phase.rs     bubble detection (bubbles.json)
```

The golden concordance gate for this crate is [`../bench/check_golden.sh`](../bench/check_golden.sh).
