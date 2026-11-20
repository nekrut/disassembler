# logan-walker

The Rust walker of [disassembler](../README.md). For each SRA accession it streams the Logan
unitig graph from S3, anchors a reference coding sequence by shared k-mers, walks the local
subgraph with GraphAligner, and stitches the walked path into a coding haplotype: one
majority-consensus primary per accession, plus the intra-host alternatives it had to choose
between. What it does and how it behaves is in [`../docs/walker.md`](../docs/walker.md); every
flag, output and exit code is in
[`../docs/walker_reference.md`](../docs/walker_reference.md).

## Build

```bash
cargo build --release            # target/release/logan-walker
cargo test --release             # offline tests; network and data tests are #[ignore]
cargo clippy --all-targets --release -- -D warnings
```

GraphAligner 1.0.17b must be on `PATH` at run time (see [`../environment.yml`](../environment.yml)
and [`../docs/notes/graphaligner_end_clipping.md`](../docs/notes/graphaligner_end_clipping.md)).
Rust installs via [rustup](https://rustup.rs). Six runtime crates and no dev-dependencies:
`ureq` (HTTP), `zstd` and `flate2` (decompression), `sha2` (provenance digests), `libc` and
`signal-hook` (the allocator pin, process groups and signal handling). The packaged build is on
Bioconda as `logan-walker`.

## Run

```bash
logan-walker -a SRR8209082 -q ../queries/query_hiv1_pol.fasta -o out -f out/reconstructed.fasta
```

`disassembler walk` passes its arguments straight to this binary, and `disassembler
reconstruct` runs it as its first stage.

## Source map

```text
src/
  main.rs        CLI, config, IN -> pipeline -> OUT
  pipeline.rs    two-pool threading and stage wiring   <- read this for the flow
  types.rs       what moves between stages             <- read this for the shapes
  config.rs      tunables
  error.rs       error type
  scratch.rs     the run's private scratch directory, its lock, and signal cleanup
  term.rs        terminal detection and progress output
  io/
    input.rs     accession list and query in
    output.rs    FASTA, VCF, SNV table, bubbles, summary and crash journal out
  stages/
    stream.rs    S3 fetch and streaming zstd decode
    parse.rs     BCALM2 records -> CSR graph and 2-bit sequences, in two passes
    anchor.rs    canonical k-mer anchoring
    extract.rs   hop-limited BFS, sequence collection, GFA emit, prune
    align.rs     GraphAligner wrapper
    stitch.rs    path -> haplotype, primary selection in multi-anchor mode
    variants.rs  prefer-majority and intra-host alternative recovery (alts, VCF, SNV table)
    phase.rs     bubble detection (bubbles.json)
tests/
  cli.rs         end-to-end runs of the built binary, with a stand-in aligner
```

`parse.rs` reads each accession twice. The first pass builds the adjacency and finds the
anchors as the bases go by, spilling the adjacency to disk once it outgrows 4 MiB; the second
keeps sequence and abundance for the extracted neighbourhood alone. That is what keeps memory a
function of the gene's neighbourhood rather than of the accession.

The golden concordance gate for this crate is
[`../bench/check_golden.sh`](../bench/check_golden.sh).
