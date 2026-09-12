<p align="center">
  <img src="assets/logo.jpg" alt="logan-disassembler" width="640">
</p>

# disassembler

Stream [Logan](https://github.com/IndexThePlanet/Logan) cDBG unitigs from S3, anchor a
reference coding sequence against them on the fly, walk the local subgraph, and reconstruct
the full-length coding haplotypes present in the sample without collapsing intra-host
variation.

This repository is a decoupled extraction of the path-walking component from
[`nekrut/logan-hypyaeon`](https://github.com/nekrut/logan-hypyaeon). It carries only the
pieces needed to go from an SRA accession list plus a query gene to reconstructed haplotypes.
The planetary-scan orchestration, HyphAeon neural selection suite, Galaxy/kmindex candidate
triage, and Roland PLM branch are intentionally left behind.

## What it does

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

## Layout

```
walker/     logan-walker: the native Rust engine (the product)
python/     disassemble.py (end-to-end driver) + build_walker_codon_msa.py (codon MSA + FastTree)
queries/    reference coding sequences (pfcrt, hiv1_pol, 25 planetary markers, master panel)
bench/      golden fixtures, the check_golden.sh concordance gate, k-mer sensitivity sim (archive/ holds retired Python-walker tooling)
docs/       codon MSA/phylogeny step, VCF example, contigs-vs-unitigs, the variance-loss investigation, upstream context
```

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
