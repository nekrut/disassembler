# Codon alignment and phylogeny

The walker reconstructs one primary coding haplotype per accession. To compare them, for a tree
or as input to a selection analysis, they have to be aligned **in frame**: every column a codon
position, and gaps whole codons rather than frameshifts. `disassembler msa` does this, and the
`reconstruct`, `from-kmindex` and `from-lexicmap` composites run it as their last stage.

## What it produces

From a FASTA of reconstructed haplotypes plus the reference CDS:

- **a codon MSA** (`--out-aln`): an in-frame alignment, reference first.
- **a tree** (`--out-tree`): a FastTree GTR maximum-likelihood tree (`FastTree -nt -gtr`) over
  that alignment, in Newick.

The composites name them `<gene>_codon_msa.fasta` and `<gene>_fasttree.nwk` in the output
directory, where `<gene>` is the reference file's stem.

## How it stays in frame

- **MAFFT `--keeplength --addfragments`.** The reference CDS is the fixed backbone, and each
  haplotype is added against it. `--keeplength` holds the alignment to the reference's length,
  so reference coordinates `1..L` are preserved and an insertion relative to the reference is
  trimmed rather than shifting the frame. A haplotype one base short of the reference becomes
  a single terminal gap, not a frameshift.
- **The reference is the first row**, named `Reference`, because downstream tools (FastTree,
  HyPhy) expect it there.
- **Record IDs are sanitised** to letters, digits and underscores, and made unique, so Newick
  and HyPhy parsers accept them.
- **Coverage gate.** A haplotype with fewer A/C/G/T bases than `--msa-min-cov` (default
  `0.70`) of the reference length is dropped before alignment rather than aligned as a
  fragment. Gaps and N do not count, so a LexicMap row, which is reference-length with `-`
  wherever no hit reached, is measured by what it covers.

`msa` prints QC on the way through: how many sequences passed the coverage gate, the final taxa
count, the alignment length and the codon count. For a coding reference the alignment length is
an exact multiple of three. MAFFT and FastTree are found on `PATH`.

**Orientation.** MAFFT does not reverse-complement its inputs. Every haplotype must be on the
same strand as the reference CDS. Walker output always is, because the walker reports
haplotypes in query orientation; a FASTA assembled some other way (for example a CDS annotated
on the minus strand of a circular genome) must be oriented before it reaches `msa`.

## Running it

End to end, from accessions:

```bash
disassembler reconstruct -a SRR6251967 -q queries/markers/zika_ns5.fasta -o out/zika
#   out/zika/reconstructed.fasta        the haplotypes
#   out/zika/zika_ns5_codon_msa.fasta   the codon alignment
#   out/zika/zika_ns5_fasttree.nwk      the tree
```

Or the MSA stage alone, on an existing FASTA:

```bash
disassembler msa \
  --ref queries/markers/zika_ns5.fasta \
  --input out/zika/reconstructed.fasta \
  --out-aln out/zika/zika_ns5_codon_msa.fasta \
  --out-tree out/zika/zika_ns5_fasttree.nwk
```

`msa` reads only the first record of `--ref`.

## `disassembler reconstruct`

`reconstruct` chains `walker`, an optional `sieve` (see [`sieve.md`](sieve.md)) and `msa`. It
needs `GraphAligner`, `MAFFT` and `FastTree` on `PATH`.

| flag | default | meaning |
|---|---|---|
| `-a` / `-l` | | one accession, or a file with one per line (exactly one of the two) |
| `-q, --query <FASTA>` | | reference CDS; repeatable, and more than one record switches the walker to multi-anchor mode (see [`walker.md`](walker.md#multi-anchor-mode)) |
| `--msa-reference <FASTA>` | first `-q` | reference for the MSA stage in multi-anchor mode |
| `-o, --out-dir <DIR>` | | output directory; the walker writes `reconstructed.fasta` here |
| `-w, --workers <N>` | `8` | walker compute workers |
| `--fetch-workers <N>` | `8` | walker S3 fetch workers |
| `--emit-vcf`, `--emit-snv-table`, `--no-emit-alts` | | walker variant outputs |
| `--prune-frac`, `--no-prefer-majority`, `--min-abund`, `--min-vaf`, `--max-size-mb` | walker defaults | passed to the walker |
| `--unitigs-cache <DIR>`, `--offline` | | walker cache (see [`walker.md`](walker.md#offline-cache)) |
| `--walker-bin <PATH>` | discovered | path to `logan-walker` |
| `--sieve` and `--sieve-*` | off | optional diversity thinning |
| `--mafft-threads <N>` | `8` | MAFFT threads |
| `--msa-min-cov <F>` | `0.70` | MSA entry gate |
| `--skip-tree` | off | stop after the walker; no sieve, MSA or tree |

Any other walker flag goes after a literal `--`, for example
`disassembler reconstruct ... -- --primary-metric pident --hops 7`. `reconstruct` sets the
walker's `-f` (so later stages can find the FASTA) and keeps `-t` at 1 (GraphAligner is not
reproducible with more threads).

`--msa-min-cov` gates only the MSA stage. The walker has its own `--min-cov` (default `0.85`)
that decides which reconstructions exist at all, and `reconstruct` never sets it: lowering
`--msa-min-cov` cannot admit haplotypes the walker already discarded. To move the walker's gate,
pass it through with `-- --min-cov <F>`. The same split applies to `--mafft-threads` (MAFFT)
versus the walker's `-t` (GraphAligner).

## LexicMap cohorts

`disassembler lexicmap` writes its clean cohort on the reference's coordinate grid already
(see [`lexicmap.md`](lexicmap.md#outputs)). With `--require-full-length`, every row is
reference-length and its only gaps are in-frame deletions (whole codons of `-`). Rows without a
deletion come through MAFFT unchanged: on three cohorts of 27 to 14,594 sequences, re-running
MAFFT reproduced them byte for byte. With the default `--min-cov 0.80`, rows can also carry `-`
at uncovered reference positions, and MAFFT may place gaps differently. `from-lexicmap` runs `msa` in both cases so that
both routes produce the same output shape.

In `from-lexicmap`, `--msa-min-cov` defaults to `0.70` or the cohort's coverage gate minus
0.05, whichever is lower (`0.45` with `--include-partials`), so the MSA keeps the rows lexicmap
admitted. Before 0.1.1 the gate measured row length, which every LexicMap row passes, so it
never filtered this route at all.

## A worked cohort

On the nine-accession Zika NS5 cohort, this step produced a 10-taxon codon alignment (reference
plus nine primaries) of 903 codons with zero internal stops, and FastTree collapsed the
identical sequences into zero-length clusters, the expected shape for a consensus-clonal set.
The intra-host SNP side of the same cohort is in [`vcf_example.md`](vcf_example.md).

## What it is not

`msa` aligns the **primary** haplotype of each accession. The intra-host alternatives
(`alt_haplotypes.fasta`) are deliberately kept out: each would enter as an independent taxon and
inflate apparent diversity, and two records from one sample are not two independent
observations. The intra-host signal has its own outputs (`--emit-vcf`, `--emit-snv-table`); see
[`vcf_example.md`](vcf_example.md).
