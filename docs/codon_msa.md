# Codon alignment and phylogeny

The walker reconstructs one coding haplotype per accession. To compare them — for a tree, or as
input to a selection analysis — they have to be aligned **in frame**, so that every column is a
codon position and gaps are whole-codon indels rather than frameshifts. That is what this step
does, in `python/build_walker_codon_msa.py`, and `python/disassemble.py` runs it straight after
reconstruction.

## What it produces

From the reconstructed haplotypes (`reconstructed.fasta`) plus the reference CDS:

- **`<gene>_codon_msa.fasta`** — an in-frame codon multiple sequence alignment, reference first.
- **`<gene>_fasttree.nwk`** — a FastTree GTR maximum-likelihood tree over that alignment.

## How it stays in frame

- **MAFFT `--keeplength --addfragments`.** The reference CDS is the fixed backbone; each
  reconstructed haplotype is added against it. `--keeplength` holds the alignment to the
  reference's length, so reference coordinates `1..L` are preserved and any insertion relative
  to the reference is trimmed rather than shifting the frame. A haplotype one base short of the
  reference (a common `path_end` trim) becomes a single terminal gap, not a frameshift.
- **Reference is sequence #1.** Downstream tools (FastTree here, HyPhy in the parent project)
  read the alignment expecting the reference row first; the builder guarantees it.
- **Header sanitisation.** Record ids are reduced to alphanumeric/underscore so newick and
  HyPhy parsers do not choke on them.
- **Coverage gate.** A haplotype below `--min-cov` (default 0.70) of the reference length is
  dropped before alignment rather than aligned as a fragment.

The builder prints QC on the way through: how many sequences passed the coverage filter, the
final taxa count, alignment length, and codon count. A clean run reports zero frameshifts and,
for a coding reference, the alignment length is an exact multiple of three.

## Running it

End to end, from accessions:

```bash
python/disassemble.py -a SRR6251967 -q queries/markers/zika_ns5.fasta -o out/zika
#   out/zika/reconstructed.fasta        the haplotypes
#   out/zika/zika_ns5_codon_msa.fasta   the codon alignment
#   out/zika/zika_ns5_fasttree.nwk      the tree
```

Or the MSA step alone, on an existing `reconstructed.fasta`:

```bash
python/build_walker_codon_msa.py \
  --ref queries/markers/zika_ns5.fasta \
  --input out/zika/reconstructed.fasta \
  --out-aln out/zika/zika_ns5_codon_msa.fasta \
  --out-tree out/zika/zika_ns5_fasttree.nwk
```

Needs `MAFFT`, `FastTree`, and biopython — all in the `disassembler` conda env
(`environment.yml`). `disassemble.py` also needs `GraphAligner` for the walker half.

## A worked cohort

On the nine-accession Zika NS5 cohort, this step produced a 10-taxon codon alignment
(reference + nine primaries) of 903 codons with zero internal stops, and the FastTree tree
collapsed the identical sequences into zero-length clusters — the expected shape for a
consensus-clonal set. The full walkthrough, including the intra-host SNP side of the same
cohort, is in [`vcf_example.md`](vcf_example.md).

## What it is not

This step aligns and trees the **primary** consensus per accession. The intra-host alternatives
(`alt_haplotypes.fasta`) are deliberately kept out of it: each would enter as an independent
taxon and inflate apparent diversity, and two records from one sample are not two independent
observations. The intra-host signal has its own outputs (`--emit-vcf`, `--emit-snv-table`); see
[`vcf_example.md`](vcf_example.md).
