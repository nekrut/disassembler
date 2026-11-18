# Documentation

`disassembler` reconstructs one gene's coding sequence from public sequencing data, either by
walking Logan's unitig graphs or by tiling a LexicMap search, and then aligns the cohort in
frame. These pages document what each stage does, what it writes, and the measurements the
defaults rest on. The project overview is in the [top-level README](../README.md).

## Start here

| if you want to | read |
|---|---|
| reconstruct haplotypes from SRA accessions | [`walker.md`](walker.md) |
| look up a walker flag, output file or exit code | [`walker_reference.md`](walker_reference.md) |
| find candidate accessions with kmindex first | [`kmindex.md`](kmindex.md) |
| work from a LexicMap search instead of unitig graphs | [`lexicmap.md`](lexicmap.md) |
| thin a cohort of near-duplicates before the tree | [`sieve.md`](sieve.md) |
| align and build a tree, or run the whole pipeline | [`codon_msa.md`](codon_msa.md) |
| parse any stage's output programmatically | [`summaries.md`](summaries.md) |
| work out why a run failed | [`troubleshooting.md`](troubleshooting.md) |

A first pass through the suite reads naturally as `walker.md`, then `summaries.md`, then
`codon_msa.md`. Everything else is reference.

## The stages

- [`walker.md`](walker.md): what `logan-walker` does, how to run it, multi-anchor mode, the
  offline cache, and what its output can and cannot be read as.
- [`walker_reference.md`](walker_reference.md): every flag, every output file, the exit
  statuses, and the memory and disk model.
- [`kmindex.md`](kmindex.md): `harvest`, `stratify` and `from-kmindex`.
- [`lexicmap.md`](lexicmap.md): the LexicMap route, its six tiers and `from-lexicmap`.
- [`sieve.md`](sieve.md): diversity thinning, and how to choose its threshold.
- [`codon_msa.md`](codon_msa.md): the codon MSA and tree, and `reconstruct`.
- [`summaries.md`](summaries.md): the summary JSON contract every stage writes, and its
  schema versioning.

## Operations

- [`troubleshooting.md`](troubleshooting.md): failure modes by symptom.
- [`galaxy.md`](galaxy.md): the Galaxy tool wrappers.
- [`packaging.md`](packaging.md): the three Bioconda packages, the pinned tools, and the
  release checklist.

## Examples

- [`examples/vcf_example.md`](examples/vcf_example.md): intra-host SNVs on a Zika NS5 cohort,
  from the walker's VCFs to a site-by-accession matrix, with the data under
  [`examples/zika_ns5/`](examples/zika_ns5/).

## Notes

Investigations that explain a default or a pin. Each one opens with what it found, when, and
whether the behaviour it describes is still current. They are evidence, not instructions.

- [`notes/contigs_vs_unitigs.md`](notes/contigs_vs_unitigs.md): why the walker streams unitigs
  rather than Logan's assembled contigs.
- [`notes/intra_accession_variance.md`](notes/intra_accession_variance.md): why path walking
  can take the reference-matching minority allele, and what the pruning and majority-preference
  defaults do about it.
- [`notes/graphaligner_end_clipping.md`](notes/graphaligner_end_clipping.md): why GraphAligner
  is pinned to 1.0.17b.
- [`notes/sieve_calibration.md`](notes/sieve_calibration.md): why the sieve's wrapper defaults
  differ from the binary's own.
- [`notes/galaxy_plan.md`](notes/galaxy_plan.md): the Galaxy workflows that are planned but not
  built.

Related documentation elsewhere in the repository: [`../bench/README.md`](../bench/README.md)
for the golden gate and the benchmarks, [`../walker/README.md`](../walker/README.md) for the
walker's source layout, [`../queries/README.md`](../queries/README.md) for the shipped marker
CDSs, and [`../panels/README.md`](../panels/README.md) for the multi-anchor panels.
