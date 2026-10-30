# The LexicMap route

`disassembler lexicmap` turns a [LexicMap](https://github.com/IndexThePlanet/LexicMap) search of
a reference CDS against Logan into a codon-checked cohort: one tiled sequence per accession,
sorted into four tiers, with the clean full-length sequences collapsed into unique haplotypes.
`disassembler from-lexicmap` adds the optional sieve and the codon MSA and tree.

This route is fast, because it never downloads unitig graphs, but it reconstructs a projection
of each accession onto the reference rather than walking a haplotype. For haplotypes, intra-host
variants and VCFs, use the walker route ([`walker.md`](walker.md)).

## Why tiling is needed

LexicMap indexes Logan's **contigs**, not the unitig graph the walker uses. On internally diverse
samples the assembler breaks a gene at every bubble it cannot pop, so the gene is spread across
many short contigs: on Zika NS5 it takes about 50 contig fragments, the longest under 100 nt, to
cover the gene (see [`contigs_vs_unitigs.md`](contigs_vs_unitigs.md)). LexicMap reports each
contig as its own HSP, so a gene split 55/45 across two contigs fails any single-row filter such
as `qcovHSP >= 80%`, even though together they cover all of it.

## What it does

1. **Tiles HSPs onto the reference.** Rows are grouped by accession (`sgenome`), and every HSP
   at or above `--min-pident` is projected onto reference positions `1..L` using its `qstart`
   and aligned `qseq`/`sseq`. Where HSPs overlap, the base from the HSP with the higher percent
   identity wins. Bases the reference does not have (insertions in the sample) have no position
   on the grid and are not kept; reference positions that no HSP covers, or that an HSP aligns
   to a gap, are written as `-`.
2. **Checks the reading frame.** The tiled sequence is translated with genetic code 11, and the
   coding body (everything before the final stop codon) is scanned for premature stops. Codons
   that contain a `-` are not checked.
3. **Sorts accessions into four tiers** (`status` column of the ledger):

   | tier | status | rule |
   |---|---|---|
   | 1 | `CLEAN_FULL_LENGTH` | coverage at least `--min-cov` (100% under `--require-full-length`) and at most `--max-internal-stops` premature stops |
   | 2 | `CLEAN_PARTIAL` | coverage at least `--min-cov-partial`, same stop rule |
   | 3 | `FLAGGED_PREMATURE_STOPS` | more than `--max-internal-stops` premature stops, with their codon positions |
   | 4 | `LOW_COVERAGE_FRAGMENT` | coverage below `--min-cov-partial` |

   `--allow-frameshift-alts` switches the stop rule off, so accessions are tiered by coverage
   alone.
4. **Collapses haplotypes.** Tier 1 sequences are collapsed into unique haplotypes as they
   stream (memory grows with the number of distinct haplotypes, not the number of rows). For
   each haplotype it records the count and frequency, nucleotide (`198G>A`) and amino acid
   (`V66V`, `I61V`) changes against the reference, a mutation class (`WT`, `WT_partial`,
   `Synonymous`, `Missense`, `Nonsense`, and `NonCoding_WT_partial` / `NonCoding_Variant` for a
   non-coding reference), and up to `--samples-per-haplotype` representative accessions. It also
   reports the Shannon entropy of the haplotype distribution.
5. **Diagnoses the regime.** Coverage and identity quantiles (p10, median, p90) and stop-codon
   rates go into the summary. If the median identity is below 85% or the median coverage below
   0.60, the target is called **divergent**: reference-guided projection fragments there, and
   the log recommends the walker route instead.

## Input

A LexicMap search table, run with all columns on (`lexicmap search --all`, or *Output all
columns* in Galaxy) so that `qseq` and `sseq` are present. It can come from:

- a local file (`-f`), plain or gzipped (gzip is detected from the content, not the name);
- standard input (`-f -`);
- a Galaxy dataset streamed over the API (`-d <dataset id>`, with `--base-url`). The API key
  comes from `--api-key`, `--api-key-file`, `$GALAXY_API_KEY` or `/tmp/gxy.txt`.

All rows for one accession must be contiguous in the input; `lexicmap` stops with an error
naming the accession if one reappears later.

Only the first record of `-q` is read, and its length defines the coordinate grid.

## Running it

```bash
# From a Galaxy dataset
disassembler lexicmap -q queries/markers/zika_ns5.fasta -d <DATASET_ID> -o out/zika_lexicmap

# From a local table (plain or gzipped)
disassembler lexicmap -q queries/query_hiv1_pol.fasta -f hiv1_pol_lexicmap.tsv.gz -o out/hiv1_pol

# From standard input
zcat hiv1_pol_lexicmap.tsv.gz | disassembler lexicmap -q queries/query_hiv1_pol.fasta -f - -o out/hiv1_pol

# Composite: lexicmap -> optional sieve -> codon MSA + tree
disassembler from-lexicmap -q queries/query_hiv1_pol.fasta -f hiv1_pol_lexicmap.tsv.gz \
    -o out/hiv1_pol --sieve
```

`-o` is an output prefix: `-o out/hiv1_pol` writes `out/hiv1_pol.*`. An existing directory, or
a path ending in `/`, writes `<dir>/<reference stem>.*`. `from-lexicmap -o DIR` writes the lexicmap files as
`DIR/lexicmap.*` and the MSA and tree alongside.

## Outputs

| file | contents |
|---|---|
| `<prefix>.clean.accessions.fasta` | one row per tier 1 accession |
| `<prefix>.clean_expanded.accessions.fasta` | tiers 1 and 2 together |
| `<prefix>.clean.msa.fasta` | one row per collapsed tier 1 haplotype, with rank, count and frequency |
| `<prefix>.clean.haplotypes.tsv` | the haplotype catalogue: changes, class, representative accessions |
| `<prefix>.clean.a2m` | the reference plus the collapsed haplotypes in A2M format, for protein language models and EVcouplings (`--no-a2m` skips it) |
| `<prefix>.flagged.accessions.fasta` | tiers 3 and 4 |
| `<prefix>.flagged.tsv` | tiers 3 and 4, with stop codon positions |
| `<prefix>.cohort_ledger.tsv` | every accession: coverage, identity, status, tier, stops, HSPs merged |
| `<prefix>.summary.json` | quantiles, tier counts, multi-HSP rescue counts, entropy, regime |

`<prefix>.msa.fasta`, `<prefix>.haplotypes.tsv` and `<prefix>.a2m` are compatibility copies of
the `clean.*` files.

Every row in the FASTA outputs is on the reference grid, so the files are already
column-aligned. Under `--require-full-length`, tier 1 rows are also gap-free. With the default
`--min-cov 0.80`, a tier 1 row can be up to 20% `-`, and a tier 2 row more.

For plots weighted by how many accessions carry each variant, use `clean.accessions.fasta` (one
row per accession), not `clean.msa.fasta` (one row per haplotype, whatever its count).

## `from-lexicmap`

Runs `lexicmap`, then (with `--sieve`) the sieve, then `msa`, on
`lexicmap.clean.accessions.fasta`, or on `lexicmap.clean_expanded.accessions.fasta` with
`--include-partials`. It takes every `lexicmap` option, the `--sieve-*` options (see
[`sieve.md`](sieve.md)), and `--mafft-threads` / `--msa-min-cov` for the MSA stage (see
[`codon_msa.md`](codon_msa.md#lexicmap-cohorts)).

## Options

| flag | default | meaning |
|---|---|---|
| `-q, -r, --query, --ref <FASTA>` | | reference in-frame CDS (first record only) |
| `-d, --dataset <ID>` | | Galaxy dataset to stream; exclusive with `-f` |
| `-f, -i, --file, --input <FILE>` | | local table, gzipped table, or `-` for stdin |
| `-o, --out-prefix, --out-dir <PATH>` | `msa/<reference stem>` | output prefix or directory |
| `--base-url <URL>` | `https://usegalaxy.org` | Galaxy server for `-d` |
| `--api-key <KEY>` / `--api-key-file <PATH>` | | Galaxy API key, literal or from a file |
| `--min-cov <F>` | `0.80` | tier 1 coverage gate |
| `--min-cov-partial <F>` | `0.50` | tier 2 coverage gate; must not exceed `--min-cov` |
| `--require-full-length` | off | tier 1 requires 100% coverage instead of `--min-cov` |
| `--min-pident <F>` | `60.0` | minimum identity for an HSP to be tiled |
| `--max-internal-stops <N>` | `0` | premature stops tolerated in the clean tiers |
| `--allow-frameshift-alts` | off | ignore the stop rule when tiering |
| `--samples-per-haplotype <N>` | `10` | representative accessions per haplotype in the TSV |
| `--max-records <N>` | | stop after N input rows (for previews) |
| `--no-a2m` | off | skip the A2M output |
| `--allow-empty` | off | exit 0 when the clean cohort is empty (non-zero otherwise) |

For a strict clean cohort, use `--require-full-length --min-pident 95 --max-internal-stops 0`.
