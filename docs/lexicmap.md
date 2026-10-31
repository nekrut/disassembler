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
   identity wins. Reference positions that no HSP covers are written as `-`.
2. **Reads the indels.** A run of `-` in `qseq` is an insertion in the sample, and a run of `-`
   in `sseq` a deletion. A deletion stays in the row as `-` and counts as covered; an insertion
   has no position on the grid, so it is recorded instead (see [Outputs](#outputs)). An indel
   counts only where its own HSP supplied the row, and one within 10 aligned bases of either end
   of its HSP is reported as `(edge)` but not trusted: aligners can close an alignment with a
   forced gap there. Indels LexicMap did not align through (a contig break, or an insertion too
   long to bridge) split the alignment into separate HSPs and cannot be seen.
3. **Checks the reading frame.** The trusted indels are walked in reference order, tracking how
   far out of frame the sample is (their net length, mod 3). A **frameshift** is a stretch that
   stays out of frame for at least 30 reference bases, or to the end of the gene. Indels that
   restore the frame sooner (for example +1 and then -1 five bases later, a common way for an
   aligner to write a divergent stretch) are treated like an in-frame indel. For the stop scan,
   every indel that is not part of a frameshift is applied to the tiled sequence (deletions
   removed, insertions put back), which is then translated with genetic code 11 and its coding
   body (everything before the final stop codon) scanned for premature stops, so a stop formed
   across an indel, or inside a briefly shifted stretch, is found. Positions are reported as
   reference codon numbers.
4. **Sorts accessions into tiers** (`status` and `tier` columns of the ledger):

   | tier | status | rule |
   |---|---|---|
   | 1 | `CLEAN_FULL_LENGTH` | coverage at least `--min-cov` (100% under `--require-full-length`), at most `--max-internal-stops` premature stops, no frameshift |
   | 2 | `CLEAN_PARTIAL` | coverage at least `--min-cov-partial`, same stop and frameshift rules |
   | 5 | `FLAGGED_FRAMESHIFT` | a frameshift (out of frame for 30 or more bases); the indels are listed in the `indels` column |
   | 3 | `FLAGGED_PREMATURE_STOPS` | more than `--max-internal-stops` premature stops, with their codon positions |
   | 4 | `LOW_COVERAGE_FRAGMENT` | coverage below `--min-cov-partial` |

   The rules are applied in that order, so an accession with a frameshift is filed under it
   whatever its coverage. `--allow-frameshift-alts` waives the frameshift rule; the stop rule
   still applies. Frameshifts are common in practice and mostly not biology: hits from
   divergent relatives, 1-nt deletions in homopolymers in Nanopore and PacBio runs, and, more
   rarely, real frameshift alleles. `flagged.tsv` lists them with their positions. The frame is
   judged against the reference, so a reference that is itself a pseudogene makes every intact
   sample look frameshifted; use an intact allele as the reference.
5. **Collapses haplotypes.** Tier 1 sequences are collapsed into unique haplotypes as they
   stream (memory grows with the number of distinct haplotypes, not the number of rows). For
   each haplotype it records the count and frequency, nucleotide (`198G>A`) and amino acid
   (`V66V`, `I61V`) changes against the reference, its indels (`685_686insT`, `2316_2321del`), a
   mutation class (`WT`, `WT_partial`, `Synonymous`, `Missense`, `Nonsense`, `InFrameIndel`,
   `Frameshift`, and `NonCoding_WT_partial` / `NonCoding_Variant` for a non-coding reference;
   the class is the most severe change), and up to `--samples-per-haplotype` representative
   accessions. A haplotype is its row plus its trusted indels, so an accession that differs
   from the wild type only by an insertion is its own haplotype even though its row is
   identical. It also reports the Shannon entropy of the haplotype distribution.
6. **Diagnoses the regime.** Coverage and identity quantiles (p10, median, p90) and stop-codon
   rates go into the summary. If the median identity is below 85% or the median coverage below
   0.60, the target is called **divergent**: reference-guided projection fragments there, and
   the log recommends the walker route instead.

## Input

A LexicMap search table, run with all columns on (`lexicmap search --all`, or *Output all
columns* in Galaxy) so that `qseq` and `sseq` are present. It can come from:

- a local file or a pipe (`-f`, including `-f <(zcat hits.tsv.gz)`);
- standard input (`-f -`);
- a Galaxy dataset streamed over the API (`-d <dataset id>`, with `--base-url`). The API key
  comes from `--api-key`, then `--api-key-file`, then `$GALAXY_API_KEY`, then `/tmp/gxy.txt`.
  Only a dataset in state `ok` is read. One still being produced (`queued`, `running`, `new`,
  `paused` and the like) exits 75 (`EX_TEMPFAIL`), so a script can try again later; one in
  any other state (`error`, for example) exits 1.

Any of them can be plain or gzip-compressed; gzip is detected from the content, not the name.
A zstd-compressed table is refused with a message, so decompress it first
(`zstd -dc hits.tsv.zst | disassembler lexicmap -f - ...`). Before 0.1.1, gzip was only
recognised in a local file, a pipe lost its contents to that check, and a queued Galaxy dataset
exited 0 having written nothing.

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
a path ending in `/`, writes `<dir>/<reference stem>.*`, as `--out-dir <dir>` does whether or not
the directory exists yet. (Before 0.1.1, `--out-dir` was another spelling of `-o`, so a directory
that did not exist yet was taken as a prefix.) `from-lexicmap -o DIR` writes the lexicmap files as
`DIR/lexicmap.*` and the MSA and tree alongside.

## Outputs

| file | contents |
|---|---|
| `<prefix>.clean.accessions.fasta` | one row per tier 1 accession |
| `<prefix>.clean_expanded.accessions.fasta` | tiers 1 and 2 together |
| `<prefix>.clean.msa.fasta` | one row per collapsed tier 1 haplotype, with rank, count and frequency |
| `<prefix>.clean.haplotypes.tsv` | the haplotype catalogue: changes, class, representative accessions |
| `<prefix>.clean.a2m` | the reference plus the collapsed haplotypes in A2M format, for protein language models and EVcouplings (`--no-a2m` skips it) |
| `<prefix>.flagged.accessions.fasta` | tiers 3, 4 and 5 |
| `<prefix>.flagged.tsv` | tiers 3, 4 and 5, with stop codon positions and indels |
| `<prefix>.cohort_ledger.tsv` | every accession: coverage, identity, status, tier, stops, HSPs merged, indels |
| `<prefix>.summary.json` | quantiles, tier counts, multi-HSP rescue counts, indel counts, entropy, regime |

`<prefix>.msa.fasta`, `<prefix>.haplotypes.tsv` and `<prefix>.a2m` are compatibility copies of
the `clean.*` files.

Every row in the FASTA outputs is on the reference grid, so the files are already
column-aligned. Under `--require-full-length`, the only `-` in a tier 1 row are in-frame
deletions. With the default `--min-cov 0.80`, a tier 1 row can also be up to 20% uncovered,
and a tier 2 row more.

Indels are written in a 1-based, HGVS-like notation: `685_686insT` (T inserted between
reference positions 685 and 686), `2316_2321del`, `709del`, with `(edge)` appended to an
untrusted one. They appear in the `indels` column of the ledger and `flagged.tsv`, as
`indels=...` in the FASTA headers, and among the `nt_mutations` of the haplotype table. Because
an insertion has no column, two haplotypes in `clean.msa.fasta` can have identical rows; their
headers and the haplotype table tell them apart.

For plots weighted by how many accessions carry each variant, use `clean.accessions.fasta` (one
row per accession), not `clean.msa.fasta` (one row per haplotype, whatever its count).

## `from-lexicmap`

Runs `lexicmap`, then (with `--sieve`) the sieve, then `msa`, on
`lexicmap.clean.accessions.fasta`, or on `lexicmap.clean_expanded.accessions.fasta` with
`--include-partials`. It takes every `lexicmap` option, the `--sieve-*` options (see
[`sieve.md`](sieve.md)), and `--mafft-threads` / `--msa-min-cov` for the MSA stage (see
[`codon_msa.md`](codon_msa.md#lexicmap-cohorts)). The sieve's gap limit and the MSA's coverage
gate default to values derived from the cohort's own coverage gate, so neither stage drops rows
for coverage that lexicmap already accepted.

## Options

| flag | default | meaning |
|---|---|---|
| `-q, -r, --query, --ref <FASTA>` | | reference in-frame CDS, plain or gzipped (first record only) |
| `-d, --dataset <ID>` | | Galaxy dataset to stream; exclusive with `-f` |
| `-f, -i, --file, --input <FILE>` | | local table, gzipped table, or `-` for stdin |
| `-o, --out-prefix <PATH>` | `msa/<reference stem>` | output prefix, or an existing directory |
| `--out-dir <DIR>` | | output directory, created if needed; exclusive with `-o` |
| `--base-url <URL>` | `https://usegalaxy.org` | Galaxy server for `-d` |
| `--api-key <KEY>` / `--api-key-file <PATH>` | | Galaxy API key, literal or from a file |
| `--min-cov <F>` | `0.80` | tier 1 coverage gate |
| `--min-cov-partial <F>` | `0.50` | tier 2 coverage gate; must not exceed `--min-cov` |
| `--require-full-length` | off | tier 1 requires 100% coverage instead of `--min-cov` |
| `--min-pident <F>` | `60.0` | minimum identity for an HSP to be tiled |
| `--max-internal-stops <N>` | `0` | premature stops tolerated in the clean tiers |
| `--allow-frameshift-alts` | off | let accessions with a frameshift into the clean tiers instead of `FLAGGED_FRAMESHIFT` |
| `--samples-per-haplotype <N>` | `10` | representative accessions per haplotype in the TSV |
| `--max-records <N>` | | stop after N input rows (for previews) |
| `--no-a2m` | off | skip the A2M output |
| `--allow-empty` | off | exit 0 when the clean cohort is empty (non-zero otherwise) |

For a strict clean cohort, use `--require-full-length --min-pident 95 --max-internal-stops 0`.
