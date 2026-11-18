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
   and aligned `qseq`/`sseq`. Where HSPs overlap, each position is decided on its own, by
   depth: the base comes from the HSP the position lies furthest inside, measured to the nearer
   of that HSP's two ends, with percent identity only as a tiebreak. Where two HSPs place
   different bases at a position at least 10 bases inside both, the position is recorded as a
   **conflict** and the accession is flagged (tier 6). Reference positions that no HSP covers
   are written as `-`.

   Through 0.1.1 the rule was whole-HSP: the HSP with the higher identity *to the reference*
   won every position it touched. Logan leaves an unpopped SNP bubble as two separate contigs,
   both of which clear LexicMap's default `-l 50 -i 70`, and the reference branch is 100%
   identical while the alternative one is not, so the reference branch overrode the real
   variant regardless of which allele the sample carries. A deletion on the losing HSP was then
   dropped as "not owned", so a premature stop or frameshift disappeared with it and the
   accession came out `CLEAN_FULL_LENGTH` wild type with nothing in any column saying so. Depth
   is the right tiebreak because a bubble branch is short: every one of its positions is near
   one of its ends, while the same positions sit deep inside the contig that spans the gene.
2. **Reads the indels.** A run of `-` in `qseq` is an insertion in the sample, and a run of `-`
   in `sseq` a deletion. A deletion stays in the row as `-` and counts as covered; an insertion
   has no position on the grid, so it is recorded instead (see [Outputs](#outputs)). An indel
   counts only where its own HSP supplied the row, and one within 10 aligned bases of the end
   of **every** HSP that carries it is reported as `(edge)` but not trusted: aligners can close
   an alignment with a forced gap there. Adjacent Logan contigs overlap by k-1 = 30 nt, so an
   indel can sit at one contig's end and deep inside the next; judged from the first contig
   alone a real frameshift was marked `(edge)` and the frame rule then ignored it. An edge
   deletion also never hides a base another HSP supplied: its `-` is replaced by that base and
   the indel is dropped, because an alignment artefact is not a statement about the sample.
   Indels LexicMap did not align through (a contig break, or an insertion too long to bridge)
   split the alignment into separate HSPs and cannot be seen.
3. **Checks the reading frame.** The trusted indels are walked in reference order, tracking how
   far out of frame the sample is (their net length, mod 3). A **frameshift** is a stretch that
   stays out of frame for at least 30 reference bases, or to the end of the gene; a
   deletion's own bases do not count towards that span, because the sample does not have them
   and the frame shifts only after them. Indels that
   restore the frame sooner (for example +1 and then -1 five bases later, a common way for an
   aligner to write a divergent stretch) are treated like an in-frame indel. For the stop scan,
   every indel that is not part of a frameshift is applied to the tiled sequence (deletions
   removed, insertions put back), which is then translated with genetic code 11 and its coding
   body scanned for premature stops, so a stop formed across an indel, or inside a briefly
   shifted stretch, is found. Positions are reported as reference codon numbers.

   The genetic code is `--genetic-code` (default 11, bacterial and plant plastid). The shipped
   `queries/markers/candida_auris_erg11.fasta` is from a CTG-clade yeast and needs
   `--genetic-code 12`, where CTG is Ser and not Leu. The coding body is everything before the
   final stop codon **when the reference ends in one**; six shipped references are polyprotein
   segments that end in a sense codon (`dengue_1_ns5` to `dengue_4_ns5`, `westnile_ns5`,
   `zika_ns5`), and for those the last codon is scanned like any other.

   **The reference is checked before anything else runs.** A reference with its own in-frame
   stop under the selected code is refused, naming the codons: its stops would otherwise be
   counted against every sample, so even perfect matches were flagged and the run ended
   reporting an empty clean cohort, which names the wrong cause. A reference whose length is
   not a multiple of 3 is treated as non-coding, with a warning; pass `--non-coding` for a
   non-coding target whose length happens to be one, which is one target in three. **Both the
   stop rule and the frame rule apply only to a coding reference**, so for an rRNA, ITS or UTR
   marker an indel allele reaches the clean cohort as `NonCoding_Variant`.
4. **Sorts accessions into tiers** (`status` and `tier` columns of the ledger):

   | tier | status | rule |
   |---|---|---|
   | 1 | `CLEAN_FULL_LENGTH` | coverage at least `--min-cov` (100% under `--require-full-length`), at most `--max-internal-stops` premature stops, no frameshift, no conflicting HSPs |
   | 2 | `CLEAN_PARTIAL` | coverage at least `--min-cov-partial`, same stop, frameshift and conflict rules |
   | 4 | `LOW_COVERAGE_FRAGMENT` | coverage below `--min-cov-partial` |
   | 5 | `FLAGGED_FRAMESHIFT` | covered, but out of frame for 30 or more bases; the indels are listed in the `indels` column |
   | 3 | `FLAGGED_PREMATURE_STOPS` | covered and in frame, but more than `--max-internal-stops` premature stops, with their codon positions |
   | 6 | `FLAGGED_MIXED_HSPS` | covered, in frame and without premature stops, but the accession's own HSPs disagree about at least one base; the positions are in the `conflict_sites` column |

   The rules are applied in that order, so coverage decides first: an accession with too little
   of the gene tiled is a fragment whatever its frame and stop codons look like, and the two
   flagged tiers mean "covered well enough to judge, and broken". Through 0.1.1-dev the frame
   and stop rules came first, which made both flagged tiers mostly fragments: on the
   SARS-CoV-2 strict pass, 61,897 of the 64,034 accessions reported as frameshifted were
   fragments at a median coverage of 0.18, where a real shift cannot be told from a gap in the
   tiling. The conflict rule is judged last, so an accession with both a premature stop and
   conflicting HSPs is reported as the stop. `--allow-frameshift-alts` waives the frameshift
   rule and `--allow-mixed` the conflict rule, each waiving only its own tier; the stop rule
   still applies. Frameshifts are common in practice and mostly not biology: hits from
   divergent relatives, 1-nt deletions in homopolymers in Nanopore and PacBio runs, and, more
   rarely, real frameshift alleles. `flagged.tsv` lists them with their positions. The frame is
   judged against the reference, so a reference that is itself a pseudogene makes every intact
   sample look frameshifted; use an intact allele as the reference.
5. **Collapses haplotypes.** Tier 1 sequences are collapsed into unique haplotypes as they
   stream. Memory grows with the number of distinct haplotypes and, far more slowly, with the
   number of accessions: roughly 100 bytes each for the contiguity check, plus one sequence
   per distinct clean haplotype. It does **not** grow with the number of rows, which is the
   unbounded quantity. For
   each haplotype it records the count and frequency, nucleotide (`198G>A`) and amino acid
   (`V66V`, `I61V`) changes against the reference, its indels (`685_686insT`, `2316_2321del`), a
   mutation class (`WT`, `WT_partial`, `Synonymous`, `Missense`, `Nonsense`, `InFrameIndel`,
   `Frameshift`, and `NonCoding_WT_partial` / `NonCoding_Variant` for a non-coding reference;
   the class is the most severe change), and up to `--samples-per-haplotype` representative
   accessions. A haplotype is its row plus its trusted indels, so an accession that differs
   from the wild type only by an insertion is its own haplotype even though its row is
   identical. It also reports the Shannon entropy of the haplotype distribution.

   Amino-acid changes are read from the sample's own sequence, not from the reference grid, so
   a trusted in-frame indel is described: `P11del`, `K10_P11insA`, `R101_E105delinsKRVQQ`.
   Each cluster of nearby indels is described on its own, so two indels in different parts of
   the gene are two changes rather than one long `delins`. 3' shifting in a repeat is not
   applied, so a call inside a homopolymer can differ from what another tool prints. An `N` or
   other ambiguity code is missing information rather than a change and is not reported as
   one.

   **Without `--require-full-length` these counts depend on coverage.** A haplotype is the
   whole reference-length row, which carries `-` at every untiled position, and with the
   default `--min-cov 0.80` a tier 1 row can be a fifth `-`, so two accessions that agree
   wherever both are covered but whose contigs stop in different places are two haplotypes.
   `clean_unique_haplotypes`, `shannon_entropy_bits` and `wildtype_fraction` therefore measure
   assembly fragmentation as well as sequence diversity, and track sequencing depth. From
   0.2.0 the summary carries coverage-independent figures beside them (the `*_full_coverage`
   and `*_folded` fields below); the original fields are unchanged, so runs already done stay
   comparable.
6. **Diagnoses the regime.** Coverage and identity quantiles (p10, median, p90) and stop-codon
   rates go into the summary. Quantiles are interpolated linearly between the two neighbouring
   values (numpy's "linear"), so an even cohort's median sits between its two middle values
   rather than being one of them. If the median identity is below 85% or the median coverage
   below 0.60, the target is called **divergent**: reference-guided projection fragments there,
   and the log recommends the walker route instead.

   A search that found nothing is `no_data`, with null quantiles and no warning, rather than
   `divergent` with an invented median of 0.0%. A cohort whose HSPs all fell below
   `--min-pident` is divergent and says so in those words.

## Input

A LexicMap search table, run with all columns on (`lexicmap search --all`, or *Output all
columns* in Galaxy) so that `qseq` and `sseq` are present. It can come from:

- a local file or a pipe (`-f`, including `-f <(zcat hits.tsv.gz)`);
- standard input (`-f -`);
- a Galaxy dataset streamed over the API (`-d <dataset id>`, with `--base-url`). The API key
  comes from `--api-key`, then `--api-key-file`, then `$GALAXY_API_KEY`, then `/tmp/gxy.txt`.
  Only a dataset in state `ok` is read. One still being produced (`new`, `upload`, `queued`,
  `running`, `setting_metadata`: Galaxy's own non-ready states) exits 75 (`EX_TEMPFAIL`), so a
  script can try again later. A `deferred` or `paused` dataset exits 1 instead, naming what to
  do: neither becomes `ok` on its own, so retrying on 75 would loop until the caller's own
  limit. Any other state (`error`, for example) exits 1.

  A download that ends before the size the server declared exits 1; rerun it. So does a
  truncated or corrupt gzip, naming the file.

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
`DIR/lexicmap.*` and the MSA and tree alongside. `DIR/lexicmap` must therefore not itself be a
directory: lexicmap would write inside it as `DIR/lexicmap/<reference stem>.*` while the
composite went on reading `DIR/lexicmap.clean.accessions.fasta`, which with an earlier run's
file at that path meant aligning the old cohort. `from-lexicmap` refuses that case, and removes
the cohort files it is about to read before lexicmap runs.

## Outputs

| file | contents |
|---|---|
| `<prefix>.clean.accessions.fasta` | one row per tier 1 accession |
| `<prefix>.clean_expanded.accessions.fasta` | tiers 1 and 2 together |
| `<prefix>.clean.msa.fasta` | one row per collapsed tier 1 haplotype, with rank, count and frequency |
| `<prefix>.clean.haplotypes.tsv` | the haplotype catalogue: changes, class, representative accessions |
| `<prefix>.clean.a2m` | the reference plus the collapsed haplotypes in A2M format, **nucleotide, not protein**: translate it before using it with a protein MSA model or EVcouplings, which would otherwise accept a 4-letter "protein" of three times the real length and score it (`--no-a2m` skips it) |
| `<prefix>.flagged.accessions.fasta` | tiers 3, 4, 5 and 6 |
| `<prefix>.cohort_ledger.tsv` | **`--cohort-ledger` only.** Every accession: coverage, identity, status, tier, stops, HSPs merged, indels, conflict sites. A tabular rendering of the summary's `records[]` |
| `<prefix>.summary.json` | the run: provenance, settings, discards, cohort statistics, and one record per accession |

`<prefix>.cohort_ledger.tsv` is **off by default** from 0.2.0, because the summary's
`records[]` carries the same row for every accession and writing both put each one on disk
twice. Pass `--cohort-ledger` when something downstream wants tabular, such as a Galaxy
dataset or a spreadsheet; `disassembler from-lexicmap` takes the same flag and forwards it.

`<prefix>.flagged.tsv` was removed in 0.2.0. It held the flagged tiers with the same columns as
the ledger, so once every accession is in the summary's `records[]` with its status it was the
same rows written a third time. Filter `records[]` or the ledger on `status` instead.

The summary follows the contract in [`summaries.md`](summaries.md), which the walker's own
summary follows too, so one parser reads both. `tool`, `inputs` and `settings` are new in
0.2.0: before that this file recorded no provenance at all, not its own version, not which
dataset or file it read, and not one of the eight flags that change its answer, so two
summaries produced under different thresholds were indistinguishable on disk. `discards.rows`
is also new and carries `malformed`, `other_query` and `qlen_mismatch`, which were counted and
printed but never recorded. A summary is now written even when the run fails, so the absence of
the file is no longer the only signal; the one exception is a Galaxy dataset that is not ready
yet, which exits `EX_TEMPFAIL` without producing anything, because no run was attempted.

`cohort` also carries, from 0.2.0, figures that do not move when an accession's contigs stop
short:

| field | meaning |
|---|---|
| `full_coverage_rows` | tier 1 rows covering the whole reference (bases placed, plus the reference bases a trusted deletion removes) |
| `clean_unique_haplotypes_full_coverage` | distinct haplotypes among those rows |
| `shannon_entropy_bits_full_coverage` | entropy over those rows |
| `partial_rows_folded` | rows that stop short and match exactly one full-coverage haplotype at every position they do cover, carrying the same trusted indels |
| `partial_rows_unresolved` | rows that match none, or more than one: never guessed, and left out of the folded figures |
| `clean_unique_haplotypes_folded` | distinct haplotypes after folding |
| `shannon_entropy_bits_folded` | entropy after folding |
| `wildtype_compatible_count` / `wildtype_fraction_folded` | wild-type rows after folding, and their share |

A partial row is never folded on its variants alone, because a row whose variant lies in its
own uncovered region would then be called wild type.

`records[]` gains `n_conflicts` and `conflict_sites` (the first 10 positions, 1-based; the
count is always exact).

`<prefix>.msa.fasta`, `<prefix>.haplotypes.tsv` and `<prefix>.a2m` are compatibility copies of
the `clean.*` files.

**The outputs appear together or not at all.** Each is written to `<name>.part` and renamed in
only once the run has finished, and anything an earlier run left under a name this run did not
write is removed, so a `--no-a2m` rerun cannot leave the previous run's A2M beside a new MSA.
`summary.json` is the exception: a failed run still writes one, saying so, which is how a
reader is meant to tell.

Whitespace inside a sequence line of the reference FASTA is ignored, and nothing may precede
the first `>` header. Path separators in the record id become `_` when the id is used to name
output files; the id itself is unchanged, since it has to keep matching the search's query
column.

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
gate default to values derived from the cohort's own coverage gate, and from the cohort's worst
row, so neither stage drops rows for coverage that lexicmap already accepted. The MSA stage
runs with `--aligned`, because this cohort is already on the reference grid: the rows are
written into the alignment as they are, rather than being re-placed by MAFFT, which could move
a base across an uncovered stretch (see [`codon_msa.md`](codon_msa.md#lexicmap-cohorts)). That second floor
is what makes the promise true for a long in-frame deletion: lexicmap counts a trusted deletion
as covered with no limit, while the row keeps `-` there.

Each run removes the sieved cohorts an earlier run left in `-o`, since they no longer match the
cohorts lexicmap has just rewritten; with `--sieve`, this run's is written again. An empty cohort
fails the run unless `--allow-empty` is given. With it, `from-lexicmap` exits 0 without aligning
and also removes an earlier run's MSA and tree.

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
| `--min-cov <F>` | `0.80` | tier 1 coverage gate; ignored under `--require-full-length` |
| `--min-cov-partial <F>` | `0.50` | tier 2 coverage gate; must not exceed `--min-cov` |
| `--require-full-length` | off | tier 1 requires 100% coverage instead of `--min-cov` |
| `--min-pident <F>` | `60.0` | minimum identity for an HSP to be tiled |
| `--genetic-code <N>` | `11` | NCBI genetic code id for translation and the stop scan (the shipped C. auris ERG11 marker needs `12`) |
| `--non-coding` | off | treat the reference as non-coding: no stop scan, no frame rule |
| `--max-internal-stops <N>` | `0` | premature stops tolerated in the clean tiers |
| `--allow-frameshift-alts` | off | let accessions with a frameshift into the clean tiers instead of `FLAGGED_FRAMESHIFT` |
| `--allow-mixed` | off | let accessions whose own HSPs disagree about a base into the clean tiers instead of `FLAGGED_MIXED_HSPS` |
| `--samples-per-haplotype <N>` | `10` | representative accessions per haplotype in the TSV |
| `--max-records <N>` | | stop after N input rows (for previews); the accession in progress is finished, so the last row of a preview is not mis-tiered |
| `--no-a2m` | off | skip the A2M output |
| `--allow-empty` | off | exit 0 when the clean cohort is empty (non-zero otherwise) |

For a strict clean cohort, use `--require-full-length --min-pident 95 --max-internal-stops 0`.
