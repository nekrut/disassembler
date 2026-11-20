# Codon alignment and phylogeny

The walker reconstructs one primary coding haplotype per accession. To compare them, for a tree
or as input to a selection analysis, they have to be put on a common coordinate system.
`disassembler msa` does this, and the `reconstruct`, `from-kmindex` and `from-lexicmap`
composites run it as their last stage.

**What the alignment is, precisely.** It is a nucleotide alignment anchored on the reference:
MAFFT runs with `--auto --keeplength --addfragments`, so every row comes out the reference's
length, every column is a reference position, and the reference's own frame therefore applies
to the whole file. That is not the same as a codon-aware alignment. `--keeplength` deletes any
base a row has beyond the reference, insertions included, and it puts gaps where nucleotide
scoring wants them, which need not be on codon boundaries. A row carrying a 1 nt indel can
come back looking intact, or with a gap that splits a codon. From 0.2.0 `msa` measures this
and says so; it does not change the alignment. See [codon QC](#codon-qc) below.

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
  and HyPhy parsers accept them. They are **not truncated**: ids used to be cut to 40
  characters, which is where a GISAID or NCBI header carries the accession that tells two
  records apart. An id that is already clean and unique keeps itself whatever order it
  arrives in. When any id changed, or any row was reverse-complemented,
  `<out-aln>.ids.tsv` records `aligned_id`, `original_header` and `reverse_complemented`.
- **Coverage gate, twice.** A haplotype with fewer A/C/G/T bases than `--msa-min-cov`
  (default `0.70`) of the reference length is dropped before alignment rather than aligned as
  a fragment. Gaps and N do not count, so a tiler row, which is reference-length with `-`
  wherever no hit reached, is measured by what it covers. It is **measured again on the
  aligned row**, with a two-point tolerance, because a contig carrying part of the gene plus
  flanking sequence passes the first gate on its raw base count and `--keeplength` then trims
  the flank away, leaving a fragment in the alignment and the tree with nothing said.

- **Orientation.** MAFFT does not reverse-complement its inputs, and `--keeplength` fits any
  input to the reference's length, so a haplotype on the other strand used to align as noise
  and exit 0. Each haplotype is now compared with the reference and its reverse complement by
  shared 15-mers, and reverse-complemented if it matches the reverse strand. Walker output is
  always in query orientation already; the case this catches is a FASTA assembled another way,
  such as a CDS annotated on the minus strand of a circular genome. On four real cohorts every
  row shared at least 896 15-mers with its reference's strand and at most 23 with the other.
  Exact 15-mers find nothing in common once a row is about 30% divergent; when they tie, all
  11-mers are compared instead and the call is made only with a clear margin, so a row that
  really cannot be placed is still reported as such rather than aligned as noise.
- **Identity check.** After MAFFT, each row's identity to the reference is measured over the
  columns where both have a base. Same-gene cohorts measured 0.88 or more (HIV-1 *pol*,
  SARS-CoV-2 spike, C. auris FKS1, H. pylori *vacA*); a row on the wrong strand, or from another
  gene, aligns at 0.40 to 0.50. Rows under `--min-identity` (default `0.60`) are named in a
  warning, and if the median is under it, `msa` fails without writing the alignment. `0`
  turns the check off.

`msa` prints QC on the way through: how many sequences passed the coverage gate, how many were
reverse-complemented, the median and lowest identity to the reference, the final taxa count,
the alignment length and the codon count. For a coding reference the alignment length is an
exact multiple of three. MAFFT and FastTree are found on `PATH`.

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

`msa` aligns against the first record of `--ref`, and warns when the file holds more. Both
`--ref` and `--input` can be plain or gzip-compressed (recognised from the content).

## `disassembler reconstruct`

`reconstruct` chains `walk`, an optional `sieve` (see [`sieve.md`](sieve.md)) and `msa`. It
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
| `--skip-tree` | off | stop after the walker; no sieve, MSA or tree (and an MSA or tree an earlier run left in `-o` is removed) |

Any other walker flag goes after a literal `--`, for example
`disassembler reconstruct ... -- --primary-metric pident --hops 7`. `reconstruct` sets the
walker's `-f` (so later stages can find the FASTA) and keeps `-t` at 1 (GraphAligner is not
reproducible with more threads).

Five flags are refused there, because each makes the walker exit 0 without writing a
reconstruction for the next stage to align: `--fetch-only`, `-h`, `--help`, `-V` and
`--version`. Use `disassembler fetch` or `disassembler walk --fetch-only` to warm a cache,
and `disassembler walk -h` or `--version` for the walker's own help. The composite also
removes the reconstruction before the walker runs, so a walker that writes nothing fails the
run rather than leaving an earlier run's cohort to be aligned as this one's.

With more than one reference record, from several `-q` or one multi-record file, the walker
runs in multi-anchor mode, and `reconstruct` says which record the MSA aligns against (the
first record of the first `-q`, unless `--msa-reference` names another). The sieve's default
`--min-len` is taken from the shortest reference record, because each primary is reconstructed
against its own anchor; see [`sieve.md`](sieve.md).

Each run removes `reconstructed.sieved.fasta` from an earlier run, since it no longer matches the
reconstruction the walker has just rewritten; with `--sieve`, this run's is written again.
`from-kmindex` ends in `reconstruct`, so the same applies there.

### When some accessions error

The walker exits 4 when some accessions errored and others reconstructed (see
[`walker_reference.md`](walker_reference.md#exit-status)). `reconstruct` does not stop there.
The haplotypes that were produced are real, so the sieve, MSA and tree all run on them, and it
exits 4 itself. Discarding a whole cohort's alignment and tree because one accession in
thousands hit a network error would contradict the walker's own rule that one bad accession
does not stop the run, and the run still reports failure, which is the part a scheduler reads.

The accessions that failed are named in the walker's `<query stem>.skipped.tsv`, and
`reconstruct` prints its path when it continues. Any other non-zero walker exit does stop the
run at stage 1, because 1 means nothing was reconstructed and 3 means nothing worked at all,
so in both cases there is nothing to carry forward.

When a stage does stop the run, `reconstruct` exits with **that stage's own code**, not a
flat 1: the walker's 3 and 1 say different things, `msa` returns 2 for a usage error as
argparse does, and a missing binary is 127. The one exception is a code that is not positive,
which means the stage did not normalise a signal death; that becomes 1, because the shell
would otherwise mask it to an unrelated byte.

`--msa-min-cov` gates only the MSA stage. The walker has its own `--min-cov` (default `0.85`)
that decides which reconstructions exist at all, and `reconstruct` never sets it: lowering
`--msa-min-cov` cannot admit haplotypes the walker already discarded. To move the walker's gate,
pass it through with `-- --min-cov <F>`. The same split applies to `--mafft-threads` (MAFFT)
versus the walker's `-t` (GraphAligner).

## Tiler cohorts

`disassembler tile` writes its clean cohort on the reference's coordinate grid already
(see [`tiler.md`](tiler.md#outputs)): every base sits at the reference position its HSP
gives it, and `-` marks a position no HSP reached. **That cohort is an alignment, so `msa`
takes it as one.** `from-lexicmap` passes `--aligned`, and `msa --aligned` writes the rows out
unchanged rather than calling MAFFT at all. The run is about twice as fast, and
`--mafft-threads` is inert there (the run says so).

Running MAFFT over those rows was wrong, not merely redundant. It discards the input's gaps and
places the bases again, and wherever several placements score the same it can slide bases
across an uncovered stretch or a deletion. On a reference HIV-1 pol with two gap-free HSPs
around an internal hole, a `C1416G` call came out of the MSA as `T1424C`, a different codon and
a different amino acid, while the haplotype table still reported `C1416G`, with exit 0 and the
identity check passing. In synthetic 200-row cohorts, 27 to 39% of rows with an internal hole
changed; rows truncated only at their ends never did.

A row that is not reference-length is an error under `--aligned`, not something to re-place
silently: the mode's whole premise is that the input grid is the output grid. Every other check
is unchanged, including the coverage gate, the identity check, the codon QC and the tree.

With `--require-full-length` every row is reference-length and its only gaps are in-frame
deletions (whole codons of `-`). Such rows did come through MAFFT unchanged, character for
character, on three cohorts of 27 to 14,594 sequences, which is why this went unnoticed: the
cohorts where it bites are the ones with an internal hole, which is what tiling is for.

In `from-lexicmap`, `--msa-min-cov` defaults to `0.70` or the cohort's coverage gate minus
0.05, whichever is lower (`0.45` with `--include-partials`), so the MSA keeps the rows the tiler
admitted. The gate measures covered bases rather than row length, which every tiler row
passes by construction.

## A worked cohort

On the nine-accession Zika NS5 cohort, this step produced a 10-taxon codon alignment (reference
plus nine primaries) of 903 codons with zero internal stops, and FastTree collapsed the
identical sequences into zero-length clusters, the expected shape for a consensus-clonal set.
The intra-host SNP side of the same cohort is in
[`examples/vcf_example.md`](examples/vcf_example.md).

## What it is not

`msa` aligns the **primary** haplotype of each accession. The intra-host alternatives
(`alt_haplotypes.fasta`) are deliberately kept out: each would enter as an independent taxon and
inflate apparent diversity, and two records from one sample are not two independent
observations. The intra-host signal has its own outputs (`--emit-vcf`, `--emit-snv-table`); see
[`examples/vcf_example.md`](examples/vcf_example.md).

## Summary output

`disassembler msa` writes `<out-aln stem>.summary.json` beside the alignment, in the shared
format of [`summaries.md`](summaries.md). It is new in 0.2.0; before it this stage produced no
machine-readable output, which was awkward for the one stage that can fail on a tool version
(MAFFT 7.525 mishandles `--addfragments`).

`discards.rows` carries `below_coverage`, `unoriented` and `lost_in_alignment`, `cohort`
carries the rows in, the rows retained, the taxa aligned and how many were reverse
complemented, and `records[]` names each aligned row by id, length and sequence hash. It is
written on every path, including the failures, so a caller can tell a run that aligned nothing
from one that never started.


## Codon QC

`--keeplength` makes every row the reference's length, which is the only thing the old check
tested, so a row that had been changed in ways that matter still passed. After alignment each
row is now compared with what went in, and anything of these kinds is reported on stderr and
written to `<out-aln>.qc.tsv`:

| column | meaning |
|---|---|
| `trimmed_bases` | bases `--keeplength` deleted: an insertion, or overhang past the reference ends |
| `internal_gap_runs` | the lengths of gap runs that touch neither end of the row |
| `split_codons` | codons where the row has one or two gaps, not three |
| `internal_stops` | in-frame stop codons before the last codon |
| `frameshift` | the trimmed count, or an internal gap run, is not a multiple of 3 |

The alignment is unchanged either way. A frameshifted row is a real finding on walker output,
where the walker gates frameshifts only for alternative alleles, so the primary of an accession
carrying a 1 nt indel reaches `msa` as it is; `from-lexicmap` filters those into its own tier 5
before this stage sees them.

## Tools and versions

The MAFFT and FastTree versions are printed at the start of a run and recorded in the summary.
MAFFT 7.525 is refused outright: it segfaults on the `--addfragments` path this stage uses, and
the published recipe leaves MAFFT unpinned, so an environment can resolve it. Install
`mafft>=7.526`.

`--no-tree` writes only the alignment, and FastTree is then neither looked up nor run, so an
alignment can be produced on a machine that does not have it. Exactly one of `--out-tree` and
`--no-tree` is required.

Duplicate rows are aligned once and copied back: under `--addfragments` each row is aligned to
the reference alone, so every copy of a sequence gets the same row. The identity median is still
weighted by copy number. A cohort of 600 rows drawn from 30 variants took 129.7 s before and
8.4 s after, for a byte-identical alignment.

A `--out-aln` or `--out-tree` that is not a plain file (a FIFO, a device such as `/dev/null`, or
a symlink) is **written into**, not replaced. A rename swaps a directory entry, so these outputs
used to turn into regular files: as root, which is the default inside a container, `--out-tree
/dev/null` left the null device holding a Newick tree for the rest of that container's life.
