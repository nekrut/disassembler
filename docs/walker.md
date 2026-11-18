# logan-walker

`logan-walker` reconstructs the coding haplotypes of one gene in each of a list of SRA runs, from
Logan's unitig graphs. It is a Rust binary (`walker/`), packaged on Bioconda as `logan-walker`.
`disassembler walker` passes its arguments straight to it, and `disassembler reconstruct` and
`from-kmindex` run it as a stage.

## How it works

For each accession, given a reference CDS:

1. **Stream** the compressed unitig file,
   `https://logan-pub.s3.amazonaws.com/u/<ACC>/<ACC>.unitigs.fa.zst`, onto disk as it arrives:
   into `--unitigs-cache` when one is set, otherwise into an unlinked file under
   `--spool-dir`. It is never held in memory. An accession whose `Content-Length` exceeds
   `--max-size-mb` is refused before anything is downloaded.
2. **Parse** it in one streaming pass into a compact graph: adjacency, `ka:f` abundance and
   2-bit sequence for every unitig. The decompressed text is never held, so a 1.28M-unitig
   accession takes about 43 MB.
3. **Anchor** the reference by canonical k-mer containment (`-k 21`): a unitig sharing at least
   `--min-shared-kmers` k-mers with it is an anchor.
4. **Extract** the `--hops` (default 5) neighbourhood around the anchors as a local GFA.
5. **Align** the reference through that subgraph with `GraphAligner -x dbg`.
6. **Stitch** the walked path into the primary haplotype. By default the walker then prunes
   low-depth branches and re-aligns (`--prune-frac`), puts the majority allele on the primary
   (`--prefer-majority`, at two-allele sites and, since 0.2.0, wherever the aligner walked a
   minority branch whose other side holds an outright majority), and writes the remaining
   intra-host alleles as alternative haplotypes (`--emit-alts`). Why these exist is in
   [`issue_intra_accession_variance_lost.md`](issue_intra_accession_variance_lost.md).
7. **Write** the primaries, the alternatives, `bubbles.json`, a summary JSON, and optionally
   per-accession VCFs and a cohort SNV table.

Accessions are fetched and processed in parallel (`--fetch-workers` download, `-w` compute), and
outputs are written in input order whatever order the work finishes in.

## Running it

```bash
logan-walker -a SRR8209082 -q queries/query_hiv1_pol.fasta -o out -f out/hiv1_pol.fasta
disassembler walker -l accessions.txt -q queries/markers/zika_ns5.fasta -o out/zika --emit-snv-table
```

GraphAligner 1.0.17b must be on `PATH`, or named with `--graphaligner <PATH>` (see
[`issue_graphaligner_end_clipping.md`](issue_graphaligner_end_clipping.md) for why that version).
The walker checks it once at startup, before downloading anything, logs its version, and warns
if it is not 1.0.17b.
`disassembler` finds the walker through `--walker-bin`, then `DISASSEMBLER_WALKER_BIN`, then
`PATH`, then `walker/target/release/logan-walker` in a source checkout.

To reproduce the original one-haplotype-per-accession behaviour (the golden fixtures' setting),
use `--prune-frac 0.0 --no-prefer-majority --no-emit-alts`.

### The reference the aligner sees

The query file is read once, at startup. What GraphAligner is handed is the record the walker
parsed, written as `>{record id}` followed by one line of uppercased sequence, never a copy of
the file. The record id is the header's first token, whitespace inside the sequence is gone and
case is folded, so the aligner and the anchoring stage see the same bases by construction.

Until 0.2.0 a query file holding exactly one record was copied through byte for byte instead,
which cost three things, all of them quiet:

- a query that can only be read once, such as `-q <(zcat ref.fa.gz)` or `-q /dev/stdin`, was
  empty from the second accession onward, so GraphAligner got a 0-byte FASTA and every
  accession after the first came back `no_path`. A plumbing failure, reported as a clean
  biological negative
- the header line reached GraphAligner intact. 1.0.17b takes the whole line as the read name
  and writes it unescaped into GAF column 1, so a tab in `>gene<TAB>description` shifted every
  column and the row was dropped, and a non-UTF-8 byte anywhere in the header failed the
  accession outright
- whitespace inside sequence lines reached it too. 1.0.17b strips only trailing whitespace, so
  codon-spaced `ATG GCT AAA` left no usable seed and `no_path` everywhere, while GenBank-style
  10-base blocks reached the alignment's fallback branch, which calls `abort()`: SIGABRT on
  every accession

A query with more than one record was always written out the normalised way, so the last two
showed only where a single record arrived on its own: a one-record file, or the prune pass
realigning the multi-anchor winner out of a pool built from one-record files.

The first had a second cause, in startup rather than in the aligner. The walker read a single
`-q` file to count its records and, finding more than one, read it again to load the whole
reference set. On a pipe that second read saw an empty stream, so a multi-record query stopped
the run with `no FASTA records found`, naming a file that was perfectly good. Both reads are
now one, so `cat panel.fa | logan-walker -q /dev/stdin` works for any record count and
`README.md`'s "so does a pipe" is true of the query as it is of the other inputs.

### The alignment is checked against the graph

GraphAligner is handed a GFA the walker wrote, and the row that comes back is checked against
that same graph before anything is stitched from it: every segment must be one the subgraph
holds, every consecutive pair must be joined by an edge in the orientations the row gives, the
steps must stitch to the `path_len` the row declares, and the alignment window must lie inside
it. A row that fails is an `error` for that accession, not a primary.

This matters because stitching cannot notice. The walk is concatenation with the k - 1 overlap
removed, so a row naming two branches of the same fork in series produces a clean sequence that
exists in no organism, and it is published carrying the identity the row itself claimed: the
review built one with 89 differences from the truth at `pident=100.0`. The causes are a
different GraphAligner build, a wrapper that renames segments, or a truncated GAF. On the
second pass the same check applies to the pruned graph, where a failure refuses the prune and
keeps the first pass rather than losing the accession.

The path column is read strictly for the same reason: `>` or `<` then a unitig id, or GAF's `*`
for "no alignment". Anything else used to be skipped silently, which turned a row the walker
could not read into a shorter walk rather than into a question.

## What the primary is

`<ACC>_primary` is the **majority-consensus path** through the accession's subgraph, not a
phased haplotype, and the difference matters whenever a sample holds more than one strain.

Requiring a connected path is what keeps out joins that no k-mer supports, which is the defect
that makes concatenated BLAST hits unusable (87% chimeric artifacts, README section 4). It does
not phase anything. In a k = 31 graph, two allele sites more than 30 bases apart are joined
through unitigs both strains carry, so every combination of their alleles is a connected path,
and the walk chooses at each site independently. A sample carrying strain A at 60% over the 5'
half and strain B at 60% over the 3' half yields a primary that matches neither.

So the walker reports what the route passed through, in the FASTA header where the sequence is:

- `min_vaf=` is the lowest frequency of any allele this sequence carries, measured at its own
  fork. `vaf=` beside it is the constant 1.0 that primaries have always carried and is kept
  for one release; `min_vaf=` is the measured one.
- `mixed=1` appears when the route crosses **two or more** two-allele sites, which is the point
  at which the published sequence is a combination nothing in the data links. One such site is
  a variant, not a phasing problem, and is not flagged.
- A site counts only when its second allele is a **substitution** that rejoins the route. A
  length-changing branch in a de Bruijn graph is usually an artifact rather than a second
  population: 190 of 203 length-changing calls across the project's panel are plus or minus
  one base in a homopolymer, and the large ones are the route being sent through a different
  stretch of graph. `ERR2028820`, the clonal control gated in `bench/golden/`, has four
  branches clearing both frequency floors and all four are length-changing or never rejoin, so
  a looser rule would have labelled it a mixture.
- `summary.json` carries the same numbers per record as `variants.two_allele_sites` and
  `variants.min_site_vaf`, and `cohort.mixed` counts the accessions flagged.

Both are measured whatever the output flags say, including under `--no-emit-alts`: a mixture
flag that is only computed when alternatives were requested would read `mixed=0` on a sequence
nothing examined.

What to do with a flagged record depends on the question. For a selection analysis, excess
substitutions in a mixed record are not one organism's divergence, and
`alt_haplotypes.fasta` and the VCF carry the alleles it had to choose between.

## Multi-anchor mode

Give the walker more than one reference, by repeating `-q` or with a FASTA that holds several
records, and it runs in multi-anchor mode. This is meant for a panel of strains of the same gene,
such as `panels/sarscov2_spike_11strains.fa`:

- the anchor pool is the union of all references' k-mers, so one subgraph is extracted per
  accession;
- the reference set is aligned through that subgraph, giving one candidate haplotype per
  reference;
- one candidate becomes the primary, chosen by `--primary-metric`, and its header carries
  `anchor=<reference id>`. Every candidate's numbers are in the summary JSON under
  `"candidates"`, so a downstream tool can re-rank.

The choice matters for more than the tag: GraphAligner walks the graph differently against
each reference, so the primary's sequence depends on which reference wins.

**One gene per pool.** Multi-anchor mode ranks every reference against every other and keeps
one primary per accession, so a pool holding two different genes reports whichever scored
higher and leaves the other as an unselected candidate: an accession carrying both `bla_KPC`
and `bla_NDM` comes back with one of them and exit 0. The walker checks the pool at startup by
grouping references that share at least 1% of their canonical 21-mers, and warns when there is
more than one group, naming them. It warns rather than refuses, because a deliberately
divergent panel also splits and only the user knows which was meant. `queries/master_planetary_targets.fasta`
is a catalogue, not a panel: its 142 records fall into 31 groups, and it should be used one
gene at a time.

**Ambiguous reference bases never count as matched.** GraphAligner treats `N` as matching any
base, and each IUPAC code as matching every base it stands for, so a reference padded with `N`
would score as though it had matched bases it does not have, beat the genuinely nearest
reference, and then drive pruning, majority correction, the alternatives and every VCF
coordinate. From 0.2.0 the walker subtracts the ambiguous bases in the aligned window from the
match count before anything ranks on it. The correction is exact for an ambiguous base aligned
as a match and conservative otherwise, so it never overstates identity. This matters for the
shipped panel: in `queries/master_planetary_targets.fasta`, `ebola_gp_var2` is 42% ambiguous
and `zika_ns5_var3` is 24%.

**Tied references are decided by depth, not by panel order.** When a sample holds two strains
that each match a panel record, both candidates score the same and before 0.2.0 the record
listed first in the file won, so a 5% strain could be published as the primary while the 95%
strain left no trace: `--prefer-majority` cannot repair it, because it skips a fork whose
taken branch is below `--min-vaf`. References scoring within `--primary-tie-eps` (default
0.002, measured as a fraction of the metric) of the best are now treated as tied, and the
deepest route among them wins. Equal depths keep the old order, so a true tie still goes to
the first reference in the pool. `--primary-tie-eps 0` restores the pre-0.2.0 behaviour.

`--primary-metric` choices:

| metric | ranks by | trade-off |
|---|---|---|
| `matched` (default) | the fraction of the reference's own bases the alignment matched, then identity | never above 1.0; a reference loses for every base of its own the walk misses or mismatches, **and for every ambiguous base it carries** (see below) |
| `cov` | walked length / reference length, then identity | not capped at 1.0, so a reference shorter than the sample's gene (one with a deletion the sample lacks) beats an exact match |
| `pident` | identity, then coverage | can prefer a truncated walk that matches its reference cleanly over a full-length one |
| `cov_x_pident` | coverage x identity | inherits `cov`'s bias toward shorter references |

`matched` became the default in 0.1.1. With `cov`, the earlier default, a 20,000-accession
SARS-CoV-2 run on the 11-strain spike panel tagged 62% of primaries Omicron XBB.1.5 and 21%
JN.1, the two shortest references; `matched` gives 43% Delta and 26% BA.2, the cohort's real
mix. To reproduce output from a run made before 0.1.1, pass `--primary-metric cov`.

Every reference record needs a distinct id: the walker stops if two records, or two `-q`
files, share one.

VCF and `snvs.tsv` positions are in the coordinates of the reference each accession's primary
was reconstructed against: `CHROM` (VCF) and `gene` (`snvs.tsv`) name that reference. In a
multi-anchor run, compare sites across accessions only within one reference, or map them to a
common one first. The codon MSA needs a single reference; `reconstruct` uses the first `-q`
unless `--msa-reference` is given.

Two bugs in 0.1.0 depended on the reference FASTA's headers. GraphAligner names each result by
the whole header line, and the walker looked results up by the record id alone, so a header
with a description (`>pfcrt chloroquine resistance transporter`) never matched. A single-record
query with such a header came back `no_path` on every accession; in multi-anchor mode, the
re-alignment after pruning never matched for a panel with descriptions, such as
`panels/sarscov2_spike_11strains.fa`, so every accession kept its unpruned result. That
re-alignment also aligned the whole panel again rather than just the chosen reference. On 30
accessions of the 20,000-accession SARS-CoV-2 run, applying the pruning changed 17 primaries:
mostly by one or two bases, as intended, but one lost 425 nt. That one is why a pruned result
must now span as much of the reference as the unpruned one (see `--prune-frac` below).

## Outputs

| file | when | contents |
|---|---|---|
| `-f` (default `reconstructed_haplotypes.fasta`) | always | one primary per reconstructed accession, `>ACC_primary cov=… pident=… nodes=… ka=… vaf=… min_vaf=…` (plus ` mixed=1` when the route crosses two or more two-allele sites, and `anchor=` in multi-anchor mode) |
| `<out-dir>/alt_haplotypes.fasta` | unless `--no-emit-alts` | intra-host alternatives, `>ACC_altN …` with their VAF and fork, at most 32 per accession |
| `<out-dir>/bubbles.json` | always | two-allele sites in the local subgraph, one record per allele set with the flanks it was found between |
| `<out-dir>/<query stem>.summary.json` | always | the run's provenance, settings, totals, discards and cohort statistics, then one record per accession (below) |
| `<out-dir>/<query stem>.skipped.tsv` | always | one row per accession that produced no haplotype: `accession`, `outcome`, `bytes`, `detail`. Written even when empty, so it can be read without first testing whether it exists |
| `<out-dir>/<ACC>.vcf` | `--emit-vcf` | intra-host SNVs of one accession, written only when it has at least one (a VCF an earlier run left for an accession with none is removed) |
| `<out-dir>/snvs.tsv` | `--emit-snv-table` | one row per (accession, site, allele) across the run. Rows sharing a `block` come from one alternative haplotype, so they do co-occur |

Allele frequencies at a fork are computed over the branches that **return to the route**. An
exit that leaves the locus for good, a repeat, a paralog or a contaminant sharing the fork's
k-mer, is not an allele of the gene and is excluded from the denominator; counting it deflated
every frequency at the site and could suppress an allele entirely. `bubbles.json` has always
done this, so the two now agree.

An alternative carrying several substitutions is called site by site rather than discarded.
In a k=31 graph two substitutions within 30 nt share one bubble and arrive as one alternative,
and two changes in one codon are the common case; until 0.2.0 both sites vanished from the VCF
and the SNV table with no message. Those rows share a `BLOCK` (VCF) or `block` (TSV), which
says they **are** on one molecule, unlike two rows that merely share a position. Past
8 linked substitutions the alternative is a rerouted stretch of graph rather than a few linked
changes, so it stays in `alt_haplotypes.fasta` and is counted on stderr.

There is no length limit on the SNV mapping. A 20 kb cap used to drop every SNV, the VCF and
the SNV table for any longer gene (SARS-CoV-2 ORF1ab is 21,290 nt) while the run still exited
0, so the absence read as "no variants". The mapping is now bounded by memory instead, at
256 MiB, which refuses the case that is actually expensive, a short but highly divergent pair,
and says so.

Each run claims its scratch directory with a `.lock` it holds open, and the next run sweeps a
directory only when that lock names this machine and is free. Liveness used to be
`kill(pid, 0)`, which answers only about the caller's PID namespace on the caller's host, so a
run in another container or on another node sharing a `--scratch-dir` looked dead and had its
in-flight files deleted; the victim then recreated the directory silently and reported those
accessions as `no_path`, the biological answer, with exit 0. A vanished GFA or query is now an
error naming the cause. **Directories left by 0.1.x carry no lock and are never swept**, so
clear them out once by hand.

A signal the walker inherited as ignored stays ignored, so `nohup logan-walker ... &` survives
the SSH session dropping. Until 0.2.0 the walker registered a SIGHUP handler regardless, so
bash's hangup killed the run with exit 129 and, because outputs are written at the end, lost
all of it. GraphAligner children inherit the ignore too, which they could not while the signal
was caught.

The summary is always written, even when another output fails: every writer is attempted and
the first error is returned at the end. Each file appears at its final path only once it is
complete, so an interrupted run leaves no truncated output. A rerun with a writer switched off
removes that writer's file from the previous run rather than leaving it beside fresh results,
and a `<ACC>.vcf` this walker did not write is never removed or overwritten.

Reference record ids must be usable as VCF contig names (letters, digits and
`!#$%&+./:;?@^_|~-`, plus `*` and `=` after the first character); a run refuses rather than
writing a VCF no parser will read. A non-ACGT reference base is written as `N` in REF.

`--min-cov` and `--min-vaf` take a finite fraction in [0, 1] and `--min-abund` a finite number
at least 0. `inf` and `NaN` used to be accepted and then written into the summary as bare
tokens no JSON parser reads.

`cov` is the walked length over the reference length, and `pident` the identity of the
alignment. When a haplotype's identity was inherited from the route before an allele was spliced
in, the header says `inherited_pident=` instead.

**`--min-cov` gates on `ref_cov`, not on `cov`.** They answer different questions and an indel
separates them: `cov` says how much sequence came out, and `ref_cov` says how much of the
reference was found, which is what the flag has always been documented as meaning. Both are in
each record in the summary, with `ref_span`, the reference interval the alignment covers.

Until 0.2.0 the gate used `cov`, so an indel in the sample moved the decision in the wrong
direction. On a 300 nt gene at the default 0.85, a sample carrying the whole gene with a 60 nt
in-frame deletion walked 240 bases, scored 0.80 and was reported `no_path`, "the gene is not
here", which biases presence calls against exactly the alleles a deletion study looks for.
An alignment that reached only reference base 240 of 300 but walked a 20 nt insertion scored
0.867 and was published as a primary. Those two now score 1.00 and 0.80 and are decided the
right way round.

`ref_span` is published as the interval rather than only the fraction because a walk that
stops 60 bases from the 3' end is a different finding from one missing 60 bases in the middle,
and the fraction alone cannot tell them apart.

The summary JSON follows the contract in [`summaries.md`](summaries.md), which every stage of
the suite shares. Eight blocks in a fixed order: `schema`, `tool` (this walker, the
`disassembler` that invoked it if any, and GraphAligner's version and path), `inputs` (each
reference's id, length and sha256, and whether the unitigs came from S3 or a cache),
`settings` (every option that changes which accessions are processed or what they produce),
`run` (status, wall clock, per-stage seconds, totals), `discards`, `cohort` and `records`.

`discards` counts everything the run rejected, in two groups: `accessions` (`never_reached`,
which covers duplicates dropped from the list, plus `no_anchors`, `no_path`, `no_links`,
`no_logan_data`, `too_large`, `subgraph_too_large` and `error`) and `alleles` (`no_rejoin`, `not_substitution`, `over_indel_cap`,
`frameshifting`, `window_not_shiftable`, `identical`). Every key is written even at zero, so a
reader can tell "this did not happen" from "this version did not measure it". Before 0.2.0
every one of these numbers existed only as a line on stderr.

`cohort` describes the run as a population rather than a list: how many accessions were
reconstructed, how many are provisional, how many were majority-corrected and at how many
sites in total, how many **distinct** primaries there are and how many of those are
singletons, the alts and bubbles recovered, and cov and pident quantiles.

Each record names its haplotype by `seq_sha256` and `length` rather than carrying the
sequence, which is already in the FASTA. `--summary-seqs` restores the literal `seq` field that
0.1.x wrote unconditionally; it was between 50% and 71% of the file on real runs. A record also
carries `majority_corrections`, the number of sites where the route was moved onto the majority
allele, `minority_rescues`, how many of those the pre-0.2.0 gate would have refused to touch
(WV-02), and `variants`, what the alt search considered and refused. Where a fork is still on
the rarer allele, `minority_left` says so and why: `alt_below_floor`, `not_outright_majority`,
`no_rejoin`, `length_refused`, `window_not_shiftable`, `coverage_gate` or `identical`. The run
totals are in `discards.minority_sites` and `cohort.accessions_with_minority_sites`.

Each accession also records two decisions that earlier versions only printed:

- `subgraph_fraction` is the extracted neighbourhood's share of the graph's unitigs, reported
  for information. `provisional` is true when the neighbourhood holds more than 3 unitigs per
  reference base, the condition the walker warns about: that is how tangled the graph is around
  the gene, which is what makes a walk ambiguous. Through 0.1.1-dev the guard was
  `subgraph_fraction` above 1% instead, which only separates anything when the target is a
  sliver of the sample. On targeted libraries the target is most of the graph, and on the
  20,000-accession SARS-CoV-2 run that guard flagged 9,440 of 9,448 successes while the 8 it
  cleared were no better.
- `prune` is how the pruning pass ended (`null` when it did not run): `outcome` is `applied`,
  `nothing_to_remove`, `would_empty`, `refused_shorter_span`, `refused_no_haplotype`,
  `refused_trimmed` or `refused_error`, with the gene's depth, the floor, the unitigs kept,
  and the reference span of each alignment.
- `bias_uncorrected` is true when the prune ran and its result was thrown away, so the
  first-pass route stands. That route is the one the aligner produced on the unpruned graph,
  where its pull toward the reference is strongest, and `prefer_majority` skips sites whose
  alternative is below the variant floor, so nothing else corrects them either. It is derived
  from `prune.outcome` rather than stored separately, and `cohort.bias_uncorrected` counts it
  across the run, which is the form that is readable at 5,000 accessions.

Two of those outcomes are new in 0.2.0. **`refused_trimmed`** is a realignment that is the
first pass's haplotype with up to three bases off an end and nothing else changed: the
three-base tolerance below exists so that one terminal codon cannot keep a reference-biased
primary, and it is worth paying only when the prune changed something. A shorter copy of the
same sequence is pure loss, and a lost terminal base can be a truncated stop codon.
**`refused_error`** is a second pass that failed outright: GraphAligner aborted, the scratch
filesystem filled, the pruned GFA could not be written. Until 0.2.0 that was not a prune
outcome at all, it was the accession's: an error propagated and a result that had already been
reconstructed became `error`, absent from the FASTA, with `prune: null`.

The run's last lines count the prune outcomes and the provisional results.

VCF positions are 1-based positions in the reference CDS (`CHROM` is the gene). In a
single-reference run every accession's VCF therefore shares one coordinate frame; in a
multi-anchor run each uses its own winning reference (see
[Multi-anchor mode](#multi-anchor-mode)). Frequencies are `ka:f` abundance ratios at the fork, a
proxy for allele frequency, per site and not phased. For a cohort, `snvs.tsv` is the output
that scales: it is sparse, concatenates across shards and loads straight into DuckDB, while a
merged multi-sample VCF becomes a huge, nearly empty matrix. A worked example is in
[`vcf_example.md`](vcf_example.md).

`bubbles.json` is computed from the graph, not from the alignment, so it is written even for
accessions with no path, and most of its bubbles are off the walked path. The VCF, the SNV table
and the alternatives come from forks on the walked path, not from `bubbles.json`.

Three things about how it is counted, all changed in 0.2.0:

- **One record per allele set, not per pair of flanks.** The same two alleles can be reachable
  between several forks and merges, and each pair used to be its own record: `SRR8209082`
  reported 3,632 records for 317 distinct allele sets, one of them sixteen times, and two SNPs
  exactly k apart form a bow-tie that counted as four. A record now carries a `flanks` list,
  whose first entry is the `fork`/`merge` it still reports as its own fields. `n_bubbles`
  therefore counts sites, which is what anyone comparing it across samples was assuming.
- **The floor is the gene's own depth, not an absolute 5x.** A branch must clear the strictest
  of `--prune-frac` of depth, `--min-abund-frac` of depth, and `--min-abund`, and the stage now
  runs after the prune, on the graph the prune left. Before 0.2.0 it ran on the neighbourhood as
  extracted with the flat 5x floor, so on `SRR8209082` every one of those 3,632 bubbles had its
  strongest allele at 66x or less against a 15,000x gene, a median of 11.9x. A default run on
  that accession now reports **none**, which is the honest answer: nothing in its neighbourhood
  is a second population by its own depth. An accession that reconstructed nothing keeps the
  absolute floor, because without a route there is no depth to be a fraction of, so a `no_path`
  record's bubble count is still the old, permissive one.
- **The two alleles must be two unitigs.** An inverted repeat, where one unitig end links to
  both strands of a neighbour, used to produce a bubble whose alleles were one sequence and its
  reverse complement at 0.5 each, and a self-looping unitig could be both an allele and a flank.

A site can be in the VCF and not in `bubbles.json`. The two stages look for different things:
this one finds diamonds, where the alternatives rejoin in a single step, while the variant
stages follow an alternative for up to `--alt-max-detour` unitigs. `SRR8209082`'s real second
population is one of those: its allele rejoins after **31** unitigs (`detour=31` on the alt
record), so it is in the VCF, the SNV table and `alt_haplotypes.fasta`, and it was never in
`bubbles.json`, before this change or after it.

## Offline cache

Every run downloads each accession from S3, which is the biggest cost when iterating on one
cohort. `--unitigs-cache DIR` reads `<DIR>/<ACC>.unitigs.fa.zst` when present and otherwise
downloads it and writes it there (atomically, so a killed run never leaves a partial file).
`--offline` refuses the S3 fallback and turns a cache miss into an error naming the accession.
Cached files are byte-identical to the S3 objects.

A cache hit is opened, not read: the decoder streams from the file handle, so a cached
accession costs no RAM of its own. That makes `--unitigs-cache --offline` the cheapest way to
run, not just the fastest. See [Memory](#memory).

The cache is the only thing that writes unitigs to disk, and it keeps every file it downloads:
the 20,000-accession SARS-CoV-2 run left about 100 GB. Without it, each file is held in memory
until it is parsed and then dropped, so a run needs RAM (see [Memory](#memory)) but almost no
disk.

```bash
# Fill the cache, no analysis (wraps logan-walker --fetch-only)
disassembler fetch -l accessions.txt --unitigs-cache ./unitigs/

# Then run as often as needed with no network
disassembler walker -l accessions.txt -q queries/markers/zika_ns5.fasta -o out/run1 \
    --unitigs-cache ./unitigs/ --offline
```

`fetch` skips accessions already in the cache; delete a file to refetch it. The composites take
the same two flags.

### Accessions with nothing to cache

Some accessions have no object to write: Logan holds no unitigs for them, or they are over
`--max-size-mb`. Both are properties of the accession rather than failures of the run, and
online the walker reports them as `no_logan_data` and `too_large` and carries on. Until 0.2.0
the follow-up `--offline` run could not tell: it saw an ordinary cache miss and called it an
`error`, so the same list with the same flags gave exit 0 online and exit 1, or 3, offline.
Since 0.2.0 any error also makes the run non-zero, which turned that into a shard of 404s
failing outright.

The fetch now leaves a one-line note, `<DIR>/<ACC>.unitigs.absent`, saying which negative it
was (`no_logan_data`, or `too_large <bytes>`), and the offline run reads it back and reports
the same status it would have reported online. Three things follow from that:

- **Deleting a note forces the next online run to ask the bucket again.** A successful fetch
  removes it by itself, because Logan does add data for accessions it had none for.
- **A `too_large` note stops applying once `--max-size-mb` is raised past the size it
  records**, so filling a cache at one cap and walking at a higher one does not keep refusing.
  The accession is then an ordinary cache miss, because the object really is not there.
- **A damaged note is no note at all.** It degrades to the plain cache-miss error rather than
  inventing an outcome.

### When a cache entry is not a Logan object

A cache entry used to be accepted on "not empty". A 200 response that was not a Logan object
at all, such as a proxy's HTML error page, was therefore written into the cache under the
accession's name, where it failed that accession with `parse failed: Unknown frame descriptor`
on every later run, online runs included, until someone found and deleted the file by hand.

The walker now checks the first four bytes for the zstd frame magic, both on the body it is
about to cache and on an entry it is about to read. A body that is not zstd is never cached
and the attempt is retried. An entry that is not zstd is removed and downloaded again; under
`--offline` there is nothing to replace it with, so it is refused with a message naming the
file, and the file is left alone for you to inspect.

### What a cached object is a copy of

Beside each cached object the walker writes `<ACC>.unitigs.fa.zst.provenance`, a few lines of
`key value` recording the object's size, its sha256, and the `ETag`, `Last-Modified` and fetch
time of the S3 object it came from. The digest is in there so the record cannot drift from the
file: if it does not match the bytes beside it, something replaced the object and the record is
ignored with a warning rather than believed.

Without it a warm run knows the digest of the bytes it holds and nothing about where they came
from, because the `ETag` is on the wire only on the run that downloaded them. An entry written
before 0.2.0 has no sidecar and is read exactly as it always was.

### `--revalidate`

Logan's `u/` objects are mutable and unversioned. A cache that was correct when it was filled
can be serving bytes the bucket has since replaced, and nothing in any output says so: two runs
of the same command analyse different data and agree in every field that records what was
asked for.

`--revalidate` makes each cache hit ask. It sends a conditional request against the `ETag` in
the sidecar, so the confirming answer is a 304 with no body: a cache that is still current
costs one round trip per accession and no transfer. It needs `--unitigs-cache`, and it is
refused with `--offline`, which forbids the network it depends on.

What it records, per accession, in `payload.revalidated`:

| value | meaning |
|---|---|
| `current` | the bucket confirmed the cached object |
| `replaced` | the bucket serves different bytes now. The cache has been refreshed and **this run used the new object**, and an earlier run of the same command did not |
| `refetched` | nothing recorded what the cached bytes were a copy of, so they were downloaded again rather than trusted |
| `unavailable` | the check could not be made, and the cached bytes were used anyway |

`replaced` is the one to alert on, and it is printed to stderr as well. `unavailable` is
deliberately not fatal: a cached object is still a usable analysis, and turning a network blip
into a lost cohort would be the worse answer. The field is absent when the flag was off, so
"not asked" and "asked, and fine" stay different answers.

It is off by default because it is a round trip per accession that a run which trusts its
cache should not have to pay. Reach for it when a result has to be defensible, not on every
iteration.

`fetch` also writes `fetch.summary.json` into `--out-dir`, listing every accession with the
size and sha256 of the bytes that landed in the cache, plus the `ETag` and `Last-Modified` of
each object it downloaded. That is what makes a filled cache a frozen input you can cite rather
than a directory of files: see [Input provenance](#input-provenance). The file is written even
when the fetch fails, so a partial cache is still accounted for.

## Options

`logan-walker` needs one of `-a`/`-l` and at least one `-q` (`--fetch-only` needs no `-q`).
Options that are on by default are turned off with a `--no-…` form.

### Input and output

| flag | default | meaning |
|---|---|---|
| `-a, --accession <ACC>` | | one SRA accession. At most 64 characters: the value becomes `<out-dir>/<ACC>.vcf`, and a longer one exceeds the filesystem's name limit. An SRA run accession given in lowercase is read as the accession it spells and the run says so; Logan's bucket is case-sensitive, so `srr8209082` fetched as written is a 404, which would be reported as `no_logan_data` for a run Logan has |
| `-l, --accession-list <FILE>` | | file with one accession per line, plain or gzipped; a repeated accession is processed once, and lowercase SRA accessions are read as above, counted and warned about once for the file. Give `-a` or `-l`, not both |
| `-q, --query <FASTA>` | | reference CDS, plain or gzipped; repeatable, and more than one record switches on multi-anchor mode |
| `-o, --out-dir <DIR>` | `walker_output` | output directory |
| `-f, --out-fasta <FILE>` | `reconstructed_haplotypes.fasta` | primary haplotypes |
| `--scratch-dir <DIR>` | `/dev/shm`, else the system temp directory | where the run's private scratch directory is made (see [Scratch files](#scratch-files)) |
| `--allow-empty` | off | exit 0 when no accession was reconstructed because nothing was found; errors still fail the run, beside successes as well as alone (see [Exit status](#exit-status)) |

### Anchoring, fetching and walking

| flag | default | meaning |
|---|---|---|
| `-k, --k <N>` | `21` | anchor k-mer size, 1 to 32 (the Logan graph's own k stays 31) |
| `--min-shared-kmers <N>` | `3` | shared k-mers a unitig needs to be an anchor, 1 to 5 |
| `--hops <N>` | `5` | hop radius of the subgraph around the anchors |
| `--max-size-mb <N>` | `4096` | refuse accessions above this compressed size. See [Choosing the size cap](#choosing-the-size-cap) |
| `--fetch-workers <N>` | `8` | parallel S3 downloads (bandwidth saturates around 8); at least 1 |
| `--fetch-attempts <N>` | `6` | failed attempts per accession before it is given up. 408, 429, 5xx and transport errors are retried with a jittered exponential backoff (1, 2, 4 ... seconds, capped at 60) and `Retry-After` is honoured; an attempt that transferred real data before the connection dropped resumes and does not count as a failure |
| `--fetch-timeout <SECS>` | `300` | budget for one fetch attempt, connect and body included. An interrupted transfer resumes from where it stopped, so this bounds how long a stalled connection goes unnoticed, not how long an accession may take |
| `-w, --workers <N>` | `8` | compute workers (see [Memory](#memory) before raising); at least 1 |
| `-t, --threads <N>` | `1` | threads per GraphAligner call; above 1 is not reproducible |
| `--graphaligner <PATH>` | `GraphAligner` on `PATH` | the GraphAligner executable |
| `--graphaligner-timeout <SECS>` | `3600` | kill a single GraphAligner call after this many seconds and fail that accession only; `0` waits forever. A stall detector, not a budget: the largest subgraph measured here (34,582 nodes) aligns in 0.47 s |
| `--min-cov <F>` | `0.85` | minimum fraction of the reference a primary must cover, measured as the reference span of the alignment (`ref_cov`) and not as the length of what was walked. See [Outputs](#outputs) |
| `--primary-metric <M>` | `matched` | multi-anchor primary selection: `matched`, `cov`, `pident` or `cov_x_pident` |
| `--primary-tie-eps <F>` | `0.002` | references scoring within this of the best are treated as tied and the deepest route wins; `0` restores pre-0.2.0 panel-order behaviour (see [Multi-anchor mode](#multi-anchor-mode)) |
| `--prune-frac <F>` | `0.01` | drop unitigs below this fraction of the gene's depth and re-align; `0.0` disables. The re-aligned result is kept only if it spans as much of the reference as the first, give or take 3 bases, **and** only if it differs from the first by more than those bases. Unitigs on the first-pass route with no deeper alternative are never dropped: see [The prune pass](#the-prune-pass) |
| `--prune-min-abund <F>` | `5.0` | floor under `--prune-frac`: a unitig at or above this `ka:f` survives the prune whatever the fraction works out to. On a gene shallower than `prune-min-abund / prune-frac` it is the whole rule, so at 180x and the default 1% the prune happens at 5.0 and not at 1.8 |
| `--no-prefer-majority` | | keep the reference-matching allele on the primary at two-allele sites |

### Intra-host variants

| flag | default | meaning |
|---|---|---|
| `--min-abund <F>` | `5.0` | minimum `ka:f` for a branch to count as an allele. **The allele floor only**: until 0.2.0 it was also the prune's lower bound, so raising it moved the floor unitigs are deleted at and could change the published primary. That is `--prune-min-abund` now. When the prune is refused, the variant stages inherit the prune's floor instead, since the sub-floor branches it meant to remove are still in the graph |
| `--min-vaf <F>` | `0.10` | minimum variant allele frequency for an alternative |
| `--min-abund-frac <F>` | `0.0` (off) | raise the floor to this fraction of the gene's depth (median `ka:f` on the path) |
| `--alt-max-detour <N>` | `100` | unitigs an alternative may traverse before rejoining the path |
| `--no-indel-alts` | | substitutions only |
| `--allow-frameshift-alts` | off | also accept indels whose length is not a multiple of 3 |
| `--max-alts <N>` | `32` | most alternative alleles per accession. Credible alleles beyond it are counted and warned about, not dropped silently, and the ones omitted are the 3' sites |
| `--max-alt-indel <N>` | `30` | largest length change treated as an indel rather than a rerouting. For an alternative haplotype this is per allele; for the consensus correction it is also a budget for the **whole** set of corrections on one accession |
| `--no-emit-alts` | | do not write `alt_haplotypes.fasta` |
| `--emit-vcf` | off | write `<out-dir>/<ACC>.vcf` |
| `--emit-snv-table` | off | write `<out-dir>/snvs.tsv` |

### Offline

| flag | meaning |
|---|---|
| `--unitigs-cache <DIR>` | read-through cache of compressed unitig files |
| `--offline` | refuse the S3 fallback; requires `--unitigs-cache` |
| `--revalidate` | confirm each cache hit against the bucket before using it; requires `--unitigs-cache`, refused with `--offline` (see [`--revalidate`](#--revalidate)) |
| `--fetch-only` | fill `--unitigs-cache` and exit before any analysis |

## Memory

A compressed unitig file used to be held by a fetch worker waiting to hand it over
(`--fetch-workers`) and by a compute worker for the whole parse (`-w`). At the defaults that is
up to 8 + 8 = 16 files, which at the then 1 GB `--max-size-mb` cap was **16 GB**. That term
is gone.

Three things removed it, all new in 0.2.0:

- **No compressed object is ever held.** A cache hit is opened, not read, and the zstd decoder
  streams from the handle. A download is written to disk as it arrives: into the cache entry
  when `--unitigs-cache` is set, and into an unlinked file under `--spool-dir` otherwise. The
  only path that still builds a `Vec` is a caller who configured neither, which no CLI run can
  be, because the walker always resolves a spool directory. Measured with `--fetch-only` on
  one accession, peak RSS went from 50.5 MiB to **5.0 MiB** on a 48.5 MB object and from
  350.5 MiB to **4.8 MiB** on a 364 MB one: the old figure tracked the object, the new one does
  not move. In exchange, an in-flight accession costs one 128 KiB copy buffer, so 16 of them
  come to 2 MiB whatever `--max-size-mb` says. This is what makes the cap raisable: before it,
  a 28 GB cap would have promised 437 GB of RAM and meant it.
- **The queue between fetch and compute no longer holds anything.** It was four slots deep,
  which is four more capped buffers in the bound for no throughput: a fetch worker with bytes
  in hand and nowhere to put them is already a held buffer. It is now a rendezvous.
- **A cached accession costs a file handle.** Measured over 20 cached 34 MB accessions at
  `-w 1 --fetch-workers 8`, peak RSS fell from 516.4 MiB to **49.9 MiB**, a 10.4x reduction,
  with the wall clock unchanged and every output byte-identical. At `-w 8 --fetch-workers 8`
  the floor is the eight resident graphs, 366.9 MiB.

The object's SHA-256 is computed as the bytes stream past rather than afterwards, so the
provenance digest costs nothing on a download. It is still a separate pass over a cache hit,
whose bytes the walker never saw arrive; on an 11 GB file the nn1 probe measured that pass at
91 s.

On top of that, each compute worker holds a parsed graph, which the walker budgets at 1.7
times the compressed size of its file.

From 0.2.0 that graph is mostly adjacency. The file is decoded **twice**: the first pass builds
adjacency and finds the anchors as the bases stream past, and the second keeps sequence and
abundance for the extracted neighbourhood alone. `ERR2028820` has 1,277,362 unitigs and a
747-node neighbourhood for PfCRT, so the per-unitig data the old parser built was almost
entirely thrown away. Over the same 20 cached accessions, peak RSS fell from 49.9 MiB to
**24.1 MiB** at `-w 1 --fetch-workers 8`, and from 366.9 MiB to **161.4 MiB** at `-w 8`. Against
0.1.1 that is 516.4 MiB to 24.1, about 21x.

From 0.2.0 the adjacency itself leaves RAM. It is the last part of the walker that scaled with
the *file* rather than the gene (12 bytes per unitig, 15.4 MB for `ERR2028820`), and it exists
only from pass 1 until `extract::expand` has finished walking it. Past 4 MiB it spills to two
files under `--spool-dir`, read back with `pread`; `parse::collect` then keeps the
neighbourhood's rows and drops the spool. The files are unlinked the moment they are created,
so they need no cleanup and survive no crash, and a spool directory that cannot be written is a
warning rather than a failed accession: the run simply uses the memory it used to use.

With that, peak RSS over the same 20 cached accessions is **12.9 MiB** at `-w 1
--fetch-workers 8` and **77.8 MiB** at `-w 8`, against 516.4 MiB at 0.1.1: about **40x**, with
no change to any output.

The second decode is not free, but most of what it used to cost was work whose result was
discarded. `parse::collect` originally called the same header reader pass 1 uses, which parses
`ka:f:` for every record and all three components of every `L:`. Pass 2 has the adjacency
already and skips most records, so on `ERR2028820` that was 1,277,362 abundances and 2,582,206
links parsed per accession for nothing. It now uses `header_light`, which reads the id, takes
the abundance only for a record it is keeping, and never looks at a link. For scale, zstd
decompresses that accession's 146 MB of text in 0.09 s, so the header walk, not the inflate,
was the larger half of pass 2.

Measured over the same 20 cached accessions at `-w 1 --fetch-workers 8`, pass 2 costs about
**0.15 s per accession**, down from 0.21 s, and the two-pass parser is now **1.4% slower** than
the single-pass one it replaced (18.27 s against 18.02 s over the 20), where before the change
it was 10.4% slower. On a 48-accession cohort spanning 18 MB to 695 MB compressed, 9.76 GB in
all, at `-w 8`: pass 2 fell from 94.9 s to 62.3 s (**34%**), pass 1 was unchanged at 234.5 s
against 237.2 s, and wall went from 49.6 s to 44.3 s. Every output byte was identical, and peak
RSS did not move (83.8 MiB against 82.5 MiB), which is the point: this one is time only.

Note what the remaining cost tracks: pass 2 stops once the last wanted unitig has gone
past, so it is governed by **where the neighbourhood sits in the file**, not by how large it is.
A neighbourhood scattered to the end costs a full second decode however few nodes it has.

What pass 2 stopped validating, pass 1 still validates. A malformed link, or an abundance that
is unparseable, negative or not finite, fails in `scan`, over the same bytes from the same
`Compressed` handle, so `collect` can never be the pass that meets one first. The id rules stay
in both readers, because a disagreement between the two decodes has to stay an error rather
than becoming a wrong answer.

The 1.7 figure was calibrated before that split, on the single-pass graph, measured one
accession at a time as peak memory minus the process baseline and the compressed file. It is
no longer the budget (see below); it survives as the fallback for a run with no spool directory
or no neighbourhood cap, and as the reason a two-pass graph sits so far under what a
single-pass one needed:

| accession | compressed | single-pass graph at peak |
|---|---|---|
| `ERR6174245` | 4.8 MB | 1.6x |
| `ERR2028820` | 32.4 MB | 1.25x |
| `SRR25605386` | 278.6 MB | 1.24x |
| `SRR15477033` | 289.1 MB | 1.62x |
| `ERR9831160` | 347.6 MB | 1.43x |

### What the budget actually is

A graph's budget is **`--max-subgraph-nodes` x 150 bytes**, plus the flush buffers the parser
keeps on the way to the spool. At the defaults that is 0.15 GB per compute worker, so `-w 8`
is **1.18 GB, and the same number at a 4 GB cap and at a 28 GB one**. Nothing in the walk
scales with accession size: the compressed object is on disk, the adjacency is on disk, and
what is left in RAM is the gene's neighbourhood.

That budget is both what the startup line prints and what the parser enforces, which was not
true before 0.2.0's last change. The line used to quote `--max-size-mb` x 1.7 and a second line
then corrected it by a factor of about fifty; whichever of the two you believed, one of them
was wrong. They are now one expression, so the promise and the refusal cannot drift apart.

It is capped by the work available: a run of n accessions holds at most min(n, `-w`) graphs, so
a single-accession run prints 0.15 GB rather than 1.18, and `--fetch-only` prints no graphs at
all because it never parses. (Before 0.1.1 the walker counted a fixed 43 MB per graph; before
0.2.0 the compressed-buffer term was 20 GB, because the fetch-to-compute queue held four more,
then 16 GB, then nothing.)

**Without a spool directory, or with `--max-subgraph-nodes 0`, nothing bounds the graph by the
gene** and the budget falls back to `--max-size-mb` x 1.7. That is the only remaining place the
size cap touches memory, and no CLI run reaches it: the walker always resolves a spool
directory, and the node cap defaults to 1,000,000.

The fallback matters in one more case. If adjacency cannot spill, because the spool directory
is unwritable, the parser keeps it in RAM, and it is still held to the budget above: the
accession is refused as `too_large`, next to the warning naming the directory that failed,
rather than quietly growing into fifty times the memory the run announced. Spilling trades RAM
for disk, and that disk is a resource the run states and enforces. See [Disk](#disk).

That 1.7 figure came from measuring real accessions one at a time, and two things make it hold
over a long list rather than drift:

- **The allocator keeps its promise.** glibc raises its own mmap threshold each time a large
  block is freed, so from the second accession on the graph arrays grow inside the heap and
  freed chunks stay resident. The walker pins the threshold at 1 MiB on glibc, which returns
  each graph to the OS when it is dropped. Worth 46 MiB of the measurement above.
- **The parser enforces the budget.** 1.7x is a property of the accessions measured, not of the
  file format: a dense enough file can build a graph many times its compressed size, and
  nothing used to stop it. Both passes measure what the graph holds and refuse the accession as
  `too_large` once it passes the budget, so neither a link-dense file in pass 1 nor a
  neighbourhood that turns out to be most of the graph in pass 2 can exceed the printed bound.
  The arrays grow by doubling and are measured periodically, so the refusal lands within one
  doubling of the budget rather than exactly on it, which is why the bound says "up to".

This is a worst case: it assumes every accession sits at the cap. `ERR2028820` is 33 MB, and the
largest *C. auris* ERG11 accession tested was 343 MB. To contain a run, lower `--max-size-mb` and
the two worker counts together, since they multiply. `-w 16 --fetch-workers 16`, for example,
allows 32 buffers and 16 graphs. An accession whose parsed graph would exceed the budget above
is skipped as `too_large`, and so is one over the cap, wherever its bytes came from (see
[Choosing the size cap](#choosing-the-size-cap)).

## Disk

The spool is a far more predictable quantity than the graph it replaced. It is two files of
`u32` and nothing else: one offset per unitig plus a terminator, and one entry per edge. No
sequences, no abundances, and no slack from an allocator growing by doubling, because it is
written out rather than grown.

| accession | compressed | spool | spool / compressed | edges per unitig |
|---|---|---|---|---|
| `ERR2028820` | 32.4 MB | 15.4 MB | 0.477 | 2.02 |
| `DRR240330` | 11.02 GB | 5.24 GiB | 0.510 | 2.64 |
| `DRR212963` | 10.99 GB | 5.41 GiB | 0.529 | 2.64 |

An 11% spread over a 340x size range, against 1.24x to 1.62x for the memory constant above.
The walker uses **0.55**, the largest measured ratio rounded up rather than an average,
because the term that drifts is edges per unitig and it drifts upward with size: larger
assemblies are branchier, so a constant taken from small accessions under-predicts the large
ones this exists to bound.

From it the walker prints a disk bound at startup alongside the memory ones, with the
consumers named separately rather than summed, because they do not share a shape:

- **The downloaded objects are transient**, when there is no cache for them to live in. One
  per in-flight accession, so `-w 8 --fetch-workers 8` at the 4 GB cap is 64 GB. This term is
  precisely the RAM that streaming moved onto disk, and it is the same size: that is the
  trade, stated rather than hidden.
- **The spool is transient.** It dies with its accession, and the ceiling is one spool per
  compute worker, so `-w 8` at the 4 GB cap is 17.6 GB. Lowering `-w` lowers it.
- **The cache is persistent.** `--unitigs-cache` outlives the run and grows with the accession
  list, not with the worker count, so 50 accessions at the 4 GB cap is 200 GB whatever `-w` is.
  With a cache set, a download lands in it and is counted here instead of as a transient
  object, never in both.

Summing them would suggest that lowering the worker counts helps with the largest of them, and
at any realistic cohort size it does not.

`--spool-dir` defaults to `--unitigs-cache` when set, otherwise the system temp directory, and
deliberately not `--scratch-dir`, which is `/dev/shm` by default: spilling into a tmpfs moves
bytes from one part of RAM to another and measures as a saving while being none. The walker
warns if you point it there anyway.

Three layers keep a full disk from turning into a failed run:

1. **The fetch refuses what will not fit.** An accession's spool is predicted from its real
   size, taken from `Content-Length` or from the cached file, never from `--max-size-mb`. The
   cap is a worst case that assumes every accession sits at it, and deciding on it here would
   refuse runs needing a fiftieth of what it predicts. The prediction is multiplied by the
   number of spools the run can hold at once and compared with the free space under
   `--spool-dir`; over it, the accession is a `too_large` negative, decided before any body is
   transferred when the server sent a `Content-Length`. An accession too small to spill at all
   is never gated, and a filesystem that will not answer disables the check rather than
   refusing everything on it.
2. **A download that runs out of space is a negative.** A write that fails with `ENOSPC` or a
   quota refusal is reported as `too_large`, the same as the size cap. Any other write failure
   is a permanent fetch error rather than a transient one: retrying a read-only directory
   three times only downloads the object three times. A cache entry that cannot be *opened*
   falls through to the spool directory with a warning, and a spool file that cannot be opened
   falls through to RAM, because a storage problem should cost performance rather than the
   accession.
3. **A spool write that runs out of space is a negative too.** The free space seen at the gate
   cannot account for growth that has not happened yet, so a disk can still fill mid-parse.
   `ENOSPC` and a quota refusal are reported as `too_large`, not as a parse error, so one
   oversized accession filling the disk does not make the run tell a scheduler that the walker
   failed. Any other write failure stays an error. A spool that cannot be **created** is
   gentler still: a warning, and that accession keeps its adjacency in RAM.

## Choosing the size cap

`--max-size-mb` defaults to **4096**, raised from 1024 in 0.2.0. It no longer guards memory,
which does not scale with accession size at all any more; it guards disk, time, and relevance.

The number comes from the bucket. A 57,000-accession sample of Logan object sizes, taken with
S3 `ListObjectsV2` (which returns the size with the key, so 57 requests cover 57,000
accessions) spread across 57 prefixes:

| percentile | size | | cap | refuses | worst-case transient disk at the defaults |
|---|---|---|---|---|---|
| p50 | 3.4 MB | | 350 MB | 5.57% | 7.0 GiB |
| p75 | 26.2 MB | | 1024 MB | 2.26% | 20.4 GiB |
| p90 | 168.3 MB | | 2048 MB | 0.90% | 40.8 GiB |
| p95 | 395.1 MB | | **4096 MB** | **0.16%** | **81.6 GiB** |
| p99 | 1,857.5 MB | | 8192 MB | 0.07% | 163.2 GiB |
| p99.9 | 5,796.3 MB | | 28000 MB | 0.00% | 557 GiB |
| max | 19,689.1 MB | | | | |

The distribution is heavily skewed: the median accession is 3.4 MB and the mean is 94 MB.
**4 GB is the knee.** It cuts refusals fourteenfold against 1 GB and recovers the 2.10% of
accessions in the 1 to 4 GB band (mean 1.81 GB). Past it the gains are negligible while the
disk cost keeps doubling.

The disk column is `16 x cap + 8 x 0.55 x cap` and assumes every one of the sixteen in-flight
accessions sits exactly at the cap, which no real cohort does: at p95 the same run uses about
8 GiB. Setting `--unitigs-cache` moves the download term into the persistent cache term rather
than adding to it. The startup disk line prints both the bound and the free space where it
will land, so check that line rather than this table.

### The cap applies wherever the bytes came from

A cached file over `--max-size-mb` is skipped as `too_large`, exactly as a download of the same
object would be. That used to be a memory argument, and it is not one any more: a cache hit
costs a file handle, and since 0.2.0 nothing in the walk scales with accession size.

It stays because the flag names a property of the accession, not of the transport. It is
threaded from the `disassembler` pipelines into every walk, so if the cache ignored it, a warm
run and a cold run of the same command over the same list would include different accessions,
and `summary.json` records the cap but not the cache's history, so nothing in the output could
explain the difference. `disassembler fetch` has its own `--max-size-mb`, which makes filling a
cache at a high cap and walking at a low one a supported thing to do; this is what makes the
walk's cap the one that decides.

The practical consequence: if you populated a cache with `disassembler fetch --max-size-mb
20000` and the walk skips those accessions, raise the walk's cap to match. The skip names the
flag and the size that tripped it.

### Not the structural ceiling

The cap is deliberately **not** set at the walker's structural ceiling. Unitig ids are packed
into 30 bits, which runs out near 29 GB compressed, and the largest object in the sample is
19.7 GB, so a cap there would refuse nothing. But the cap is also a relevance filter: the two
11 GB accessions measured took about 16 minutes each and both returned `no_path`, which is
what a single-gene walk over an environmental metagenome tends to do. Raise it per run when
you want those accessions, with the disk to match.

## Scratch files

Each accession's subgraph GFA, alignment GAF and a temporary copy of the reference are written to
a scratch directory while it is aligned. Every run makes its own, `logan-walker.<pid>.<suffix>`
with mode 0700, under `--scratch-dir`. The default is `/dev/shm`, which is tmpfs, so the files use
RAM; where it does not exist (macOS), the system temp directory is used instead. The files are
removed as each accession finishes, and the directory is removed when the run ends, on any of
the signals that end a run: SIGINT, SIGTERM, SIGHUP, SIGPIPE, SIGQUIT, SIGUSR1, SIGUSR2,
SIGXCPU and SIGXFSZ. A run killed outright (SIGKILL, or the out-of-memory killer) cannot clean
up, so the next run with the same `--scratch-dir` removes any directory whose process no
longer exists. Set `LOGAN_KEEP_SCRATCH=1` to keep the directory for debugging.

Until 0.2.0 the list was SIGINT, SIGTERM and SIGHUP only, so every other ending left the
directory in `/dev/shm` holding RAM until a later run swept it. The ones that matter in
practice are `logan-walker ... | head` (SIGPIPE, which is what makes a pipeline exit 141),
SIGUSR1 and SIGUSR2, which is how SGE and LSF announce a soft wall-time limit, and SIGXCPU
from a CPU ulimit or a cgroup. A signal the walker inherited as ignored stays ignored, so
`nohup` still works.

SIGALRM, SIGVTALRM and SIGPROF are deliberately not in the list. Nothing sends those to
announce a limit; what does send them is a timer or a sampling profiler belonging to the
process itself, and catching them would end a profiled run rather than clean up after one.

**Any of those signals also stops the GraphAligner processes the run started.** They are the only
processes the walker starts, and these scratch files are what they are working on. Until
0.2.0 nothing tracked them, so `kill <walker pid>` or `timeout 3600 logan-walker` left up to
`-w` of them running, reparented to init, aligning against files that had just been deleted.
Ctrl-C in a terminal never showed this, because a terminal signals the whole foreground
process group, and neither do schedulers that kill a process group or a cgroup; a plain `kill`
to the walker alone does. Each aligner now runs in its own process group, which is what lets
one signal reach a wrapper script's children as well, and on Linux it is additionally marked
to die with the walker, so even SIGKILL and the out-of-memory killer take the aligners with
them.

Before 0.1.1 every run shared `/dev/shm/logan-walker`: created by one user, it failed other
users' runs with a bare "permission denied", and an interrupted run left its files there,
holding RAM until reboot.

`/dev/shm` is only 64 MB in a default Docker container; point `--scratch-dir` elsewhere there.

## Input provenance

The Logan objects are unversioned and mutable. They live at a stable URL, `u/{ACC}/`, and the
bucket has no versioning, so a replaced object is served under the same address with nothing in
the response to say it changed. They do get replaced: a 5,000-key sample across five key ranges
on 2026-11-16 found 3.1% of objects rewritten since the 2024-03 bulk load, in waves dated
2025-12-10, 2025-12-14 and 2026-02-08.

Every accession's record in the summary therefore carries a `payload` block naming the bytes it
was built from:

```json
"payload": {"source": "s3", "bytes": 1861058,
            "sha256": "4f2a...", "etag": "d2dd288b238c3c88b46acdcd78860dce",
            "last_modified": "2024-03-31T17:38:36Z", "fetched_at": "2026-11-16T16:47:50Z"}
```

- `sha256` is of the compressed bytes exactly as read, which are byte-identical to the S3
  object and to the cache file. It is the member that works on every path, because a cache hit
  and an `--offline` run have no HTTP response to take headers from. Two runs consumed the
  same input if and only if these match.
- `source` says which path produced the bytes: `s3`, `cache`, or `none` when none were read.
  A run with a cache and no `--offline` reports `cache+s3` at the run level while each
  accession came from exactly one of the two, and this is where that is recorded.
- `etag` is the object's own identifier. Below S3's 5 GB single-PUT limit it is the MD5 of the
  content, which `--max-size-mb` (default 4096) keeps every ingested object inside, so it
  is directly comparable with what the bucket serves today. Above that limit S3 appends
  `-<parts>` and the value is opaque.
- `fetched_at` is absent for a cache hit, whose bytes may predate the run by any amount.
  `run.started` is the only honest clock there.

A `too_large` skip carries a `payload` too, with the size that decided it, taken from
`Content-Length` before any body was transferred, or from the cached file. The skip's `detail`
says which limit was reached, since `too_large` is not only `--max-size-mb`.

To check whether a cohort's inputs have moved since a run, compare each record's `etag` against
a `HEAD` of its object. The bucket also honours conditional requests, so
`If-None-Match: "<etag>"` returns `304` when the object is unchanged and the full body when it
is not.

### When a fetch goes wrong

A transfer that is interrupted **resumes**. The next attempt asks for the rest with a `Range`
header and an `If-Match` on the object's `ETag`, so a dropped connection costs the time it was
up, not the bytes it transferred. The `If-Match` is what makes that safe: the `u/` objects are
mutable, and without it a resume could staple the tail of a replaced object onto the head of
the one the transfer started with. A server that cannot honour the range answers 200 or 412
and the transfer starts again from the beginning.

Before 0.2.0 every attempt reissued a plain `GET` from byte 0 against a fixed 900 s
whole-request timeout with no flag, which made that timeout a hard floor on transfer rate: an
object whose per-stream rate fell below size/900s could not be fetched **at all**, however many
attempts it was given, because every attempt timed out in the same place. At the 4 GB cap that
floor is 4.7 MB/s per stream.

Retries used to be three attempts at a fixed `500 ms x attempt`, which is 1.5 s in total, and
408 and 429 were not retried at all. A 60-second network, proxy or DNS interruption therefore
discarded one accession per fetch worker every 1.5 s: about 320 accessions marked `error` at
the default eight workers. The default is now six failed attempts with a jittered exponential
backoff, which covers about half a minute, and `--fetch-attempts 10` covers about four. The
jitter matters because fetch workers fail together, for the same reason, and would otherwise
all return together.

### The prune pass

The floor is one number for the whole gene, a fraction of the gene's own sequencing depth, and
the pass removes every unitig in the neighbourhood below it and re-aligns. That is right for a
reference-matching error branch sitting beside a deeper true allele, which is what the pass
exists to remove.

It is wrong for a stretch of the gene that is genuinely covered below the floor: amplicon
dropout, a 3' fall-off, any uneven coverage. Removing it cuts the route in two, the
realignment cannot span the reference, the prune is refused, and the first-pass route stands
with the aligner's reference pull still in it. On ARTIC-style amplicon data that is routine,
not a corner case.

So a unitig on the first-pass route is kept, whatever the floor says, when the route had **no
choice** at that point: no sibling at or above the floor leaving the predecessor by the end the
route leaves by, and none entering the successor by the end it enters by. If there is a deeper
sibling the aligner could have gone that way, and the unitig is a droppable alternative. If
there is not, the unitig is the gene, thinly covered. The ends of a route count as having no
sibling, because there is nothing there to choose between.

When the prune is refused anyway, the record says so as `bias_uncorrected`.

The floor is `max(--prune-frac x depth, --prune-min-abund)`. The two settings are separate,
and `--min-abund` is neither of them: it gates which branches count as alleles, and a run that
raises it to be stricter about variants no longer deletes unitigs at the higher number as a
side effect. The lower bound matters more often than it looks, because on anything shallower
than about 500x the default 1% works out below 5.0, so the prune happens at 5.0: `ERR2028820`
at 180x prunes at 5.0 and loses 268 of its 747 neighbourhood unitigs, all of them off the
route, with a byte-identical haplotype.

The depth the floor is a fraction of is the **median `ka:f` over the route's steps**, not
weighted by unitig length. A k-mer-weighted median was tried and measured: on `SRR8209082`,
whose depth is independently known to be about 14,500x, the step median gives 15,000x and the
weighted one gives 91.7x, which collapses the floor to `--min-abund` and stops the prune
removing anything. A tangled neighbourhood's route threads long, low-coverage unitigs that are
not the gene, while the gene itself is cut into short high-coverage pieces by its own
variation, so weighting by sequence weights the noise.

## Reproducibility

- **Input bytes.** See [Input provenance](#input-provenance). Reproducing a run means
  reproducing its inputs first, and the summary's per-accession `sha256` is what makes that
  checkable rather than assumed.
- **GraphAligner threads.** GraphAligner is not deterministic with more than one thread: on
  `SRR8209082`, `-t 1` scored 2691 on every run, while `-t 2` gave 2686 or 2691 and five runs
  produced two different primaries. Keep `-t 1` and get parallelism from `-w`.
- **GraphAligner version.** Different versions can take different routes through a large
  subgraph, and 1.0.19 and later drop the final query base. The output and the golden fixtures
  are pinned to 1.0.17b.
- **The golden gate.** `bench/check_golden.sh` runs the walker under the original-behaviour flags
  and checks the primaries of `ERR2028820` (PfCRT) and `SRR12363134` (*C. auris* ERG11) byte for
  byte. `SRR8209082` (HIV-1 *pol*) is reported but not gated: its 5-hop subgraph is about 10% of
  the whole graph (against 0.033% or less for the gated accessions), and its unpruned route is
  not reproducible across GraphAligner versions. The walker warns whenever a subgraph exceeds 1%
  of the graph. See [`../bench/README.md`](../bench/README.md).

## Exit status

Each accession ends with one status, printed as it finishes and recorded in the summary JSON:

| status | meaning | kind |
|---|---|---|
| `success` | a primary haplotype was reconstructed | |
| `no_anchors` | no unitig shares enough k-mers with the reference: the gene is not there | nothing found |
| `no_path` | GraphAligner ran, wrote its GAF, and either found no alignment or none spanning `--min-cov` of the reference | nothing found |
| `no_links` | nothing was reconstructed **and** the unitig file carried no `L:` adjacency at all, so no walk could span more than one unitig. A fact about the Logan object, not about the sample: Logan publishes such objects (0 links in 15,864 records for `ERR3988483`). An accession whose gene fits inside a single unitig still reconstructs and is `success` | nothing found |
| `no_logan_data` | Logan has no unitigs for the accession: the bucket returns 404, or it serves an object that decompresses to nothing (a 55-byte empty zstd frame, which about 0.4% of accessions have) | nothing found |
| `too_large` | over `--max-size-mb`, so skipped before downloading; or its parsed graph would exceed the per-accession graph budget; or its spooled adjacency will not fit under `--spool-dir` (see [Disk](#disk)) | nothing found |
| `subgraph_too_large` | the gene's neighbourhood is over `--max-subgraph-nodes` | nothing found |
| `error` | anything else: a network or cache failure, a corrupt file, a cached object that is truncated (a link naming a record the file does not hold, which is what a copy cut at a zstd frame boundary looks like), a GAF row that is not a walk of the graph it was aligned against, GraphAligner crashing, one GraphAligner call passing `--graphaligner-timeout`, or a program on `--graphaligner` that exits 0 without writing a GAF at all (1.0.17b creates the file before it aligns, so that is not an alignment result) | error |

One bad accession does not stop the run: the others are still walked, and whatever was
reconstructed is still written. **It does make the run fail, though.** An accession that
errored is a failure of the run, not a thing that was not found, so from 0.2.0 any error at
all puts the exit code above zero:

| exit | when |
|---|---|
| `0` | no accession errored, and either something was reconstructed or `--allow-empty` is set |
| `1` | nothing was reconstructed: either nothing was found and `--allow-empty` is not set, or some found nothing and the rest errored |
| `3` | every accession errored. The run did not work, so `--allow-empty` does not apply |
| `4` | something was reconstructed and at least one accession errored |

`2` is deliberately absent. The `disassembler` CLI in front of the walker is a Python
argparse program, and argparse exits 2 when it rejects a command line, so a 2 always means
"those arguments were not understood" and never names an outcome of the walk. No walker exit
code is assigned to it.

Exit 4 is the new one, and it is separate from 1 and 3 on purpose: the run produced real
output, so a wrapper that wants to keep `reconstructed.fasta` and carry on can tell this case
from a run that produced nothing. The accessions that failed are listed in
`<query stem>.skipped.tsv` with their `error` outcome and detail.

Until 0.2.0 this case exited 0, because one success was enough for the whole run to pass. A
5,000-accession run on 2026-11-17 in which 8 accessions crashed reported success to the
scheduler that launched it, and the only way to see the failures was to read
`records[].status` out of the summary by hand.

Usage and I/O errors also exit 1. Until 0.1.1, `--allow-empty` turned every one of these into 0,
including a run in which GraphAligner failed on every accession, and `no_logan_data` and
`too_large` were reported as `error`.

With `--fetch-only`, success means every accession is now in the cache. `--allow-empty` accepts
accessions with no Logan data or over the size cap; any other failure exits 1, or 3 when every
accession failed.

If the walker's standard output is closed early (`logan-walker ... | head`), the walker stops
on SIGPIPE, as other Unix tools do, and exits 141. Before 0.1.1 it panicked and then hung.
Through 0.1.1 it ended on the spot, mid-write, with nothing run on the way out; from 0.2.0 it
cleans up first (see [Scratch files](#scratch-files)).

## Logan unitig files

Logan's unitigs are BCALM2 output from a k = 31 compacted de Bruijn graph. A header looks like

```
>ERR2028820_2933 ka:f:5.3 L:-:443893:+ L:+:1239391:-
```

- The unitig id is the text after the last `_`. Ids run from 0 to n-1 with no gaps.
- `ka:f` is the unitig's mean k-mer abundance; the walker uses it as depth. Logan's own format
  notes say that in some large accessions this field is spelled `km:f` instead and that the two
  carry the same information, so the walker reads either. A record carrying neither is not read
  as 0x coverage: those records are counted, warned about once per accession, and reported as
  `abundance_missing` in `summary.json` (per record, and as a count of affected accessions in
  the `cohort` block). The field is written only when the count is non-zero. It matters because
  an abundance the walker cannot read leaves every depth-based stage inert: the prune floor
  collapses, no fork has a majority to correct towards, no branch clears `--min-abund`, and the
  VCF comes back empty, none of which is an error.
- `ka:f` **is clipped at 25,000.** Measured on three cached accessions: `ERR2028820` has 664 unitigs at exactly
  25000.0 against 4 at the next value down, `SRR8209082` has 802 against 24, and
  `SRR12363134` has one at 25000.0 whose next-highest neighbour is 14500.0. Where a branch is
  clipped, the fork total is a lower bound and every frequency drawn from it is **overstated**:
  a 3,000x allele beside a true 250,000x one is 1.2%, but against a clipped 25,000 it reports
  10.7% and clears the default `--min-vaf`. The true depth is not recoverable from the file, so
  the walker flags rather than drops: `KaSat` in the VCF FILTER column, a `ka_saturated` column
  in `snvs.tsv`, `ka_sat=1` on the alt header, `"ka_saturated"` on each bubble, and `ka_max`
  with `depth_saturated` per accession in the summary. `bcftools view -f PASS` now drops these
  rows, which is the intended behaviour and a change from 0.1.x, where they passed unmarked.
- `L:s1:t:s2` in unitig `u`'s header is the traversal from `u` on strand `s1` to `t` on strand
  `s2`. Every edge appears in both endpoints' headers, so the graph is closed under reversal.
  Some objects carry no `L:` field at all, in any record. That is not an error and the walker
  parses them, but nothing can be walked across unitigs, so an accession that reconstructs
  nothing from one is reported as `no_links` rather than `no_path` (see
  [Exit status](#exit-status)).
- Consecutive unitigs overlap by k - 1 = 30 bases, which the walker removes when it stitches a
  path into a sequence.
- Sequences are plain uppercase `ACGT`; anything else is an error.

A 404 from S3 means the accession has no Logan unitigs and is not retried; server errors are
retried three times. The objects are unversioned and mutable, so each run records what it read:
see [Input provenance](#input-provenance).

[`contigs_vs_unitigs.md`](contigs_vs_unitigs.md) explains why the walker uses unitigs rather
than Logan's assembled contigs. The default `-k 21` comes from the simulation in
[`../bench/simulated_kmer_sensitivity/`](../bench/simulated_kmer_sensitivity/README.md).
