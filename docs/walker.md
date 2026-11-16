# logan-walker

`logan-walker` reconstructs the coding haplotypes of one gene in each of a list of SRA runs, from
Logan's unitig graphs. It is a Rust binary (`walker/`), packaged on Bioconda as `logan-walker`.
`disassembler walker` passes its arguments straight to it, and `disassembler reconstruct` and
`from-kmindex` run it as a stage.

## How it works

For each accession, given a reference CDS:

1. **Stream** the compressed unitig file,
   `https://logan-pub.s3.amazonaws.com/u/<ACC>/<ACC>.unitigs.fa.zst`, into memory. An accession
   whose `Content-Length` exceeds `--max-size-mb` is refused before anything is downloaded.
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

`--primary-metric` choices:

| metric | ranks by | trade-off |
|---|---|---|
| `matched` (default) | the fraction of the reference's own bases the alignment matched, then identity | never above 1.0; a reference loses for every base of its own the walk misses or mismatches |
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
| `-f` (default `reconstructed_haplotypes.fasta`) | always | one primary per reconstructed accession, `>ACC_primary cov=… pident=… nodes=…` (plus `anchor=` in multi-anchor mode) |
| `<out-dir>/alt_haplotypes.fasta` | unless `--no-emit-alts` | intra-host alternatives, `>ACC_altN …` with their VAF and fork, at most 32 per accession |
| `<out-dir>/bubbles.json` | always | bubbles in the local subgraph |
| `<out-dir>/<query stem>.summary.json` | always | the run's provenance, settings, totals, discards and cohort statistics, then one record per accession (below) |
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

The summary JSON follows the contract in [`summaries.md`](summaries.md), which every stage of
the suite shares. Eight blocks in a fixed order: `schema`, `tool` (this walker, the
`disassembler` that invoked it if any, and GraphAligner's version and path), `inputs` (each
reference's id, length and sha256, and whether the unitigs came from S3 or a cache),
`settings` (every option that changes which accessions are processed or what they produce),
`run` (status, wall clock, per-stage seconds, totals), `discards`, `cohort` and `records`.

`discards` counts everything the run rejected, in two groups: `accessions` (`never_reached`,
which covers duplicates dropped from the list, plus `no_anchors`, `no_path`, `no_logan_data`,
`too_large` and `error`) and `alleles` (`no_rejoin`, `not_substitution`, `over_indel_cap`,
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
  `nothing_to_remove`, `would_empty`, `refused_shorter_span` or `refused_no_haplotype`, with the
  gene's depth, the floor, the unitigs kept, and the reference span of each alignment.

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
| `-a, --accession <ACC>` | | one SRA accession. At most 64 characters: the value becomes `<out-dir>/<ACC>.vcf`, and a longer one exceeds the filesystem's name limit |
| `-l, --accession-list <FILE>` | | file with one accession per line, plain or gzipped; a repeated accession is processed once. Give `-a` or `-l`, not both |
| `-q, --query <FASTA>` | | reference CDS, plain or gzipped; repeatable, and more than one record switches on multi-anchor mode |
| `-o, --out-dir <DIR>` | `walker_output` | output directory |
| `-f, --out-fasta <FILE>` | `reconstructed_haplotypes.fasta` | primary haplotypes |
| `--scratch-dir <DIR>` | `/dev/shm`, else the system temp directory | where the run's private scratch directory is made (see [Scratch files](#scratch-files)) |
| `--allow-empty` | off | exit 0 when no accession was reconstructed because nothing was found; errors still fail (see [Exit status](#exit-status)) |

### Anchoring, fetching and walking

| flag | default | meaning |
|---|---|---|
| `-k, --k <N>` | `21` | anchor k-mer size, 1 to 32 (the Logan graph's own k stays 31) |
| `--min-shared-kmers <N>` | `3` | shared k-mers a unitig needs to be an anchor, 1 to 5 |
| `--hops <N>` | `5` | hop radius of the subgraph around the anchors |
| `--max-size-mb <N>` | `1024` | refuse accessions above this compressed size |
| `--fetch-workers <N>` | `8` | parallel S3 downloads (bandwidth saturates around 8); at least 1 |
| `-w, --workers <N>` | `8` | compute workers (see [Memory](#memory) before raising); at least 1 |
| `-t, --threads <N>` | `1` | threads per GraphAligner call; above 1 is not reproducible |
| `--graphaligner <PATH>` | `GraphAligner` on `PATH` | the GraphAligner executable |
| `--min-cov <F>` | `0.85` | minimum fraction of the reference a primary must cover |
| `--primary-metric <M>` | `matched` | multi-anchor primary selection: `matched`, `cov`, `pident` or `cov_x_pident` |
| `--prune-frac <F>` | `0.01` | drop unitigs below this fraction of the gene's depth and re-align; `0.0` disables. The re-aligned result is kept only if it spans as much of the reference as the first, give or take 3 bases: a strict rule refused a prune over one terminal base |
| `--no-prefer-majority` | | keep the reference-matching allele on the primary at two-allele sites |

### Intra-host variants

| flag | default | meaning |
|---|---|---|
| `--min-abund <F>` | `5.0` | minimum `ka:f` for a branch to count as an allele. When the prune is refused, the variant stages inherit the prune's floor instead, since the sub-floor branches it meant to remove are still in the graph |
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
| `--fetch-only` | fill `--unitigs-cache` and exit before any analysis |

## Memory

A compressed unitig file can be held by a fetch worker waiting to hand it over
(`--fetch-workers`) and by a compute worker for the whole parse (`-w`). At the defaults that is
up to 8 + 8 = 16 files, so **16 GB** at the 1 GB `--max-size-mb` cap.

Two things keep that figure from being reached in practice, both new in 0.2.0:

- **A cached accession costs a file handle, not its own size.** A cache hit opens the file and
  the zstd decoder streams from the handle, so the bytes are never resident. A download with
  `--unitigs-cache` set is written to the cache and the buffer dropped for the handle, which is
  what removes the compute worker's copy, the one held longest. Only a download with no cache
  configured is still held in RAM. Measured over 20 cached 34 MB accessions at
  `-w 1 --fetch-workers 8`, peak RSS fell from 516.4 MiB to **49.9 MiB**, a 10.4x reduction,
  with the wall clock unchanged and every output byte-identical. At `-w 8 --fetch-workers 8`
  the floor is the eight resident graphs, 366.9 MiB.
- **The queue between fetch and compute no longer holds anything.** It was four slots deep,
  which is four more capped buffers in the bound for no throughput: a fetch worker with bytes
  in hand and nowhere to put them is already a held buffer. It is now a rendezvous.

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

The 1.7 budget was calibrated before that split, on the single-pass graph, measured one
accession at a time as peak memory minus the process baseline and the compressed file. It is
kept as the budget because it is the figure the printed bound and the parser both work from,
and because it is now a ceiling the two-pass graph sits well under rather than a prediction of
what it will use:

| accession | compressed | single-pass graph at peak |
|---|---|---|
| `ERR6174245` | 4.8 MB | 1.6x |
| `ERR2028820` | 32.4 MB | 1.25x |
| `SRR25605386` | 278.6 MB | 1.24x |
| `SRR15477033` | 289.1 MB | 1.62x |
| `ERR9831160` | 347.6 MB | 1.43x |

The walker counts each graph at 1.7 times the cap, so the default worst case is 16 GB plus
8 x 1.7 GB, **29.6 GB**, and it prints the bound for the settings and the accession list it was
actually given at startup. (Before 0.1.1 it counted a fixed 43 MB per graph; before 0.2.0 the
buffer term was 20 GB, because the fetch-to-compute queue held four more.)

The printed bound is capped by the work available: a run of n accessions can hold at most n
buffers and at most min(n, `-w`) graphs, so a single-accession run prints 2.7 GB rather than
29.6. `--fetch-only` prints one buffer per fetch worker and no graphs, because it never parses.

That bound assumes the graph is in RAM and scales with `--max-size-mb`. With adjacency spooled
and `--max-subgraph-nodes` set, neither holds, and the walker prints a second line with the
ceiling that actually applies: about 150 bytes per retained node, so `-w 8` with the default
1,000,000-node cap is **1.12 GB regardless of accession size**. That figure does not mention
the file size at all, which is the point. It is also what `--max-size-mb` can eventually be
relaxed against, since the compressed cap exists to bound memory and memory is now bounded by
the gene instead.

Spilling trades RAM for disk: up to roughly the compressed size again per in-flight accession,
so `-w 8 --fetch-workers 8` at the 1 GB cap can want several GB of scratch under
`--spool-dir`. It defaults to `--unitigs-cache` when set, otherwise the system temp directory,
and deliberately not `--scratch-dir`, which is `/dev/shm` by default: spilling into a tmpfs
moves bytes from one part of RAM to another and measures as a saving while being none. The
walker warns if you point it there anyway.

That 1.7 figure came from measuring real accessions one at a time, and two things make it hold
over a long list rather than drift:

- **The allocator keeps its promise.** glibc raises its own mmap threshold each time a large
  block is freed, so from the second accession on the graph arrays grow inside the heap and
  freed chunks stay resident. The walker pins the threshold at 1 MiB on glibc, which returns
  each graph to the OS when it is dropped. Worth 46 MiB of the measurement above.
- **The parser enforces the budget.** 1.7x is a property of the accessions measured, not of the
  file format: a dense enough file can build a graph many times its compressed size, and
  nothing used to stop it. Both passes measure what the graph holds and refuse the accession as
  `too_large` once it passes `--max-size-mb` x 1.7, so neither a link-dense file in pass 1 nor a
  neighbourhood that turns out to be most of the graph in pass 2 can exceed the printed bound.
  The arrays grow by doubling and are measured periodically, so the refusal lands within one
  doubling of the budget rather than exactly on it, which is why the bound says "up to".

This is a worst case: it assumes every accession sits at the cap. `ERR2028820` is 33 MB, and the
largest *C. auris* ERG11 accession tested was 343 MB. To contain a run, lower `--max-size-mb` and
the two worker counts together, since they multiply. `-w 16 --fetch-workers 16`, for example,
allows 32 buffers and 16 graphs. The cap applies to files read from `--unitigs-cache` too: a
cached file over it is skipped as `too_large`, like a download, and so is one whose parsed graph
would exceed the budget above.

## Scratch files

Each accession's subgraph GFA, alignment GAF and a temporary copy of the reference are written to
a scratch directory while it is aligned. Every run makes its own, `logan-walker.<pid>.<suffix>`
with mode 0700, under `--scratch-dir`. The default is `/dev/shm`, which is tmpfs, so the files use
RAM; where it does not exist (macOS), the system temp directory is used instead. The files are
removed as each accession finishes, and the directory is removed when the run ends, including
on SIGINT, SIGTERM or SIGHUP. A run killed outright (SIGKILL, or the out-of-memory killer)
cannot clean up, so the next run with the same `--scratch-dir` removes any directory whose
process no longer exists. Set `LOGAN_KEEP_SCRATCH=1` to keep the directory for debugging.

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
  content, which `--max-size-mb` (default 1024) keeps every ingested object well inside, so it
  is directly comparable with what the bucket serves today. Above that limit S3 appends
  `-<parts>` and the value is opaque.
- `fetched_at` is absent for a cache hit, whose bytes may predate the run by any amount.
  `run.started` is the only honest clock there.

A `too_large` skip carries a `payload` too, with the size that tripped `--max-size-mb`, taken
from `Content-Length` before any body was transferred.

To check whether a cohort's inputs have moved since a run, compare each record's `etag` against
a `HEAD` of its object. The bucket also honours conditional requests, so
`If-None-Match: "<etag>"` returns `304` when the object is unchanged and the full body when it
is not.

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
| `no_path` | anchors were found, but no alignment through the subgraph covers `--min-cov` of the reference | nothing found |
| `no_logan_data` | Logan has no unitigs for the accession (HTTP 404) | nothing found |
| `too_large` | over `--max-size-mb`, so skipped before downloading, or its parsed graph would exceed `--max-size-mb` x 1.7 | nothing found |
| `subgraph_too_large` | the gene's neighbourhood is over `--max-subgraph-nodes` | nothing found |
| `error` | anything else: a network or cache failure, a corrupt file, GraphAligner crashing | error |

One bad accession does not stop the run, and one success is enough for it to exit 0. When none
succeeds:

| exit | when |
|---|---|
| `0` | every accession found nothing, and `--allow-empty` is set |
| `1` | every accession found nothing, without `--allow-empty`; or some found nothing and the rest errored |
| `3` | every accession errored. The run did not work, so `--allow-empty` does not apply |

Usage and I/O errors also exit 1. Until 0.1.1, `--allow-empty` turned every one of these into 0,
including a run in which GraphAligner failed on every accession, and `no_logan_data` and
`too_large` were reported as `error`.

With `--fetch-only`, success means every accession is now in the cache. `--allow-empty` accepts
accessions with no Logan data or over the size cap; any other failure exits 1, or 3 when every
accession failed.

If the walker's standard output is closed early (`logan-walker ... | head`), the walker stops
on SIGPIPE, as other Unix tools do. Before 0.1.1 it panicked and then hung.

## Logan unitig files

Logan's unitigs are BCALM2 output from a k = 31 compacted de Bruijn graph. A header looks like

```
>ERR2028820_2933 ka:f:5.3 L:-:443893:+ L:+:1239391:-
```

- The unitig id is the text after the last `_`. Ids run from 0 to n-1 with no gaps.
- `ka:f` is the unitig's mean k-mer abundance; the walker uses it as depth. **It is clipped at
  25,000.** Measured on three cached accessions: `ERR2028820` has 664 unitigs at exactly
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
- Consecutive unitigs overlap by k - 1 = 30 bases, which the walker removes when it stitches a
  path into a sequence.
- Sequences are plain uppercase `ACGT`; anything else is an error.

A 404 from S3 means the accession has no Logan unitigs and is not retried; server errors are
retried three times. The objects are unversioned and mutable, so each run records what it read:
see [Input provenance](#input-provenance).

[`contigs_vs_unitigs.md`](contigs_vs_unitigs.md) explains why the walker uses unitigs rather
than Logan's assembled contigs. The default `-k 21` comes from the simulation in
[`../bench/simulated_kmer_sensitivity/`](../bench/simulated_kmer_sensitivity/README.md).
