# logan-walker

`logan-walker` reconstructs the coding sequence of one gene in each of a list of SRA runs, from
Logan's unitig graphs. It is a Rust binary (`walker/`), packaged on Bioconda as `logan-walker`.
`disassembler walk` (`disassembler walker` through 0.1.x, which still works) passes its
arguments straight to it, and `disassembler reconstruct` and
`from-kmindex` run it as a stage.

This page is the guide: what the walker does, how to run it, and how it behaves. Every flag,
every output file, the exit statuses and the resource model are in
[`walker_reference.md`](walker_reference.md).

- [How it works](#how-it-works)
- [Running it](#running-it)
- [What the primary is](#what-the-primary-is)
- [The prune pass](#the-prune-pass)
- [Multi-anchor mode](#multi-anchor-mode)
- [Offline cache](#offline-cache)
- [Choosing the size cap](#choosing-the-size-cap)
- [Input provenance](#input-provenance)
- [Reproducibility](#reproducibility)

## How it works

For each accession, given a reference CDS:

1. **Stream** the compressed unitig file,
   `https://logan-pub.s3.amazonaws.com/u/<ACC>/<ACC>.unitigs.fa.zst`, onto disk as it arrives:
   into `--unitigs-cache` when one is set, otherwise into an unlinked file under
   `--spool-dir`. It is never held in memory. An accession whose `Content-Length` exceeds
   `--max-size-mb` is refused before anything is downloaded.
2. **Parse** it in two streaming passes. The first builds the adjacency and finds the anchors
   as the bases go by, spilling the adjacency to disk once it outgrows 4 MiB; the second keeps
   sequence and abundance for the extracted neighbourhood alone. The decompressed text is
   never held, and what stays in RAM is the gene's neighbourhood rather than the accession.
3. **Anchor** the reference by canonical k-mer containment (`-k 21`): a unitig sharing at least
   `--min-shared-kmers` k-mers with it is an anchor.
4. **Extract** the `--hops` (default 5) neighbourhood around the anchors as a local GFA.
5. **Align** the reference through that subgraph with `GraphAligner -x dbg`.
6. **Stitch** the walked path into the primary haplotype. By default the walker then prunes
   low-depth branches and re-aligns (`--prune-frac`), puts the majority allele on the primary
   (`--prefer-majority`, at two-allele sites), and writes the remaining intra-host alleles as
   alternative haplotypes (`--emit-alts`). Why these exist is in
   [`notes/intra_accession_variance.md`](notes/intra_accession_variance.md).
7. **Write** the primaries, the alternatives, `bubbles.json`, a summary JSON, a `skipped.tsv`
   naming every accession that produced nothing, and optionally per-accession VCFs and a
   cohort SNV table. A journal of the run's progress is written as it goes, so a run that is
   killed still leaves what it had reconstructed.

Accessions are fetched and processed in parallel (`--fetch-workers` download, `-w` compute), and
outputs are written in input order whatever order the work finishes in.

## Running it

```bash
logan-walker -a SRR8209082 -q queries/query_hiv1_pol.fasta -o out -f out/hiv1_pol.fasta
disassembler walk -l accessions.txt -q queries/markers/zika_ns5.fasta -o out/zika --emit-snv-table
```

GraphAligner 1.0.17b must be on `PATH`, or named with `--graphaligner <PATH>` (see
[`notes/graphaligner_end_clipping.md`](notes/graphaligner_end_clipping.md) for why that version).
The walker checks it once at startup, before downloading anything, logs its version, and warns
if it is not 1.0.17b.
`disassembler` finds the walker through `--walker-bin`, then `DISASSEMBLER_WALKER_BIN`, then
`PATH`, then `walker/target/release/logan-walker` in a source checkout.

To reproduce the original one-haplotype-per-accession behaviour (the golden fixtures' setting),
use `--prune-frac 0.0 --no-prefer-majority --no-emit-alts`.

### The reference the aligner sees

The query file is read once, at startup, and GraphAligner is handed the record the walker
parsed rather than a copy of the file: `>{record id}` followed by one line of uppercased
sequence, where the record id is the header's first token. Whitespace inside the sequence is
gone and case is folded, so the aligner and the anchoring stage see the same bases by
construction, and a header that GraphAligner would mis-parse never reaches it.

That also makes a query readable once enough. `-q <(zcat ref.fa.gz)`, `-q /dev/stdin` and
`cat panel.fa | logan-walker -q /dev/stdin` all work, for any number of records.

### The alignment is checked against the graph

GraphAligner is handed a GFA the walker wrote, and the row that comes back is checked against
that same graph before anything is stitched from it: every segment must be one the subgraph
holds, every consecutive pair must be joined by an edge in the orientations the row gives, the
steps must stitch to the `path_len` the row declares, and the alignment window must lie inside
it. A row that fails is an `error` for that accession, not a primary.

This matters because stitching cannot notice. The walk is concatenation with the k - 1 overlap
removed, so a row naming two branches of the same fork in series produces a clean sequence that
exists in no organism, and it is published carrying the identity the row itself claimed. An
inconsistent row built for the test suite yields a primary with 89 differences from the truth
at `pident=100.0`. The causes are a different GraphAligner build, a wrapper that renames
segments, or a truncated GAF. On the second pass the same check applies to the pruned graph,
where a failure refuses the prune and keeps the first pass rather than losing the accession.

The path column is read strictly for the same reason: `>` or `<` then a unitig id, or GAF's
`*` for "no alignment". Anything else fails the accession rather than being skipped, because a
row the walker cannot read must not become a shorter walk.

## What the primary is

`<ACC>_primary` is the **majority-consensus path** through the accession's subgraph, not a
phased haplotype, and the difference matters whenever a sample holds more than one strain.

Requiring a connected path is what keeps out joins that no k-mer supports, which is the defect
that makes concatenated BLAST hits unusable (87% chimeric artifacts, measured on this
project's own data). It does not phase anything. In a k = 31 graph, two allele sites more than
30 bases apart are joined through unitigs both strains carry, so every combination of their
alleles is a connected path,
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

## The prune pass

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

The depth the floor is a fraction of is the **median `ka:f` over the route's steps**,
deliberately not weighted by unitig length. Weighting it by sequence weights the noise: a
tangled neighbourhood's route threads long, low-coverage unitigs that are not the gene, while
the gene is cut into short high-coverage pieces by its own variation. Measured on
`SRR8209082`, whose depth is independently known to be about 14,500x, the step median gives
15,000x and a k-mer-weighted median gives 91.7x, which collapses the floor to `--min-abund`
and stops the prune removing anything at all.

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
divergent panel also splits and only the user knows which was meant.
`queries/master_planetary_targets.fasta` is a catalogue, not a panel: its 142 records fall into
31 groups, and it should be used one gene at a time.

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

## Offline cache

Every run downloads each accession from S3, which is the biggest cost when iterating on one
cohort. `--unitigs-cache DIR` reads `<DIR>/<ACC>.unitigs.fa.zst` when present and otherwise
downloads it and writes it there (atomically, so a killed run never leaves a partial file).
`--offline` refuses the S3 fallback and turns a cache miss into an error naming the accession.
Cached files are byte-identical to the S3 objects.

A cache hit is opened, not read: the decoder streams from the file handle, so a cached
accession costs no RAM of its own. That makes `--unitigs-cache --offline` the cheapest way to
run, not just the fastest. See [Memory](walker_reference.md#memory).

What the cache changes is persistence, not whether unitigs reach the disk. Every run spools a
download as it arrives; without a cache that file is unlinked once it has been parsed, and with
one it stays. So the cache trades disk for network, and it keeps everything it downloads: the
20,000-accession SARS-CoV-2 run left about 100 GB behind.

```bash
# Fill the cache, no analysis (wraps logan-walker --fetch-only)
disassembler fetch -l accessions.txt --unitigs-cache ./unitigs/

# Then run as often as needed with no network
disassembler walk -l accessions.txt -q queries/markers/zika_ns5.fasta -o out/run1 \
    --unitigs-cache ./unitigs/ --offline
```

`fetch` skips accessions already in the cache; delete a file to refetch it. The composites take
the same two flags.

### Accessions with nothing to cache

Some accessions have no object to write: Logan holds no unitigs for them, or they are over
`--max-size-mb`. Both are properties of the accession rather than failures of the run, and
online the walker reports them as `no_logan_data` and `too_large` and carries on. An offline
run has no way to see that by itself: all it has is a cache miss, which is an error, so the
same list with the same flags would succeed online and fail offline.

So the fetch leaves a one-line note, `<DIR>/<ACC>.unitigs.absent`, saying which negative it
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

A 200 response is not proof of a Logan object. A proxy's HTML error page is a body like any
other, and written into the cache under an accession's name it fails that accession with
`parse failed: Unknown frame descriptor` on every later run, online runs included, until
someone finds and deletes the file by hand.

So the walker checks the first four bytes for the zstd frame magic, both on the body it is
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
object would be. That is not a memory argument: a cache hit costs a file handle, and nothing in
the walk scales with accession size.

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

## Input provenance

The Logan objects are unversioned and mutable. They live at a stable URL, `u/{ACC}/`, and the
bucket has no versioning, so a replaced object is served under the same address with nothing in
the response to say it changed. They do get replaced: a 5,000-key sample across five key ranges
on 2026-11-16 found 3.1% of objects rewritten since the 2024-03 bulk load, in waves dated
2025-12-10, 2025-12-14 and 2026-02-08.

So every accession's record in the summary carries a `payload` block naming the bytes it was
built from: where they came from, how many there were, their SHA-256, and the object's `ETag`,
`Last-Modified` and fetch clock where an HTTP response supplied them.
[`summaries.md`](summaries.md) defines the fields; two of them decide what you can do with
them:

- **`sha256`** is of the compressed bytes exactly as read, which are byte-identical to the S3
  object and to the cache file. It is the member that works on every path, because a cache hit
  and an `--offline` run have no HTTP response to take headers from. Two runs consumed the same
  input if and only if these match.
- **`etag`** is the object's own identifier, and below S3's 5 GB single-PUT limit it is the MD5
  of the content, which `--max-size-mb` keeps every ingested object inside. So it is directly
  comparable with what the bucket serves today.

To check whether a cohort's inputs have moved since a run, compare each record's `etag` against
a `HEAD` of its object. The bucket honours conditional requests, so `If-None-Match: "<etag>"`
returns `304` when the object is unchanged and the full body when it is not. `--revalidate`
does this for you on every cache hit.

A `too_large` skip carries a `payload` too, with the size that decided it, taken from
`Content-Length` before any body was transferred or from the cached file, and its `detail` says
which limit was reached.

### When a fetch goes wrong

A transfer that is interrupted **resumes**. The next attempt asks for the rest with a `Range`
header and an `If-Match` on the object's `ETag`, so a dropped connection costs the time it was
up, not the bytes it transferred. The `If-Match` is what makes that safe: the `u/` objects are
mutable, and without it a resume could staple the tail of a replaced object onto the head of
the one the transfer started with. A server that cannot honour the range answers 200 or 412 and
the transfer starts again from the beginning.

`--fetch-attempts` (default 6) bounds the attempts and `--fetch-timeout` (default 300 s) bounds
one attempt rather than the accession, since a resumed transfer keeps what it had. Backoff is
exponential and jittered, which covers about half a minute at the default and about four
minutes at `--fetch-attempts 10`. The jitter matters because fetch workers fail together, for
the same reason, and would otherwise all return together. 408, 429, 5xx and transport errors
are retried and `Retry-After` is honoured; a 404 is not, because it means Logan has no unitigs
for that accession.

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
- **The golden gate.** `bench/check_golden.sh` runs two sets. The default-flag set gates every
  output of all three fixture accessions byte for byte. The old-behaviour set runs the walker
  under `--prune-frac 0 --no-prefer-majority --no-emit-alts` and gates only the bare primary,
  and there `SRR8209082` (HIV-1 *pol*) is reported rather than gated: its 5-hop neighbourhood
  is 42% of the whole graph, against 0.06% and 0.02% for the other two, so its **unpruned**
  route is not reproducible across GraphAligner versions. The default pruning is what settles
  it. See [`../bench/README.md`](../bench/README.md).
