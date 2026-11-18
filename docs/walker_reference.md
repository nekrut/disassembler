# logan-walker reference

Every option, every file the walker writes, what each exit status means, and what the run
costs in memory and disk. The guide, which explains what the walker does and how it behaves,
is [`walker.md`](walker.md).

- [Options](#options)
- [Outputs](#outputs)
- [Exit status](#exit-status)
- [Memory](#memory)
- [Disk](#disk)
- [Scratch files](#scratch-files)
- [Logan unitig files](#logan-unitig-files)

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
| `--spool-dir <DIR>` | `--unitigs-cache` when set, else the system temp directory | where a download, and the adjacency past 4 MiB, go. Deliberately not `--scratch-dir`, which is tmpfs by default, where spilling moves bytes from one part of RAM to another. An accession whose spool will not fit the free space here is skipped as `too_large`, from its real size rather than from `--max-size-mb` (see [Disk](#disk)) |
| `--summary-seqs` | off | also write each haplotype's literal sequence into the summary, as 0.1.x did. The summary otherwise carries `seq_sha256` and `length`, and the sequence itself is in the FASTA |
| `--allow-empty` | off | exit 0 when no accession was reconstructed because nothing was found; errors still fail the run, beside successes as well as alone (see [Exit status](#exit-status)) |

### Anchoring, fetching and walking

| flag | default | meaning |
|---|---|---|
| `-k, --k <N>` | `21` | anchor k-mer size, 1 to 32 (the Logan graph's own k stays 31) |
| `--min-shared-kmers <N>` | `3` | shared k-mers a unitig needs to be an anchor, 1 to 5 |
| `--hops <N>` | `5` | hop radius of the subgraph around the anchors |
| `--max-size-mb <N>` | `4096` | refuse accessions above this compressed size. See [Choosing the size cap](walker.md#choosing-the-size-cap) |
| `--max-subgraph-nodes <N>` | `1000000` (`0` = off) | refuse an accession whose anchor neighbourhood exceeds this many unitigs, as `subgraph_too_large`. Everything after extraction is proportional to the neighbourhood, so this, and not `--max-size-mb`, is what bounds memory once the gene is found (see [Memory](#memory)) |
| `--fetch-workers <N>` | `8` | parallel S3 downloads (bandwidth saturates around 8); at least 1 |
| `--fetch-attempts <N>` | `6` | failed attempts per accession before it is given up. 408, 429, 5xx and transport errors are retried with a jittered exponential backoff (1, 2, 4 ... seconds, capped at 60) and `Retry-After` is honoured; an attempt that transferred real data before the connection dropped resumes and does not count as a failure |
| `--fetch-timeout <SECS>` | `300` | budget for one fetch attempt, connect and body included. An interrupted transfer resumes from where it stopped, so this bounds how long a stalled connection goes unnoticed, not how long an accession may take |
| `-w, --workers <N>` | `8` | compute workers (see [Memory](#memory) before raising); at least 1 |
| `-t, --threads <N>` | `1` | threads per GraphAligner call; above 1 is not reproducible |
| `--graphaligner <PATH>` | `GraphAligner` on `PATH` | the GraphAligner executable |
| `--graphaligner-timeout <SECS>` | `3600` | kill a single GraphAligner call after this many seconds and fail that accession only; `0` waits forever. A stall detector, not a budget: the largest subgraph measured here (34,582 nodes) aligns in 0.47 s |
| `--min-cov <F>` | `0.85` | minimum fraction of the reference a primary must cover, measured as the reference span of the alignment (`ref_cov`) and not as the length of what was walked. See [Outputs](#outputs) |
| `--primary-metric <M>` | `matched` | multi-anchor primary selection: `matched`, `cov`, `pident` or `cov_x_pident` |
| `--primary-tie-eps <F>` | `0.002` | references scoring within this of the best are treated as tied and the deepest route wins; `0` decides a tie by panel order instead (see [Multi-anchor mode](walker.md#multi-anchor-mode)) |
| `--prune-frac <F>` | `0.01` | drop unitigs below this fraction of the gene's depth and re-align; `0.0` disables. The re-aligned result is kept only if it spans as much of the reference as the first, give or take 3 bases, **and** only if it differs from the first by more than those bases. Unitigs on the first-pass route with no deeper alternative are never dropped: see [The prune pass](walker.md#the-prune-pass) |
| `--prune-min-abund <F>` | `5.0` | floor under `--prune-frac`: a unitig at or above this `ka:f` survives the prune whatever the fraction works out to. On a gene shallower than `prune-min-abund / prune-frac` it is the whole rule, so at 180x and the default 1% the prune happens at 5.0 and not at 1.8 |
| `--no-prefer-majority` | | keep the reference-matching allele on the primary at two-allele sites |

### Intra-host variants

| flag | default | meaning |
|---|---|---|
| `--min-abund <F>` | `5.0` | minimum `ka:f` for a branch to count as an allele. **The allele floor only**; the prune's lower bound is `--prune-min-abund`. When the prune is refused, the variant stages inherit the prune's floor instead, since the sub-floor branches it meant to remove are still in the graph |
| `--min-vaf <F>` | `0.10` | minimum variant allele frequency for an alternative |
| `--min-abund-frac <F>` | `0.0` (off) | raise the floor to this fraction of the gene's depth (median `ka:f` on the path) |
| `--rescue-minority` | off | also move the route off a branch the aligner walked that is itself a minority of its fork, onto an outright majority. **Off by default**: on a 5,000-accession slice the correction loop did not converge with it on, making 97,068 corrections against 470 and 7.8 per net substitution. The sites are reported either way, in `minority_left` |
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
| `--revalidate` | confirm each cache hit against the bucket before using it; requires `--unitigs-cache`, refused with `--offline` (see [`--revalidate`](walker.md#--revalidate)) |
| `--fetch-only` | fill `--unitigs-cache` and exit before any analysis |

## Outputs

| file | when | contents |
|---|---|---|
| `-f` (default `reconstructed_haplotypes.fasta`) | always | one primary per reconstructed accession, `>ACC_primary cov=.. pident=.. nodes=.. ka=.. vaf=.. min_vaf=..` (plus ` mixed=1` when the route crosses two or more two-allele sites, and `anchor=` in multi-anchor mode) |
| `<out-dir>/alt_haplotypes.fasta` | unless `--no-emit-alts` | intra-host alternatives, `>ACC_altN ..` with their VAF and fork, at most `--max-alts` per accession |
| `<out-dir>/bubbles.json` | always | two-allele sites in the local subgraph, one record per allele set with the flanks it was found between |
| `<out-dir>/<query stem>.summary.json` | always | the run's provenance, settings, totals, discards and cohort statistics, then one record per accession |
| `<out-dir>/<query stem>.skipped.tsv` | always | one row per accession that produced no haplotype: `accession`, `outcome`, `bytes`, `detail`. Written even when empty, so it can be read without first testing whether it exists |
| `<out-dir>/<ACC>.vcf` | `--emit-vcf` | intra-host SNVs of one accession, written only when it has at least one (a VCF an earlier run left for an accession with none is removed) |
| `<out-dir>/snvs.tsv` | `--emit-snv-table` | one row per (accession, site, allele) across the run. Rows sharing a `block` come from one alternative haplotype, so they do co-occur |
| `<out-dir>/<query stem>.partial.fasta` | while a run is in flight | the haplotypes so far, in completion order, byte-identical to the records the final FASTA will hold |
| `<out-dir>/<query stem>.journal.jsonl` | while a run is in flight | one line per finished accession: `index`, `acc`, `status`, `haplotypes`, `length`, `seq_sha256`, `error` |

Every file is written at its final path only once it is complete, so an interrupted run leaves
no truncated output. A rerun with a writer switched off removes that writer's file from the
previous run rather than leaving it beside fresh results, and a `<ACC>.vcf` this walker did not
write is never removed or overwritten. The summary is written even when another output fails:
every writer is attempted and the first error is returned at the end.

The output locations are checked before the cohort starts, so an `--out-dir` that cannot be
made or an `-f` whose parent does not exist fails in a fraction of a second rather than after
every accession has been fetched and aligned. The probe leaves the filesystem as it found it.

### What a killed run leaves

Every file above except the last two is written after the last accession finishes, so a
walltime SIGTERM, an OOM kill or a dropped connection part way through a long shard would
otherwise throw away every reconstruction.

The two journal files are written and flushed as each accession lands, so a run that dies
leaves the sequences it had already made and a per-accession record of how far it got. They are
a salvage, not an output: the alternatives, the VCFs, `bubbles.json` and `summary.json` come
only from a run that reaches the end, and the journal's `acc` column is the list of accessions
a rerun can leave out. On a clean finish both are removed, so their presence means a run is in
flight or one died here.

### The primary's header

`cov` is the walked length over the reference length and `pident` the identity of the
alignment. When a haplotype's identity was inherited from the route before an allele was
spliced in, the header says `inherited_pident=` instead, so a consumer reading `pident=` never
gets a number measured on a different sequence.

`vaf=` on a primary is the constant 1.0 and is kept for one release. `min_vaf=` is the measured
counterpart: the lowest frequency of any allele this sequence carries. ` mixed=1` is added from
two two-allele sites up, where the published sequence is a combination nothing in the data
links. See [What the primary is](walker.md#what-the-primary-is).

**`--min-cov` gates on `ref_cov`, not on `cov`.** They answer different questions and an indel
separates them: `cov` says how much sequence came out, `ref_cov` says how much of the reference
was found, which is what the flag has always meant. On a 300 nt gene at the default 0.85, a
sample carrying the whole gene with a 60 nt in-frame deletion walks 240 bases; gated on `cov`
it scores 0.80 and is reported `no_path`, "the gene is not here", which biases presence calls
against exactly the alleles a deletion study looks for. Both numbers are in each summary
record, with `ref_span`, the reference interval the alignment covers, because a walk that stops
60 bases from the 3' end is a different finding from one missing 60 bases in the middle.

### summary.json

The summary follows the contract in [`summaries.md`](summaries.md), which every stage of the
suite shares and which defines the eight blocks, their order and the stability rules. Keys are
not repeated here, so that the two cannot drift apart; what follows is only what is particular
to the walker.

Each record names its haplotype by `seq_sha256` and `length` rather than carrying the sequence,
which is already in the FASTA. `--summary-seqs` restores the literal `seq` field.

Three per-accession decisions are worth knowing before reading a record:

- `provisional` is true when the neighbourhood holds more than
  `--subgraph-nodes-per-base-guard` unitigs per reference base. That is how tangled the graph
  is around the gene, which is what makes a walk ambiguous. `subgraph_fraction`, the
  neighbourhood's share of the graph, is reported beside it for information only.
- `prune` is how the pruning pass ended, or `null` when it did not run: `applied`,
  `nothing_to_remove`, `would_empty`, `refused_shorter_span`, `refused_no_haplotype`,
  `refused_trimmed` or `refused_error`, with the gene's depth, the floor, the unitigs kept and
  the reference span of each alignment. `refused_trimmed` is a realignment that is the first
  pass's haplotype with up to three bases off an end and nothing else changed, which is pure
  loss and can be a truncated stop codon. `refused_error` is a second pass that failed
  outright, which refuses the prune rather than failing the accession.
- `bias_uncorrected` is true when the prune ran and its result was thrown away, so the
  first-pass route stands. That route is the one the aligner produced on the unpruned graph,
  where its pull toward the reference is strongest, and `prefer_majority` skips sites whose
  alternative is below the variant floor, so nothing else corrects them either.
  `cohort.bias_uncorrected` counts it across the run, which is the form that is readable at
  5,000 accessions.

Where a fork on the published route still carries the rarer allele, `minority_left` says so and
counts the reason, and `cohort.accessions_with_minority_sites` carries the run total.

### bubbles.json

Computed from the graph rather than from the alignment, so it is written even for an accession
with no path, and most of its records are off the walked route. The VCF, the SNV table and the
alternatives come from forks **on** the route, not from this file.

- **One record per allele set.** The same two alleles can be reachable between several forks
  and merges; each record carries a `flanks` list, whose first entry is the `fork`/`merge` it
  also reports as its own fields. `n_bubbles` therefore counts sites.
- **The floor is the gene's own depth.** A branch must clear the strictest of `--prune-frac`
  of depth, `--min-abund-frac` of depth, and `--min-abund`, and the stage runs after the prune,
  on the graph the prune left. An accession that reconstructed nothing keeps the absolute
  floor, because without a route there is no depth to take a fraction of.
- **The two alleles must be two distinct unitigs**, neither of which is a flank, so an
  inverted repeat does not become a bubble whose two alleles are one sequence and its reverse
  complement.

A site can be in the VCF and not here. The two stages look for different things: this one finds
diamonds, where the alternatives rejoin in a single step, while the variant stages follow an
alternative for up to `--alt-max-detour` unitigs. On `SRR8209082` the real second population
rejoins after 31 unitigs, so it reaches the VCF, the SNV table and `alt_haplotypes.fasta` and
has never been in `bubbles.json`.

### VCF and the SNV table

Positions are 1-based in the reference CDS, and `CHROM` is the gene. In a single-reference run
every accession's VCF shares one coordinate frame; in a multi-anchor run each uses its own
winning reference. Reference record ids must be usable as VCF contig names (letters, digits and
`!#$%&+./:;?@^_|~-`, plus `*` and `=` after the first character), and a run refuses rather than
writing a VCF no parser will read. A non-ACGT reference base is written as `N` in REF.

Frequencies are `ka:f` abundance ratios at the fork, a proxy for allele frequency, per site and
not phased. **The denominator counts only branches that return to the route**: an exit that
leaves the locus for good, a repeat, a paralog or a contaminant sharing the fork's k-mer, is
not an allele of this gene, and counting it deflates every frequency at the site.

An alternative carrying several substitutions is called site by site rather than discarded. In
a k=31 graph two substitutions within 30 nt share one bubble and arrive as one alternative, and
two changes in one codon are the common case. Those rows share a `BLOCK` (VCF) or `block`
(TSV), which says they **are** on one molecule, unlike two rows that merely share a position.
Past 8 linked substitutions the alternative is a rerouted stretch of graph rather than a few
linked changes, so it stays in `alt_haplotypes.fasta` and is counted on stderr.

There is no length limit on the SNV mapping; it is bounded by memory, at 256 MiB, which refuses
the case that is actually expensive, a short but highly divergent pair, and says so.

For a cohort, `snvs.tsv` is the output that scales: it is sparse, concatenates across shards
and loads straight into DuckDB, while a merged multi-sample VCF becomes a huge, nearly empty
matrix. A worked example is in [`examples/vcf_example.md`](examples/vcf_example.md).

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
| `error` | anything else, including an internal panic in a fetch or compute worker (the panic ends that accession and names itself in `detail`, rather than ending the run): a network or cache failure, a corrupt file, a cached object that is truncated (a link naming a record the file does not hold, which is what a copy cut at a zstd frame boundary looks like), a GAF row that is not a walk of the graph it was aligned against, GraphAligner crashing, one GraphAligner call passing `--graphaligner-timeout`, or a program on `--graphaligner` that exits 0 without writing a GAF at all (1.0.17b creates the file before it aligns, so that is not an alignment result) | error |

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

4 is separate from 1 and 3 on purpose: the run produced real output, so a wrapper that wants
to keep `reconstructed.fasta` and carry on can tell this case from a run that produced nothing.
The accessions that failed are listed in `<query stem>.skipped.tsv` with their `error` outcome
and detail. Through 0.1.x a single success was enough for the whole run to exit 0, so a wrapper
written against that behaviour needs to handle 4.

Usage and I/O errors also exit 1, and `--allow-empty` does not cover them.

With `--fetch-only`, success means every accession is now in the cache. `--allow-empty` accepts
accessions with no Logan data or over the size cap; any other failure exits 1, or 3 when every
accession failed.

If the walker's standard output is closed early (`logan-walker ... | head`), the walker stops
on SIGPIPE, as other Unix tools do, and exits 141, removing its scratch directory on the way
out (see [Scratch files](#scratch-files)).

## Memory

Nothing in the walk scales with the size of an accession. The compressed object streams to
disk as it arrives, the adjacency spills to disk once it outgrows 4 MiB, and what stays in RAM
is the gene's neighbourhood.

So the bound is a product of two settings, neither of which is `--max-size-mb`:

```text
per compute worker   --max-subgraph-nodes x 150 bytes, plus the parser's flush buffers
at the defaults      0.16 GB each, so -w 8 is 1.24 GB
in-flight accessions one 128 KiB stream buffer each, so 16 of them is 2 MiB
```

The same numbers hold at a 4 GB cap and at a 28 GB one. The walker prints the bound at startup,
capped by the work available: a run of n accessions holds at most min(n, `-w`) graphs, so a
single-accession run prints 0.16 GB, and `--fetch-only` prints no graphs at all because it
never parses.

**That printed bound is also what the parser enforces.** Both passes measure what the graph
holds and refuse the accession as `too_large` once it crosses the budget, so a file dense
enough to build a graph many times its compressed size cannot exceed what the run announced.
The arrays grow by doubling and are measured periodically, so the refusal lands within one
doubling of the budget rather than exactly on it, which is why the line says "up to".

To contain a run, lower `-w` and `--fetch-workers`; they multiply. Excluded from the bound:
GraphAligner, the scratch GFA and the output buffers.

### When the gene does not bound it

Two settings turn the bound back into a function of accession size, and the fallback is
`--max-size-mb` x 1.7, a ratio measured across accessions from 4.8 MB to 348 MB:

- `--max-subgraph-nodes 0`, which removes the neighbourhood cap;
- no spool directory, which no CLI run can be, because the walker always resolves one.

If adjacency *cannot* spill, because the spool directory is unwritable, the parser keeps it in
RAM and still holds it to the budget: the accession is refused as `too_large`, beside a warning
naming the directory that failed, rather than quietly growing into fifty times the announced
memory. Spilling trades RAM for disk, and the disk is bounded too. See [Disk](#disk).

One allocator note, because it is load-bearing on long lists: glibc raises its own mmap
threshold each time a large block is freed, which would leave every graph resident from the
second accession onward. The walker pins the threshold at 1 MiB, so each graph returns to the
OS when it is dropped.

## Disk

What the spool costs is predictable. It is two files of `u32` and nothing else: one offset per
unitig plus a terminator, and one entry per edge. No sequences, no abundances, and no slack
from an allocator growing by doubling, because it is written out rather than grown.

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

## Scratch files

Each accession's subgraph GFA, alignment GAF and a temporary copy of the reference are written
to a scratch directory while it is aligned. Every run makes its own,
`logan-walker.<pid>.<suffix>` with mode 0700, under `--scratch-dir`. The default is `/dev/shm`,
which is tmpfs, so the files use RAM; where it does not exist (macOS), the system temp
directory is used. Files are removed as each accession finishes and the directory when the run
ends. Set `LOGAN_KEEP_SCRATCH=1` to keep it for debugging.

`/dev/shm` is only 64 MB in a default Docker container, so point `--scratch-dir` elsewhere
there.

### Signals

The directory is removed on any signal that ends a run: SIGINT, SIGTERM, SIGHUP, SIGPIPE,
SIGQUIT, SIGUSR1, SIGUSR2, SIGXCPU and SIGXFSZ. That list is what a scheduler actually sends:
SIGUSR1 and SIGUSR2 are how SGE and LSF announce a soft wall-time limit, SIGXCPU comes from a
CPU ulimit or a cgroup, and SIGPIPE is what `logan-walker ... | head` produces, which is why a
pipeline exits 141.

SIGALRM, SIGVTALRM and SIGPROF are deliberately not in the list. Nothing sends those to
announce a limit; what does send them is a timer or a sampling profiler belonging to the
process itself, and catching them would end a profiled run rather than clean up after one.

A signal the walker inherited as ignored stays ignored, so `nohup logan-walker ... &` survives
the SSH session dropping, and GraphAligner children inherit the ignore too.

**Any of those signals also stops the GraphAligner processes the run started.** They are the
only processes the walker starts, and the scratch files are what they are working on. Each
aligner runs in its own process group, which is what lets one signal reach a wrapper script's
children as well, and on Linux it is additionally marked to die with the walker, so even
SIGKILL and the out-of-memory killer take the aligners with them.

### Sweeping another run's directory

A run killed outright cannot clean up, so the next run with the same `--scratch-dir` removes
any directory whose process no longer exists. Each run claims its directory with a `.lock` it
holds open, and a directory is swept only when that lock names this machine and is free.
Liveness by PID alone would answer only about the caller's PID namespace on the caller's host,
so a run in another container or on another node sharing a `--scratch-dir` would look dead and
have its in-flight files deleted. **Directories left by 0.1.x carry no lock and are never
swept**, so clear those out once by hand.

## Logan unitig files

Logan's unitigs are BCALM2 output from a k = 31 compacted de Bruijn graph. A header looks like

```text
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
- `ka:f` **is clipped at 25,000.** Measured on three cached accessions: `ERR2028820` has 664
  unitigs at exactly 25000.0 against 4 at the next value down, `SRR8209082` has 802 against 24, and
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

A 404 from S3 means the accession has no Logan unitigs and is not retried; 408, 429, 5xx and
transport errors are retried up to `--fetch-attempts` times (default 6). The objects are
unversioned and mutable, so each run records what it read: see
[Input provenance](walker.md#input-provenance).

[`notes/contigs_vs_unitigs.md`](notes/contigs_vs_unitigs.md) explains why the walker uses
unitigs rather than Logan's assembled contigs. The default `-k 21` comes from the simulation in
[`../bench/simulated_kmer_sensitivity/`](../bench/simulated_kmer_sensitivity/README.md).
