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
   low-depth branches and re-aligns (`--prune-frac`), puts the majority allele on the primary at
   two-allele sites (`--prefer-majority`), and writes the remaining intra-host alleles as
   alternative haplotypes (`--emit-alts`). Why these exist is in
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

GraphAligner 1.0.17b must be on `PATH` (see
[`issue_graphaligner_end_clipping.md`](issue_graphaligner_end_clipping.md) for why that version).
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

## Outputs

| file | when | contents |
|---|---|---|
| `-f` (default `reconstructed_haplotypes.fasta`) | always | one primary per reconstructed accession, `>ACC_primary cov=… pident=… nodes=…` (plus `anchor=` in multi-anchor mode) |
| `<out-dir>/alt_haplotypes.fasta` | unless `--no-emit-alts` | intra-host alternatives, `>ACC_altN …` with their VAF and fork, at most 32 per accession |
| `<out-dir>/bubbles.json` | always | bubbles in the local subgraph |
| `<out-dir>/<query stem>_walker_summary.json` | always | per-accession status, metrics and candidates |
| `<out-dir>/<ACC>.vcf` | `--emit-vcf` | intra-host SNVs of one accession, written only when it has at least one |
| `<out-dir>/snvs.tsv` | `--emit-snv-table` | one row per (accession, site, allele) across the run |

`cov` is the walked length over the reference length, and `pident` the identity of the
alignment. When a haplotype's identity was inherited from the route before an allele was spliced
in, the header says `inherited_pident=` instead.

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

```bash
# Fill the cache, no analysis (wraps logan-walker --fetch-only)
disassembler fetch -l accessions.txt --unitigs-cache ./unitigs/

# Then run as often as needed with no network
disassembler walker -l accessions.txt -q queries/markers/zika_ns5.fasta -o out/run1 \
    --unitigs-cache ./unitigs/ --offline
```

`fetch` skips accessions already in the cache; delete a file to refetch it. The composites take
the same two flags.

## Options

`logan-walker` needs one of `-a`/`-l` and at least one `-q` (`--fetch-only` needs no `-q`).
Options that are on by default are turned off with a `--no-…` form.

### Input and output

| flag | default | meaning |
|---|---|---|
| `-a, --accession <ACC>` | | one SRA accession |
| `-l, --accession-list <FILE>` | | file with one accession per line |
| `-q, --query <FASTA>` | | reference CDS; repeatable, and more than one record switches on multi-anchor mode |
| `-o, --out-dir <DIR>` | `walker_output` | output directory |
| `-f, --out-fasta <FILE>` | `reconstructed_haplotypes.fasta` | primary haplotypes |
| `--scratch-dir <DIR>` | `/dev/shm/logan-walker` | per-accession GFA and GAF files (not shown by `--help`) |
| `--allow-empty` | off | exit 0 even if no accession was reconstructed |

### Anchoring, fetching and walking

| flag | default | meaning |
|---|---|---|
| `-k, --k <N>` | `21` | anchor k-mer size, 1 to 32 (the Logan graph's own k stays 31) |
| `--min-shared-kmers <N>` | `3` | shared k-mers a unitig needs to be an anchor, 1 to 5 |
| `--hops <N>` | `5` | hop radius of the subgraph around the anchors |
| `--max-size-mb <N>` | `1024` | refuse accessions above this compressed size |
| `--fetch-workers <N>` | `8` | parallel S3 downloads (bandwidth saturates around 8) |
| `-w, --workers <N>` | `8` | compute workers (see [Memory](#memory) before raising) |
| `-t, --threads <N>` | `1` | threads per GraphAligner call; above 1 is not reproducible |
| `--min-cov <F>` | `0.85` | minimum fraction of the reference a primary must cover |
| `--primary-metric <M>` | `matched` | multi-anchor primary selection: `matched`, `cov`, `pident` or `cov_x_pident` |
| `--prune-frac <F>` | `0.01` | drop unitigs below this fraction of the gene's depth and re-align; `0.0` disables |
| `--no-prefer-majority` | | keep the reference-matching allele on the primary at two-allele sites |

### Intra-host variants

| flag | default | meaning |
|---|---|---|
| `--min-abund <F>` | `5.0` | minimum `ka:f` for a branch to count as an allele |
| `--min-vaf <F>` | `0.10` | minimum variant allele frequency for an alternative |
| `--min-abund-frac <F>` | `0.0` (off) | raise the floor to this fraction of the gene's depth (median `ka:f` on the path) |
| `--alt-max-detour <N>` | `100` | unitigs an alternative may traverse before rejoining the path |
| `--no-indel-alts` | | substitutions only |
| `--allow-frameshift-alts` | off | also accept indels whose length is not a multiple of 3 |
| `--max-alt-indel <N>` | `30` | largest length change treated as an indel rather than a rerouting |
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

A compressed unitig file is held in three places at once: queued in a bounded channel (4 slots,
not a flag), held by a fetch worker waiting to hand it over (`--fetch-workers`), and held by a
compute worker for the whole parse (`-w`). At the defaults that is up to 4 + 8 + 8 = 20 files, so
**20 GB** at the 1 GB `--max-size-mb` cap. On top of that, each compute worker holds its parsed
graph, roughly 1.3 times the compressed size of the file it is working on.

This is a worst case: it assumes every accession sits at the cap. `ERR2028820` is 33 MB, and the
largest *C. auris* ERG11 accession tested was 343 MB. The walker prints its buffer bound at
startup. To contain a run, lower `--max-size-mb` and the two worker counts together, since they
multiply. `-w 16 --fetch-workers 16`, for example, allows 36 buffers.

## Scratch files

Each accession's subgraph GFA, alignment GAF and a temporary copy of the reference go to
`--scratch-dir`, default `/dev/shm/logan-walker`. That is tmpfs, so the files use RAM. Their names
carry the process ID and a per-task counter, so concurrent runs can share the directory, and they
are removed as each accession finishes. Set `LOGAN_KEEP_SCRATCH=1` to keep them for debugging.

On a shared machine or inside a container, set `--scratch-dir` to a private directory: `/dev/shm`
is shared by every user on the host and is only 64 MB in a default Docker container.

## Reproducibility

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

A run that reconstructs no accession at all exits non-zero, so a scheduler cannot mistake an
all-failed shard for a good one; `--allow-empty` makes it exit 0. Individual accessions that fail
(no Logan record, no anchors, no path) are reported in the summary and do not stop the run.

## Logan unitig files

Logan's unitigs are BCALM2 output from a k = 31 compacted de Bruijn graph. A header looks like

```
>ERR2028820_2933 ka:f:5.3 L:-:443893:+ L:+:1239391:-
```

- The unitig id is the text after the last `_`. Ids run from 0 to n-1 with no gaps.
- `ka:f` is the unitig's mean k-mer abundance; the walker uses it as depth.
- `L:s1:t:s2` in unitig `u`'s header is the traversal from `u` on strand `s1` to `t` on strand
  `s2`. Every edge appears in both endpoints' headers, so the graph is closed under reversal.
- Consecutive unitigs overlap by k - 1 = 30 bases, which the walker removes when it stitches a
  path into a sequence.
- Sequences are plain uppercase `ACGT`; anything else is an error.

A 404 from S3 means the accession has no Logan unitigs and is not retried; server errors are
retried three times.

[`contigs_vs_unitigs.md`](contigs_vs_unitigs.md) explains why the walker uses unitigs rather
than Logan's assembled contigs. The default `-k 21` comes from the simulation in
[`../bench/simulated_kmer_sensitivity/`](../bench/simulated_kmer_sensitivity/README.md).
