# disassembler on Galaxy

The disassembler suite is being wrapped as Galaxy tools so that a user can go from a kmindex or
LexicMap search to a codon MSA and tree without leaving Galaxy. The six tool wrappers are under
review in [galaxyproject/tools-iuc#8800](https://github.com/galaxyproject/tools-iuc/pull/8800)
and will publish to the Tool Shed as `suite_disassembler`. The workflows are not written yet.

Each tool is a thin wrapper around the same subcommand a terminal user runs, so the two stay in
sync, and the tools install the Bioconda packages described in [`packaging.md`](packaging.md).
Options that only make sense outside Galaxy (local caches, S3 fallbacks, API keys, binary path
overrides, output paths) are not exposed; they stay on the command line.

## Tools

| tool | wraps | inputs | notes |
|---|---|---|---|
| `disassembler_harvest` | `harvest --shards-dir` | the shard collection from the IUC `kmindex_query` tool, in either `json` or `json_vec` format | outputs per-target accession lists and score files as collections |
| `disassembler_stratify` | `stratify` | harvest's per-target scores collection | score bins as a repeat; the upper bound defaults to 1.01 so a perfect score is included |
| `disassembler_walker` | `logan-walker` | an accession list dataset or a single accession, plus the reference CDS | see below |
| `disassembler_lexicmap` | `lexicmap -f` | a LexicMap search table with all columns, plus the reference CDS | |
| `disassembler_sieve` | `sieve` | a FASTA, plus a reference or an explicit length gate | gene-scale defaults come from the subcommand |
| `disassembler_msa` | `msa` | haplotypes FASTA, plus the reference CDS | pins MAFFT 7.526 (7.525 has a bug in `--addfragments`, the path `msa` uses) |

Two subcommands are not tools:

- **`fetch`** fills a local cache for `--unitigs-cache`. Galaxy already keeps tool outputs as
  datasets, and a collection of unitig files (median about 50 MB each) would count against the
  user's storage quota. It could return as a data manager on a dedicated Logan Galaxy instance.
- **The composites** (`reconstruct`, `from-kmindex`, `from-lexicmap`) become workflows.

### The walker tool

The walker takes an accession list dataset, so how the list is split before it reaches the tool
decides the scatter strategy (below). The wrapper:

- streams unitigs from the public Logan bucket inside the job, so the Galaxy server needs
  outbound network access;
- sets `--scratch-dir` inside the job directory, because the walker's default `/dev/shm` is
  shared by every job on a node and is only 64 MB inside a default Docker container;
- sets `-w` from `GALAXY_SLOTS`, and trims `--fetch-workers` so the walker's worst-case buffer
  memory fits `GALAXY_MEMORY_MB` when the runner sets it. That worst case is a good deal
  smaller from 0.2.0: the fetch-to-compute queue no longer holds buffers, and an accession read
  from `--unitigs-cache` is streamed from its file handle rather than held, so the wrapper's
  formula can be recomputed against the bound the walker now prints at startup;
- keeps `-t` at 1, because GraphAligner is not reproducible with more threads;
- links the reference under its dataset name, because the walker names its summary and
  multi-anchor tags after the reference file's stem.

A job that reconstructs nothing fails unless `--allow-empty` is set, so a workflow decides
whether an empty shard is a soft or a hard failure. `--allow-empty` accepts only accessions where
nothing was found (no Logan data, over the size cap, no anchors, no path); a job in which every
accession errored fails with exit code 3 regardless (see
[`walker_reference.md`](walker_reference.md#exit-status)).

### Scatter

Galaxy's **Split file** tool splits the accession list, and the split decides the strategy:

1. **One accession per job:** per-accession retry and provenance, but scheduler overhead
   becomes a problem above roughly 10,000 accessions.
2. **N accessions per job:** each job runs the walker's own list loop. Practical well past
   100,000 accessions, and the only realistic option at SARS-CoV-2 scale (10^5 to 10^6).
3. **One job:** the whole list in one walker run. No per-accession isolation.

Downstream, the IUC **Collapse Collection** tool concatenates the per-job FASTAs before the sieve
or MSA. With one accession per job, either set `--allow-empty` or filter failed datasets before
collapsing.

The workflows that will tie these tools together, and the questions still open about them, are
in [`notes/galaxy_plan.md`](notes/galaxy_plan.md).
