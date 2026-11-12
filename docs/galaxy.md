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
| `disassembler_harvest` | `harvest --shards-dir` | the shard collection from the IUC `kmindex_query` tool | outputs per-target accession lists and score files as collections |
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
  memory fits `GALAXY_MEMORY_MB` when the runner sets it;
- keeps `-t` at 1, because GraphAligner is not reproducible with more threads;
- links the reference under its dataset name, because the walker names its summary and
  multi-anchor tags after the reference file's stem.

A job that reconstructs nothing fails unless `--allow-empty` is set, so a workflow decides
whether an empty shard is a soft or a hard failure. `--allow-empty` accepts only accessions where
nothing was found (no Logan data, over the size cap, no anchors, no path); a job in which every
accession errored fails with exit code 3 regardless (see
[`walker.md`](walker.md#exit-status)).

## Workflows

Three workflows are planned, one per composite:

```
from-kmindex:   kmindex_query -> harvest -> stratify (optional) -> walker (scattered)
                -> sieve (optional) -> msa
from-lexicmap:  LexicMap search -> lexicmap -> sieve (optional) -> msa
reconstruct:    accession list -> walker (scattered) -> sieve (optional) -> msa
```

Command-line flags become workflow choices: `--sieve` is whether the sieve step is connected,
and `--skip-tree` is where the workflow stops. `--skip-harvest` becomes a choice of entry point
(start from an accession list instead of a kmindex collection). For cohorts, `snvs.tsv` is the
scalable variant output; a merged multi-sample VCF can be built with the IUC `bcftools merge`
tool from the walker's per-accession VCFs.

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

## Open questions

- **Splitting `msa`.** The disassembler-specific part of `msa` is the reference-guided,
  in-frame preparation; MAFFT and FastTree already have well-maintained IUC tools. Splitting it
  into a preparation tool plus the IUC tools would make each step inspectable and turn
  `--skip-tree` into simply not connecting FastTree. That needs a preparation-only mode in the
  CLI first, so `disassembler_msa` is one tool for now.
- **A dedicated Logan Galaxy instance.** usegalaxy.org charges user storage and rate-limits
  API calls. A dedicated instance would make `fetch` worth wrapping and could allow a larger
  default `--max-size-mb`.
- **Metadata.** At SARS-CoV-2 scale a random subset of millions of candidates is not very
  informative; joining SRA metadata (date, platform, country) before stratifying would make
  it so. That belongs in the workflow, not in disassembler.
