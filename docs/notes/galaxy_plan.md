# Galaxy workflows: the plan

**Status:** written 2026-11-18, alongside the six tool wrappers. Nothing here is built yet. What
the wrappers do today is in [`../galaxy.md`](../galaxy.md).

## The three workflows

Three workflows are planned, one per composite:

```text
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

## Questions still open

- **A dedicated Logan Galaxy instance.** usegalaxy.org charges user storage and rate-limits
  API calls. A dedicated instance would make `fetch` worth wrapping and could allow a larger
  default `--max-size-mb`.
- **Metadata.** At SARS-CoV-2 scale a random subset of millions of candidates is not very
  informative; joining SRA metadata (date, platform, country) before stratifying would make
  it so. That belongs in the workflow, not in disassembler.
