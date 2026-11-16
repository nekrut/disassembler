# Summary outputs: one contract for every stage

Every stage of the suite writes a summary. Until 0.2.0 each one was written independently, so
two summaries from the same run shared no keys, no file-naming convention and no idea of what
belongs in a summary at all. This document is the contract they now share.

It is normative. A new stage implements this; it does not invent a shape.

## Why this exists

Three findings from the 0.2.0 audit, each verified in code:

**The machine-readable summary was a strict subset of stdout, and the part dropped was the part
about what went wrong.** `lexicmap` counts `n_malformed`, `n_other_query` and
`n_qlen_mismatch`, prints all three and records none. `harvest` counts `n_failed`, `n_empty`
and `n_bad_target`, prints all three and records none. The walker prints a full variant
rejection breakdown and a majority-correction count and records neither. A run that silently
discarded most of its input and a clean run produced summaries that looked the same.

**The two pathways disagreed arbitrarily.** The walker summary carried full provenance
(versions, references, 19 settings) and no population statistics. The LexicMap summary carried
quantiles, tiers and entropy and no provenance at all: not its own version, not the input
identity, not one of the eight flags that change its answer. Neither difference was a
considered decision.

**Summaries carried unbounded data.** Measured on three real runs, between 50.1% and 71.4% of a
walker summary was embedded haplotype sequence, a verbatim second copy of
`reconstructed.fasta`. Extrapolated to a 20,000-accession run that is roughly 130 MB of JSON,
about 67 MB of it duplicate FASTA.

## The four rules

1. **One shape.** The top-level blocks below appear in every summary, in this order, with these
   names. A difference between two summaries means a real difference between the stages.
2. **Bounded by construction.** A summary holds run-level facts plus one small row per unit of
   work. It never holds sequence, and never holds anything that scales with variants, unitigs,
   HSPs or bubbles. The ceiling is **1 KB per record**, measured, with sub-objects on one line
   and the ones that did not run omitted entirely.
3. **Everything printed is also recorded.** If a line goes to stdout or stderr describing a
   decision, a rejection or a discard, it has a field in the summary. A reader of the summary
   alone must be able to reconstruct what the operator saw.
4. **A summary is always written**, including on failure, with a status naming what failed.
   Absence of the file is never the error signal. One exception: a stage that exits
   `EX_TEMPFAIL` because its input is not ready yet writes nothing, because no run was
   attempted and the caller is expected to retry.

## The schema

```jsonc
{
  "schema": 1,

  "tool": {
    "name": "logan-walker",
    "version": "0.2.0",
    "disassembler": "0.2.0",          // null when the stage runs standalone
    "external": [                      // every external binary whose version changes results
      {"name": "GraphAligner", "version": "1.0.17b", "path": "/opt/conda/bin/GraphAligner"}
    ]
  },

  "inputs": {
    "references": [                    // null for a stage with no reference
      {"id": "PfCRT", "length": 1275, "sha256": "..."}
    ],
    "source": {
      "kind": "s3",                    // s3 | file | galaxy | stdin | cache
      "id": "s3://logan-pub/u/",       // URI, path, or Galaxy dataset id
      "n_requested": 50                // units of work asked for, before any filtering
    }
  },

  "settings": { },                     // every flag that changes the answer, nothing else

  "run": {
    "status": "ok",                    // ok | partial | failed
    "exit_code": 0,
    "started": "2026-11-15T09:12:44Z",
    "wall_seconds": 182.437,
    "stage_seconds": { },              // the stage's own internal phases
    "n_records": 50,
    "n_ok": 45
  },

  "discards": { },                     // every rejection, keyed by reason, counts only

  "cohort": { },                       // population statistics; omitted where meaningless

  "records": [ ]                       // one small row per unit of work
}
```

### `tool`

What produced this file. `external` carries any binary whose version changes the result:
GraphAligner for the walker, MAFFT and FastTree for `msa`, `l1_pipeline_serial` for `sieve`.
Record the resolved path as well as the version, because the version alone does not say which
of several installed copies ran.

### `inputs`

What was consumed. `source.n_requested` is the count asked for, before any filtering, so
`n_requested` minus `run.n_records` is itself a discard that has to be explained in `discards`.
Reference sequences are identified by sha256, not by file path, because a path does not travel
between machines and a renamed or edited reference must not look identical.

### `settings`

Every flag that changes the answer, and nothing that does not. Workers, timeouts, output paths
and verbosity are excluded by design: they change cost, not results. The walker's existing
19-key block is the model. A new flag that changes results and is not added here is a bug.

### `run`

`status` is `ok` when every unit of work succeeded, `partial` when some did and some did not,
and `failed` when the stage could not do its job at all. A `failed` summary still carries
`tool`, `inputs` and `settings`, plus an `error` string saying what failed, so a failed run is
diagnosable from the summary alone.

### `discards`

Counts only, keyed by reason, flat where possible. Every reason the stage can reject something
for appears, including the ones that are zero, so that a consumer can tell "none of these
happened" from "this version did not measure it".

### `cohort`

Population statistics over whatever the stage produced: quantiles, counts of distinct outputs,
diversity, a top-N preview. Omitted by a stage with no population to describe.

### `records`

One row per unit of work: per accession, per target, per sequence. **400 bytes is the ceiling.**
Sub-objects are written on one line each. A record names its outputs rather than carrying them:
a sequence appears as `seq_sha256` and `length`, never as `seq`.

## The size budget

| what | rule |
|---|---|
| run-level blocks | O(1), about 2 KB |
| `records[]` | O(units of work), 1 KB each |
| `cohort.top_n` | bounded by `--summary-top-n`, default 10 |
| anything O(variants), O(unitigs), O(HSPs), O(bubbles) | **not in the summary** |

Measured on a real walker record (ERR2028820 against PfCRT, with the variant, prune and timing
sub-objects all populated): **1,132 bytes**, against roughly 4.9 KB for the same record in
0.1.1. A record that produced nothing is far smaller, because a sub-object that did not run is
omitted rather than written as null, and in a cohort like `bench/hiv1_50.txt` most records are
of that kind. A 20,000-accession walker run therefore lands near 23 MB against 130 MB.

Sequences are identified by `seq_sha256` (lowercase hex, sha256 of the uppercase sequence with
no header and no line breaks). This answers every question the old `seq` field was used for:
equality between runs, counting distinct outputs, and joining a record to its FASTA entry.
`--summary-seqs` restores the literal sequence for a caller that genuinely needs it.

## What stays a separate file

A file, not a summary field, when it scales with something other than the unit of work:

| file | why |
|---|---|
| `reconstructed.fasta`, `alt_haplotypes.fasta` | sequence |
| `<prefix>.clean.haplotypes.tsv` | carries a sequence per haplotype |
| `bubbles.json` | 3,632 bubbles for a single accession, about 1 MB |
| `snvs.tsv`, per-accession VCFs | scales with variants times accessions |
| `<target>.scores.json` | scales with accessions times targets |

Each of these has a **count** in `records[]`. The file holds the data; the summary holds the
shape of it.

Conversely, a per-unit table small enough to obey the size rule belongs **in** the summary,
not beside it. `<prefix>.cohort_ledger.tsv` is nine short columns per accession, so it becomes
`records[]`, and the TSV itself is demoted to an opt-in rendering behind `--cohort-ledger`:
writing it by default put every accession on disk twice. `<prefix>.flagged.tsv` was a strict
subset of that ledger and is removed outright.

That is the general rule for a file the summary now subsumes. If something downstream wants
the tabular form, it stays available behind a flag and must carry exactly the same rows, which
is a thing to test; it does not stay on by default.

## File naming

`<stem>.summary.json`, beside the stage's other outputs. This retires four conventions:
`summary.json`, `bins_summary.json`, `<prefix>.summary.json` and
`<query stem>_walker_summary.json`.

## Stability

`schema` is the contract. Within a major schema version, fields are only ever added, never
renamed and never removed, and a consumer that does not recognise a field ignores it. A change
that breaks a reader increments `schema`.

Summaries are byte-comparable between runs once the volatile fields are masked:
`run.started`, `run.wall_seconds` and `run.stage_seconds`, plus `tool.external[].path`.
`bench/check_golden.sh` masks exactly those and gates the rest.

## Per-stage instantiation

| stage | `records[]` is | `cohort` carries | `discards` carries |
|---|---|---|---|
| `harvest` | one target | score quantiles | `failed`, `empty`, `bad_target` |
| `stratify` | one target | selection rate | under-filled bins |
| `walker` | one accession | cov and pident quantiles, distinct primaries, singleton rate | fetch failures, `no_logan_data`, `too_large`, and the full variant rejection breakdown |
| `lexicmap` | one accession (the former ledger row) | tiers, quantiles, entropy (raw, full-coverage and folded), top N | `malformed`, `other_query`, `qlen_mismatch`, and the flagged tiers |
| `sieve` | one output haplotype | dedup ratio, admission rate | dedup drops, completeness-gate drops, table saturation |
| `msa` | one row | rows in and out, length agreement | rows dropped or lengthened by MAFFT |
| composites | one stage | none | which stage failed |

### Composite manifests

`reconstruct`, `from-kmindex` and `from-lexicmap` write a summary whose `records[]` is one entry
per stage, each naming that stage's own summary file:

```jsonc
"records": [
  {"stage": "harvest", "status": "ok", "wall_seconds": 31.2,
   "summary": "cohort/harvest.summary.json"},
  {"stage": "walker",  "status": "partial", "wall_seconds": 1821.9,
   "summary": "out/pfcrt.summary.json"}
]
```

A composite's `run.status` is the worst of its stages'. This is the only place a summary points
at another summary, and it points by relative path so the set can be moved together.

## Migration from the 0.1.x summaries

- **`seq` becomes `seq_sha256`.** Scripts that compared sequences between runs compare hashes
  instead; scripts that took `len(p["seq"])` use the existing `length` field. `--summary-seqs`
  restores the old field. Known consumers live in the `disassembler-dump` analyses.
- **`<query stem>_walker_summary.json` becomes `<query stem>.summary.json`.** The Galaxy
  wrappers name this file, so tools-iuc #8800 moves in step.
- **`<prefix>.flagged.tsv` is removed**, since `records[]` carries every flagged accession with
  its reason.
- The walker's old top-level `logan_walker`, `graphaligner`, `references`, `settings`, `total`,
  `success`, `results` and `bubble_events` keys move into the blocks above. `results` becomes
  `records`, `total` and `success` become `run.n_records` and `run.n_ok`, and `bubble_events`
  becomes a per-record `bubbles` count, which it already duplicated.
