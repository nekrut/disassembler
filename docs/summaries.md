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
  "schema": 2,

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
    "stage_seconds": { },              // the stage's own internal phases; a phase that no
                                       // longer exists reports 0 rather than disappearing,
                                       // so two schema-1 summaries have the same keys
    "n_records": 50,
    "n_ok": 45
  },

  "discards": { },                     // every rejection, keyed by reason, counts only

  "cohort": { },                       // population statistics; omitted where meaningless

  "records": [                       // one small row per unit of work
    {
      "acc": "SRR8209082",
      "status": "success",
      "payload": {                     // what this unit of work was built from
        "source": "s3",                // s3 | cache | none
        "bytes": 1861058,
        "sha256": "...",               // of the bytes as read
        "etag": "d2dd288b238c...",     // live fetches only
        "last_modified": "2024-03-31T17:38:36Z",
        "fetched_at": "2026-11-16T16:47:50Z"
      }
    }
  ]
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

`source` describes the run. It cannot describe a single unit of work, because a run with a
cache and no `--offline` reports `cache+s3` while each accession came from exactly one of the
two. That is what `records[].payload` is for.

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
diversity. Omitted by a stage with no population to describe.

### `records`

One row per unit of work: per accession, per target, per sequence. Sub-objects are written on
one line each, and a record stays inside the per-record budget below. A record names its
outputs rather than carrying them: a sequence appears as `seq_sha256` and `length`, never as
`seq`.

#### `records[].payload`

What the unit of work was built from, as opposed to what the stage made of it. Always present,
including as `{"source": "none"}` when no bytes were read, so that "nothing to record" and
"this version does not record it" stay different answers. Absent members are omitted rather
than written as null.

| member | meaning |
|---|---|
| `source` | `s3`, `cache`, or `none` when no bytes were read |
| `bytes` | size as read. Also set on a `too_large` skip: from `Content-Length` when the object was refused before transfer, and from the object itself when the refusal came later, such as a parsed graph over its memory budget |
| `sha256` | digest of the bytes exactly as read |
| `etag` | the object's `ETag`, quotes stripped, when a live fetch supplied one |
| `last_modified` | the object's `Last-Modified`, normalised to `YYYY-MM-DDTHH:MM:SSZ` |
| `fetched_at` | when this run read the object. Absent for a cache hit, whose bytes may predate the run by any amount |

`sha256` is the member that works on every path, because it comes from the bytes rather than
from the transport: a cache hit and an `--offline` run have no response to take headers from.
`etag` and `last_modified` are recorded in addition when they exist, because they are
comparable with what the bucket serves today.

This matters because the Logan objects at `logan-pub.s3.amazonaws.com/u/` are unversioned and
mutable, and they do get replaced. A 5,000-key sample across five key ranges on 2026-11-16
found 3.1% of objects rewritten since the 2024-03 bulk load, in waves dated 2025-12-10,
2025-12-14 and 2026-02-08. Without a recorded digest, a rerun after an upstream regeneration
produces different results with no detectable difference in any output.

## The size budget

| what | rule |
|---|---|
| run-level blocks | O(1), about 2 KB |
| `records[]` | O(units of work), 1 KB each |
| anything O(variants), O(unitigs), O(HSPs), O(bubbles) | **not in the summary** |

Measured on a real walker record (ERR2028820 against PfCRT, with the variant, prune, timing
and payload sub-objects all populated): **about 1.3 KB**, against roughly 4.9 KB for the same
record in 0.1.1. `payload` accounts for about 170 bytes of that when every member is present,
and about 30 when the bytes came from a cache. A record that produced nothing is far smaller,
because a sub-object that did not run is omitted rather than written as null, and in a cohort
like `bench/hiv1_50.txt` most records are of that kind. A 20,000-accession walker run therefore
lands near 23 MB against 130 MB.

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

**0.2.0 writes schema 2.** Almost everything it changed was additive, but two changes are not,
and the version exists to say so:

- the walker's `n_bubbles`, and the per-record `bubbles` count, now count **allele sets**
  rather than the (fork, merge) pairs an allele set was found between, so the same key returns
  a different number for the same accession (`SRR8209082`: 3,632 under schema 1, 0 under
  schema 2, after the floor change that went with it);
- `no_links` is a new value of `status`, and a reader that switches over a closed set of
  statuses does not know it.

A schema-1 reader therefore keeps parsing a schema-2 summary and draws the wrong conclusion
from it, which is the case the version number is for.

Summaries are byte-comparable between runs once the volatile fields are masked.
`bench/mask_summary.py`, which `bench/check_golden.sh` runs before diffing, replaces:

- `tool.version`, `tool.disassembler` and each `tool.external[].path`;
- `inputs.source.id` and `inputs.source.kind`;
- `run.started`, `run.wall_seconds` and `run.stage_seconds`;
- per record, `timing` and `wall_seconds`, and in `payload` the `source`, `etag`,
  `last_modified` and `fetched_at`.

Everything else is gated. `payload.bytes` and `payload.sha256` are deliberately not masked,
because they are the input the run consumed and must not move; neither is
`payload.revalidated`. Masking `tool.version` is what makes the fixtures survive a version
bump, so the gate does not catch a version regression on its own; the version floor in
`src/disassembler/resolve.py` is what does.

## Per-stage instantiation

| stage | `records[]` is | `cohort` carries | `discards` carries |
|---|---|---|---|
| `harvest` | one target | score quantiles | `failed`, `empty`, `bad_target` |
| `stratify` | one target | selection rate | under-filled bins |
| `walker` | one accession | cov and pident quantiles, distinct primaries, singleton rate | three groups: `accessions` (eight outcomes), `alleles` (seven reasons an alternative was rejected) and `minority_sites` (nine reasons the primary still carries the rarer allele) |
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

## Migration from schema 1 (0.2.0)

- **`n_bubbles` and `records[].bubbles` count allele sets.** One two-allele site is one
  record however many pairs of flanks it sits between, and each record lists them under
  `flanks`. The gate also became depth-relative and moved after the prune, so the counts fall
  on deep samples whatever a reader does with them. Comparisons of bubble counts across
  schema versions are meaningless; within a version they are now a diversity signal rather
  than a property of local topology.
- **`status` can be `no_links`**: the unitig object carried no `L:` adjacency at all, so
  nothing could be walked across unitigs. It is a negative like `no_path`, not an error, and
  it is counted in `discards.accessions.no_links`.
- Additive, listed here because consumers ask: `records[].variants` gained
  `two_allele_sites` and `min_site_vaf`; `cohort` gained `mixed` and `bias_uncorrected`;
  `records[]` gained `abundance_missing` (written only when non-zero), `majority_corrections`
  and `minority_rescues`; haplotypes gained `ref_cov` and `ref_span`; the walker's `settings`
  gained `prune_min_abund` and the lexicmap `settings` gained `min_conflict_hsp_bp`.

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
