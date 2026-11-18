# Golden-truth fixtures

Concordance baselines for `logan-walker`, checked by
[`../check_golden.sh`](../check_golden.sh). There are two sets, because the walker's defaults
changed after the first one was captured and both behaviours are worth gating:

- **old-behaviour** (this directory): the primary haplotype each accession must reproduce byte
  for byte under `--prune-frac 0.0 --no-prefer-majority --no-emit-alts`.
- **default-flag** (`default/`): every output the walk produces with no flags but
  `--emit-vcf --emit-snv-table`, which is what users actually run. A writer that is off by
  default produces nothing to gate.

All fixtures were captured under GraphAligner 1.0.17b, which is recorded in
`GRAPHALIGNER_VERSION`; the gate warns when it runs under a different one.

## The old-behaviour set

| accession | gene | query | primary | reference | gated |
|---|---|---|---|---|---|
| `ERR2028820` | PfCRT | `queries/pfcrt_only.fasta` | 1,275 nt | 1,275 nt | byte-exact |
| `SRR12363134` | *C. auris* ERG11 | `queries/markers/candida_auris_erg11.fasta` | 1,575 nt | 1,575 nt | byte-exact |
| `SRR8209082` | HIV-1 *pol* | `queries/query_hiv1_pol.fasta` | 2,740 nt | 2,739 nt | reported only |

Each `*.primary.fasta` is the walker's FASTA output verbatim. `SRR8209082` is reported rather
than gated here because its **unpruned** route is not reproducible across GraphAligner versions:
its 34,582-node neighbourhood is 42% of the whole graph, and 1.0.17b and 1.0.20 take different
routes through it. The default pruning is what stabilises it, which is why the other set can
gate it.

## The default-flag set

| accession | primary | alts | bubbles | SNV rows | VCF | summary |
|---|---|---|---|---|---|---|
| `ERR2028820` | 1,275 nt | none | 1 | none | none | gated |
| `SRR12363134` | 1,575 nt | none | 0 | none | none | gated |
| `SRR8209082` | 2,738 nt | 1 (`vaf=0.310`) | 0 | 1 | 1 row | gated |

Three things worth knowing before trusting a green run:

1. **`SRR8209082` is gated here, unlike above**, on the evidence of three independent default
   runs that agreed byte for byte across all five outputs. That covers one machine at one
   thread; GraphAligner is not deterministic on a graph this size with more than one thread. If
   it flaps in CI, demote it to the script's informational section rather than removing the
   gate.
2. **Only `SRR8209082`'s sequence changes between the two sets**, 2,738 nt against 2,740 nt. For
   the other two the default walk reproduces the same primary despite pruning being active (479
   of 747 and 135 of 296 unitigs kept), so there these fixtures gate the header statistics,
   bubbles, alts and SNV table rather than the sequence.
3. **None of the three exercises the majority correction.** `SRR8209082`'s walk already took
   the majority branch (`min_vaf=0.690` against the alt's `vaf=0.310`), `ERR2028820`'s four
   candidate alleles are all rejected for unrelated reasons, and `SRR12363134` has no on-path
   fork. A green run here says nothing about `prefer_majority`; the only coverage of that path
   is `walker/tests/cli.rs`.

`alts.fasta` is empty for the two clonal isolates, and a **missing** VCF fixture is itself an
assertion: those two recover no SNVs, so a VCF appearing is a real change.

Why `SRR8209082` has no bubbles while carrying a real second population: its one genuine
alternative rejoins the route after 31 unitigs (`detour=31` on its alt record), and this stage
finds only diamonds, where the alternatives rejoin in a single step. It is therefore in the
VCF, the SNV table and the alt FASTA, and has never been in `bubbles.json`.

## What the gate asserts

The old-behaviour set byte-compares the primary sequence and nothing else. The default-flag set
compares the full primary record, the alt haplotypes, the SNV table, the VCF, `bubbles.json` and
the summary.

Summaries are compared through [`../mask_summary.py`](../mask_summary.py), which replaces what
cannot reproduce: the clock (`run.started`, `run.wall_seconds`, `run.stage_seconds`, each
record's `timing` and `wall_seconds`), the machine (`tool.external[].path`, `inputs.source.id`
and `kind`, since a cached CI run and a cold local one differ there), the versions under test
(`tool.version`, `tool.disassembler`), and each record's `payload.source`, `payload.etag`,
`payload.last_modified` and `payload.fetched_at`. A warm cache has no HTTP response behind it,
so the middle two are absent on a cached run and present on a cold one; the mask inserts all
four whether or not they were there, so the two forms compare equal.

`payload.bytes` and `payload.sha256` are **gated**, not masked. They are identical on every
ingest path, because a cache file is byte-identical to the S3 object, so gating them makes the
run fail when Logan regenerates a gated accession. The `u/` objects are unversioned and mutable
and a 5,000-key sample on 2026-11-16 found 3.1% rewritten since the 2024-03 load, so that is a
real event rather than a hypothetical one. It is a different failure from a walk regression, and
the diff says which.

Everything else in the summary is gated: every count in `discards` and `cohort`, every
`seq_sha256`, every setting, and the block order. No normalisation is applied anywhere. The VCF
carries no `##fileDate` and no command line, and the FASTA headers hold only walk statistics, so
every file is byte-comparable as written.

A clean run leaves no `.partial.fasta` or `.journal.jsonl`: the walker writes both as it goes
and removes them once the real outputs are on disk, so finding one in an output directory means
a run is in flight or one died there.

## Regenerating a fixture

Run the gate with `LOGAN_GOLDEN_WORK` set so its outputs survive, then copy what changed:

```bash
LOGAN_GOLDEN_WORK=/tmp/golden bench/check_golden.sh
cp /tmp/golden/def-<tag>/<query stem>.summary.json default/<tag>.summary.json
```

The other default-flag fixtures come from the same directory: `primary.fasta`,
`alt_haplotypes.fasta`, `snvs.tsv`, `<acc>.vcf` and `bubbles.json` become
`default/<tag>.{primary.fasta,alts.fasta,snvs.tsv,vcf,bubbles.json}`. Without
`LOGAN_GOLDEN_WORK` the gate works in a temp directory it deletes on exit and there is nothing
to copy.

For the old-behaviour set, run the walker with `--prune-frac 0.0 --no-prefer-majority
--no-emit-alts` under the GraphAligner version in `environment.yml`, save the primary, and
update `GRAPHALIGNER_VERSION`.

Regenerate from a **cold** run when the thing being re-captured is a summary, so that
`payload` is populated from a real HTTP response rather than from a cache.

## Why these accessions

They are the two named in the original concordance criterion, plus `SRR12363134`. Neither of
the first two appears in `hiv1_50.txt` or `pfcrt_50.txt`, so the shipped benchmark harness does
not exercise them.

`SRR8209082` is worth treating with care as a performance baseline. Its 81,534 unitigs are
roughly 37x below the median of a 50-accession ERG11 run, and its 5-hop subgraph pulls in
34,582 nodes, 42% of its entire graph. It is fast and tangled, which is not the common case.

## History

The `*.golden.json` files in this directory are from the original capture, 2026-09-10, by the
Python prototype walker under GraphAligner 1.0.20. Nothing in the current tree reads them; they
are kept as a record, along with the graph statistics they carry (`ERR2028820`: 1,277,362
unitigs, 76 anchors, 356 subgraph nodes; `SRR8209082`: 81,534 unitigs, 4,102 anchors, 26,323
subgraph nodes at the time). The verifier that read them, `verify_golden.py`, was retired with
that walker.

| when | what changed |
|---|---|
| 2026-10-18 | all three fixtures recaptured from the Rust walker under GraphAligner 1.0.17b. The two gated primaries gained their final base, the stop codon 1.0.19 and 1.0.20 clip |
| 2026-11-13 | the default-flag set added, because the old gate switches off three of the four steps that shape a default primary and compares only the bare sequence |
| 2026-11-18 | the summaries moved to `"schema": 2` and gained `prune_min_abund`, `abundance_missing`, `no_links`, `mixed`, `two_allele_sites` and `min_site_vaf`; the primary headers gained `min_vaf=`. The majority correction, which had been discarded unless the aligner's row spanned the whole path, now runs under default flags, so `SRR8209082`'s route moved from 2,506 to 2,479 unitigs. Its `bubbles.json` went from 3,632 bubbles to none, under a depth-relative floor: every one of the 3,632 had its strongest allele at 66x or less against a 15,000x gene, a median of 11.9x. The file is now eight lines and committed rather than gated by a checksum |

Each of those recaptures was verified to contain only its intended diff before the fixtures were
replaced.
