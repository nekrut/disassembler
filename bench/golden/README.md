# Golden-truth fixtures

Concordance baselines for the path walker, covering success criteria 2 and 3 of the
optimization issue (`docs/issue_path_walker_optimization.md`).

## Current fixtures (recaptured 2026-10-18)

| Accession | Gene | Query | Primary | Reference | Gate |
|---|---|---|---|---|---|
| `ERR2028820` | PfCRT | `queries/pfcrt_only.fasta` | 1,275 nt | 1,275 nt | byte-exact |
| `SRR12363134` | *C. auris* ERG11 | `queries/markers/candida_auris_erg11.fasta` | 1,575 nt | 1,575 nt | byte-exact |
| `SRR8209082` | HIV-1 pol | `queries/query_hiv1_pol.fasta` | 2,740 nt | 2,739 nt | reported only |

Captured from the Rust walker (`logan-walker` 0.1.0, commit `f655edd`) under GraphAligner
1.0.17b, with `--prune-frac 0.0 --no-prefer-majority --no-emit-alts`. Each `*.primary.fasta` is
the walker's FASTA output verbatim. The aligner version is recorded in `GRAPHALIGNER_VERSION`, and
`check_golden.sh` warns when it runs under a different one.

The two gated primaries are the previous fixtures plus their final base, which restores the stop
codon that GraphAligner 1.0.19 and 1.0.20 clip (`docs/issue_graphaligner_end_clipping.md`). Under
1.0.20 the gate now fails on exactly that base. `SRR8209082` takes a different unpruned route under
1.0.17b than in the earlier capture; it is reported, not gated (see `check_golden.sh`).

## Original capture (2026-09-10)

| Accession | Gene | Query | Unitigs | Anchors | Subgraph | Primary | Bubbles |
|---|---|---|---|---|---|---|---|
| `ERR2028820` | PfCRT | `queries/pfcrt_only.fasta` | 1,277,362 | 76 | 356 | 1,274 nt | 202 |
| `SRR8209082` | HIV-1 pol | `queries/query_hiv1_pol.fasta` | 81,534 | 4,102 | 26,323 | 2,738 nt | 11,825 |

Captured on the unmodified Python walker at commit `81fc1fd`, GraphAligner 1.0.20, BLAST 2.17.0+.
The graph statistics above and the `*.golden.json` files are from this capture.

## Usage

```bash
bench/check_golden.sh              # exit 0 = concordant, 1 = regression
```

`check_golden.sh` runs the Rust `logan-walker` under old-behaviour flags and diffs the primary
byte-for-byte against the `*.primary.fasta` fixtures here (needs `GraphAligner` on `PATH` and S3
access). Until 2026-10-18 the `ERR2028820.pfcrt` and `SRR8209082.hiv1_pol` fixtures were the
original Python walker's output at commit `81fc1fd`, and `SRR12363134.candida_auris_erg11` was
captured from the Rust walker and cross-validated against the Python walker (identical 1,574 nt)
before that walker was removed. To regenerate a fixture, run the Rust walker with
`--prune-frac 0.0 --no-prefer-majority --no-emit-alts` under the GraphAligner version in
`environment.yml`, save the primary, and update `GRAPHALIGNER_VERSION`. The retired
Python-based `capture_golden.py`/`verify_golden.py` are under `../archive/`.

## What is asserted, and what is not

**What the live gate actually checks.** `check_golden.sh` byte-compares the primary sequence
against the `*.primary.fasta` fixtures and nothing else. The richer field list below describes
`verify_golden.py`, which moved to `../archive/` with the Python walker and is not runnable. The
`*.golden.json` files are therefore not read by anything in the current tree; they are kept as a
record of the original capture.

Every asserted field was verified byte-stable across two `PYTHONHASHSEED` values (1 and 7):

- graph shape (`n_unitigs`, `n_anchors`, `subgraph_nodes`)
- the primary haplotype, by SHA-256 of the sequence, plus `cov` and `pident`
- `n_bubbles`, and every bubble's fork node with its full branch list (node, `ka:f`, VAF)

The bubble `merge` node is deliberately **not** asserted. `phase_superbubbles` picks it with
`reconverged[0]` (in the original Python walker) from a `Counter` populated by iterating a `set`,
so it depends on string hash randomisation. Measured across seeds 1 and 7, it differed on
10,922 of 11,825 fork/merge pairs for `SRR8209082` and 18 of 202 for `ERR2028820`. The branch
calls themselves, which carry the biology, were identical in every case.

Timings are excluded. They are not reproducible and belong in benchmark output, not a fixture.

## Why these two accessions

They are the two named in the issue's concordance criterion. Note that neither appears in
`hiv1_50.txt` or `pfcrt_50.txt`, so the shipped benchmark harness does not exercise them.

`SRR8209082` is also worth treating with care as a performance baseline. Its 81,534 unitigs are
roughly 37x below the median of a 50-accession ERG11 run (2,986,072) measured during the
optimization work (that cohort list is not committed), and its 5-hop subgraph pulls in 26,323
nodes, about a third of its entire graph. It is fast and tangled, which
is not the common case.

## Regression coverage (verified)

The verifier was checked against two deliberate regressions:

1. A fixture with one minor allele removed. Caught, reporting `n_bubbles,
   quasispecies_branches`.
2. The `int`/`str` key migration hazard, injected by making `local_nodes` integer-keyed while
   `seqs` stays string-keyed. This empties the GFA and GraphAligner aborts with SIGABRT. The
   capture exits non-zero and the verifier reports a failure.

Hazard 2 is the one to keep in mind when bounding the subgraph index. In the full walker the
exception was caught per accession in the original Python walker and recorded as
`status: "error: ..."`, so the process still exits 0 and simply reconstructs nothing. With no
log redirect in `run_benchmark.sh`, that is easy to mistake for a working run.
