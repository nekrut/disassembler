# Golden-truth fixtures

Concordance baselines for the path walker, covering success criteria 2 and 3 of the
optimization issue (`docs/issue_path_walker_optimization.md`).

| Accession | Gene | Query | Unitigs | Anchors | Subgraph | Primary | Bubbles |
|---|---|---|---|---|---|---|---|
| `ERR2028820` | PfCRT | `cds/pfcrt_only.fasta` | 1,277,362 | 76 | 356 | 1,274 nt | 202 |
| `SRR8209082` | HIV-1 pol | `viral/query_hiv1_pol.fasta` | 81,534 | 4,102 | 26,323 | 2,738 nt | 11,825 |

Captured on the unmodified walker at commit `81fc1fd`, GraphAligner 1.0.20, BLAST 2.17.0+.

## Usage

```bash
python3 bench/verify_golden.py     # exit 0 = concordant, 1 = regression
```

Regenerate after an intentional, reviewed behaviour change:

```bash
python3 bench/capture_golden.py --accession ERR2028820 \
    --query queries/pfcrt_only.fasta --workdir /tmp/cap
```

Both scripts need `blastn`, `makeblastdb`, and `GraphAligner` on `PATH`.

## What is asserted, and what is not

Every asserted field was verified byte-stable across two `PYTHONHASHSEED` values (1 and 7):

- graph shape (`n_unitigs`, `n_anchors`, `subgraph_nodes`)
- the primary haplotype, by SHA-256 of the sequence, plus `cov` and `pident`
- `n_bubbles`, and every bubble's fork node with its full branch list (node, `ka:f`, VAF)

The bubble `merge` node is deliberately **not** asserted. `phase_superbubbles` picks it with
`reconverged[0]` (`python/logan_path_walker.py:259`) from a `Counter` populated by iterating a `set`,
so it depends on string hash randomisation. Measured across seeds 1 and 7, it differed on
10,922 of 11,825 fork/merge pairs for `SRR8209082` and 18 of 202 for `ERR2028820`. The branch
calls themselves, which carry the biology, were identical in every case.

Timings are excluded. They are not reproducible and belong in benchmark output, not a fixture.

## Why these two accessions

They are the two named in the issue's concordance criterion. Note that neither appears in
`hiv1_50.txt` or `pfcrt_50.txt`, so the shipped benchmark harness does not exercise them.

`SRR8209082` is also worth treating with care as a performance baseline. Its 81,534 unitigs are
roughly 37x below the median of the committed 50-accession ERG11 run (2,986,072), and its 5-hop
subgraph pulls in 26,323 nodes, about a third of its entire graph. It is fast and tangled, which
is not the common case.

## Regression coverage (verified)

The verifier was checked against two deliberate regressions:

1. A fixture with one minor allele removed. Caught, reporting `n_bubbles,
   quasispecies_branches`.
2. The `int`/`str` key migration hazard, injected by making `local_nodes` integer-keyed while
   `seqs` stays string-keyed. This empties the GFA and GraphAligner aborts with SIGABRT. The
   capture exits non-zero and the verifier reports a failure.

Hazard 2 is the one to keep in mind when bounding the subgraph index. In the full walker the
exception is caught per accession at `python/logan_path_walker.py:392` and recorded as
`status: "error: ..."`, so the process still exits 0 and simply reconstructs nothing. With no
log redirect in `run_benchmark.sh`, that is easy to mistake for a working run.
