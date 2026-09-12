
# Path walker benchmarks and measurements

| script | purpose |
|---|---|
| `capture_golden.py` | capture a golden-truth fixture for one accession |
| `verify_golden.py` | check a walker build against `golden/` (exit 0 = concordant) |
| `measure_fork_loss.py` | measure the haplotypes discarded at fork choice |
| `run_benchmark.sh` | timing runs over the accession lists in this directory |

`golden/README.md` documents the fixtures, including which fields are stable enough to assert.

---

## Fork-loss measurement

`measure_fork_loss.py` reproduces every figure in
`work/issue_intra_accession_variance_lost.md`: the haplotypes discarded when the walk takes
the reference-matching branch at a fork instead of the abundant one.

```bash
python3 work/benchmark/measure_fork_loss.py --all
python3 work/benchmark/measure_fork_loss.py --all --json
python3 work/benchmark/measure_fork_loss.py -a SRR8209082 -q work/viral/query_hiv1_pol.fasta
```

Needs `blastn`, `makeblastdb`, `GraphAligner`, `curl` and `zstd`. Binaries resolve the same way
`logan_path_walker.find_bin` resolves them, with `.condaenv-walker/bin` added as a fallback.
Decompressed unitigs are cached under `--workdir` (default `/tmp/measure_fork_loss`), so the
first run of an accession streams from S3 and later runs take seconds.

### Why it is a copy

`logan_path_walker.py` deletes its GFA and GAF in a `finally` block (`:327-330`), so the
evidence is destroyed on every run. Stages 1 to 5 here are copied from the walker and keep
those files. Phasing is reimplemented with the orientation fix from
`rust-walker/src/stages/phase.rs`, which the Python walker does not have.

Because it is a copy, it can drift:

```bash
python3 work/benchmark/measure_fork_loss.py --all --self-check   # exit 1 on drift
```

`--self-check` fails unless the copied stages still reproduce the golden graph shape
(`n_unitigs`, `n_anchors`, `subgraph_nodes` from `golden/`) and the orientation check still
passes on every path step.

### What is asserted, and what is not

The orientation check is a precondition, not a result. Reading `L:s1:t:s2` in `u`'s header as
the traversal `(u,s1) -> (t,s2)` must make every consecutive GAF path step a real oriented
edge; the flipped-target reading must match none. Measured: 1278/1278 and 0/1278 on
`SRR8209082`, 114/114 and 0/114 on `ERR2028820`. Every claim about which nodes lie on the
walked path depends on this, so it is checked before anything else is reported.

All reported figures are byte-stable across `PYTHONHASHSEED` 1 and 7 and across repeat runs.
Branch and neighbour iteration is sorted for this reason: unsorted, the BFS walks Python sets
of string tuples and picks a different rejoin point per process, which gave 48, 50 and 51
clean alternatives on three runs of the same input.

Two things are deliberately **not** stable and are not asserted:

- **GraphAligner's row count.** `SRR8209082` yields one or two rows across runs (measured 1,
  1, 1, 2, 1). The best row by residue matches was byte-identical in all five, so the primary
  haplotype is stable and the count is reported rather than asserted.
- **The alternative haplotypes themselves**, which no fixture pins. They are constructed by
  splicing a detour back into the walked path, then scored with `blastn` against the reference,
  so they are validated on each run rather than compared against a stored expectation.
