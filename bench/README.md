# Benchmarks and concordance

| path | purpose |
|---|---|
| `check_golden.sh` | golden gate: runs `logan-walker` under old-behaviour flags and checks the primary byte for byte against `golden/` |
| `golden/` | golden-truth fixtures (see `golden/README.md`) |
| `simulated_kmer_sensitivity/` | the 121,500-trial simulation behind the anchor default `k = 21` |
| `*_50.txt`, `smoke.txt` | accession lists |

## Golden concordance gate

```bash
bench/check_golden.sh          # exit 0 = concordant, 1 = regression
```

Byte-exact on two clonal, deterministic accessions, `ERR2028820` (PfCRT) and `SRR12363134`
(*C. auris* ERG11). `SRR8209082` (HIV-1 *pol*) is reported but not gated: its subgraph is a
large share of the whole graph, and its unpruned route is not reproducible across GraphAligner
versions. The gate needs GraphAligner 1.0.17b on `PATH` and network access to the Logan S3
bucket, so it is the manual `golden` CI job, not the every-push gate. The header of the script
documents why each accession is or isn't gated.

## The reference-bias finding

Why the walk can take a reference-matching minority allele, how that was measured, and what the
walker's pruning and majority-preference defaults change, is written up in
[`../docs/issue_intra_accession_variance_lost.md`](../docs/issue_intra_accession_variance_lost.md).
