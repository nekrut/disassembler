# Benchmarks and concordance

| path | purpose |
|---|---|
| `check_golden.sh` | Rust-native golden gate: runs `logan-walker` under old-behaviour flags and checks the primary byte-for-byte against `golden/` |
| `golden/` | golden-truth fixtures (see `golden/README.md`) |
| `simulated_kmer_sensitivity/` | the 121,500-trial simulation establishing k=21 |
| `*_50.txt`, `smoke.txt` | accession lists |
| `archive/` | retired tooling that targeted the removed Python walker (see `archive/README.md`) |

## Golden concordance gate

```bash
bench/check_golden.sh          # exit 0 = concordant, 1 = regression
```

Byte-exact on two clonal, deterministic accessions — `ERR2028820` (PfCRT) and `SRR12363134`
(*C. auris* ERG11) — whose primaries reproduce byte for byte. `SRR8209082` (HIV-1 *pol*) is
reported but not gated: its subgraph is a third of the whole graph and GraphAligner is not
deterministic on it when unpruned. Needs `GraphAligner` on `PATH` and network access to the
Logan S3 bucket, so it is the manual `golden` CI job, not the every-push gate. The header of the
script documents why each accession is or isn't gated.

## The reference-bias finding

The measurement of the haplotypes discarded at fork choice — why the walk can take a
reference-matching minority allele — is written up in full in
[`../docs/issue_intra_accession_variance_lost.md`](../docs/issue_intra_accession_variance_lost.md).
The script that produced those figures, `measure_fork_loss.py`, is under `archive/` because it
was a copy of the removed Python walker's stages.
