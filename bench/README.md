# Benchmarks and concordance

| path | purpose |
|---|---|
| `check_golden.sh` | the golden gate: runs `logan-walker` twice and checks its outputs against `golden/` |
| `mask_summary.py` | replaces the volatile fields of a summary JSON so two runs can be diffed |
| `golden/` | the golden-truth fixtures (see [`golden/README.md`](golden/README.md)) |
| `simulated_kmer_sensitivity/` | the 121,500-trial simulation behind the anchor default `k = 21` |
| `*_50.txt`, `smoke.txt` | accession lists |

## The golden gate

```bash
bench/check_golden.sh          # exit 0 = concordant, 1 = regression
```

It runs two sets over the same three accessions, `ERR2028820` (PfCRT), `SRR12363134`
(*C. auris* ERG11) and `SRR8209082` (HIV-1 *pol*):

| set | flags | what is gated |
|---|---|---|
| old-behaviour | `--prune-frac 0 --no-prefer-majority --no-emit-alts` | the primary haplotype of `ERR2028820` and `SRR12363134`, byte for byte. `SRR8209082` is reported but not gated |
| default-flag | none | every output of all three accessions, byte for byte, including the summaries |

`SRR8209082` is ungated in the first set only. Its 34,582-node neighbourhood is 42% of the whole
graph, against 0.06% and 0.02% for the other two, so its **unpruned** route depends on the
GraphAligner version and is not reproducible across them. The default pruning is what settles
it, which is why the second set can gate it. The header of the script records why each
accession is or is not gated in each set.

Summaries are compared through `mask_summary.py`, which replaces the fields that cannot
reproduce (timings, paths, the tool version, the fetch metadata) and leaves everything else to
be gated, `payload.sha256` included. The list of masked fields is in
[`../docs/summaries.md`](../docs/summaries.md).

The gate needs GraphAligner 1.0.17b on `PATH` and network access to the Logan S3 bucket, so it
runs in its own workflow (`.github/workflows/golden.yml`) rather than in the every-push
`ci.yml`: on pushes to `main`, on `v*` and `logan-walker-v*` tags, on pull requests touching
`walker/`, weekly, and on manual dispatch. Set `LOGAN_UNITIGS_CACHE` to keep the downloaded
objects between runs, and `LOGAN_GOLDEN_WORK` to keep the outputs a run produced.

## The reference-bias finding

Why the walk can take a reference-matching minority allele, how that was measured, and what the
walker's pruning and majority-preference defaults change, is written up in
[`../docs/notes/intra_accession_variance.md`](../docs/notes/intra_accession_variance.md).
