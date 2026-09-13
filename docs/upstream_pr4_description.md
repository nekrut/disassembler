> **Archival: verbatim text of upstream `nekrut/logan-hypyaeon` PR #4.** The body below is left
> unedited as a record of what was proposed, so it names flags and paths that no longer exist here.
> Read none of it as current behaviour; for that see [`../README.md`](../README.md).
>
> | named below | status in this repo |
> |---|---|
> | `--alt-engine detour\|variants\|off` | Removed. The `detour` engine is retired and `stages::variants` is the sole intra-host alt engine, so there is no longer an engine to select. |
> | `-H many` | Removed, along with `--haplotypes`. Alternatives are now controlled by `--emit-alts` (on by default) and written to `alt_haplotypes.fasta`, **not** into `reconstructed.fasta`, so they no longer enter the codon MSA as extra taxa. That change is the fix for the problem this PR text describes. |
> | `work/issue_intra_accession_variance_lost.md` | Now [`issue_intra_accession_variance_lost.md`](issue_intra_accession_variance_lost.md), in this directory. |
> | `work/benchmark/measure_fork_loss.py` | Now `bench/archive/measure_fork_loss.py`, and no longer runnable: it depends on the Python reference walker, which has been removed. The figures stand as measured; they cannot be reproduced from that script. |
>
> The GraphAligner threading result in section 1 still holds, and `-t 1` is still the default.

From Claude:

This adds intra-accession haplotype recovery and the fix for what was blocking it. Your detour stage is kept intact and selectable alongside it (--alt-engine detour|variants|off), so both can be compared on the same data.

  Two things worth a look before the rest.

  1. GraphAligner is not reproducible above one thread. Same graph, same query, five invocations, best-alignment score:

  -t 1   2691 2691 2691 2691 2691    deterministic, and the best score
  -t 2   2686 2691 2691 2691 2691
  -t 8   2686 2686 2686 2686 2686    stable, but consistently worse

  Not a tie-break: a losing run fails to find the better alignment. This reaches the output, five identical walker runs produced two different primary haplotypes. Default is now -t 1, which costs 0.47s against 0.41s on the largest graph here, and throughput comes from compute_workers regardless.

  2. What -H many emits, measured across the full 200-accession panel. With pruning on, both stages find the same number of sites but classify them differently:

                                       detour         variants
  alternatives                         269            267
    genuinely two alleles (VAF < 0.5)  131 (48.7%)    267 (100%)
    route took a <10% branch            98 (36.4%)      0
  adds a gap its own primary lacked    114 (42.4%)      6 (2.2%)

  The records above 0.5 VAF are not co-existing alleles. They are sites where GraphAligner took the rarer branch because it resembled the reference, so the majority allele is what gets emitted as the "alternative". Since -H many writes these into reconstructed.fasta, each becomes an independent taxon in the codon MSA.

  The mean_ka ranking cannot separate the two cases: on SRR8209082 the gap between _major and the runner-up was 0.036%, while swapping a single unitig out of ~1,300 moves that mean by 0.08%. The comparison that does work is between the two branches at the fork, which is 8,000x against 3,600x.

  Flagging it because dengue_4_ns5 quotes 2,023% yield with -H many, and that ratio shapes the BUSTED and MEME input.

  Root cause, fix and panel validation are in work/issue_intra_accession_variance_lost.md. Every figure above reproduces via work/benchmark/measure_fork_loss.py.
