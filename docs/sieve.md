# Diversity thinning with the k-mer sieve

`disassembler sieve`, and `--sieve` in the composites, thin a cohort of near-duplicate
haplotypes before the codon MSA and tree. The MSA and tree are `O(N^2)`, so on a large cohort of
nearly identical lineages, thinning first can be a large speedup with little loss of
phylogenetic signal. The sieve is **off by default**; without `--sieve`, the composites behave
exactly as if it did not exist.

The sieve is the `l1_pipeline_serial` binary from
[kmer-diversity-sieve](https://github.com/stevenweaver/kmer-diversity-sieve) (vendored as the
`extern/kmer-diversity-sieve` submodule, and packaged on Bioconda as `kmer-diversity-sieve`).
It reads sequences in order and promotes one only if it is novel enough against everything
already promoted, using MinHash sketches and LSH.

## What it is not

The sieve is a **lineage-level** reducer. What decides whether a sequence survives is its
absolute divergence from an already-promoted representative, not whether it came from a
different accession. It is **not** an intra-host SNV tool: sequences that differ by a base or
two are below the resolution any usable sketch can distinguish. To keep intra-host single-SNV
variation, use the walker's `--emit-vcf` / `--emit-snv-table` instead (see
[`examples/vcf_example.md`](examples/vcf_example.md)).

## Running it

```bash
# As a stage of a composite
disassembler reconstruct -l accessions.txt -q queries/markers/zika_ns5.fasta -o out/zika --sieve

# On its own
disassembler sieve -i out/zika/reconstructed.fasta -o out/zika/sieved.fasta \
    --ref queries/markers/zika_ns5.fasta
```

| composite flag | atomic flag | default | meaning |
|---|---|---|---|
| `--sieve` | | off | enable the stage |
| `--sieve-threshold <F>` | `--threshold` | `0.005` | novelty threshold tau. A sequence is promoted only if `1 - maxJaccard >= tau` against the promoted set, so a **higher** tau collapses **more** |
| `--sieve-sketch <M>` | `--sketch` | `256` | MinHash sketch size. Novelty finer than `1/M` cannot be resolved |
| `--sieve-min-len <N>` | `--min-len` | a fraction of the shortest reference record (below) | minimum sequence length |
| `--sieve-max-gap-frac <F>` (`from-lexicmap` only) | `--max-gap-frac` | `0.05`; derived in `from-lexicmap` (below) | the largest share of a sequence that may be gaps (`-`) or other non-ACGT characters |
| `--sieve-bin <PATH>` | `--sieve-bin` | discovered | path to `l1_pipeline_serial` (from `PATH`, or a build in `extern/kmer-diversity-sieve`) |

The reference is used only to derive `--min-len`, and which fraction of it depends on what the
cohort came from:

| route | default `--min-len` |
|---|---|
| `reconstruct`, `from-kmindex` | `min(0.85, --min-cov)` of the shortest `-q` record, following the walker's own coverage gate |
| `disassembler sieve`, `from-lexicmap` | 95% of the shortest `--ref` record |

The composites follow the walker because anything the walker published is a reconstruction this
step is not meant to re-judge: it thins near-duplicates, and a flat 95% drops every haplotype
between the walker's gate and 95% although the walker and the MSA both accept it, which is
exactly the walk that stops short of a terminus or the sample with a real in-frame deletion. It
is the shortest record because in multi-anchor mode each primary is reconstructed against its
own reference, and there is no margin on top because `cov` is measured against the winning
anchor, which is at least that long. The atomic `sieve` and `from-lexicmap` have no walker gate
to follow, so pass `--min-len` explicitly on a cohort the walker produced.

The input and `--ref` can be gzip-compressed; the wrapper decompresses the input into the
deduplicated copy it hands the binary, which reads only plain FASTA.

### Before the sieve runs

**Duplicates are removed.** Every route removes byte-identical sequences first, keeping the
first copy, and logs how many it removed. A later copy could never be promoted, but it would
still feed the binary's admission check, which rejects a sequence as an error when it has far
more never-seen k-mers than the running mean. Each duplicate adds a zero to that mean, so once
the binary's 2,000-record warmup is over, a genuinely new lineage looks like an outlier. On four
C. auris cohorts, the raw input lost 26 to 273 sequences this way; the deduplicated input lost
none.

**The completeness gate.** The binary drops a sequence shorter than `--min-len`, or one in
which more than `--max-gap-frac` of the characters are not A, C, G or T. Gaps (`-`) count.
Walker haplotypes have no gaps, so for them this is a check on N. LexicMap rows are on the
reference grid, with `-` wherever no hit reached, so the same limit acts as a second coverage
gate: at 5% it dropped 22% of a default C. auris FKS1 cohort, and a third of the distinct
haplotypes in a SARS-CoV-2 spike cohort, all rows lexicmap had already accepted.
`from-lexicmap` therefore derives the limit from the cohort's own coverage gate, as
`1 - gate + 0.05`:

| cohort fed to the sieve | coverage gate | default `--sieve-max-gap-frac` |
|---|---|---|
| default | `--min-cov` (`0.80`) | `0.25` |
| `--require-full-length` | `1.0` | `0.05` |
| `--include-partials` | `--min-cov-partial` (`0.50`) | `0.55` |

The extra 0.05 leaves room for in-frame deletions, which count as covered but are `-` in the
row, **and the derived value is never set below the cohort's own worst row**. lexicmap counts a
trusted deletion as covered with no limit, so a clean, full-length row whose deletion is longer
than 5% of the gene (45 nt in 900) used to be dropped here, by a step documented as
near-duplicate thinning, with exit 0 and a `6 -> 5 haplotypes` line. Indel alleles are often
what a cohort is for. The run prints the worst row's share beside the gate it derived. To run
the atomic `sieve` on a LexicMap cohort, set `--max-gap-frac` the same way. Either way, the log
says how many sequences the gate dropped, and for which reason.

Keeping partial rows has a cost. The sieve compares k-mer sets, so a partial row looks novel
next to a full-length one it otherwise matches, and is rarely thinned. On the FKS1 cohort above,
the sieve promoted the same 391 full-length sequences at either limit, plus 1,284 partial rows
at 0.25; 154 of those matched a promoted full-length row wherever they had a base. To thin
harder, use `--require-full-length`, or pass a lower `--sieve-max-gap-frac`.

To build the binary from the submodule (a single C source, so only a C compiler is needed):

```bash
git submodule update --init && (cd extern/kmer-diversity-sieve && make serial)
```

## Choosing the threshold

The defaults above are **gene-scale** and deliberately differ from the sieve binary's own
defaults (`tau 0.15`, sketch `64`, `--min-len 27000`), which were calibrated for deduplicating
millions of near-identical SARS-CoV-2 genomes. At the binary's defaults, a gene-scale cohort of
one virus collapses to about one representative.

The divergence a threshold resolves is roughly

```text
d ~ tau / (2k)
```

at the sieve's k-mer size `k = 15`, independent of gene length. The binary's `tau 0.15` sits at
about 0.5% divergence; the wrapper's `0.005` sits near 0.017%. An intra-species gene cohort is
often well below 0.5% (the Zika NS5 test cohort has a median pairwise divergence of 0.18%).

The sketch size sets a separate floor that tau cannot cross: with `M = 64`, novelty below about
0.05% divergence is invisible at any tau, which is why the wrapper uses `256`. A bigger sketch
is not strictly better, because it also widens the LSH bands so that near-duplicates are
compared less often.

Measured at the wrapper defaults:

| cohort | reconstructions | unique | promoted | notes |
|---|---|---|---|---|
| Zika NS5 (36 accessions, PRJNA362728) | 36 | 12 | 9 | topology-identical (RF = 0) on the 9 tips, including the long-branch lineage the binary's defaults drop |
| HIV-1 pol (17 accessions near SRR8209082) | 17 | 17 | 17 | pairwise divergence is 2 to 3%, well above the boundary, so nothing collapses |
| 33 intra-host alternatives of one HIV-1 pol accession | 33 | 33 | 28 | single-SNV alternatives sit near the resolution floor |

**Rule of thumb:** look at your cohort's median pairwise divergence. Below about 0.5%, the sieve
can help; between 0.5% and 1%, it might; above 1%, it is either a no-op or a signal-losing
collapse. Sweep tau on your own data before trusting a setting, and expect the sieve to pay off
only when the cohort is large enough for the MSA and tree to dominate the run time.

The full calibration study, with the measured transition curves and the reproduction commands,
is in [`notes/sieve_calibration.md`](notes/sieve_calibration.md).

## Outputs

| file | contents |
|---|---|
| `<output>.fasta` | the promoted haplotypes |
| `<output>.summary.json` | the run, in the shared format of [`summaries.md`](summaries.md) |

`discards.records` carries `exact_duplicates`, `shorter_than_min_len`, `too_gappy` and
`not_promoted`, `cohort` carries the admission rate, and `records[]` names each promoted
haplotype by id and sequence hash. The summary is written on every exit path, including the
failures, and is new in 0.2.0: before it, this stage produced no machine-readable output at all.

## Limits and behaviour worth knowing

**The binary's admission gate is disarmed.** It has a second gate besides novelty: after a
warm-up it rejects a record whose count of never-seen k-mers is a statistical outlier, and a
rejected record is never sketched, promoted or added to the k-mer table, so the rest of its
clade goes with it. The wrapper does not arm it, because the cut-off it imposes is a property
of the file rather than of the biology: the band is `mean + 5*sd` over the running singleton
counts, and the first record alone contributes about fifty times what a later one does, so
there is no divergence a user can stay under. Measured on one 2,100-record cohort, a lineage 3%
divergent placed after record 2,000 was dropped and the same lineage placed first was kept.
Deduplication is therefore a speed-up and a memory saving, not a correctness fix.

**Headers are not shown to the binary.** It buckets novelty by a month parsed from `|YYYY-MM`
in the header, so twelve copies of one sequence dated in twelve different months all survived,
and it selects what to emit by the header's first 255 characters, so two records sharing a
header were both emitted when only one was promoted and a long header lost its tail. The
wrapper gives it `>r0`, `>r1`, ... and puts the original headers back afterwards.

**Two upstream limits the wrapper guards rather than fixes.** These are in
`kmer-diversity-sieve` itself (github.com/stevenweaver/kmer-diversity-sieve), not in this repo:

- the promoted set's LSH index is a fixed 2^20-slot table that is never grown, and lookup
  stops only at an empty slot, so past about 65,536 admitted records the binary can loop
  forever at 100% CPU. An input that would admit more than that is refused, with the count,
  rather than being handed over. Chunking and re-sieving is not a way round it, since the
  union can exceed the limit again.
- only the first ~65,550 k-mers of a record are sketched, so two records identical over their
  first 65.5 kb are judged near-duplicates whatever follows. A longer record earns a warning.

**`-march=native`.** The Makefile builds for the machine that built it, so a binary copied to
another CPU dies with SIGILL. That is now reported by name, with the rebuild command, rather
than as "exited -4".
