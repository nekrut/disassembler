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
[`vcf_example.md`](vcf_example.md)).

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
| `--sieve-min-len <N>` | `--min-len` | 95% of the reference length | minimum sequence length. The 5% margin is there because reconstructions are often a few bases shorter than the reference |
| `--sieve-bin <PATH>` | `--sieve-bin` | discovered | path to `l1_pipeline_serial` (from `PATH`, or a build in `extern/kmer-diversity-sieve`) |

The atomic `sieve` takes the reference through `--ref`, used only to derive `--min-len` from its
first record. The wrapper always passes the sieve binary `--max-n-frac 0.05`, so a sequence
with more than 5% non-ACGT characters is dropped as low quality. In `from-lexicmap`, byte-identical
duplicates are removed (keeping the first) before the sieve runs.

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

```
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
is in [`issue_sieve_calibration.md`](issue_sieve_calibration.md).
