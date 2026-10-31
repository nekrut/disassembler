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
| `--sieve-min-len <N>` | `--min-len` | 95% of the shortest reference record | minimum sequence length. The 5% margin is there because reconstructions are often a few bases shorter than the reference; the shortest record, because in multi-anchor mode each primary is reconstructed against its own reference |
| `--sieve-max-gap-frac <F>` (`from-lexicmap` only) | `--max-gap-frac` | `0.05`; derived in `from-lexicmap` (below) | the largest share of a sequence that may be gaps (`-`) or other non-ACGT characters |
| `--sieve-bin <PATH>` | `--sieve-bin` | discovered | path to `l1_pipeline_serial` (from `PATH`, or a build in `extern/kmer-diversity-sieve`) |

The atomic `sieve` takes the reference through `--ref`, used only to derive `--min-len` from its
shortest record. The input and `--ref` can be gzip-compressed; the wrapper decompresses the input
into the deduplicated copy it hands the binary, which reads only plain FASTA. `reconstruct` and `from-kmindex` derive it from every `-q` record.

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
row. To run the atomic `sieve` on a LexicMap cohort, set `--max-gap-frac` the same way. Either
way, the log says how many sequences the gate dropped, and for which reason.

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
