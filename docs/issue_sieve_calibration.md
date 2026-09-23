# The k-mer diversity sieve is calibrated for archives, not for gene cohorts

**Status:** characterised, defaults unchanged. Investigated 2026-09-21.

`--sieve` collapses a nine-accession Zika NS5 cohort to a single haplotype at its default
threshold, discarding the only lineage in the set that carries a non-negligible branch length.
This is not a defect in the sieve. It is operating exactly as specified, and its specification
targets a regime roughly ten times more divergent than an intra-species gene cohort.

Two separate things go wrong, and only one of them is reachable from `disassemble.py`.

---

## 1. The short version

| | |
|---|---|
| Symptom | a cohort of distinct haplotypes collapses to one at `--sieve-threshold 0.15` |
| Cause 1 | `tau = 0.15` corresponds to **0.5% nucleotide divergence**; the cohort's median pairwise divergence is 0.18% |
| Cause 2 | the default 64-position MinHash sketch cannot resolve novelty finer than `1/64`, i.e. ~0.05% divergence, whatever `tau` is set to |
| Governing relation | **`d ~ tau / (2k)`**, independent of gene length |
| Working setting here | `--threshold 0.005 --sketch 256`, which is lossless on the test cohort |
| Blocker | `disassemble.py` does not expose `--sketch`, so only half the fix is reachable |
| Not a defect in | the sieve's correctness, its cost, or its behaviour on the archive-scale data it was built for |

---

## 2. How the sieve decides

`is_novel` (`scripts/l1_pipeline.c`) promotes a sequence if, for every already-promoted sequence
it collides with under LSH, the estimated Jaccard similarity is at or below `1 - tau`:

```c
double need = 1.0 - threshold;
for (int band = 0; band < B; band++) {
    bucket_t *b = bucket_find(p, band_key(mbucket, band, sk, R));
    if (!b) continue;
    for (int t = 0; t < b->n; t++)
        if (jaccard(sk, &p->sk[b->ids[t] * p->m], p->m) > need) return 0;
}
return 1;
```

Jaccard is estimated by counting matching positions in an `m`-position MinHash sketch:

```c
static double jaccard(const uint64_t *a, const uint64_t *b, int m) {
    int e = 0; for (int i = 0; i < m; i++) if (a[i] == b[i]) e++;
    return (double)e / m;
}
```

So `tau` is a **novelty floor**: promote iff `1 - maxJaccard >= tau`. Higher `tau` collapses more.

Defaults, from `l1_pipeline.c`:

| constant | value |
|---|---|
| `DEF_K` | 15 |
| `DEF_M` (sketch) | 64 |
| `DEF_THRESHOLD` | 0.15 |
| `DEF_WARMUP` | 2000 |
| `LSH_BANDS` | 16, so `R = m / 16` rows per band |
| default `--min-len` | 27000, i.e. a whole coronavirus genome |

---

## 3. What `tau` means in nucleotides

An isolated SNV destroys up to `k` k-mers in each sequence. For a sequence of `K` k-mers with `s`
substitutions:

```
J = (K - k*s) / (K + k*s)          novelty = 1 - J  ~  2*k*d
```

where `d = s / L` is fractional divergence. Inverting:

> **`d ~ tau / (2k)`**

At `k = 15` and `tau = 0.15` that is `d = 0.005`, **0.5% nucleotide divergence**.

`L` cancels, so the threshold is **independent of gene length**. One `tau` does work across genes;
it is simply pointed at the wrong divergence.

---

## 4. Measured: the boundary is correctly centred but fuzzy

Reference plus one mutant at a controlled number of SNVs, 20 independent mutants per cell,
reporting the fraction of trials in which the two were kept apart. `tau = 0.15`.

```
 sketch |     5    10    15    20    25    30    40    60   SNVs
        |  0.2%  0.4%  0.6%  0.7%  0.9%  1.1%  1.5%  2.2%   divergence
     64 |  0.00  0.00  0.45  0.65  1.00  0.95  1.00  1.00
    256 |  0.00  0.00  0.55  1.00  1.00  1.00  1.00  1.00
```

The 50% point sits at 15 SNVs in 2,709 nt, where `novelty = 0.154` crosses `tau = 0.15`. The
estimator is therefore unbiased and centred exactly where section 3 predicts.

What is wrong is the **width**. At `sketch = 64` the decision is close to a coin flip across a
two-fold range of divergence, because 64 samples give a standard error near 0.045 on the Jaccard
estimate. Near the boundary the outcome depends on which k-mers happened to hash where, not on the
biology. Raising the sketch narrows the zone.

Length independence, same settings, fraction kept apart:

```
 gene nt |   0.2   0.3   0.4   0.5   0.6   0.8   1.0   1.5   % divergence
     900 |  0.00  0.17  0.42  0.42  0.50  1.00  1.00  1.00
    2709 |  0.00  0.00  0.08  0.17  0.75  0.92  1.00  1.00
    6000 |  0.00  0.00  0.17  0.42  0.67  1.00  1.00  1.00
```

All three transition in the same 0.3 to 0.8% band, as the algebra requires.

---

## 5. The test cohort is less diverse than it looks

Nine reconstructed Zika NS5 primaries (`SRR6251963`, `6251966`, `6251967`, `6251969`, `6251973`,
`6251979`, `6251984`, `6251992`, `6251993`), 2,708 nt:

```
9 haplotypes -> 4 UNIQUE sequences
  4 accessions   identical        SRR6251973 / 79 / 92 / 93
  3 accessions   identical        SRR6251963 / 66 / 69
  1 accession     1 SNV apart     SRR6251967
  1 accession    16 SNVs apart    SRR6251984
```

```
pairwise SNVs  min 0, median 5, max 16
pairwise %div  min 0.00%, median 0.18%, max 0.59%
pairs at or above the 0.5% boundary   4 of 36
```

So "9 collapsed to 1" overstates the loss. **The lossless answer is 4, not 9.** Only one pair
separation in the cohort reaches the sieve's boundary at all, and it lands in the fuzzy part of it.

---

## 6. Which lineage dies, and at what threshold

```
tau = 0.005  ->  SRR6251963, SRR6251973, SRR6251984      3 of 4
tau = 0.02   ->  SRR6251963, SRR6251973, SRR6251984      3 of 4
tau = 0.04   ->  SRR6251963, SRR6251984                  2 of 4
tau = 0.08   ->  SRR6251963, SRR6251984                  2 of 4
tau = 0.15   ->  SRR6251963                              1 of 4   <- default
```

At the default the survivor is the first record in file order, and the loss is **`SRR6251984`**,
the 16-SNV lineage at 0.59% divergence. On a FastTree GTR tree of this cohort it is the only
haplotype tip with a non-negligible branch length (0.0041, against approximately 0 for every
other tip), so the default discards precisely the sequence the phylogeny depends on.

Note that lowering `tau` alone never reaches 4 of 4. That is the second cause.

---

## 7. The sketch sets a resolution floor that `tau` cannot cross

Jaccard is estimated from `m` samples, so novelty cannot be resolved more finely than `1/m`. At
`m = 64` that floor is 0.0156 novelty, about 0.05% divergence by section 3. The 1-SNV separation
in this cohort is 0.037% divergence, which is **below the floor**, so no value of `tau` recovers
it.

Raising the sketch, at `tau = 0.005`:

```
sketch   64  ->  3 of 4 unique    (floor 0.0156 novelty)
sketch  256  ->  4 of 4           (floor 0.0039)
sketch 1024  ->  4 of 4           (floor 0.0010)
sketch 4096  ->  4 of 4
```

`--threshold 0.005 --sketch 256` is **lossless on this cohort**: it recovers all four unique
sequences, including both the 1-SNV and the 16-SNV lineages.

Varying `k` from 9 to 31 changes nothing at that setting, so `k` is not the lever.

### The sketch is not a free dial

`R = m / LSH_BANDS`, so a larger sketch makes each LSH band wider, band collisions rarer, and
near-duplicates stop being compared at all. On a synthetic 1,000-sequence cohort built from 8
lineages:

```
sketch   64  ->  kept 8     correct
sketch  256  ->  kept 8     correct
sketch 1024  ->  kept 41    leaking near-duplicates
```

256 has headroom on both sides. Larger is not better.

---

## 8. The defaults are right for the job the sieve was built for

From the sieve's own README, validated on a SARS-CoV-2 consensus archive:

| Metric | Value |
|---|---|
| Input | 489 GiB, ~17.6M genomes |
| Promoted | 1,443 |
| Overall reduction | 12,209x |
| Wall time | ~20 min on 128 MPI ranks |

That is the design target: an archive in which most records are near-exact duplicates, reduced
four orders of magnitude to roughly one representative per broad lineage. The default `--min-len`
of 27000 says the same thing. Keeping one genome per lineage out of millions is a different
problem from keeping every distinct haplotype in a cohort of nine.

The cost is negligible and the behaviour is correct when the cohort genuinely has lineage
structure. Synthetic cohort, 8 lineages separated by more than 1% divergence, `tau = 0.15`:

```
     N   sieve s   kept   MSA+tree full   MSA+tree sieved       net
   100      0.05      8            2.72            1.08     +1.59
   400      0.11      8           24.78            1.08    +23.60
  1000      0.24      8           26.39            1.05    +25.10
```

It recovers exactly the 8 real lineages, costs a fraction of a second, and saves 25 s at N = 1000.
A larger sketch stays cheap: 0.39 s at `--sketch 256` on the same 1,000 sequences.

So the sieve is worth keeping. It is the gene-scale regime that sits outside its range.

---

## 9. Only half the fix is reachable from disassemble.py

`python/disassemble.py` builds:

```python
scmd = [sieve_bin, "--input", str(recon),
        "--min-len", str(sieve_min_len), "--max-n-frac", "0.05",
        "--threshold", str(args.sieve_threshold), "--warmup", "2000",
        "--emit-promoted", str(emit_prefix)]
```

The sieve accepts eleven flags:

```
--emit-promoted --input --k --limit --max-n-frac --min-len
--outlier-sd --sketch --threads --threshold --warmup
```

Four are passed, two are hardcoded (`--max-n-frac 0.05`, `--warmup 2000`), and five are
unreachable: **`--sketch`**, `--k`, `--threads`, `--limit`, `--outlier-sd`.

The consequence is that `--sieve-threshold` is the only available knob, and on its own it cannot
fix a gene-scale cohort: lowering `tau` past the sketch's resolution floor moves the boundary
somewhere the estimator cannot see. That is why the sieve currently reads as unusable here rather
than merely mis-defaulted.

`--threads` being unreachable also means the serial binary is always what runs, which is
irrelevant at these sizes but worth knowing.

---

## 10. Options

1. **Expose `--sieve-sketch`** in `disassemble.py`, defaulting to 64 so nothing changes unless
   asked, and put `d ~ tau / (2k)` in the `--sieve-threshold` help so the threshold is choosable
   rather than guessable. This is the smallest change that makes the sieve usable at gene scale.
2. **Leave the sieve default-off**, as it is, and document that it targets archive-scale
   redundancy. The `--sieve` help already warns it is a lineage-level reducer.
3. **Derive `tau` from the cohort** rather than accepting a constant: measure pairwise divergence
   on a sample of the reconstructions and set `tau ~ 2*k*d_min`. More machinery than the problem
   currently justifies.

Nothing here has been changed. The sieve remains default-off, and its defaults are untouched.

---

## 11. Reproducing this investigation

Requires the submodule built (`git submodule update --init && (cd extern/kmer-diversity-sieve &&
make serial)`).

```bash
SV=extern/kmer-diversity-sieve/l1_pipeline_serial

# section 6: which lineage dies, on a real reconstructed cohort
for tau in 0.005 0.02 0.04 0.08 0.15; do
  rm -f /tmp/z.*
  $SV --input out/zika/reconstructed.fasta --min-len 2500 --max-n-frac 0.05 \
      --threshold $tau --warmup 2000 --emit-promoted /tmp/z >/dev/null 2>&1
  echo "tau=$tau -> $(cat /tmp/z.rank*.fasta | grep -c '^>')"
done

# section 7: the resolution floor
for m in 64 256 1024; do
  rm -f /tmp/z.*
  $SV --input out/zika/reconstructed.fasta --min-len 2500 --max-n-frac 0.05 \
      --threshold 0.005 --sketch $m --warmup 2000 --emit-promoted /tmp/z >/dev/null 2>&1
  echo "sketch=$m -> $(cat /tmp/z.rank*.fasta | grep -c '^>') of 4 unique"
done
```

Sections 4 and 8 use synthetic cohorts: a random 2,709 nt reference, mutants at controlled SNV
counts, and for section 8 a cohort built from 8 lineages with within-lineage noise. The nine Zika
accessions in section 5 reconstruct in about 30 s at `-w 8` against
`queries/markers/zika_ns5.fasta`.
