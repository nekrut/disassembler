# Path walking selects the reference-matching minority allele

**Status:** diagnosed 2026-09-11 and fixed. The fix has been on by default since 2026-09-12, when
the competing `detour` engine was also retired.

Why the walker has these three options, what bias they correct, and what they change on real data.
Every flag and default is in the CLI reference, [walker.md](walker.md).

| flag | default | effect |
|---|---|---|
| `--prune-frac <F>` | `0.01` | align once, drop unitigs below F of the gene's depth, align again |
| `--prefer-majority` / `--no-prefer-majority` | on | at a two-allele site, the primary takes the more abundant allele |
| `--emit-alts` / `--no-emit-alts` | on | each remaining second allele goes to `alt_haplotypes.fasta` |

Related: `--no-indel-alts`, `--max-alt-indel` (default 30), `--allow-frameshift-alts` (off),
`--alt-max-detour` (default 100), `--min-abund-frac` (default 0.0, off). The pre-fix behaviour is
`--prune-frac 0.0 --no-prefer-majority --no-emit-alts`, the flags
[bench/check_golden.sh](../bench/check_golden.sh) runs under.

## 1. The finding

The walker emits one primary sequence per accession, `<ACC>_primary`, by walking the single path
GraphAligner returns. GraphAligner scores a path by residue matches to the reference. At a fork
where the sample's dominant allele differs from the reference and a rare allele matches it, the
walk takes the rare allele. On `SRR8209082` (HIV-1 *pol*), before the fix, **49 of 50**
recoverable alternative haplotypes carried the majority allele that the walk discarded, with a
median of **778x** more coverage than the branch walked. None of these sites was reported as a
bubble.

## 2. How it was measured

The original prototype walker deleted its GFA and GAF after every run, so the measurement used a
replica of its stages 1 to 5 that kept them and reproduces both golden fixtures exactly
(`ERR2028820`: 356 subgraph nodes; `SRR8209082`: 26,323).

The orientation convention was checked against the path itself: reading `L:s1:t:s2` in unitig
`u`'s header as the traversal `(u,s1) -> (t,s2)` makes **1278 of 1278** consecutive GAF path steps
a real oriented edge, and flipping the target strand scores **0 of 1278**. Alternatives were built
by splicing each untaken branch back into the path and scored with `blastn` against the reference
CDS, the test the pipeline applies to its own output. Every figure is byte-stable across repeat
runs and two hash seeds (see [bench/golden/README.md](../bench/golden/README.md) on hash order).

## 3. Where the loss happens

The walked path for `SRR8209082` is 1,279 oriented nodes, 1,216 of them forks. At each fork the
branch taken was compared with those not taken, under the walker's own gates (`ka:f >= 5.0` and
VAF >= 10%, the `--min-abund` and `--min-vaf` defaults). That finds 58 credible untaken branches,
57 of which rejoin the path within 8 hops, and 50 clean single-allele swaps (coverage >= 0.85,
<= 50 nt). All 50 are full length, 2,738 nt against a 2,739 nt reference, with coverage 1.000 and
`blastn` identity 98.028% to 98.174% against the primary's 98.174%. They pass every gate the
primary passes. What separates them is coverage, in the wrong direction:

```
discarded branch better supported than the branch walked        49 / 50
discarded branch is the majority allele (VAF >= 0.5)            49 / 50
discarded branch carries >= 90% of local coverage               49 / 50
coverage ratio discarded : walked      median 778x, max 10,417x
walked branch < 10x while discarded branch >= 10,000x           16 cases
alternatives LESS identical to the reference than the primary   40 / 50
alternatives MORE identical                                      0 / 50
```

The path is not low-coverage overall: `ka:f` along the walked nodes has a median of 15,000, and
only 12 of 1,279 fall below the 5.0x floor. The walk runs almost entirely through the dominant
population and, at roughly 4% of its steps, drops into a near-absent allele because that allele
matches the reference. For selection analysis this matters: the divergence being suppressed is
the divergence dN/dS is computed from.

One case in full, path index 1236 near the 3' end of *pol*:

```
fork (54954,-), reached from (54872,-)
  branch (13182,+)  ka:f     18.9  VAF 0.002   len 35
  branch (54790,+)  ka:f  12000.0  VAF 0.997   len 31
  branch (54791,+)  ka:f     10.0  VAF 0.001   len 31
  branch (54792,+)  ka:f      2.0  VAF 0.0004  len 61   <- walked

pos 2639   primary T   alternative A    TAATAGTG[T/A]CATAAAAG
pos 2640   primary C   alternative A    AATAGTGT[C/A]ATAAAAGT
```

The discarded branch rejoins the path two nodes later, so the swap is clean and local. The emitted
primary carried `TC` at a site where 99.7% of the sample's k-mer coverage says `AA`.

**Why bubble phasing does not catch it.** A bubble needs two branches that clear both gates and
reconverge one step later. Of the 58 credible forks, 8 fail the VAF gate (the walked branch holds
well under 1%), 50 reconverge further out, and 0 are reported. Conversely, none of the 14 bubbles
found on `SRR8209082` has its fork on the walked path, since most of a 5-hop neighbourhood is not
the gene. For intra-host variants, read `alt_haplotypes.fasta` or the `--emit-vcf` output
([vcf_example.md](vcf_example.md)), not `bubbles.json`.

## 4. It depends on the sample

`ERR2028820`, PfCRT in *P. falciparum*, shows none of the effect. Its walked path is 115 nodes with
63 forks and a median `ka:f` of 180, against 1,279, 1,216 and 15,000 for `SRR8209082`. It has 4
credible untaken branches against 58, and 0 clean alternatives against 50: 3 of the 4 rejoin the
path and none survives, one rewriting 827 nt rather than swapping an allele and two falling to
0.669 and 0.704 coverage, below the 0.85 gate. The effect belongs to high-diversity samples.

## 5. The fix

### 5.1 `--prune-frac`: take the error branches away from the aligner

Pruning deletes low-coverage branches before GraphAligner sees them, which needs the gene's depth.
That number does not exist before an alignment and cannot be taken from the anchors:
`SRR8209082`'s anchors have a median `ka:f` of 70x against a true 14,500x, because
`--min-shared-kmers 3` admits thousands of shallow incidental matches. So the walker aligns once
to measure depth, drops unitigs below `prune_frac` of it, and aligns again. If the second
alignment fails, the first pass's haplotype stands.

```
SRR8209082 (HIV-1 pol, 14,500x)   nodes  aligns    cov   min ka:f  under 1% of depth
  unpruned                        34582     3    0.9996      2.0        118
  1% floor (145x)                  6930     1    0.9993    147.5          2
  5% floor (725x)                  5343     1    0.9993   7000.0          0
```

The three unpruned alignments differed at 33 positions: three routes through the same noise, each
spending 8 to 11% of its unitigs below 1% of the gene's depth. Pruning is also what makes variant
calls trustworthy. Unpruned, the variant stage finds 12 apparently biallelic sites on
`SRR8209082`, every one between 1x and 90x in a 14,500x sample; pruned, it finds one, at 8,000x
against 3,600x.

### 5.2 `--prefer-majority`: put the consensus on the abundant allele

Pruning cannot settle a site where both alleles are real: no coverage floor separates a 35% allele
from a 65% one without also deleting every genuine minority variant, so the aligner still chooses
by resemblance to the reference. `--prefer-majority` compares the two branches at the fork itself
(on `SRR8209082`, 8,000x against 3,600x) and moves the primary onto the more abundant one. That
comparison is local and decisive, where a whole-route average is not (section 7).

### 5.3 `--emit-alts`: report the second allele

At each fork on the route, both branches must clear the abundance and frequency gates; that "both"
separates two real populations from the rare-branch sites of section 3. Each surviving site becomes
one record, the primary with that one allele swapped, in `alt_haplotypes.fasta` rather than the
main FASTA, because downstream every main-FASTA sequence is an independent taxon. Never more than
one site per record: a de Bruijn graph cannot say which alternatives share a molecule, so combining
them would invent sequences ([contigs_vs_unitigs.md](contigs_vs_unitigs.md) on chimeras). Headers
carry `vaf`, `alt_kaf`, `primary_vaf`, `dlen` and `span`; `inherited_pident` is the primary's.

**Indel alleles.** In-frame indels of up to `--max-alt-indel` (30) bases are accepted, with the
alignment window shifted by the length change; a splice that straddles either edge of the window is
refused. Larger changes are the route rewriting the gene rather than editing it. Frameshifting
indels need `--allow-frameshift-alts`: 190 of the 203 credible length-changing calls on the panel
are this case, with a median magnitude of one base, and they cannot be told apart from the
homopolymer-length artifacts a de Bruijn graph produces. `--prefer-majority` applies the same
length rules, so an indel site can be corrected too. `--no-indel-alts` restores substitutions only.

### 5.4 What it recovers

```
                          depth   primary                       alternatives recovered
ERR2028820  PfCRT          180x   1274 nt, byte-identical                 0
ERR899743   ERG11          230x   1574 nt, byte-identical                 0
SRR8209082  HIV-1 pol   14,500x   2737 nt, cov 0.9993                     1   VAF 0.31, 3,600x
SRR6251967  Zika NS5    25,000x   2708 nt, cov 0.9996, 99.82% id         11   VAF 0.10-0.29, 3,000-10,000x
SRR6007473  Ebola GP    25,000x   2030 nt, cov 0.9995, 100.00% id        13   VAF 0.10-0.34, 3,000-13,000x
```

Clonal isolates are untouched; every viral alternative is full length and gapless. Pruning also
improved the consensus where the bias had been worst: Zika coverage 0.9974 to 0.9996 and identity
99.08% to 99.82%, Ebola to a 100.00% match. Lengths one base short of the reference are GraphAligner
clipping the final base ([issue_graphaligner_end_clipping.md](issue_graphaligner_end_clipping.md)).

## 6. Panel run: 200 accessions, three arms

200 accessions over 16 targets, run three ways on one build: no pruning, pruning at 1% of depth,
and pruning plus `--prefer-majority`.

**Pruning neither rescues nor breaks accessions.** 187 of 200 reconstruct in all three arms. The 52
gained over the original prototype's 135 (six targets going from zero to working) come from k=21
anchoring and the Rust port, not pruning; the arms agreeing accession for accession shows it.

**Pruning removes the noise calls.** Alternatives fall from 540 across 91 accessions to 311 across
69. Measured against each sample's own gene depth:

```
support of each alternative     <1% of depth    1-5%       >=5%
  no prune                        94 (17.4%)   19 (3.5%)  427 (79.1%)
  pruned                          10 ( 3.2%)    9 (2.9%)  292 (93.9%)
```

Removals concentrate where sequencing is deepest (`zika_ns5` 165 to 60, `h5n1_pb2` 46 to 9); clonal
targets are untouched (`candida_auris_erg11` 51 and `m_oryzae_avrpik` 7 in both arms).

**Pruning cleans the consensus.** It changes the primary on 84 of 187 accessions. Gaps against the
reference CDS are the objective test, since an indel in a coding sequence is a hard defect for a
codon-aware alignment: 52 of the 84 were gap-free before and 72 after, with fewer gaps on 26 and
more on 5. Identity is higher on 44 and lower on 38, which is what removing a reference bias should
look like: the point is to come off the error branches, not to move toward the reference.

**What `--prefer-majority` changes.** The alternative is the more abundant allele at 59 of 311
sites (19.0%) with pruning alone, and at 1 of 311 (0.3%) with `--prefer-majority`. The maximum
alternative VAF falls from 0.883 to 0.499; staying below 0.5 means the primary carries the majority
allele at every site reported. It changed 26 primaries over 78 sites and nothing else. The
corrections are real: five independent `p_infestans_avr3a` isolates were corrected at the same two
positions, 119 and 132, and all five had matched the reference exactly beforehand.

**One accession not to trust.** `h5n1_ha SRR19878648` is the outlier of the 26: 27 corrections
moved 91 positions in 11 blocks, one 339 nt wide, and put four gaps into an alignment that had
none. Its route dips to 13.5x in a 230x sample, so those majority calls rest on very few reads. It
may be a second influenza clade or a cascade; it is not resolved. The walker warns above 12
corrections, since the other 25 sit between 1 and 11; treat a warned accession as provisional and
compare it with `--no-prefer-majority`.

**Indel alleles, panel-wide.** A later run at the then-current defaults emitted 298 alternatives, 31
of them in-frame indels, and refused 618 frameshifting candidates and 66 over the 30 nt cap. On all
31 the flanks outside the edited stretch are byte-identical to the primary, so nothing downstream of
the indel is shifted. With the length rules shared, 0 of 298 alternatives are above 0.5 VAF.

That run and section 7 reconstruct 173 accessions rather than 187 (the indel run lost a further 8 to
transient DNS failures). The difference is entirely the compressed-size cap then in force, 350 MB,
which refused 14 accessions that otherwise reconstruct. The default is now 1 GB (`--max-size-mb`).

### 6.1 Why these are the defaults

- `--prune-frac 0.01` (was 0.0): no accession lost, 94% of the weakly supported alternatives gone,
  twenty more sequences free of indels.
- `--prefer-majority` (was off): a consensus carrying the rarer base is wrong, and the corrections
  replicate across samples. It changes the sequence feeding a selection analysis, so use
  `--no-prefer-majority` when reproducing an older result.
- `--emit-alts` (was off): it writes its own file and changes nothing downstream.
- `-t 1` (was 2): at `-t 2`, five identical runs produced two different primaries. Use `-w` instead.
- `--min-abund-frac` stays 0.0: it never had a panel arm, and pruning leaves it little to target.

## 7. Why the `detour` engine was retired

A second implementation of the same recovery, the `detour` engine, was run against the current
variant stage over the whole panel at equal cost. With pruning on:

```
                                        variants        detour
alternatives emitted                    267 / 64 acc    269 / 59 acc
  of those, two alleles genuinely       267 (100%)      131 (48.7%)
adds a gap its own primary lacked         6 (2.2%)      114 (42.4%)
primaries gap-free                      160 / 173       151 / 173
```

Equal sensitivity, and a nineteen-fold difference in how often a call introduces an indel. `detour`
also did not tell "two alleles" from "the route took the wrong one": 51% of its calls were the
second, emitted as the first. Its correction ranked candidates on mean `ka:f` over the whole route,
which cannot carry the decision: on `SRR8209082` the chosen `_major` beat the runner-up by
**0.036%**, while one swapped unitig out of ~1,300 moves the mean by 0.08%. Unpruned, `detour` did
find more genuine sites (819 against 469), 203 of them length-changing; that capability became the
in-frame indel support in section 5.3, and the engine was removed.

## 8. Not addressed

Weighting downstream (no selection tool can say "this sequence is 31% of that sample", so that
stays the analyst's call), phasing across sites, and frameshift alleles, which stay opt-in.
