# Path walking selects the reference-matching minority allele

**Status:** diagnosed, then fixed. Sections 1 to 13 are the diagnosis and stand as written.
Section 15 records the fix and what it recovers. The Rust walker gained two opt-in flags,
`--prune-frac` and `--emit-alts`; defaults are unchanged.
The Rust walker under `work/rust-walker/` reproduces the existing one-haplotype-per-accession
behaviour deliberately, because the task at hand is speed, not method.

Rewritten 2026-09-11 after a measurement pass. An earlier version of this document, written
2026-09-10 against commit `81fc1fd`, attributed the loss to bubble phasing and to discarded
GraphAligner rows. Both were wrong, and section 3 below says exactly how. The loss is real,
it is larger than first described, and it happens somewhere else.

---

## 1. The short version

The walker emits one sequence per accession, named `<ACC>_primary`. It is chosen by walking
the single path GraphAligner returns, and GraphAligner scores paths by residue matches to the
reference. At forks where the sample's dominant allele differs from the reference and a rare
allele matches it, the walk takes the rare allele.

Measured on `SRR8209082` (HIV-1 *pol*): of the 1,216 forks the walked path passes through,
**50** have an alternative branch that leaves the path and rejoins it a few nodes later,
yielding a full-length alternative haplotype. In **49 of those 50 the discarded branch is the
majority allele**, carrying at least 90% of local coverage. The median coverage ratio between
the branch discarded and the branch walked is **778x**, and the maximum is **10,417x**.

None of these are reported as bubbles. The quasispecies phaser does not see any of them.

---

## 2. How this was measured

`logan_path_walker.py` deletes its GFA and GAF in a `finally` block (`:327-330`), so the
evidence is destroyed on every run. The measurement used a replica of stages 1 to 5, copied
verbatim from the walker, that keeps those files. It reproduces both golden fixtures exactly:

| accession | unitigs | anchors | subgraph nodes | matches `work/benchmark/golden/` |
|---|---|---|---|---|
| `ERR2028820` | 1,277,362 | 76 | 356 | yes |
| `SRR8209082` | 81,534 | 4,102 | 26,323 | yes |

Phasing was reimplemented with the orientation fix described in `rust-walker/src/stages/phase.rs`.
Before trusting any claim about which nodes lie on the walked path, the orientation convention
was validated against the path itself: reading `L:s1:t:s2` in `u`'s header as the traversal
`(u,s1) -> (t,s2)` makes **1278 of 1278** consecutive GAF path steps a real oriented edge. The
alternative reading, flipping the target strand, scores **0 of 1278**.

Alternative haplotypes were built by splicing the detour back into the path and were then
scored with `blastn` against the reference CDS, so "valid" means valid by the same test the
pipeline applies to its own output.

Everything below is reproduced by `work/benchmark/measure_fork_loss.py`:

```bash
python3 work/benchmark/measure_fork_loss.py --all              # both accessions
python3 work/benchmark/measure_fork_loss.py --all --self-check # guards the copied stages
```

`--self-check` exits non-zero unless the copied stages still reproduce the golden graph shape
and the orientation check still passes, so the copy cannot drift from the walker unnoticed.

On determinism, following the convention in `work/benchmark/golden/README.md`: every figure
quoted here is byte-stable across `PYTHONHASHSEED` 1 and 7 and across repeat runs. Reaching
that required sorting the branch and neighbour iteration in the measurement code, which
otherwise walks Python sets of string tuples and picks a different rejoin point per process.
An earlier unsorted version of this measurement gave 48, 50 and 51 clean alternatives on three
runs of the same input, which is the same defect class the golden README documents for
`reconverged[0]`.

---

## 3. What the previous version of this document got wrong

**Wrong: "GraphAligner reports several alignments, the walker keeps one."** It cited
`logan_path_walker.py:193-194`, `rows.sort(...)` then `best = rows[0]`, and presented that as
the site where haplotypes are lost. With the GAF retained, `ERR2028820` produces one row every
time. `SRR8209082` produces one or two, varying between runs of the same input: five repeat
invocations of GraphAligner gave rows 1, 1, 1, 2, 1.

When the second row appears it is not a competing haplotype. It is the same locus scored
identically, and it is flagged ambiguous:

```
row 1   q 0-2737   1279 path nodes   matches 2688   mapq 60
row 2   q 0-2739   1277 path nodes   matches 2688   mapq 0
```

Across all five runs the best row by residue matches was byte-identical, so the emitted
primary is stable. `best = rows[0]` therefore discards a row occasionally and a distinct
haplotype never. It is not where the loss happens.

**Wrong: `res["reconstructed"][0]` at `:379` discards haplotypes.** `align_and_reconstruct`
appends at most one element (`:226`), so the list has length 1 by construction. That index is
not a discard site.

**Wrong: the reported bubbles are discarded alleles of the gene.** Of the 14 bubbles the
corrected phasing finds on `SRR8209082`, **0** have their fork on the walked path. The worked
example the previous version led with, fork 17604 to merge 58278 with branches 43677 at
`ka:f` 12000 and 43767 at 2475, is off the path in every part:

```
unitig 17604: in subgraph=yes  on walked path=no  ka:f=   11.5
unitig 43677: in subgraph=yes  on walked path=no  ka:f=12000.0
unitig 43767: in subgraph=yes  on walked path=no  ka:f= 2475.0
unitig 58278: in subgraph=yes  on walked path=no  ka:f=    4.4
```

The subgraph is a 5-hop neighbourhood, so most of it is not the gene. A bubble in that
neighbourhood is not by itself an allele of the reconstructed CDS, and this one is not.

---

## 4. Where the loss actually happens

The walked path for `SRR8209082` is 1,279 oriented nodes, and 1,216 of them are forks. For
each fork, the branch the path took was compared with the branches it did not, under the
pipeline's own gates (`ka:f >= 5.0`, VAF >= 10%).

```
credible untaken branches at forks on the walked path        58
  ... that rejoin the walked path within 8 hops              57
  ... clean single-allele swaps (cov >= 0.85, <= 50 nt)      50
```

Every one of the 50 is a full-length haplotype, 2,738 nt against a 2,739 nt reference, and
each clears the walker's own acceptance gates:

```
blastn identity to the reference CDS
  primary                98.174%
  50 alternatives        98.028% to 98.174%
coverage of the reference
  primary and all 50     1.000
```

They are not fragments, not chimeras, and not low-quality. By every criterion the pipeline
uses to decide that `<ACC>_primary` is worth keeping, these are worth keeping too.

What separates them from the primary is coverage, in the wrong direction:

```
discarded branch better supported than the branch walked     49 / 50
discarded branch is the majority allele (VAF >= 0.5)         49 / 50
discarded branch carries >= 90% of local coverage            49 / 50
coverage ratio discarded : walked      median 778x, max 10,417x
walked branch < 10x while discarded branch >= 10,000x        16 cases
```

This is not a case of the path being low-coverage overall. It is not:

```
ka:f along the 1,279 walked nodes   min 2.0   median 15,000   max 25,000
ka:f over all 26,323 subgraph nodes min 1.0   median    9.3   max 25,000
walked nodes below the pipeline's own 5.0x floor    12 / 1,279
```

The walk is almost entirely through the dominant population. At roughly 4% of its steps it
drops into a near-absent allele, because that allele matches the reference.

---

## 5. One case in full

Path index 1236, near the 3' end of *pol*:

```
fork (54954,-), reached from (54872,-)
  branch (13182,+)  ka:f     18.9  VAF 0.002   len 35
  branch (54790,+)  ka:f  12000.0  VAF 0.997   len 31
  branch (54791,+)  ka:f     10.0  VAF 0.001   len 31
  branch (54792,+)  ka:f      2.0  VAF 0.0004  len 61   <- walked
```

The discarded branch rejoins the walked path two nodes later, via
`(54790,+) -> (11952,+) -> (66696,+)`, so the substitution is clean and local. The two
resulting haplotypes are both 2,738 nt and differ at two adjacent positions:

```
pos 2639   primary T   alternative A    TAATAGTG[T/A]CATAAAAG
pos 2640   primary C   alternative A    AATAGTGT[C/A]ATAAAAGT
```

The sequence that goes forward into MAFFT, FastTree and HyPhy carries `TC` at a site where
99.7% of the sample's k-mer coverage says `AA`.

---

## 6. Why the phaser cannot see any of this

`phase_superbubbles` requires at least two branches to clear both gates
(`logan_path_walker.py:263-269`). Classifying the same 58 forks:

```
no reconvergence at the depth the phaser tests                50
bubble-shaped, but suppressed by the VAF gate                  8
reported as a bubble                                           0
```

The suppression is self-inflicted. A bubble is only reported when two branches each hold at
least 10% of local coverage. At these forks one branch holds 99% and the branch the walker
took holds well under 1%, so `len(valid_branches) < 2` and nothing is reported. The gate
filters out the site precisely because the walk chose the rare allele there.

The remaining 50 do not reconverge within one step, which is all `phase_superbubbles` tests
(`:253-258`). They reconverge at depths 2 to 4, or rejoin the path further along. A one-step
reconvergence test finds clean SNP bubbles and little else.

---

## 7. The reference pull, previously flagged as unverified

The earlier version of this document raised, in a section explicitly labelled unverified,
that sorting GAF rows by column 10 selects the path most similar to the reference rather than
the most abundant. The underlying selection happens inside GraphAligner rather than in that
sort, but the effect is now measured:

```
alternatives LESS identical to the reference than the primary   40 / 50
alternatives MORE identical                                      0 / 50
```

Not one alternative beats the primary on identity to the reference. That is what a
reference-maximising objective produces, and it is the shape that matters for selection
analysis: the divergence being suppressed is exactly the divergence dN/dS is computed from.
The bias is conservative, and it is systematic rather than occasional.

---

## 8. This is sample dependent

The previous version generalised from one accession to every accession. That does not hold.
`ERR2028820`, PfCRT in *P. falciparum*, shows none of the effect:

```
                                  ERR2028820      SRR8209082
walked path nodes                        115           1,279
forks on the walked path                  63           1,216
credible untaken branches                  4              58
clean alternative haplotypes               0              50
bubbles after the orientation fix          0              14
ka:f along the path, median              180          15,000
walked nodes below 10x                     0              25
```

Of `ERR2028820`'s 4 credible untaken branches, 3 rejoin the path and none survives: one
rewrites 827 nt of the sequence rather than swapping an allele, and the other two fall to
0.669 and 0.704 coverage, below the 0.85 `min_cov` gate. A clonal parasite isolate at a few hundred x coverage does not have the
structure that a deep viral sample has. The effect belongs to high-diversity samples.

---

## 9. The bubble detail is still thrown away at the source

Independent of everything above, and unchanged from the earlier version of this document:
`phase_superbubbles` computes each fork, each branch, its `ka:f` and its VAF, and then only
the count survives.

```
317:      "bubbles": len(bubbles),
385:          bubble_reports.append({"acc": acc, "bubbles": res["bubbles"]})
```

Tracing every consumer:

| consumer | reads |
|---|---|
| `build_walker_codon_msa.py` | the reconstructed FASTA. Bubbles are not an input. |
| `run_panel_analysis.py` | no reference to bubbles anywhere in the file |
| MAFFT, FastTree, HyphAeon | the codon MSA, built only from the FASTA |
| `summarize_all_targets.py:164` | `total_bubbles += b.get("bubbles", 0)` |
| `generate_report_figures.py:133` | `bubbles = [r['bubbles'] for r in walker_data['results']]` |

The only thing done with a bubble is to add its count to a total and colour a scatter plot by
it. A reported quasispecies figure cannot be audited by anyone, including the person who
produced it. The Rust walker already fixes this much: `Bubble` carries fork, merge, and every
branch with its orientation, `ka:f` and VAF, and `io/output.rs` writes `bubbles.json`.

Note that those counts are also inflated. The committed *C. auris* ERG11 run reports 20,256
bubbles across 50 accessions, and the report figure of 1.57M inherits the same artifact: the
strandless adjacency in `phase_superbubbles` reads a chain as a fork and a bidirected edge as
a reconvergence, giving `fork == merge` on 93.1% of `ERR2028820`'s 202 bubbles and 58.5% of
`SRR8209082`'s 11,825. Corrected, the same subgraphs hold 0 and 14.

---

## 10. What the committed output looks like

Confirmed against `panel_output/candida_auris_erg11/`:

```
haplotypes per successful accession -> {1: 50}
all ids end in _primary               -> True (50 records)
total bubbles reported                -> 20256
codon MSA: 51 records, 6 unique sequences
```

51 is 50 accessions plus the reference, and `candida_auris_erg11.busted.json` accordingly
reports `taxa: 6`. Statistical power is six sequences.

---

## 11. The documentation disagrees with itself

`HANDOVER.md` states the limitation plainly:

> line 323: Logan = assembled unitigs -> within-sample variation collapsed; frequencies are
> per-haplotype not population; **not phased**.

`README.md` reads the other way:

> line 109: executes orientation-aware bubble phasing to reconstruct primary haplotypes
> **alongside genuine intra-host quasispecies**

> line 480: reconstruct the full-length CDS, and **phase intra-host quasispecies without
> consensus loss**

> line 524: **Orientation-aware superbubble phasing** ... to phase genuine within-host
> quasispecies

"Alongside", "without consensus loss" and "phase" all imply the variants reach the output.
They do not. The `HANDOVER.md` note is the accurate one. On the evidence in section 4,
"without consensus loss" is the wrong way round: at the sites where the sample disagrees with
the reference, the emitted sequence is not the consensus.

---

## 12. What is still genuinely preserved

**Cross-sample diversity.** Selection inference operates on differences between the sequences
in the alignment, placed on a tree. One haplotype per accession across thousands of runs is
real population-scale variation. Section 7 means it is biased toward the reference, not that
it is absent.

**Chimera avoidance.** Requiring a connected path through the graph is what stops fragments of
different strains being stitched into sequences that never existed. README section 4 measures
naive coordinate stitching at 87% chimeric artifacts. That benefit is fully retained, and the
alternatives in section 4 are connected paths too, so emitting them would not forfeit it.

**Sensitivity to divergent strains.** Walking a graph finds a gene whose sequence differs
substantially from the reference, where a collapsed contig consensus might not.

---

## 13. What changing this would involve

Two separable changes, in increasing order of difficulty.

**Pick the abundant branch, not the reference-matching one.** This is the section 7 bias and
it does not require emitting anything new. The walk already has `ka:f` for every branch. A
path that breaks ties, or overrides GraphAligner, by coverage would emit the actual consensus
of the sample. On the 50 sites measured here it would change the primary haplotype at 49 of
them. This is the change with the clearest justification: the current output is not the
consensus of anything, and no downstream consumer wants the 2x allele over the 12,000x one.

**Emit the alternatives as additional records.** The machinery is most of the way there. The
detour is local, it rejoins within a few nodes, and the spliced sequence is already known to
be valid. A second record such as `<ACC>_alt1` carrying its VAF in the header is mechanical.
Two genuine design problems come with it, which is very likely why it was not done.

*Weighting in the tree.* Forty near-identical sequences from one outbreak already enter as
forty independent observations. HyPhy has no way to express "this sequence is 17% of that
sample", so a naive addition inflates apparent diversity rather than describing it.

*Which combinations to emit.* With *n* independent forks there are up to 2^n haplotypes, and a
de Bruijn graph cannot say which alternatives sit on the same molecule. That is the "not
phased" in the `HANDOVER.md` note. Emitting every combination would invent sequences that may
not exist. Emitting one variant at a time is defensible but arbitrary.

Neither problem applies to the first change, and neither is a reason to leave the per-bubble
detail unexported.

---

## 14. Related

* `work/benchmark/measure_fork_loss.py`, which reproduces every figure in this document.
* `work/plan_rust_path_walker.md` and `work/rust_walker_as_built.md`, the speed work this was
  found during.
* `work/benchmark/golden/`, whose bubble counts pin the pre-fix strandless numbers and will
  need regenerating alongside any phasing change.
* Worth noting for whoever does that regeneration: `phase_real` in
  `rust-walker/src/stages/phase.rs` builds its subgraph from the Rust `anchor::find` rather
  than from `blastn`, and gets 422 and 8,466 subgraph nodes where the Python walker gets 356
  and 26,323. The two anchoring paths do not agree, so bubble counts measured through that
  test are not directly comparable with the golden fixtures. The Rust test reports 1 and 11
  bubbles on its subgraphs; the corrected phasing on the Python subgraphs gives 0 and 14. The
  conclusions in this document are all measured on the Python subgraphs, which are the ones
  the pipeline actually uses.

---

## 15. Resolution

Two changes, both off by default, in `work/rust-walker/`.

### 15.1 `--prune-frac`: remove the reference bias

The defect in section 4 is that GraphAligner scores a route by resemblance to the reference,
so where the sample's real sequence differs from the reference and an error branch matches
it, the error wins. The fix is to delete those branches from the search space before the
aligner sees them, which needs the gene's sequencing depth.

That number is not available before an alignment, and cannot be taken from the anchor set:
`SRR8209082`'s anchors have a median `ka:f` of 70x against a true 14,500x, because
`min_shared_kmers = 3` admits thousands of shallow incidental matches. So the walker now
aligns once to measure depth, drops unitigs below `prune_frac` of it, and aligns again. If
the second alignment fails, the first pass's haplotype stands, so this can improve an
accession and cannot break one.

```
SRR8209082 (HIV-1 pol, 14,500x)   nodes  aligns    cov   min ka:f  under 1% of depth
  unpruned                        34582     3    0.9996      2.0        118
  1% floor (145x)                  6930     1    0.9993    147.5          2
  5% floor (725x)                  5343     1    0.9993   7000.0          0
```

The three alignments collapsing to one is the intended result, not a loss. They differed at
33 positions and were three routes through the same noise, each spending 8 to 11% of its
unitigs below 1% of the gene's depth.

### 15.2 `--emit-alts`: recover the second allele

At each fork the route passes through, both branches must clear the abundance and frequency
gates. That "both" is the entire filter, and it is what separates two real populations from
the 54 sites in section 4 where the aligner took a rare branch. Each surviving site is
emitted as one extra record: the consensus with exactly one substitution, written to
`alt_haplotypes.fasta` rather than the main FASTA, because every sequence in the main FASTA
becomes an independent taxon downstream and two records from one sample are not two
independent observations.

Never more than one substitution per record. With *n* sites there are up to 2^n haplotypes
and the graph cannot say which alternatives share a molecule, so combining them would invent
sequences, which is the chimera problem the walker exists to prevent.

### 15.3 What it recovers

Pruning is what makes the variant calls trustworthy. On `SRR8209082` unpruned, the stage
finds 12 apparently biallelic sites, every one between 1x and 90x in a 14,500x sample. Pruned,
it finds one, at 8,000x against 3,600x.

```
                          depth   primary                       alternatives recovered
ERR2028820  PfCRT          180x   1274 nt, byte-identical                 0
ERR899743   ERG11          230x   1574 nt, byte-identical                 0
SRR8209082  HIV-1 pol   14,500x   2737 nt, cov 0.9993                     1   VAF 0.31, 3,600x
SRR6251967  Zika NS5    25,000x   2708 nt, cov 0.9996, 99.82% id         11   VAF 0.10-0.29, 3,000-10,000x
SRR6007473  Ebola GP    25,000x   2030 nt, cov 0.9995, 100.00% id        13   VAF 0.10-0.34, 3,000-13,000x
```

Clonal isolates are untouched: identical haplotype with and without pruning, no alternatives.
Every alternative on the three viral samples is full length, gapless against the reference,
and differs from the consensus at exactly one position, and on `SRR6251967` the eleven sites
are eleven distinct positions spread across the gene. That is the intra-host quasispecies the
pipeline was documented as capturing and did not.

Pruning also improved the consensus itself where the reference bias had been worst: Zika
coverage 0.9974 to 0.9996 and identity 99.08% to 99.82%, Ebola to a 100.00% match.

### 15.4 Not addressed

* **Indel alleles.** An alternative that changes the sequence length is rejected, because
  `path_start`/`path_end` are offsets into the stitched sequence and a shifted one would
  silently reframe the CDS. On `ERR2028820` this rejects three candidates at +1, +61 and
  +125 nt; the +125 would have rewritten a chunk of the gene. Supporting them means adjusting
  those offsets, and the failure mode is a frameshifted coding sequence entering a selection
  analysis, so it should be done deliberately.
* **Weighting downstream.** HyPhy has no way to express "this sequence is 31% of that
  sample". The separate output file leaves that decision open rather than making it silently.
* **`min_abund_frac`.** A depth-relative floor for the variant gates exists and is off. With
  pruning on it has little left to do, since pruning has already removed what it targeted.

---

## 16. Panel run, 200 accessions, three arms

Every accession in `work/panel_output/` (16 targets, 200 accessions), run three ways on the
same build: no pruning, pruning at 1% of depth, and pruning plus `--prefer-majority`.

### 16.1 Pruning does not rescue accessions, and does not break any

```
reconstructed        no prune 187 / 200
                     pruned   187 / 200
                     pruned + prefer-majority   187 / 200
```

Identical in all three arms. The panel does reconstruct 187 against the committed Python
walker's 135, 52 more, with six targets going from zero to working (`dengue_1`, `dengue_2`,
`measles_h`, `rabies_g`, `rsv_a_f`, `westnile_ns5`). **None of that is pruning.** It is the
k=21 anchoring from upstream's `3aabf8f` plus the Rust port, and the three arms agreeing to
the accession is what establishes that.

### 16.2 What pruning does: it removes the noise calls

```
alternatives emitted      no prune 540 across 91 accessions
                          pruned   311 across 69 accessions

support of each alternative, against its own sample's gene depth
                          <1% of depth    1-5%     >=5%
  no prune                  94 (17.4%)   19 (3.5%)  427 (79.1%)
  pruned                    10 ( 3.2%)    9 (2.9%)  292 (93.9%)
```

Pruning removes 229 alternatives and 94% of the sub-1% ones. The concentration is where the
sequencing is deepest: `zika_ns5` 165 to 60, `h5n1_pb2` 46 to 9, `mpox_f13l` 39 to 13,
`ebola_gp` 59 to 35. Clonal targets are untouched, `candida_auris_erg11` 51 in both arms and
`m_oryzae_avrpik` 7 in both.

### 16.3 What pruning does to the consensus

It changes the emitted haplotype on 84 of 187 accessions, and the objective test is gaps
against the reference CDS, since an indel in a coding sequence is a hard defect for a
codon-aware alignment:

```
gap-free before pruning   52 of 84
gap-free after pruning    72 of 84
fewer gaps after          26
more gaps after            5
identity higher after     44
identity lower after      38
```

A net twenty sequences stop carrying an indel. Identity to the reference is near evenly
split, which is what removing a reference bias should look like: it is not supposed to move
the consensus toward the reference, only off the error branches.

The largest position-level differences are indel offsets rather than substitutions. On
`SRR6007473` the unpruned haplotype aligns at 99.115% with five gaps and the pruned one at
100.000% with none; comparing them base by base counts 945 differences, all of it the frame
offset that the gap created.

### 16.4 What `--prefer-majority` does

```
sites where the alternative is the MORE abundant allele
  pruned                      59 of 311  (19.0%)
  pruned + prefer-majority     1 of 311  ( 0.3%)

alternative VAF, maximum      0.883  ->  0.499
```

A maximum below 0.5 is the invariant: after correction the primary carries the majority
allele at every site it reports. It changed 26 primaries over 78 sites, and nothing else:
every difference between the two arms traces to a reported correction, none unexplained.

Evidence the corrections are real rather than noise: five independent `p_infestans_avr3a`
isolates were corrected at the **same two positions**, 119 and 132, with a sixth at the same
sites offset by its alignment start. All five had matched the reference exactly beforehand.
Agreement across five independent samples at identical coordinates is not something a noise
process produces, and a perfect reference match is the symptom this exists to remove.

### 16.5 One accession not to trust

`h5n1_ha SRR19878648` is the single outlier of the 26. It took 27 corrections that moved 91
positions in 11 blocks across the gene, one of them 339 nt wide, and put four gaps into an
alignment that had none. Its route also dips to 13.5x in a 230x sample, so the majority call
at those sites rests on very few reads. It may be a second influenza clade in the sample or a
cascade; it is not resolved here. `prefer_majority` now warns above 12 corrections, since the
other 25 sit between 1 and 11.

### 16.6 Defaults

Set from the panel evidence above:

```
prune_frac            0.01     was 0.0
prefer_majority       true     was false
emit_alts             true     was false
graphaligner_threads  1        was 2
min_abund_frac        0.0      unchanged, no panel evidence either way
```

`--prune-frac 0.01` because it never changed whether an accession reconstructed, 187 of 200
in all three arms, while removing 94% of the weakly supported alternatives and leaving twenty
more sequences free of indels.

`prefer_majority` because a consensus carrying the rarer base at a site is simply wrong, and
because the corrections replicate across independent samples (16.4). It is the one default
here that changes the sequence feeding the selection analysis, so `--no-prefer-majority`
restores the old behaviour when reproducing an older result.

`emit_alts` because it writes to its own file and changes nothing downstream.

`min_abund_frac` stays at 0.0. It was never given a panel arm, and with pruning on it has
little left to target, so defaulting it would be guessing.

Old behaviour in full: `--prune-frac 0.0 --no-prefer-majority --no-emit-alts`.

---

## 17. variants.rs against detour.rs

Two implementations of the same recovery landed in this tree, upstream's `stages::detour`
and the local `stages::variants`. Both were run over the whole 200-accession panel twice,
with pruning off and on, on the same build. Cost is not a discriminator: the haplotype stage
has a median of 0.001s in both, against a 1.7s median for the accession as a whole.

### 17.1 With pruning on, which is the default

```
                                        variants        detour
alternatives emitted                    267 / 64 acc    269 / 59 acc
  of those, two alleles genuinely       267 (100%)      131 (48.7%)
  route had taken the weaker branch       0             40  (14.9%)
  route had taken a <10% branch           0             98  (36.4%)
full length (>= 0.95 of the reference)  266 / 267       268 / 269
adds a gap its own primary lacked         6 (2.2%)      114 (42.4%)
primaries gap-free                      160 / 173       151 / 173
extras on clonal targets                 67              74
```

Equal sensitivity in sites found, and a nineteen-fold difference in the rate at which a call
introduces an indel against the reference. An indel in a coding sequence is a hard defect for
a codon-aware alignment, so that gap rate is the single most important number here.

The VAF split is the other one. Both stages are looking at the same forks, but `detour` does
not distinguish "this sample carries two alleles" from "the route took the wrong one": 51% of
its calls are the second case, emitted as though they were the first. `variants` resolves
those into the consensus first (`--prefer-majority`) and then reports only what is left, so
every record it emits has the primary on the majority allele and the alternative below 50%.

### 17.2 With pruning off

```
                                        variants        detour
alternatives emitted                    469             2137
  VAF 0.10-0.50, two real alleles       469 (100%)       819 (38.3%)
  VAF >= 0.90, route took a <10% branch   0             1175 (55.0%)
full length                             469 / 469       2100 / 2137
adds a gap its own primary lacked        40 (8.5%)       465 (21.8%)
plasmid_mcr1 (clonal, 10 isolates)        0               5
```

Unpruned, `detour` finds more genuine sites, 819 against 469. The gap is two things: 203 of
them are length-changing alleles, which `variants` rejects deliberately, and the rest are
mostly the cap at `max_alts = 32`. That is a real sensitivity advantage and it is the part
worth taking.

### 17.3 Why detour's own correction does not work

`detour` handles the reference-bias case by ranking every candidate on `mean_ka`, the average
coverage over the whole route, and naming the top one `_major`. The statistic cannot carry
that decision. On `SRR8209082` it produced 55 candidates whose mean coverage spans 1.65%, and
the gap between the chosen `_major` and the runner-up was **0.036%**. One swapped unitig out
of ~1,300 moves the mean by 0.08%, so the ranking is decided well below the resolution of the
thing being ranked.

Comparing the two branches *at the fork*, which is what `--prefer-majority` does, is a
comparison between 8,000x and 3,600x. Local and decisive, where the global average is not.

Two smaller defects in the same output: `pident` is copied from the primary onto every
alternative, so all 55 records claim 98.2% however much they differ; and `vaf` is the
frequency at the defining fork rather than the haplotype's own share, so `_major` can carry
`vaf=0.97` while `minor29` carries `vaf=1.00`.

### 17.4 What to do

Keep `variants` as the implementation and take one thing from `detour`: alleles that change
the sequence length. Those are 203 of its credible calls and `variants` currently rejects
every one, because `path_start`/`path_end` are offsets into the stitched sequence and a
shifted one silently reframes the CDS. The fix is to adjust those offsets by the length
delta, which is mechanical, and the reason it has not been done is that the failure mode is a
frameshifted coding sequence reaching a selection analysis. It needs its own validation pass,
not a flag flip.

Not recommended: adopting `detour` as-is, which puts a 42% gap-introduction rate and 51%
mislabelled calls into `reconstructed.fasta` and therefore into the codon alignment; or
running both, which emits the same alleles twice under two schemes.

---

## 18. Indel alleles

Section 17.4 named the one capability worth taking from `detour`: alleles that change the
sequence length. `variants` refused every one of them, because `path_start` and `path_end`
are offsets into the stitched sequence and reusing them after a splice that adds or removes
bases slices the wrong window. `variants::shifted_window` now does the arithmetic, and the
shift is only applied when the splice sits entirely inside the alignment window, so a splice
straddling either edge is still refused rather than guessed at.

### 18.1 Why most of them are still refused

Of the 203 length-changing calls `detour` makes in the credible band across the panel:

```
in-frame, delta divisible by 3       13   ( 6.4%)
frameshifting                       190   (93.6%)
most common deltas   +1 (66)  -1 (42)  +100 (32)  -2 (10)  +142 (10)
|delta|              median 1, max 360; 64 calls exceed 60 nt
```

So the class is not one thing. The +100, +142 and -282 tail is the route being sent through a
different stretch of graph, which rewrites the gene rather than editing it, and
`max_alt_indel = 30` excludes it. The ±1 and ±2 bulk shifts the reading frame of everything
downstream, and cannot be told apart here from the homopolymer-length artifacts a de Bruijn
graph produces. A frameshifted record does more damage to a codon-aware selection analysis
than a missing one, so `allow_frameshift_alts` is off by default.

What that leaves is the in-frame subset, which is a real class: a codon inserted or deleted.

### 18.2 What it recovers, over the whole panel

200 accessions, current defaults:

```
primaries                           165      (173 less 8 transient DNS failures)
alternatives                        298
  of those, indel alleles            31
refused, frameshifting              618      twenty times what it admits
refused, over the 30 nt cap          66
refused, window not shiftable         0
refused, no rejoin                  476

deltas  -18(1) -3(16) +3(5) +6(2) +9(1) +15(1) +18(1) +21(3) +24(1)
target  mpox_f13l 15, ebola_gp 14, h5n1_pb2 1, rsv_b_f 1
```

Every check holds:

```
in-frame, delta divisible by 3                      31 / 31
length invariant, allele == primary + delta          0 violations
full length against the reference (>= 0.95)         30 / 31
clean local edit, flanks byte-identical             31 / 31
```

The last one is the real frame test and it is exact: outside the edited stretch, the allele
and its primary are byte-identical, so nothing downstream of the indel is shifted or
scrambled. Translated stop counts move on 12 of 31, which is expected rather than alarming,
because an inserted or deleted codon can itself be a stop; the flank test is what rules out a
frameshift, and it passes on all 31.

The window guard never fired across 200 accessions, so the coordinate arithmetic held
wherever it was reached.

### 18.3 Two things the panel exposed

**The two stages have to share the length rules.** `find_alts` accepted in-frame indels while
`prefer_majority` still demanded equal length, so at an indel site the consensus could not be
corrected and the allele was emitted holding the majority. Measured on `ebola_gp`: 1 of 48
alternatives, an indel. With the rules shared, the panel-wide count of alternatives above 0.5
VAF is **0 of 298**. Sharing them also means the baseline length has to move when a
length-changing swap lands, or later swaps are measured against the wrong route.

**The cap bounds the net change, not the edit.** `max_alt_indel` limits how much an allele
adds or removes, which is not the same as how much sequence it replaces, and the two come
apart: one allele nets +21 nt while replacing 194.

```
replaced span      substitution alleles (267)    indel alleles (31)
  1-3 nt                209                          4
  4-30 nt                44                         22
  31-100 nt              14                          3
  >100 nt                 0                          1
  median                  1 nt                       9 nt
```

No cap is applied on the span, because one would discard the phased blocks along with the
reroutings: an allele carrying several linked substitutions across 50 nt is a haplotype block
worth having, and that is exactly what a span cap would remove. The span is written into each
record's header instead, so the distinction is a reader's to make rather than a filter's.

### 18.4 Settings

```
allow_indel_alts        true    substitutions and in-frame indels
max_alt_indel           30      beyond this it is a rerouting, not an allele
allow_frameshift_alts   false   opt in with --allow-frameshift-alts
```

`--no-indel-alts` restores the substitutions-only behaviour that held before this.

---

## 19. The 350 MB size cap, measured

Upstream's `8b692ac` lowered `max_compressed_bytes` from 2 GB to 350 MB, so accessions above that
are refused on Content-Length before anything is downloaded. That is a sound idea for environmental
metagenomes, but it is not free on this panel.

```
accessions over 350 MB                          18 of 200
  ...that DID reconstruct under the 2 GB cap    14   <- these are now refused
  ...that failed for other reasons anyway        4
accessions over 2 GB                             3 of 200
largest: SRR12016917 3,982 MB, SRR20458000 2,991 MB, SRR12487681 2,303 MB
```

So the cap costs **14 of the 187 successful reconstructions, 7.5%**, including `SRR6252011`
(1,054 MB, zika), `SRR15530513` (1,642 MB) and `SRR24066904` (1,441 MB). It also happens to drop
`SRR19878648` (767 MB), the H5N1 accession section 16.5 flags as untrustworthy.

Only 3 of 200 exceed the old 2 GB cap, so a value between the two, around 1.5 GB, would refuse the
genuine metagenomes and keep all 14. The cap is currently left at upstream's 350 MB; this is recorded
so the trade is visible rather than discovered later as a drop in yield.

A related reading note: a panel run at the current default reports 173 reconstructions where earlier
sections report 187. The difference is entirely this cap, not a regression in the walker.
