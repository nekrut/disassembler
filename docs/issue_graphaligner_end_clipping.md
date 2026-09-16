# GraphAligner 1.0.19 and later clip the last base of the query

**Status:** diagnosed, not yet fixed. Investigated 2026-09-16.

`logan-walker` emits haplotypes exactly one base shorter than the reference CDS. The missing base
is the last one, which for an in-frame CDS is the final base of the stop codon. This is not a
walker defect. GraphAligner 1.0.19 and 1.0.20 fail to align the final base of a query even when
the query is present in the graph verbatim; 1.0.18 and earlier align it correctly.

The repo's own `environment.yml` pins `graphaligner=1.0.20`, so the documented setup is the
affected one, and the golden fixtures under `bench/golden/` were captured under an affected build.

---

## 1. The short version

| | |
|---|---|
| Symptom | every reconstructed haplotype is `len(reference) - 1` |
| Cause | GraphAligner stops the alignment one query base early |
| First affected version | 1.0.19 |
| Last good version | 1.0.18 |
| Walker involvement | none, the slice arithmetic is correct |
| Is the base recoverable | yes, it is present in the graph |

---

## 2. Minimal reproduction

No Logan data and no walker are involved. Build a one-segment GFA whose sequence *is* the query,
then align the query to it. A perfect self-alignment should be end to end.

```bash
python3 - <<'PY'
seq = ''.join(l.strip() for l in open('queries/pfcrt_only.fasta') if not l.startswith('>'))
open('q.fa', 'w').write('>q\n' + seq + '\n')
open('single.gfa', 'w').write('H\tVN:Z:1.0\nS\t1\t' + seq + '\n')
PY

GraphAligner -g single.gfa -f q.fa -a out.gaf -x dbg -t 1
```

Under 1.0.20 this prints `End-to-end alignments: 0 (0bp)` and writes a single GAF row:

```
q	1275	0	1274	+	>1	1275	0	1274	1274	1274	60	NM:i:0	AS:f:1274	dv:f:0	id:f:1	cg:Z:1274=
```

Reading the columns: query length 1,275, query span `[0, 1274)`, path length 1,275, path span
`[0, 1274)`, 1,274 residue matches, `NM:i:0`, `id:f:1`, CIGAR `1274=`.

So the query is 1,275 bases, the graph holds those same 1,275 bases, the alignment is a flawless
run of 1,274 matches with zero mismatches, and the last base is simply left off. The aligner
declines to align a sequence end to end against an identical copy of itself.

---

## 3. Version bisect

Same GFA, same query, same command, four bioconda builds:

| GraphAligner | `qend` / `qlen` | bases lost | reported |
|---|---|---|---|
| 1.0.17b | 1275 / 1275 | 0 | `End-to-end alignments: 1 (1275bp)` |
| 1.0.18 | 1275 / 1275 | 0 | `End-to-end alignments: 1 (1275bp)` |
| **1.0.19** | 1274 / 1275 | **1** | `End-to-end alignments: 0 (0bp)` |
| **1.0.20** | 1274 / 1275 | **1** | `End-to-end alignments: 0 (0bp)` |

The regression lands between 1.0.18 and 1.0.19.

---

## 4. What it is not

**Not the preset.** `-x dbg` and `-x vg` give identical results on every version tested.

**Not a tunable.** `--X-drop 1000000`, `-b 1000`, `-b 5000` and `-C -1` all still lose the base.
`--precise-clipping` (default 0.66) only makes it worse when raised: 0.9 loses 8 bases, 0.99 loses
43. Values below the default abort.

**Not the introduction of `--precise-clipping`.** The flag is present in 1.0.17b and 1.0.18 too,
so its arrival is not the change.

**Not length-specific.** Reproduced on random 100, 300, 1,000 and 3,000 bp sequences. Usually one
base is lost; the 1,000 bp case lost two, so the size of the loss is not fixed, only its position
at the 3' end.

**Not the walker.** `walker/src/stages/stitch.rs` slices `stitched[path_start..path_end]`, which
is the correct half-open GAF span. The removed Python walker did the identical
`stitched[ps:pe]`. Both faithfully return what the aligner gave them.

---

## 5. On real Logan graphs

`SRR12363134` against `queries/markers/candida_auris_erg11.fasta` (1,575 nt), with
`LOGAN_KEEP_SCRATCH=1` to retain the subgraph and GAF:

```
qlen=1575  qstart=0  qend=1574      query span 1574 of 1575
plen=1594  pstart=9  pend=1583      path span  1574, with 11 bases of path left over
```

The path continues eleven bases past where the alignment stops, so the sequence was in the graph
the whole time. Rebuilding the stitched path by hand from the retained GFA:

```
stitched[9:1583]  = 1,574 nt, what the walker emits today
stitched[9:1584]  = 1,575 nt, byte-identical to the previously committed panel sequence
```

The recovered base is `A`, completing the terminal `TAA`.

---

## 6. It explains a panel that looked unreproducible

The 50-accession *C. auris* ERG11 panel in `nekrut/logan-hypyaeon`
(`work/panel_output/candida_auris_erg11/`) holds 1,575 nt sequences with `cov: 1.0`. Reruns here
produced 1,575 nt under GraphAligner 1.0.17b and old-behaviour flags:

```bash
logan-walker -l erg11_5.txt -q queries/markers/candida_auris_erg11.fasta -o out \
  --prune-frac 0.0 --no-prefer-majority --no-emit-alts
```

```
ERR899743  SRR12363134  SRR14253779  SRR19738638  SRR19738642
5/5 byte-identical to the committed panel
```

Under the pinned 1.0.20 the same command gives 0/50 byte-identical, with all 50 exactly one base
short. Sixteen of the fifty additionally differ by 1 to 3 substitutions, and all sixteen reduce to
the trailing base alone once `--prune-frac 0.0 --no-prefer-majority` is applied, so those
substitutions are the current defaults working as designed and are not part of this issue.

---

## 7. Consequences

**The golden fixtures encode the bug.** All three are `len(reference) - 1`:

| fixture | reference | golden |
|---|---|---|
| `ERR2028820.pfcrt.primary.fasta` | 1,275 | 1,274 |
| `SRR12363134.candida_auris_erg11.primary.fasta` | 1,575 | 1,574 |
| `SRR8209082.hiv1_pol.primary.fasta` | 2,739 | 2,738 |

They were captured under an affected build, so `bench/check_golden.sh` passes on 1.0.20 and would
**fail** on 1.0.18, which produces the biologically correct full-length sequence. The gate is
pinned to the regression and cannot detect it.

**`environment.yml` pins the affected version.** Anyone following the documented setup gets
truncated output.

**Outputs from different eras are not comparable at the 3' end.** Panels produced on the cluster
with an older GraphAligner are full length; anything produced through the pinned env is one base
shorter.

**Downstream impact is small but non-zero.** MAFFT `--keeplength` turns the missing base into a gap
in the final alignment column, so the codon MSA still reports the full reference length and 525
codons for ERG11. Reported coverage becomes 1574/1575 = 0.99937 rather than 1.0, which only matters
against a `--min-cov` of exactly 1.0. Selection analyses that strip the stop codon are unaffected.

---

## 8. Options

1. **Report upstream and pin 1.0.18 meanwhile.** A self-alignment that is not end to end is an
   aligner bug, and 1.0.18 is a drop-in that does not have it. The golden fixtures have to be
   recaptured, since every one of them gains a base.
2. **Pin 1.0.18 and recapture the goldens, without an upstream report.** Same work, leaves the bug
   in place for everyone else.
3. **Compensate in the walker.** Extend the slice by `qlen - qend` at the 3' end, bounded by
   `path_len - path_end`. Verified above to recover exactly the right base. It emits sequence the
   aligner did not confirm, which is a real change in what a reconstructed haplotype means, so it
   should not be done quietly if it is done at all.

Option 1 is preferred. Whichever is chosen, `bench/golden/` must be recaptured and the capture
should record the GraphAligner version alongside the fixtures, so a future version change is
visible rather than silent.

---

## 9. Reproducing this investigation

```bash
conda create -n ga118 -c conda-forge -c bioconda graphaligner=1.0.18 -y
conda create -n ga120 -c conda-forge -c bioconda graphaligner=1.0.20 -y
# then run the section 2 repro under each
```

Related: [`rust_walker_as_built.md`](rust_walker_as_built.md) section 7 on why `SRR8209082` is
reported but not gated, and [`issue_intra_accession_variance_lost.md`](issue_intra_accession_variance_lost.md)
for the pruning and majority-preference defaults that account for the substitutions in section 6.
