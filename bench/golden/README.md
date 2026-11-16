# Golden-truth fixtures

Concordance baselines for `logan-walker`, in two sets, both checked by `bench/check_golden.sh`:

- **old-behaviour** (this directory): the primary haplotype each accession must reproduce byte
  for byte under `--prune-frac 0.0 --no-prefer-majority --no-emit-alts`.
- **default-flag** (`default/`): every output the walk produces with no flags at all, which is
  what users actually run.

## Current fixtures (recaptured 2026-10-18)

| Accession | Gene | Query | Primary | Reference | Gate |
|---|---|---|---|---|---|
| `ERR2028820` | PfCRT | `queries/pfcrt_only.fasta` | 1,275 nt | 1,275 nt | byte-exact |
| `SRR12363134` | *C. auris* ERG11 | `queries/markers/candida_auris_erg11.fasta` | 1,575 nt | 1,575 nt | byte-exact |
| `SRR8209082` | HIV-1 pol | `queries/query_hiv1_pol.fasta` | 2,740 nt | 2,739 nt | reported only |

Captured from the Rust walker (`logan-walker` 0.1.0, commit `f655edd`) under GraphAligner
1.0.17b, with `--prune-frac 0.0 --no-prefer-majority --no-emit-alts`. Each `*.primary.fasta` is
the walker's FASTA output verbatim. The aligner version is recorded in `GRAPHALIGNER_VERSION`, and
`check_golden.sh` warns when it runs under a different one.

The two gated primaries are the previous fixtures plus their final base, which restores the stop
codon that GraphAligner 1.0.19 and 1.0.20 clip (`docs/issue_graphaligner_end_clipping.md`). Under
1.0.20 the gate now fails on exactly that base. `SRR8209082` takes a different unpruned route under
1.0.17b than in the earlier capture; it is reported, not gated (see `check_golden.sh`).

## Default-flag fixtures (`default/`, captured 2026-11-13)

Captured from `logan-walker` 0.1.1 at commit `bf7560c` under GraphAligner 1.0.17b with
`graphaligner_threads` at its default of 1, with no flags but `--emit-vcf --emit-snv-table`
(a writer that is off by default produces nothing to gate).

| Accession | Primary | Alts | Bubbles | SNVs | VCF |
|---|---|---|---|---|---|
| `ERR2028820` | 1,275 nt | none | 1 | none | none |
| `SRR12363134` | 1,575 nt | none | 0 | none | none |
| `SRR8209082` | 2,738 nt | 1 (`vaf=0.310`) | 3,632 | 1 | 1 row |

Why this set exists: the old-behaviour gate switches off three of the four steps that shape a
default primary, and compares only the bare sequence. So the prune/realign acceptance,
`variants::prefer_majority`, and the alt, VCF and SNV-table writers had no end-to-end check
anywhere in the tree (review finding TC-02).

Three things worth knowing before trusting a green run:

1. **`SRR8209082` is gated here, unlike above.** Its route is stabilised by the default pruning
   that the old gate turns off, and three independent default runs on 2026-11-13 agreed byte for
   byte across all five outputs. That evidence covers one machine at one thread only;
   GraphAligner is not deterministic on a graph this size with more than one thread. If it flaps
   in CI, demote it to the informational section of the script rather than removing the gate.
2. **Its default primary is 2,738 nt against the old-behaviour fixture's 2,740 nt.** It is the
   one accession here whose answer the defaults actually change, which is why leaving it ungated
   left the largest hole. For the other two the default walk reproduces the same primary
   sequence despite pruning being active (479 of 747 and 135 of 296 unitigs kept), so there
   these fixtures gate the header statistics, bubbles, alts and SNV table rather than the
   sequence.
3. **None of the three exercises the majority correction.** No run printed a `consensus:` line:
   `SRR8209082`'s walk already took the majority branch (`primary_vaf=0.690` against the alt's
   `vaf=0.310`), `ERR2028820`'s four biallelic candidates were all rejected for unrelated
   reasons, and `SRR12363134` has no on-path fork. So a green run here says nothing about
   `prefer_majority`; the only coverage of that path is `walker/tests/cli.rs`.

### `SRR8209082` was recaptured for WV-02 (2026-11-15)

Its default-flag fixtures are the only ones WV-02 changed, and they changed on purpose.
`prefer_majority` now moves the route off a branch the aligner walked that is itself a minority
of its fork, when the other side holds an outright majority. On this accession that is 20
corrections, all of them rescues the pre-0.2.0 gate refused to make.

The new primary was checked before the fixture was replaced, because 20 corrections is the
shape the stage warns about:

- same length, 2,738 nt, and **no new stop codon** (0 before, 0 after);
- **10 substitutions, every one a transition**, spread across the gene;
- **7 synonymous to 3 non-synonymous**, and all three amino-acid changes conservative
  (R422K, I737V, I858V);
- the corrected route's mean `ka:f` **rises from 16,987 to 17,758**, so it sits on
  better-covered unitigs than the one it replaced.

That is recovered allelic variation, not a rerouted walk. `alt_haplotypes.fasta` gains a second
record and the VCF and SNV table each gain a row, which is the other half of WV-02: those
alleles were previously neither corrected nor reported. `ERR2028820` and `SRR12363134` are
byte-identical through the change; only their summaries moved, for the new fields.

### The summary is gated too (new in 0.2.0)

Each default-flag accession also has a `<tag>.summary.json` fixture, compared with its
volatile fields masked by `bench/mask_summary.py`: the clock (`run.started`,
`run.wall_seconds`, `run.stage_seconds`, each record's `timing`), the machine
(`tool.external[].path`, `inputs.source.id` and `kind`, since a cached CI run and a cold local
one differ there), and the versions under test (`tool.version`, `tool.disassembler`).

Everything else is gated, which is the point. Every count in `discards` and `cohort`,
every `seq_sha256`, every setting and the block order are all compared byte for byte. Those
numbers existed only on stderr until 0.2.0, so nothing could regress them visibly: a walk that
started refusing every candidate allele, or silently stopped majority-correcting, produced the
same gated outputs as one that did not.

To regenerate: run the gate, then copy each `<query stem>.summary.json` to
`default/<tag>.summary.json`.

`alts.fasta` is empty for the two clonal isolates, and a **missing** VCF fixture is itself an
assertion: those two recover no SNPs, so a VCF appearing is a real change. `SRR8209082`'s
`bubbles.json` is gated by sha256 rather than committed, because 3,632 bubbles make a 1 MB file
that no one could read as a diff and that is regenerable in a minute.

No normalisation is applied. The VCF carries no `##fileDate` and no command line, and the FASTA
headers hold only walk statistics, so every file is byte-comparable as written.

To regenerate: run the walker as above, then copy `primary.fasta`, `alt_haplotypes.fasta`,
`snvs.tsv`, `<acc>.vcf` and `bubbles.json` to `default/<tag>.{primary.fasta,alts.fasta,snvs.tsv,
vcf,bubbles.json}`, or record its sha256 where the file is large.

## Original capture (2026-09-10)

| Accession | Gene | Query | Unitigs | Anchors | Subgraph | Primary | Bubbles |
|---|---|---|---|---|---|---|---|
| `ERR2028820` | PfCRT | `queries/pfcrt_only.fasta` | 1,277,362 | 76 | 356 | 1,274 nt | 202 |
| `SRR8209082` | HIV-1 pol | `queries/query_hiv1_pol.fasta` | 81,534 | 4,102 | 26,323 | 2,738 nt | 11,825 |

Captured from the original Python prototype walker, under GraphAligner 1.0.20 and BLAST 2.17.0+.
The graph statistics above and the `*.golden.json` files are from this capture.

## Usage

```bash
bench/check_golden.sh              # exit 0 = concordant, 1 = regression
```

`check_golden.sh` runs the Rust `logan-walker` under old-behaviour flags and diffs the primary
byte-for-byte against the `*.primary.fasta` fixtures here (needs `GraphAligner` on `PATH` and S3
access). Until 2026-10-18 the `ERR2028820.pfcrt` and `SRR8209082.hiv1_pol` fixtures were the
original Python walker's output, and `SRR12363134.candida_auris_erg11` was
captured from the Rust walker and cross-validated against the Python walker (identical 1,574 nt)
before that walker was removed. To regenerate a fixture, run the Rust walker with
`--prune-frac 0.0 --no-prefer-majority --no-emit-alts` under the GraphAligner version in
`environment.yml`, save the primary, and update `GRAPHALIGNER_VERSION`. The retired
Python-based `capture_golden.py`/`verify_golden.py` are under `../archive/`.

## What is asserted, and what is not

**What the live gate actually checks.** In the old-behaviour set, `check_golden.sh`
byte-compares the primary sequence against the `*.primary.fasta` fixtures and nothing else. The
`default/` set above additionally compares the full primary record, the alt haplotypes, the SNV
table, the VCF and the bubbles. The richer field list below describes
`verify_golden.py`, which moved to `../archive/` with the Python walker and is not runnable. The
`*.golden.json` files are therefore not read by anything in the current tree; they are kept as a
record of the original capture.

Every asserted field was verified byte-stable across two `PYTHONHASHSEED` values (1 and 7):

- graph shape (`n_unitigs`, `n_anchors`, `subgraph_nodes`)
- the primary haplotype, by SHA-256 of the sequence, plus `cov` and `pident`
- `n_bubbles`, and every bubble's fork node with its full branch list (node, `ka:f`, VAF)

The bubble `merge` node is deliberately **not** asserted. `phase_superbubbles` picks it with
`reconverged[0]` (in the original Python walker) from a `Counter` populated by iterating a `set`,
so it depends on string hash randomisation. Measured across seeds 1 and 7, it differed on
10,922 of 11,825 fork/merge pairs for `SRR8209082` and 18 of 202 for `ERR2028820`. The branch
calls themselves, which carry the biology, were identical in every case.

Timings are excluded. They are not reproducible and belong in benchmark output, not a fixture.

## Why these two accessions

They are the two named in the issue's concordance criterion. Note that neither appears in
`hiv1_50.txt` or `pfcrt_50.txt`, so the shipped benchmark harness does not exercise them.

`SRR8209082` is also worth treating with care as a performance baseline. Its 81,534 unitigs are
roughly 37x below the median of a 50-accession ERG11 run (2,986,072) measured during the
optimization work (that cohort list is not committed), and its 5-hop subgraph pulls in 26,323
nodes, about a third of its entire graph. It is fast and tangled, which
is not the common case.

## Regression coverage (verified)

The verifier was checked against two deliberate regressions:

1. A fixture with one minor allele removed. Caught, reporting `n_bubbles,
   quasispecies_branches`.
2. The `int`/`str` key migration hazard, injected by making `local_nodes` integer-keyed while
   `seqs` stays string-keyed. This empties the GFA and GraphAligner aborts with SIGABRT. The
   capture exits non-zero and the verifier reports a failure.

Hazard 2 is the one to keep in mind when bounding the subgraph index. In the full walker the
exception was caught per accession in the original Python walker and recorded as
`status: "error: ..."`, so the process still exits 0 and simply reconstructs nothing. With no
log redirect in `run_benchmark.sh`, that is easy to mistake for a working run.
