# The Rust path walker, as built

Companion to [`plan_rust_path_walker.md`](plan_rust_path_walker.md), which is deliberately
left as written so it still represents the original plan. This document records what the
implementation actually turned out to be, including the places where building it contradicted
the plan.

Crate: `work/rust-walker/`, binary `logan-walker`. Written 2026-09-10 against commit `94ce920`.

**Status: the pipeline runs start to finish.** Every stage from accession list to codon MSA is
implemented and validated. It is not production-ready; see section 9.

> **Also read [`issue_intra_accession_variance_lost.md`](issue_intra_accession_variance_lost.md).**
> While porting the phasing stage it became clear that the pipeline detects real within-sample
> variants and then discards them, emitting one arbitrary haplotype per accession. The Rust
> port reproduces that behaviour deliberately, because the task was speed. That document
> explains what is being thrown away and what it would cost to keep.

---

## 1. What exists

```
src/
  main.rs        CLI, config, explicit IN -> pipeline -> OUT
  pipeline.rs    threading structure and stage wiring
  types.rs       what moves between stages
  config.rs      tunables
  error.rs       error type
  io/input.rs    accession list, reference CDS
  io/output.rs   reconstructed.fasta, bubbles.json, summary.json
  stages/stream.rs   S3 fetch + streaming zstd decode
  stages/parse.rs    BCALM2 records -> CSR graph + 2-bit sequences
  stages/anchor.rs   canonical k-mer anchoring (replaces BLAST)
  stages/extract.rs  5-hop BFS, sequence collection, GFA emit
  stages/align.rs    GraphAligner wrapper (external C++ binary)
  stages/stitch.rs   path -> haplotype
  stages/phase.rs    bubble detection and phasing
```

4,995 lines, 111 test functions, 3 dependencies (`ureq`, `zstd`, and `sha2` for tests only),
3.0 MB static binary, zero warnings.

Run the fixed smoke batch with `cargo smoke` from the crate directory.

---

## 2. Measured performance

`ERR2028820` against PfCRT, single accession, isolated, same machine:

| stage | Python | Rust | note |
|---|---|---|---|
| `stream` | 3.8 s | 4.8 s | network bound, not comparable run to run |
| `blast` -> `anchor` | 5.8 s | **0.203 s** | 28x, and no BLAST database written to disk |
| `subgraph` -> `parse` + `extract` | 4.8 s | **0.599 s** | 8x, and now also extracts every sequence |
| `walk` -> `align` | 0.04 s | 0.032 s | GraphAligner, unchanged |
| `stitch`, `phase` | included above | under 1 ms each | |
| **compute total** | **10.6 s** | **0.835 s** | **12.7x** |

Compute per accession drops from about 10.6 s to 0.835 s. The `stream` figure varies with the
link and is excluded from the ratio; on this hardware it dominates either way, which is the
point section 10 of the plan makes.

Memory, `ERR2028820`: `UnitigGraph` is **40.9 MB** resident, holding adjacency, abundance and
all 68,745,081 bases 2-bit packed, against a measured ~897 MB for the Python dictionaries.
Whole-process peak for a single accession is **84 MB**, and 235 MB for a four-worker run.

Four-accession ERG11 batch, 211 MB downloaded, 4 of 4 reconstructed: **13.7 s wall**, about
1.2 s of compute per accession.

### A measurement correction

Peak RSS figures quoted during development as 360 to 400 MB were wrong. They came from
`/usr/bin/time cargo run ...`, which measures rustc compiling, not the walker. Running the
binary directly gives the numbers above.

---

## 3. Bugs found and fixed

### 3.1 Strand-less adjacency in bubble phasing

`phase_superbubbles` builds its adjacency while discarding the `s1` and `s2` strand fields
(`logan_path_walker.py:244-247`). Two consequences follow.

A unitig has two ends, so a neighbour at each is an ordinary chain. Without orientation that
reads as a two-way fork. And because Logan writes every edge into both endpoints' headers,
`u -> v` and `v -> u` are the same edge seen twice, so both branches appear to lead back to
the fork, which is then reported as its own merge point.

Measured: `fork == merge` on **93.1%** of `ERR2028820`'s 202 bubbles and **58.5%** of
`SRR8209082`'s 11,825.

The Rust version works over oriented nodes, `(unitig, strand)`. An edge in `u`'s header is the
traversal `(u, s1) -> (t, s2)`, and the reverse traversal arrives independently from `t`'s own
header, so the graph is closed under reversal with each direction recorded once.

Effect on real data:

| accession | subgraph nodes | strand-less "forks" | oriented forks | real bubbles |
|---|---|---|---|---|
| `ERR2028820` | 422 | 254 | 232 | **1** |
| `SRR8209082` | 8,466 | 6,873 | 9,079 | **11** |

Note the fork counts are comparable. The collapse is not over-filtering: it is that the
reconvergence test, which the artifact had been satisfying trivially, became a real condition.

The 1.57M superbubble figure in the committed report is a sum of the uncorrected counts.

### 3.2 Hash-order dependent output

`c_node = reconverged[0]` (`:259`) selects from a `Counter` populated by iterating a `set`, so
the reported merge node depends on `PYTHONHASHSEED`. Measured across two seeds it differed on
10,922 of 11,825 fork/merge pairs for `SRR8209082`. The Rust version sorts and takes the
smallest, and the whole stage is order-independent.

### 3.3 Silent failure on an empty subgraph

An empty GFA makes GraphAligner abort with SIGABRT. In the Python walker the resulting
exception is caught per accession at `:392`, so the process exits 0 having reconstructed
nothing, which with no log redirect in `run_benchmark.sh` reads as a successful run.

Rust keeps the cases apart: `Ok(None)` is no alignment found, `Err` is the tool crashing, and
a nodeless subgraph short-circuits before spawning.

### 3.4 Scratch files leaked into RAM

Python unlinks its intermediates in a `finally` block. Early Rust drafts did not, leaving GFAs
in `/dev/shm`, which is RAM, so the growth would not appear in any process's RSS. `Subgraph`
now implements `Drop`, which runs on every exit path including errors. `LOGAN_KEEP_SCRATCH=1`
retains them for debugging.

---

## 4. Things only visible once it was built

### 4.1 The guard the plan proposed does not work

The plan recommended "if anchors exceed roughly 1% of unitigs, fall back to BLAST", citing a
500x separation. That figure is `SRR8209082`'s **BLAST** anchor fraction. The guard would run
on **k-mer** anchors, which at the default `min_shared_kmers = 3` are 92 of 81,534 unitigs, or
0.113%. The threshold sits 9x above the value it must trip on, and stays silent on the one
accession it exists to catch. Anchor fraction separates the good and bad cases by about 1.2x.

The signal that works is the neighbourhood size after expansion:

| accession | subgraph fraction | reconstruction |
|---|---|---|
| `SRR21578072` | 0.008% | matches BLAST |
| `SRR14253779` | 0.009% | matches BLAST |
| `SRR12363134` | 0.013% | matches BLAST |
| `ERR899743` | 0.019% | matches BLAST |
| `ERR2028820` | 0.033% | matches BLAST |
| `SRR8209082` | **10.383%** | **differs from BLAST** |

A 315x gap with nothing in between. The guard lives in `extract`, after expansion, and warns
above 1%.

### 4.2 Decoding twice saves nothing

The plan weighed holding the decompressed buffer (145 MB per worker) against retaining the
compressed buffer and decoding twice, preferring the second. Both peak at about 199 MB per
worker: decoding twice only shortens the time at peak, because the second pass rebuilds the
same 145 MB.

What was built instead is a single streaming decode in which `parse::build` extracts
adjacency, abundance and 2-bit packed sequences as the bytes go past, so the decompressed text
is never materialised. Sequences cost 17.2 MB packed against 145 MB as text. Peak per worker
becomes roughly `compressed + graph`, about 76 MB, and the three file traversals the Python
walker makes collapse to one.

### 4.3 Structural facts about Logan output

Each was checked against real files because the implementation depends on it, and each is
enforced at runtime so a violation is a loud error rather than a wrong graph.

* **Unitig ids are dense and strictly sequential**, 0 to n-1 with no gaps. This is what makes
  CSR construction trivial: a record's `L:` tokens are exactly the edges incident to that
  unitig, so appending in record order already produces CSR layout, with no edge sort and no
  temporary edge list. The plan had budgeted a ~20 MB transient for this.
* **The edge relation is symmetric**: 2,578,430 link tokens, zero missing reciprocals. So
  following CSR targets forward equals Python's explicitly symmetrised adjacency.
* **1,156 self-loops** exist and must not spin the BFS.
* **Sequences are pure uppercase ACGT**, verified over 850,835,755 bp across 10 accessions
  spanning bacterial, fungal, plant-pathogen, avian-flu and viral datasets. This is what makes
  2-bit packing lossless, and the absence of lowercase removes a case-folding branch from a
  loop that runs once per base.
* **Ids need 30 bits, not the 21 the plan assumed.** 21 was sized from `ERR2028820`'s 1.28M
  unitigs; `DRR029119` has about 13M and needs 24.

### 4.4 GAF columns 8 and 9 are path coordinates

The skeleton named them `query_start` and `query_end`. Reading Python's stitching closely
(`logan_path_walker.py:217-220`) they index the **stitched path**, not the query. Left
misnamed, `stitch` would have sliced with the wrong coordinates and produced a haplotype of
plausible length that was wrong, and only the hash comparison against the golden fixture would
have caught it.

### 4.5 `ureq` truncates bodies at 10 MiB

`Body::read_to_vec()` caps at 10 MiB by default. Every real accession exceeds it
(`ERR2028820` is 33 MB, the ERG11 panel reaches 343 MB). The limit is raised explicitly. Left
at the default this presents as a Logan or network fault.

### 4.6 Frontier 2 is not Python-impossible

An earlier assumption recorded during planning was that k-mer anchoring could not beat BLAST
in CPython. Measured in the same run it is roughly 2x faster on large accessions:
`ERR899743` 14.9 s BLAST against 6.9 s k-mer. The Rust gain is much larger, but the frontier
is not compiled-only.

---

## 5. Concordance

Every stage is validated against the Python implementation's own numbers, not against itself.

| check | result |
|---|---|
| zstd decode, 145 MB | byte-identical to the `zstd` CLI |
| every sequence round-trip | 1,277,362 sequences byte-identical to the original text |
| graph shape | 1,277,362 unitigs, 2,582,206 edges, 68,745,081 bases, all exact |
| canonical 31-mer counts | 1,245 / 2,709 / 1,545 for PfCRT / HIV-1 pol / ERG11, exact |
| anchor counts, 5 thresholds | 104, 83, 77, 59, 36, exact |
| subgraph sizes, 5 thresholds | 452, 422, 422, 386, 348, exact |
| GFA output | **byte-identical**, 44,859 bytes, same 77 anchors |
| GAF fields | identical: 115 nodes, 68 forward, 47 reverse, plen 1288, 1268/1274 |
| final haplotype | **SHA-256 identical** to `work/benchmark/golden/` |

The last line is the one that matters: fetch through stitch in Rust lands on the same 1,274
base haplotype Python produced, having replaced BLAST with k-mer anchoring and built a
different subgraph (422 nodes against 356) along the way.

---

## 6. Deliberate divergences from Python

1. **Bubbles are correct rather than compatible.** Section 3.1. The golden fixtures pin the
   pre-fix counts and need regenerating, with the diff reviewed.
2. **Phasing runs even when no path is found.** It reads the graph, not the GAF. A `no_path`
   accession reports bubbles here and reports none in Python.
3. **Per-bubble detail is written.** `bubbles.json` has no Python counterpart; there only
   `len(bubbles)` survives. See the linked variance document.
4. **Timing field names differ.** Python reports `stream`/`blast`/`subgraph`/`walk`, which do
   not map one to one onto seven stages. Everything else in the summary keeps Python's field
   names so the two files stay diffable.

---

## 7. A reproducibility problem that is not ours

The committed ERG11 panel cannot be reproduced exactly in this environment, independent of
the port. For `ERR899743`:

| run | anchors | subgraph | haplotype |
|---|---|---|---|
| Rust, k-mer M=3 | 159 | 533 | 1,574 nt |
| Python, k-mer M=3 | 159 | 533 | 1,574 nt, identical GFA and GAF |
| Python, BLAST, re-run here | 216 | 542 | 1,574 nt |
| committed run, BLAST | 216 | 542 | **1,575 nt** |

Rust and Python agree exactly. Re-running Python with BLAST here also gives 1,574. The
committed output differs by one trailing base, with the same anchors and the same subgraph
size, so the variable is the environment rather than the walker. The most likely cause is the
GraphAligner build, since the difference is where the alignment starts (`ps=1` against `ps=0`).
Our sequence is an exact prefix of the committed one.

This is why the golden fixtures were captured locally rather than taken from committed output.
`ERR2028820` matches its fixture byte for byte, including the hash, because both were produced
with the same GraphAligner.

---

## 8. End to end

Verified on four *C. auris* ERG11 accessions:

```
[1/4] SRR12363134: DONE in 1.1s (nodes=253, cov=1.00, pid=100.0%, bubbles=0)
[2/4] SRR14253779: DONE in 1.0s (nodes=150, cov=1.00, pid=100.0%, bubbles=0)
[3/4] SRR21578072: DONE in 1.2s (nodes=153, cov=1.00, pid=100.0%, bubbles=10)
[4/4] ERR899743:   DONE in 1.5s (nodes=533, cov=1.00, pid=99.9%,  bubbles=0)
4/4 accessions reconstructed, wall 13.7 s
```

Feeding `reconstructed.fasta` to the unmodified `build_walker_codon_msa.py`:

```
4/4 sequences passed coverage filter
MAFFT aligned 5 records
Codon alignment complete: 5 taxa x 1575 nt (525 codons)
```

All records in frame, zero internal stop codons, three of four haplotypes identical, which is
consistent with *C. auris* ERG11 being highly clonal. FastTree produced a phylogeny.

HyphAeon was not run: it needs a CUDA GPU and the private `veg/HyphAeon` repository.

---

## 9. Not done

* **`run_panel_analysis.py` still points at the Python walker by default.** Set
  `LOGAN_WALKER_BIN=work/rust-walker/target/release/logan-walker` to route step 1 through the
  Rust binary. The seam is in place and unset by default, so nothing changes without opting in.
* **Golden fixtures need regenerating** for the corrected phasing, as a separate reviewed
  change.
* **Concurrency is threads over a bounded channel**, which is correct and memory-safe but has
  not been tuned. `crossbeam-channel` would remove the `Arc<Mutex<Receiver>>` the standard
  library's single-consumer channel forces.
* **No BLAST fallback.** The guard warns; it does not yet fall back for accessions above the
  threshold.
* **The subgraph-fraction guard rests on one positive example.** The 315x separation is clean
  across 10 accessions, but all the discordant evidence is `SRR8209082`. Widen the panel before
  trusting the exact threshold.
* **Bubble reconvergence is one hop only**, matching Python. A variant whose alternative paths
  run several unitigs before meeting is missed. Pre-existing, not a regression.
* **The hardware ceiling still applies.** Measured S3 throughput here is 6.9 MB/s single stream
  and 17.4 MB/s saturated, so a median accession needs about 4.1 s of transfer. Compute is now
  well under that, meaning the walker is network bound on this machine and the remaining gains
  are only realised next to the data.
