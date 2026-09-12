> **Historical — original design plan.** Written before the port; paths like `work/rust-walker/`
> refer to the upstream layout (here the crate is `walker/`). Kept as-is as the record of the
> plan. For current behaviour see [`../README.md`](../README.md) and [`rust_walker_as_built.md`](rust_walker_as_built.md).

# Plan: a targeted Rust rewrite of the path walker

Status: in progress. Five of seven stages implemented in `work/rust-walker/`; see section 12.
Addresses the five frontiers in [`issue_path_walker_optimization.md`](issue_path_walker_optimization.md).
All measurements below were taken on commit `81fc1fd` with GraphAligner 1.0.20 and BLAST 2.17.0+,
each run contained with `systemd-run --user --scope -p MemoryMax=<N>G -p MemorySwapMax=0`.

---

## 1. Summary

Replace the unitigs-to-haplotype span of `logan_path_walker.py` with a single Rust binary, keeping
GraphAligner as a subprocess and leaving the rest of the pipeline untouched.

The case for this is not the elimination of interprocess overhead, which measurement shows is
negligible. It is that the memory model becomes bounded by construction rather than by hope, that
the parse and scan stages are limited by interpreter speed rather than by algorithm, and that
fetch concurrency and compute concurrency can finally be sized independently.

One expectation needs setting up front. On a workstation this rewrite buys almost no throughput,
because the S3 download is already the limiting factor. The gain is realised when the walker runs
next to the data, in `us-east-1`. See section 10.

---

## 2. What the measurements actually say

Single accession, `ERR2028820`, 1,277,362 unitigs, `--workers 1`, total 14.6 s, peak RSS 897 MB:

| stage | seconds | share | frontier |
|---|---|---|---|
| `blast` | 5.8 | 40% | F2 |
| `subgraph` | 4.8 | 33% | F3 |
| `stream` | 3.8 | 26% | F1 |
| `walk` | 0.04 | 0.3% | F4 |

Three corrections to the issue's premises follow from this.

**BLAST is the largest stage, not subgraph extraction.** The issue implies the reverse.

**The "62.3% of runtime" figure for `extract_local_subgraph` in the committed ERG11 run is a thread
contention artifact.** That run used `--workers 16` against a `ThreadPoolExecutor`, so its per
accession numbers are contended wall times, not costs. Isolated cost is 3.76 microseconds per
unitig, and scaling is linear.

**The issue's stated baseline of 2 to 15 seconds per accession comes from unrepresentative
accessions.** `SRR8209082` measures 2.0 s and holds 81,534 unitigs, roughly 37x below the ERG11
median of 2,986,072.

Process spawn costs about 5 to 20 ms. Eliminating all three subprocesses would save on the order of
50 ms out of 14,600. What is expensive about BLAST is not the fork, it is that its interface
requires writing a 145 MB on-disk index to answer a question about 1,275 bp.

---

## 3. Scope

```
IN   accession list + reference CDS
     -> stream, anchor, subgraph, GFA emit, [GraphAligner], stitch, phase
     -> reconstructed.fasta + bubbles.json + summary.json

OUT  kmindex discovery, MAFFT, FastTree, HyphAeon, reporting  (all stay in Python)
```

The binary is a drop-in replacement for `logan_path_walker.py`: same CLI surface, same output files.
`run_panel_analysis.py` and `build_walker_codon_msa.py` do not change, which also means the new
walker can ship behind a flag and be A/B tested against the Python one on real cohorts.

---

## 4. Architecture

```
[ N async fetch tasks ]  --bounded channel-->  [ rayon compute pool ]
  reqwest + zstd                                 parse -> anchor -> BFS
  sized to bandwidth                             -> GFA -> align -> stitch -> phase
```

Stage by stage:

- **fetch**: streaming HTTPS with connection pooling, zstd decoded in memory.
- **parse**: byte slice tokeniser over the BCALM2 headers. No regex, no allocation per record.
- **anchor**: rolling canonical 31-mer scan against a hash set built from the query, folded into
  the same pass as the parse.
- **subgraph**: 5-hop BFS over a CSR adjacency structure.
- **emit**: GFA written to `/dev/shm`.
- **align**: GraphAligner subprocess, unchanged.
- **stitch**: path string to haplotype, skipping the k-1 overlap.
- **phase**: bubble detection, with the strand bug fixed (section 7).

### The two-pass problem, and why there is now only one pass

The BFS needs the full edge set before it can expand, and sequence extraction needs the surviving
nodes afterwards. Today this is solved by writing 145 MB to disk and reading it twice, which with
the record count scan in `stream_unitigs:62-67` makes three full traversals per accession.

*Resolved 2026-09-10, differently from the first draft of this section.* That draft weighed holding
the decompressed buffer (145 MB per worker) against retaining the compressed buffer (33 MB) and
decoding twice, and preferred the second. Working the numbers, **both peak at about 199 MB per
worker**: decoding twice merely shortens the time at peak, because the second pass rebuilds the same
145 MB. Neither helps the tail, where an accession decompresses to ~1.44 GB.

The design actually adopted is a **single streaming decode** in which `parse::build` extracts
adjacency, abundance and 2-bit packed sequences as the bytes go past, so the decompressed text is
never materialised at all. Sequences cost 17.2 MB packed against 145 MB as text, an 8.4x saving that
also removes the second decode. Peak per worker becomes roughly `compressed + graph`, about 76 MB.
Three file traversals become one.

---

## 5. Memory model

Unitig ids are dense and contiguous, and strictly sequential: `ERR2028820` runs 0 to 1,277,361 with
no gaps. Python cannot exploit this because it treats ids as opaque strings and reaches for dicts.
Rust indexes plain arrays directly, which removes every hash map from the hot path.

Measured shape of that accession: 2,582,206 `L:` tokens (2.02 per unitig), 68,745,081 bp of sequence,
mean unitig length 53.8 bp.

| structure | type | bytes |
|---|---|---|
| `offsets` | `u32[n+1]` | 5.1 MB |
| `targets` | `u32[m]`, 30-bit id plus 2 strand bits packed in | 10.3 MB |
| `abundance` | `f32[n]` | 5.1 MB |
| `seq_2bit` | 4 bases per byte | 17.2 MB |
| `seq_offset` | `u32[n+1]` | 5.1 MB |
| total | | **~43 MB** (40.9 MB measured) |

Against a measured 897 MB in Python for adjacency alone, and against 145 MB for the decompressed
text. The Python cost works out to about 700 bytes per unitig against roughly 114 bytes on disk, and
that constant is what makes a few concurrent large accessions reach the 21.3 GB that once froze a
workstation.

Note the id field is 30 bits, not the 21 an earlier draft proposed. 21 was sized from this
accession's 1.28M unitigs; `DRR029119` has ~13M and needs 24.

Holding sequences here is what lets `extract::collect_seqs` be O(nodes kept) rather than O(all
unitigs), and it removes the second pass over the file entirely. At ~43 MB per accession, 16
concurrent workers is about 690 MB. The out-of-memory failure stops being a risk to manage and
becomes a state the program cannot reach.

The strand bits packed into the spare bits of the target word are exactly the fields
`phase_superbubbles` currently discards, so the data layout fixes that bug for free.

---

## 6. Concurrency

Fetch and compute want different widths, and conflating them is part of why `--workers 16` was so
destructive. Measured S3 throughput saturates around 8 parallel streams at 17.4 MB/s aggregate,
while compute wants one task per core.

Two pools joined by a bounded channel. Capacity multiplied by maximum buffer size is a hard,
provable ceiling on resident memory. `ThreadPoolExecutor` never had one, which is precisely why a
single accession could grow without limit and take the machine with it.

Never move to process or task parallelism before the index is bounded. Threads sharing one heap is
what capped the original damage at 21.3 GB. Independent heaps at today's memory cost would be worse.

---

## 7. Frontier mapping

| frontier | verdict | notes |
|---|---|---|
| F1 in-memory streaming | do | real prize is collapsing three file traversals into one, not the disk write |
| F2 k-mer anchoring | do, validated | see section 8 |
| F3 bounded index | do, first | primarily a stability fix; 33% of runtime |
| F4 replace GraphAligner | **decline** | 0.3% of runtime |
| F5 true concurrency | do, last | only safe after F3 |

### Why F4 is declined

The walk is 0.04 s out of 14.6 s. The issue's premise that GraphAligner traverses "the compact
100-node DAG" is wrong twice over. Measured subgraphs are 356 nodes for `ERR2028820`, a median of 330
across the ERG11 panel, and 26,323 for `SRR8209082`, which is the issue's own HIV-1 benchmark. It is
also not a DAG but a bidirected graph with cycles, which is why real paths carry reverse complement
steps: the `ERR2028820` path has 68 forward and 47 reverse among its 115 nodes.

There is no mature sequence-to-graph aligner crate in Rust. `rust-bio` is pairwise, and the
`spoa`/`abpoa` bindings do partial order alignment over DAGs. This would be written from scratch,
and matching GraphAligner's sensitivity on divergent strains is the part that takes months.

Writing the 38 KB GFA to `/dev/shm` turns F4's "intermediate disk I/O" complaint into a RAM write and
captures essentially all of its available gain.

### The superbubble bug, fixed in passing

`phase_superbubbles` builds its adjacency while discarding the strand fields `s1` and `s2`
(`logan_path_walker.py:244-247`), so a bidirected edge recorded from both endpoints looks like a
reconvergence. Measured consequence: `fork == merge` on 93.1% of `ERR2028820`'s 202 bubbles and 58.5%
of `SRR8209082`'s 11,825. What is being detected is `A->B, A->C, B->A, C->A`, which is two edges read
backwards rather than a diamond.

Sequence this carefully. Port to byte-identical behaviour first and verify against the fixtures, then
fix the strand handling as a separate reviewable change with a deliberate fixture regeneration.

---

## 8. Validation

`work/benchmark/golden/` holds concordance fixtures for `ERR2028820` and `SRR8209082`, the two
accessions the issue names. `verify_golden.py` exits 0 when concordant. Asserted fields were verified
byte-stable across two `PYTHONHASHSEED` values; the bubble `merge` node is deliberately excluded
because it is hash-order dependent.

The verifier has been checked against two deliberate regressions, so it is known to be capable of
failing: a dropped minor allele, and the int/str key migration hazard.

**Watch for that second one during F3.** `seqs` is keyed by `str` at `:143` and `:151`, the frontier
test is `current_uid in local_nodes` at `:146`, and edges are filtered by `if tgt in seqs` at `:160`.
Convert the index to integers without converting that side and every lookup silently misses, the GFA
comes out empty, GraphAligner aborts with SIGABRT, the worker catches it at `:392`, and the run exits
0 having reconstructed nothing. With no log redirect in `run_benchmark.sh`, that reads as success.

### Anchor recall, measured

Frontier 2 was validated across 10 accessions, 3 organisms and 3 genes, by running the full
downstream on each anchor set and comparing haplotype SHA-256.

- **5 of the 6 reconstructing accessions produce byte-identical haplotypes** for a minimum shared
  k-mer count M in [1,5]. Above that it breaks, with two accessions shifting haplotype at M=10.
- **Anchor recall is a poor proxy for correctness.** Output was identical at recall as low as 0.375,
  because the 5-hop BFS bridges missing anchors. Gate on the reconstructed sequence, not on overlap.
- **An apparent recall collapse to 0.018 was an artifact of `-evalue 1e-10`.** A bare 31 bp exact
  match scores below that threshold, so BLAST never reports it. Relaxing to `1e-5` puts 91.5% of the
  k-mer anchors inside BLAST's set.
- **Do not build the consecutive-run refinement.** Requiring R consecutive matching k-mers was
  expected to track BLAST's 35 bp floor better than a plain count. Measured, the two agree in 22 of
  24 cells and the exceptions differ by one anchor, because a unitig is contiguous and its matching
  k-mers form one block. Use the plain count; it needs one counter instead of three.

### The one failure, and the guard

`SRR8209082` never reproduces its BLAST haplotype at any threshold. Its graph is saturated with the
target, 4,102 anchors out of 81,534 unitigs and a 26,323-node subgraph, so many near-equal paths
exist and the haplotype is under-determined. This is not a parameter to tune.

It is cheap to detect, but not the way this document originally proposed. **Corrected 2026-09-10
during implementation:** an "anchors exceed 1% of unitigs" rule does not work. The 5.0% figure for
`SRR8209082` is its BLAST anchor fraction, whereas the guard runs on k-mer anchors, which at the
default `min_shared_kmers = 3` are only 92 of 81,534 unitigs, or 0.113%. The threshold would sit 9x
above the value it must trip on and stay silent on the one accession it exists to catch.

The signal that works is the size of the 5-hop neighbourhood after expansion, as a fraction of the
whole graph. Measured across all 10 accessions at the default threshold, every accession whose k-mer
reconstruction matched BLAST sits at or below **0.033%**, while the single discordant one sits at
**10.383%**. That is a 315x gap with nothing in between, so a 1% threshold has an order of magnitude
of margin on both sides, and it measures the actual mechanism: a neighbourhood swallowing a tenth of
the graph offers many near-equal paths, which is what leaves the haplotype under-determined.

The guard therefore lives in `extract`, after expansion, not in `anchor`. Note it catches over-large
neighbourhoods only; the four `no_path` accessions sit at 0.002 to 0.004%, so too little signal is a
separate failure mode needing a separate check.

### Benchmark inputs need attention

Neither golden accession appears in `hiv1_50.txt` or `pfcrt_50.txt`, so the shipped harness never
exercised them. Both lists also carry duds: `DRR000252`, `DRR028850` and `DRR028861` yield 1 to 3
anchors and no path, and `DRR006432` and `DRR006471` also fail to reconstruct. Any throughput number
denominated in "50 accessions" is conflating speed with reconstruction yield. Curate the lists, or
report yield alongside throughput.

---

## 9. Milestones

1. **Port with no behaviour change.** Fetch, parse, CSR, BFS, GFA, GraphAligner subprocess, stitch,
   phase. Target: `verify_golden.py` passes. Keep BLAST in this step.
2. **Swap in k-mer anchoring** behind a flag, with the subgraph-fraction guard. Re-verify.
3. **Fix the strand bug** in phasing, regenerate fixtures deliberately, review the diff.
4. **Concurrency**, once memory is bounded: bounded channel, separate fetch and compute widths.
5. **Widen the validation panel** to more targets before trusting the guard threshold.

Rough effort for milestones 1 through 4 is one to two weeks for someone comfortable in Rust. The
anchor validation in milestone 5 is the item most likely to extend that.

---

## 10. Throughput expectations

Measured S3 bandwidth from a workstation: 6.9 MB/s on a single stream, 17.4 MB/s aggregate across 8
parallel range requests, which is where the link saturates. A median ERG11 accession is about 71 MB
compressed, so roughly 4.1 s of unavoidable transfer.

| | compute per accession | network per accession | bound |
|---|---|---|---|
| Python today | 14.6 s | 4.1 s | compute |
| Python with F2 and F3 | ~12 s | 4.1 s | compute |
| Rust | ~0.5 s | 4.1 s | **network** |

On this hardware both the fixed Python and the Rust version land near **0.245 accessions per second**,
because the uplink caps them. Sixteen cores would sit idle.

Next to the bucket in `us-east-1`, at roughly 1.2 GB/s, transfer falls to about 0.06 s per accession
and compute becomes the whole story. Rust reaches tens of accessions per second there. Python does
not reach the issue's target of 10 even with unlimited bandwidth, because its floor is set by
interpreter speed over 145 MB of text.

**The honest deliverable is therefore twofold: make compute cease to be the bottleneck, and state the
bandwidth requirement.** The 10 accessions per second target is an infrastructure decision as much as
a code one.

---

## 11. Open questions

- Does the subgraph-fraction guard hold across more targets? The 315x separation is clean on 10
  accessions, but all the discordant evidence is a single accession.
- Are dense contiguous unitig ids guaranteed by BCALM2, or only usual? Validate at parse time and
  fall back to a hash map rather than trusting it silently.
- Should the walker emit the full per-bubble detail? Today `phase_superbubbles` computes fork, merge,
  branch abundances and VAFs, and only `len(bubbles)` survives into the summary, so the quasispecies
  claim cannot be audited from the output.
- Is `SRR8209082`-class path ambiguity worth reporting explicitly, for instance by emitting the top
  N paths with scores rather than a single haplotype?

---

## 12. Implementation progress

Crate lives at `work/rust-walker/`. Run the fixed smoke batch with `cargo smoke` from that
directory. Every stage is cross-validated against the Python implementation's recorded numbers
rather than against itself.

| stage | file | status | measured |
|---|---|---|---|
| fetch | `stages/stream.rs` | done | pooled client, retry classification, size guard |
| decode | `stages/stream.rs` | done | 145 MB byte-identical to the `zstd` CLI |
| parse | `stages/parse.rs` | done | 0.61 s, graph 40.9 MB |
| anchor | `stages/anchor.rs` | done | 0.206 s vs BLAST 5.8 s, 28x |
| expand | `stages/extract.rs` | done | under 1 ms |
| collect_seqs + GFA | `stages/extract.rs` | not started | |
| align (GraphAligner) | `stages/align.rs` | not started | |
| stitch | `stages/stitch.rs` | not started | |
| phase | `stages/phase.rs` | not started | |
| output | `io/output.rs` | not started | |

Concordance with Python so far, all exact: canonical 31-mer counts (1,245 / 2,709 / 1,545 for
PfCRT / HIV-1 pol / ERG11); anchor counts on `ERR2028820` at five thresholds (104, 83, 77, 59, 36);
subgraph sizes at those thresholds (452, 422, 422, 386, 348), plus `SRR8209082` at 24,331 and 8,466.

### Structural facts verified, and enforced at runtime

- unitig ids are dense and **strictly sequential**, which is what makes CSR construction trivial
- every sequence is on a single line, though multi-line is handled anyway
- the edge relation is **symmetric**: 2,578,430 link tokens, zero missing reciprocals, so following
  CSR targets forward equals Python's explicitly symmetrised adjacency
- 1,156 self-loops exist and must not spin the BFS
- sequences are pure uppercase ACGT, so 2-bit packing is lossless (850 Mbp scanned, zero violations)

### Gotcha

`ureq`'s `Body::read_to_vec()` caps at 10 MiB by default. Every real accession exceeds it, so the
limit is raised explicitly. Left at the default this presents as a Logan or network fault.
