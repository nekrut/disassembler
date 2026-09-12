## Task: Optimization & Performance Acceleration of Multi-Accession Path Walking (`logan_path_walker.py`)

**Assignee:** Andrew (@as042)  
**Objective:** Scale multi-accession Logan path walking and orientation-aware quasispecies extraction from ~0.2–0.5 accessions/second to **$>10–50$ accessions/second**, eliminating I/O bottlenecks and subshell overhead across planetary cohorts ($N > 10{,}000$).

---

### 1. Context & Architectural Challenge

Path traversal through localized de Bruijn subgraphs is blazingly fast, yet end-to-end multi-accession throughput remains throttled by disk I/O, subshell orchestration, and sequential text parsing. Inside [`work/logan_path_walker.py`](work/logan_path_walker.py), GraphAligner navigates an extracted 5-hop local subgraph GFA in a mere **0.07 to 0.33 seconds per sample**. However, total wall time per accession currently spans **2 to 15 seconds**, dominated by:
1. Writing compressed AWS S3 unitigs to disk and invoking external `curl | zstd` pipelines.
2. Generating temporary BLAST nucleotide database indices on disk (`makeblastdb`) to identify matching anchor segments.
3. Parsing large unitig FASTA files line-by-line in pure Python to construct adjacency graphs from BCALM2 `L:` tokens.
4. Spawning external `GraphAligner -x dbg` subprocesses that read and write intermediate GFA and GAF text files to disk.

When scaling to planetary surveillance cohorts comprising 10,000 to 100,000 candidate runs (e.g. *Candida auris* *ERG11*, Avian Flu H5N1, or HIV-1 *pol*), processing accessions at 0.3/second requires days of continuous compute. We need to decouple path walking from intermediate disk I/O and process-spawning overhead to deliver true high-throughput parallel execution.

---

### 2. Environment Setup for You and Your Agent

To establish an identical, isolated testing and benchmarking sandbox:

```bash
# 1. Clone the repository and navigate to root
git clone https://github.com/nekrut/logan-hypyaeon.git
cd logan-hypyaeon

# 2. Configure Python virtualenv (uv recommended)
cd work
uv venv .venv
source .venv/bin/activate
uv pip install biopython requests duckdb pandas pyarrow zstandard

# 3. Configure bioinformatics binaries via micromamba / conda
micromamba create -n codon -c conda-forge -c bioconda graphaligner blast fasttree zstd curl
# Ensure codon binaries are accessible on PATH:
export PATH="$HOME/micromamba/envs/codon/bin:$PATH"
```

---

### 3. Standardized Test Datasets & Benchmark Harness

Pre-curated candidate accession lists and benchmark targets are available in the repository:

- **Candidate Accession Lists (in `work/benchmark/`):**
  - `work/benchmark/hiv1_50.txt`: 50 diverse HIV-1 *pol* candidate accessions.
  - `work/benchmark/pfcrt_50.txt`: 50 *Plasmodium falciparum* chloroquine resistance candidates.
  - `work/benchmark/sarscov2_50.txt`: 50 SARS-CoV-2 candidate accessions.
- **Reference In-Frame Coding Sequences (in `work/cds/individual/`):**
  - `work/cds/individual/candida_auris_erg11.fasta` (*C. auris* *ERG11*; 1,575 nt / 525 codons).
  - `work/viral/query_hiv1_pol.fasta` (HIV-1 Group M *pol*; 2,739 nt / 913 codons).
  - `work/cds/pfcrt_only.fasta` (*P. falciparum* PfCRT; 1,275 nt / 425 codons).

#### Executing the Standardized Baseline Benchmark:
We have provided an automated benchmark script [`work/benchmark/run_benchmark.sh`](work/benchmark/run_benchmark.sh):
```bash
./work/benchmark/run_benchmark.sh
```
This script measures wall-clock time, throughput (accessions/second), and successful reconstruction yield across 50 accessions. Use this harness to evaluate the speedup of each optimization iteration.

---

### 4. Five Optimization Frontiers to Explore with Your Agent

#### Frontier 1: In-Memory S3 Streaming & Decompression (Zero Intermediate Disk Writes)
- **Current Bottleneck:** `stream_unitigs()` launches subshell processes (`curl ... | zstd -dc`) to write uncompressed unitigs to disk (`{acc}.unitigs.fa`), reading them back immediately.
- **Opportunity:** Stream HTTP chunks directly into Python memory using `urllib3`, `requests`, or `httpx` with persistent connection pooling, decompressing in-memory via `zstandard.ZstdDecompressor().read()`. This eliminates all intermediate `.unitigs.fa` disk writes.

#### Frontier 2: Replacing `makeblastdb` with In-Memory $k$-mer Anchor Hashing
- **Current Bottleneck:** `find_anchors_blast()` invokes `makeblastdb` on disk, followed by `blastn -task megablast` to identify matching unitig IDs. This consumes 0.5–2.0 seconds per accession in process spin-up and index writing.
- **Opportunity:** The query CDS is known and short (1–3 kb; ~1,000 $k$-mers for $k=31$). An in-memory hash set of the query 31-mers can scan streamed unitig sequences on-the-fly. A unitig is flagged as an anchor if it shares $\ge M$ matching 31-mers with the query. This replaces BLAST entirely with a 5-millisecond pure-memory lookup.

#### Frontier 3: High-Performance Topological Subgraph Extraction
- **Current Bottleneck:** `extract_local_subgraph()` parses multi-megabyte unitig files line-by-line in pure Python to construct adjacency dictionaries and expand the 5-hop neighborhood.
- **Opportunity:** Implement the neighborhood BFS expansion in a compiled C/C++, Rust, or Cython extension, or optimize line parsing using fast C-level string methods. Compacting the 5-hop BFS will drop graph parsing from 1–2 seconds to $<50 \text{ ms}$.

#### Frontier 4: Memory-Piped GraphAligner Execution
- **Current Bottleneck:** Subgraphs are dumped to disk (`{acc}.sub.gfa`), GraphAligner is invoked via subprocess, and GAF alignments are read from disk (`{acc}.sub.gaf`).
- **Opportunity:** Explore passing the GFA directly to GraphAligner via stdin pipes (`/dev/stdin`) and capturing GAF output via stdout, or explore whether the compact 100-node DAG can be traversed directly using a fast dynamic programming / DAG shortest-path implementation in memory.

#### Frontier 5: True Multi-Process Concurrency (Bypassing Python GIL)
- **Current Bottleneck:** `ThreadPoolExecutor` suffers from Python Global Interpreter Lock (GIL) contention when multiple threads execute CPU-bound string parsing simultaneously.
- **Opportunity:** Transition the parallel worker orchestrator to `ProcessPoolExecutor`, Ray, or multiprocessing with shared memory buffers to achieve linear scaling across 16–64 CPU cores.

---

### 5. Success Criteria & Deliverables

1. **Throughput Target:** Exceed **$\ge 10 \text{ accessions/second}$** on a 16-core machine (a $\ge 20\times$ acceleration over the current baseline).
2. **Concordance Invariant:** 100% sequence identity concordance on reconstructed primary haplotypes against the golden-truth outputs for `ERR2028820` (PfCRT) and `SRR8209082` (HIV-1 *pol*).
3. **Quasispecies Preservation:** Retain orientation-aware superbubble parsing ($ka:f \ge 5.0$, $\text{VAF} \ge 10\%$) without dropping minor alleles.
4. **Pull Request:** Submit optimizations to `main` with before/after benchmark numbers emitted by `run_benchmark.sh`.
