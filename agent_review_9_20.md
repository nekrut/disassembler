# Code review findings, 2026-09-20

Outstanding problems found during a three-part review of this repository (Python and tests, the
Rust walker, documentation and hygiene). The documentation and hygiene items from that review
have already been fixed; what follows is everything still open.

Each finding gives a file and line, the defect, and a concrete failure scenario. Status is
**CONFIRMED** (reproduced, or traced precisely through the code) or **PLAUSIBLE** (reasoned from
the code but not demonstrated). Items marked *(verified directly)* were re-checked by hand after
the review rather than taken on the reviewer's word.

---

## 0. Cross-cutting: multi-record query FASTA is mishandled in three places

Three independent reviews landed on this. The repo ships `queries/master_planetary_targets.fasta`
with **160 records**, and every entry point below does something different and wrong with it.
Worth fixing as one coordinated change.

| location | behavior |
|---|---|
| `walker/src/stages/align.rs:33-46,144-162` | `io/input.rs:128-136` warns "using only the first" and returns record 1, but the *whole file* is handed to GraphAligner and the best-scoring row across all records wins, ignoring GAF column 1 (query name). `stitch` then slices with that row's coordinates and divides by record 1's length. |
| `python/disassemble.py:108` | sums residues across every record to derive the sieve's auto `--min-len`. |
| `python/lexicmap_streamer.py:150-172` | keeps record 1's header but concatenates every record's sequence into one pseudo-gene. |

**CONFIRMED.** Reviewer reproduced the walker case on `SRR8209082`: a 1-record query gives
`cov=0.997, pident=99.0, 299 nt`; adding a second record gives `cov=2.097, pident=94.0, 629 nt`,
the wrong record's path under the first record's name. A `cov` above 1.0 is impossible and would
make a cheap assertion.

*(verified directly)* `disassemble.py:108` yields `min-len 2573` for the single-record Zika query
but **297,645** for the 160-record master panel, which rejects every haplotype and then aborts
with an error blaming `--sieve-threshold`, the wrong lever.

**CONFIRMED** for the streamer: `load_reference` on the master panel returns `ref_len=313311`
under the name `h5n1_pb2_reference`. Every accession then lands in tier 4 at `cov` around 0.009
and the whole cohort is discarded, at exit 0.

---

## 1. Silent wrong results

### 1.1 Concurrent runs corrupt each other's scratch files
`walker/src/stages/extract.rs:264`, with `align.rs:99-131` and `pipeline.rs:340`.

The GFA and GAF paths are keyed on the accession alone, with no PID, worker id or run component,
and `--scratch-dir` defaults to `/dev/shm/logan-walker` for every process. `io/input.rs:19-53`
deliberately keeps duplicate accessions, and `Subgraph::drop` (`types.rs:196-203`) unlinks the
file on every exit path.

**CONFIRMED, reproduced.** Four copies of `SRR8209082` at `-w 4`, three consecutive runs:

```
run 1: 1 ERROR (GraphAligner "non supported overlap: 3", a half-written GFA), 2 SIGABRT, 1 DONE
run 2: 3 no_path, 1 SIGABRT   -> 0 records, exit status 0
run 3: 3 no_path, 1 DONE
```

Single-worker control: `DONE, cov=1.00`. Two mechanisms compound. One worker truncates or unlinks
a GFA another is aligning against, and `align.rs:126` maps a missing GAF (`ErrorKind::NotFound`)
to `Ok(None)`, so a file deleted by another worker becomes `Outcome::NoPath`, indistinguishable
from a real biological negative.

The likelier real trigger is two `logan-walker` processes on different genes side by side, since
they share the same default scratch dir. `io/output.rs:284` has the same problem for VCFs
(`<out_dir>/<acc>.vcf`), so duplicate accessions overwrite each other's output.

*(verified directly)* the path construction at `extract.rs:264`.

Fix: include a per-run unique component (PID, or a run UUID) in the scratch filenames, and
distinguish "GAF absent because alignment found nothing" from "GAF absent because the file
vanished".

### 1.2 The process exits 0 when every accession fails
`walker/src/main.rs:263-267, 308-312`.

Per-accession failures become `AccessionResult::failed` and never reach `real_main`'s error path;
`real_main` returns `Ok(())` regardless of `summary.success`.

**CONFIRMED, reproduced.** A list of two 404 accessions reports both as ERROR, prints
`0/2 accessions reconstructed`, writes a 0-byte `reconstructed.fasta`, and exits **0**. Any
wrapper or scheduler reads that as a successful shard. The same happens when GraphAligner is
missing entirely.

This is the failure mode `docs/rust_walker_as_built.md` section 3.3 claims the Rust rewrite keeps
apart. Suggested: exit non-zero when `reconstructed == 0`, or add an explicit
`--allow-empty` for callers that genuinely tolerate it.

### 1.3 Output record order is nondeterministic, and it silently decides what `--sieve` keeps
`walker/src/pipeline.rs:132-140`; consumed by `io/output.rs:65, 137, 232, 335`.

Results are collected in completion order, so `reconstructed.fasta`, `summary.json`,
`bubbles.json` and `snvs.tsv` are not byte-reproducible across runs.

*(verified directly)* on our own nine-accession Zika run:

```
input order : 963 966 967 969 973 979 984 992 993
output order: 979 992 973 966 984 993 969 963 967
```

The compounding problem, which neither reviewer noted: the sieve promotes the **first** record it
sees as the representative for a near-duplicate group. `SRR6251979_primary` survived the sieve and
is exactly the first line of `reconstructed.fasta`. So **which haplotype survives `--sieve` is
decided by S3 race outcomes rather than by biology**, and rerunning the identical command can
yield a different representative and a different tree.

Fix: sort results by input index before writing. That fixes reproducibility and this at once.

### 1.4 `mean_ka` is always 0.0 but published in every header
`walker/src/stages/stitch.rs:144` is the only writer (`mean_ka: 0.0`); printed at
`io/output.rs:70` in the primary FASTA header and `:189` in `summary.json`. `types.rs:248`
documents it as "Mean k-mer coverage (`ka:f`) across the unitig path".

**CONFIRMED** *(verified directly)*: `grep` finds no other assignment, and every header this repo
produces reads `ka=0.0`. Either compute it (`extract::path_depth` already does the equivalent) or
remove the field.

### 1.5 Alt haplotypes report the primary's `pident`
`walker/src/stages/variants.rs:245-250, 465-470`, via `copy_coords` (`:569`).

`num_matches` and `aln_len` are carried over unchanged, so every record in `alt_haplotypes.fasta`
reports the primary's identity, and after `prefer_majority` swaps, the primary's own `pident`
still describes the pre-correction alignment. **CONFIRMED.** The labels are stale, not merely
approximate.

### 1.6 Total prune empties the subgraph while logging "nothing to remove"
`walker/src/pipeline.rs:406-412`.

The condition is `dropped == 0 || kept_nodes.is_empty()`. In the second case `dropped == before`,
yet the code prints "nothing to remove" and assigns the empty node and sequence sets back.

**CONFIRMED by reading.** An accession whose gene sits entirely below `min_abund` still aligns and
stitches a primary from pass 1, then `variants::view` returns `None` (`variants.rs:85`) and
`find_alts` returns early (`variants.rs:313-318`): no majority correction, no alts, no VCF rows,
no SNV rows, accompanied by a log line asserting nothing was pruned.

### 1.7 Harvest reports success after fetching nothing
`python/harvest_kmindex.py:49-60`, merge loop `:99-126`, final print `:148`.

`fetch_shard` retries five times then returns `{}`, indistinguishable from a legitimately empty
shard. Nothing counts or reports failures.

**CONFIRMED by trace.** An expired API key makes every shard fetch fail; the run prints
`Harvest complete: 0 unique SRA runs across 0 targets`, writes a `UNION.accessions.txt` containing
a single newline and a `summary.json` of `{"_union": {"n": 0}}`, and exits 0. A downstream
`--diversity` or `disassemble.py -l` then runs on an empty cohort.

Secondary: the backoff at `:59` sleeps even after the final attempt, wasting 2.5 s per shard. At
10,000 shards and 32 workers that is roughly 40 minutes of pure sleeping on a path that produces
nothing.

---

## 2. Correctness in bubble and variant calling

### 2.1 Mirror-dedup can drop a genuine bubble in both directions
`walker/src/stages/phase.rs:141-151`.

`merge = min(reconverged)` is an arbitrary pick, and the dedup `if fork > (merge ^ 1) { continue }`
assumes the mirror direction selects the mirrored merge. When a fork has more than one
reconvergence point, the two directions pick different merges and both can be suppressed.

**CONFIRMED by faithful simulation** of the algorithm (oriented adjacency, `min` merge, `^1`
dedup). Minimal counterexample, 4 unitigs, edges `(u, s1, t, s2)`:

```
(0,-,2,+) (1,+,3,-) (2,-,1,+) (2,+,1,+) (3,-,0,+) (3,-,2,-) (3,-,2,+) (3,+,1,+)
```

Oriented node `3-` forks to `2-` and `2+`, which both reconverge at `1+`, a real diamond.
Direction `3-` picks `merge=1+` (id 3) and is skipped because `6 > 3^1 = 2`. Its mirror fork `1-`
picks `merge=0+` (id 1) and is skipped because `2 > 1^1 = 0`. The diamond appears in neither
direction of `bubbles.json`.

This is a simulation rather than an observed failure on real data, so reproduce against a real
graph before acting.

### 2.2 Branches that never reach the merge are reported as bubble alleles
`walker/src/stages/phase.rs:155-174`.

`total` and `kept` iterate `branches`, meaning every exit from the oriented fork, not only those
whose successor set contains the chosen `merge`. `variants::find_alts` does this correctly, proving
a rejoin per candidate with `rejoin()`, so there is a working model in-tree.

**CONFIRMED by reading.** A fork with exits b1, b2 (meeting at m) and b3 (going elsewhere at 400x)
emits `Bubble{fork, merge: m, branches: [b1, b2, b3]}`. b3 is published as an allele of a bubble it
is not part of, and b1/b2's VAFs are diluted by b3's coverage, which can push a real minority
allele under `min_vaf` and suppress the bubble entirely via `kept.len() < 2`.

---

## 3. Robustness and crashes

### 3.1 FastTree is never preflighted; failure lands after the MSA and leaves a 0-byte tree
`python/build_walker_codon_msa.py:106` and `:136-140`; `python/disassemble.py:99-100` preflights
only `mafft`.

`find_bin` returns the bare name `"fasttree"` when it is not found, and the `subprocess.run` is
unwrapped, so it raises `FileNotFoundError` instead of reaching the `sys.exit("FastTree failed: ...")`
path immediately below.

**CONFIRMED, reproduced.** With a stub mafft and no fasttree:

```
FileNotFoundError: [Errno 2] No such file or directory: 'fasttree'
a.fasta  5511 bytes
t.nwk       0 bytes      <- looks like a real output to anything downstream
```

Failure scenario: a 2,000-accession cohort in an env with mafft but no fasttree passes the
preflight, burns hours through the walk and MAFFT, then dies on an unhandled traceback and leaves
an empty `.nwk`. The same class applies to `mafft` if `find_bin` resolves it to a non-executable
path, since `find_bin` tests `candidate.exists()` rather than `os.access(..., X_OK)`.

### 3.2 A stray micromamba env silently replaces the pinned toolchain
`python/build_walker_codon_msa.py:25-35`.

```python
MAMBA_PREFIX = Path(os.environ.get("MAMBA_ROOT_PREFIX", str(Path.home() / "micromamba")))
CODON_ENV_BIN = MAMBA_PREFIX / "envs" / "codon" / "bin"
```

`find_bin` checks `CODON_ENV_BIN` **before** `shutil.which`. Flagged independently by two
reviewers.

**CONFIRMED, reproduced** with `MAMBA_ROOT_PREFIX` pointed at a temp prefix holding
`envs/codon/bin/mafft`: `which mafft` resolves to the `disassembler` env (what `disassemble.py`
preflighted) while `find_bin('mafft')` resolves to the codon env (what actually ran). On any
machine with a micromamba env named `codon`, the pipeline silently uses different MAFFT and
FastTree versions, which makes `environment.yml`'s claim that the pins match the golden fixtures
untrue.

Leftover from the upstream layout. Either drop the special case or make it opt-in and log which
binary was chosen.

### 3.3 `ZeroDivisionError` for `--diversity` below 5
`python/harvest_kmindex.py:151-156` (`step = len(items) / k` at `:155`), reached from `:195-197`.

`quota_canon = int(total_target * 0.20)` is 0 for any `--diversity` under 5, and
`stratified_sample` only guards `len(items) <= k`, which is false when the bracket is non-empty and
`k == 0`.

**CONFIRMED, reproduced.** With one scores file containing `{0.995, 0.991, 0.85, 0.6}`:
`--diversity 4` raises `ZeroDivisionError`; `--diversity 10` and `2500` are fine. Same crash for
`--diversity 1/2/3` whenever the `>= 0.99` bracket is non-empty, and for `1/2` whenever the
0.80 to 0.99 bracket is non-empty (`quota_high = int(n * 0.40)`).

### 3.4 Untrusted target names used as dict keys and filenames
`python/harvest_kmindex.py:119` and `:133-134`. **PLAUSIBLE.**

`target_name.split()[0]` raises `IndexError` on an empty or whitespace-only target name. The only
`try/except` in that loop wraps `fut.result()` (`:107-110`), so this escapes `as_completed`, tears
down the executor and kills a multi-hour harvest with **nothing written**, because all merging is
in memory and files are written only at `:130-146`.

`out_dir / f"{target}.accessions.txt"` uses the value verbatim, so a target containing `/` writes
outside `out_dir` or raises, and one containing `..` escapes it.

### 3.5 One malformed field aborts the stream after outputs are partly written
`python/lexicmap_streamer.py:499-503`, bare `float()` and `int()` on `pident`, `qcov`, `qstart`,
`qend`.

**CONFIRMED, reproduced** with a row whose `pident` is `NA`:
`ValueError: could not convert string to float: 'NA'`. The five cohort handles opened at `:353-357`
are closed by the `finally` at `:532`, but `.clean.accessions.fasta`, `.flagged.tsv` and
`.cohort_ledger.tsv` are left holding a truncated cohort while `.summary.json` is never written, so
the output directory looks half-plausible. Realistic triggers: concatenated LexicMap outputs with
an embedded second header row, `NA` or `-` in a numeric column, a `qseq` containing a tab.

### 3.6 Empty or header-only reference gives `ZeroDivisionError`
`python/lexicmap_streamer.py:424` (`cov = covered_bases / ref_len`). `load_reference`'s only guard
is `path.exists()`, and nothing checks `ref_len > 0`. **CONFIRMED** with a zero-byte reference. A
reference path that is a *directory* also passes `path.exists()` and raises `IsADirectoryError`
from `open()` at `:158` rather than the intended `FileNotFoundError`.

### 3.7 `-a` bypasses accession validation
`walker/src/main.rs:278` versus `io/input.rs:57-78`. The `--accession-list` path calls
`validate_accession`; the single-accession path does not. The value is interpolated into a URL and
two filesystem paths (`extract.rs:264`, `output.rs:284`), so `-a ../../foo` escapes both `out_dir`
and the scratch dir. Low severity, since it is the user's own command line, but the validator
already exists. *(verified directly: `validate_accession` is called only from `input.rs:32`.)*

---

## 4. Correctness gaps in the streamer

### 4.1 `qlen` and the `query` column are parsed but never used
`python/lexicmap_streamer.py:325-346`; `grep -n "qlen\|idx_query"` returns nothing further.

Two consequences. There is **no reference/query length agreement check**, though the TSV carries
`qlen` per row, so a reference FASTA of a different length or gene version mis-projects every HSP
onto the grid with no warning. And **multi-target TSVs are merged onto one reference**: a LexicMap
search against the 160-target master panel emits rows for all targets in one file, distinguished
only by the `query` column, but grouping is by `sgenome` alone (`:497`, `:513-519`), so one
accession hit by several genes has HSPs from *different genes* merged into a single grid. If the
file is sorted by query rather than accession it instead trips the `ValueError` at `:411`. Rows
should be filtered to `query == ref_name`.

### 4.2 Dead branch in mutation classification
`python/lexicmap_streamer.py:662-663`. Both `not aa_muts` cases are consumed by the first two arms,
so `aa_muts` is always truthy by the fourth and `"Complex"` can never be emitted.

Related: a haplotype differing from the reference only at positions it does not cover (all `-`)
produces empty `nt_muts` and `aa_muts` and is labelled `WT` in the TSV, while `summary.json`'s
`wildtype_count` (`clean_full_haps.get(ref_seq, 0)`) excludes it, so the two outputs disagree on
how many wild-type sequences there are.

---

## 5. Performance

### 5.1 Memory ceiling is understated roughly 5x
`walker/src/pipeline.rs:54-58` claims "Capacity times the largest retained buffer is a provable
ceiling on resident memory." It is not. Concurrently live compressed buffers are
`channel_capacity` (4, `pipeline.rs:90`) **plus** one per in-flight fetch worker (8, blocked in
`stream.rs:117-122`) **plus** one per compute worker, because `run_one` holds `fetched` through the
whole parse and drops it only at `pipeline.rs:284`. That is 20 buffers at the defaults, and
`max_compressed_bytes` is now 1 GB (`config.rs:322`).

**CONFIRMED by reading.** A list of 343 MB ERG11-scale accessions implies roughly 6.9 GB of
compressed buffers alone, plus 8 CSR graphs and `ureq`'s transient 2x growth during
`read_to_vec`. Worst case at the default cap is around 20 GB. This is precisely the OOM the
two-pool design exists to prevent, and the documented bound is wrong by the factor that matters.

### 5.2 Oriented adjacency rebuilt 3+ times per accession, including inside a loop
`walker/src/variants.rs:170-171` calls `view()` once per swap iteration (up to `max_alts * 2` = 64),
and `view()` (`variants.rs:83-103`) rebuilds `local` (a `HashMap` over all subgraph nodes) and `out`
(`2n` separately allocated `Vec`s) from scratch each time, although both depend only on
`graph` and `subgraph` and never on the evolving path. `find_alts` (`variants.rs:320-339`)
duplicates the same construction inline instead of calling `view()`, and `phase` (`phase.rs:90-110`)
builds it a third time.

**CONFIRMED by reading.** On `SRR8209082`'s 34,582-node subgraph that is roughly 69k `Vec`
allocations plus a sort per build. Hoist `out` and `local` out of the loop and have `find_alts`
reuse `view()`.

### 5.3 The VCF/SNV coordinate map runs a full O(n*m) Needleman-Wunsch twice per accession
`walker/src/io/output.rs:503-556` allocates a full `(n+1)*(m+1)` byte traceback matrix, and
`snv_sites` is called once from `write_vcfs` (`:279`) and again from `write_snv_table` (`:336`),
recomputing the identical alignment. The guard at `:377` permits 20,000 x 20,000, which is a
**400 MB** single allocation.

**CONFIRMED by reading.** With both `--emit-vcf` and `--emit-snv-table` on a 200-accession panel and
a 10 kb query that is 400 full 100M-cell DPs run serially after the pipeline. Compute the map once
per result and reuse it; a banded DP would make the guard unnecessary.

### 5.4 Every streamer output is written twice, and mutation calling runs twice
`python/lexicmap_streamer.py:605`, `:617`, `:624`. The `.clean.*` files and the undocumented
`.msa.fasta` / `.haplotypes.tsv` / `.a2m` compat aliases are byte-identical siblings. For the TSV
the duplication is not only I/O: the loop body re-runs `translate_nt(seq)` and the per-position
nucleotide and amino-acid diff for every haplotype on the second pass. On a cohort with ~10^5
unique 3 kb haplotypes that is a second full O(haplotypes x ref_len) translation pass and a second
multi-GB write for no new information.

Relatedly, the module docstring's "O(1) streaming memory" claim for the collapser is wrong:
`clean_full_haps`, `clean_full_pident_sum` and `clean_full_samples` (`:360-362`) are keyed by the
full merged sequence and grow with the number of *distinct* haplotypes, which on a divergent target
approaches the accession count.

### 5.5 `rejoin`'s BFS clones the whole trail per enqueued node
`walker/src/stages/variants.rs:662-686`, `let mut t = trail.clone(); t.push(next)` per edge, with
`trail` up to `alt_max_detour` (default 100) long. **PLAUSIBLE**, not measured. The BFS stops at the
first rejoin, normally 1 to 2 hops, so the common case is cheap; the expensive case is the
`no_rejoin` outcome the stage's own log counts, where the BFS exhausts the reachable set, up to
69k nodes on `SRR8209082`, each copying up to 100 `u32`. Parent-pointer BFS with a single
reconstruct at the end removes it.

---

## 6. Test and CI defects

The two sieve tests are now wired into CI, but they carry the defects below.

### 6.1 `test_sieve.sh` misreports "zero promoted" as "did not exit 0"
`tests/test_sieve.sh:81` and `:95`. `grep -c` exits 1 when the count is 0, and it is the last
command in `run_sieve`, so the function returns 1 even though the sieve exited 0. Line 85's
`|| fail "sieve run 1 did not exit 0"` fires, and the dedicated assertion at `:95`
(`run1 promoted $N1 < 1`) is unreachable.

**CONFIRMED, reproduced** with `MIN_LEN=99999` to force zero promotions: the script prints
`FAIL: sieve run 1 did not exit 0` when the sieve in fact exited 0 and promoted 0. A regression
that makes the sieve promote nothing would be diagnosed as a crash. The same path triggers if
`--emit-promoted` ever changes shard naming, since an unexpanded glob also yields count 0.
*(verified directly: `grep -c` on a file with no matches exits 1.)*

### 6.2 `test_full_length_only_requires_complete_coverage` does not test the flag
`python/test_lexicmap_streamer.py:158-172`. SRR010 covers 100% and SRR011 covers 50%, and the test
asserts `clean_full_length == 1`. With the flag **off** the threshold is `min_coverage = 0.80`,
which SRR010 also clears and SRR011 also misses, so the assertion holds either way.

**CONFIRMED** by running the fixture both ways: both give `clean_full_length=1`. Deleting
`full_cov = 1.0 if full_length_only else min_coverage` (`lexicmap_streamer.py:441`) entirely would
leave this test green. A discriminating fixture needs an accession in the 0.80 to 0.99 band, for
example `qend=11` giving cov 0.917.

### 6.3 `test_sieve_integration.py` reuses one record file across both runs
`tests/test_sieve_integration.py:141-148`. The off-path run's return code is discarded and it
reuses the same `work/stage2_input.txt` the on-path run wrote, so the
`else: check(False, "off-path: stage 2 was invoked")` branch at `:148` is unreachable, and a second
run that dies before stage 2 is misreported as "stage 2 got the wrong `--input`" rather than
"stage 2 never ran". Use a distinct record file per run.

### 6.4 `check_golden.sh` discards the walker's exit status and stderr
`bench/check_golden.sh:37-39` runs the walker with `>/dev/null 2>&1` and never inspects the return
code. If the Logan S3 fetch times out, `recon.fasta` is absent, `.got` is empty, and the gate
prints `FAIL ERR2028820 (got 0 nt vs golden 1476 nt)` and exits 1, reporting an infrastructure
outage as a walk regression with every diagnostic thrown away. Capture the rc and echo the log tail
on non-zero.

### 6.5 `tests/test_sieve.sh:63` bound is decoupled from the fixture
An unexpected fixture size only prints a note, while `MAX_PROMOTED=3` stays enforced regardless, so
growing the fixture silently changes what the bound means.

---

## 7. Smaller items

- `python/disassemble.py:140, 169, 195`: all three `subprocess.run(..., check=True)` calls surface
  tool failure as a raw `CalledProcessError` traceback, inconsistent with the `die()` style used
  for every other error in the file.
- `python/disassemble.py:125, 127, 131, 133`: the `is not None` guards are always true, since the
  defaults are non-`None` numbers. Harmless today only because `extra` is appended last and the
  walker's parser is last-wins on duplicates, which is what makes the docstring's
  `-- --prune-frac 0.0` example work.
- `python/disassemble.py:96`: `Path(walker).is_file()` accepts a non-executable file, giving an
  uncaught `PermissionError` at `:140` instead of the intended `die()`.
- `python/disassemble.py:108`: the file handle from `query.open()` is never closed, and the read
  happens unconditionally even when `--sieve` is off.
- `python/disassemble.py:152`: `--sieve --skip-tree` is a silent no-op. The design choice is
  reasonable and commented, but nothing tells the user and the `--sieve` help text does not mention
  it. A user asking for a thinned FASTA without a tree gets the unthinned one and exit 0.
- `python/build_walker_codon_msa.py:83, 101`: temp files are `out_aln.parent/.tmp_<input name>`, so
  two concurrent runs writing into one output directory clobber each other's MAFFT input, and the
  temps leak if MAFFT is interrupted.
- `python/build_walker_codon_msa.py:118-125`: `out_aln` is opened for writing before MAFFT runs, so
  a MAFFT failure replaces a previously good alignment with a truncated or empty file before
  `sys.exit`.
- `python/harvest_kmindex.py:49`: `fetch_shard`'s `api_key` parameter is unused, since auth comes
  from the session headers set at `:91`; the `requests.Session` is never closed.
- `walker/src/pipeline.rs:26` and `config.rs:236-237` say the compute pool is sized to the core
  count; `compute_workers` is a hardcoded 8 and never queries the machine.
- `assets/candida_auris_venn.png` (310 KB) and `.svg` (43 KB) are tracked but referenced by no
  markdown file. Added by commit `a3999eb` with no accompanying write-up. Either the figure lost
  its document or the document was never written. Left in place deliberately, since deleting
  committed assets should be a human call.

---

## 8. Checked and found correct

Recorded so nobody re-derives them.

**Rust:** 2-bit pack/unpack and mid-byte unitig boundaries; canonical k-mer rolling and the
`j+1>=k` window bound in both `index_query` and `find`; CSR construction and the density check;
`graph_k` versus `k` separation, with overlap correctly using `graph_k` everywhere; `stitch`'s
overlap trim, revcomp and `path_start`/`path_end` clamping; the oriented-edge encoding `2*i +
is_forward`, consistent across `phase`, `variants` and the GAF `>`/`<` parse; `shifted_window`'s
prefix and replaced-span arithmetic; `snv_sites`' triallelic coalescing and 1-based POS; the NW
traceback's boundary rows.

**Python:** the `argparse.REMAINDER` `--` handling in `disassemble.py:136`, where typo'd flags
still error with exit 2 rather than being swallowed; the sieve stale-shard cleanup and `rank*`
concat at `:162-174`; `find_internal_stops`' terminal-codon exclusion and `translate_nt`'s bounds;
`merge_accession_hsps`' query-gap coordinate advance and the non-gap-only `mean_pid` numerator;
`sanitize_id` collision handling.

**Build and CI:** `cargo clippy --all-targets --release -- -D warnings` exits 0; `cargo test
--release` gives 133 passed, 0 failed, 17 ignored; all pinned GitHub Actions versions exist; the
submodule is cleanly pinned at `39d90c0` and `make serial` produces `l1_pipeline_serial` as
documented.

---

## 9. Suggested order

1. **1.1** scratch-path collision and **1.2** exit code. Both silently produce wrong results that
   look like success, and both are small changes.
2. **1.3** sort results by input index. Fixes reproducibility and the nondeterministic sieve
   representative together.
3. **Section 0**, the multi-record FASTA handling, as one coordinated change across the three call
   sites. Add a `cov <= 1.0` assertion while there.
4. **3.1** and **3.2**, the two toolchain-resolution problems, which are cheap and currently make
   the environment pins unreliable.
5. **5.1** correct the documented memory bound, then **5.2** and **5.3** if throughput matters.
6. The test defects in section 6, which currently mask the regressions they exist to catch.
