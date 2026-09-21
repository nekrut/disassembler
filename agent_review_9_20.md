# Code review findings, 2026-09-20

Outstanding problems found during a three-part review of this repository (Python and tests, the
Rust walker, documentation and hygiene). The documentation and hygiene items from that review
have already been fixed; what follows is everything still open.

Each finding gives a file and line, the defect, and a concrete failure scenario. Status is
**CONFIRMED** (reproduced, or traced precisely through the code) or **PLAUSIBLE** (reasoned from
the code but not demonstrated). Items marked *(verified directly)* were re-checked by hand after
the review rather than taken on the reviewer's word.

**Re-verified 2026-09-21.** Every finding below was checked against the tree a second time. All
cited addresses resolve, the reproductions in sections 3 and 6 were re-run, and section 8's build
results were reproduced exactly. Two findings did not survive that pass as written and have been
corrected in place: **2.1**, whose counterexample was not a graph this stage can see, and **4.1**,
whose columns turn out not to be parsed at all. The defects both describe are real; only the
mechanism was wrong.

**Status 2026-09-21.** Everything in sections 0 through 4, 5.1, 5.3, 5.4, 6 and 7 is fixed, each
with its own "fix as applied" note and, where the defect was reachable from a test, a regression
test that fails when the fix is reverted.

**Deliberately not fixed: 5.2 and 5.5.** Both were measured first, and the numbers did not justify
the change; section 5 opens with the measurements and the reasoning. In short, the walker is I/O
bound by about an order of magnitude, and the accession whose subgraph makes those two expensive
is also the one that takes under two seconds to fetch. Also open by design: the unreferenced
`assets/candida_auris_venn.*`, which is a human call, not a code change.

Suite state at that point: `cargo clippy --all-targets --release -- -D warnings` exits 0,
`cargo test --release` gives **155 passed, 0 failed, 17 ignored** (133 before this work),
`python/test_lexicmap_streamer.py` **25 passed** (7 before), the new
`python/test_harvest_kmindex.py` **17 passed**, and both sieve tests pass against the real
binary.

---

## 0. Cross-cutting: multi-record query FASTA is mishandled in three places

> **FIXED 2026-09-21.** All three call sites now select record 1 consistently, and the walker acts
> on that selection instead of undoing it. The description and line numbers below record what was
> wrong and no longer describe the tree; see **Fix as applied** at the end of the section.

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

### Fix as applied, 2026-09-21

| file | change |
|---|---|
| `walker/src/types.rs` | `Query` gains `n_records`, so a stage can tell that the file held more than the record it was handed. |
| `walker/src/stages/align.rs` | `run` writes the selected record to `<scratch>/<acc>.sub.query.fa` and aligns that when the query is multi-record, so GraphAligner never sees the other genes and does not spend 160x the alignment work on them. `command` takes the FASTA path now, not the `Query`. |
| `walker/src/stages/align.rs` | `parse_gaf` takes the expected query name and filters on GAF column 1 before ranking by matches, so a higher-scoring row for another gene cannot win. |
| `walker/src/stages/stitch.rs` | Warns when `cov` exceeds 1.02, as a tripwire in case another path reintroduces this. |
| `python/disassemble.py` | `first_record_len()` replaces the whole-file sum. It reads record 1 only, closes its handle, and runs only when `--sieve` is on. |
| `python/lexicmap_streamer.py` | `load_reference` stops at the second header and warns, instead of concatenating every record. |

Single-record queries, which is every golden fixture and every well-formed run, take exactly the
path they took before: no temp file is written and GraphAligner is handed `query.path` untouched.

**Verified.** On `queries/master_planetary_targets.fasta` (160 records), `load_reference` returns
`h5n1_pb2_reference` at 2,280 nt, matching its own header, where it used to return 313,311; and the
sieve's auto `--min-len` is 2,166, where it used to be 297,645. `cargo clippy --all-targets
--release -- -D warnings` exits 0 and `cargo test --release` gives **135 passed, 0 failed, 17
ignored** (133 before; two new tests in `align.rs` cover the query-name filter and the
unknown-name fallback). The streamer suite is **8 passed** (7 before; one new test for a
multi-record reference).

**Deviation from the suggestion above.** The `cov` check is a warning, not an assertion. A genuine
alignment can exceed 1.0 slightly, since insertions in the path push `path_end` past
`path_start + qlen`, so a hard assert would fire on valid data. The threshold is 1.02, well under
the 2.097 this bug produced.

**Deliberately not done.** A multi-record query is still accepted rather than rejected.
`io/input.rs` has carried an explicit "use the first, warn" decision since the port, citing
tools-iuc #8388, and all three sites now honour it rather than contradict it. Making it a hard
error is a separate call, and would be the stronger fix if no caller legitimately passes a panel.

---

## 1. Silent wrong results

### 1.1 Concurrent runs corrupt each other's scratch files
`walker/src/stages/extract.rs:264`, with `align.rs:99-131` and `pipeline.rs:340`.

> **FIXED 2026-09-21.** Scratch names now carry a PID and a per-task counter. See
> **Fix as applied** below; the line numbers in this finding predate the change.

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

**Fix as applied, 2026-09-21.** `extract::scratch_stem` mints `<acc>.<pid>.<n>`, where `n` comes
from a process-wide atomic counter, and the GFA is written to `<stem>.sub.gfa`. The GAF
(`pipeline.rs`) and the single-record query FASTA (`align.rs`, added for section 0) are both
derived with `with_extension`, so they inherit the uniqueness; a test pins that, because a stem
full of dots would otherwise let `with_extension` eat the wrong suffix. The PID separates
concurrent processes, which is the likelier trigger given `--scratch-dir` defaults to
`/dev/shm/logan-walker` for every process, and the counter separates duplicate accessions inside
one process.

The second half needed no code. With names unique per task, a missing GAF can no longer be another
worker's doing, so `align.rs`'s `NotFound => Ok(None)` is now sound rather than a conflation; the
comment there says so and names what it depends on.

`io/output.rs:284` is left writing `<out_dir>/<acc>.vcf`. That name is documented in two README
rows and in `docs/vcf_example.md:15`, and the committed examples use it, so renaming it would
break more than it fixes. Duplicate accessions produce results for the same accession, so the
overwrite is benign; `write_vcfs` now warns when it sees one rather than overwriting silently.

**Verified** on a real duplicate at `-w 2`: `SRR8209082` twice against `query_hiv1_pol.fasta`
wrote `SRR8209082.11257.0.sub.gfa` and `SRR8209082.11257.1.sub.gfa`, two files where the old code
wrote one path twice. Both fetched from S3 and extracted a 34,582-node subgraph without touching
each other.

### 1.2 The process exits 0 when every accession fails
`walker/src/main.rs:263-267, 308-312`.

> **FIXED 2026-09-21.** A run that reconstructs nothing now exits non-zero unless
> `--allow-empty` is passed.

Per-accession failures become `AccessionResult::failed` and never reach `real_main`'s error path;
`real_main` returns `Ok(())` regardless of `summary.success`.

**CONFIRMED, reproduced.** A list of two 404 accessions reports both as ERROR, prints
`0/2 accessions reconstructed`, writes a 0-byte `reconstructed.fasta`, and exits **0**. Any
wrapper or scheduler reads that as a successful shard. The same happens when GraphAligner is
missing entirely.

`docs/rust_walker_as_built.md:122-129` names this Python behaviour, exiting 0 having reconstructed
nothing, as one the rewrite improves on. What section 3.3 actually claims to fix is the
per-accession `Ok(None)` versus `Err` distinction, and the walker does do that, so this finding is
an adjacent gap rather than a broken promise. The claim section 3.3 does break is **1.1**, where a
GAF deleted by another process becomes `Ok(None)` and is reported as a biological negative.

Suggested: exit non-zero when `reconstructed == 0`, or add an explicit `--allow-empty` for callers
that genuinely tolerate it.

**Fix as applied, 2026-09-21.** Both. `real_main` returns a new
`Error::NothingReconstructed { total }` when `success == 0 && total > 0`, and `--allow-empty` opts
out for callers that genuinely tolerate an empty shard. The error names the count and points at the
per-accession log rather than reading as a usage mistake. Outputs are still written before the
check, so an empty `reconstructed.fasta` and the summary stay on disk for inspection; only the exit
status changes. A partially successful run is untouched, since one bad accession must not sink a
10,000-accession shard.

**Verified**, reproducing the finding's own case: two 404 accessions print
`0/2 accessions reconstructed` and now exit **1**, and exit **0** with `--allow-empty`. The
"GraphAligner missing entirely" variant the finding mentions also exits 1 now, confirmed on a real
`SRR8209082` run that reached stage 5 with no GraphAligner installed.

### 1.3 Output record order is nondeterministic, and it silently decides what `--sieve` keeps
`walker/src/pipeline.rs:132-140`; consumed by `io/output.rs:65, 137, 232, 335`.

> **FIXED 2026-09-21.** Results are sorted by input index before they leave `pipeline::run`.

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

**Fix as applied, 2026-09-21.** The input position is carried end to end rather than reconstructed
at the end, because duplicate accessions make "position of this accession in the list" ambiguous.
The work queue sends `(index, accession)`; `FetchedUnitigs` and `AccessionResult` each carry an
`index`; and `pipeline::run` sorts on it before building the `Summary`. Every failure path had to
be covered too, since a result can be born in the fetch worker (fetch failed), in `run_one` (no
anchors, or a stage error), or in the compute worker's `catch_unwind` (a stage panicked), and a run
with some failures would otherwise come back partly ordered. A test pins that.

Progress reporting is deliberately left in completion order: printing each accession as it lands is
what makes a stalled run visible, and it is what the Python walker did. Only the *stored* order
changed.

**Verified** on six accessions across two runs. Completion order differed between them
(`96 91 92 93 95 97`, then `96 95 91 97 92 ...`), and both runs stored `95 91 97 93 92 96`, which
is the input order exactly. All four outputs are byte-identical between the two runs:
`reconstructed.fasta`, `bubbles.json`, `alt_haplotypes.fasta` and the summary JSON.

Since the sieve promotes the first record it sees, its representative is now chosen by the list the
user supplied rather than by a race. That makes the choice *deterministic*, not *principled*: it is
still "whichever haplotype came from the first accession in the list", which is worth knowing given
the sieve's measured tendency to collapse a cohort to one survivor.

### 1.4 `mean_ka` is always 0.0 but published in every header
`walker/src/stages/stitch.rs:144` is the only writer (`mean_ka: 0.0`); printed at
`io/output.rs:70` in the primary FASTA header and `:189` in `summary.json`. `types.rs:248`
documents it as "Mean k-mer coverage (`ka:f`) across the unitig path".

**CONFIRMED** *(verified directly)*: `grep` finds no other assignment, and every header this repo
produces reads `ka=0.0`. Either compute it (`extract::path_depth` already does the equivalent) or
remove the field.

**FIXED 2026-09-21: computed.** `extract::path_mean_ka` is new alongside `path_depth`, and
`stitch` fills the field from it. `stitch` gained a `&UnitigGraph` parameter, since abundance
lives there and not in the `Subgraph`; threading it in rather than annotating afterwards means a
new call site cannot forget.

Deliberately the arithmetic **mean**, not `path_depth`'s median, even though the review points at
`path_depth` as the equivalent: the field is named `mean_ka` and is published as `ka=`, so it has
to be what it says. The median remains the right statistic for *deciding* a prune floor, where one
repeat-collapsed unitig must not drag the estimate, and `path_depth` stays exactly that. Both are
documented in terms of each other. Unweighted per step, so a repeated unitig counts twice.

**Verified** by four assertions across two tests: the mean of a two-node path, that it tracks the
graph rather than any fixed value, that a repeated step counts twice, and that on a 1/2/300 path
`mean_ka` is 101.0 while `path_depth` stays 2.0, so the two statistics cannot silently converge.

### 1.5 Alt haplotypes report the primary's `pident`
`walker/src/stages/variants.rs:245-250, 465-470`, via `copy_coords` (`:569`).

`num_matches` and `aln_len` are carried over unchanged, so every record in `alt_haplotypes.fasta`
reports the primary's identity, and after `prefer_majority` swaps, the primary's own `pident`
still describes the pre-correction alignment. **CONFIRMED.** The labels are stale, not merely
approximate.

**FIXED 2026-09-21: labelled, not recomputed.** A spliced path has no GAF of its own, and
measuring identity properly needs an alignment per haplotype: one Needleman-Wunsch per alt, up to
`max_alts` = 32 per accession, which is the cost **5.3** already objects to paying twice. Rather
than fabricate a number or add that cost, the published key now says which number it is.

`Haplotype` gains `pident_inherited`. `copy_coords` still carries the coordinates forward, because
they are needed, but the haplotype built from a spliced path is marked: `pipeline` sets it where
`prefer_majority` swaps the route, and `write_alts` labels every alt unconditionally, since an alt
is always spliced. The FASTA key is `pident=` only when the number was measured on that sequence's
own alignment, and `inherited_pident=` otherwise; `summary.json` gains a `pident_inherited` flag
next to `pident`. A consumer reading `pident=` therefore never gets a number about a different
sequence, which is the defect.

The prune-and-realign path is unaffected: it produces a real second GAF, so its numbers are
measured and it keeps a plain `pident=`.

**Verified** by two tests asserting the key that appears and, just as importantly, that a bare
`pident=` never appears on a spliced record. No doc or in-repo consumer parses these headers, and
the golden gate strips them, so the rename breaks nothing.

### 1.6 Total prune empties the subgraph while logging "nothing to remove"
`walker/src/pipeline.rs:406-412`.

The condition is `dropped == 0 || kept_nodes.is_empty()`. In the second case `dropped == before`,
yet the code prints "nothing to remove" and assigns the empty node and sequence sets back.

**CONFIRMED by reading.** An accession whose gene sits entirely below `min_abund` still aligns and
stitches a primary from pass 1, then `variants::view` returns `None` (`variants.rs:85`) and
`find_alts` returns early (`variants.rs:313-318`): no majority correction, no alts, no VCF rows,
no SNV rows, accompanied by a log line asserting nothing was pruned.

**FIXED 2026-09-21.** The survivor count is taken **before** anything is consumed, because
`prune_nodes` takes ownership of the node and sequence vectors and a prune that turns out to be
unusable could not otherwise be undone. The three cases are now a `PruneOutcome` enum rather than
a two-armed `if`:

* `WouldEmpty` (nothing clears the floor): skip the prune and keep the unpruned subgraph. This is
  the same "never a regression" rule the failed-realign branch below already applies, and it is
  what the old code should have done instead of assigning the empty vectors back.
* `NothingToRemove` (everything clears it): unchanged.
* `Prune`: unchanged.

Each logs a different line, so "would remove all N unitigs; skipping" can no longer be read as
"nothing to remove".

**Verified** by two tests: one pinning the three cases apart, including the `(0, 0)` and `(1, 1)`
edges; one showing the case is reachable on real data, with a three-unitig gene at 1-3x under the
5.0 `min_abund` floor, where the counting predicate and `prune_nodes` agree the subgraph would be
emptied.

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

**FIXED 2026-09-21.** `fetch_shard` returns `None` on exhausted retries and `{}` only for a shard
that genuinely matched nothing, so the two are distinguishable; it logs the last error with the
shard name, and no longer sleeps after the final attempt. The merge loop counts failures and
empties separately, the summary line reports both and switches its marker from `✓` to `!` when
anything failed, and a run in which **every** shard failed raises rather than returning, on the
same reasoning as **1.2**: a stage that produced no cohort must not look like one that produced an
empty cohort.

A partial failure still returns what was harvested, with a loud `the cohort below is INCOMPLETE`
warning, since discarding good shards would be worse than an incomplete cohort the operator is
told about.

**Verified** by a new `python/test_harvest_kmindex.py`, 8 tests: a failed fetch is `None` and an
empty shard is `{}`; a good response is not retried; exactly `ATTEMPTS - 1` sleeps totalling 5.0 s
rather than 7.5 s; a late success still returns data; all-failed raises naming the count;
all-empty is **not** an error and still writes `summary.json`; and a partial failure returns the
shards that worked. Wired into CI, skipping cleanly where `requests` is absent, which it is on the
current runner.

---

## 2. Correctness in bubble and variant calling

### 2.1 A genuine bubble can be dropped in both directions
`walker/src/stages/phase.rs:142-151`. **Rewritten 2026-09-21**; the original counterexample was
not reversal-closed and did not hold up, see below. **Mechanism corrected again while fixing it**,
see "one cause, not two".

**Only one reconvergence point per fork is ever considered.** `reconverged.sort_unstable()` then
`let merge = reconverged[0]` (`:141-142`) keeps the smallest and discards the rest, so a fork with
several reconvergence points reports at most one bubble no matter which direction it is seen from.

**The dedup then compares forks across different merges.** `if fork > (merge ^ 1) { continue }`
(`:149`) is sound when both directions are considering the same merge. Because of the pick above
they need not be, so both copies can be suppressed and the bubble appears in neither direction of
`bubbles.json`.

#### One cause, not two

An earlier version of this section (written 2026-09-21, before the fix) called these two
independent defects. They are not. **The dedup was never wrong on its own.** Once every
reconvergence point is enumerated, `f <= m^1` and `m^1 <= f` cannot both fail, so the fork-only
test emits each bubble exactly once. Verified two ways: algebraically, and by running both dedup
forms over 300,000 random reversal-closed graphs with zero disagreements. Caught by mutation
testing the fix, where reverting the dedup alone killed no test.

So the single root cause is `reconverged[0]`. The dedup only looked defective because that pick
left the two directions comparing different merges.

**CONFIRMED by simulation** of the algorithm (oriented adjacency, `min` merge, `^1` dedup), run
over edge sets **closed under reversal**. The closure is what makes a simulation meaningful:
`phase.rs:69-72` records that Logan writes every edge into both endpoints' headers, so
`(u,s1) -> (t,s2)` and `(t,!s2) -> (u,!s1)` both arrive independently and `out` is always closed
under reversal. An edge set that is not closed is not a graph this stage can be handed.

Counterexample, 4 unitigs, oriented edges, reversal-closed:

```
0-->1-  0-->1+  0-->2+  0-->3-  0+->2-  0+->2+  1-->0+  1-->3-
1+->0+  1+->2+  2-->0-  2-->0+  2-->1-  2+->0-  3+->0+  3+->1+
```

Oriented node `2-` exits to `0-`, `0+` and `1-`, and reconverges at both `2+` and `3-`, so
`(2-, 3-)` is a real diamond. It is then lost once to each mechanism. Forward, `2-` keeps
`merge = 2+` as the smaller candidate, emits that bubble, and never reports `3-`. In the mirror
direction `3+` reconverges only at `2+`, picks `merge = 2+`, and is suppressed by the dedup
(`7 > 5^1 = 4`). The algorithm emits `(0-,0+)`, `(0+,0-)` and `(2-,2+)`; `(2-,3-)` appears in
neither direction.

Found by brute-force search over random reversal-closed graphs (n = 3 to 6, no self-loops) for a
genuine diamond emitted in neither direction. Still a simulation rather than an observed failure
on real data, so reproduce against a real graph before acting.

**FIXED 2026-09-21.** The loop now iterates **every** reconvergence point rather than
`reconverged[0]`, so a fork contributes one bubble per merge. `bubbles.sort_by_key` gained the
merge, since a fork can now contribute several and the order would otherwise depend on discovery.

The dedup is now a comparison of the whole `(fork, merge)` pair against its mirror. That is
equivalent to the old fork-only test given the fix above, as established in "one cause, not two";
it is kept because it does not depend on that invariant, so a future change that again reports one
merge per fork cannot silently turn the dedup into a bubble-loser.

**Verified** by two tests plus mutation testing. `a_fork_with_two_merges_reports_both` uses a fork
with two genuine diamonds and asserts both merges are reported.
`a_bubble_whose_directions_disagree_on_the_merge_survives` is the counterexample above, written as
a BCALM2 fixture. Both fixtures assert their own reversal closure first, via a helper that walks
the parsed edges and checks each mirror is present, because an unclosed fixture is not a graph this
stage can be handed, which is what invalidated the original counterexample. Reverting the
`reconverged[0]` pick fails both tests.

### 2.2 Branches that never reach the merge are reported as bubble alleles
`walker/src/stages/phase.rs:155-174`.

`total` and `kept` iterate `branches`, meaning every exit from the oriented fork, not only those
whose successor set contains the chosen `merge`. `variants::find_alts` does this correctly, proving
a rejoin per candidate with `rejoin()`, so there is a working model in-tree.

**CONFIRMED by reading.** A fork with exits b1, b2 (meeting at m) and b3 (going elsewhere at 400x)
emits `Bubble{fork, merge: m, branches: [b1, b2, b3]}`. b3 is published as an allele of a bubble it
is not part of, and b1/b2's VAFs are diluted by b3's coverage, which can push a real minority
allele under `min_vaf` and suppress the bubble entirely via `kept.len() < 2`.

**FIXED 2026-09-21.** `total` and `kept` now iterate only the exits whose successor set contains
the chosen merge, found by binary search since `out[b]` is already sorted and deduped. The existing
ordering is preserved: a branch removed by `min_abund` or `min_vaf` still contributed to `total`,
which is what matches the Python behaviour; only branches that are not part of the bubble at all
are excluded.

**Verified** by `a_branch_that_misses_the_merge_is_not_an_allele`, built as the finding's own
scenario: exits at 80x and 20x meeting at the merge, plus a 400x exit going elsewhere. The bubble
reports exactly the two real alleles at VAF 0.8 and 0.2. Counting the third would put the minority
allele at 20/500 = 0.04, under the 0.10 floor, and `kept.len() < 2` would then drop the bubble
entirely, so the fixture exercises the suppression path and not just the labelling. Reverting the
filter fails it.

---

## 3. Robustness and crashes

### 3.1 FastTree is never preflighted; failure lands after the MSA and leaves a 0-byte tree
`python/build_walker_codon_msa.py:106` and `:136-140`; `python/disassemble.py:99-100` preflights
only `mafft`.

> **FIXED 2026-09-21.** Both tools are resolved before any work, and the tree is written to a
> temp and renamed on success.

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

**Fix as applied, 2026-09-21.** `find_bin` raises `FileNotFoundError` naming the tool instead of
returning the bare name, and `build_walker_codon_msa.main` resolves **both** binaries immediately
after argument parsing, before reading a single FASTA. `disassemble.py` preflights both too. The
executability problem is gone with it: resolution is `shutil.which` throughout, which tests
`X_OK`, rather than `Path.exists()`.

The tree is now written to `.tmp_<name>.<pid>` and `os.replace`d into position only after FastTree
exits 0, so a failure leaves no 0-byte `.nwk` and does not destroy a previous good tree. The
alignment got the same treatment, which closes the section 7 item at `:118-125`.

Both spellings of FastTree are accepted (`fasttree`, then `FastTree`), because bioconda ships it
under either depending on the build and resolving early would otherwise turn a wrong guess into a
blocked run. `disassemble.py`'s preflight uses the same list, so it cannot reject a run the MSA
step would have accepted.

**Verified.** With a stub mafft and no fasttree, `build_walker_codon_msa` exits
`Error: 'fasttree' or 'FastTree' not found on PATH or in ...` and writes **neither** output file,
where before it produced a 0-byte `t.nwk`. `disassemble.py` exits 1 at preflight with
`fasttree not on PATH (needed for the ML tree)`. With only a capital-F `FastTree` on PATH, both
resolve it and proceed.

### 3.2 A stray micromamba env silently replaces the pinned toolchain
`python/build_walker_codon_msa.py:25-35`.

> **FIXED 2026-09-21.** PATH now wins; the `codon` env is a fallback that announces itself.

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

**Fix as applied, 2026-09-21.** The order is inverted: `shutil.which(name)` is tried first, and
`CODON_ENV_BIN` only if the tool is not on PATH at all. So `environment.yml`'s pins describe what
actually runs whenever the env is active, and `find_bin` agrees with the `shutil.which` that
`disassemble.py` preflights.

The special case is kept rather than dropped, and logs when it fires. It is not dead weight:
`docs/issue_path_walker_optimization.md:40` documents creating exactly that env
(`micromamba create -n codon ... fasttree ...`), so dropping it would break a setup this repo's
own docs tell you to build. Keeping it as a fallback means such a machine still works, and the
stderr line says which binary was chosen and that the pins do not describe it.

**Verified**, reproducing the finding's setup. With `MAMBA_ROOT_PREFIX` pointed at a prefix
holding `envs/codon/bin/mafft` *and* a different mafft on PATH, `find_bin('mafft')` and
`shutil.which('mafft')` now return the same path, where before they disagreed. With mafft removed
from PATH, the fallback still resolves and prints
`[*] mafft: not on PATH, falling back to ... (environment.yml pins do not describe this binary)`.

### 3.3 `ZeroDivisionError` for `--diversity` below 5
`python/harvest_kmindex.py:151-156` (`step = len(items) / k` at `:155`), reached from `:195-197`.

`quota_canon = int(total_target * 0.20)` is 0 for any `--diversity` under 5, and
`stratified_sample` only guards `len(items) <= k`, which is false when the bracket is non-empty and
`k == 0`.

**CONFIRMED, reproduced.** With one scores file containing `{0.995, 0.991, 0.85, 0.6}`:
`--diversity 4` raises `ZeroDivisionError`; `--diversity 10` and `2500` are fine. Same crash for
`--diversity 1/2/3` whenever the `>= 0.99` bracket is non-empty, and for `1/2` whenever the
0.80 to 0.99 bracket is non-empty (`quota_high = int(n * 0.40)`).

**FIXED 2026-09-21.** `stratified_sample` returns `[]` for a quota of zero or less. A zero quota
is legitimate rather than an error: contributing nothing from that bracket is the right answer at
that size, and `stratify_diversity`'s existing top-up (`if len(selected) < total_target`) then
fills the cohort from the pooled remainder, so the requested count is still produced.

**Verified** on the finding's own fixture: every `--diversity` from 1 to 4 now writes a cohort of
`min(k, 4)` accessions, and 10 and 2500 are unchanged. Mutation-tested: removing the guard fails
five tests.

### 3.4 Untrusted target names used as dict keys and filenames
`python/harvest_kmindex.py:119` and `:133-134`. **PLAUSIBLE.**

`target_name.split()[0]` raises `IndexError` on an empty or whitespace-only target name. The only
`try/except` in that loop wraps `fut.result()` (`:107-110`), so this escapes `as_completed`, tears
down the executor and kills a multi-hour harvest with **nothing written**, because all merging is
in memory and files are written only at `:130-146`.

`out_dir / f"{target}.accessions.txt"` uses the value verbatim, so a target containing `/` writes
outside `out_dir` or raises, and one containing `..` escapes it.

**FIXED 2026-09-21.** A `safe_target_name` helper replaces the inline `split()[0]`: it returns
`None` for anything unusable instead of raising, strips path separators, and refuses `.` and `..`.
The merge loop counts and skips those, and the whole per-shard payload merge is now wrapped so
that no malformed shard can tear down the executor, which was the part that turned a data defect
into total data loss.

Worth recording: **the empty-name case is reachable without a malformed shard.** The existing
`phiX174` strip turns a target named exactly `phiX174` into an empty string, which then wrote
dotfiles. That is the same defect arriving through ordinary data rather than a hostile name, and
it is now rejected with the rest.

**Verified** by four tests covering blank and whitespace names, the `phiX174` collapse, separators
and traversal, and that ordinary names such as `candida_auris_erg11` survive intact.
Mutation-tested: removing the blank guard fails the suite.

### 3.5 One malformed field aborts the stream after outputs are partly written
`python/lexicmap_streamer.py:499-503`, bare `float()` and `int()` on `pident`, `qcov`, `qstart`,
`qend`.

**CONFIRMED, reproduced** with a row whose `pident` is `NA`:
`ValueError: could not convert string to float: 'NA'`. The five cohort handles opened at `:353-357`
are closed by the `finally` at `:532`, but `.clean.accessions.fasta`, `.flagged.tsv` and
`.cohort_ledger.tsv` are left holding a truncated cohort while `.summary.json` is never written, so
the output directory looks half-plausible. Realistic triggers: concatenated LexicMap outputs with
an embedded second header row, `NA` or `-` in a numeric column, a `qseq` containing a tab.

**FIXED 2026-09-21.** The four numeric parses are wrapped: a row that will not parse is skipped
and counted, the first five are named individually on stderr, and the total is reported at the
end. The run therefore completes and `.summary.json` is written, which is what made the old output
directory misleading.

Skipping silently would trade a crash for quiet data loss, so there are two guards against that:
the count is always reported, and a file in which **every** data row failed raises instead, on the
same reasoning as **1.2** and **1.7**. All rows unparseable is a schema mismatch, not dirty data,
and the cohort would be empty rather than incomplete.

**Verified** by three tests: the `NA` case completes with the good row kept and `.summary.json`
written; an embedded second header row is skipped and both real accessions survive; an all-bad
file raises naming the row count. Mutation-tested: narrowing the `except` fails three tests.

### 3.6 Empty or header-only reference gives `ZeroDivisionError`
`python/lexicmap_streamer.py:424` (`cov = covered_bases / ref_len`). `load_reference`'s only guard
is `path.exists()`, and nothing checks `ref_len > 0`. **CONFIRMED** with a zero-byte reference. A
reference path that is a *directory* also passes `path.exists()` and raises `IsADirectoryError`
from `open()` at `:158` rather than the intended `FileNotFoundError`.

**FIXED 2026-09-21.** Both in `load_reference`, the single choke point every coverage figure
divides through: the existence check is now `path.is_file()`, so a directory raises the intended
`FileNotFoundError`, and an empty or header-only file raises a `ValueError` naming the path, the
record count and the first header. Previously it sailed through to `cov = covered_bases / ref_len`
and died there with a bare `ZeroDivisionError` that named neither the file nor the reason.

**Verified** by four tests: empty file, header-only file, a directory, and that a real reference
still loads. Mutation-tested: removing the empty guard fails two tests.

### 3.7 `-a` bypasses accession validation
`walker/src/main.rs:278` versus `io/input.rs:57-78`. The `--accession-list` path calls
`validate_accession`; the single-accession path does not. The value is interpolated into a URL and
two filesystem paths (`extract.rs:264`, `output.rs:284`), so `-a ../../foo` escapes both `out_dir`
and the scratch dir. Low severity, since it is the user's own command line, but the validator
already exists. *(verified directly: `validate_accession` is called only from `input.rs:32`.)*

**FIXED 2026-09-21.** `validate_accession` now takes a `source` string (`"--accession"`, or
`"list.txt:12"`) instead of a path and a line number, which is what had confined it to the list
path, and `main.rs` calls it on the `-a` value. It also rejects an empty accession, which the old
character check accepted.

**Verified** end to end: `-a '../../foo'`, `-a 'SRR8209082/../evil'` and `-a ''` are now usage
errors naming `--accession`, while `-a SRR99999991` proceeds to the fetch as before. Three unit
tests cover traversal and separators, real accessions of every SRA prefix, and that the message
names its source.

---

## 4. Correctness gaps in the streamer

### 4.1 `qlen` and the `query` column are never read
`python/lexicmap_streamer.py:325-346`. **Corrected 2026-09-21**: an earlier draft said these were
parsed and then ignored. They are not parsed. `find_col` resolves seven columns (`sgenome`,
`pident`, `qcovHSP`, `qstart`, `qend`, `qseq`, `sseq`) and neither `qlen` nor `query` is among
them; `grep -n "qlen\|idx_query"` returns nothing anywhere in the file.

Two consequences. There is **no reference/query length agreement check**, though the TSV carries
`qlen` per row, so a reference FASTA of a different length or gene version mis-projects every HSP
onto the grid with no warning. And **multi-target TSVs are merged onto one reference**: a LexicMap
search against the 160-target master panel emits rows for all targets in one file, distinguished
only by the `query` column, but grouping is by `sgenome` alone (`:497`, `:513-519`), so one
accession hit by several genes has HSPs from *different genes* merged into a single grid. If the
file is sorted by query rather than accession it instead trips the `ValueError` at `:411`. Rows
should be filtered to `query == ref_name`, which means adding that column to `find_col` first and
deciding what to do with a file that does not carry it.

**FIXED 2026-09-21.** Both columns are now resolved, and both are **optional**: older or
hand-made TSVs may not carry them, and requiring them would reject working inputs. Each absence
prints a line saying which check is therefore not running.

With a `query` column, rows naming another query are skipped and counted. With a `qlen` column,
a row whose `qlen` disagrees with the loaded reference **raises**, rather than warning: every HSP
is projected onto a grid of `len(ref_seq)`, so a different gene version mis-places all of them and
the output would be uniformly wrong rather than partly wrong. The message names both lengths and
says to use the reference the search used.

Filtering everything out also raises, on the same reasoning as **1.7** and **3.5**: every row
naming a different query means the reference record name does not match the search's query name,
not that there were no hits, and an empty cohort would look plausible.

**Verified** by six tests, three of which were sharpened after they turned out not to
discriminate. The multi-target case now gives the two genes **different halves** of the grid,
because overlapping spans hid the defect: overlap resolution dropped the intruder anyway, so the
first version of the test passed with the filter removed. It now asserts the accession tiles to
`ATGGCT------` and that geneB's bases never appear. The sorted-by-query case was likewise rebuilt
so the accession genuinely reappears after its block closes, which is what trips the `ValueError`.
With the filter removed, all three fail.

### 4.2 Dead branch in mutation classification
`python/lexicmap_streamer.py:662-663`. Both `not aa_muts` cases are consumed by the first two arms,
so `aa_muts` is always truthy by the fourth and `"Complex"` can never be emitted.

Related: a haplotype differing from the reference only at positions it does not cover (all `-`)
produces empty `nt_muts` and `aa_muts` and is labelled `WT` in the TSV, while `summary.json`'s
`wildtype_count` (`clean_full_haps.get(ref_seq, 0)`) excludes it, so the two outputs disagree on
how many wild-type sequences there are.

**FIXED 2026-09-21.** `Complex` is removed rather than made reachable: the arms above it already
cover every case, and inventing a meaning for it would be inventing semantics the caller never
asked for. The final arm is now plainly `Missense`, with a comment recording why there is nothing
below it.

The `WT` disagreement is resolved by making the label mean the same thing in both outputs. A
sequence that matches the reference everywhere it is covered but is not equal to it is now
`WT_partial`, not `WT`, so counting `mutation_class == "WT"` in the TSV equals `wildtype_count` in
`summary.json`. The non-CDS path had the identical problem and gets the parallel
`NonCoding_WT_partial`. `README.md`'s class list is updated.

**Verified** by five tests covering Missense, Synonymous, Nonsense, an exact-match `WT`, and a
partially covered match now classified `WT_partial`. Mutation-tested: reverting `WT_partial` to
`WT`, or restoring the dead `Complex` pair, each fails.

Worth recording for whoever reads the class list: **`Nonsense` is nearly unreachable in the clean
table by design.** It requires an amino-acid change to a stop, which is an internal stop, which
the premature-stop screen flags out of the clean cohort. The test has to pass
`max_internal_stops=1` to observe it. That is the screen working, not a defect, but it means a
cohort of clean haplotypes will essentially never show the class.

---

## 5. Performance

### Measured 2026-09-21, before deciding what to fix

Everything in this section was timed against real data rather than estimated, because the
answer turned out to depend on it. Stage timings, single stream, this link:

| accession | compressed | stream | parse | anchor | extract+phase | subgraph nodes |
|---|---|---|---|---|---|---|
| `SRR8209082` | 1.8 MB | 0.39-0.53 s | 0.043 s | 0.012 s | 0.024 s | **34,582** |
| `ERR2028820` | 32.4 MB | **5.41 s** | 0.636 s | 0.277 s | ~0.000 s | 747 |
| `SRR12363134` | 46.3 MB | **6.25 s** | 0.882 s | 0.352 s | ~0.000 s | 296 |

Two things follow, and they decided which findings were worth acting on.

**The walker is I/O bound by roughly an order of magnitude.** Streaming is 85-90% of pre-align
time on the large accessions, at 6.0-7.4 MB/s. Against the 17.4 MB/s ceiling recorded in the
optimisation baseline, and with 8 compute workers, compute occupies about 6% of the wall-clock
budget. Shaving CPU out of the walker does not move throughput.

**The expensive-subgraph accession is the cheap-to-fetch one.** `SRR8209082` has 34,582 subgraph
nodes and 1.8 MB to download; the two goldens have 296 and 747 nodes and 30-46 MB. 5.2 and 5.5
scale with subgraph size and the I/O scales with file size, and in this data they are
anti-correlated, so the stages those findings target are cheapest exactly where the accession is
slowest.

| finding | before | after | measured cost | measured saving | status |
|---|---|---|---|---|---|
| **5.2** view rebuild | O((2+k)(n+E)) | O(n+E) | 3.7 ms per build at n=34,582 | k=3: 15 ms; k=64: 243 ms | **skipped** |
| **5.3** NW | 2 x O(n*m) | O(n*b) | 12.5 ms, 1.6 MB (pfcrt) | **22-313x, see below** | **FIXED** |
| **5.4** double write | 2 x O(H*L) | O(H*L) | 1.90 s per 2,000 x 3 kb | **24%** | **FIXED** |
| **5.5** trail-cloning BFS | O(V*d) | O(V+E) | 6.5 ms exhaustive BFS | 1.6 ms | **skipped** |

**Why 5.2 was skipped.** The saving is real but small at realistic `k`. `k` is the number of
sites where the route took the minority allele, and the ERG11 panel evidence (16 of 50 accessions
differing by 1 to 3 substitutions) puts it at 1-3, not the 64 the cap allows: 15 ms, against a
~2.6 s wall-clock budget per accession, on the one accession shape where the subgraph is large
enough to matter. It is a constant-factor win that removes a multiplier of `k`, not an asymptotic
one, and it touches `view`, `find_alts` and `phase` at once, which is more churn than 0.5% of
compute justifies. Worth doing if the walker ever becomes CPU-bound, for instance behind a local
unitig cache.

**Why 5.5 was skipped.** It is the only genuine asymptotic improvement in this section, O(V*d) to
O(V+E) with `d = alt_max_detour = 100`, and it still measured **6.5 ms against 4.9 ms** on an
exhaustive BFS over the 34,582-node subgraph (54,722 nodes enqueued, 142,179 oriented edges). The
bound is loose in practice: most enqueued nodes sit at shallow depth, so the average cloned trail
is nowhere near `d`. 1.6 ms per no-rejoin BFS does not pay for rewriting a correctness-sensitive
traversal. The asymptotic argument would start to bite on a subgraph far deeper than any observed
here.

### 5.1 Memory ceiling is understated roughly 5x
> **FIXED 2026-09-21.** The claim is corrected, the real bound is a tested function, and the
> binary prints it at startup.

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

**Re-derived 2026-09-21 before fixing, and the count holds.** A buffer is owned by the fetch
worker from the moment `read_to_vec` returns until its blocking `send` completes, then by the
channel, then by the compute worker until `drop(fetched)` once the CSR exists. Those are three
disjoint holders, so `4 + 8 + 8 = 20` is right, and `max_compressed_bytes` is 1 GB.

**Fix as applied.** The wrong sentence is replaced with the derivation and the arithmetic, in
both places it appeared: the module doc, and `run_one`'s "Memory shape" section, which told the
reader to multiply by `compute_workers` for "the real ceiling" and so understated it the same
way. The diagram annotation that called capacity "the hard memory ceiling" is corrected too, and
its stale `~20 MB` CSR figure is now `~43 MB`, matching `types.rs`, `parse.rs` and the README.

Rather than leave the number in prose where it went stale once already, it is now
`Config::peak_compressed_buffers()` and `Config::max_compressed_resident_bytes()`, with tests
pinning 20 buffers and 20 GB at the defaults, ~6.9 GB at ERG11 panel size, that each of the three
terms moves the bound, and that an absurd `--max-size-mb` saturates rather than wrapping. Every
term is a flag or a config field, so the binary prints the computed figure at startup:

```
[logan-walker] memory bound: up to 20 compressed buffers resident (4 queued + 8 fetching +
8 parsing) x 1024 MB cap = 20.0 GB worst case, plus ~43 MB per compute worker for the CSR graph
```

**Verified** live across three configurations: defaults give 20 buffers / 20.0 GB;
`-w 16 --fetch-workers 16` gives 36 buffers / 36.0 GB; and `--max-size-mb 350 -w 4
--fetch-workers 4` gives 12 buffers / 4.1 GB.

That middle case matters: the README's own cohort example recommends `-w 16 --fetch-workers 16`,
which nearly doubles the bound, so the README now says so next to the example and carries a
**Memory** section with the derivation. No behaviour changed; the accounting is now stated where
it can be checked instead of asserted where it was wrong.

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

**FIXED 2026-09-21, both halves.**

*Computed once.* `snv_sites_all` runs the alignment once per result and both writers take the
rows, so enabling both flags no longer doubles the work. `main` computes it when either flag is
on.

*Banded.* The two sequences are the same gene at high identity, so the optimal path hugs the main
diagonal. `align_primary_to_query` now runs a banded DP whose width starts at
`max(32, |n-m| + 16)`. Memory becomes O(n*b) rather than O(n*m).

The band is a speedup and never a different answer: the traceback reports whether it touched a
band edge, and if it did the band is doubled and the alignment redone, falling back to the exact
full-matrix `align_full` at full width. Worst case is therefore exactly the old behaviour.

**Measured**, banded against exact:

| query | full DP | banded | speedup | matrix, full -> banded |
|---|---|---|---|---|
| pfcrt 1275 | 12.5 ms | 0.57 ms | **21.9x** | 1.6 MB -> 0.08 MB |
| ERG11 1575 | 13.4 ms | 0.43 ms | **31.2x** | 2.4 MB -> 0.10 MB |
| hiv1_pol 2739 | 14.9 ms | 0.39 ms | **38.1x** | 7.2 MB -> 0.17 MB |
| 10 kb | 221.2 ms | 1.46 ms | **151x** | 95.4 MB -> 0.62 MB |
| 20 kb (guard cap) | 918.9 ms | 2.94 ms | **313x** | **381.5 MB -> 1.24 MB** |

The last row is the reason this was worth doing ahead of the rest of section 5: the time was never
the problem, a 381 MB single allocation on top of the compressed-buffer ceiling in **5.1** was.

**Verified exact** by three tests, including 400 random sequence pairs carrying substitutions and
multi-base indels, each asserting the banded result equals `align_full`'s, plus a deliberately
divergent pair that forces the band fully open and degenerate length shapes.

### 5.4 Every streamer output is written twice, and mutation calling runs twice
`python/lexicmap_streamer.py:605`, `:617`, `:624`. The `.clean.*` files and the undocumented
`.msa.fasta` / `.haplotypes.tsv` / `.a2m` compat aliases are byte-identical siblings. For the TSV
the duplication is not only I/O: the loop body re-runs `translate_nt(seq)` and the per-position
nucleotide and amino-acid diff for every haplotype on the second pass. On a cohort with ~10^5
unique 3 kb haplotypes that is a second full O(haplotypes x ref_len) translation pass and a second
multi-GB write for no new information.

**FIXED 2026-09-21.** Both files are still written, because they are documented outputs and
dropping one is a separate decision; what stopped is regenerating the second. Each record is now
formatted once and written to both handles, for the MSA FASTA and for the haplotype TSV, whose
loop body is the expensive one (`translate_nt` plus the per-position nucleotide and amino-acid
diff per haplotype).

**Measured** on 2,000 distinct 3 kb haplotypes: **1.90 s to 1.45 s, 24% faster**, with all seven
output files byte-identical before and after. Output bytes are unchanged at 48.9 MB, which is the
point: an earlier probe that deleted the aliases entirely ran in 1.50 s, so essentially the whole
saving was recomputation rather than I/O, and keeping the aliases costs nothing.

Scaling linearly in haplotype count, the review's 10^5-haplotype cohort saves roughly 22 s per
target. Unlike the walker findings this is not gated on S3: the streamer never touches unitigs.

Relatedly, the module docstring's "O(1) streaming memory" claim for the collapser is wrong:
`clean_full_haps`, `clean_full_pident_sum` and `clean_full_samples` (`:369-371`) are keyed by the
full merged sequence and grow with the number of *distinct* haplotypes, which on a divergent target
approaches the accession count.

**The docstring half is FIXED 2026-09-21; the double write is not.** The false claim is a
documentation defect of the same kind as **5.1**, so it is corrected here rather than held behind
the performance work: the module docstring and the README now say memory is O(distinct
haplotypes) and that it is the row count that streams, which is the unbounded quantity. The
duplicate writes and the second translation pass are untouched and remain open.

### 5.5 `rejoin`'s BFS clones the whole trail per enqueued node
`walker/src/stages/variants.rs:662-686`, `let mut t = trail.clone(); t.push(next)` per edge, with
`trail` up to `alt_max_detour` (default 100) long. **PLAUSIBLE**, not measured. The BFS stops at the
first rejoin, normally 1 to 2 hops, so the common case is cheap; the expensive case is the
`no_rejoin` outcome the stage's own log counts, where the BFS exhausts the reachable set, up to
69k nodes on `SRR8209082`, each copying up to 100 `u32`. Parent-pointer BFS with a single
reconstruct at the end removes it.

---

## 6. Test and CI defects

> **ALL FIXED 2026-09-21.** Each item carries its fix below. Verified by running the suites
> against the real sieve binary (submodule initialised and built with `make serial`), and the
> golden gate against stubbed walkers covering pass, regression and infrastructure-failure.

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

**Fixed.** The count and the exit status are separated: the count goes into a local, and
`run_sieve` returns 0 explicitly when the sieve exited 0. **Verified** with a stub sieve that
exits 0 and promotes nothing: the script now reports `FAIL: run1 promoted 0 < 1 (nothing
promoted)` instead of `FAIL: sieve run 1 did not exit 0`, so the assertion at `:95` is reachable.
A stub that genuinely exits 3 still reports `sieve run 1 did not exit 0`, and the real sieve still
passes (1/33 promoted, deterministic).

### 6.2 `test_full_length_only_requires_complete_coverage` does not test the flag
`python/test_lexicmap_streamer.py:158-172`. SRR010 covers 100% and SRR011 covers 50%, and the test
asserts `clean_full_length == 1`. With the flag **off** the threshold is `min_coverage = 0.80`,
which SRR010 also clears and SRR011 also misses, so the assertion holds either way.

**CONFIRMED** by running the fixture both ways: both give `clean_full_length=1`. Deleting
`full_cov = 1.0 if full_length_only else min_coverage` (`lexicmap_streamer.py:441`) entirely would
leave this test green. A discriminating fixture needs an accession in the 0.80 to 0.99 band, for
example `qend=11` giving cov 0.917.

**Fixed** with exactly that: `SRR012` at `qend=11`, cov 0.917, clean-full with the flag off and
not with it on. The test now runs the fixture **both** ways in one case and asserts the two
results differ, so the discrimination cannot quietly decay again.

**Verified by mutation**, using the review's own criterion: replacing that line with
`full_cov = min_coverage` now fails the suite with `AssertionError: 2 != 1`, where before the
mutant survived.

### 6.3 `test_sieve_integration.py` reuses one record file across both runs
`tests/test_sieve_integration.py:141-148`. The off-path run's return code is discarded and it
reuses the same `work/stage2_input.txt` the on-path run wrote, so the
`else: check(False, "off-path: stage 2 was invoked")` branch at `:148` is unreachable, and a second
run that dies before stage 2 is misreported as "stage 2 got the wrong `--input`" rather than
"stage 2 never ran". Use a distinct record file per run.

**Fixed.** The stage-2 stub records its `--input` next to the alignment it was asked to write
(`Path(aln).with_name('stage2_input.txt')`), so each run writes into its own `-o` directory, and
the off-path run's return code is now captured and checked.

**Verified**, and the reachable branch immediately earned its keep: while the `fasttree` preflight
added in **3.1** was unstubbed, the off-path run died before stage 2 and the test correctly
reported `FAIL: off-path: stage 2 was invoked`, which is precisely the diagnosis this branch
could not previously give. With both tools stubbed, all seven checks pass.

### 6.4 `check_golden.sh` discards the walker's exit status and stderr
`bench/check_golden.sh:37-39` runs the walker with `>/dev/null 2>&1` and never inspects the return
code. If the Logan S3 fetch times out, `recon.fasta` is absent, `.got` is empty, and the gate
prints `FAIL ERR2028820 (got 0 nt vs golden 1476 nt)` and exits 1, reporting an infrastructure
outage as a walk regression with every diagnostic thrown away. Capture the rc and echo the log tail
on non-zero.

**Fixed.** `run_primary` keeps stderr in `$WORK/<acc>.log` and returns the walker's exit status;
a non-zero status prints `ERROR <acc> (walker exited non-zero: a run failure, NOT a walk
regression)` followed by the log tail. This works precisely because **1.2** made the walker exit
non-zero when nothing reconstructed, which is the S3-timeout case; before that it exited 0 and
there was no status to catch.

Writing this surfaced a latent bug in the fix itself, worth recording: `local acc="$1"
log="$WORK/$acc.log"` expands `$acc` **before** the assignment takes effect, so the log was named
from whatever the caller last left in the global `acc`. The two gates masked it by coincidence
(same name), and only the ungated `SRR8209082` run exposed it, by writing its log over
`SRR12363134.log`. Split into two statements, with a comment.

The gate also prints the GraphAligner version and path it ran with, because the fixtures are
version-dependent: 1.0.19+ clips the query's final base, so these goldens PASS on the pinned 1.0.20
and would FAIL on 1.0.18. A green gate now says which behaviour it green-lit. See
`docs/issue_graphaligner_end_clipping.md`.

**Verified** against three stub walkers: one that exits non-zero (reports ERROR with the log tail,
for all three accessions), one that emits the committed goldens (`PASS`, `golden: concordant`), and
one that exits 0 with a wrong sequence (still `FAIL ... got 8 nt vs golden 1274 nt`).

### 6.5 `tests/test_sieve.sh:63` bound is decoupled from the fixture
An unexpected fixture size only prints a note, while `MAX_PROMOTED=3` stays enforced regardless, so
growing the fixture silently changes what the bound means.

**Fixed.** The expected size is now a named constant next to the bound it calibrates, and a
mismatch is a hard `fail` telling the maintainer to re-derive `MAX_PROMOTED` and update both.
**Verified** by truncating the fixture to 2 records: the script exits 1 with
`fixture has 2 haplotypes, expected 33. MAX_PROMOTED=3 is calibrated against that fixture...`.

---

## 7. Smaller items

- ~~`python/disassemble.py:140, 169, 195`: all three `subprocess.run(..., check=True)` calls surface
  tool failure as a raw `CalledProcessError` traceback, inconsistent with the `die()` style used
  for every other error in the file.~~ **FIXED 2026-09-21.** A `run_step()` helper reports the
  stage and exit code via `die()`. This became urgent once **1.2** made logan-walker exit non-zero
  on an all-failed shard, which is a routine outcome and was surfacing as a traceback.
- ~~`python/disassemble.py:125, 127, 131, 133`: the `is not None` guards are always true, since the
  defaults are non-`None` numbers. Harmless today only because `extra` is appended last and the
  walker's parser is last-wins on duplicates, which is what makes the docstring's
  `-- --prune-frac 0.0` example work.~~ **FIXED 2026-09-21.** The four argparse defaults are now
  `None`, so the guards mean what they say and the walker's own defaults apply unless the user
  asks otherwise. Each help string names the walker default for reference instead of duplicating
  it as behaviour. **`-w` and `--fetch-workers` had the identical problem one line above** and are
  fixed with them. Verified: with no flags the walker is invoked with none of the six; with
  `--prune-frac 0.0 --min-vaf 0.25 -w 4` exactly those three appear.
- ~~`python/disassemble.py:96`: `Path(walker).is_file()` accepts a non-executable file, giving an
  uncaught `PermissionError` at `:140` instead of the intended `die()`.~~ **FIXED 2026-09-21.** A
  `resolve_executable()` helper applies `os.access(..., X_OK)` to an explicit path, which is the
  test `shutil.which` already applies to PATH entries, and `--sieve-bin` gets the same treatment.
  Verified on a mode-644 file: `logan-walker not found, or not executable, at ...`, exit 1.
- ~~`python/disassemble.py:108`: the file handle from `query.open()` is never closed, and the read
  happens unconditionally even when `--sieve` is off.~~ **FIXED 2026-09-21** as part of section 0;
  `first_record_len()` uses a `with` block and is called only on the `--sieve` path.
- ~~`python/disassemble.py:152`: `--sieve --skip-tree` is a silent no-op. The design choice is
  reasonable and commented, but nothing tells the user and the `--sieve` help text does not mention
  it. A user asking for a thinned FASTA without a tree gets the unthinned one and exit 0.~~
  **FIXED 2026-09-21** as the review frames it, by telling the user rather than changing the
  design: a warning naming what was ignored, why, and what to drop to get the thinned file, plus a
  line in the `--sieve` help. The alternative reading, that `--sieve --skip-tree` should run the
  sieve and write `reconstructed.sieved.fasta` without a tree, is a behaviour change and is left
  as a decision rather than made quietly.
- ~~`python/build_walker_codon_msa.py:83, 101`: temp files are `out_aln.parent/.tmp_<input name>`, so
  two concurrent runs writing into one output directory clobber each other's MAFFT input, and the
  temps leak if MAFFT is interrupted.~~ **FIXED 2026-09-21.** Every temp in that function now
  carries the PID, matching the scratch fix in **1.1**. The leak on interruption is unchanged for
  the two MAFFT inputs; the two output temps are removed by a `finally`.
- ~~`python/build_walker_codon_msa.py:118-125`: `out_aln` is opened for writing before MAFFT runs, so
  a MAFFT failure replaces a previously good alignment with a truncated or empty file before
  `sys.exit`.~~ **FIXED 2026-09-21** alongside **3.1**: temp file, then `os.replace` on success.
- ~~`python/harvest_kmindex.py:49`: `fetch_shard`'s `api_key` parameter is unused, since auth comes
  from the session headers set at `:91`; the `requests.Session` is never closed.~~
  **FIXED 2026-09-21.** The dead parameter is gone from the signature and both call sites, and the
  executor block is wrapped in `try`/`finally` so the session closes on every path, which matters
  at 32 pooled connections against one host. Two tests cover the close, on the success path and
  when the harvest raises.
- ~~`walker/src/pipeline.rs:26` and `config.rs:236-237` say the compute pool is sized to the core
  count; `compute_workers` is a hardcoded 8 and never queries the machine.~~ **FIXED 2026-09-21 by
  correcting the docs, not the code.** Auto-sizing would be the wrong fix: each compute worker
  holds a compressed buffer through the whole parse, so `available_parallelism()` on a 32-core
  machine means 4 + 8 + 32 = 44 buffers, 44 GB at the default `--max-size-mb` (see **5.1**). The
  diagram, the `config.rs` doc comment and the README row now say it is a fixed default of 8 that
  is deliberately not auto-sized, and point at the memory bound before raising it.
- `assets/candida_auris_venn.png` (310 KB) and `.svg` (43 KB) are tracked but referenced by no
  markdown file other than this one. Added by commit `a3999eb` with no accompanying write-up. Either the figure lost
  its document or the document was never written. Left in place deliberately, since deleting
  committed assets should be a human call. **Still open 2026-09-21, and still a human call:** the
  choice is between writing the document the figure belongs to and deleting both files, and
  neither is a decision to make on someone else's behalf.

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

1. ~~**1.1** scratch-path collision and **1.2** exit code. Both silently produce wrong results that
   look like success, and both are small changes.~~ **DONE 2026-09-21.**
2. ~~**1.3** sort results by input index. Fixes reproducibility and the nondeterministic sieve
   representative together.~~ **DONE 2026-09-21.**
3. ~~**Section 0**, the multi-record FASTA handling, as one coordinated change across the three call
   sites. Add a `cov <= 1.0` assertion while there.~~ **DONE 2026-09-21.**
4. ~~**3.1** and **3.2**, the two toolchain-resolution problems, which are cheap and currently make
   the environment pins unreliable.~~ **DONE 2026-09-21.**
5. ~~**5.1** correct the documented memory bound~~ **DONE 2026-09-21**. **5.3 also DONE** (the
   banded DP, which was about the 381 MB allocation rather than the time), and **5.4 DONE**.
   **5.2 and 5.5 measured and declined**, see the measurements at the top of section 5.
6. ~~The test defects in section 6, which currently mask the regressions they exist to catch.~~
   **DONE 2026-09-21.**
