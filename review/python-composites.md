# Python: composites and plumbing (cmd/pipe_*.py, cli.py, stage.py, resolve.py, validators.py, textio.py, cmd/walker.py, cmd/fetch.py)

The composites (`reconstruct`, `from-kmindex`, `from-lexicmap`) are already guarded in several places: `stage_args` catches renamed destinations, the `--` passthrough refuses the walker flags the composite owns, stale sieved copies are removed before each run, and `from-lexicmap` derives its sieve and MSA gates from lexicmap's own coverage gate. The weak points sit at the seams between stages. Handoff files are found by a fixed name and never checked to be this run's (PC-04, PC-05, PC-11). Gates in one stage are hard-coded rather than derived from the stage before, so rows that the walker or lexicmap accepted are silently dropped downstream (PC-02, PC-03), and `from-lexicmap` re-aligns rows that are already on the reference grid (PC-01). Binary resolution silently replaces a binary the user named (PC-06, PC-07). Most fixes share one theme: a composite should check what the previous stage actually wrote and inherit its settings, instead of copying literals and option declarations by hand (PC-19, PC-20).

| ID | Severity | Where | Title | Fix effort |
|---|---|---|---|---|
| PC-01 | medium | `src/disassembler/cmd/pipe_from_lexicmap.py:372` | from-lexicmap re-aligns reference-grid rows with MAFFT, moving bases off their LexicMap coordinates | small |
| PC-02 | medium | `src/disassembler/cmd/pipe_reconstruct.py:414` | `reconstruct --sieve` drops walker primaries at 85-95% of reference length | small |
| PC-03 | medium | `src/disassembler/cmd/pipe_from_lexicmap.py:339` | `from-lexicmap --sieve` drops clean rows whose in-frame deletion exceeds the 5% margin | small |
| PC-04 | medium | `src/disassembler/cmd/pipe_reconstruct.py:392` | Composites align a stale handoff FASTA from an earlier run when this run did not write it | small |
| PC-05 | medium | `src/disassembler/cmd/pipe_from_kmindex.py:342` | from-kmindex reconstructs an earlier harvest's cohort when the new harvest lacks `--target` | small |
| PC-06 | medium | `src/disassembler/resolve.py:47` | Unusable `--walker-bin`/`--sieve-bin` or `DISASSEMBLER_*_BIN` silently falls back to another binary | small |
| PC-07 | medium | `src/disassembler/resolve.py:37` | Relative `--walker-bin`/`--sieve-bin` (`./name`) is re-resolved through PATH by subprocess | trivial |
| PC-08 | medium | `src/disassembler/cmd/pipe_from_lexicmap.py:363` | Composites cannot set msa's `--min-identity`; `from-lexicmap --min-pident` below 60 fails at stage 3 | small |
| PC-09 | medium | `src/disassembler/cmd/walker.py:71` | Child-tool wrappers do not forward SIGTERM/SIGHUP, orphaning walker/MAFFT and misreporting signal exits | small |
| PC-10 | low | `src/disassembler/cmd/pipe_reconstruct.py:390` | reconstruct/from-kmindex turn walker exit 3 into 1 and have no working `--allow-empty` | small |
| PC-11 | low | `src/disassembler/cmd/pipe_reconstruct.py:432` | Stale-output cleanup keyed on this run's reference stem misses an earlier run's MSA and tree | small |
| PC-12 | low | `src/disassembler/cmd/pipe_from_lexicmap.py:206` | `-d` and `-f` are documented as mutually exclusive but both are accepted; Galaxy wins | trivial |
| PC-13 | low | `src/disassembler/cmd/pipe_from_kmindex.py:258` | from-kmindex and from-lexicmap reject bad inputs and missing tools only after stage 1 | small |
| PC-14 | low | `src/disassembler/cmd/pipe_reconstruct.py:198` | A stray positional token sends it and all later options to the walker passthrough | small |
| PC-15 | low | `src/disassembler/cmd/pipe_reconstruct.py:213` | Passthrough and validator ranges disagree with the walker (`-t`, `--prune-frac`, `--min-abund`, repeated flags) | small |
| PC-16 | low | `src/disassembler/cmd/pipe_reconstruct.py:361` | Empty `-a ""` crashes reconstruct and fetch with a TypeError traceback | trivial |
| PC-17 | low | `src/disassembler/textio.py:31` | textio's gzip reader never closes the file it opened; truncated gzip ends in a traceback | small |
| PC-18 | low | `src/disassembler/textio.py:48` | Inputs are read whole with no decompressed-size cap; small gzip inputs can cost GBs of RAM | small |
| PC-19 | low | `src/disassembler/cmd/pipe_reconstruct.py:251` | Sieve/MSA defaults copied as literals into composites and inert-knob checks | small |
| PC-20 | low | `src/disassembler/cmd/pipe_from_kmindex.py:235` | Composites re-declare about 40 stage options by hand; help text has already drifted | medium |
| PC-21 | low | `src/disassembler/cli.py:75` | `cli.main` has no top-level error handling: OS errors, truncated gzip, a closed stdout and Ctrl-C end in tracebacks | small |
| PC-22 | low | `src/disassembler/cmd/pipe_from_kmindex.py:236` | `--skip-tree` help says the sieve still runs; `--sieve --skip-tree` still requires a sieve binary | trivial |

### PC-01 · medium · from-lexicmap re-aligns reference-grid rows with MAFFT, moving bases off their LexicMap coordinates

`src/disassembler/cmd/pipe_from_lexicmap.py:372` · correctness · confirmed (reproduced) · fix effort: small · sources: PA-01, SC2-08

**Problem.** lexicmap writes every base at the reference position its HSP's `qstart` gives and `-` where no HSP reached, so `lexicmap.clean.accessions.fasta` is already an alignment on the reference grid. `from-lexicmap` still hands it to `msa_cmd.run(msa_ns)` (`pipe_from_lexicmap.py:372`), which always runs `mafft --auto --keeplength --addfragments` (`src/disassembler/cmd/msa.py:279-289`). MAFFT discards the input gaps and places the bases again. Where several placements score the same, it slides bases across an uncovered stretch or a deletion. The module docstring (`pipe_from_lexicmap.py:8-9`) says only that "MAFFT may place those gaps differently".

**Impact.** The codon MSA and the tree can disagree with lexicmap's own haplotype table and ledger, with exit 0 and the identity check passing. The trigger is the default `--min-cov 0.80` on genes covered by several HSPs with an internal hole, which is the case tiling exists for. Rated medium, not high: when HSP ends are trimmed the way a local aligner ends them (on matches, at least 50 nt), about 40% of affected rows still shift, but the base state changed in at most 3 of 300 rows at 0-1.5% divergence and in none at 0-8%. Most shifts are equivalent placements of the same gap.

**Evidence.** Reference HIV-1 pol (2739 nt), MAFFT 7.505; accession SRR0000001 carries C1416G and two gap-free HSPs (1-1418, 1425-2739). `disassembler from-lexicmap -q ref.fasta -f hits.tsv -o out` exits 0 and the haplotype table reports `1416C>G`, but the MSA row moved the SNV:

```
ref  ATTGGCAAGCCACCTGGATTCCTGA
lex  ATTGGCAAGCGAC------TCCTGA      SNVs [C1416G]
msa  ATTGGCAAGC------GACTCCTGA      SNVs [T1424C]  (codon 475 ATT>ACT, I>T)
```

In synthetic 200-row cohorts with HSP boundaries at arbitrary points, 27-39% of rows with an internal hole changed and 1-5 rows per cohort lost or gained a SNV. Rows truncated only at the ends never changed (0/200). The `--sieve` route behaves the same.

**Proposed fix.** Add an explicit `msa --aligned` mode and set it from `from-lexicmap`. Do not detect grid input by row length, because walker rows can be reference-length by chance.

```python
# msa.py
p.add_argument("--aligned", action="store_true",
               help="input rows are already on the reference grid; skip MAFFT")
mafft_bin = None if args.aligned else find_bin("mafft")
# row loop under --aligned: no reverse-complement; collect rows with len(seq) != ref_len
# after the loop: any off-grid rows -> error naming them, exit 1
# instead of MAFFT: write >Reference + rows; then the existing length check,
# identity_ok, os.replace and FastTree, unchanged
# pipe_from_lexicmap.py
msa_ns = stage_args(msa_cmd, ..., aligned=True)
```

Keep the coverage gate, ID sanitising, `Reference` as the first row and the identity check. With the fix, every MSA row is byte-identical to lexicmap's cohort and the run is about twice as fast. `from-lexicmap --mafft-threads` then does nothing: warn when it is set, as `pipe_from_lexicmap.py:249-266` already does for inert sieve options.

**Tests and docs.** New test: the codon MSA rows equal `lexicmap.clean.accessions.fasta`; C1416G is kept and 1419-1424 stay `-`; MAFFT is not invoked; `msa --aligned` with a row of the wrong length exits 1 and writes nothing. Update `pipe_from_lexicmap.py:3-10`, the `msa.py` docstring and help, `docs/codon_msa.md:121-129`, `docs/lexicmap.md:160-169`. The Galaxy `disassembler_msa.xml` (tools-iuc#8800) needs an "already aligned" boolean mapped to `--aligned`, and its claim that rows are "gap-free" at 100% coverage is wrong for rows with deletions.

### PC-02 · medium · `reconstruct --sieve` drops walker primaries at 85-95% of reference length

`src/disassembler/cmd/pipe_reconstruct.py:414` · correctness · confirmed (reproduced) · fix effort: small · sources: IC-10, PE-05, PC-05

**Problem.** Under `--sieve`, `reconstruct` sets the sieve's `--min-len` to `int(length * 0.95)` of the shortest reference (`pipe_reconstruct.py:414`). The walker emits a primary whenever `cov = target.len() / query.seq.len()` (`walker/src/stages/stitch.rs:198`) is at least `min_cov`, default 0.85 (`walker/src/config.rs:455`), and the MSA stage admits rows with 70% of the reference in ACGT. Every haplotype at 85-95% of reference length, which both neighbouring stages accept, is removed at a step documented as near-duplicate thinning. `from-lexicmap` derives its gates from the upstream gate for exactly this reason; `reconstruct` does not.

**Impact.** Walks that do not reach a terminus, and samples with a genuine in-frame deletion, vanish from the sieved FASTA, MSA and tree with exit 0. These are the divergent lineages thinning is meant to keep. The sieve log prints `completeness gate drops N`, but the final `[reconstruct] done` summary does not. Passing `-- --min-cov 0.6` to the walker widens the gap.

**Evidence.** Offline cache with two accessions against a 63-nt query, one carrying a 6-nt deletion (walker reports `cov=0.90`). `disassembler reconstruct -l accs.txt -q ref.fasta -o out --unitigs-cache cache --offline` puts both primaries in the MSA; adding `--sieve` prints `[sieve] completeness gate drops 1 of 2: 1 shorter than 59 nt (--min-len)` and the MSA holds only the reference and one primary.

**Proposed fix.** Derive the default from the walker gate in force:

```python
WALKER_DEFAULT_MIN_COV = 0.85            # walker/src/config.rs
def _walker_min_cov(extra):              # last `--min-cov X` in the passthrough wins, as in main.rs
    ...
frac = min(0.95, _walker_min_cov(raw_extra))
min_len = int(length * frac)             # print it always, not only in multi-anchor mode
```

No extra margin is needed: cov is measured against the winning anchor, which is at least as long as the shortest reference. The walker has only the `--min-cov X` form, no `=` form. Put the fraction in a shared constant (see PC-19). Longer term, gate failures could pass through the sieve unthinned instead of being dropped, which also removes PC-03.

**Tests and docs.** `tests/test_composites.py:277` pins `int(9 * 0.95)`; change it to `int(9 * 0.85)` and add a case with `-- --min-cov 0.6`. New test: a stub walker writing a 0.90-length primary reaches the MSA input under `--sieve`. Update `docs/sieve.md:40` and the `--sieve-min-len` help.

### PC-03 · medium · `from-lexicmap --sieve` drops clean rows whose in-frame deletion exceeds the 5% margin

`src/disassembler/cmd/pipe_from_lexicmap.py:339` · correctness · confirmed (reproduced) · fix effort: small · sources: PA-05, PC-05

**Problem.** lexicmap counts a trusted deletion as covered with no limit: `covered = len(based) + sum(e.length for e in indels if e.kind == "del" and not e.edge)` (`src/disassembler/cmd/lexicmap.py:420`), while the row keeps `-` there. `from-lexicmap` then sets the sieve's non-ACGT limit to `round(min(1.0, 1.0 - gate + GATE_MARGIN), 6)` (`pipe_from_lexicmap.py:339`), a fixed 0.05 allowance (`GATE_MARGIN`, `:46`) for deletions plus Ns. A clean row whose deletions exceed 5% of the reference is dropped by the sieve's completeness gate. This contradicts the module docstring (`:16-17`) and `docs/lexicmap.md:165-167` ("neither stage drops rows for coverage that lexicmap already accepted"). The same margin sets the MSA gate (`:362`), which can bite under `--include-partials`.

**Impact.** `from-lexicmap --require-full-length --sieve` on a cohort with an in-frame deletion longer than about 5% of the gene (45 nt in 900 nt) loses that row from the sieved cohort, MSA and tree with exit 0. Indel alleles are often what a study is after. The `[sieve] 6 -> 5 haplotypes` summary reads as thinning.

**Evidence.** A LexicMap table for a 900-nt ORF with five SNV samples and one sample with a 60-nt in-frame deletion (ledger: `CLEAN_FULL_LENGTH`, coverage 1.0). Without `--sieve` the MSA has 7 rows. With `--sieve`: `--max-gap-frac 0.05 (cohort coverage gate 1.0)`, `completeness gate drops 1 of 6: ... 1 with more than 5% gaps`, and the deletion sample is absent from the MSA.

**Proposed fix.** Never set the derived limit below the cohort's own worst row:

```python
max_gap_frac = min(1.0, max(1.0 - gate + GATE_MARGIN, _max_non_acgt_share(cohort)))
max_gap_frac = math.ceil(max_gap_frac * 1e6) / 1e6   # round up, never below a row's share
```

`_max_non_acgt_share` is one FASTA pass, done where `n_cohort` is counted (`:311`). Use the same floor for `msa_min_cov` (`:362`). A simpler equivalent is to default the limit to 1.0 in `from-lexicmap`, since lexicmap has already gated coverage. See PC-02 for the same defect in `reconstruct`.

**Tests and docs.** `tests/test_sieve_integration.py:392-399` (`TestFromLexicmapGates`) pins 0.25/0.05; update it and add a full-length cohort with a long in-frame deletion that must reach the MSA input. Update `docs/sieve.md:64-75` and the `--sieve-max-gap-frac` help.

### PC-04 · medium · Composites align a stale handoff FASTA from an earlier run when this run did not write it

`src/disassembler/cmd/pipe_reconstruct.py:392` · correctness · confirmed (reproduced) · fix effort: small · sources: PC-07, P4-04, PE-02, P4-13

**Problem.** `reconstruct` (and `from-kmindex`, which ends in it) finds the walker output by a fixed name and never checks that this run wrote it: `n_recon = sum(1 for line in recon.open() if line.startswith(">")) if recon.is_file() else 0` (`pipe_reconstruct.py:392`). The passthrough scrub (`:213-224`) refuses only `-f/-o/-a/-l/-q`, so `-- --fetch-only`, `-- --help` and `-- --version` reach the walker, which returns success without writing `-f` (`walker/src/main.rs:340-347`, `:515-522`). `from-lexicmap` has the same hole: it passes `out_prefix=<out-dir>/lexicmap` (`pipe_from_lexicmap.py:241`), `lexicmap` turns an existing directory prefix into `<dir>/<ref>.*` (`lexicmap.py:576-578`), and the composite then reads `<out-dir>/lexicmap.clean.accessions.fasta` (`pipe_from_lexicmap.py:294`).

**Impact.** In an `--out-dir` holding an earlier run, the user gets "reconstructed N haplotype(s)", an MSA and a tree of the old cohort, and exit 0, all presented as this run's. `docs/codon_msa.md:100` invites any other walker flag after `--`, and `--fetch-only` is a natural cache warm-up. For `from-lexicmap`, the trigger is a directory `<out-dir>/lexicmap` (for example from `disassembler lexicmap --out-dir results/lexicmap` followed by `from-lexicmap -o results`); without an earlier composite run it fails with "expected lexicmap cohort not found", with one it silently uses the old cohort.

**Evidence.** Seed `out/reconstructed.fasta` with two records, then run `disassembler reconstruct -a SRR0000001 -q ref.fasta -o out --unitigs-cache cache --offline -- --fetch-only`: the walker prints `1 fetched ...`, the composite reports `reconstructed 2 haplotype(s)`, and the MSA holds the two seeded records, exit 0. `-- --version` and `-- --help` behave the same. For `from-lexicmap`: run once into `out`, `mkdir out/lexicmap`, rerun with a different hits table; the new cohort lands in `out/lexicmap/<ref>.*` and the MSA is built from the old one, exit 0.

**Proposed fix.**

```python
# pipe_reconstruct.run, after the preflights
recon.unlink(missing_ok=True)
...
if not recon.is_file():
    return _die(f"logan-walker exited 0 but wrote no {recon}; nothing to align.")
# _scrub_walker_passthrough: also refuse --fetch-only (-> `disassembler fetch`),
# -h/--help (-> `disassembler walker -h`), -V/--version, with their own wording

# pipe_from_lexicmap.run, after prefix = out_dir / "lexicmap"
if prefix.is_dir():
    print(f"from-lexicmap: {prefix} is a directory; lexicmap would write inside it, "
          f"not at {prefix}.*. Remove it or choose another --out-dir.", file=sys.stderr)
    return 1
for sfx in (".clean.accessions.fasta", ".clean_expanded.accessions.fasta"):
    Path(f"{prefix}{sfx}").unlink(missing_ok=True)
```

The unlink is the real fix; the scrub gives a clearer message. Longer term, move the directory heuristic out of `stream_lexicmap_to_msa` into `lexicmap.run`'s argument handling so programmatic callers get the exact prefix they pass. Related: PC-05 (stale harvest files), PC-11 (stale MSA/tree).

**Tests and docs.** In `tests/test_composites.py`: a stale `reconstructed.fasta` plus a stub walker that only exits 0 gives rc 1 and no MSA; `-- --fetch-only`, `-h`, `--help`, `-V`, `--version` are refused; an `out/lexicmap/` directory beside a stale cohort gives rc 1. Update `docs/codon_msa.md:100-103` and the passthrough help (`pipe_reconstruct.py:186-196`) to list refused flags; `docs/lexicmap.md:124-125` should say `<out-dir>/lexicmap` must not be a directory.

### PC-05 · medium · from-kmindex reconstructs an earlier harvest's cohort when the new harvest lacks `--target`

`src/disassembler/cmd/pipe_from_kmindex.py:342` · correctness · confirmed (reproduced) · fix effort: small · sources: PE-03, P4-01, PC-04

**Problem.** After a fresh harvest, `from-kmindex` picks the accession list by name only: `acc_list = cohort_dir / f"{args.target}.accessions.txt"` then `if not acc_list.is_file():` (`pipe_from_kmindex.py:341-342`). Harvest writes per-target files only for this run's targets (`src/disassembler/cmd/harvest.py:179-186`) and never removes those of an earlier harvest into the same directory. With `--score-bin`, stratify scans the whole cohort directory (`pipe_from_kmindex.py:326`, `stratify.py:192-201`), so a stale `<target>.scores.json` is binned too. Nothing compares `--target` with what this harvest produced.

**Impact.** Re-running `from-kmindex` into the same `-o` (or a shared `--cohort-dir`) with a collection in which `--target` has no hits runs the walker, MSA and tree on the previous collection's accessions and exits 0. The output looks like a valid result for the new collection. In a fresh `-o` the same command correctly fails. Galaxy jobs get fresh directories, so this is a CLI re-run hazard. The same class was already fixed for `from-lexicmap` (`tests/test_composites.py:344`).

**Evidence.** Shard set A has target `ns5` (SRR_OLD1, SRR_OLD2); shard set B has only `HA`. `from-kmindex --target ns5 -o o --shards-dir shA`, then the same with `--shards-dir shB`: the second run exits 0, `summary.json` lists only `HA`, the walker gets `-l o/cohort/ns5.accessions.txt`, and the MSA rows are SRR_OLD1 and SRR_OLD2. With `--score-bin 0.0:1.0:10` the same happens via `ns5.bins.txt`.

**Proposed fix.**
1. After a successful (not skipped) harvest, read the just-written `cohort_dir/summary.json`; if `args.target` is not among its keys (other than `_union`), print `target 'ns5' has no hits in the harvest in <dir> (targets: HA)` and return 1. Do not apply this under `--skip-harvest`, where hand-made cohorts must still work.
2. After a fresh harvest, pass `scores=[cohort_dir / f"{args.target}.scores.json"]` to stratify; keep the directory scan under `--skip-harvest` (as `docs/kmindex.md:122` promises).
3. In `harvest._write_outputs`, read the previous `summary.json` and unlink `<t>.accessions.txt` / `<t>.scores.json` for each earlier target `t` not in this run, only when `safe_target_name(t) == t`.

**Tests and docs.** `TestFromKmindex`: harvest `gene`, then harvest a collection with only `HA` into the same `-o`, with and without `--score-bin`; assert rc 1, the message, that the walker did not run, and that the `gene.*` files are gone. A harvest unit test that a second harvest removes the dropped target's files. `docs/kmindex.md` (around line 106): `--target` must be among the targets of the harvest just run.

### PC-06 · medium · Unusable `--walker-bin`/`--sieve-bin` or `DISASSEMBLER_*_BIN` silently falls back to another binary

`src/disassembler/resolve.py:47` · provenance · confirmed (reproduced) · fix effort: small · sources: PE-06, BL2-04, IC-08, P2-04, P4-02, PA-03, PC-03, PX-09

**Problem.** Both resolvers are one `or` chain: `_executable(explicit) or _executable(os.environ.get("DISASSEMBLER_WALKER_BIN")) or _from_path("logan-walker") or _executable(<repo>/walker/target/release/logan-walker)` (`resolve.py:47-52`; `resolve_sieve` at `:57-64` has the same shape). `_executable` returns None for a missing or non-executable path (`:32-38`), so a choice the user made explicitly is treated as "not given" and the next source wins. No caller prints which binary was resolved, and `summary.json` records the crate version, not the path.

**Impact.** A typo, a moved build, a missing `chmod +x` or a stale `DISASSEMBLER_WALKER_BIN`, on a host that also has `logan-walker` on PATH (conda) or an in-tree build, runs the other binary with exit 0. A comparison of two builds compares one build with itself, and same-version builds cannot be told apart afterwards. The same holds for the sieve, and the composites' preflights report success.

**Evidence.** `disassembler walker --walker-bin /nonexistent/logan-walker --version` prints `logan-walker 0.1.1` and exits 0 (the in-tree build answered). `DISASSEMBLER_WALKER_BIN=/nonexistent/lw disassembler walker --version` does the same. With a stub `logan-walker` on PATH, the stub runs.

**Proposed fix.**

```python
def _resolve(explicit, env_var, name, in_tree):
    if explicit:
        return _executable(explicit)          # no fall-through
    env = os.environ.get(env_var)
    if env:
        return _executable(env)               # no fall-through
    return shutil.which(name) or _executable(in_tree)
```

`resolve_walker`/`resolve_sieve` call `_resolve(...)`. Add `unresolved_reason(explicit, env_var)` (for example `"--walker-bin '/x' does not exist"` or `"is not an executable file"`) and use it in the five not-found messages (`cmd/walker.py:52`, `cmd/sieve.py:155`, `cmd/pipe_reconstruct.py:322`, `:346`, `cmd/pipe_from_lexicmap.py:244`). Print the resolved path once per run (`[reconstruct] walker: /abs/path`). Ship together with PC-07, which changes the same function. Optionally probe the in-tree path only when `_repo_root()/"pyproject.toml"` exists, since for an installed package that root is `<prefix>/lib`.

**Tests and docs.** New `tests/test_resolve.py`: a missing explicit path, a non-executable explicit path and a bad env value each return None with a decoy on PATH; PATH is used when nothing explicit is given. Update the `resolve.py` docstring, `docs/walker.md:44`, `README.md:53`, `docs/sieve.md:42`: an explicit path or env value is final, and an unusable one is an error.

### PC-07 · medium · Relative `--walker-bin`/`--sieve-bin` (`./name`) is re-resolved through PATH by subprocess

`src/disassembler/resolve.py:37` · correctness · confirmed (reproduced) · fix effort: trivial · sources: PB-02, P4-05

**Problem.** `_executable` validates the path relative to the cwd but returns `str(q)` (`resolve.py:37`), and `str(Path("./logan-walker"))` is `"logan-walker"`. A name without a slash is then looked up on PATH by `subprocess.run([resolved, *argv])` (`cmd/walker.py:71`, used by `fetch` and every composite) and by `subprocess.run(cmd)` in `cmd/sieve.py:219`. The atomic `walker` uses `os.execv` (`cmd/walker.py:91`), which does not search PATH, so the atomic and composite paths run different files. The same applies to a bare name present in the cwd and to `DISASSEMBLER_*_BIN=./name`.

**Impact.** `cd walker/target/release && disassembler reconstruct ... --walker-bin ./logan-walker` runs a same-named binary from PATH (for example a conda install) silently; with none on PATH, the preflight reports the binary found and the stage then dies with an uncaught `FileNotFoundError` traceback. Likely for developers testing a local build.

**Evidence.** `python -c "from disassembler.resolve import resolve_walker; print(repr(resolve_walker('./logan-walker')))"` in a directory holding an executable `logan-walker` prints `'logan-walker'`. With a stub on PATH, `disassembler fetch -a SRR1 --unitigs-cache c --walker-bin ./logan-walker` runs the PATH stub; `disassembler walker --walker-bin ./logan-walker --version` runs the local one.

**Proposed fix.**

```diff
     if q.is_file() and os.access(q, os.X_OK):
-        return str(q)
+        return os.path.abspath(q)   # never let subprocess re-search PATH
```

`abspath` rather than `resolve()` keeps a symlinked build's name in logs. Ship with PC-06.

**Tests and docs.** `tests/test_resolve.py::test_relative_explicit_is_absolute`: chdir into a temp dir holding an executable `logan-walker`, put a decoy on PATH, assert `resolve_walker("./logan-walker") == str(cwd / "logan-walker")`. No doc change.

### PC-08 · medium · Composites cannot set msa's `--min-identity`; `from-lexicmap --min-pident` below 60 fails at stage 3

`src/disassembler/cmd/pipe_from_lexicmap.py:363` · usability · confirmed (reproduced) · fix effort: small · sources: BL2-03, P4-07, PB-05, PE-08

**Problem.** 0.1.1 added a median-identity gate to `msa` (`--min-identity`, default 0.60, `src/disassembler/cmd/msa.py:166`); on failure it refuses to write the alignment (`msa.py:129-134`). No composite exposes or forwards it: the msa namespaces at `pipe_from_lexicmap.py:363-371` and `pipe_reconstruct.py:440-448` set only ref, input, out_aln, out_tree, mafft_threads and msa_min_cov, so `stage_args` fills in 0.60. `from-lexicmap` accepts `--min-pident` from 0 to 100 and derives its other downstream gates from the cohort gate, but not this one. The error message names `--min-identity` and `--ref`, neither of which exists on the composite.

**Impact.** `from-lexicmap --min-pident 50` on divergent but genuine homologues (53-59% identity) admits them as clean, then exits 1 at stage 3 with no MSA or tree; adding `--min-identity 0` gives "unrecognized arguments". In `reconstruct`/`from-kmindex` the same failure comes after the whole walker stage. The workaround is to run `disassembler msa ... --min-identity 0` by hand on the kept cohort FASTA. This is a 0.1.0 -> 0.1.1 regression.

**Evidence.** Six 600-nt hits at 53-59% identity: `disassembler from-lexicmap -q ref.fasta -f hits.tsv -o out --min-pident 50` prints `median identity to the reference is 0.558, under --min-identity 0.60 ... not written` and exits 1. `disassembler msa -r ref.fasta -i out/lexicmap.clean.accessions.fasta -a m.fa -t m.nwk --min-identity 0` exits 0. `disassembler reconstruct -h | grep -c identity` gives 0.

**Proposed fix.** Name the flag `--msa-min-identity`, beside `--msa-min-cov` (a bare `--min-identity` in `from-lexicmap` would read as lexicmap's `--min-pident`).

```python
# reconstruct, from-kmindex
g_msa.add_argument("--msa-min-identity", type=fraction, default=0.60,
                   help="msa's --min-identity: median identity to the reference below "
                        "which no alignment is written; 0 disables")
# reconstruct msa stage_args: min_identity=args.msa_min_identity
# from-kmindex reconstruct stage_args: msa_min_identity=args.msa_min_identity
# _warn_inert_knobs: inert under --skip-tree
# from-lexicmap: default=None, then
msa_min_identity = args.msa_min_identity
if msa_min_identity is None:
    msa_min_identity = round(max(0.0, min(0.60, args.min_pident / 100 - GATE_MARGIN)), 6)
```

Print the derived value on the stage-3 line. The derived default at `--min-pident 60` is 0.55; keep 0.60 when `--min-pident` is left at its default if exact 0.1.1 behaviour matters. Optionally let `identity_ok` take the flag name so the message names the composite's flag. PC-20 would add this option once for all composites.

**Tests and docs.** Capture tests: `from-lexicmap --min-pident 50` forwards `min_identity == 0.45`, the default forwards 0.55 (or 0.60), an explicit value passes through; `reconstruct`/`from-kmindex` forward `--msa-min-identity`; `--skip-tree` warns that it is inert. Update `docs/codon_msa.md:45`, the option tables in `docs/lexicmap.md` and the reconstruct/from-kmindex docs. The Galaxy msa wrapper needs a `min_identity` parameter at the 0.1.1 bump.

### PC-09 · medium · Child-tool wrappers do not forward SIGTERM/SIGHUP, orphaning walker/MAFFT and misreporting signal exits

`src/disassembler/cmd/walker.py:71` · reliability · confirmed (reproduced) · fix effort: small · sources: PE-07, P4-03

**Problem.** Every composite, `fetch`, `sieve` and `msa` start their child with a bare `subprocess.run`, e.g. `return subprocess.run([resolved, *argv]).returncode` (`cmd/walker.py:71`; also `cmd/sieve.py:219`, `cmd/msa.py:289`, `:319`). Python installs no SIGTERM/SIGHUP handler, so a signal sent to the disassembler PID kills only Python: the walker (or MAFFT, FastTree, sieve) is reparented to init and keeps running, and Python's `finally` blocks (msa's `.tmp_*` cleanup, sieve's shard cleanup) never run. A child killed by a signal returns a negative code that `sys.exit` turns into `256-N` (SIGKILL gives 247, not 137).

**Impact.** Triggered by signals aimed at the Python PID only: `kill <pid>`, `Popen.terminate()` from an orchestrator, `timeout --foreground`. The walker keeps up to tens of GB and its `/dev/shm` scratch until it finishes, and writes `reconstructed.fasta` into an output directory whose run already failed with 143. Not triggered by terminal Ctrl-C, plain `timeout`, Galaxy's `killpg`, systemd or Slurm, which signal the whole group. Exit 247 instead of 137 hides OOM kills from scripts that test for 137.

**Evidence.** With a stand-in walker that sleeps 3 s then writes a marker, `disassembler fetch -a SRR1 --unitigs-cache cache --walker-bin stub.sh &` followed by `kill -TERM <wrapper pid>`: the wrapper exits 143, the child's parent becomes PID 1, and the marker is written 3 s later. A stub that SIGKILLs itself gives exit 247 through `fetch` and 137 through `walker`.

**Proposed fix.** One helper, `src/disassembler/childproc.py::run_child(cmd, name=..., **popen_kw)`, used at all five call sites:

```python
proc = subprocess.Popen(cmd, process_group=0, **popen_kw)    # Python >= 3.11 is the floor
# while it runs: SIGINT/SIGTERM/SIGHUP handlers call os.killpg(proc.pid, sig)
proc.communicate()
# restore handlers; rc < 0 -> 128 + N with "walker: logan-walker was killed by SIGKILL";
# if a signal was received: raise SystemExit(rc or 128 + sig) so finally blocks run
```

`fetch.run` may instead `os.execv` as `walker.run` does. `pipe_reconstruct` should return the walker's rc when it is 3 or at least 128 (see PC-10). In its own process group the child no longer gets tty SIGQUIT or SIGTSTP directly; forward SIGTSTP/SIGCONT or accept that. Handlers are installed only in the main thread.

**Tests and docs.** New `tests/test_childproc.py`: run `sh -c 'sleep 5'` via `run_child` in a subprocess, send SIGTERM, assert exit 143 and that the grandchild is gone; a child running `kill -9 $$` gives 137. `tests/test_msa.py:127` must patch `run_child` instead of `msa.subprocess.run`. Update the `walker.py:61-67` docstring and the "Exit status" section of `docs/walker.md`.

### PC-10 · low · reconstruct/from-kmindex turn walker exit 3 into 1 and have no working `--allow-empty`

`src/disassembler/cmd/pipe_reconstruct.py:390` · usability · confirmed (reproduced) · fix effort: small · sources: BL2-09, P4-08, PC-08, PS-04, PE-08

**Problem.** Any non-zero walker status becomes 1 (`return _die(f"logan-walker (stage 1) failed (exit {rc}); ...")`, `pipe_reconstruct.py:390`), so "every accession errored" (exit 3) looks like "nothing found". `if n_recon == 0: return _die("no haplotypes reconstructed; nothing to align.")` (`:394-395`) runs before the `--skip-tree` return, and neither composite has `--allow-empty`; `-- --allow-empty` makes the walker exit 0 but the composite still exits 1, although the walker's hint tells users to pass it. Batch schedulers cannot retry only failed shards.

**Proposed fix.** Add `--allow-empty` to both parsers (from-kmindex forwards it); treat `"--allow-empty" in raw_extra` the same; append it to `walker_argv` once. On walker failure print the message and `return rc if rc > 0 else 128 - rc`. Move the `sieved/aln/tree` setup above the empty check; when `n_recon == 0` and `allow_empty`, unlink those three and return 0.

**Tests and docs.** A fake walker exiting 3 gives rc 3; `--allow-empty` (flag or via `--`) gives 0, removes stale outputs and forwards the flag once. Document exit statuses 0/1/2/3 in `docs/kmindex.md` and `docs/codon_msa.md`.

### PC-11 · low · Stale-output cleanup keyed on this run's reference stem misses an earlier run's MSA and tree

`src/disassembler/cmd/pipe_reconstruct.py:432` · correctness · confirmed (reproduced) · fix effort: small · sources: SWP4-02

**Problem.** Output names come from the current reference (`stem = msa_ref.stem`, `pipe_reconstruct.py:398-400`), and the `--skip-tree` cleanup removes only `for stale in (aln, tree)` (`:432-433`). After `reconstruct -q gene.fasta -o out` then `reconstruct -q geneB.fasta -o out --skip-tree`, `out/gene_codon_msa.fasta` and `gene_fasttree.nwk` from the other cohort remain, contrary to `docs/codon_msa.md:98`. `from-lexicmap --allow-empty` (`pipe_from_lexicmap.py:307-321`) has the same key.

**Proposed fix.** Record derived outputs in a manifest in `-o` and clear it once the walker has rewritten the reconstruction; keep the stem unlink as a fallback; accept only bare names from the manifest. Do not glob `*_codon_msa.fasta`, which would delete user files.

```python
DERIVED_MANIFEST = ".disassembler_derived"
def clear_derived(out): ...            # unlink each bare name listed, then the manifest
def record_derived(out, *paths): ...   # write p.name per line
```

**Tests and docs.** Two runs with different `-q`, the second with `--skip-tree`, leave no `*_codon_msa.fasta`; a from-lexicmap twin. Mention the manifest in `docs/codon_msa.md` and `docs/lexicmap.md`.

### PC-12 · low · `-d` and `-f` are documented as mutually exclusive but both are accepted; Galaxy wins

`src/disassembler/cmd/pipe_from_lexicmap.py:206` · correctness · confirmed (reproduced) · fix effort: trivial · sources: P4-13

**Problem.** `from-lexicmap` and `lexicmap` check only that one of `-d`/`-f` is given (`pipe_from_lexicmap.py:206`, `src/disassembler/cmd/lexicmap.py:1279`), although the help says "Mutually exclusive". With both, the dataset branch wins (`lexicmap.py:540`, `:552`): the local file is silently ignored and Galaxy is contacted, giving an error or another dataset's data.

**Proposed fix.** In both `run`s: `if args.dataset and args.file: print("... -d and -f are mutually exclusive", file=sys.stderr); return 2`, or an argparse mutually exclusive group.

**Tests and docs.** `-d` with `-f` exits 2 in `lexicmap` and `from-lexicmap`. Scripts passing both today will now fail.

### PC-13 · low · from-kmindex and from-lexicmap reject bad inputs and missing tools only after stage 1

`src/disassembler/cmd/pipe_from_kmindex.py:258` · usability · confirmed (reproduced) · fix effort: small · sources: P4-09, PC-09, PS-02

**Problem.** The comment at `pipe_from_kmindex.py:258` promises to "Fail fast ... BEFORE spending a Galaxy round-trip on harvest", but only the cache, passthrough, harvest mode and API key are checked. A typo in `-q`, a missing `--msa-reference`, a bad `--score-bin`, `--prune-frac 1.0` or a missing tool is reported only after the harvest, which can take hours. `reconstruct` checks `--msa-reference` only with `is_file()` (`pipe_reconstruct.py:297`), so a reference without records fails after the walker. `from-lexicmap` checks for MAFFT/FastTree only after streaming the dataset. No output is wrong; `--skip-harvest` recovers, but the error does not say so.

**Proposed fix.** Move `reconstruct`'s checks (query and MSA reference exist and have records, cache, scrub, walker, MAFFT/FastTree, sieve) into `pipe_reconstruct.preflight(args) -> str | None`; call it, plus `stratify_cmd.assert_bins_disjoint([parse_bin_spec(b) for b in args.score_bins])`, before `cohort_dir.mkdir` in `from-kmindex`; call `msa_cmd.find_bin("mafft")` and `find_bin("fasttree", "FastTree")` before stage 1 of `from-lexicmap`.

**Tests and docs.** For each bad input, `from-kmindex --shards-dir` exits non-zero with no `*.accessions.txt` written; `from-lexicmap` without MAFFT writes no `lexicmap.*`.

### PC-14 · low · A stray positional token sends it and all later options to the walker passthrough

`src/disassembler/cmd/pipe_reconstruct.py:198` · usability · confirmed (reproduced) · fix effort: small · sources: SWP4-01

**Problem.** The passthrough is `g_pass.add_argument("walker_args", nargs=argparse.REMAINDER, ...)` (`pipe_reconstruct.py:198`; same in `pipe_from_kmindex.py:252`), and `run` accepts it whether or not it began with `--` (`pipe_reconstruct.py:315`, `pipe_from_kmindex.py:265-267`). A stray word, typically a second FASTA after `-q`, swallows every later option without a `--`. The walker then fails with `unknown flag geneB.fasta`; in `from-kmindex` only after the harvest, and if `-o` came later argparse reports a misleading "-o/--out-dir is required".

**Proposed fix.** Split at `--` in `cli.main` before argparse runs, as it already does for `walker`; drop the REMAINDER positional and `set_defaults(walker_args=[])`. In `stage.stage_args`, apply parser defaults before the unknown-option check and before `values`. A weaker alternative: fail in `run` when `walker_args and walker_args[0] != "--"`.

**Tests and docs.** `reconstruct -a SRR1 -q ref ref -o out` gives rc 2, "unrecognized arguments", no out dir; same for `from-kmindex`.

### PC-15 · low · Passthrough and validator ranges disagree with the walker (`-t`, `--prune-frac`, `--min-abund`, repeated flags)

`src/disassembler/cmd/pipe_reconstruct.py:213` · usability · confirmed (reproduced) · fix effort: small · sources: PE-08, SWD-07, P4-09

**Problem.** Three mismatches. (1) Help (`pipe_reconstruct.py:73`, `pipe_from_kmindex.py:135`) and docs (`docs/codon_msa.md:102`, `docs/galaxy.md:42`) say `-t` is "kept at 1 for GraphAligner reproducibility", but the composite neither passes `-t 1` nor lists `-t` in `COMPOSITE_OWNED_WALKER_FLAGS` (`:213-224`); `-- -t 8` goes through. (2) `--prune-frac` uses `fraction`, closed [0,1] (`:112`), while the walker requires [0,1) (`walker/src/main.rs:300`); `--min-abund` uses `positive_float` (`:121`) while the walker accepts 0. (3) Surfaced flags (`-w`, `--fetch-workers`, `--unitigs-cache`, `--offline`, `--max-size-mb`) may be repeated after `--`, and the walker's last value overrides the one the composite validated (`-w 2 -- -w 64` runs 64 workers).

**Proposed fix.** Either append `-t 1` and add `-t`/`--threads` to the owned table, or reword help and docs to "walker default 1; settable via `--`". Add `half_open_fraction` and `non_negative_float` validators for `--prune-frac` and `--min-abund` (in `reconstruct`, `from-kmindex` and `cmd/walker.py`). Add the surfaced flags to the owned table with "use the composite's own X".

**Tests and docs.** Scrub tests for `-- -t 8` and `-- -w 64`; validator tests at 1.0 and 0; update `docs/codon_msa.md:102` and `docs/galaxy.md:42`.

### PC-16 · low · Empty `-a ""` crashes reconstruct and fetch with a TypeError traceback

`src/disassembler/cmd/pipe_reconstruct.py:361` · correctness · confirmed (reproduced) · fix effort: trivial · sources: P4-06

**Problem.** `-a ""` satisfies argparse's required group but is falsy, so `if args.accession:` (`pipe_reconstruct.py:361`; `cmd/fetch.py:59`) takes the `-l` branch with `accession_list=None`. `reconstruct` dies at `' '.join(walker_argv)` (`:387`) and `fetch` inside `subprocess.run`, each with a traceback instead of a usage error. Triggered by an unset shell variable or an empty wrapper parameter.

**Proposed fix.**

```python
# validators.py
def non_empty(s: str) -> str:
    if not s.strip():
        raise argparse.ArgumentTypeError("must not be empty")
    return s
# pipe_reconstruct.py and fetch.py: type=non_empty on -a and -l; `if args.accession is not None:`
```

**Tests and docs.** `main(["reconstruct", "-a", "", ...])` and `main(["fetch", "-a", "", ...])` raise `SystemExit(2)`; a unit test for `non_empty`.

### PC-17 · low · textio's gzip reader never closes the file it opened; truncated gzip ends in a traceback

`src/disassembler/textio.py:31` · reliability · confirmed (reproduced) · fix effort: small · sources: PE-10, P4-11, PB-08, T1-10

**Problem.** `text_stream` wraps gzip as `buf = gzip.GzipFile(fileobj=buf)` (`textio.py:31`); a `GzipFile` given `fileobj=` does not close it, and `open_text` (`:35-42`) returns only the text wrapper, so `with open_text(p)` leaves the raw handle open until garbage collection. This is the root cause of the known unclosed handle at `stratify.py:195`. A truncated `.gz` raises `EOFError`, which the one-line error handlers (`msa.py:192-200`, `sieve.py:193`, `pipe_reconstruct.py:287`) do not catch, so `msa` prints a 20-line traceback. On CPython, refcounting frees descriptors quickly; no output is wrong. Four `sum(1 for line in X.open() ...)` counters (for example `pipe_reconstruct.py:392`, `pipe_from_lexicmap.py:311`) also leak.

**Proposed fix.**

```python
class _GzipReader(gzip.GzipFile):
    def __init__(self, fileobj, owns=False):
        super().__init__(fileobj=fileobj, mode="rb")
        self._owned = fileobj if owns else None
    def _wrap(self, fn, *a):
        try: return fn(*a)
        except (EOFError, zlib.error) as e:
            raise OSError(f"truncated or corrupt gzip data ({e})") from None
    def read(self, size=-1):  return self._wrap(super().read, size)
    def read1(self, size=-1): return self._wrap(super().read1, size)
    def readinto(self, b):    return self._wrap(super().readinto, b)
    def close(self):
        try: super().close()
        finally:
            if self._owned is not None: self._owned.close()
# text_stream(raw, owns=False); open_text passes owns=True
```

Rewrite the counters as `with X.open() as fh: n = sum(...)`.

**Tests and docs.** Assert the raw handle is closed after `with open_text(gz)`; a truncated gzip to `msa` returns 1 with one stderr line.

### PC-18 · low · Inputs are read whole with no decompressed-size cap; small gzip inputs can cost GBs of RAM

`src/disassembler/textio.py:48` · robustness · confirmed (reproduced) · fix effort: small · sources: PX-08

**Problem.** `read_text` returns `fh.read()` (`textio.py:48`) after transparent gunzip, `decode_bytes` (`:53`) does the same for HTTP bodies, `harvest` then `json.loads` the whole text (`cmd/harvest.py:254`), and LexicMap reads with an unbounded `readline()` (`cmd/lexicmap.py:584`, `:810`). Peak memory is about twice the decompressed size, which gzip reaches from a file about 1000 times smaller: a 97 KB gzipped shard peaked at 219 MB RSS in `disassembler harvest --shards-dir`. `timeout=30` on `sess.get` is per read, not a total deadline. Inputs are the user's own files or chosen Galaxy server, so this is amplification, not a security hole.

**Proposed fix.** `read_text(path, max_bytes=...)` reads in chunks and raises `ValueError(f"{path}: inflates past {max_bytes >> 20} MB")`, with a generous default (about 2 GB; measure real shards first) and a `--max-input-mb` option on `harvest`; stream `fetch_shard` with `stream=True`, a `Content-Length` check and a wall-clock deadline; use `readline(MAX_LINE)` in LexicMap.

**Tests and docs.** A gzip of 2 MB of spaces with `max_bytes=1 MB` raises; document the cap in `docs/kmindex.md` and `docs/lexicmap.md`.

### PC-19 · low · Sieve/MSA defaults copied as literals into composites and inert-knob checks

`src/disassembler/cmd/pipe_reconstruct.py:251` · reuse · confirmed (reproduced) · fix effort: small · sources: PS-08

**Problem.** Sieve tau 0.005, sketch 256, MAFFT threads 8 and MSA min-cov 0.70 are literals in the atomics, in every composite parser and in two inert-knob checks (`if args.sieve_threshold != 0.005:`, `pipe_reconstruct.py:251`; `pipe_from_lexicmap.py:252-262`). Composites always pass their own copy, so recalibrating tau in `sieve.py` changes `disassembler sieve` only, and changing a composite default without its check makes every run print a false "set but --sieve is off" warning. An explicit `--sieve-threshold 0.005` without `--sieve` gets no warning today. `reconstruct` also hard-codes the 0.95 `--min-len` rule (`:414`) instead of calling `sieve.auto_min_len`.

**Proposed fix.** Constants `DEFAULT_THRESHOLD`, `DEFAULT_SKETCH`, `AUTO_MIN_LEN_FRACTION`, `min_len_for(records)` in `sieve.py`, and `DEFAULT_MAFFT_THREADS`, `DEFAULT_MIN_COV` in `msa.py`. Composite parsers use `default=None` with help built from the constants; `stage.given(**kv)` drops None values so stages apply their own defaults; inert checks become `is not None`.

**Tests and docs.** `--sieve` without sieve flags logs `--threshold {sieve.DEFAULT_THRESHOLD}`; an explicit `--sieve-sketch 256` without `--sieve` warns; no flags, no warning.

### PC-20 · low · Composites re-declare about 40 stage options by hand; help text has already drifted

`src/disassembler/cmd/pipe_from_kmindex.py:235` · reuse · confirmed (reproduced) · fix effort: medium · sources: PS-07

**Problem.** `from-kmindex` re-declares all 25 of `reconstruct`'s walker, sieve and MSA options and forwards each by name (`pipe_from_kmindex.py:143-237`, `:352-380`); `from-lexicmap` copies 16 of `lexicmap`'s filters (`pipe_from_lexicmap.py:77-139`, `:270-289`). Help differs for 9 of 25 and 8 of 16. One is misleading: `from-kmindex --skip-tree` says "Stop after the walker (and sieve, if enabled)" (`:235-236`), but the sieve never runs under `--skip-tree` (`pipe_reconstruct.py:408`), as `reconstruct`'s help (`:182-183`) says (PC-22). Every option change must be made in two or three places; PC-08 is one result.

**Proposed fix.** In `pipe_reconstruct.py`, add `add_walker_options(g)`, `add_sieve_options(g)`, `add_msa_options(g)`, each returning the tuple of dests it added, with an optional `note=` for composite-specific wording; forward with `stage_args(pipe_reconstruct, accession_list=..., **{d: getattr(args, d) for d in FORWARDED_DESTS})`. Add `lexicmap.add_filter_options(parser)` returning `FILTER_DESTS`. Keep explicit tuples rather than intersecting `_actions`, so `stage_args` still catches a renamed dest.

**Tests and docs.** A parity test (same default, type and help for every shared dest) and a forwarding test with non-default sentinels.

### PC-21 · low · `cli.main` has no top-level error handling: OS errors, truncated gzip, a closed stdout and Ctrl-C end in tracebacks

`src/disassembler/cli.py:75` · reliability · confirmed (reproduced) · fix effort: small · sources: P4-10

**Problem.** `cli.main` dispatches with a bare `rc = args.func(args)` (`cli.py:75`; the `walker` shortcut at `:70` is the same), and nothing in `src/` catches `KeyboardInterrupt`, `BrokenPipeError` or `EOFError`. The composites create `-o` with `out.mkdir(parents=True, exist_ok=True)` (`pipe_reconstruct.py:280`, `pipe_from_lexicmap.py:239`) before any validation, so an `-o` that names a regular file raises `FileExistsError`. The first FASTA read is guarded by `except (OSError, ValueError)` (`pipe_reconstruct.py:287`, `pipe_from_lexicmap.py:222`), which misses the `EOFError` gzip raises for a truncated `-q` (PC-17 has the root cause in `textio`). Progress lines go to stdout, so `disassembler reconstruct ... | head -1` raises `BrokenPipeError` after stage 1: exit 1, `reconstructed.fasta` written, no MSA or tree. Ctrl-C prints a `KeyboardInterrupt` traceback. Every case exits non-zero, so no output is silently wrong; the cost is tracebacks on the terminal and in Galaxy's stderr.

**Proposed fix.** Wrap both dispatch points in `cli.main` (the `walker` shortcut has no `args.cmd`; use the literal `walker` there):

```python
try:
    rc = args.func(args)
except KeyboardInterrupt:
    print(f"disassembler {args.cmd}: interrupted", file=sys.stderr); return 130
except BrokenPipeError:
    os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno()); return 141
except (OSError, EOFError) as e:
    print(f"disassembler {args.cmd}: error: {e}", file=sys.stderr); return 1
```

Re-raise when an environment variable such as `DISASSEMBLER_DEBUG=1` is set, so bug reports can still carry a traceback. Add `EOFError` and `zlib.error` to the `except (OSError, ValueError)` clauses around `fasta_has_record`, or let PC-17's reader turn them into `OSError`. Optionally send progress lines to stderr so a stdout consumer cannot break a run. Applied to a copy of the package, the three reproduced cases print one line each (`disassembler reconstruct: error: [Errno 17] File exists: 'afile'`, rc 1; silent rc 141; `...: error: Compressed file ended before the end-of-stream marker was reached`, rc 1).

**Tests and docs.** `tests/test_cli.py`: `cli.main` with `-o` pointing at a file and with a truncated gzip `-q` returns 1 with a one-line stderr message; a command whose `run` is patched to raise `KeyboardInterrupt` or `BrokenPipeError` returns 130 or 141. `tests/test_msa.py:128` calls `msa.run` directly and is unaffected. List the exit codes wherever exit statuses are documented (see PC-10).

### PC-22 · low · `--skip-tree` help says the sieve still runs; `--sieve --skip-tree` still requires a sieve binary

`src/disassembler/cmd/pipe_from_kmindex.py:236` · docs · confirmed (reproduced) · fix effort: trivial · sources: P4-12, D2-14

**Problem.** Under `--skip-tree` the composites stop right after the walker; the sieve is skipped (`if args.sieve and not args.skip_tree:`, `pipe_reconstruct.py:408`), as `docs/codon_msa.md:98` says. The `from-kmindex` help says otherwise: "Stop after the walker (and sieve, if enabled). Do not build the MSA or tree." (`pipe_from_kmindex.py:236-237`). The `reconstruct` help (`pipe_reconstruct.py:182-183`) contradicts itself, adding "--sieve is also ignored under this flag." to the same sentence. The sieve preflight `if args.sieve and not resolve_sieve(args.sieve_bin):` (`pipe_reconstruct.py:344`) ignores `--skip-tree`, so a stage that will not run can abort the command. With a sieve binary available, `--sieve --skip-tree` exits 0 with only `reconstructed.fasta`, not the `reconstructed.sieved.fasta` the `from-kmindex` help promises (a stderr warning explains it). Without one (a pip or Bioconda install with no `l1_pipeline_serial` on PATH), the run prints the "--sieve is ignored" warning, then dies with "kmer-diversity-sieve not found", rc 1, before the walker starts. No output is wrong.

**Proposed fix.**

```diff
-    if args.sieve and not resolve_sieve(args.sieve_bin):
+    if args.sieve and not args.skip_tree and not resolve_sieve(args.sieve_bin):
```

Use one help string in both parsers: `"Stop after the walker: no sieve, MSA or tree (--sieve is ignored). An MSA, tree or sieved FASTA an earlier run left in -o is removed."` The removal matches `pipe_reconstruct.py:404` and `:432-433`. Sharing the option declarations (PC-20) would stop the two help texts from drifting again.

**Tests and docs.** `tests/test_composites.py`: `--sieve --skip-tree --sieve-bin /nonexistent` with `DISASSEMBLER_SIEVE_BIN` unset gives rc 0, the warning in the log and no `reconstructed.sieved.fasta`; neither parser's help contains "and sieve, if enabled". Add "sieved FASTA" to the list of removed files in `docs/codon_msa.md:98`.
