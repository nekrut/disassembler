# Python: codon MSA and diversity sieve (cmd/msa.py, cmd/sieve.py, extern/kmer-diversity-sieve)

On the inputs the composites normally produce (walker and LexicMap rows with short, unique ids, clean upper-case FASTA, cohorts of a few thousand rows), `msa` and `sieve` give correct output, and the plumbing around the regular-file case (temp file, then `os.replace`) is sound. The fragile parts are scale and edges. The vendored sieve binary has fixed-size tables and an admission gate that the wrapper leaves armed, so large cohorts silently lose divergent lineages or hang, and its header-keyed logic (month buckets, 255-character emit keys) makes the result depend on header text. In `msa`, the only check after MAFFT is a row-length test that `--keeplength` guarantees, so frameshifts, a wrong record set and lower-cased output all pass with exit 0, and temp files are named after user files inside user directories. Most fixes share two patterns: run each external tool in a private temp directory with neutral names and publish only a verified result atomically, and check what the tool reported (REJECT count, record count, codon QC, tool versions) instead of trusting exit status 0.

| ID | Severity | Where | Title | Fix effort |
|---|---|---|---|---|
| PM-01 | high | `src/disassembler/cmd/sieve.py:214` | Sieve wrapper leaves the binary's admission gate armed; divergent records after record 2,000 are dropped | small |
| PM-02 | high | `extern/kmer-diversity-sieve/scripts/l1_pipeline.c:155` | Sieve binary hangs forever once its fixed 2^20-slot LSH table fills (more than 65,536 divergent records) | small |
| PM-03 | medium | `src/disassembler/cmd/msa.py:282` | msa silently hides frameshifts and trimmed insertions; docs promise codon-wise gaps it does not make | small |
| PM-04 | medium | `extern/kmer-diversity-sieve/scripts/l1_pipeline.c:217` | Sieve judges novelty only within a month parsed from `\|YYYY-MM` headers; undocumented | small |
| PM-05 | medium | `src/disassembler/cmd/sieve.py:182` | sieve's shard glob deletes look-alike user files and misses shards for names with glob characters | small |
| PM-06 | medium | `src/disassembler/cmd/msa.py:229` | msa temp names built from user basenames collide; MAFFT fails or a 2-row alignment replaces `--out-aln` | small |
| PM-07 | medium | `src/disassembler/cmd/sieve.py:88` | sieve passes FASTA through unnormalised: a BOM drops record 1, whitespace counts as bad bases | trivial |
| PM-08 | medium | `extern/kmer-diversity-sieve/scripts/l1_pipeline.c:264` | Sieve sketches only the first 65,536 k-mers of a record, silently merging long loci | trivial |
| PM-09 | medium | `src/disassembler/cmd/msa.py:326` | msa renames temp files over special-file outputs: `/dev/null` replaced as root, traceback otherwise | small |
| PM-10 | medium | `src/disassembler/cmd/msa.py:251` | msa aligns every exact-duplicate row with MAFFT; aligning distinct rows once gives identical output | small |
| PM-11 | low | `src/disassembler/cmd/msa.py:279` | Codon MSA is written in lower case; docs claim rows come through "byte for byte" | trivial |
| PM-12 | low | `src/disassembler/cmd/msa.py:248` | msa keeps a minus-strand row it cannot orient; it is aligned as noise with only a warning | small |
| PM-13 | low | `extern/kmer-diversity-sieve/scripts/l1_pipeline.c:226` | Sieve output re-selected by 255-character header: same-header records over-emitted, long headers truncated | small |
| PM-14 | low | `src/disassembler/cmd/sieve.py:224` | A sieve run that promotes nothing truncates the previous output and blames `--threshold` | small |
| PM-15 | low | `src/disassembler/cmd/sieve.py:199` | Sieve dedup records no member ids and counts empty records as exact duplicates | small |
| PM-16 | low | `src/disassembler/cmd/msa.py:289` | Python stages record no tool versions; msa discards the MAFFT and FastTree version output | small |
| PM-17 | low | `src/disassembler/cmd/msa.py:186` | msa runs whatever MAFFT is on PATH; no guard against the known-bad 7.525 | small |
| PM-18 | low | `src/disassembler/cmd/msa.py:160` | msa cannot write an alignment without FastTree; composite `--skip-tree` help misstates the sieve | medium |
| PM-19 | low | `src/disassembler/cmd/sieve.py:219` | Sieve binary killed by a signal is reported as "exited -4" with exit status 252 | trivial |
| PM-20 | low | `src/disassembler/cmd/sieve.py:194` | sieve reports an unwritable output directory as "cannot read the input FASTA" | trivial |
| PM-21 | low | `src/disassembler/cmd/msa.py:297` | msa and other subcommands crash on a non-UTF-8 stdout because of the U+2713 marker | trivial |
| PM-22 | low | `src/disassembler/cmd/msa.py:57` | msa `sanitize_id` is quadratic when many ids share their first 40 characters | trivial |
| PM-23 | low | `src/disassembler/cmd/msa.py:197` | msa holds cohort and alignment in memory and checks rows with per-character Python loops | small |
| PM-24 | low | `src/disassembler/cmd/sieve.py:71` | Sieve dedup holds every distinct sequence in memory; the wrapper re-reads files it could stream | small |
| PM-25 | low | `src/disassembler/cmd/msa.py:242` | msa coverage gate counts raw bases before alignment; flanked partial-gene rows enter MSA and tree | small |
| PM-26 | low | `src/disassembler/cmd/msa.py:54` | msa truncates record ids to 40 characters and renames clean ids order-dependently, with no name map | small |
| PM-27 | low | `src/disassembler/cmd/msa.py:229` | msa and sieve scratch files use predictable PID names and follow planted symlinks | small |

### PM-01 · high · Sieve wrapper leaves the binary's admission gate armed; divergent records after record 2,000 are dropped

`src/disassembler/cmd/sieve.py:214` · science · confirmed (reproduced) · fix effort: small · sources: PE-01, P2-01, PA-02, SC2-02

**Problem.** `sieve_run` always passes `"--warmup", "2000"` (`sieve.py:214`) and never `--outlier-sd`. After 2,000 records the binary rejects a record as an "error" when its count of never-seen k-mers exceeds `mean + outlier_sd*sqrt(var) + 0.5` (`extern/kmer-diversity-sieve/scripts/l1_pipeline.c:209-210`). A rejected record is never sketched, promoted or added to the k-mer count table (`:212`), so the rest of its clade is rejected too. In a cohort of near-identical unique haplotypes the cut-off is about 1/(3·√n) novel divergence: 0.7% at n = 2,100, 0.5% at n = 5,000, whatever the gene length. The module docstring (`sieve.py:63-69`) and `docs/sieve.md:50-56` present byte-exact dedup as the cure, but dedup does nothing when the records are unique. The binary prints `REJECT (error) : N` on stdout; the wrapper's summary line and exit status ignore it.

**Impact.** `disassembler sieve` and `reconstruct` / `from-kmindex` / `from-lexicmap --sieve` at defaults, on more than 2,000 unique records: a divergent lineage placed after record 2,000 vanishes from the MSA and tree with exit 0, and which lineages survive depends on file order (`harvest.py:180` sorts by score, and the walker keeps input order). Cohorts of up to 2,000 records are unaffected. The sieve is opt-in, so this is high, not critical.

**Evidence.** Sieve built with `make serial` from the submodule. 2,100 unique 1,500-nt haplotypes followed by one accession 0.67% away: `REJECT 1`, accession not kept, rc 0; the same accession placed first is kept. A 5-member lineage 2% away after 2,100 records: `[sieve] 2105 -> 2003 haplotypes`, 0 of 5 kept (5 of 5 when placed first). An independent run with a 10-member clade 0.83% away gave `REJECT (error) : 10`, 0 of 10 kept, rc 0.

**Proposed fix.** Never arm the gate (`--warmup` is read with `atoi`, and the test needs `n_mine > warmup`), and check the binary's report:

```diff
+ADMISSION_OFF_WARMUP = 2**31 - 1
+_REJECT_RE = re.compile(r"REJECT \(error\)\s*:\s*(\d+)")
-            "--warmup", "2000",
+            "--warmup", str(ADMISSION_OFF_WARMUP),
-        rc = subprocess.run(cmd).returncode
+        proc = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
+        sys.stdout.write(proc.stdout); rc = proc.returncode
+        m = _REJECT_RE.search(proc.stdout)   # warn if missing; return 1 if count > 0
```

Do not use `--outlier-sd 1e9` instead: it still rejects when the variance is 0. With this change every reproduction keeps the lineage in either order, `REJECT 0`, output for cohorts of up to 2,000 records is byte-identical, and the 18 existing sieve and dedup tests pass. A direct-binary Galaxy sieve tool should pass `--warmup 2147483647` too. The same `subprocess.run` line is changed by PM-19; apply both together.

**Tests and docs.** Add `test_a_late_divergent_lineage_survives_a_large_unique_cohort` to `tests/test_sieve_integration.py` (2,100 unique haplotypes, then a divergent one; fails at HEAD). Make `tests/test_dedup_sieve_equivalence.py:73` pass the wrapper's flags. Rewrite `sieve.py:63-69` and `docs/sieve.md:50-56` ("admission is disabled; dedup is only a speed-up"), and correct `docs/issue_sieve_calibration.md:252-261`.

### PM-02 · high · Sieve binary hangs forever once its fixed 2^20-slot LSH table fills (more than 65,536 divergent records)

`extern/kmer-diversity-sieve/scripts/l1_pipeline.c:155` · reliability · confirmed (reproduced) · fix effort: small · sources: SWP2-01

**Problem.** The promoted set's LSH index is an open-addressing hash table allocated once with 2^20 slots (`pset_init(&ps,m,20)`, `l1_pipeline.c:252`, also `:375` and the threads/MPI merge table at `:496`) and never grown or counted. Each promoted record inserts 16 band keys (`LSH_BANDS`, `:168`). Lookup stops only at an empty slot: `while(p->buk[idx].ids){if(p->buk[idx].key==key)return &p->buk[idx];idx=(idx+1)&p->mask;}` (`:155`, same loop in `bucket_add` at `:151`). Once 2^20 distinct keys are stored, the next lookup of an absent key loops forever.

**Impact.** Any sieve input in which more than 65,536 records pass the completeness gate and are novel with distinct bands: cohorts of about 2% pairwise divergence (HIV-1 pol-like) or many LexicMap partial rows. Closer cohorts store 2.4 to 7 keys per promoted record and fill at roughly 150k to 430k promoted. `disassembler sieve` and every `--sieve` composite print the command line, then sit at 100% CPU and never exit; a Galaxy job runs to its walltime. `docs/galaxy.md:74-80` plans cohorts of 10^5 to 10^6 accessions concatenated before the sieve. The published Bioconda `kmer-diversity-sieve` has the same `l1_pipeline.c`.

**Evidence.** With 65,537 random 200-nt records the run never finishes; with 65,536 it exits 0 in 6 s (`PROMOTED 65536`). gdb on the hung process: `bucket_find ... l1_pipeline.c:155`, `p->n=65536`, `nbuk=1048576`, `idx` cycling. An 80,000-record star phylogeny (1,000 nt, 12 SNVs per branch) also hangs; the threads build hangs at the merge.

```sh
python3 -c "import random
for i in range(65537): print(f'>r{i}'); print(''.join(random.choices('ACGT', k=200)))" > rand.fa
disassembler sieve -i rand.fa -o out.fa --min-len 150 --sieve-bin <make-serial build>/l1_pipeline_serial
```

**Proposed fix.** Upstream, grow the table at half load. Bucket id lists move intact, so results do not change:

```c
typedef struct { /* ... */ bucket_t*buk; uint64_t nbuk,mask,nused; } pset_t;   /* nused=0 in pset_init/pset_free */
static void pset_grow(pset_t*p){ uint64_t nn=p->nbuk*2; bucket_t*nb=calloc(nn,sizeof(bucket_t));
    if(!nb){fprintf(stderr,"pset grow alloc fail\n");abort();}
    for(uint64_t b=0;b<p->nbuk;b++) if(p->buk[b].ids){ uint64_t k=p->buk[b].key,i=(k^(k>>29))&(nn-1);
        while(nb[i].ids) i=(i+1)&(nn-1); nb[i]=p->buk[b]; }
    free(p->buk); p->buk=nb; p->nbuk=nn; p->mask=nn-1; }
/* bucket_add: when the probe ends on an empty slot and 2*(nused+1) > nbuk,
   call pset_grow, re-probe, then nused++ */
```

With this patch the 80,000-record cohort finishes in 72 s and the promoted FASTA is byte-identical to the unpatched binary on inputs that do not hang. Then tag a release, bump Bioconda, re-pin the submodule and raise the recipe floor. Until then, in `sieve_run` after `completeness_failures`, refuse with a one-line message when `n_unique - short - gappy > 65536` (only admitted records reach the table). Chunking and re-sieving is not a fix, since the union can exceed the limit again.

**Tests and docs.** Upstream `tests/unit_test.c`: `test_pset_grows` (4-slot table, 1,000 keys, all found, absent key returns NULL; hangs today). Here: a sieve test with 65,537 random 200-nt records under `subprocess.run(..., timeout=120)`. `docs/sieve.md`: state the limit until the fixed version is pinned; `docs/packaging.md:42`: sieve version floor.

### PM-03 · medium · msa silently hides frameshifts and trimmed insertions; docs promise codon-wise gaps it does not make

`src/disassembler/cmd/msa.py:282` · science · confirmed (reproduced) · fix effort: small · sources: SC2-03, P2-03, PA-06

**Problem.** `msa` runs nucleotide MAFFT with `--auto --keeplength --addfragments` (`msa.py:279-286`). `--keeplength` deletes every base a row has beyond the reference and places gaps where nucleotide scoring puts them, not on codon boundaries. The only check afterwards is `if len(r.seq) != ref_len:` (`msa.py:299`), which `--keeplength` guarantees, under the banner "Validating codon lengths..." (`msa.py:297`). Nothing checks for codon-splitting gaps, removed bases or internal stops. `docs/codon_msa.md:4-5` says every column is a codon position and gaps are "whole codons rather than frameshifts", and `:24-25` says a haplotype one base short "becomes a single terminal gap", which holds only when the missing base is at the end.

**Impact.** Walker primaries carrying a 1-nt indel (assembly homopolymer error or real frameshift allele) reach `msa` from `reconstruct` and `from-kmindex`; the walker gates frameshifts only for alternative alleles (`walker/src/stages/variants.rs:451`). A +1 insertion row comes out identical to the reference, a −1 deletion shifts the rest of the row into the reference frame, in-frame insertions vanish from MSA and tree, and the run ends "Codon alignment complete", exit 0. dN/dS on this file treats frameshifted alleles as intact codons. `from-lexicmap` is mostly protected by its tier-5 frameshift filter.

**Evidence.** A 600-nt CDS with 7 edited rows, `disassembler msa -r ref.fa -i in.fa -a out.aln.fa -t out.nwk` (MAFFT 7.505): rc 0, no warning. `ins1_at300`, `ins3_at300` and `ins12_at150` are identical to the reference (tree branch length 0.0); `del1_at200` has one gap at column 200, splitting codon 67; `del3_at301` splits codons 101-102; an internal stop passes.

**Proposed fix.** Warnings only; the alignment is unchanged. While writing the MAFFT input, record each row's A/C/G/T count by sanitised id. After `identity_ok`, run `codon_qc(row, input_bases)` per row:

```python
trimmed = input_acgt - aligned_acgt                    # bases --keeplength removed
internal_runs = [len(m) for m in re.finditer(r"(?<=[^-])-+(?=[^-])", row)]  # not touching either end
split_codons = [c for c in range(0, len(row), 3) if 0 < row[c:c+3].count("-") < 3]
stops = [c for c in range(0, len(row) - 3, 3) if row[c:c+3].upper() in ("TAA", "TAG", "TGA")]  # last codon excluded
frameshift = trimmed % 3 != 0 or any(r % 3 for r in internal_runs)
```

Write `<out-aln>.qc.tsv` (temp file, then rename) and print one `[warn]` line per class with a count and the first ids. Replace the tautological length loop with an assert or drop it. Optionally add `--drop-frameshifts`. Prototyped: it flags exactly the edited rows above, the alignment is byte-identical, and `test_msa`, `test_composites` and `test_msa_real_tools` pass. The trimmed count also includes overhang past the reference ends; say so in the warning.

**Tests and docs.** A `tests/test_msa_real_tools.py` case with the 7-row fixture asserting the QC flags; unit tests of `codon_qc` on hand-built rows. Rewrite `docs/codon_msa.md:4-5`, `:12`, `:21-25` (nucleotide, reference-anchored, insertions deleted, gaps not forced to codon boundaries) and the `msa.py:3-4` docstring; list the QC file.

### PM-04 · medium · Sieve judges novelty only within a month parsed from `|YYYY-MM` headers; undocumented

`extern/kmer-diversity-sieve/scripts/l1_pipeline.c:217` · science · confirmed (reproduced) · fix effort: small · sources: P2-06, D2-07

**Problem.** The binary folds a month bucket parsed from the header into every LSH band key: `int mb=month_bucket(hdr);` then `if(is_novel(&ps,mb,sk,threshold,R,B)) pset_add(...)` (`l1_pipeline.c:217-218`, key at `:148`). `month_bucket` (`:183-190`) reads the text after the first `|` and returns `year*12+month` if it starts `YYYY-MM`, else −1. A record is therefore compared only with promoted records of the same month. The binary has no flag to turn this off, and the wrapper hands it the user's headers verbatim (`sieve.py:84`). `docs/sieve.md:12-13` says a sequence is promoted only if novel "against everything already promoted"; `disassembler sieve -h` says nothing about headers.

**Impact.** `disassembler sieve` on exports with `acc|2021-03-15|...` headers (NCBI Virus, GISAID, Nextstrain style) keeps near-identical sequences from different months, so the thinning and the MSA/tree size depend on header text. Walker and LexicMap headers contain no `|`; composite routes are unaffected.

**Evidence.** 12 copies of one random 30-kb sequence, each with one distinct SNV, `disassembler sieve -i X.fasta -o out.fasta --min-len 1000`: headers `>s0` give `12 -> 1`; `>s0|2021-03-01` (one month) `12 -> 1`; `>s0|2021-01-15` … `>s11|2021-12-15` `12 -> 12`; `>s0|chr1|x` `12 -> 1`.

**Proposed fix.** Give the binary neutral ids and restore the headers afterwards:

```diff
-def byte_exact_dedup(src, dst) -> tuple[int, int]:
+def byte_exact_dedup(src, dst, headers: list[str] | None = None) -> tuple[int, int]:
-            fout.write(header)
+            if headers is None:
+                fout.write(header)
+            else:
+                fout.write(f">r{len(headers)}\n"); headers.append(header)
 # sieve_run
-            n_in, n_unique = byte_exact_dedup(input_fasta, deduped)
+            orig_headers: list[str] = []
+            n_in, n_unique = byte_exact_dedup(input_fasta, deduped, orig_headers)
-                fh.write(rf.read_text())
+                for line in rf.open():
+                    if line.startswith(">"):
+                        line = orig_headers[int(line[2:].split()[0])]
+                    fh.write(line)
```

All four inputs then give `12 -> 1` with the original header in the output; the 18 existing sieve tests pass. If month partitioning is wanted, add an explicit `--by-month` flag that appends the date field to the neutral id. This same change fixes PM-13; PM-13's variant keeps the date field by default, so pick one policy for both.

**Tests and docs.** `tests/test_sieve_integration.py::test_dated_headers_do_not_partition_novelty` (6 one-SNV variants with `|2021-MM-15 desc` headers give exactly `["s0|2021-01-15 desc"]`; fails today). State in `docs/sieve.md` and the module docstring that header dates are ignored (or document `--by-month`).

### PM-05 · medium · sieve's shard glob deletes look-alike user files and misses shards for names with glob characters

`src/disassembler/cmd/sieve.py:182` · reliability · confirmed (reproduced) · fix effort: small · sources: PX-04, P2-07, PA-04, PE-09

**Problem.** `sieve_run` runs the binary with `--emit-promoted <output without suffix>` in the user's output directory and finds the shards with `sorted(output_fasta.parent.glob(f"{emit_prefix.name}.rank*.fasta"))` (`sieve.py:182`). It unlinks every match before the run, even before the input is read (`:185-186`), concatenates every match into the output (`:225-226`), and unlinks every match again in `finally` (`:231-232`). `*` matches any text, and the stem is not escaped, so `[...]` becomes a character class.

**Impact.** (a) `-o d/cohort.fasta` silently deletes `d/cohort.ranked_by_abundance.fasta`, `d/cohort.rank1_manual.fasta` or a shard from a manual `--emit-promoted d/cohort` run, even when the run then fails on a missing input; a file matching at concatenation time is merged into the output. (b) `-o 'd/cohort[v2].fasta'`: the glob misses the shard, the wrapper writes a 0-byte output, says "promoted 0 haplotypes; try a lower --threshold", exits 1 and leaves the shard behind. (c) Two concurrent runs whose outputs share a stem in one directory sweep and merge each other's shards. Composites use fixed names, so exposure is direct CLI use.

**Evidence.** With `cohort.notes.fasta`, `cohort.rank1_manual.fasta` and `cohort.ranked_by_abundance.fasta` in `a/`, `disassembler sieve -i in.fa -o a/cohort.fasta --min-len 1000` exits 0 and leaves only `cohort.fasta` and `cohort.notes.fasta`. With `-i /no/such.fa` the look-alike is deleted and the run fails. `-o 'b/cohort[v2].fasta'` gives `30 -> 0 haplotypes`, rc 1, a 0-byte output and a 45-kB leftover `cohort[v2].rank0.fasta`.

**Proposed fix.** Run the binary in a private directory so no user text enters a pattern; replace `sieve.py:179-232` with:

```python
work = Path(tempfile.mkdtemp(prefix=".sieve_", dir=output_fasta.parent))
try:
    deduped = work / "in.fasta"; emit_prefix = work / "promoted"
    ...  # dedup, gate report, cmd: unchanged, using these paths
    shards = sorted(work.glob("promoted.rank*.fasta"),
                    key=lambda p: int(p.name[len("promoted.rank"):-len(".fasta")]))
    staged = work / "out.fasta"
    with staged.open("w") as fh:
        for rf in shards: fh.write(rf.read_text())
    os.replace(staged, output_fasta)
finally:
    shutil.rmtree(work, ignore_errors=True)
```

Numeric order also fixes `rank10` sorting before `rank2` in multi-shard builds. Prototyped: siblings survive, the bracket name gives all 30 records with rc 0, nothing is left behind, and the 19 existing sieve tests pass. `glob.escape` plus a `rank\d+` filter is the weaker fix: it still sweeps the user's directory. PM-14 needs the same staged-then-replace write.

**Tests and docs.** `TestSieveRun.test_sibling_files_and_glob_characters_in_the_name_are_safe` in `tests/test_sieve_integration.py` (fails today). Update the comments at `sieve.py:184`, `:188`, `:228`.

### PM-06 · medium · msa temp names built from user basenames collide; MAFFT fails or a 2-row alignment replaces `--out-aln`

`src/disassembler/cmd/msa.py:229` · correctness · confirmed (reproduced) · fix effort: small · sources: PE-04, P2-02, PB-07

**Problem.** `msa` names its three scratch files after the user's files, all in `out_aln.parent`: `.tmp_{in_p.name}.{pid}` (`msa.py:229`), `.tmp_ref_{ref_p.name}.{pid}` (`:230`) and `.tmp_{out_aln.name}.{pid}` (`:234`). Ordinary name choices make two of them the same path, and nothing after MAFFT checks that the output holds the reference plus every retained row (`:296-307` check only length and identity).

**Impact.** Three collisions. (1) `-i raw/cohort.fasta -a aligned/cohort.fasta`: writing MAFFT's output (`:288`) truncates its input; the run ends `MAFFT failed: Check out/.tmp_cohort.fasta.<pid>`, exit 1, with no hint of the cause. (2) `-r pol.fasta -i ref_pol.fasta`: the reference overwrites the cleaned input, MAFFT aligns the reference to itself, the log says `MAFFT aligned 2 records` and `median 1.000`, and a 2-row alignment (`>Reference` twice) replaces any earlier good `--out-aln`; only FastTree's duplicate-name error makes it exit 1. (3) `-a o3/ref_pol.fasta` with `-r pol.fasta`: MAFFT's output truncates the reference, exit 1. Composites use fixed names and are unaffected.

**Evidence.** 900-nt reference, 5 mutated haplotypes, `disassembler msa -r pol.fasta ...`: the control writes 6 taxa, rc 0; the three cases give the outcomes above, and after case 2 the output file holds `>Reference`, `>Reference`.

**Proposed fix.** A private temp directory with fixed names, and a record-count check before publishing:

```diff
+import tempfile
-    clean_input_p = out_aln.parent / f".tmp_{in_p.name}.{os.getpid()}"
-    clean_ref_p = out_aln.parent / f".tmp_ref_{ref_p.name}.{os.getpid()}"
+    tmp_dir = Path(tempfile.mkdtemp(prefix=".tmp_msa_", dir=out_aln.parent))
+    clean_input_p = tmp_dir / "input.fasta"
+    clean_ref_p = tmp_dir / "reference.fasta"
-    tmp_aln = out_aln.with_name(f".tmp_{out_aln.name}.{os.getpid()}")
+    tmp_aln = tmp_dir / "aligned.fasta"
 ...
+        if len(aligned_recs) != retained + 1 or aligned_recs[0].id != "Reference":
+            print(f"Error: MAFFT returned {len(aligned_recs)} records, expected the "
+                  f"reference plus {retained} sequences", file=sys.stderr); return 1
 ...
-        for temp in (clean_input_p, clean_ref_p, tmp_aln):
-            temp.unlink(missing_ok=True)
+        shutil.rmtree(tmp_dir, ignore_errors=True)
```

Keeping `tmp_dir` beside `--out-aln` keeps `os.replace` on one filesystem. Prototyped: all three cases finish with rc 0 and 6 taxa, no `.tmp_*` left, 16 `test_msa*` tests pass. PM-09 moves the same temps for special-file outputs; implement both in one change.

**Tests and docs.** In `tests/test_msa.py` with stub tools: cases (1)-(3) give rc 0 and retained+1 rows; a stub MAFFT that drops a record gives rc 1 and leaves an existing `--out-aln` untouched.

### PM-07 · medium · sieve passes FASTA through unnormalised: a BOM drops record 1, whitespace counts as bad bases

`src/disassembler/cmd/sieve.py:88` · correctness · confirmed (reproduced) · fix effort: trivial · sources: SWP2-02

**Problem.** `textio.text_stream` decodes with `encoding="utf-8"` (`src/disassembler/textio.py:32`), so a UTF-8 byte-order mark survives as U+FEFF, the first header fails `line.startswith(">")` (`sieve.py:90`), and its lines are dropped while `header is None` (`:79-80`): the record vanishes and is not counted in `n_in`. Separately, `_emit` writes lines verbatim (`fout.write(header); fout.writelines(seq_lines)`, `:84-85`). The binary strips only `\n`/`\r` (`l1_pipeline.c:276`) and `n_fraction` counts every other byte as N (`:132`), while the wrapper's own predictor `completeness_failures` strips whitespace (`sieve.py:127`) and predicts no loss.

**Impact.** A BOM-prefixed FASTA loses record 1 with `dedup: 10 -> 10` logged for an 11-record file, rc 0. Wrapped lines with trailing spaces push records with 4% N past the 5% gate; they are dropped with rc 0 and no "completeness gate drops" line, so the loss reads as thinning. Walker and LexicMap output is clean; this hits hand-made inputs to `disassembler sieve` and the Galaxy tool.

**Evidence.** 900-nt records on 60-nt lines, `--min-len 800`: BOM file gives `dedup: 10 -> 10` and the first record is absent from the output; with one trailing space per line, 5 records at 4% N give `LOW-QUAL (partial) : 5`, `10 -> 5`, rc 0; the same file without spaces gives `10 -> 10`.

**Proposed fix.** `textio.py:32`: `encoding="utf-8-sig"` (same result without a BOM; also fixes `msa`, `auto_min_len` and `stratify`). In `byte_exact_dedup`, normalise what the binary sees:

```python
seq = "".join("".join(seq_lines).split())
key = seq.upper()
if seq and key not in seen:
    seen.add(key)
    fout.write(header.rstrip("\r\n") + "\n" + seq + "\n")
...
elif header is None and line.strip():
    raise ValueError(f"{src}: text before the first '>' header")
```

The binary unwraps sequences on emit, so output is unchanged for clean input. Prototyped: BOM file `11 -> 11` with record 1 present; whitespace file `10 -> 10`. Rejecting leading text is a behaviour change for `;` comment lines, which Biopython in `msa` already rejects.

**Tests and docs.** `tests/test_sieve_integration.py`: BOM input keeps and counts record 1; trailing-space input gives the clean file's promoted set; text before the first header raises. Update the `byte_exact_dedup` docstring.

### PM-08 · medium · Sieve sketches only the first 65,536 k-mers of a record, silently merging long loci

`extern/kmer-diversity-sieve/scripts/l1_pipeline.c:264` · correctness · confirmed (reproduced) · fix effort: trivial · sources: P2-05

**Problem.** The per-shard k-mer buffer has fixed capacity, `size_t seqcap=1<<20,hdrcap=1<<12,kcap=1<<16;` (`l1_pipeline.c:264`, MPI path `:386`), and is never grown, while the sequence buffer is (`:277`). `genome_kmers` stores with `if(nout<cap)out[nout++]=can;` (`:116`, `:118`), so k-mers past position ~65,550 are dropped without a message. Admission and novelty see only the first ~65.5 kb of each record; the upstream README (line 17) claims no genome-length assumption.

**Impact.** Records over ~65.5 kb (very long genes such as TTN, large viral genomes, operon-scale queries) that are identical over their first 65.5 kb are judged near-duplicates; one is dropped, rc 0, nothing on stderr. Typical CDS inputs are a few kb, so this is edge-case input.

**Evidence.** `disassembler sieve -i long.fasta -o out.fasta --min-len 1000`: two 100-kb records sharing the first 70 kb give `2 -> 1`; two 70-kb records differing only after 65,600 give `2 -> 1`; two 50-kb records sharing the first 20 kb give `2 -> 2`. With the fix all three give `2 -> 2`; upstream `make unit` passes 23/23.

**Proposed fix.** Upstream, in `PIPE_PROCESS` before the `genome_kmers` call (`:204`), covering serial, threads and MPI:

```c
+    if(seqlen>kcap){ kcap=seqlen; kbuf=realloc(kbuf,kcap*sizeof(uint64_t)); \
+        if(!kbuf){fprintf(stderr,"kbuf alloc fail\n");L1_ABORT();} } \
     size_t nu=genome_kmers(seq,seqlen,k,kbuf,kcap); \
```

Then re-pin the submodule. Until then, have `completeness_failures` also return the longest length and refuse (rc 2) or warn when it exceeds 65,550 nt.

**Tests and docs.** Upstream end-to-end case: two 70-kb records differing only after 66 kb give 2 promoted. Here: the same fixture in `tests/test_sieve_integration.py` (skipped without the binary) and a unit test for the length warning. State the limit in `docs/sieve.md` until the fixed pin lands.

### PM-09 · medium · msa renames temp files over special-file outputs: `/dev/null` replaced as root, traceback otherwise

`src/disassembler/cmd/msa.py:326` · reliability · confirmed (reproduced) · fix effort: small · sources: BL2-01

**Problem.** Every output is written to `.tmp_<name>.<pid>` beside the target and renamed over it: `tmp_tree = out_tree.with_name(...)` (`msa.py:316`) then `os.replace(tmp_tree, out_tree)` (`:326`); the alignment likewise (`:234`, `:307`). A rename swaps the directory entry, so a device, FIFO or symlink output (`/dev/null`, `/dev/stdout`, a named pipe) becomes a regular file. Where the target's directory cannot take new files (`/dev` for a normal user, `/dev/fd` for bash `>(...)`), the temp `open()` raises an uncaught exception after MAFFT has run. FastTree also re-reads `out_aln` (`:318`).

**Impact.** `--out-tree` is required (PM-18), so `--out-tree /dev/null` is a natural way to keep only the alignment. As root (default in a container): exit 0, and `/dev/null` becomes a regular file holding the Newick tree for the rest of that container. As a normal user: `PermissionError` traceback, exit 1, after the alignment is published. `--out-tree >(gzip > t.gz)`: `FileNotFoundError`, exit 1. A FIFO output is replaced and its reader blocks forever. No shipped command line reaches this: composites and the Galaxy wrapper pass regular paths, Galaxy and Singularity run as non-root, and progress already goes to stdout, so `--out-aln /dev/stdout` was never usable. Same code in 0.1.0.

**Evidence.** On a scratch character device (`mknod c 1 3`), a FIFO and a symlink to stdout: device and FIFO outputs become regular files with rc 0 (the FIFO reader times out with 0 bytes); `>(cat > f)` gives `FileNotFoundError: '/dev/fd/.tmp_63.<pid>'`; a device in a read-only directory without root capabilities gives `PermissionError` after MAFFT.

**Proposed fix.** Rename only onto a regular file or a new path; copy into anything else, keeping its temps in `$TMPDIR`. Use `os.stat` and `S_ISREG`: `os.path.realpath('/dev/stdout')` on a pipe resolves to a non-existent `/proc/PID/fd/pipe:[N]`.

```python
def is_special(path):   # symlink, device, FIFO, socket
    try: return path.is_symlink() or not stat.S_ISREG(os.stat(path).st_mode)
    except FileNotFoundError: return False
def publish(tmp, dest):
    if is_special(dest):
        with open(tmp, "rb") as s, open(dest, "wb") as d: shutil.copyfileobj(s, d)
    else: os.replace(tmp, dest)
```

Put the MAFFT input temps, `tmp_aln` and `tmp_tree` in `tempfile.gettempdir()` when the matching output is special; call `publish()` at `:307` and `:326`; feed FastTree from `tmp_aln` (a special `out_aln` cannot be read back); remove all temps in one `finally`. Prototyped (~45-60 lines): all cases exit 0, devices, FIFO and link survive, the FIFO reader gets the tree, and the Python suite passes. Combine with PM-06's temp directory.

**Tests and docs.** `tests/test_msa.py`: `test_a_fifo_output_is_written_into_not_replaced` and `test_a_symlinked_output_keeps_its_link` (both fail today). `docs/codon_msa.md` "What it produces": `--out-tree /dev/null` discards the tree; pipe and device outputs are written into.

### PM-10 · medium · msa aligns every exact-duplicate row with MAFFT; aligning distinct rows once gives identical output

`src/disassembler/cmd/msa.py:251` · efficiency · confirmed (reproduced) · fix effort: small · sources: PF-01

**Problem.** `msa.run` writes every retained row to MAFFT's `--addfragments` file (`fh.write(f">{clean_id}\n{seq}\n")`, `msa.py:251`) without checking for repeats. Under `--addfragments` each fragment is aligned to the reference alone, so every copy of a sequence gets the same row, and MAFFT's cost grows linearly with rows. `from-lexicmap` sends `clean.accessions.fasta` (one row per accession) by default (`pipe_from_lexicmap.py:296`) and `reconstruct` one haplotype per accession (`pipe_reconstruct.py:407`); only the opt-in sieve deduplicates.

**Impact.** Correct output, but N MAFFT runs instead of U for clonal or conserved genes (the docs describe the Zika cohort as "consensus-clonal"). At the measured rate a 14,594-row cohort costs about half an hour of MAFFT time, almost all on copies.

**Evidence.** 3,822-nt reference, 600 rows drawn from 30 variants: 129.7 s user CPU, `MAFFT aligned 601 records`; the 30 distinct rows alone: 8.4 s; each of the 600 rows equals the distinct-run row for its sequence. With the patch: identical alignment (`cmp`), 9.3 s. A 240-row set with 1-6 nt indels, trimmed partials and 15% reverse-complemented copies: alignment and tree identical, 51.6 s to 6.0 s; 43 msa, composite and sieve tests pass.

**Proposed fix.** Key each row on its oriented, upper-cased sequence; send only the first copy to MAFFT; expand back to every id in input order:

```python
rep = rep_of.setdefault(seq, clean_id); rows.append((clean_id, rep))
if rep == clean_id: fh.write(f">{clean_id}\n{seq}\n")
...
distinct = list(SeqIO.parse(tmp_aln, "fasta")); by_id = {r.id: r for r in distinct}
# error if len(distinct) != len(rep_of) + 1 or a rep id is missing
aligned_recs = [distinct[0]] + [copy of by_id[rep] renamed clean_id for clean_id, rep in rows]
# length and identity checks on the expanded rows; rewrite tmp_aln with SeqIO "fasta" if expanded
```

Run the identity check on the expanded rows so the `--min-identity` median keeps weighting by copy number. Biopython's `"fasta"` writer matches MAFFT's 60-column output byte for byte; `"fasta-2line"` does not. Log "U distinct of N". Above 500,000 rows MAFFT `--auto` picks its strategy by row count, so deduplication could change the strategy there.

**Tests and docs.** `tests/test_msa.py`: a fake `mafft` that counts its records; 5 rows with 2 distinct sequences (one reverse-complemented) send 3 records and give 6 output rows in input order. `docs/codon_msa.md` (MAFFT bullet near line 21): duplicates are aligned once and copied back.

### PM-11 · low · Codon MSA is written in lower case; docs claim rows come through "byte for byte"

`src/disassembler/cmd/msa.py:279` · docs · confirmed (reproduced) · fix effort: trivial · sources: SWP2-05

**Problem.** Every stage upper-cases sequence (`msa.py:239`, `:276`), but the MAFFT argv (`msa.py:279-286`) lacks `--preservecase`, and MAFFT lower-cases its output by default. The final `<gene>_codon_msa.fasta` of every composite is the only lower-case file in the pipeline. `docs/codon_msa.md:125-127` says unchanged rows come through "byte for byte"; no row does (case, and 60-column wrapping). The tree is unaffected, but case-sensitive comparisons fail and soft-masking tools see a fully masked alignment. Reproduced with MAFFT 7.505: `head -2 aln.fa` shows `atgcag...`.

**Proposed fix.**

```diff
             "--keeplength",
+            "--preservecase",
             "--thread", str(args.mafft_threads),
```

Reword the docs to "the same sequence, character for character (MAFFT wraps lines at 60 columns)", and keep `msa.py:148-149` consistent.

**Tests and docs.** `tests/test_msa_real_tools.py::test_rows_keep_the_input_case` (fails today); the fake `mafft` in `tests/test_msa.py` only `cat`s, so assert `--preservecase` in its recorded argv.

### PM-12 · low · msa keeps a minus-strand row it cannot orient; it is aligned as noise with only a warning

`src/disassembler/cmd/msa.py:248` · correctness · confirmed (reproduced) · fix effort: small · sources: BL2-08

**Problem.** Orientation uses exact 15-mers sampled every 4th position (`msa.py:72-74`). A row sharing none with either strand ties and is written as given (`elif strand == 0: unoriented += 1`, `:248-249`). MAFFT does not reverse-complement, so a divergent minus-strand row aligns at ~0.43 identity; `identity_ok` warns per row and fails only on the cohort median (`:123-129`). The result is one noise taxon in alignment and tree, exit 0, and a warning ("Check they are the same gene") that names the wrong cause. Needs a non-walker FASTA with rows over ~30% divergent: 7% of reverse rows are unoriented at 30%, 72% at 40%.

**Proposed fix.** On a 15-mer tie, fall back to all 11-mers against precomputed reference 11-mer sets and decide only with a clear margin:

```python
if n_fwd == n_rev and fallback is not None:
    kmers = kmer_set(seq, 1, 11)
    a, b = len(kmers & fallback[0]), len(kmers & fallback[1])
    return 0 if max(a, b) < 2*min(a, b) + 3 else (a > b) - (b > a)
```

Prototyped: the test row is flipped (lowest identity 0.892); unoriented rate falls to 0/300 at 30%, 23/300 at 40%. For rows still tied, optionally align both orientations and keep the better; reword the warning to "possibly the other strand".

**Tests and docs.** Unit test on `orientation` with a reverse-complemented row mutated every 12 nt; mention the fallback in `docs/codon_msa.md:36-41`.

### PM-13 · low · Sieve output re-selected by 255-character header: same-header records over-emitted, long headers truncated

`extern/kmer-diversity-sieve/scripts/l1_pipeline.c:226` · correctness · confirmed (reproduced) · fix effort: small · sources: P2-08

**Problem.** The binary chooses what to emit by header: `FLUSH_EMIT` copies the first 255 characters into `key` and emits every record whose key is in `promoted_hdrs` (`l1_pipeline.c:226-227`, `HW 256` at `:193`), and writes that truncated copy as the output header (`:301`). The wrapper dedups on sequence only (`sieve.py:82`). Two records sharing a header (two concatenated `reconstructed.fasta` runs) are both emitted though only one was promoted (`PROMOTED 3` but `[sieve] 4 -> 4`); `from-lexicmap` headers with long `indels=` lists lose their tail.

**Proposed fix.** The neutral-id rewrite in PM-04 fixes this as well: the binary sees `>r<i>`, and the wrapper streams the shards, writes the saved original header for each `>r<i>`, and counts records while copying. If month buckets are kept (`--by-month`), append the date field after `|` (at most 32 characters) to the neutral id so `month_bucket()` sees the same text. Prototyped: duplicate-header input gives `4 -> 3`; a 311-character header survives intact; 39 sieve, dedup and composite tests pass.

**Tests and docs.** `TestSieveRun`: output count equals the binary's PROMOTED count with duplicate headers; a 300-character header survives. `docs/sieve.md` "Duplicates are removed": the binary sees synthetic ids.

### PM-14 · low · A sieve run that promotes nothing truncates the previous output and blames `--threshold`

`src/disassembler/cmd/sieve.py:224` · usability · confirmed (reproduced) · fix effort: small · sources: SWP2-03

**Problem.** `sieve_run` opens `output_fasta` for writing (`sieve.py:224`) before counting records (`:234`), so a run that promotes nothing exits 1 leaving a 0-byte file over any earlier result. The hint "try a lower --threshold, ..." (`:238`) leads with an impossible cause: `is_novel` on an empty promoted set returns 1 (`l1_pipeline.c:156-162`), so zero promoted means every record failed the completeness gate or has no A/C/G/T 15-mer. Reproduced: a rerun with `--min-len 5000` on 900-nt records empties a 9,060-byte output; `--threshold` 0, 0.99 and 1 all give `10 -> 10`.

**Proposed fix.** Write to a temp beside the output, count `>` lines while copying shards, and `os.replace` only when `n_out > 0` (PM-05's staged write does this). For zero, report from `short` and `gappy`: "promoted 0 of N; <output> left unchanged. X shorter than --min-len M, Y over --max-gap-frac F", and drop the `--threshold` suggestion. This also removes the unclosed `output_fasta.open()` at `:234`.

**Tests and docs.** Sentinel output, `min_len` above every record: rc 1, sentinel unchanged, no `.tmp_*`, no "--threshold" in stderr.

### PM-15 · low · Sieve dedup records no member ids and counts empty records as exact duplicates

`src/disassembler/cmd/sieve.py:199` · provenance · confirmed (reproduced) · fix effort: small · sources: P2-12

**Problem.** `byte_exact_dedup` skips a record when its sequence is empty or already seen (`sieve.py:82`) and keeps no record of either. `sieve_run` logs `({n_in - n_unique} exact duplicates removed, first copy kept)` (`:199-200`), so header-only records count as duplicates (`4 -> 2 records (2 exact duplicates removed)` for two real and two empty records). No output says which accessions were folded into which representative or which unique records were not promoted, so in `reconstruct --sieve` the 24 accessions sharing a haplotype with the kept one are untraceable.

**Proposed fix.** Count empty records separately (`N exact duplicates removed, M empty records dropped`). Keep later duplicate ids per kept record, and after mapping promoted records back (PM-13), write an optional `--clusters FILE` TSV: `representative<TAB>member<TAB>reason`, reason `exact_duplicate` or `not_promoted`. Prototyped as opt-in; existing tests pass. Making it always-on changes the composites' output set (`tests/test_sieve_integration.py:275`, Galaxy outputs).

**Tests and docs.** Empty-record log line; TSV contents for a 5-copy cohort. Document in `docs/sieve.md:50-51`.

### PM-16 · low · Python stages record no tool versions; msa discards the MAFFT and FastTree version output

`src/disassembler/cmd/msa.py:289` · provenance · confirmed (reproduced) · fix effort: small · sources: SWK-01

**Problem.** Only the walker records versions (`walker/src/io/output.rs:201`). `msa` captures MAFFT's stderr and prints it only on failure (`msa.py:289`) and runs FastTree with `-quiet` (`:320`); `sieve` logs the binary path but no version (`sieve.py:218`); `harvest` and `lexicmap` summary JSONs carry no disassembler version. MAFFT is unpinned in the published recipe, so two installs can give different `aln.fasta`/`tree.nwk` with nothing in either output naming the tools. A default `msa` run's stdout and stderr contain no version string.

**Proposed fix.** Add `src/disassembler/provenance.py` with `mafft_version(bin)` (first line of `mafft --version`), `fasttree_version(bin)` (first line of `FastTree -expert`), `sha256(path)` and an atomic JSON `write()`. In `msa.run` print `[*] tools: disassembler X; MAFFT v…; FastTree …` after `find_bin` and write `<out_aln>.provenance.json` (versions, command lines, input hashes). `sieve` and composites write `run_info.json` (disassembler version, sieve path and sha256, walker `--version`); `harvest` and `lexicmap` add `"disassembler": __version__`. Shares the `--version` call with PM-17.

**Tests and docs.** Stub `mafft`/`FastTree` answering version queries; assert the JSON. Document in `docs/codon_msa.md`, `docs/sieve.md`.

### PM-17 · low · msa runs whatever MAFFT is on PATH; no guard against the known-bad 7.525

`src/disassembler/cmd/msa.py:186` · reliability · confirmed (reproduced) · fix effort: small · sources: PK2-06

**Problem.** `mafft_bin = find_bin("mafft")` (`msa.py:186`) is a `shutil.which`; nothing runs `mafft --version`, and the composites' preflight (`pipe_reconstruct.py:340`) checks presence only. `docs/packaging.md:39` records that 7.525 segfaults on the `--addfragments` path `msa` uses, and the published recipe leaves MAFFT unpinned. An environment that resolves 7.525 fails with a raw `MAFFT failed:` log and no hint, in `reconstruct` only after the walk. (The 7.525 crash itself was not run here; a shim reporting 7.525 showed `--version` is never called.)

**Proposed fix.** A deny-list, called after `find_bin` and in the composites' preflight:

```python
BAD_MAFFT = {(7, 525): "segfaults on the --addfragments path msa uses; fixed in 7.526"}
def check_mafft(mafft_bin):
    try:
        r = subprocess.run([mafft_bin, "--version"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return
    m = re.search(r"v(\d+)\.(\d+)", r.stderr + r.stdout)   # MAFFT prints to stderr
    if m and (int(m[1]), int(m[2])) in BAD_MAFFT:
        raise RuntimeError(f"{mafft_bin} is MAFFT {m[1]}.{m[2]}, which "
                           f"{BAD_MAFFT[(int(m[1]), int(m[2]))]}. Install mafft>=7.526.")
```

Unparsable output passes, so stubs and 7.505 keep working.

**Tests and docs.** Stub printing `v7.525` raises; `v7.526` and garbage pass; `reconstruct` dies before walking. Mention in `docs/packaging.md:39`.

### PM-18 · low · msa cannot write an alignment without FastTree; composite `--skip-tree` help misstates the sieve

`src/disassembler/cmd/msa.py:160` · feature-gap/docs · confirmed (reproduced) · fix effort: medium · sources: BL2-02

**Problem.** `--out-tree` is `required=True` (`msa.py:160`) and FastTree is looked up before any work (`:187`), so an alignment-only run needs FastTree installed, and there is no prep-only mode to hand the gated, oriented rows to the IUC MAFFT and FastTree tools (blocking the Galaxy split in `docs/galaxy.md:86-89`). The composites' `--skip-tree` help says "Stop after the walker (and sieve, if enabled)" (`pipe_from_kmindex.py:236`; `pipe_reconstruct.py:182` adds "--sieve is also ignored under this flag", contradicting itself), while `pipe_reconstruct.py:408` (`if args.sieve and not args.skip_tree:`) skips the sieve. `from-lexicmap` has no `--skip-tree`. No wrong output.

**Proposed fix.** (1) `--out-tree` optional plus `--no-tree`, exactly one required; look up FastTree only when a tree is wanted. (2) `--prep-only --out-prep-ref R --out-prep-seqs S`: write `>Reference` and the gated rows, print the equivalent `mafft --auto --keeplength --addfragments S R`, need no tool. (3) `--aligned ALN`: check `Reference` is first and equals `--ref` ungapped, rerun length and identity checks, write `--out-aln`, optional tree. (4) Composite help: "Stop after the walker; the sieve, MSA and tree are skipped". The help fix alone is trivial.

**Tests and docs.** `--no-tree` without FastTree on PATH; `--prep-only` with neither tool; prep, external MAFFT, `--aligned` equals the fused run; update `tests/test_composites.py:95-98`, `docs/codon_msa.md`, `docs/galaxy.md:64,88`.

### PM-19 · low · Sieve binary killed by a signal is reported as "exited -4" with exit status 252

`src/disassembler/cmd/sieve.py:219` · reliability · confirmed (reproduced) · fix effort: trivial · sources: P2-13

**Problem.** `rc = subprocess.run(cmd).returncode` (`sieve.py:219`) is negative for a signal death; it is printed as `sieve: kmer-diversity-sieve exited -4` and returned to `sys.exit`, so the shell sees 252 (SIGILL) or 250 (SIGABRT), and `reconstruct` says `sieve failed (exit -4)`. The sieve Makefile builds with `-march=native`, so a binary moved to another CPU dies with SIGILL. Same pattern for the walker (`walker.py:71`, `:93`; an OOM kill reads `exited -9`). Reproduced with a `--sieve-bin` script that runs `kill -ILL $$`.

**Proposed fix.** One helper used at all three sites:

```python
def report_rc(rc: int, what: str, prog: str) -> int:
    if rc < 0:
        sig = signal.Signals(-rc).name
        hint = (" (built with -march=native on another CPU? rebuild with `make serial` here)"
                if -rc == signal.SIGILL and prog == "sieve" else "")
        print(f"{prog}: {what} was killed by {sig}{hint}", file=sys.stderr)
        return 128 - rc
    if rc:
        print(f"{prog}: {what} exited {rc}", file=sys.stderr)
    return rc
```

Apply together with PM-01, which rewrites the same line.

**Tests and docs.** `kill -ILL $$` stub returns 132 with `SIGILL` in stderr; walker stub with `kill -KILL $$` returns 137. One line in `docs/sieve.md` on `-march=native`.

### PM-20 · low · sieve reports an unwritable output directory as "cannot read the input FASTA"

`src/disassembler/cmd/sieve.py:194` · usability · confirmed (reproduced) · fix effort: trivial · sources: SWP2-04

**Problem.** The handler around `byte_exact_dedup` (`sieve.py:190-195`) prints `sieve: cannot read the input FASTA: {e}` for any `OSError`, but that function also creates the dedup copy beside the output (`dst.open("w")`, `:88`), so a read-only or full output directory is blamed on the input. In `msa` the temp writes (`msa.py:237`, `:275`, `:288`) end in a traceback, and `output_fasta.parent.mkdir` (`sieve.py:177`) raises `FileExistsError` when a path component is a regular file.

**Proposed fix.**

```python
except (OSError, ValueError) as e:
    if isinstance(e, OSError) and e.filename == str(deduped):
        print(f"sieve: cannot write beside the output {output_fasta}: {e}", file=sys.stderr)
    else:
        print(f"sieve: cannot read the input FASTA: {e}", file=sys.stderr)
```

Use a neutral message when `filename is None` (ENOSPC), guard the `mkdir`, and in `msa` wrap the temp writes: `except OSError as e: print(f"Error: cannot write beside --out-aln {out_aln}: {e}"); return 1`.

**Tests and docs.** `-o <tmp>/afile/out.fasta` with `afile` a regular file gives a one-line error and rc 1.

### PM-21 · low · msa and other subcommands crash on a non-UTF-8 stdout because of the U+2713 marker

`src/disassembler/cmd/msa.py:297` · reliability · confirmed (reproduced) · fix effort: trivial · sources: P2-14

**Problem.** Progress lines print `✓` (`msa.py:297`, `:331`; also `harvest.py:191`, `:209`, `stratify.py:285`, `:294`, `lexicmap.py:949`, `:1171-1177`) and `cli.main` sets no error handler. With `PYTHONIOENCODING=ascii`, a legacy single-byte locale, or `LC_ALL=C PYTHONUTF8=0 PYTHONCOERCECLOCALE=0`, the print raises `UnicodeEncodeError` after MAFFT finishes and before `os.replace` (`:307`); the `finally` deletes the alignment, rc 1. Plain `LC_ALL=C` is safe on Python 3.7 and later. Reproduce: `PYTHONIOENCODING=ascii disassembler msa -r tests/fixture_query.fasta -i tests/fixture_reconstructed.fasta -a p.aln -t p.nwk`.

**Proposed fix.** At the top of `cli.main`:

```python
for stream in (sys.stdout, sys.stderr):
    try:
        stream.reconfigure(errors="backslashreplace")
    except (AttributeError, ValueError):
        pass
```

Optionally use ASCII `[ok]` instead of the check mark.

**Tests and docs.** `tests/test_cli.py`: subprocess run with `PYTHONIOENCODING=ascii` exits 0 and writes the output.

### PM-22 · low · msa `sanitize_id` is quadratic when many ids share their first 40 characters

`src/disassembler/cmd/msa.py:57` · efficiency · confirmed (reproduced) · fix effort: trivial · sources: PF-05, PX-10

**Problem.** `sanitize_id` truncates to 40 characters (`msa.py:54`) and on a collision restarts the suffix search at 1 for every record (`while candidate in used:`, `:57`), so n ids with the same 40-character prefix cost about n²/2 probes: 10,000 such ids took 31.8 s in renaming alone. Walker and LexicMap ids are short and unique; this needs user-supplied long ids (e.g. `Severe_acute_respiratory_syndrome_coronavirus_2_isolate_<n>`).

**Proposed fix.** Remember where each base's search stopped; `used` only grows, so the names are identical:

```diff
-def sanitize_id(s: str, used: set) -> str:
+def sanitize_id(s: str, used: set, next_idx: dict | None = None) -> str:
-    idx = 1
+    idx = next_idx.get(base, 1) if next_idx is not None else 1
     while candidate in used: ...
+    if next_idx is not None:
+        next_idx[base] = idx
```

and pass a `next_suffix: dict[str, int] = {}` at `:250`. 10,000 ids: 31.8 s to 0.05 s, same names. Consider also writing an original-to-sanitised id TSV, since truncation is not reversible (PM-26).

**Tests and docs.** 20,000 identical long ids give 20,000 distinct names in under 1 s; equivalence with the old algorithm on a mixed list.

### PM-23 · low · msa holds cohort and alignment in memory and checks rows with per-character Python loops

`src/disassembler/cmd/msa.py:197` · efficiency · confirmed (reproduced) · fix effort: small · sources: PF-02

**Problem.** `input_records = list(SeqIO.parse(fh, "fasta"))` (`msa.py:197`) and `aligned_recs = list(SeqIO.parse(tmp_aln, "fasta"))` (`:296`) materialise everything though each is iterated once; orientation builds a `set` per k-mer (`:87-88`) for a filter that changes nothing; identity loops per column and re-uppercases the reference per row (`:102`). With stub tools, 20,000 rows of 3,822 nt take 21.4 s and 222 MB in Python alone (a streaming version: 5.4 s, 46 MB, identical alignment). Small next to MAFFT itself, but it inflates the memory a Galaxy job must request.

**Proposed fix.** Stream input with `SimpleFastaParser` (id = `title.split(None, 1)[0]`), keeping the `(OSError, ValueError)` handler and the "no records" error before MAFFT; stream `tmp_aln` keeping only `(id, identity)`; build orientation k-mers with a set comprehension; compute identity with the reference precomputed and `bytes.translate`-based counting.

**Tests and docs.** Identity and orientation equivalence on gaps, N, IUPAC and lower case; non-FASTA and empty inputs still fail before MAFFT.

### PM-24 · low · Sieve dedup holds every distinct sequence in memory; the wrapper re-reads files it could stream

`src/disassembler/cmd/sieve.py:71` · efficiency · confirmed (reproduced) · fix effort: small · sources: PF-07

**Problem.** `seen: set[str]` (`sieve.py:71`) holds every distinct upper-cased sequence (about 1 GB at 250,000 gene-length rows, transient before the binary starts); `completeness_failures` re-reads the dedup file (`:201`); shards are loaded whole (`fh.write(rf.read_text())`, `:226`); and the output is re-read through an unclosed handle (`:234`, a ResourceWarning). 50,000 × 3,822-nt records: peak RSS 198 MB, 20 MB with digests, identical output. The binary's own run time dwarfs the extra passes.

**Proposed fix.** Key `seen` on `hashlib.blake2b(seq.encode(), digest_size=16).digest()`; compute short and gappy counts in the same pass and return them; stream shard concatenation with `shutil.copyfileobj` and count `>` lines while copying (shared with PM-14). Keep first-copy-wins order.

**Tests and docs.** The dedup pass returns the same short/gappy counts as `completeness_failures`; the existing equivalence test pins the output.

### PM-25 · low · msa coverage gate counts raw bases before alignment; flanked partial-gene rows enter MSA and tree

`src/disassembler/cmd/msa.py:242` · science · confirmed (reproduced) · fix effort: small · sources: P2-09

**Problem.** The only coverage check runs on the unaligned row: `if real_bases(seq) < ref_len * args.msa_min_cov: continue` (`msa.py:242`). A row that carries part of the gene plus flanking sequence (an assembly contig, a Prokka or Bakta extract) passes on its raw base count, and nothing measures coverage after MAFFT, although `docs/codon_msa.md:30-32` says such a haplotype "is dropped before alignment rather than aligned as a fragment". The outcome depends on the flank side. Flank outside the gene: `--keeplength` trims it and the row enters the MSA and tree as a 20-25% fragment with no message. Flank on the gene side: it is packed into reference columns (coverage 0.99, identity about 0.5); `identity_ok` only warns (`msa.py:123-128`), so the junk row stays in the tree. Reproduced with MAFFT 7.505: contigs holding the last 25% and the first 20% of the gene log `6/6 sequences passed coverage filter` and the tree has 7 taxa, rc 0; with inward flanks the same rows align at identity 0.53 and 0.52, are warned about and kept. Walker haplotypes span the gene and LexicMap rows sit on the reference grid, so the composites are not affected; this is atomic `msa` on other input.

**Proposed fix.** In `run()`, after `identity_ok` and before `os.replace` (`msa.py:307`), measure coverage again on each aligned row and drop the short ones:

```python
kept = [aligned_recs[0]] + [r for r in aligned_recs[1:]
                            if real_bases(str(r.seq).upper()) >= ref_len * args.msa_min_cov]
# log each dropped id with its aligned coverage on stderr; return 1 if only the reference is left
# if anything was dropped: SeqIO.write(kept, tmp_aln, "fasta")
```

Add an opt-in `--drop-low-identity` that also drops rows under `--min-identity`; it is the only check that catches the inward-flank case. Composites build a `Namespace`, so read the flag with `getattr(args, "drop_low_identity", False)`. `--keeplength` deletes insertions, so a row just over the threshold can fall under it after alignment (walker rows with insertions; `from-lexicmap` sets `msa_min_cov = gate - GATE_MARGIN`, `pipe_from_lexicmap.py:362`); give the second gate a small tolerance. Prototyped: the outward-flank input loses 2 rows (`covers 0.25/0.20 after alignment`) and gives 5 taxa; the inward-flank input with `--drop-low-identity` loses 3 rows and gives 4 taxa; `tests/test_msa.py` passes. PM-03's per-row codon QC runs at the same point and can share the pass.

**Tests and docs.** `tests/test_msa_real_tools.py`: a flanked 25% contig is absent from alignment and tree; a unit test of the post-alignment filter on a hand-made aligned FASTA. Correct `docs/codon_msa.md:30-32` and the `--msa-min-cov` help (`msa.py:162-165`: "measured again on the aligned row"); document `--drop-low-identity` and add it to the Galaxy msa wrapper.

### PM-26 · low · msa truncates record ids to 40 characters and renames clean ids order-dependently, with no name map

`src/disassembler/cmd/msa.py:54` · provenance · confirmed (reproduced) · fix effort: small · sources: P2-11

**Problem.** `sanitize_id` maps every character outside `[A-Za-z0-9_]` to `_` (`msa.py:51`), cuts the result to 40 characters (`base = clean[:40]`, `:54`) and resolves collisions with `_1`, `_2` in input order (`:57-58`). The original header is written nowhere. Two headers that differ only after character 40 become `X` and `X_1` and lose the part that told them apart, and an id that was already clean can be renamed while another record takes its name. `docs/codon_msa.md:28-29` says ids are sanitised "and made unique" but not that they are truncated. Reproduced: input ids `SRR1.1`, `SRR1_1` and two `hCoV-19/USA/CA-CDC-QDX12345678/2021|EPI_ISL_100000{1,2}|x` records give rows `SRR1_1` (really `SRR1.1`), `SRR1_1_1` (really `SRR1_1`), `hCoV_19_USA_CA_CDC_QDX12345678_2021_EPI_` and `..._EPI__1` (EPI accessions gone), rc 0. Tree tips can then be traced to input records only by file order. Walker and LexicMap ids are short and clean, so composites are not affected; this hits atomic `msa` on GISAID- or NCBI-style headers.

**Proposed fix.** Replace `sanitize_id` with a two-pass `assign_ids(raw_ids)` over the rows that survive the coverage and orientation steps:

```python
def clean_id(s): return re.sub(r"[^A-Za-z0-9_]", "_", s).strip("_") or "seq"   # no [:40]
# pass 1: an id that is already clean, unique and not "Reference" keeps itself
# pass 2: the rest get clean_id(raw), then "_N" suffixes skipping every name taken in pass 1
```

When any id changed or a row was reverse-complemented, write `<out_aln>.ids.tsv` (`aligned_id`, `original_header`, `reverse_complemented`) and log its path. Prototyped: the long ids become `..._EPI_ISL_1000001_x` and `..._EPI_ISL_1000002_x`, the TSV records `SRR1_1 <- SRR1.1` and `SRR1_1_1 <- SRR1_1`, and `tests/test_msa.py` passes. Keep PM-22's per-base suffix counter in the second pass. Longer tip labels can trouble tools with name caps (old 10-character PHYLIP); FastTree, Newick and HyPhy accept them.

**Tests and docs.** Unit tests of `assign_ids`: a 60-character id is not truncated; a clean unique id keeps its name whatever the order; `Reference` stays reserved; the TSV is written when an id changes. State the rule and the mapping file in `docs/codon_msa.md:28-29`; the Galaxy wrapper may expose the TSV as an optional output.

### PM-27 · low · msa and sieve scratch files use predictable PID names and follow planted symlinks

`src/disassembler/cmd/msa.py:229` · security · confirmed (reproduced) · fix effort: small · sources: PX-06

**Problem.** Scratch files are named from the user's file name plus `os.getpid()` beside the output: `.tmp_{in_p.name}.{pid}`, `.tmp_ref_{ref_p.name}.{pid}`, `.tmp_{out_aln.name}.{pid}` and `.tmp_{out_tree.name}.{pid}` (`msa.py:229-230`, `:234`, `:316`), and sieve's `.tmp_{stem}.deduped.{pid}.fasta` (`sieve.py:189`). They are opened with plain `open(..., "w")` (`msa.py:237`, `:275`, `:288`; `sieve.py:88`), which follows an existing symlink and truncates its target. `os.replace(tmp_aln, out_aln)` (`msa.py:307`) then moves a planted symlink itself into place as the alignment. Another account that can write the output directory can plant symlinks for a range of PIDs: the victim's run overwrites any file the victim can write with `>Reference\n...`, and `--out-aln` is left as a symlink to an attacker-owned file holding the alignment. Reproduced with links at the predictable names in the same process: rc 0, the victim's `authorized_keys` starts `>Reference\nGCTAAAG...`, and `aln.fa` points at the attacker's file. Exposure needs a group-writable shared directory (a sticky world-writable one is covered by `fs.protected_symlinks=1` on most distributions); single-user and Galaxy job directories are not exposed. A long input file name also overflows NAME_MAX in these temp names and crashes with a traceback.

**Proposed fix.** Create scratch files only under a fresh `tempfile.mkdtemp(dir=<output dir>)` (O_EXCL, random name, mode 0700), as PM-06 proposes for the msa inputs and alignment and PM-05 for the sieve's dedup copy and shards, and move `tmp_tree` there too. Files opened with plain `open` inside that directory get the usual umask mode (0644 under umask 022). If per-file `tempfile.mkstemp` is used instead, `os.fchmod(fd, 0o666 & ~umask)` before publishing, or the alignment and tree come out 0600, and create the files inside the `try:` so `finally` removes them. Prototyped with `mkstemp` and the chmod: the planted links' targets are untouched, `aln.fa` is a regular 0644 file, the MSA tests pass, and the NAME_MAX crash is gone. Implementing PM-05, PM-06 and PM-09 together with the tree temp closes this.

**Tests and docs.** `tests/test_msa.py`: monkeypatch `os.getpid`, plant a symlink at each old name, assert the targets are unchanged, `--out-aln` and `--out-tree` are regular files, and their mode equals `0o666 & ~umask`; the same for `sieve`. No docs change.
