# Documentation and reference data (README.md, docs/*.md, queries/, panels/, docs/examples/)

The prose documentation is detailed and mostly accurate about how the tools behave. The single-record marker files under `queries/markers/` are clean, and the 0.1.1 spike-panel repair together with `tests/test_panels.py` closed the defects it targeted. The weak points are the multi-record reference files: two records of the documented 11-strain SARS-CoV-2 spike panel carry the wrong lineage, and two gene groups plus seven N-laden records in `queries/master_planetary_targets.fasta` do not hold what their names say. Nothing in CI checks that a record is the gene or lineage its name claims, and nothing reads `queries/` at all. The documentation has drifted from 0.1.1 in a number of places: the VCF ALT semantics, the provisional (tangle) guard, the sieve threshold, install limits and output naming. Most fixes are data swaps or doc edits; the lasting fix shared by the data entries is a test that checks content against labels (lineage markers, homology to the group reference, ACGT-only), not just file format.

| ID | Severity | Where | Title | Fix effort |
|---|---|---|---|---|
| DD-01 | high | `panels/sarscov2_spike_11strains.fa:586` | Spike panel record `Omicron_EG51` is a 2020 B.1.2 spike (Wuhan + D614G), not EG.5.1 | small |
| DD-02 | medium | `panels/sarscov2_spike_11strains.fa:456` | Spike panel record `Omicron_BA5` is a BA.2 reporter-clone spike, not BA.5 | small |
| DD-03 | medium | `docs/vcf_example.md:28` | VCF docs call ALT the minor allele; the walker also emits the sample's majority allele as ALT | trivial |
| DD-04 | medium | `walker/src/pipeline.rs:577` | The documented Zika worked example is flagged `provisional` by the 0.1.1 tangle guard, and no doc says so | medium |
| DD-05 | medium | `queries/master_planetary_targets.fasta:61` | Master target file: `erg11` and `avr3a` var1..9 are unrelated genes from the reference's own contig | small |
| DD-06 | low | `docs/kmindex.md:36` | kmindex docs omit the required `--format json`; other formats harvest to an empty result with exit 0 | small |
| DD-07 | low | `README.md:60` | Quick start uses repository-relative `queries/` paths that a Bioconda install does not provide | trivial |
| DD-08 | low | `README.md:30` | README install line omits the linux-64-only and Python 3.11/3.12 limits of the Bioconda stack | trivial |
| DD-09 | low | `walker/src/io/output.rs:493` | Per-sample values written to VCF INFO, so the documented `bcftools merge` mixes INFO across samples | small |
| DD-10 | low | `queries/master_planetary_targets.fasta:181` | Seven master-file records contain N or IUPAC codes, two of them 32-42% N; no test reads `queries/` | small |
| DD-11 | low | `README.md:98` | README promises pipe inputs; `reconstruct`, `from-lexicmap` and `lexicmap` reject a piped reference | small |
| DD-12 | low | `docs/vcf_example.md:83` | `vcf_example.md` SQL reads an undefined table `snv`; output block and sample-group text do not match the data | trivial |
| DD-13 | low | `walker/src/main.rs:208` | `--min-shared-kmers 0` runs as 1 but the summary records 0; docs imply an enforced 1-5 range | trivial |
| DD-14 | low | `docs/sieve.md:4` | Sieve docs: wrong O(N^2) motive, misstated effective threshold, 9-vs-10 count, unreproducible tables | small |
| DD-15 | low | `docs/lexicmap.md:122` | `lexicmap.md` says directory-mode outputs use the reference file stem; the code uses the record id | trivial |
| DD-16 | low | `docs/walker.md:288` | `walker.md` Reproducibility describes the removed 1%-of-graph warning and a wrong 10% share | trivial |
| DD-17 | low | `docs/contigs_vs_unitigs.md:87` | `contigs_vs_unitigs.md`: stale 2,708 nt primaries and an understated 15-33x size ratio | trivial |
| DD-18 | low | `docs/lexicmap.md:3` | LexicMap links point to a nonexistent `IndexThePlanet/LexicMap` repository | trivial |
| DD-19 | low | `docs/walker.md:341` | Smaller doc drift: retry count, walker crate list, Galaxy tag naming | trivial |

### DD-01 · high · Spike panel record `Omicron_EG51` is a 2020 B.1.2 spike (Wuhan + D614G), not EG.5.1

`panels/sarscov2_spike_11strains.fa:586` · data · confirmed (reproduced) · fix effort: small · sources: DA-01

**Problem.** The record `>Omicron_EG51 lineage=EG.5.1 source=OR310579.1` is byte-identical to the `[gene=S]` CDS of GenBank OR310579.1. That genome is "SARS-CoV-2/human/USA/Mayo6318/2020", collected 2020-08-24 and classified by NCBI as Pango B.1.2. Its spike is 3,822 nt and differs from `Wuhan_Hu1_ref` by one base, A1841G (spike D614G). It carries none of the EG.5.1 changes (the XBB.1.5 set plus Q52H and F456L); a real EG.5.1 spike is 3,810 nt. The record is unchanged since 0.1.0, and `panels/README.md:63` claims the untouched records are correct.

**Impact.** In multi-anchor mode with this panel, which `docs/walker.md:54` documents, any D614G sample without other panel markers (most B.1.x samples from mid-2020 to early 2021) matches EG51 one base better than Wuhan. EG51 wins outright under every `--primary-metric`, not on a tie. The primary is published as `anchor=Omicron_EG51`, the VCF `CHROM` and `snvs.tsv` `gene` name EG51, and lineage tallies built from `anchor=` (as in `docs/walker.md:76-78`) count 2020 samples as a 2023 variant. Real EG.5.1 samples cannot be tagged EG.5.1 and fall to `Omicron_XBB15`. The bases of the primary are not corrupted; the label is, and users cannot detect it.

**Evidence.** In a checkout:
```
python3 -c "from Bio import SeqIO
r={x.id:str(x.seq) for x in SeqIO.parse('panels/sarscov2_spike_11strains.fa','fasta')}
w,e=r['Wuhan_Hu1_ref'],r['Omicron_EG51']
print(len(w),len(e),[(i+1,a,b) for i,(a,b) in enumerate(zip(w,e)) if a!=b])"
3822 3822 [(1841, 'A', 'G')]
```
NCBI esummary for OR310579.1 gives the 2020 Mayo title and collection date. A release-walker run with a stand-in GraphAligner on synthetic single-unitig samples tagged a Wuhan+D614G sample and a B.1.177-like sample `anchor=Omicron_EG51` (3,820 matched bases against 3,819 for Wuhan), and a real EG.5.1 spike (BS007713.1) `anchor=Omicron_XBB15`.

**Proposed fix.** Replace the record with a clean EG.5.1 spike taken byte for byte from GenBank, wrapped at 60 nt. Two candidates were checked: LC780101.1 (Sagamihara, 2023-08-12, Pango EG.5.1; 3,810 nt, ATG...stop, ACGT only, exactly the XBB.1.5 set plus Q52H and F456L) and BS007377.1 (3,810 nt, same markers, in frame, no ambiguity codes).
```
->Omicron_EG51 lineage=EG.5.1 source=OR310579.1   + 64 sequence lines
+>Omicron_EG51 lineage=EG.5.1 source=LC780101.1   + [gene=S] CDS of LC780101.1 (3,810 nt)
```
With the replacement, the same four synthetic samples are tagged Wuhan, Wuhan, Wuhan and `Omicron_EG51`, and the existing six panel tests still pass. Land together with DD-02.

**Tests and docs.** In `tests/test_panels.py` add (a) `test_records_are_distinct_anchors`: every pair of records must differ by at least 3 x 21 private 21-mers (the shipped Wuhan/EG51 pair has 21; the next closest pair has 105), and (b) a lineage-marker table of 3-4 defining spike changes per record, checked by protein alignment to Wuhan with Biopython's `PairwiseAligner` (this also catches DD-02). Update `panels/README.md:27` (source, 3,810 nt), add EG51 to the "replaced" list, and correct `:63`. Re-measure the cohort percentages in `docs/walker.md:76-78`, which were measured with this panel.

### DD-02 · medium · Spike panel record `Omicron_BA5` is a BA.2 reporter-clone spike, not BA.5

`panels/sarscov2_spike_11strains.fa:456` · data · confirmed (reproduced) · fix effort: small · sources: DA-04

**Problem.** `>Omicron_BA5 lineage=BA.5 source=PV800153.1` names a GenBank record titled "Mutant ... clone SARS-CoV-2-BA.2_nLuc, complete genome", a 2025 BA.2 nanoluciferase reporter clone. The panel copy is faithful; the source choice was wrong. The spike has none of BA.5's defining changes (del69-70, L452R, F486V), keeps BA.2's Q493R (which BA.5 reverted), is 3,813 nt (BA.2 length; BA.5 is 3,807) and differs from the panel's `Omicron_BA2` by 5 nt. Present since 0.1.0; `panels/README.md:25` repeats the label. Related: DD-01.

**Impact.** Multi-anchor runs with the documented panel have two BA.2 anchors and no BA.5 anchor. BA.2 samples split between `Omicron_BA2` and `Omicron_BA5` on a 5-nt margin, and BA.5 samples are tagged almost at random. `anchor=`, VCF `CHROM` and `snvs.tsv` `gene` carry a BA.5 name for BA.2 coordinates; the "26% BA.2" in `docs/walker.md:78` cannot be trusted. Single-reference runs are unaffected.

**Evidence.** The panel `Omicron_BA5` sequence equals the `[gene=S]` CDS fetched for PV800153.1. A protein alignment to Wuhan shows T19I, del24-26, G142D, V213G and Q493R, and no H69-, V70-, L452R or F486V; `Omicron_BA2` and `Omicron_BA5` are both 3,813 nt with 5 mismatches.

**Proposed fix.** Replace the record with the `[gene=S]` CDS of OR326389.1 (IHUCOVID-095226, BA.5.2.1, France 2022-09; 3,807 nt, one terminal stop, no ambiguity codes, the canonical BA.5 set with no private changes):
```
->Omicron_BA5 lineage=BA.5 source=PV800153.1   + 64 sequence lines
+>Omicron_BA5 lineage=BA.5 source=OR326389.1   + [gene=S] CDS of OR326389.1 (3,807 nt)
```
After the swap the six panel tests pass and BA2 vs BA5 differ by exactly H69-, V70-, L452R, F486V and Q493R.

**Tests and docs.** The marker table from DD-01 (BA5 must carry H69-, V70-, L452R, F486V and not Q493R; BA2 must carry Q493R and not L452R) fails on the shipped panel and passes after the swap. Note that a k-mer near-duplicate threshold does not catch this pair. Update `panels/README.md:25`, the "replaced" list, `docs/walker.md:76-78`, and add a changelog note that `anchor=Omicron_BA5` results are not comparable across the change.

### DD-03 · medium · VCF docs call ALT the minor allele; the walker also emits the sample's majority allele as ALT

`docs/vcf_example.md:28` · docs · confirmed (reproduced) · fix effort: trivial · sources: D1-03, IC-12

**Problem.** `docs/vcf_example.md:28` says "`ALT` the minor allele" and `:111` "Each cell is a minor-allele frequency"; `:88` and `:114` say sites "carry a minor allele". The writer emits every sample allele that differs from the reference base:
```rust
if minor != ref_base { alts.push((minor, a.alt_vaf, minor_ka)); }
if major != ref_base { alts.push((major, a.taken_vaf, major_ka)); }   // walker/src/io/output.rs:628
```
Where the reference base is the sample's minority, ALT is the majority allele with AF > 0.5. `snvs.tsv` has the same rows. The VCF header ("Allele frequency") is correct.

**Impact.** This is common when the query comes from a different strain than the cohort, and `--prefer-majority` (default) produces these rows by design. A user following the page averages majority and minority frequencies in the documented `mean_af` query, or loses the site with an `af < 0.5` filter. The committed Zika example never shows it (max AF 0.2857), so the page looks consistent.

**Evidence.** A one-bubble fixture (reference branch ka 30, alternative branch ka 100) run with `--emit-vcf --emit-snv-table` produced `gene 106 . T A . PASS AF=0.7692;AO=100.0;RO=30.0;FORK=0`; with the depths swapped it produced `AF=0.2308`.

**Proposed fix.** Docs only:
- `:28`: "`ALT` lists each allele seen in the sample that differs from the reference base, `AF` its `ka:f` frequency. It is the sample's minority allele when the majority matches the reference; AF above 0.5 means the sample's consensus differs from the reference at this site (RO is then the minority). At a triallelic site both sample alleles can be ALTs."
- `:111` "minor-allele frequency" -> "non-reference allele frequency"; `:88`, `:114` "carry a minor allele" -> "carry a non-reference allele".
- One sentence with the same semantics in `docs/walker.md` near `:136`.

Optional: a `role` column (`minor`/`major`) in `snvs.tsv` and an INFO flag `MAJ`; this changes the schema and needs a Galaxy datatype update.

**Tests and docs.** Add a unit test in `walker/src/io/output.rs`: a fork whose minor allele equals the reference yields one ALT equal to the major allele with AF = `taken_vaf`.

### DD-04 · medium · The documented Zika worked example is flagged `provisional` by the 0.1.1 tangle guard, and no doc says so

`walker/src/pipeline.rs:577` · docs · confirmed (reproduced) · fix effort: medium · sources: D1-04

**Problem.** 0.1.1 marks a result `provisional` when the 5-hop neighbourhood holds more than 3.0 unitigs per reference base, counted on all unitigs before pruning (`walker/src/stages/extract.rs:631`, guard at `walker/src/config.rs:469`). The flagship examples (Zika NS5 in `docs/vcf_example.md:11`, `docs/contigs_vs_unitigs.md:47-61`, `docs/issue_intra_accession_variance_lost.md:161`) sit far above that line yet are presented as clean reconstructions, while `docs/walker.md:121-128` and `walker/src/config.rs:292-296` say easy cohorts sit far below the cut.

**Impact.** Running the documented command on SRR6251967 prints the WARNING and sets `"provisional": true`. Deep viral samples (above roughly 1,000x) are likely all flagged, because sequencing-error unitigs dominate the neighbourhood. Users either distrust sound results or learn to ignore the flag. No wrong sequence is produced.

**Evidence.** On the real SRR6251967 unitigs, the walker printed `the 5-hop neighbourhood holds 203.5 unitigs per reference base (551400 unitigs over 2709 bases, 21.62% of the graph), above the 3.0 guard` and wrote `"provisional": true`. An instrumented count gave 19.8 unitigs per base above the 1% prune floor (ka >= 250) and 2.17 per base above `min_vaf x depth` (ka >= 2,500).

**Proposed fix.**
1. Docs now: state in `vcf_example.md`, `contigs_vs_unitigs.md` and `issue_intra_accession_variance_lost.md` §5 that these runs are `provisional: true` and why (depth-driven error tips); say in `walker.md` what to do with a provisional result; correct `walker.md:121-128` and `config.rs:292-296`.
2. Recalibrate: counting after pruning is not enough (19.8 per base here). Count neighbourhood unitigs with `ka >= max(min_abund, min_vaf * route_depth)`, using the route depth that `path_depth` already computes; this gives 2.17 per base here but moves the guard after alignment. Re-validate on the SARS-CoV-2 calibration set (89.7% vs 51.3% complete-ORF split) plus Zika, Ebola and HIV accessions before changing the default.

**Tests and docs.** Unit test in `extract.rs`: many low-ka tips around a high-depth path are not flagged; a tangle of high-ka branches is. Related: DD-16, DD-17.

### DD-05 · medium · Master target file: `erg11` and `avr3a` var1..9 are unrelated genes from the reference's own contig

`queries/master_planetary_targets.fasta:61` · data · confirmed (reproduced) · fix effort: small · sources: DA-02, DA-03

**Problem.** Records follow the `<gene>_reference` + `<gene>_varN` strain-panel convention of `panels/README.md`, but two groups are not alleles. `candida_auris_erg11_var1..9` (`:61-77`) are other CDSs of NC_140807.1 (an F-box protein, a tRNA adenylyltransferase, an NADH kinase, Bud22p and others); ERG11 itself is B9J08_03698. `p_infestans_avr3a_var1..9` (`:85-101`) are PITG_14293-14328 from NW_003303725.1; AVR3a is PITG_14371. None shares a 21-mer with its reference; protein relatedness is 0.01-0.06 against 0.26 or more in every other group.

**Impact.** A user who cuts out a group as a multi-anchor panel (`grep -A1 '^>candida_auris_erg11_'`) gets the unrelated gene as primary for exactly the samples the group exists to find: any ERG11 or AVR3a differing from the reference (azole-resistance Y132F, K143R; AVR3a-KI) loses to an exact-matching unrelated gene under the default `matched` metric. Output is `anchor=candida_auris_erg11_var1` with no warning; `msa` only warns. Medium rather than high because no command, test, wrapper or package uses the file; `queries/markers/` is correct.

**Evidence.** A release-walker run with a stand-in ungapped aligner on synthetic unitigs: a sample with ERG11 Y132F produced `>SRRFAKE1_primary cov=1.000 pident=100.0 ... anchor=candida_auris_erg11_var1` (1,626 nt, the B9J08_03025 CDS); an AVR3a-KI sample produced `anchor=p_infestans_avr3a_var1`. The record headers themselves show `[locus_tag=B9J08_03025]`, `[protein=F-box protein]` and so on.

**Proposed fix.** Minimal: delete `:61-78` and `:85-102`, leaving each group reference-only (the file drops to 142 records). Alternative: replace the 18 records with real alleles copied from GenBank `fasta_cds_na` with accessions in the header, e.g. ERG11 PV630305.1 (Y132F), PV630302.1 (K143R), OK564647.1; AVR3a KF154428.1 (KI), KP317588.1; a maintainer must curate the choice. Optional: in multi-anchor mode, warn about a reference that shares no k-mer with any other reference.

**Tests and docs.** Add `tests/test_master_panel_homology.py`: each `<g>_varN` protein must have at least 0.25 identity with `<g>_reference` (Biopython `PairwiseAligner`, BLOSUM62); it fails on exactly these 18 records, and the mcr-3..10 homologs pass at 0.33-0.37. Record the reason for the threshold in the docstring. Describe the file in `README.md:75-76` and fix the record count, also quoted at `walker/src/stages/align.rs:128`. Related: DD-10.

### DD-06 · low · kmindex docs omit the required `--format json`; other formats harvest to an empty result with exit 0

`docs/kmindex.md:36` · docs · confirmed (reproduced) · fix effort: small · sources: P3-14

**Problem.** `docs/kmindex.md` never says which kmindex output `harvest` reads. `_merge_payload` (`src/disassembler/cmd/harvest.py:150-166`) understands only `--format json` (`{index: {query: {sample: score}}}`); `json_vec` and single-query `jsonl` entries are skipped silently, and `harvest` prints `[✓] Harvest complete: 0 unique SRA runs across 0 targets.` with exit 0. Also: `:36` "every accession the query returned" is wrong, because `harvest.py:165` never stores a 0.0 score; `:41-42` omits the name rewriting at `harvest.py:75-77`; `:60` and `:129` disagree (5.2 vs 5.1 million).

**Evidence.** A `json_vec` shard `{"idx1":{"zika_ns5":{"SRR1":{"R":0.9,"P":[1,1,0]}}}}` passed to `disassembler harvest --shards-dir D -o O` gave the empty "✓" line and exit 0.

**Proposed fix.** Document "Shards must be kmindex `--format json`; json_vec, jsonl and matrix are not read", ":36 every accession with a nonzero score", the name rules, and 5.2 million. In `_merge_payload`, count entries skipped for shape and entries merged in this payload; if `skipped > 0 and merged == 0`, raise `ValueError("not kmindex --format json ...")`, which callers already count as a failed shard, so an all-wrong-format run exits 1.

**Tests and docs.** In `tests/test_harvest.py`: a `json_vec` shard and a one-line `jsonl` shard must fail with exit 1; a zero score is absent from `scores.json`.

### DD-07 · low · Quick start uses repository-relative `queries/` paths that a Bioconda install does not provide

`README.md:60` · docs · confirmed (reproduced) · fix effort: trivial · sources: SWK-04

**Problem.** `pyproject.toml:31-32` packages only `src/`, and the Bioconda recipe runs only `pip install .`, so the package holds the Python modules and no `queries/` or `panels/`. Every Quick start command (`README.md:60`, `:64`, `:69`, `:73`), `README.md:75` and the `lexicmap --help` epilog use `queries/...` paths without saying they need a source checkout. After `conda install`, the first Quick start line fails with `cannot read query queries/query_hiv1_pol.fasta: No such file or directory`, exit 1, and leaves an empty `out/`.

**Proposed fix.** Add at the top of Quick start:
```markdown
Paths under `queries/` are relative to a source checkout. After a Bioconda install, fetch them first:

    curl -L https://github.com/nekrut/disassembler/archive/refs/tags/v0.1.1.tar.gz \
      | tar -xz --strip-components=1 disassembler-0.1.1/queries disassembler-0.1.1/panels
```
This yields `queries/` and `panels/` byte-identical to the tag. Optionally ship them in the recipe under `$PREFIX/share/disassembler/`.

**Tests and docs.** README, `docs/packaging.md`, the `lexicmap.py` epilog; add the tag bump to the release checklist.

### DD-08 · low · README install line omits the linux-64-only and Python 3.11/3.12 limits of the Bioconda stack

`README.md:30` · docs · confirmed (by reading) · fix effort: trivial · sources: D1-13

**Problem.** `README.md:30-33` says "**Bioconda** (Linux)" and `conda install -c conda-forge -c bioconda disassembler`. logan-walker 0.1.1 is linux-64 only and depends on `graphaligner 1.0.17b`, whose builds need `libzlib <1.3`, while conda-forge Python 3.13 (and current 3.12 builds) need `libzlib >=1.3.1`. disassembler declares only `python >=3.11`. Installing into a Python 3.13 environment or on linux-aarch64 fails to solve, or proposes downgrades, with no hint. Package metadata on api.anaconda.org shows these pins.

**Proposed fix.**
```diff
-**Bioconda** (Linux):
+**Bioconda** (x86-64 Linux only; Python 3.11 or 3.12):
-conda install -c conda-forge -c bioconda disassembler
+conda create -n disassembler -c conda-forge -c bioconda disassembler
```
Add one sentence on the GraphAligner 1.0.17b pin; set `python >=3.11,<3.13` in the next recipe bump.

**Tests and docs.** Same note in `docs/packaging.md` "Platforms"; optionally a CI `conda create --dry-run` job.

### DD-09 · low · Per-sample values written to VCF INFO, so the documented `bcftools merge` mixes INFO across samples

`walker/src/io/output.rs:493` · data · confirmed (by reading) · fix effort: small · sources: D1-12

**Problem.** Each per-accession VCF line carries the sample's AF, AO, RO and FORK in INFO as well as in FORMAT (`output.rs:493`); FORK is INFO-only. `bcftools merge` (the recipe at `docs/vcf_example.md:99` and the Galaxy route at `docs/galaxy.md:66`) keeps one INFO value per key from different inputs. In the shipped `docs/examples/zika_ns5/zika_ns5.cohort.merged.vcf`, all 12 multi-sample sites mix samples; at POS 517, `FORK` is SRR6251967's unitig while `AF=0.1327;AO=3900` are SRR6251993's. FORMAT columns and `snvs.tsv` are correct. The committed file's `##bcftools_mergeCommand` (line 14) also embeds the author's local home and scratch paths.

**Proposed fix.**
```diff
-"...PASS\tAF={};AO={};RO={:.1};FORK={}\tGT:AF:AO:RO\t.:{}:{}:{:.1}"
+"...PASS\t.\tGT:AF:AO:RO:FORK\t.:{}:{}:{:.1}:{}"
```
Replace the INFO header lines with `##FORMAT=<ID=FORK,Number=1,Type=Integer,...>`. Regenerate the Zika example with `bcftools merge --no-version`. If INFO must stay, add `bcftools annotate -x INFO` to the recipe and a "read FORMAT only" note.

**Tests and docs.** Unit test in `output.rs` `vcf_tests` (INFO is `.`, FORMAT is `GT:AF:AO:RO:FORK`); update `vcf_example.md:97-105`, `galaxy.md:66`; changelog note for the format change.

### DD-10 · low · Seven master-file records contain N or IUPAC codes, two of them 32-42% N; no test reads `queries/`

`queries/master_planetary_targets.fasta:181` · data · confirmed (reproduced) · fix effort: small · sources: BL2-14, DA-05

**Problem.** Seven of the 160 records have non-ACGT bases: `zika_ns5_var3` (line 181, 659 N in 2,709 nt), `ebola_gp_var2` (line 287, 861 N + 1 Y in 2,031 nt), `rabies_g_var4` (3 N), `rabies_g_var2` (1 N), and single IUPAC codes in `h5n1_ha_var6`, `dengue_1_ns5_var8` and `westnile_ns5_var5`. The 0.1.1 regression test reads only `panels/` (`PANEL_DIR = REPO_ROOT / "panels"`, `tests/test_panels.py:23`). The walker skips k-mers touching ambiguity codes and warns, so the cost is lost anchoring power; under `--primary-metric pident` a half-length walk over the ACGT part of `ebola_gp_var2` can win. Related: DD-05.

**Proposed fix.** Drop `zika_ns5_var3` and `ebola_gp_var2`; replace the five single-code records with unambiguous source CDSs or drop them. Do not edit bases by hand (the 0.1.0 Delta frame shift came from that).

**Tests and docs.** In `tests/test_panels.py`, loop over `QUERY_DIR.rglob("*.fa")` and `*.fasta`: non-empty, `len % 3 == 0`, ACGT only, unique ids, no internal stop, non-empty glob (no ATG/stop check: many records are mat_peptides). Describe the file in `README.md:75-76`.

### DD-11 · low · README promises pipe inputs; `reconstruct`, `from-lexicmap` and `lexicmap` reject a piped reference

`README.md:98` · correctness · confirmed (reproduced) · fix effort: small · sources: SWD-05, SWR1-02

**Problem.** `README.md:96-99` says every input, FASTA references included, may be a pipe such as `<(zcat hits.tsv.gz)`. `pipe_reconstruct.py:283`, `pipe_from_lexicmap.py:217` and `lexicmap.py:204` test `is_file()`, so `disassembler reconstruct -q <(cat q.fasta) ...` prints `query not found: /dev/fd/63`, exit 1. Relaxing the test alone fails differently, because each composite reads the reference several times and a pipe can be read once. `msa -r/-i` handles pipes. The walker's own multi-read of `-q` is a separate defect.

**Proposed fix.** Materialise non-regular references once, before anything reads them:
```python
def materialize(path, dest_dir, name) -> Path:      # textio.py
    p = Path(path)
    if p.is_file(): return p
    if not p.exists(): raise FileNotFoundError(f"not found: {p}")
    if p.is_dir(): raise IsADirectoryError(f"is a directory, not a file: {p}")
    dest = Path(dest_dir); dest.mkdir(parents=True, exist_ok=True)
    with open(p, "rb") as s, open(dest / name, "wb") as d: shutil.copyfileobj(s, d)
    return dest / name
```
Use it in `pipe_reconstruct.py:281-284` and `:296-297`, and `pipe_from_lexicmap.py:217`; in `lexicmap.py:204` use `not path.exists() or path.is_dir()`.

**Tests and docs.** FIFO tests (`os.mkfifo` fed from a thread) in `tests/test_composites.py`. Until the walker side is fixed, narrow `README.md:96-99` to inputs read once.

### DD-12 · low · `vcf_example.md` SQL reads an undefined table `snv`; output block and sample-group text do not match the data

`docs/vcf_example.md:83` · docs · confirmed (reproduced) · fix effort: trivial · sources: SWD-08

**Problem.** The co-occurrence query reads `FROM snv` (`:83-84`), a table the page never creates (the other query reads `read_csv_auto('snvs.tsv', ...)` at `:66`), so DuckDB reports a missing table. The output block at `:71-76` shows columns the query does not return (it returns pos, ref, alt, n_samples, mean_af, min_af, max_af). `:117-118` calls the second site group "largely disjoint", but SRR6251967, 973 and 979 carry sites of both groups (recomputed from `docs/examples/zika_ns5/zika_ns5.snvs.tsv`).

**Proposed fix.** Add `CREATE TABLE snv AS SELECT * FROM read_csv_auto('snvs.tsv', sep='\t');` before both queries (a `CREATE VIEW` over `shard_*/snvs.tsv` for `:57`); regenerate the output block; reword `:117-118`: "Sites 517, 871, 1762 and 2153 are in SRR6251967, 6251992 and 6251993, and 2064 in 6251973, 979, 992 and 993; SRR6251967, 973 and 979 also carry 2225/2525, so the groups overlap." Recheck `:112`.

### DD-13 · low · `--min-shared-kmers 0` runs as 1 but the summary records 0; docs imply an enforced 1-5 range

`walker/src/main.rs:208` · cli-validation · confirmed (reproduced) · fix effort: trivial · sources: D1-09, IC-12

**Problem.** `main.rs:208` accepts any `u32`; `walker/src/stages/anchor.rs:263` clamps silently (`cfg.min_shared_kmers.max(1)`), and the summary writes the raw value (`walker/src/io/output.rs:219`), so a run at 0 is recorded as 0 but ran at 1. `--help` "validated 1..=5" (`main.rs:75`) and `docs/walker.md:195` "1 to 5" read as an enforced range, though "validated" means tested; values above 5 are legitimate and accepted.

**Proposed fix.**
```diff
-  args.cfg.min_shared_kmers = value()?.parse().map_err(...)?;
+  let n = at_least_one("--min-shared-kmers", &value()?)?;
+  args.cfg.min_shared_kmers = u32::try_from(n).map_err(|e| Usage(format!("--min-shared-kmers: {e}")))?;
-  --min-shared-kmers <N>    Anchor threshold, validated 1..=5 [default: 3]
+  --min-shared-kmers <N>    Anchor threshold, at least 1 [default: 3; outputs tested stable over 1..=5]
```
Do not cap at 5.

**Tests and docs.** `parse_args` test: `0` is a usage error, `7` is accepted. `docs/walker.md:195`: "at least 1; haplotypes were identical across 1 to 5 in validation, and some shift at 10".

### DD-14 · low · Sieve docs: wrong O(N^2) motive, misstated effective threshold, 9-vs-10 count, unreproducible tables

`docs/sieve.md:4` · docs · confirmed (reproduced) · fix effort: small · sources: SC2-09, D2-11

**Problem.** (1) `:4` calls the MSA and tree O(N^2); `msa` uses `mafft --addfragments`, which is linear in N (measured 93 s at N=150, 187 s at N=300). (2) `:103` "`0.005` sits near 0.017%": the binary promotes a sequence only when it differs in `ceil(tau*m)` of m sketch slots (`l1_pipeline.c:131,161`), so at m=256 tau 0.005 acts as 0.0078 (about 0.026%); runs at tau 0.004, 0.005 and 0.0078 gave identical FASTAs. (3) `:115` reports 9 promoted for Zika, `docs/issue_sieve_calibration.md:293` reports 10. (4) RF = 0 compares against the pruned full tree and does not measure loss. (5) The reproduction inputs are untracked. (6) The `:46` sentence about `-q` records belongs with `--min-len`.

**Proposed fix.** Reword `:4` ("MSA time grows linearly with N ..."); state the effective threshold `ceil(tau*M)/M` and optionally log it from `cmd/sieve.py`:
```python
eff = math.ceil(threshold * sketch - 1e-9) / sketch
print(f"[*] effective novelty threshold {eff:.4f}")
```
Re-run Zika and fix 9 or 10; commit the accession list and scripts under `bench/sieve_calibration/` or mark the tables as not reproducible; move the `:46` sentence.

**Tests and docs.** If the log line is added, a test that tau 0.005 with sketch 256 reports 0.0078.

### DD-15 · low · `lexicmap.md` says directory-mode outputs use the reference file stem; the code uses the record id

`docs/lexicmap.md:122` · docs · confirmed (reproduced) · fix effort: trivial · sources: D2-05, IC-12

**Problem.** `docs/lexicmap.md:122` ("writes `<dir>/<reference stem>.*`") and `:181` (default `msa/<reference stem>`) disagree with the code, which uses the first word of the header (`ref_name = header.split()[0] if header else path.stem`, `src/disassembler/cmd/lexicmap.py:237`). `-q queries/query_hiv1_pol.fasta -o out/` writes `out/hiv1_pol.*`, not `out/query_hiv1_pol.*`. `--help` is already correct; the docs drifted.

**Proposed fix.** Fix the docs, not the code (changing the code would rename every 0.1.x user's outputs):
```diff
-a path ending in `/`, writes `<dir>/<reference stem>.*`, ...
+a path ending in `/`, writes `<dir>/<record id>.*`, where the record id is the first word of the
+reference's header (`hiv1_pol` for `queries/query_hiv1_pol.fasta`), ...
-| `-o, --out-prefix <PATH>` | `msa/<reference stem>` | ...
+| `-o, --out-prefix <PATH>` | `msa/<record id>` | ...
```
If the separate fix that sanitises `/` and `\` in record ids lands, mention it here.

**Tests and docs.** Optional test: `-o DIR/` with `query_hiv1_pol.fasta` (header `>hiv1_pol`) produces `DIR/hiv1_pol.clean.msa.fasta`.

### DD-16 · low · `walker.md` Reproducibility describes the removed 1%-of-graph warning and a wrong 10% share

`docs/walker.md:288` · docs · confirmed (reproduced) · fix effort: trivial · sources: D1-06, IC-12

**Problem.** `docs/walker.md:286-289` says SRR8209082's subgraph is "about 10% of the whole graph" and "The walker warns whenever a subgraph exceeds 1% of the graph". 0.1.1 warns only on the unitigs-per-base guard (`walker/src/stages/extract.rs:631`), as the same page says at `:121-128`; a fixture whose subgraph is 100% of the graph ran without a warning. The golden fixture records 26,323 of 81,534 unitigs (a third), and 42% unpruned.

**Proposed fix.** Replace `:286-289` with:
```markdown
byte. `SRR8209082` (HIV-1 *pol*) is reported but not gated: its 5-hop neighbourhood is
about 42% of the graph unpruned (34,582 of 81,534 unitigs; the golden fixture's pruned
subgraph is 26,323, a third), against 0.033% or less for the gated accessions, and its
unpruned route is not reproducible across GraphAligner versions. The walker marks a result
`provisional` and warns when the neighbourhood holds more than 3 unitigs per reference base
(see Outputs). See [`../bench/README.md`](../bench/README.md).
```

### DD-17 · low · `contigs_vs_unitigs.md`: stale 2,708 nt primaries and an understated 15-33x size ratio

`docs/contigs_vs_unitigs.md:87` · docs · confirmed (reproduced) · fix effort: trivial · sources: D1-11

**Problem.** The table (`:47-55`) and `:60` ("the full 2,708 nt") describe a 2,709 nt reference (`queries/markers/zika_ns5.fasta`); 2,708 is the GraphAligner 1.0.19/1.0.20 length, which drops the final base, while the pinned 1.0.17b gives 2,709. `:87` "Contigs are 15-33x smaller ... (2-7 MB vs hundreds of MB here)" does not fit the nine Zika samples: `curl -sSI` Content-Length on `logan-pub` `u/` vs `c/` gives 0.4-1.3 MiB vs 28-68 MiB, 39-83x; 13-20x holds only for the ERG11 controls. Related: DD-04.

**Proposed fix.** `:60` -> "reconstructed the full-length gene as one connected sequence"; footnote the table "lengths measured under GraphAligner 1.0.20; 1.0.17b gives 2,709 nt" (same on `issue_intra_accession_variance_lost.md:158-162`); `:87` -> "Compressed, the contig files are 40-80x smaller than the unitig files for these Zika samples (0.4-1.3 MiB against 28-68 MiB) and 13-20x for the ERG11 controls."

### DD-18 · low · LexicMap links point to a nonexistent `IndexThePlanet/LexicMap` repository

`docs/lexicmap.md:3` · docs · confirmed (reproduced) · fix effort: trivial · sources: D2-13

**Problem.** `README.md:19` and `docs/lexicmap.md:3` link `https://github.com/IndexThePlanet/LexicMap`, which does not exist (`raw.githubusercontent.com/IndexThePlanet/LexicMap/main/README.md` returns 404). LexicMap lives at `shenwei356/LexicMap`, the `home` Bioconda records. Users following the link to build the index hit a 404.

**Proposed fix.**
```diff
-[LexicMap](https://github.com/IndexThePlanet/LexicMap)
+[LexicMap](https://github.com/shenwei356/LexicMap)
```
in both files. Optionally add a non-blocking link-check CI job (e.g. `lychee README.md docs`).

### DD-19 · low · Smaller doc drift: retry count, walker crate list, Galaxy tag naming

`docs/walker.md:341` · docs · confirmed (by reading) · fix effort: trivial · sources: IC-12

**Problem.**
- `docs/walker.md:340-341` "server errors are retried three times": `fetch_attempts: 3` (`walker/src/config.rs:476`) is three attempts in total (`for attempt in 1..=self.attempts`, `walker/src/stages/stream.rs:124`), i.e. two retries.
- `walker/README.md:19` "Two runtime crates, `ureq` and `zstd`": `walker/Cargo.toml:19-24` lists flate2, libc, signal-hook, ureq and zstd.
- `docs/galaxy.md:43-44` says the walker names "its summary and multi-anchor tags after the reference file's stem": `anchor=` is the record id, and the summary is `<stem>__<record_id>` for multi-record files (`walker/src/io/input.rs:292-296`).

**Proposed fix.** "are tried up to three times (two retries)"; "Runtime crates: `flate2`, `libc`, `signal-hook`, `ureq`, `zstd`"; "because the walker names its summary after the reference file's stem (`<stem>__<record id>` for a multi-record file)".
