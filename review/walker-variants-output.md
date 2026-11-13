# Walker: majority correction, variants, VCF and output writers (variants.rs, io/output.rs)

The variant stage reads the graph correctly at the level of single forks: the window arithmetic in `find_alts`, the VCF's REF/ALT semantics and the banded aligner's common case all behave as designed and are covered by unit tests. The weak point is the decision logic around those forks. `prefer_majority`, on by default, compares a stitched window against a whole-path length and so discards every correction whenever the GAF row starts or ends inside the path, as the real GraphAligner rows in the repository do (WV-01). Its gates and `find_alts`'s gates divide the fork's coverage by exits that never rejoin the route, treat Logan's clipped `ka:f` values as linear depth, and skip walked minority branches between the prune floor and `--min-vaf`. The second theme is silence: alleles are dropped by a cap, a length guard or a multi-substitution rule, and a "treat as provisional" decision is printed only to stderr, while `summary.json` and the VCF record nothing. Most fixes are small. They share two directions: put both stages on one shared set of helpers (`view()`, the fork total, the depth floor), and give every skip or override a machine-readable trace in `summary.json` or a VCF FILTER. Several entries (WV-02, WV-07, WV-11, WV-12) only become reachable on real data once WV-01 is fixed, so they should ship together with it, followed by a panel re-run.

| ID | Severity | Where | Title | Fix effort |
|---|---|---|---|---|
| WV-01 | critical | `walker/src/stages/variants.rs:259` | Default-on majority correction is discarded unless the alignment spans the whole stitched path | small |
| WV-02 | high | `walker/src/stages/variants.rs:188` | A walked allele at 1-10% of a fork is neither corrected nor reported | small |
| WV-03 | high | `walker/src/stages/variants.rs:385` | Logan `ka:f` is clipped at 25,000; VAFs and gates use it unflagged | small |
| WV-04 | medium | `walker/src/stages/variants.rs:385` | Fork VAF counts exits that never rejoin: deflated AF, lost alleles, blocked majority swap | small |
| WV-05 | medium | `walker/src/stages/variants.rs:553` | A refused prune reopens the variant gate: sub-floor branches published as PASS alleles | small |
| WV-06 | medium | `walker/src/stages/variants.rs:375` | Hard-coded cap of 32 alternatives silently truncates alt FASTA, VCF and snvs.tsv | small |
| WV-07 | medium | `walker/src/stages/variants.rs:291` | ">12 corrections, treat as provisional" is stderr-only; summary says `provisional: false` | small |
| WV-08 | medium | `walker/src/io/output.rs:612` | VCF and snvs.tsv drop every allele carrying two or more substitutions | small |
| WV-09 | medium | `walker/src/io/output.rs:581` | Stale 20 kb guard silently drops all SNV/VCF output for longer genes | small |
| WV-10 | medium | `walker/src/io/output.rs:726` | Banded alignment accepts a clipped in-band optimum; SNVs get wrong POS/REF | trivial |
| WV-11 | low | `walker/src/stages/variants.rs:230` | `--allow-frameshift-alts` also lets the majority swap reframe `_primary` | trivial |
| WV-12 | low | `walker/src/stages/variants.rs:192` | `prefer_majority` tries only the most abundant alternative at a fork | small |
| WV-13 | low | `docs/vcf_example.md:28` | VCF docs call ALT/AF the minor allele; ALT can be the majority allele | trivial |
| WV-14 | low | `walker/src/io/output.rs:615` | VCF text not validated: raw id as contig, IUPAC REF, per-sample values in INFO | small |
| WV-15 | low | `walker/src/io/output.rs:223` | Non-finite floats written as bare `inf`/`NaN`, making summary.json invalid | small |
| WV-16 | low | `walker/src/io/output.rs:475` | `--emit-vcf` deletes or overwrites any `<out-dir>/<ACC>.vcf`, walker-written or not | small |
| WV-17 | low | `walker/src/main.rs:562` | Reruns into one `--out-dir` keep stale variant outputs beside fresh primaries | small |
| WV-18 | low | `walker/src/main.rs:558` | Outputs written in place with the summary last; a late write error loses the summary | small |
| WV-19 | low | `walker/src/io/output.rs:476` | An over-long accession makes `--emit-vcf` abort before the summary is written | trivial |
| WV-20 | low | `walker/src/io/output.rs:207` | summary.json identifies references by id and length only; omits emit flags and variant decisions | small |
| WV-21 | low | `walker/src/stages/variants.rs:41` | Code docs say alt haplotypes are substitutions only and opt-in; indels are on by default | trivial |
| WV-22 | low | `walker/src/stages/variants.rs:320` | `find_alts` ignores the shared `view()`; oriented adjacency is built three times | small |
| WV-23 | low | `walker/src/stages/variants.rs:538` | variants.rs re-implements extract's median depth and relative floor | trivial |
| WV-24 | low | `walker/src/io/output.rs:67` | output.rs repeats writer setup and JSON arms; pident formula duplicated; header doc incomplete | small |

### WV-01 · critical · Default-on majority correction is discarded unless the alignment spans the whole stitched path

`walker/src/stages/variants.rs:259` · correctness · confirmed (reproduced) · fix effort: small · sources: RC-01, RS-01, SC1-01

**Problem.** `prefer_majority` rejects a swap when `if restitched.seq.len() != spliced_len { continue; }` (`variants.rs:259`). `spliced_len` is `stitched_len(cfg, &steps, subgraph)` (`variants.rs:225`), the length of the whole stitched path, but `stitch` returns only the aligned window, `stitched[start..end.max(start)]` (`walker/src/stages/stitch.rs:196`). The two are equal only when the GAF row has `path_start == 0` and `path_end >= path_len`. The comment above the guard says it "cannot fire today"; it fires on both real GAF rows the repository records (`walker/src/stages/align.rs:426`: path length 1288, window 7..1281). Every unit test uses a path with `path_start: 0` and the golden gate runs `--no-prefer-majority` (`bench/check_golden.sh:59`), so nothing exercises it. The guard was added in 257401e and is in v0.1.0 and v0.1.1.

**Impact.** With default options, any accession whose alignment starts or ends inside the path (only a tip-to-tip alignment escapes) gets no majority correction. Where the aligner took the rarer, reference-matching branch, `_primary` keeps it, no `consensus:` line is printed, and the outputs are byte-identical to a `--no-prefer-majority` run; the alternative is published with `vaf` above 0.5. The project's own panel measured the correction at 59 of 311 two-allele sites across 26 of 200 accessions (`docs/issue_intra_accession_variance_lost.md:197-200`), and the Galaxy wrapper offers the option checked by default. Each base is still a real allele, and the majority allele survives in `alt_haplotypes.fasta` with a correct VAF.

**Evidence.** A synthetic accession (300 nt gene, two 60/40 SNP sites) with a stand-in aligner: with a GAF window 0..300 the primary carries the 60% strain and stderr prints `consensus: 2 site(s) ... corrected`; with the same data and a window 40..340 the primary carries the 40% strain, no consensus line, alts `vaf=0.600 primary_vaf=0.400`, identical to `--no-prefer-majority`.

**Proposed fix.** Compare windows, as `find_alts` already does (`variants.rs:481`):

```rust
let Some(first) = stitch::stitch(cfg, graph, path, subgraph, query)? else { return Ok(None) };
let mut window = first.seq;
// ... in the loop, replacing the length check:
if restitched.seq.len() as i64 != window.len() as i64 + delta { continue; }
if restitched.seq == window { continue; } // same published sequence: not a correction
window = restitched.seq;
```

Prefer this over re-deriving the clamp arithmetic by hand, which is how the bug arose. The fix activates swap code, including indel swaps, that has not run on real rows in any release; the `SRR19878648` cascade (27 swaps) becomes reachable again.

**Tests and docs.** Add `variants.rs` tests with the window inside the path (substitution, +3 indel, no-op swap) and a `tests/cli.rs` case with a flanked GAF row (e.g. `140 10 130`); add a default-flags run to the golden gate. Re-run the 200-accession panel and refresh `config.rs:182-196`, `variants.rs:286-290` and `docs/issue_intra_accession_variance_lost.md:197-209`. Release note: 0.1.0 and 0.1.1 primaries were not majority-corrected. Land WV-07, WV-11 and WV-12 with this fix.

### WV-02 · high · A walked allele at 1-10% of a fork is neither corrected nor reported

`walker/src/stages/variants.rs:188` · science · confirmed (reproduced) · fix effort: small · sources: SC1-03

**Problem.** Pruning (`--prune-frac 0.01`) removes walked branches below 1% of gene depth (`walker/src/stages/extract.rs:157-163`). `prefer_majority` corrects a site only if the walked branch also clears `--min-vaf` (0.10): `if taken_ab < floor || (taken_ab as f64 / total) < cfg.min_vaf { continue; // ... error-branch case }` (`variants.rs:188-190`), and `find_alts` uses the same gate (`variants.rs:393-395`). A walked branch between 1% and 10% of the fork is therefore not pruned, not swapped and not reported.

**Impact.** Trigger: a sample whose majority differs from the reference, with a reference-like minority at 1-10% of the fork (a lab-strain or control contaminant, a wild-type minority beside a resistance mutation). GraphAligner walks the minority because it matches the reference. With default options `_primary` carries the minority base at every such site with `pident=100.0`, and there is no alt record, VCF or snvs.tsv row, summary field, or `provisional` flag. In the real SRR8209082 unitigs the reference allele sits continuously between high-depth flanks at three positions (T 550x vs C 20,000x; G 310x vs A 25,000x; T 244x vs C 15,000x), all above the ~145x prune floor; the route itself was not observed because GraphAligner was not run.

**Evidence.** A synthetic 900 nt gene at 10,000x with two strains five SNVs apart, run with the shipped binary and a stand-in aligner (no flanks, so WV-01 does not interfere): with the reference-like strain at 2%, 5% or 9%, the primary equals the minority strain, with no `consensus:` line and zero alts, VCF rows or snvs.tsv rows. At 15% the run prints `consensus: 5 site(s) ... corrected`.

**Proposed fix.** Let `prefer_majority` move a prune-surviving walked branch onto the fork's outright majority, with the existing rejoin and length rules:

```rust
let taken_minor = taken_ab < floor || (taken_ab as f64 / total) < cfg.min_vaf;
let prune_survivor = cfg.prune_frac > 0.0
    && taken_ab >= crate::stages::extract::depth_floor(cfg, depth);
if taken_minor && !prune_survivor { continue; }
// in the alternative loop:
if taken_minor && (ab as f64 / total) < 0.5 { continue; }
```

A second verifier found this exemption leaves error alleles in place when the prune is refused (WV-05). Dropping the exemption, so any walked minority moves onto a strict majority, fixed every tested case but fails the deliberate test `an_error_branch_site_is_not_touched`; the maintainers must choose. Either way, add a per-accession `minority_sites` list to summary.json (forks where an untaken branch holds at least half the fork and the route was not moved, with a reason: below_prune_floor, frameshift, over_indel_cap, no_rejoin, prefer_majority_off), and count these swaps separately. Requires WV-01 to take effect on real alignments.

**Tests and docs.** Add `a_prune_surviving_minority_branch_is_corrected` (fails today), `without_pruning_a_minority_branch_is_left_to_the_old_rule` and `off_a_minority_branch_only_onto_an_outright_majority`. Remove the "bias fixed" claims at `README.md:121`, `docs/issue_intra_accession_variance_lost.md:3` and `variants.rs:16-28`; update `docs/walker.md:23-24` and the `prefer_majority` doc in `config.rs`. Re-run the panel.

### WV-03 · high · Logan `ka:f` is clipped at 25,000; VAFs and gates use it unflagged

`walker/src/stages/variants.rs:385` · science · confirmed (reproduced) · fix effort: small · sources: SC1-02, R2-01, SWD-01

**Problem.** Logan caps abundances (its documentation states a 50,000 cap on (k+1)-mer counts); the `ka:f` values actually published stop at exactly 25000.0. The walker treats `ka:f` as linear depth: `let total: f64 = branches.iter().map(|&b| abund(b) as f64).sum();` (`variants.rs:385`), then `taken_ab as f64 / total` and `alt_ab as f64 / total` (`variants.rs:390`, `:402`) are gated on `--min-vaf` and published as AF/vaf. When a branch is clipped, `total` is a lower bound and the minor allele's share is overestimated. The same values drive `prefer_majority` (`variants.rs:183-198`), bubble VAFs (`walker/src/stages/phase.rs:183`) and the prune floor (`walker/src/stages/extract.rs:127`, `:157`). The VCF FILTER is always `PASS` (`walker/src/io/output.rs:493`) and `docs/walker.md:333` describes `ka:f` only as depth.

**Impact.** Deep samples (viral, amplicon) with default options. A 3,000x minor allele beside a true 250,000x major allele (1.2%) is emitted as AF 0.107, PASS, and clears `--min-vaf 0.10`; the prune floor becomes 250x instead of 2,500x. The repository's showcase is in this regime: `docs/examples/zika_ns5/SRR6251967.zika_ns5.vcf` has `RO=25000.0` in 11 of 11 rows, and 2,461 of the 2,564 NS5 reference unitigs in that run sit at 25000.0. Affected outputs: `alt_haplotypes.fasta` and `bubbles.json` (default), the VCF and snvs.tsv. When both branches are clipped, the tie leaves the minority allele on the primary.

**Evidence.** In the full SRR8209082 unitig object, 802 of 81,534 unitigs carry `ka:f:25000.0` and none is higher (next values 24500 and 24000). A one-bubble fixture with major 25,000 and minor 3,000 gives `gene 121 . C T . PASS AF=0.1071;AO=3000.0;RO=25000.0;FORK=0`; the same fixture with an uncapped major of 250,000 gives no row.

**Proposed fix.** Flag, do not drop; the true depth cannot be recovered.

```rust
// types.rs
pub const KA_CEILING: f32 = 25_000.0;   // empirical; see parse-time warning below
pub struct AltHaplotype { /* ... */ pub ka_saturated: bool }
// variants.rs, find_alts
let ka_saturated = branches.iter().any(|&b| abund(b) >= KA_CEILING);
```

In `output.rs`, write FILTER `KaSat` instead of `PASS` with a `##FILTER` header line (OR the flag when rows coalesce), append a trailing `ka_saturated` column to snvs.tsv, append ` ka_sat=1` to the alt header only when set (goldens stay byte-identical), and add `"depth_saturated"` to the summary's prune block. Also flag saturated bubbles in `phase.rs`, report `ka_max`, and warn in `parse::build` if any `ka:f` exceeds the constant.

**Tests and docs.** Add `vcf_tests::a_saturated_fork_is_flagged_not_passed` and a `variants.rs` test where a branch at 25,000 sets the flag and one at 24,999 does not; update the prune JSON assertion in `identity_label_tests`. Docs: `docs/walker.md:333`, `docs/vcf_example.md`, `docs/issue_intra_accession_variance_lost.md:161-162`, the VCF `##source` line; regenerate `docs/examples/zika_ns5/`. Changelog: `bcftools view -f PASS` now drops these rows. See WV-04 for a separate error in the same `total`.

### WV-04 · medium · Fork VAF counts exits that never rejoin: deflated AF, lost alleles, blocked majority swap

`walker/src/stages/variants.rs:385` · correctness · confirmed (reproduced) · fix effort: small · sources: RC-02, R5-06, SC1-08

**Problem.** `find_alts` and `prefer_majority` sum every out-edge of the fork (`variants.rs:385`, `:183`), including exits that never return to the route. `phase.rs:153-158` documents this as a defect it already fixed and restricts its own total to in-bubble branches, so `bubbles.json` and the VCF/alt FASTA give the same site different VAFs. In addition, `prefer_majority` picks `best` (`variants.rs:192-205`) before calling `rejoin`; if the best exit does not rejoin, the whole site is skipped.

**Impact.** A fork with a well-covered third exit that leaves the locus (repeat, paralog, contaminant sharing the 30-mer): walked 80x, allele 20x, exit 400x gives VAF 0.04, so no alt, VCF or snvs.tsv row, while `bubbles.json` reports 0.2. Walked 20x, allele 80x, exit 400x leaves the 20% base on `_primary`. In the shipped Zika VCF, POS 246 has `AF=0.1093` with `AO=4400` and `RO=25000` (two-allele ratio 0.150).

**Evidence.** Fixtures for both cases with the shipped binary: zero alts and no VCF in the first, no `consensus:` line and the minority on the primary in the second; with the fix both match the no-third-exit controls (`AF=0.2000`; `consensus: 1 site(s)`).

**Proposed fix.**

```rust
fn site_total(out, at, branches, taken, fork_index, max_detour, abund) -> f64 {
    branches.iter()
        .filter(|&&b| b == taken || rejoin(out, at, b, fork_index, max_detour).is_some())
        .map(|&b| abund(b) as f64).sum()
}
```

Use it at both sites, and in `prefer_majority` call `rejoin` inside the candidate loop so `best` holds only a rejoining exit (the total fix alone leaves the second case wrong). Optionally cache rejoin results and report excluded coverage in an `OTHER` INFO field.

**Tests and docs.** Unit tests on a three-exit fork for both cases, plus a CLI test. Update the AF description at `walker/src/io/output.rs:688` and `docs/walker.md:212`; regenerate `docs/examples/zika_ns5` and the golden fixtures.

### WV-05 · medium · A refused prune reopens the variant gate: sub-floor branches published as PASS alleles

`walker/src/stages/variants.rs:553` · science · confirmed (reproduced) · fix effort: small · sources: SWR5-01

**Problem.** `abundance_floor` returns `cfg.min_abund` (5x) unless `--min-abund-frac` is set (`variants.rs:553-559`; default `min_abund_frac: 0.0`, `walker/src/config.rs:458`). The design relies on the prune pass to remove sub-floor unitigs first (`walker/src/pipeline.rs:858`). When the prune is refused (no haplotype, shorter span, would empty), `pipeline.rs:843-845` rebuilds the unpruned subgraph and `prefer_majority`/`find_alts` run on it with the 5x gate.

**Impact.** Gene at 1000x, fork with taken 30x and untaken 8x (below the 10x prune floor), pruned realignment fails: `alt_haplotypes.fasta` gets `alt1 vaf=0.211`, the VCF a `PASS` row and snvs.tsv a row. Only `summary.json`'s `prune.outcome` hints at it.

**Evidence.** Fixture where the stand-in aligner answers only the unpruned pass: stderr shows `floor 10.0x cost the alignment; keeping the unpruned result`, then `variants: ... floor 5.0x ... 1 recovered` and `gene 43 . A C . PASS AF=0.2105;AO=8.0;RO=30.0`.

**Proposed fix.** In `pipeline.rs`, before the variant block:

```rust
let vcfg_raised;
let vcfg: &Config = match prune.as_ref() {
    Some(r) if !matches!(r.result, PruneResult::Applied | PruneResult::NothingToRemove)
        && r.floor > cfg.min_abund =>
        { vcfg_raised = Config { min_abund: r.floor, ..cfg.clone() }; &vcfg_raised }
    _ => cfg,
};
// pass vcfg to prefer_majority and find_alts
```

Applied and pruned-off runs stay byte-identical. Optionally add FILTER `UnprunedGraph`.

**Tests and docs.** A CLI test with an aligner that fails the pruned pass, asserting no alt below the prune floor. Update the `pipeline.rs:858` comment, the `min_abund_frac` doc and `docs/walker.md:211-213`.

### WV-06 · medium · Hard-coded cap of 32 alternatives silently truncates alt FASTA, VCF and snvs.tsv

`walker/src/stages/variants.rs:375` · correctness · confirmed (reproduced) · fix effort: small · sources: R5-05, SC1-04

**Problem.** `find_alts` walks 5' to 3' and stops at `if alts.len() >= cfg.max_alts { break; }` (`variants.rs:375-377`, again `:398`), with `max_alts: 32` (`walker/src/config.rs:466`) and no CLI flag. There is no warning, and the stderr counters stop with the loop. `prefer_majority` stops silently after `max_alts*2` = 64 swaps (`variants.rs:170`). `docs/contigs_vs_unitigs.md:78` reports 12-32 on-path bubbles per Zika accession, at the cap.

**Impact.** Accessions with more than 32 credible sites (quasispecies, mixed infections, long genes) get output for the 5' part only; downstream sites look invariant. The only clue is `bubbles=40` on the DONE line.

**Evidence.** 40 separate SNP bubbles: stderr `32 on-path fork(s), 32 biallelic candidate(s) -> 32 recovered`, 32 alt records, last snvs.tsv row at 1921; sites 1981-2401 absent.

**Proposed fix.** Drop the outer `break`; check the cap after the VAF/floor gate and count `beyond_cap += 1; continue`. When `beyond_cap > 0`, print `WARNING: N credible alternative allele(s) beyond --max-alts ... omitted` and add `alts_beyond_cap` to the per-accession summary. Add `--max-alts <N>` and consider a higher default. Same counter and warning for the `prefer_majority` iteration cap.

**Tests and docs.** Unit test with more than `max_alts` diamonds asserting the beyond-cap count; summary and CLI tests. Update `docs/walker.md:106` and the `config.rs:254-257` comment.

### WV-07 · medium · ">12 corrections, treat as provisional" is stderr-only; summary says `provisional: false`

`walker/src/stages/variants.rs:291` · provenance · confirmed (reproduced) · fix effort: small · sources: RC-06, RA-08, R5-13

**Problem.** When `swaps > 12` (`variants.rs:291`), `prefer_majority` prints "treat this accession as provisional", then drops the count: it returns `Option<GafPath>`, and the only producer of `AccessionResult::provisional` is `let provisional = stages::extract::exceeds_guard(cfg, &nodes, ref_len);` (`walker/src/pipeline.rs:577`). summary.json writes `"provisional": false` and records no correction count; `pident_inherited: true` is set from one correction up. `types.rs:484-486` already moved the neighbourhood guard out of stderr for this reason.

**Impact.** A rerouted primary (the panel's SRR19878648 took 27 corrections and added four gaps to the alignment) enters the FASTA, codon MSA and tree with no machine-readable mark; sharded and Galaxy runs do not read stderr. Rare today because of WV-01; routine once WV-01 is fixed.

**Evidence.** 13 two-allele bubbles with a whole-path GAF row: stderr prints the warning, summary has `"provisional": false` and no count. 12 bubbles: no warning, as intended.

**Proposed fix.** Return `Result<Option<(GafPath, usize)>>` from `prefer_majority` with `pub const MANY_CORRECTIONS: usize = 12`; add `majority_corrections: usize` to `AccessionResult`; at `pipeline.rs:866` set `majority_corrections = swaps; provisional |= swaps > MANY_CORRECTIONS;`; write `"majority_corrections"` and preferably `"provisional_reasons": ["tangled_neighbourhood", "many_majority_corrections"]`; count each reason separately in `decision_lines` (`walker/src/main.rs:389`).

**Tests and docs.** `variants.rs` test returning `(path, 13)`; summary and `decision_lines` tests; a CLI test. Define both rules in `docs/walker.md:120-133`; update `docs/issue_intra_accession_variance_lost.md:207-209`.

### WV-08 · medium · VCF and snvs.tsv drop every allele carrying two or more substitutions

`walker/src/io/output.rs:612` · science · confirmed (reproduced) · fix effort: small · sources: D1-02, R5-04, SC1-04

**Problem.** `snv_sites` sets `multiple = true` on the second differing base (`output.rs:603-607`) and then discards the allele: `let (Some(i), false) = (diff, multiple) else { continue };` (`output.rs:612`). In a k=31 graph, two substitutions within 30 nt on one haplotype share one bubble, so `find_alts` returns one alternative with both, and neither site is reported. User docs say "Substitutions only; indel alleles stay in `alt_haplotypes.fasta`" (`docs/vcf_example.md:132`), implying all substitutions are reported.

**Impact.** Any minor haplotype with two SNVs less than 31 nt apart (including two changes in one codon) is absent from `<ACC>.vcf` and snvs.tsv with no message; `types.rs:328` counts 58 of 267 panel substitution alleles spanning more than 3 nt. Per-site counts and cohort matrices are biased against dense variable regions. Opt-in outputs only.

**Evidence.** One bubble with a single substitution yields a snvs.tsv row and a VCF; adding a second substitution 4 nt away yields `1 recovered`, an alt header with `span=5`, no VCF and a header-only snvs.tsv.

**Proposed fix.**

```rust
let diffs: Vec<usize> = primary.seq.iter().zip(alt_seq).enumerate()
    .filter(|(_, (x, y))| x != y).map(|(i, _)| i).collect();
if diffs.is_empty() || diffs.len() > MAX_LINKED_SNVS { skipped_block += 1; continue; }
for i in diffs { /* existing row construction, unchanged */ }
```

Rows already share FORK; add `##INFO=<ID=BLOCK,...>` and amend the `##source` sentence "two rows do not assert co-occurrence" for rows sharing a BLOCK. Report skipped alleles (indel or over the cap, e.g. 8) on stderr and in the summary.

**Tests and docs.** `two_linked_substitutions_yield_two_rows`, a cap test and a CLI case. Update `docs/vcf_example.md:132`, the `docs/walker.md` outputs table and the `snv_sites` doc comment.

### WV-09 · medium · Stale 20 kb guard silently drops all SNV/VCF output for longer genes

`walker/src/io/output.rs:581` · correctness · confirmed (reproduced) · fix effort: small · sources: RA-07, R5-08

**Problem.** `if primary.seq.len() > 20_000 || query.seq.len() > 20_000 { eprintln!(...); return Vec::new(); }` (`output.rs:581-589`) dates from the full O(n·m) DP. Since 91ab897 the banded `align_primary_to_query` costs `(n+1)·w` bytes at high identity, so the guard rejects cheap inputs while doing nothing about divergent short pairs. Its last banded passes can also allocate about twice the full matrix.

**Impact.** `--emit-vcf`/`--emit-snv-table` against a reference over 20 kb (SARS-CoV-2 ORF1ab, 21,290 nt; whole coronavirus genomes): no `<ACC>.vcf`, header-only snvs.tsv, nothing in summary.json, exit 0. `docs/walker.md:109` says a VCF is written only when there are sites, so its absence reads as "no variants". Shipped panels are all under 20 kb.

**Evidence.** A 60/40 site at mid-gene: 19,998 nt gives a VCF row; 21,291 nt gives `SNP sites skipped: sequence longer than 20 kb`, no VCF, exit 0.

**Proposed fix.** Bound by memory, not length:

```rust
const SNV_ALIGN_BUDGET: usize = 256 << 20;
fn align_within(p: &[u8], q: &[u8], budget: usize) -> Option<Vec<Option<usize>>> {
    let (n, m) = (p.len(), q.len());
    let mut band = (n.abs_diff(m) + 16).max(32);
    while band < n.max(m) {
        let w = n.abs_diff(m) + 2 * band + 1;
        if (n + 1).saturating_mul(w) > budget { return None; }
        if let Some(map) = align_banded(p, q, band) { return Some(map); }
        band *= 2;
    }
    if (n + 1).saturating_mul(m + 1) > budget { return None; }
    Some(align_full(p, q))
}
```

Delete the length test; on `None`, warn "too divergent to map within 256 MB" and record `snv_skipped` in summary.json. Coordinate with WV-10, which changes `align_banded`.

**Tests and docs.** `long_gene_at_high_identity_is_not_skipped` (21,291 nt, one SNV) and `a_pass_over_the_budget_is_refused_not_allocated`. Update the `output.rs:709` comment, `docs/walker.md:109-110` and `docs/vcf_example.md`.

### WV-10 · medium · Banded alignment accepts a clipped in-band optimum; SNVs get wrong POS/REF

`walker/src/io/output.rs:726` · correctness · confirmed (reproduced) · fix effort: trivial · sources: R5-02

**Problem.** `align_primary_to_query` returns the first banded result whose traceback avoids the band edge (`output.rs:726-727`); the only rejection is `if d_here <= d_min || d_here >= d_max { return None; }` (`output.rs:804`). That is not an optimality certificate: when the optimum leaves the band, the best in-band path can be strictly worse and still stay off the edge. The doc comment's "never a wrong answer" (`output.rs:714-719`) is false.

**Impact.** A primary with an insertion and a deletion of opposite sign, each longer than the initial band (32 nt for equal lengths), e.g. +33 and -33 nt 900 nt apart in a 2.4 kb gene: about 940 bases map to wrong reference positions, and an SNV at true POS 761 is written at 783-828 with REF from the wrong base, silently. Rare.

**Evidence.** The three functions (`output.rs:720-882`) compiled standalone: for ±30 nt the banded and full maps agree; for ±33, ±40, ±60 nt the banded score is 1367/1336/1289 against 2233/2198/2098 full, and the SNV lands at 809/783/828 instead of 761. The repository's equivalence test uses at most three indels of at most 9 nt (`output.rs:914-916`).

**Proposed fix.** In `align_banded`, before traceback, accept only with a score certificate (any path reaching diagonal `d_max+1` or `d_min-1` needs at least `g` gaps):

```rust
let best = prev[col(n, m)] as i64;
let g_hi = (d_max + 1) + (d_max + 1 - delta);
let g_lo = (1 - d_min) + (delta - d_min + 1);
let g = g_hi.min(g_lo) as i64;
if best < (n as i64 + m as i64 - g) / 2 - 2 * g { return None; } // widen the band
```

Keep the edge check. With this, all cases match `align_full`, 400 random pairs show no difference, and a 2.7 kb 98%-identity pair is still accepted at band 32.

**Tests and docs.** Add compensating ±40 nt indels to `banded_equivalence_tests` and an `snv_sites` POS test; rewrite the comment at `output.rs:714-719`.

### WV-11 · low · `--allow-frameshift-alts` also lets the majority swap reframe `_primary`

`walker/src/stages/variants.rs:230` · correctness · confirmed (reproduced) · fix effort: trivial · sources: SWR5-02

**Problem.** `prefer_majority` reuses the side-file length rule `|| (delta % 3 != 0 && !cfg.allow_frameshift_alts))` (`variants.rs:230`), so the flag documented as admitting frameshifting alleles into `alt_haplotypes.fasta` also lets a frameshifting majority branch into `_primary`, against the stated intent at `variants.rs:254-258` and in the test `a_frameshifting_majority_is_not_swapped_in`. On the panel, 190 of 203 length-changing calls are ±1 nt homopolymer artifacts, so users of the flag (offered in the Galaxy wrapper) would get out-of-frame primaries in the codon MSA. Latent until WV-01 is fixed; medium after that.

**Evidence.** Taken 40x in frame, +1 nt branch 60x, full-span GAF: with the flag the primary grows from 63 to 64 nt and stderr prints `consensus: 1 site(s) ... corrected`.

**Proposed fix.** Replace the clause with `|| delta % 3 != 0` and reword the comment at `variants.rs:220-224` ("except that the consensus never takes a frameshift").

**Tests and docs.** Extend the existing test with `Config { allow_frameshift_alts: true, ..maj(true) }` expecting `None`. Update `main.rs:104-105`, `docs/walker.md:216` and the Galaxy help: "side file only; the primary is never reframed".

### WV-12 · low · `prefer_majority` tries only the most abundant alternative at a fork

`walker/src/stages/variants.rs:192` · correctness · confirmed (reproduced) · fix effort: small · sources: SWR5-03

**Problem.** `prefer_majority` keeps one `best` candidate (`variants.rs:192`); if it is then refused (no rejoin, length or frameshift rule, window, coverage), the loop `continue`s to the next fork, and a second credible branch that also beats the walked one is never tried. Triallelic fork, walked 25x, +1 nt artifact 45x, substitution 30x: the primary keeps the 25% allele and the alt reports `vaf=0.300 primary_vaf=0.250`, identical to `--no-prefer-majority`. Uncommon; latent until WV-01 is fixed.

**Proposed fix.** Try each better alternative in descending order:

```rust
let mut better: Vec<(Oriented, f32)> = branches.iter().copied()
    .filter(|&a| a != taken).map(|a| (a, abund(a)))
    .filter(|&(_, ab)| ab >= floor && (ab as f64 / total) >= cfg.min_vaf && ab > taken_ab)
    .collect();
better.sort_by(|a, b| b.1.total_cmp(&a.1)); // stable: ties keep branch order
for (alt, _) in better { /* existing body; its `continue`s try the next alt */ }
```

**Tests and docs.** Add `a_refused_top_alternative_falls_through_to_the_next`; update the `prefer_majority` doc. Interacts with WV-04's rejoin change.

### WV-13 · low · VCF docs call ALT/AF the minor allele; ALT can be the majority allele

`docs/vcf_example.md:28` · docs · confirmed (reproduced) · fix effort: trivial · sources: R5-14

**Problem.** The doc says "`ALT` the minor allele, `AF` its `ka:f` frequency" (`docs/vcf_example.md:28`) and "Each cell is a minor-allele frequency" (`:111`). The code lists every sample allele that differs from the reference base (`walker/src/io/output.rs:625-630`), so where the sample's majority differs from the reference, ALT is the majority allele with AF above 0.5, and where neither matches there are two ALTs with `RO=0.0`. The VCF is correct; users following the doc mix frequencies above 0.5 into minor-allele analyses. The committed Zika example has no such row.

**Proposed fix.** Line 28: "`ALT` is each sample allele that differs from the gene base (usually the minor allele; the majority allele when the sample's consensus differs from the reference; both when neither matches), `AF` that allele's frequency, `RO` the gene base's coverage (0 when absent)." Line 111: frequency of the non-reference allele; grey means not polymorphic in that sample. Optionally show `minor_af = min(af, 1-af)` in the SQL. Also edit `:88`, `:114`.

**Tests and docs.** Pin the major-allele-ALT case with an `snv_sites` unit test.

### WV-14 · low · VCF text not validated: raw id as contig, IUPAC REF, per-sample values in INFO

`walker/src/io/output.rs:615` · data · confirmed (reproduced) · fix effort: small · sources: R5-11

**Problem.** `##contig=<ID={chrom},...>` (`output.rs:687`) and CHROM take the raw first token of the reference header, or the file stem, which may contain spaces (`output.rs:562-568`); `>gene,with<odd>=chars` produces a malformed header line. REF is copied verbatim, `let ref_base = query.seq[qpos];` (`output.rs:615`), so an IUPAC code gives `REF=R`, which VCF 4.2 does not allow. AF/AO/RO/FORK are written to INFO as well as FORMAT (`output.rs:493`), so after `bcftools merge` INFO holds one sample's value (`zika_ns5.cohort.merged.vcf:19`: `AF=0.1327` against samples 0.1564, 0.132, 0.1327). No shipped reference has such characters; how htslib reacts was not observed.

**Proposed fix.** Reject reference ids outside the VCF contig grammar `[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*` with a usage error after loading queries; write any non-ACGT reference base as `N`; move AF/AO/RO/FORK to FORMAT only, or document that merged INFO is unreliable.

**Tests and docs.** Unit test for the `N` mapping; CLI test for an id with `,`; document the id rule in `docs/walker.md` and INFO versus FORMAT in `docs/vcf_example.md`.

### WV-15 · low · Non-finite floats written as bare `inf`/`NaN`, making summary.json invalid

`walker/src/io/output.rs:223` · correctness · confirmed (reproduced) · fix effort: small · sources: R5-10, RB-07, SWR2-01

**Problem.** Settings are written as `format!("{:?}", c.min_cov)` (`output.rs:223`, also `:227-229`) and `"mean_ka": {:.1}` (`:271`), producing `inf`/`NaN`, which no JSON parser accepts. `--min-abund`, `--min-vaf` and `--min-cov` are a bare `parse()` with no range check (`walker/src/main.rs:191`, `:194`, `:197`), and `ka:f:` accepts `inf`/`nan` (`walker/src/stages/parse.rs:147`). `--min-abund inf` exits 0 with a summary `json.load` rejects; `--min-cov NaN` also disables the coverage gate; `--min-vaf 5` is accepted.

**Proposed fix.** A `finite_in(flag, raw, lo, hi)` helper in `main.rs` (`--min-cov`, `--min-vaf` in [0, 1]; `--min-abund` >= 0); `.filter(|v| v.is_finite() && *v >= 0.0)` in `parse.rs`; a `json_f64` writer emitting `null` for non-finite values.

**Tests and docs.** `parse_args` and `parse.rs` rejection tests; a CLI test parsing the whole summary strictly; state ranges in `docs/walker.md:202, 211-212`.

### WV-16 · low · `--emit-vcf` deletes or overwrites any `<out-dir>/<ACC>.vcf`, walker-written or not

`walker/src/io/output.rs:475` · data · confirmed (reproduced) · fix effort: small · sources: RE-05, R5-09

**Problem.** For an accession without sites `write_vcfs` calls `std::fs::remove_file(&path)` (`output.rs:475`); otherwise `File::create` truncates it (`:484`). Neither checks that the file starts with the walker's `##source=logan-walker` header (`:683`). Pointing `-o` at a directory holding per-sample VCFs from a read-based caller with the same `<ACC>.vcf` names deletes or replaces them without a message. Observed: a `##source=bcftools` `SRRFAKE0.vcf` vanished on a run with no sites, exit 0.

**Proposed fix.** A `walker_owned(path)` check (missing, or second line starts with `##source=logan-walker`) before both calls; otherwise error with "refusing to replace <path>: not written by logan-walker; use another --out-dir". Keep the stale-VCF removal for walker-owned files.

**Tests and docs.** Extend `a_rerun_without_sites_removes_the_old_vcf` with a foreign VCF; update `docs/walker.md:109`. See WV-17 and WV-19 for the same code path.

### WV-17 · low · Reruns into one `--out-dir` keep stale variant outputs beside fresh primaries

`walker/src/main.rs:562` · reliability · confirmed (reproduced) · fix effort: small · sources: R5-09

**Problem.** Each optional writer runs only when its flag is set (`main.rs:562-577`), so a rerun without `--emit-vcf`, `--emit-snv-table` or with `--no-emit-alts` leaves the previous run's `<ACC>.vcf`, `snvs.tsv` and `alt_haplotypes.fasta` next to a fresh `r.fasta` and summary; summary.json does not record the emit flags (WV-20), so nothing shows the mismatch. The variant file names carry no gene, so two genes run into one out-dir overwrite each other's `snvs.tsv`, `bubbles.json` and alts while both summaries survive. Observed with a rerun two seconds later: the three variant files kept the first run's timestamps. Galaxy runs in fresh directories and is unaffected.

**Proposed fix.** Define `VCF_SIGNATURE = b"##fileformat=VCFv4.2\n##source=logan-walker"` and `remove_if_ours(path, sig)`; when a flag is off, remove this run's walker-written `<ACC>.vcf` files and `snvs.tsv` (signature = its header) and log each removal. For `alt_haplotypes.fasta` (no signature), always write it or prefix the variant files with the query stem like the summary; renaming is a breaking change for `docs/vcf_example.md`, `pipe_*.py` and the Galaxy wrappers.

**Tests and docs.** CLI test: run with both flags, rerun without, assert the files are gone. Update `docs/walker.md:109`.

### WV-18 · low · Outputs written in place with the summary last; a late write error loses the summary

`walker/src/main.rs:558` · reliability · confirmed (reproduced) · fix effort: small · sources: R5-12

**Problem.** Every writer opens the final path, `BufWriter::new(std::fs::File::create(path)?)` (`walker/src/io/output.rs:72` and five more), without temp-file-and-rename or `sync_all`, and `main.rs:558-581` writes the FASTA first and the summary last, each step returning on `?`. A write error after compute (full disk, a directory named `<ACC>.vcf`) exits 1 with `r.fasta` and `bubbles.json` present and no `alt_haplotypes.fasta` or summary, losing the shard's per-accession status and timings; a signal mid-write leaves a truncated file under the final name. The exit status is non-zero, so the failure is detectable. The cache writer already writes atomically (`walker/src/stages/stream.rs:249`).

**Proposed fix.** An `AtomicFile` writing `<dst>.tmp.<pid>` whose `commit()` flushes, syncs and renames and whose `Drop` removes an uncommitted temp file; in `main.rs`, attempt every writer, always write the summary, then return the first error.

**Tests and docs.** CLI test with a directory named `<ACC>.vcf`: exit 1, summary and snvs.tsv present, no `*.tmp.*`. Unit test for `AtomicFile` drop. Note in `docs/walker.md` that the summary is always written.

### WV-19 · low · An over-long accession makes `--emit-vcf` abort before the summary is written

`walker/src/io/output.rs:476` · reliability · confirmed (reproduced) · fix effort: trivial · sources: RX-06

**Problem.** `validate_accession` (`walker/src/io/input.rs:96-106`) has no length limit, and `write_vcfs` treats every removal error except NotFound as fatal: `Err(e) if e.kind() != std::io::ErrorKind::NotFound => return Err(e.into()),` (`output.rs:476`). An accession line of 252 or more characters (accessions joined without a separator) makes `<acc>.vcf` exceed NAME_MAX; the run exits 1 with `File name too long (os error 36)` and writes no snvs.tsv, alt FASTA or summary, although the good accessions reconstructed. Without `--emit-vcf` the same list exits 0. Any other removal error (EACCES, EISDIR) has the same effect.

**Proposed fix.**

```diff
 // io/input.rs validate_accession
+    if s.len() > 64 {
+        return Err(Error::Usage(format!(
+            "{source}: accession of {} characters is not usable (at most 64)", s.len())));
+    }
 // io/output.rs write_vcfs
-  Err(e) if e.kind() != NotFound => return Err(e.into()),
+  Err(e) if e.kind() != NotFound => eprintln!("[warn] cannot remove stale {}: {e}", path.display()),
```

**Tests and docs.** Unit test for 64 and 65 characters; document the limit in `docs/walker.md`. WV-18's write-the-summary-anyway change covers the general case.

### WV-20 · low · summary.json identifies references by id and length only; omits emit flags and variant decisions

`walker/src/io/output.rs:207` · provenance · confirmed (by reading) · fix effort: small · sources: R5-13

**Problem.** `write_run_info` writes each reference as `{"id": ..., "length": ...}` (`output.rs:205-211`), so a changed reference with the same id and length leaves identical provenance. The settings array (`output.rs:216`) omits `emit_alts`, `emit_vcf` and `emit_snv_table`, and `write_result` records nothing the variant stage decided: alt count, the cap (WV-06), SNV site count or the 20 kb skip (WV-09), and correction count (WV-07). Nothing in the sequences is wrong; the audit trail is incomplete.

**Proposed fix.** Add `"sha256"` of the uppercased reference sequence (move `sha2` from dev-dependencies to dependencies in `walker/Cargo.toml`; it is already in Cargo.lock); add the three emit flags to settings; carry `alts`, `alts_capped` and `snv_sites` (`None` when skipped) on `AccessionResult` and write them, appending fields at the end of each object.

**Tests and docs.** Extend the summary unit test; update the `docs/walker.md` field table and any golden summaries.

### WV-21 · low · Code docs say alt haplotypes are substitutions only and opt-in; indels are on by default

`walker/src/stages/variants.rs:41` · docs · confirmed (reproduced) · fix effort: trivial · sources: RS-12, SWR5-04

**Problem.** The module doc has "# Why substitutions only" (`variants.rs:41`) and says each record carries exactly one substitution (`:34-35`); `:424` says "Substitutions only"; `:150-151` says a swap must "leave the sequence exactly as long"; `:8` says the output is empty unless `cfg.emit_alts`; `walker/src/stages/mod.rs:11` labels the stage "(opt-in)"; `walker/src/types.rs:300`, `:466-467` repeat "one substitution". The defaults are `allow_indel_alts: true`, `max_alt_indel: 30`, `emit_alts: true` (`walker/src/config.rs:460-463`), in-frame indels are accepted at `variants.rs:432-454`, and `:313` also runs for `emit_vcf`/`emit_snv_table`. A maintainer trusting the note could treat `path_start`/`path_end` as fixed after a splice, the bug `shifted_window` prevents. `cargo test accepts_an_in_frame_indel_allele` passes on the shipped tree. User-facing help and `docs/walker.md:215` are correct.

**Proposed fix.** Replace the section with "# Which alleles are accepted" (substitutions always; indels when `allow_indel_alts`, |delta| ≤ `max_alt_indel`, in-frame unless `allow_frameshift_alts`, window moved by `shifted_window`; keep the 1,279-to-334-node measurement as the reason for the cap); fix `:8`, `:34-35`, `:150-151`, `:424`, `mod.rs:11` ("on by default") and the `types.rs` lines; move the doc line at `:582` from `tail_contrib` to `stitched_len`.

**Tests and docs.** No new behaviour to pin: `accepts_an_in_frame_indel_allele` and `an_in_frame_indel_majority_is_swapped_in` already cover it. Comments only; user docs are correct.

### WV-22 · low · `find_alts` ignores the shared `view()`; oriented adjacency is built three times

`walker/src/stages/variants.rs:320` · reuse · confirmed (reproduced) · fix effort: small · sources: RS-03, SWR5-04

**Problem.** `View` is documented as "Shared by `prefer_majority` and `find_alts` so the two cannot drift apart" (`variants.rs:72`), but only `prefer_majority` calls `view()` (`:171`); `find_alts` rebuilds the same `local`, `out`, `walk`, `at` inline (`:320-355`), and `phase::phase` builds the adjacency a third time (`walker/src/stages/phase.rs:71-93`). `View.local` is never read (hidden by `#![allow(dead_code, unused_variables)]`, `walker/src/main.rs:24`). No wrong output today; a later adjacency change landing in one copy would make bubbles.json, the corrected primary and the alts disagree about which forks exist.

**Proposed fix.**

```rust
let Some(View { out, walk, at }) = view(graph, subgraph, path) else {
    return Ok(Vec::new()); // path leaves the subgraph; stitch would bail too
};
```

Remove `local` from `View`. Optionally move the builder and the `2*i+is_forward` decoding (`orient_of`, `strand_of`, phase's `unlocal`) into one module used by `phase`, `view` and the diagnostic test. Outputs are byte-identical; do this before any adjacency fix.

**Tests and docs.** Add a test asserting that `phase` reports a bubble at every fork where `find_alts` recovers an alt, so the readings cannot diverge silently. Reword or delete the comment at `variants.rs:326`.

### WV-23 · low · variants.rs re-implements extract's median depth and relative floor

`walker/src/stages/variants.rs:538` · reuse · confirmed (by reading) · fix effort: trivial · sources: RS-06

**Problem.** `variants::path_depth` (`variants.rs:538-548`) repeats `extract::path_depth`'s median (`walker/src/stages/extract.rs:127-133`), and `abundance_floor` (`variants.rs:553-559`) is `extract::depth_floor` (`extract.rs:157-163`) with a different fraction. Identical today; a change to one copy (a NaN-safe sort, length weighting) would make the prune floor and the allele floor measure gene depth differently.

**Proposed fix.**

```rust
// extract.rs
pub fn relative_floor(frac: f64, depth: f32, min_abund: f32) -> f32 {
    if frac > 0.0 { (frac as f32 * depth).max(min_abund) } else { min_abund }
}
pub fn depth_floor(cfg: &Config, depth: f32) -> f32 { relative_floor(cfg.prune_frac, depth, cfg.min_abund) }
// variants.rs: delete path_depth; call extract::path_depth(graph, &path.steps)
fn abundance_floor(cfg: &Config, depth: f32) -> f32 {
    crate::stages::extract::relative_floor(cfg.min_abund_frac, depth, cfg.min_abund)
}
```

**Tests and docs.** Rewrite `depth_is_the_median_of_the_walked_path` to call `extract::path_depth`. No user docs change.

### WV-24 · low · output.rs repeats writer setup and JSON arms; pident formula duplicated; header doc incomplete

`walker/src/io/output.rs:67` · reuse · confirmed (by reading) · fix effort: small · sources: RS-14

**Problem.** Five writers open with the same parent-dir-and-`BufWriter` block (`output.rs:67` and four more); `write_fasta` and `write_alts` build the same `anchor` suffix; `write_result` spells out seven identical `Option`-to-`null` arms (`output.rs:294-321`). The identity formula appears in `walker/src/stages/stitch.rs:213-217` and again in `walker/src/pipeline.rs:647-653`. The `write_fasta` doc (`output.rs:61`) and `docs/walker.md:105` omit the `ka=` and `vaf=` fields that line 84 writes, and the `write_alts` doc leaves a `text` fence unclosed so prose renders as code.

**Proposed fix.** Helpers `create(path)`, `anchor_suffix(h)`, `json_opt(w, key, v)` (keep `Display`, not `Debug`, so `100` stays `100`) and `GafPath::pident()` used at both call sites; correct the header doc to `>{id} cov={:.3} pident={:.1} nodes={n} ka={:.1} vaf={:.2}` plus optional ` anchor={rid}`. Outputs were byte-identical after the refactor. Close the fence in the `write_alts` doc (`output.rs:103-111`).

**Tests and docs.** Tests for `pident()` with `aln_len == 0` and 95/100, and for the full FASTA header including `ka=`, `vaf=` and `anchor=`; update `docs/walker.md:105`.
