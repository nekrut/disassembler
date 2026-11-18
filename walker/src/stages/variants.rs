//! Step 8: recover the intra-accession haplotypes the primary walk leaves behind.
//!
//! `stitch` emits one sequence per accession. Where the sample genuinely carries two alleles
//! at a site, the second one is discarded, and nothing downstream ever sees it. This stage
//! puts it back, as one extra record per site.
//!
//! IN   UnitigGraph, Subgraph, GafPath, Query
//! OUT  `Vec<AltHaplotype>`. Written to `alt_haplotypes.fasta` unless `--no-emit-alts`, and
//!      also computed for `--emit-vcf` and `--emit-snv-table`
//!
//! # What counts as a site worth recovering
//!
//! A fork the walked path passes through, where **both** the branch taken and the branch not
//! taken clear `min_abund` and `min_vaf`. That "both" is the whole filter for what
//! `find_alts` reports, and it separates real quasispecies from the far more common failure
//! next door. Measured on `SRR8209082` (HIV-1 *pol*), over the 58 on-path forks with a
//! credible untaken branch:
//!
//! ```text
//! walked branch is itself under 10%   54   the aligner took a rare allele; not a variant
//! both branches between 10% and 90%    4   genuinely two populations
//! ```
//!
//! The 54 are the reference-pull defect: GraphAligner scores paths by residue matches to the
//! reference, so where the sample's dominant allele differs from the reference and an error
//! branch matches it, the walk takes the error branch. Emitting those as variants would
//! publish sequencing error as biology, so they are not recovered here as alternatives. They
//! are counted in `minority_left` instead, and `--rescue-minority` is what acts on them.
//! `docs/notes/intra_accession_variance.md` covers the defect itself.
//!
//! # Why one variant per record
//!
//! With *n* independent sites there are up to 2^n haplotypes, and a de Bruijn graph cannot
//! say which alternatives sit on the same molecule. Emitting combinations would invent
//! sequences that may not exist in any organism, which is the chimera problem the walker
//! exists to avoid. Every record here is the primary path with **exactly one** allele changed,
//! so it is a real connected path through the graph and asserts only the minimal thing.
//!
//! In practice the cap is not the binding constraint: real site counts are single digits.
//!
//! # Which alleles are accepted
//!
//! Substitutions always. An indel allele when `allow_indel_alts` (**on by default**), its
//! `|delta|` is at most `max_alt_indel` (30), and it is in frame unless
//! `allow_frameshift_alts`. The alignment window moves with the splice via `shifted_window`,
//! so `path_start`/`path_end` are **not** fixed across a swap and code that assumes they are
//! is wrong.
//!
//! The length rules are not bookkeeping, they are the guard against the dangerous shape: a
//! detour that rejoins far downstream is a path *shortcut* that silently deletes gene.
//! Measured while prototyping, a coverage-greedy re-walk without any such constraint cut a
//! 1,279-node path to 334 and the sequence from 2,738 to 653 nt. `prefer_majority` applies
//! `max_alt_indel` to the whole correction set for the same reason: one panel accession took
//! 64 individually legal swaps compounding into +1,080 nt.

use crate::config::Config;
use crate::error::Result;
use crate::stages::{parse, stitch};
use crate::types::{
    AltHaplotype, GafPath, MajorityReport, MinorityLeft, Orient, PathStep, Query, Strand,
    Subgraph, UnitigGraph, UnitigId, VariantStats, KA_CEILING,
};
use std::collections::HashMap;

/// Oriented node id, packed as `2 * local_index + is_forward`, matching `phase`.
type Oriented = u32;

fn orient_of(o: Oriented) -> Orient {
    if o % 2 == 1 {
        Orient::Forward
    } else {
        Orient::Reverse
    }
}

/// The subgraph as an oriented adjacency, plus where the route sits inside it.
///
/// Shared by `prefer_majority` and `find_alts` so the two cannot drift apart in how they
/// read the graph. Oriented ids are packed as `2 * local_index + is_forward`, matching
/// `phase`: an edge in u's header is the traversal `(u, s1) -> (t, s2)`, and the reverse
/// arrives independently from t's own header.
struct View {
    out: Vec<Vec<Oriented>>,
    walk: Vec<Oriented>,
    at: HashMap<Oriented, usize>,
}

fn view(graph: &UnitigGraph, subgraph: &Subgraph, path: &GafPath) -> Option<View> {
    let n = subgraph.nodes.len();
    if n == 0 || path.steps.len() < 2 {
        return None;
    }
    let mut local: HashMap<UnitigId, usize> = HashMap::with_capacity(n * 2);
    for (i, &id) in subgraph.nodes.iter().enumerate() {
        local.insert(id, i);
    }
    let mut out: Vec<Vec<Oriented>> = vec![Vec::new(); n * 2];
    for (i, &id) in subgraph.nodes.iter().enumerate() {
        for &packed in graph.edges(id) {
            let (target, src_fwd, dst_fwd) = parse::unpack_edge(packed);
            let Some(&j) = local.get(&target) else { continue };
            out[i * 2 + usize::from(src_fwd)].push((j * 2 + usize::from(dst_fwd)) as Oriented);
        }
    }
    for v in out.iter_mut() {
        v.sort_unstable();
        v.dedup();
    }
    let mut walk: Vec<Oriented> = Vec::with_capacity(path.steps.len());
    for step in &path.steps {
        let &i = local.get(&step.node)?; // route left the subgraph; stitch would bail too
        let fwd = matches!(step.orient, Orient::Forward);
        walk.push((i * 2 + usize::from(fwd)) as Oriented);
    }
    let mut at: HashMap<Oriented, usize> = HashMap::with_capacity(walk.len());
    for (i, &o) in walk.iter().enumerate() {
        // First occurrence wins; a repeated unitig makes "where does it rejoin" ambiguous,
        // and the rejoin search only ever moves forward.
        at.entry(o).or_insert(i);
    }
    Some(View { out, walk, at })
}

/// Rewrite the route so that at every two-allele site it carries the **more abundant**
/// allele.
///
/// # Why this is needed on top of pruning
///
/// `--prune-frac` removes branches that are negligible against the gene's depth, which stops
/// the aligner walking onto sequencing error. It cannot help where both alleles are real:
/// there is no coverage floor that removes a 35% allele without also removing every genuine
/// minority variant. At those sites GraphAligner still prefers whichever branch resembles
/// the reference, and roughly half the time that is the rarer one.
///
/// Measured over the first 68 two-allele sites of a panel run, after pruning:
///
/// ```text
/// alternative is the more abundant allele   32 of 68  (47%)
///   of those, 0.50-0.55  near tie            3        harmless either way
///              0.55-0.65                    11
///              0.65-0.80                    13        primary carries the minority
///              0.80+                         5
/// ```
///
/// The variance is not lost in those cases, since both alleles reach the output. What is
/// wrong is the labelling: `_primary` is the record that feeds the codon alignment and the
/// selection analysis, and at those positions it carries the rarer base.
///
/// # Why this is not the repair that failed
///
/// An earlier attempt rewrote the route at *every* fork toward better coverage. It was
/// unstable: unguarded it cut a 1,279-unitig route to 334, and guarded it drifted, 657 swaps
/// for a net loss of identity. The difference here is the precondition. A swap is considered
/// only where the alternative clearly beats what was walked, which is a handful of sites per
/// accession rather than every fork, and each one must still rejoin the route and obey the
/// length rules (equal length, or an in-frame indel inside `max_alt_indel`, counted against
/// the whole correction set).
///
/// Returns the corrected route, or `None` if nothing needed changing.
pub fn prefer_majority(
    cfg: &Config,
    graph: &UnitigGraph,
    subgraph: &Subgraph,
    path: &GafPath,
    query: &Query,
) -> Result<(Option<GafPath>, MajorityReport)> {
    let mut report = MajorityReport::default();
    if !cfg.prefer_majority {
        // Nothing is examined, so nothing is reported. `settings.prefer_majority` in the
        // summary already says the stage was off; an empty report must not read as "looked
        // and found none".
        return Ok((None, report));
    }
    let mut baseline_len = stitched_len(cfg, &path.steps, subgraph);
    let mut current = GafPath { steps: path.steps.clone(), ..copy_coords(path) };
    let mut swaps = 0usize;

    // The sequence the aligner's window actually covers, which is what gets published.
    // `stitch` returns `stitched[path_start..path_end]` (`stitch.rs`), not the whole path, so
    // this is the sequence each swap must preserve the length of.
    //
    // The guard below must compare against this and not against the whole stitched path:
    // the two are equal only for a tip-to-tip alignment, and every real GAF row has a flank,
    // so comparing against the path silently rejects every correction.
    let Some(first) = stitch::stitch(cfg, graph, path, subgraph, query)? else {
        // The route does not clear the coverage gate; nothing to correct.
        return Ok((None, report));
    };
    let mut window = first.seq;
    // Total length the corrections have moved the primary by, across every swap so far.
    //
    // `--max-alt-indel` is a budget for the whole correction set, not a per-swap cap,
    // because swaps compound: one panel accession took 64 corrections, each in frame and
    // inside the 30 nt cap, and its ebola GP primary went from 2,049 nt to 3,129 against a
    // 2,031 nt reference. That is the greedy whole-route rewrite this stage exists to avoid,
    // reached one legal step at a time.
    let mut cumulative_delta: i64 = 0;

    // Each swap fixes one site and makes that site's taken branch the majority, so it cannot
    // be chosen again. The cap is a backstop against a shape that defeats that reasoning.
    let passes = cfg.max_alts.max(1) * 2;
    let mut hit_cap = true;
    for _ in 0..passes {
        let Some(v) = view(graph, subgraph, &current) else { break };
        let abund = |o: Oriented| graph.abundance_of(subgraph.nodes[(o as usize) / 2]);
        let depth = path_depth(graph, subgraph, &v.walk);
        let floor = abundance_floor(cfg, depth);

        let mut applied = false;
        // Reset every pass. The pass that applies nothing is the one that visited every fork
        // on the final route, so its tally is the set of sites actually left behind.
        let mut left = MinorityLeft::default();
        for i in 0..v.walk.len() - 1 {
            let (fork, taken) = (v.walk[i], v.walk[i + 1]);
            let branches = &v.out[fork as usize];
            if branches.len() < 2 || !branches.contains(&taken) {
                continue;
            }
            // Exits that never return to the route are not alleles of this gene and
            // must not sit in the denominator.
            let total = site_total(&v.out, &v.at, branches, taken, i, cfg.alt_max_detour, abund);
            if total <= 0.0 {
                continue;
            }
            let taken_ab = abund(taken);
            // A fork whose walked branch is itself below `min_abund` or `min_vaf` falls
            // between the two mechanisms: pruning removes branches below a fraction of GENE
            // DEPTH, this gate is a fraction of THE FORK. That is where a reference-like
            // contaminant or a wild-type minority beside a resistance mutation sits, and
            // GraphAligner walks it precisely because it matches the reference.
            //
            // Such a site is always reported. Acting on it is `--rescue-minority`, and then
            // only onto a branch holding an outright majority of the fork: "more than half of
            // everything here" is a far stronger claim than "bigger than what we took", and
            // the weaker one becomes the greedy whole-route rewrite.
            let taken_minor = taken_ab < floor || (taken_ab as f64 / total) < cfg.min_vaf;
            // A fork the aligner walked on the rarer side. If the route is still on it when
            // the loop ends, `_primary` knowingly carries the minority allele, and the report
            // says why.
            let untaken_max = branches
                .iter()
                .copied()
                .filter(|&b| b != taken)
                .map(abund)
                .fold(0.0f32, f32::max);
            let mut why: Option<fn(&mut MinorityLeft)> = None;
            // The rescue's bar is a ratio against `total`, and where a branch sits at
            // Logan's ceiling `total` is a LOWER BOUND, so the ratio is overstated and the
            // bar is not the bar. Only the rescue path is refused for this: an ordinary
            // two-allele correction needs the alternative to BEAT the walked branch, which
            // clipping cannot manufacture, since two clipped branches compare equal.
            if taken_minor {
                // Off by default: on a 5,000-accession slice, acting here took corrections
                // from 470 to 97,068 and the median corrected accession hit the loop's own
                // backstop, with 64 swaps producing 6 bases of net change, so the swaps were
                // undoing one another. The site is still reported; only the acting is off.
                if !cfg.rescue_minority {
                    left.exempt_minority += 1;
                    continue;
                }
                // Even with rescues on, refuse where the fork is clipped: the rescue's only
                // rail is a ratio against a total that saturation makes a lower bound.
                if branches.iter().any(|&b| abund(b) >= KA_CEILING) {
                    left.saturated += 1;
                    continue;
                }
            }
            // Every credible alternative that beats what we took, best first.
            //
            // All of them, not just the most abundant: when the best is refused further down
            // (no rejoin, a length or frameshift rule, the shifted window, the coverage gate)
            // a second credible branch that also beats the walked one should still be tried.
            // The `continue`s below therefore mean "try the next alternative", not "give up
            // on this fork".
            let mut better: Vec<(Oriented, f32)> = branches
                .iter()
                .copied()
                .filter(|&a| a != taken)
                .map(|a| (a, abund(a)))
                .filter(|&(_, ab)| {
                    // Always a credible allele that beats what we took. The VAF bar is the
                    // ordinary `min_vaf`, except when the walked branch is itself a minority
                    // where only an outright majority of the fork will do.
                    ab >= floor
                        && ab > taken_ab
                        && (ab as f64 / total) >= if taken_minor { 0.5 } else { cfg.min_vaf }
                })
                .collect();
            if better.is_empty() {
                if taken_minor {
                    if untaken_max < floor {
                        // The other side is below `min_abund` too: the whole fork is noise.
                        left.alt_below_floor += 1;
                    } else {
                        // Credible, but no single branch holds half the fork. A three-way
                        // split with the route on the smallest share lands here, and moving
                        // on a plurality is exactly the greedy rewrite this must not become.
                        left.not_outright_majority += 1;
                    }
                }
                continue;
            }
            // Descending abundance. `sort_by` is stable, so ties keep branch order and the
            // choice stays deterministic.
            better.sort_by(|a, b| b.1.total_cmp(&a.1));

            for (alt, _) in better {
                let Some((rejoin_at, detour)) = rejoin(&v.out, &v.at, alt, i, cfg.alt_max_detour)
                else {
                    why = Some(|l| l.no_rejoin += 1);
                    continue;
                };
                let mut steps: Vec<PathStep> = Vec::with_capacity(current.steps.len() + detour.len());
                steps.extend_from_slice(&current.steps[..=i]);
                for &o in &detour {
                    steps.push(PathStep {
                        node: subgraph.nodes[(o as usize) / 2],
                        orient: orient_of(o),
                    });
                }
                steps.extend_from_slice(&current.steps[rejoin_at..]);

                // Almost the same length rules as `find_alts`, so the two stages agree on what
                // counts as an allele. Without that agreement an indel site could not be
                // corrected and its alternative would be emitted with a VAF above 0.5, meaning
                // the primary had kept the rarer base there. Measured before the rules were
                // shared: 1 of 48 alternatives on `ebola_gp`, every one of them an indel.
                //
                // The one deliberate difference is the frameshift rule, below.
                let spliced_len = stitched_len(cfg, &steps, subgraph);
                let delta = spliced_len as i64 - baseline_len as i64;
                if delta != 0
                    && (!cfg.allow_indel_alts
                        || delta.unsigned_abs() as usize > cfg.max_alt_indel
                        // The consensus never takes a frameshift, whatever the flags say.
                        // `--allow-frameshift-alts` admits such alleles to the side file;
                        // honouring it here too would reframe `_primary`, which is what the
                        // codon alignment reads. Deliberately not `&& !allow_frameshift_alts`.
                        || delta % 3 != 0)
                    {
                        why = Some(|l| l.length_refused += 1);
                        continue;
                    }
                // The same cap, against the running total rather than this one swap. See
                // `cumulative_delta` above for the accession that made this necessary.
                if (cumulative_delta + delta).unsigned_abs() as usize > cfg.max_alt_indel {
                    why = Some(|l| l.length_refused += 1);
                    continue;
                }
                let prefix_len = stitched_len(cfg, &current.steps[..=i], subgraph);
                let replaced_len = tail_contrib(cfg, &current.steps[i + 1..rejoin_at], subgraph);
                let Some((new_start, new_end)) = shifted_window(
                    current.path_start,
                    current.path_end,
                    prefix_len,
                    replaced_len,
                    delta,
                ) else {
                    why = Some(|l| l.window_not_shiftable += 1);
                    continue;
                };
                let candidate = GafPath {
                    steps,
                    path_start: new_start,
                    path_end: new_end,
                    ..copy_coords(path)
                };
                let Some(restitched) = stitch::stitch(cfg, graph, &candidate, subgraph, query)? else {
                    // Would fail the coverage gate; leave the site alone.
                    why = Some(|l| l.coverage_gate += 1);
                    continue;
                };
                // The window arithmetic must land exactly, the same check `find_alts` makes on
                // its own candidates. The restitched window is the sequence that would be
                // published, so if it is not the current window's length plus the delta, the
                // shift was wrong and the swap would truncate or reframe `_primary`, the record
                // the codon alignment and selection analysis read.
                //
                // Compare windows, not whole paths. `spliced_len` is the length of the entire
                // stitched path and belongs to the indel rules above, which are about how much
                // sequence the detour adds or removes; it is not what `stitch` returns.
                if restitched.seq.len() as i64 != window.len() as i64 + delta {
                    why = Some(|l| l.window_not_shiftable += 1);
                    continue;
                }
                // A different route through identical sequence is not a correction, and counting
                // it as one would print a `consensus:` line for a no-op. `find_alts` drops the
                // same case as `identical`.
                if restitched.seq == window {
                    why = Some(|l| l.identical += 1);
                    continue;
                }
                // A length-changing swap moves the baseline every later swap is measured
                // against, so it has to be updated rather than left at the original route's. The
                // window moves with it, for the same reason.
                baseline_len = spliced_len;
                window = restitched.seq;
                current = candidate;
                cumulative_delta += delta;
                swaps += 1;
                if taken_minor {
                    // Counted apart from an ordinary two-allele correction, because only
                    // this one fires under `--rescue-minority`.
                    report.minority_rescues += 1;
                }
                applied = true;
                break;
            }
            if applied {
                break;
            }
            // Record the refusal for ANY fork a correction was tried on, not only minority
            // ones. The cumulative indel cap and the rejoin and window rules fire on
            // ordinary two-allele sites too, and a refusal nobody counts is invisible.
            if let Some(bump) = why {
                bump(&mut left);
            }
        }
        if !applied {
            report.left = left;
            hit_cap = false;
            break;
        }
    }
    // Exhausting the loop means corrections were still being applied when the backstop
    // stopped us, so this route is not a few swapped alleles. Silence here is what let the
    // walker20k slice pin 1,460 accessions at the cap unnoticed.
    report.hit_pass_cap = hit_cap && swaps > 0;
    if report.hit_pass_cap {
        eprintln!(
            "[{}] WARNING: consensus correction stopped at the {} pass backstop with swaps \
             still being applied. This route was rewritten rather than having a few alleles \
             corrected; treat the primary as provisional.",
            subgraph.accession, passes
        );
    }

    report.swaps = swaps;
    if report.left.total() > 0 {
        // Not a warning: with `--prune-frac` refusing a prune, or an indel the consensus will
        // not take, leaving a site is the correct outcome. It is reported because `_primary`
        // carries the rarer allele there and nothing else in the output would say so.
        eprintln!(
            "[{}] consensus: {} fork(s) left where a correction was considered ({} alt below \
             floor, {} no outright majority, {} walked-minority exempt, {} fork \
             saturated, {} no rejoin, {} length refused, {} window, {} coverage, \
             {} identical)",
            subgraph.accession,
            report.left.total(),
            report.left.alt_below_floor,
            report.left.not_outright_majority,
            report.left.exempt_minority,
            report.left.saturated,
            report.left.no_rejoin,
            report.left.length_refused,
            report.left.window_not_shiftable,
            report.left.coverage_gate,
            report.left.identical,
        );
    }
    if swaps == 0 {
        return Ok((None, report));
    }
    eprintln!(
        "[{}] consensus: {} site(s) where the route carried the minority allele, corrected \
         ({} of them a minority the old gate refused to touch)",
        subgraph.accession, swaps, report.minority_rescues
    );
    // A handful of corrections is the expected shape. A large number means the route has
    // been rerouted through a substantially different sequence rather than having a few
    // alleles swapped, and that is worth a human look before the result is used.
    //
    // Measured over a 200-accession panel: 26 accessions were corrected, 25 of them at 1 to
    // 11 sites and cleanly. The exception, `SRR19878648` (H5N1 HA), took 27 corrections that
    // moved 91 positions in 11 blocks across the gene, one of them 339 nt wide, and put four
    // gaps into an alignment that previously had none. Its route also dips to 13.5x in a
    // 230x sample, so the majority call at those sites rests on very few reads.
    // The band covers ordinary corrections only. Minority rescues are legitimate in bulk on
    // a clean accession, so they are counted beside it rather than folded into it.
    let ordinary = swaps - report.minority_rescues;
    if ordinary > 12 {
        eprintln!(
            "[{}] WARNING: {} ordinary consensus corrections (plus {} minority rescues) \
             is far more than the 1-11 seen across the panel. The route may have been rerouted \
             through a divergent segment rather than having single alleles swapped; treat this \
             accession as provisional and compare it against the unswapped result.",
            subgraph.accession, ordinary, report.minority_rescues
        );
    }
    Ok((Some(current), report))
}

/// Recover one alternative haplotype per genuinely biallelic site on the walked path.
pub fn find_alts(
    cfg: &Config,
    graph: &UnitigGraph,
    subgraph: &Subgraph,
    path: &GafPath,
    query: &Query,
) -> Result<(Vec<AltHaplotype>, VariantStats)> {
    // The fork scan runs whatever the output flags say; only the splicing is conditional.
    // The scan decides `mixed=` and `min_vaf=` on the primary's own header, and a mixture
    // flag computed only when alternatives were asked for is worse than none:
    // `--no-emit-alts` would publish `mixed=0` for a sequence nobody looked at.
    let want_alts = cfg.emit_alts || cfg.emit_vcf || cfg.emit_snv_table;
    if path.steps.len() < 2 || subgraph.nodes.is_empty() {
        return Ok((Vec::new(), VariantStats::default()));
    }

    // Shares `view()` with `prefer_majority` rather than rebuilding it: an adjacency change
    // landing in one copy would make the corrected primary and the alts disagree about which
    // forks exist, which no output reveals. `phase` still builds its own, differently shaped.
    let Some(View { out, walk, at }) = view(graph, subgraph, path) else {
        // The path leaves the subgraph; `stitch` would bail too. Nothing was examined, so
        // the stats stay zero rather than claiming forks that were never looked at.
        return Ok((Vec::new(), VariantStats::default()));
    };

    let Some(primary) = stitch::stitch(cfg, graph, path, subgraph, query)? else {
        return Ok((Vec::new(), VariantStats::default()));
    };
    let baseline_len = stitched_len(cfg, &path.steps, subgraph);

    let abund = |o: Oriented| graph.abundance_of(subgraph.nodes[(o as usize) / 2]);
    let depth = path_depth(graph, subgraph, &walk);
    let floor = abundance_floor(cfg, depth);

    // Counters, so a zero result can be told apart from a filtered one. Calibration needs
    // to know which gate is doing the work, not just the total.
    let (mut forks, mut biallelic, mut no_rejoin, mut not_substitution, mut identical) =
        (0usize, 0usize, 0usize, 0usize, 0usize);
    // Sites, not alleles, and the weakest allele the published route carries.
    let mut two_allele_sites = 0usize;
    let mut min_site_vaf = 1.0f64;
    let (mut rejected_large, mut rejected_frameshift, mut rejected_window, mut indels) =
        (0usize, 0usize, 0usize, 0usize);
    let mut deltas: Vec<i64> = Vec::new();
    let mut alts: Vec<AltHaplotype> = Vec::new();
    // The scan runs to the end and counts what it leaves out. Breaking at the cap would
    // give output for the 5' end of the gene only, with the rest looking invariant and the
    // stderr counters stopping too, so nothing would say the result was truncated.
    let mut beyond_cap = 0usize;
    for i in 0..walk.len() - 1 {
        let fork = walk[i];
        let taken = walk[i + 1];
        let branches = &out[fork as usize];
        if branches.len() < 2 || !branches.contains(&taken) {
            continue;
        }
        forks += 1;
        // As in `prefer_majority`: only branches that rejoin count toward the total.
        let total = site_total(&out, &at, branches, taken, i, cfg.alt_max_detour, abund);
        if total <= 0.0 {
            continue;
        }
        let taken_ab = abund(taken);
        let taken_vaf = taken_ab as f64 / total;
        // The branch the walk took must itself be a credible allele. If it is not, this is
        // the reference-pull defect and not a site with two populations.
        if taken_ab < floor || taken_vaf < cfg.min_vaf {
            continue;
        }

        let mut second_allele_here = false;
        for &alt in branches {
            if alt == taken {
                continue;
            }
            let alt_ab = abund(alt);
            let alt_vaf = alt_ab as f64 / total;
            if alt_ab < floor || alt_vaf < cfg.min_vaf {
                continue;
            }
            biallelic += 1;
            let Some((rejoin_at, detour)) = rejoin(&out, &at, alt, i, cfg.alt_max_detour) else {
                no_rejoin += 1;
                continue;
            };

            // Splice: everything up to the fork, the detour, then the path from where it
            // rejoins. Every step is a real edge, so this is still a connected path.
            let mut steps: Vec<PathStep> = Vec::with_capacity(path.steps.len() + detour.len());
            steps.extend_from_slice(&path.steps[..=i]);
            for &o in &detour {
                steps.push(PathStep {
                    node: subgraph.nodes[(o as usize) / 2],
                    orient: orient_of(o),
                });
            }
            steps.extend_from_slice(&path.steps[rejoin_at..]);

            // The length rules start here. See the module note "Which alleles are accepted":
            // a substitution always, an in-frame indel inside `max_alt_indel` by default, and
            // nothing that rejoins far downstream, which would be a shortcut deleting gene.
            let spliced_len = stitched_len(cfg, &steps, subgraph);
            let delta = spliced_len as i64 - baseline_len as i64;

            // A mixture site is one where the sample holds two POPULATIONS, and in a de
            // Bruijn graph a length-changing branch is usually not one: 190 of 203
            // length-changing calls across the panel are plus or minus one base in a
            // homopolymer, and the large ones are the route being sent through a different
            // stretch of graph. Only an equal-length rejoining branch counts towards `mixed=`.
            //
            // Measured on the gated fixtures: `ERR2028820` has four branches that clear the
            // abundance and frequency floors, and all four are length-changing or never
            // rejoin (2 over the 30 nt indel cap, 1 frameshifting, 1 no rejoin). Counting
            // those would have labelled the repository's clonal control a mixture, which is
            // the measurement that set this rule.
            if delta == 0 {
                second_allele_here = true;
            }
            // Everything below builds the allele as a sequence, which only an output needs.
            if !want_alts {
                continue;
            }
            // Counted after the gates, so this is "credible alleles we could not emit", not
            // "branches we did not look at". After the rejoin test too, so it means
            // an allele the cap displaced rather than one that was never emittable.
            if alts.len() >= cfg.max_alts {
                beyond_cap += 1;
                continue;
            }

            // An allele that changes the length is an indel, and three things have to hold
            // before one is trustworthy enough to emit. See `shifted_window` for the first.
            if delta != 0 {
                if !cfg.allow_indel_alts {
                    not_substitution += 1;
                    deltas.push(delta);
                    continue;
                }
                // Large deltas are not indel alleles. Measured across the panel, the
                // length-changing calls cluster at +-1 and +-2 with a long tail at +100,
                // +142, -282: the tail is the route being sent through a different stretch
                // of graph, which rewrites the gene rather than editing it.
                if delta.unsigned_abs() as usize > cfg.max_alt_indel {
                    rejected_large += 1;
                    continue;
                }
                // A delta that is not a multiple of three shifts the reading frame of
                // everything downstream. Real frameshift alleles exist, but they are
                // indistinguishable here from the homopolymer-length artifacts a de Bruijn
                // graph produces, and 94% of the panel's length-changing calls are this
                // case with a median magnitude of one base. In a codon-aware analysis a
                // frameshifted record is worse than a missing one, so it is off by default.
                if delta % 3 != 0 && !cfg.allow_frameshift_alts {
                    rejected_frameshift += 1;
                    continue;
                }
            }

            let prefix_len = stitched_len(cfg, &path.steps[..=i], subgraph);
            let replaced_len = tail_contrib(cfg, &path.steps[i + 1..rejoin_at], subgraph);
            let Some((new_start, new_end)) =
                shifted_window(path.path_start, path.path_end, prefix_len, replaced_len, delta)
            else {
                rejected_window += 1;
                continue;
            };
            let spliced = GafPath {
                steps,
                path_start: new_start,
                path_end: new_end,
                ..copy_coords(path)
            };
            let Some(mut hap) = stitch::stitch(cfg, graph, &spliced, subgraph, query)? else {
                continue;
            };
            if hap.seq == primary.seq {
                identical += 1; // a different route through identical sequence is not a variant
                continue;
            }
            // The window arithmetic must land exactly. If the emitted sequence is not the
            // primary's length plus the delta, the shift was wrong and the record would be
            // truncated or reframed, so it is dropped rather than published.
            let expected = primary.seq.len() as i64 + delta;
            if hap.seq.len() as i64 != expected {
                rejected_window += 1;
                continue;
            }
            if delta != 0 {
                indels += 1;
            }
            // How much sequence this allele actually replaces: trim the identical head and
            // tail, and what is left between them is the edit.
            let (a, b) = (&hap.seq, &primary.seq);
            let head = a.iter().zip(b.iter()).take_while(|(x, y)| x == y).count();
            let tail = a
                .iter()
                .rev()
                .zip(b.iter().rev())
                .take_while(|(x, y)| x == y)
                .count()
                .min(a.len().min(b.len()).saturating_sub(head));
            let replaced_span = a.len().max(b.len()).saturating_sub(head + tail);
            hap.id = format!("{}_alt{}", subgraph.accession, alts.len() + 1);
            alts.push(AltHaplotype {
                // The VAFs just computed rest on a total that is a lower bound.
                ka_saturated: branches.iter().any(|&b| abund(b) >= KA_CEILING),
                haplotype: hap,
                fork: subgraph.nodes[(fork as usize) / 2],
                fork_strand: strand_of(fork),
                alt_node: subgraph.nodes[(alt as usize) / 2],
                alt_strand: strand_of(alt),
                taken_abundance: taken_ab,
                taken_vaf,
                alt_abundance: alt_ab,
                alt_vaf,
                detour_nodes: detour.len(),
                length_delta: delta,
                replaced_span,
            });
        }
        // One site, however many substitution alleles it holds, and the frequency of the one
        // the primary took: two such sites more than k - 1 apart are independent choices that
        // no read links, so the published sequence may exist in no molecule.
        if second_allele_here {
            two_allele_sites += 1;
            min_site_vaf = min_site_vaf.min(taken_vaf);
        }
    }
    eprintln!(
        "[{}] variants: depth {:.1}x, floor {:.1}x; {} on-path fork(s), \
         {} biallelic candidate(s) -> {} recovered \
         ({} indel allele(s); {} no rejoin within {}, {} length-changing refused{}, \
         {} over the {} nt indel cap, {} frameshifting, {} window not shiftable, \
         {} identical to primary)",
        subgraph.accession, depth, floor, forks, biallelic, alts.len(),
        indels, no_rejoin, cfg.alt_max_detour, not_substitution,
        if deltas.is_empty() { String::new() } else { format!(" {deltas:?}") },
        rejected_large, cfg.max_alt_indel, rejected_frameshift, rejected_window,
        identical
    );
    if !want_alts {
        eprintln!(
            "[{}] variants: no variant output requested, so alleles were counted and not \
             built; the {} two-allele site(s) on the route are still reported as mixed= and \
             min_vaf= on the primary",
            subgraph.accession, two_allele_sites
        );
    }
    // Say it, rather than leaving the caller to infer truncation from a bubble count.
    if beyond_cap > 0 {
        eprintln!(
            "[{}] WARNING: {} credible alternative allele(s) beyond --max-alts {} were \
             omitted. alt_haplotypes.fasta, the VCF and snvs.tsv are incomplete for this \
             accession, and the omitted sites are the 3' ones: raise --max-alts to recover \
             them.",
            subgraph.accession, beyond_cap, cfg.max_alts
        );
    }
    // The same numbers, for `summary.json`. They are assembled here rather than counted
    // again by the caller so the printed line and the recorded block cannot drift apart.
    let stats = VariantStats {
        depth,
        floor,
        forks,
        biallelic,
        two_allele_sites,
        min_site_vaf,
        recovered: alts.len(),
        beyond_cap,
        indels,
        no_rejoin,
        not_substitution,
        over_indel_cap: rejected_large,
        frameshifting: rejected_frameshift,
        window_not_shiftable: rejected_window,
        identical,
    };
    Ok((alts, stats))
}

/// The gene's own sequencing depth: median `ka:f` over the nodes the walk uses.
///
/// The walked path is the gene, so this measures the quantity the floor is relative to.
/// Median rather than mean because a single repeat-collapsed unitig at 25,000x should not
/// drag the estimate for an otherwise ordinary accession.
fn path_depth(graph: &UnitigGraph, subgraph: &Subgraph, walk: &[Oriented]) -> f32 {
    // The same median `extract` takes for the prune floor, over the oriented walk
    // rather than over `PathStep`s. Only the projection differs, so only it lives here.
    crate::stages::extract::median_abundance(
        walk.iter()
            .map(|&o| graph.abundance_of(subgraph.nodes[(o as usize) / 2]))
            .collect(),
    )
}

/// The abundance a branch must clear to count as an allele.
///
/// Never below `min_abund`: this can only tighten the existing gate, never loosen it.
fn abundance_floor(cfg: &Config, depth: f32) -> f32 {
    // `extract::depth_floor` with `--min-abund-frac` in place of `--prune-frac`.
    crate::stages::extract::relative_floor(cfg.min_abund_frac, depth, cfg.min_abund)
}

fn strand_of(o: Oriented) -> Strand {
    if o % 2 == 1 {
        Strand::Forward
    } else {
        Strand::Reverse
    }
}

fn copy_coords(p: &GafPath) -> GafPath {
    GafPath {
        steps: Vec::new(),
        query_start: p.query_start,
        query_end: p.query_end,
        path_len: p.path_len,
        path_start: p.path_start,
        path_end: p.path_end,
        num_matches: p.num_matches,
        aln_len: p.aln_len,
    }
}

/// Length of the sequence a path stitches to, without building it.
/// Bases a run of steps contributes when it sits *after* the start of a path.
///
/// Every step but the first contributes `len - (k-1)`, so a stretch in the middle of a route
/// has a length that does not depend on what precedes it. That is what makes a splice's
/// effect on the total computable without restitching.
fn tail_contrib(cfg: &Config, steps: &[PathStep], subgraph: &Subgraph) -> usize {
    let overlap = cfg.graph_k.saturating_sub(1);
    let mut total = 0usize;
    for step in steps {
        let Ok(idx) = subgraph.nodes.binary_search(&step.node) else {
            return usize::MAX;
        };
        total += subgraph.seqs[idx].len().saturating_sub(overlap);
    }
    total
}

/// Where the alignment window lands after a splice that changes the sequence length.
///
/// `path_start` and `path_end` are offsets into the stitched sequence, so a splice that adds
/// or removes bases moves everything after it. Reusing the old `path_end` would slice the
/// wrong window and silently reframe the coding sequence, which is why length-changing
/// alternatives were rejected outright until this existed.
///
/// The shift is only valid when the splice sits entirely inside the window: the unchanged
/// prefix must reach at least `path_start`, and the stretch being replaced must end at or
/// before `path_end`. A splice straddling either edge is not a simple shift, and is refused.
///
/// Returns the new `(start, end)`, or `None` when the shift cannot be justified.
fn shifted_window(
    path_start: usize,
    path_end: usize,
    prefix_len: usize,
    replaced_len: usize,
    delta: i64,
) -> Option<(usize, usize)> {
    if delta == 0 {
        return Some((path_start, path_end));
    }
    if prefix_len < path_start {
        return None; // the splice begins before the window; the prefix is not intact
    }
    if prefix_len + replaced_len > path_end {
        return None; // the replaced stretch runs past the window's end
    }
    let end = path_end as i64 + delta;
    if end <= path_start as i64 {
        return None;
    }
    Some((path_start, end as usize))
}

fn stitched_len(cfg: &Config, steps: &[PathStep], subgraph: &Subgraph) -> usize {
    stitch::stitched_len(cfg, steps, subgraph)
}

/// Shortest way from `alt` back onto the walked path, strictly after the fork.
///
/// Returns the path index it rejoins at and the detour, excluding the rejoin node itself
/// because the path supplies that. Breadth-first, so the detour found is the shortest, and
/// bounded so a tangled neighbourhood cannot make this expensive.
/// Total abundance at a fork, counting only the branches that come back to the route.
///
/// An exit that leaves the locus for good, a repeat, a paralog or a contaminant sharing the
/// fork's 30-mer, is not an allele of this gene, and counting it in the denominator deflates
/// every frequency at the site. With the route at 80x, the allele at 20x and an exit at 400x,
/// the allele reads 0.04 instead of 0.2 and clears no gate.
///
/// `phase` already restricts its own total to in-bubble branches and says so
/// (`phase.rs`), which is why `bubbles.json` and the VCF could report different VAFs for the
/// same site. They now agree on what belongs in the denominator.
///
/// The taken branch is always counted: it is the route, so it rejoins by definition.
fn site_total<F: Fn(Oriented) -> f32>(
    out: &[Vec<Oriented>],
    at: &HashMap<Oriented, usize>,
    branches: &[Oriented],
    taken: Oriented,
    fork_index: usize,
    max_detour: usize,
    abund: F,
) -> f64 {
    branches
        .iter()
        .filter(|&&b| b == taken || rejoin(out, at, b, fork_index, max_detour).is_some())
        .map(|&b| abund(b) as f64)
        .sum()
}

fn rejoin(
    out: &[Vec<Oriented>],
    at: &HashMap<Oriented, usize>,
    alt: Oriented,
    fork_index: usize,
    max_detour: usize,
) -> Option<(usize, Vec<Oriented>)> {
    let mut queue: std::collections::VecDeque<(Oriented, Vec<Oriented>)> =
        std::collections::VecDeque::new();
    let mut seen: std::collections::HashSet<Oriented> = std::collections::HashSet::new();
    seen.insert(alt);
    queue.push_back((alt, vec![alt]));

    while let Some((node, trail)) = queue.pop_front() {
        if let Some(&j) = at.get(&node) {
            if j > fork_index + 1 {
                let mut detour = trail;
                detour.pop(); // the rejoin node belongs to the path, not the detour
                return Some((j, detour));
            }
        }
        if trail.len() > max_detour {
            continue;
        }
        for &next in &out[node as usize] {
            if seen.insert(next) {
                let mut t = trail.clone();
                t.push(next);
                queue.push_back((next, t));
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stages::parse;
    use crate::types::PathStep;
    use std::path::PathBuf;

    // Both stages gained a second return value in 0.2.0, the counts `summary.json` now
    // records (`docs/summaries.md` rule 3). These shims shadow the real functions so the
    // tests below keep asserting on what they are about, the alleles and the route. A test
    // that is about the counts calls `super::` directly.
    fn find_alts(
        cfg: &Config,
        graph: &UnitigGraph,
        subgraph: &Subgraph,
        path: &GafPath,
        query: &Query,
    ) -> Result<Vec<AltHaplotype>> {
        super::find_alts(cfg, graph, subgraph, path, query).map(|(alts, _)| alts)
    }

    fn prefer_majority(
        cfg: &Config,
        graph: &UnitigGraph,
        subgraph: &Subgraph,
        path: &GafPath,
        query: &Query,
    ) -> Result<Option<GafPath>> {
        super::prefer_majority(cfg, graph, subgraph, path, query).map(|(path, _)| path)
    }

    /// Two diamonds in series: `0 -> {1,2} -> 3 -> {4,5} -> 6`, k=4. The walk takes 1 and 4.
    ///
    /// The two sites are joined only through node 3, which both strains carry,
    /// so every combination of their alleles is a connected path and nothing in the graph
    /// says which pairs occur together in a molecule.
    fn two_diamonds(a_taken: f32, a_alt: f32, b_taken: f32, b_alt: f32) -> UnitigGraph {
        let fa = format!(
            ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n\
             >T_1 ka:f:{a_taken}  L:-:0:-  L:+:3:+\nCGTA\n\
             >T_2 ka:f:{a_alt}  L:-:0:-  L:+:3:+\nCGTC\n\
             >T_3 ka:f:100.0  L:-:1:-  L:-:2:-  L:+:4:+  L:+:5:+\nGTAG\n\
             >T_4 ka:f:{b_taken}  L:-:3:-  L:+:6:+\nTAGA\n\
             >T_5 ka:f:{b_alt}  L:-:3:-  L:+:6:+\nTAGC\n\
             >T_6 ka:f:100.0  L:-:4:-  L:-:5:-\nAGAT\n"
        );
        parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.into_bytes())), u64::MAX).unwrap()
    }

    /// The route through both diamonds, taking branch 1 and then branch 4.
    fn walked_two() -> GafPath {
        GafPath {
            steps: vec![
                PathStep { orient: Orient::Forward, node: 0 },
                PathStep { orient: Orient::Forward, node: 1 },
                PathStep { orient: Orient::Forward, node: 3 },
                PathStep { orient: Orient::Forward, node: 4 },
                PathStep { orient: Orient::Forward, node: 6 },
            ],
            query_start: 0,
            query_end: 8,
            path_len: 8,
            path_start: 0,
            path_end: 8,
            num_matches: 99,
            aln_len: 100,
        }
    }

    /// A primary that passes two two-allele sites is a combination of alleles nothing links:
    /// the sites are joined through unitigs both strains carry, so the walk chooses at each
    /// independently and may publish a sequence carried by no molecule.
    #[test]
    fn two_sites_on_the_route_are_counted_and_the_weakest_allele_is_reported() {
        let g = two_diamonds(60.0, 40.0, 55.0, 45.0);
        let sub = subgraph_of(&g);
        let (_, stats) =
            super::find_alts(&cfg(true), &g, &sub, &walked_two(), &q("ACGTAGAT")).expect("find_alts");
        assert_eq!(stats.forks, 2, "the route passes two forks");
        assert_eq!(stats.two_allele_sites, 2, "both carry a credible second allele");
        assert!(
            (stats.min_site_vaf - 0.55).abs() < 1e-9,
            "the weakest allele the primary carries is the 55% one, got {}",
            stats.min_site_vaf
        );
        // Sites, not alleles: `biallelic` counts one per alternative and is the number the
        // alt writer needs. Here they agree, which is what makes the distinction checkable
        // in the three-branch case below.
        assert_eq!(stats.biallelic, 2);
    }

    /// A site is one fork however many alleles hang off it, and a clonal route has none.
    #[test]
    fn a_three_branch_fork_is_one_site_and_a_clonal_route_is_none() {
        // Node 0 forks three ways: 1 (taken, 50), 2 (30) and a third at 20, all rejoining.
        let fa = ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+  L:+:4:+\nACGT\n\
                  >T_1 ka:f:50.0  L:-:0:-  L:+:3:+\nCGTA\n\
                  >T_2 ka:f:30.0  L:-:0:-  L:+:3:+\nCGTC\n\
                  >T_3 ka:f:100.0  L:-:1:-  L:-:2:-  L:-:4:-\nGTAG\n\
                  >T_4 ka:f:20.0  L:-:0:-  L:+:3:+\nCGTG\n";
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX).unwrap();
        let sub = subgraph_of(&g);
        let (_, stats) =
            super::find_alts(&cfg(true), &g, &sub, &walked(), &q("ACGTAG")).expect("find_alts");
        assert_eq!(stats.two_allele_sites, 1, "one fork is one site");
        assert_eq!(stats.biallelic, 2, "two alternatives hang off it");
        assert!((stats.min_site_vaf - 0.5).abs() < 1e-9);

        // The clonal case: the alternative is sequencing error, under the floor.
        let g = diamond(100.0, 1.0);
        let sub = subgraph_of(&g);
        let (_, stats) =
            super::find_alts(&cfg(true), &g, &sub, &walked(), &q("ACGTAG")).expect("find_alts");
        assert_eq!(stats.two_allele_sites, 0);
        assert_eq!(stats.min_site_vaf, 1.0, "no site means no weakest allele, which is 1.0");
    }

    /// The mixture scan must not depend on the output flags: it decides `mixed=` on the
    /// primary's own header, so skipping it under `--no-emit-alts` would publish `mixed=0`
    /// for a sequence nothing examined.
    #[test]
    fn the_mixture_is_measured_even_when_no_variant_output_is_wanted() {
        let g = two_diamonds(60.0, 40.0, 55.0, 45.0);
        let sub = subgraph_of(&g);
        let quiet = Config { emit_alts: false, emit_vcf: false, emit_snv_table: false, ..cfg(true) };
        let (alts, stats) =
            super::find_alts(&quiet, &g, &sub, &walked_two(), &q("ACGTAGAT")).expect("find_alts");
        assert!(alts.is_empty(), "no allele file was asked for");
        assert_eq!(stats.two_allele_sites, 2, "the sites are still counted");
        assert!((stats.min_site_vaf - 0.55).abs() < 1e-9);
    }

    /// A diamond: node 0 forks to either 1 or 2, both rejoin at 3. k=4, so overlap is 3 and
    /// each branch contributes one base. The walked path takes branch 1.
    ///
    /// `ka:f` is set per test to move the site between the regimes that matter.
    fn diamond(taken: f32, alt: f32) -> UnitigGraph {
        let fa = format!(
            ">T_0 ka:f:{main}  L:+:1:+  L:+:2:+\nACGT\n\
             >T_1 ka:f:{taken}  L:-:0:-  L:+:3:+\nCGTA\n\
             >T_2 ka:f:{alt}  L:-:0:-  L:+:3:+\nCGTC\n\
             >T_3 ka:f:{main}  L:-:1:-  L:-:2:-\nGTAG\n",
            main = 100.0, taken = taken, alt = alt
        );
        parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.into_bytes())), u64::MAX).unwrap()
    }

    fn subgraph_of(g: &UnitigGraph) -> Subgraph {
        let nodes: Vec<UnitigId> = (0..g.n_unitigs as u32).collect();
        let seqs = crate::stages::extract::collect_seqs(g, &nodes).unwrap();
        Subgraph {
            accession: "ACC".into(),
            nodes,
            seqs,
            gfa_path: PathBuf::from("/tmp/logan-variants-test-nonexistent.gfa"),
            n_edges: 0,
            owns_gfa: false,
        }
    }
    /// The walked route with an explicit alignment window, so a test can put the window
    /// inside the path as a real GAF row does. `walked()` is the spanning case.
    fn walked_window(path_start: usize, path_end: usize) -> GafPath {
        GafPath {
            steps: vec![
                PathStep { orient: Orient::Forward, node: 0 },
                PathStep { orient: Orient::Forward, node: 1 },
                PathStep { orient: Orient::Forward, node: 3 },
            ],
            query_start: 0,
            query_end: 7,
            path_len: 7,
            path_start,
            path_end,
            num_matches: 99,
            aln_len: 100,
        }
    }
    /// A window spanning the whole path. Convenient, but NOT representative: see
    /// `the_majority_is_swapped_in_when_the_window_is_inside_the_path` for why that matters.
    fn walked() -> GafPath {
        walked_window(0, 7)
    }
    fn q(seq: &str) -> Query {
        Query {
            name: "q".into(),
            record_id: "q".into(),
            seq: seq.as_bytes().to_vec(),
            path: PathBuf::from("q.fa"),
            n_records: 1,
        }
    }
    fn cfg(emit: bool) -> Config {
        Config { graph_k: 4, min_cov: 0.0, emit_alts: emit, ..Config::default() }
    }

    fn maj(on: bool) -> Config {
        Config { graph_k: 4, min_cov: 0.0, prefer_majority: on, ..Config::default() }
    }

    #[test]
    fn majority_correction_can_be_switched_off() {
        assert!(Config::default().prefer_majority, "on by default since the panel run");
        let g = diamond(40.0, 60.0);
        let s = subgraph_of(&g);
        assert!(prefer_majority(&maj(false), &g, &s, &walked(), &q("ACGTAG")).unwrap().is_none());
    }

    /// The route took the 40% allele because it matched the reference better. Both alleles
    /// are real, so the primary should carry the 60% one.
    #[test]
    fn the_primary_takes_the_more_abundant_allele() {
        let g = diamond(40.0, 60.0);
        let s = subgraph_of(&g);
        let query = q("ACGTAG");
        let before = stitch::stitch(&maj(true), &g, &walked(), &s, &query).unwrap().unwrap();
        let fixed = prefer_majority(&maj(true), &g, &s, &walked(), &query)
            .unwrap()
            .expect("the site should be corrected");
        let after = stitch::stitch(&maj(true), &g, &fixed, &s, &query).unwrap().unwrap();
        assert_ne!(after.seq, before.seq);
        assert_eq!(after.seq.len(), before.seq.len(), "substitution, not indel");
        assert!(fixed.steps.iter().any(|p| p.node == 2), "must now traverse the 60% branch");
        assert!(!fixed.steps.iter().any(|p| p.node == 1), "and not the 40% one");
    }

    #[test]
    fn a_route_already_carrying_the_majority_is_left_alone() {
        let g = diamond(60.0, 40.0);
        let s = subgraph_of(&g);
        assert!(prefer_majority(&maj(true), &g, &s, &walked(), &q("ACGTAG")).unwrap().is_none());
    }

    /// A fork whose walked branch is a minority is reported and, by default, left alone.
    /// What would stop acting on it becoming a greedy whole-route rewrite is the bar on the
    /// other side: an outright majority of the fork, here 12,000 of 12,003.
    #[test]
    fn an_error_branch_site_is_reported_but_not_touched() {
        let g = diamond(3.0, 12000.0);
        let s = subgraph_of(&g);
        let (fixed, report) =
            super::prefer_majority(&maj(true), &g, &s, &walked(), &q("ACGTAG")).unwrap();
        assert!(fixed.is_none(), "the exemption is the default again");
        assert_eq!(report.swaps, 0);
        // The site is first-class data even though nothing acts on it.
        assert_eq!(report.left.exempt_minority, 1, "{:?}", report.left);
        assert_eq!(report.left.total(), 1);
    }

    /// `--rescue-minority` acts on the site instead of only reporting it. Off by default
    /// because the correction loop does not converge with it on, but it stays reachable.
    #[test]
    fn rescue_minority_moves_the_route_off_a_walked_minority_branch() {
        let g = diamond(3.0, 12000.0);
        let s = subgraph_of(&g);
        let c = Config { rescue_minority: true, ..maj(true) };
        let (fixed, report) =
            super::prefer_majority(&c, &g, &s, &walked(), &q("ACGTAG")).unwrap();
        let fixed = fixed.expect("0.02% walked against a 99.98% majority");
        assert!(fixed.steps.iter().any(|p| p.node == 2), "must traverse the majority branch");
        assert_eq!((report.swaps, report.minority_rescues), (1, 1));
        assert_eq!(report.left.total(), 0);
    }

    /// The flag is set by the ceiling exactly, not by "large".
    ///
    /// 24,999 is an ordinary high-coverage fork whose ratios mean what they say. 25,000 is
    /// Logan's clip, so the fork total is a lower bound and both VAFs on the record are
    /// overstated by an unknown amount.
    #[test]
    fn an_allele_is_flagged_saturated_only_at_the_ceiling() {
        let near = diamond(24_999.0, 20_000.0);
        let s = subgraph_of(&near);
        let (alts, _) = super::find_alts(&cfg(true), &near, &s, &walked(), &q("ACGTAG")).unwrap();
        assert_eq!(alts.len(), 1, "a two-allele site either way");
        assert!(!alts[0].ka_saturated, "24,999 is below the ceiling");

        let at = diamond(25_000.0, 20_000.0);
        let s = subgraph_of(&at);
        let (alts, _) = super::find_alts(&cfg(true), &at, &s, &walked(), &q("ACGTAG")).unwrap();
        assert_eq!(alts.len(), 1);
        assert!(alts[0].ka_saturated, "25,000 is the clip; the VAFs are lower bounds");
    }

    /// The cap is a budget for the whole correction set, not for one swap.
    ///
    /// Found on a panel, not by reasoning: one accession took 64 corrections that were each
    /// in frame and each inside the 30 nt cap, and its ebola GP primary grew from 2,049 nt to
    /// 3,129 against a 2,031 nt reference. Every step was legal; the total was the
    /// whole-route rewrite this stage exists to avoid. Here two +30 nt swaps are each legal
    /// and the second must still be refused.
    #[test]
    fn the_indel_cap_bounds_the_whole_correction_set_not_one_swap() {
        // Two diamonds in series, each with a +30 nt in-frame majority branch.
        let ins = "A".repeat(30);
        // Ordinary two-allele sites (40 against 60), not minority rescues, so this tests the
        // cumulative cap on the default path rather than behind `--rescue-minority`.
        let fa = format!(
            ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n\
             >T_1 ka:f:40.0   L:-:0:-  L:+:3:+\nCGTA\n\
             >T_2 ka:f:60.0   L:-:0:-  L:+:3:+\nCGT{ins}A\n\
             >T_3 ka:f:100.0  L:-:1:-  L:-:2:-  L:+:4:+  L:+:5:+\nGTAC\n\
             >T_4 ka:f:40.0   L:-:3:-  L:+:6:+\nTACG\n\
             >T_5 ka:f:60.0   L:-:3:-  L:+:6:+\nTAC{ins}G\n\
             >T_6 ka:f:100.0  L:-:4:-  L:-:5:-\nACGT\n",
            ins = ins
        );
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.into_bytes())), u64::MAX)
            .unwrap();
        let s = subgraph_of(&g);
        let walk = GafPath {
            steps: vec![
                PathStep { orient: Orient::Forward, node: 0 },
                PathStep { orient: Orient::Forward, node: 1 },
                PathStep { orient: Orient::Forward, node: 3 },
                PathStep { orient: Orient::Forward, node: 4 },
                PathStep { orient: Orient::Forward, node: 6 },
            ],
            ..walked()
        };
        let c = Config { min_cov: 0.0, ..maj(true) };
        let (fixed, report) =
            super::prefer_majority(&c, &g, &s, &walk, &q("ACGTACGT")).unwrap();
        assert!(fixed.is_some(), "the first +30 is inside the cap and must be taken");
        assert_eq!(report.swaps, 1, "the second would put the total at 60, over the 30 cap");
        assert_eq!(report.left.length_refused, 1, "and it is reported, not silent: {:?}",
                   report.left);
    }

    /// The case proper: a reference-like minority between the prune floor and `--min-vaf`,
    /// which neither mechanism reaches. Reproduced at 2%, 5% and 9% of the fork.
    #[test]
    fn a_walked_minority_between_the_prune_floor_and_min_vaf_is_corrected() {
        // 5% of the fork: above the 1%-of-depth prune floor, below --min-vaf 0.10.
        let g = diamond(50.0, 950.0);
        let s = subgraph_of(&g);
        let c = Config { rescue_minority: true, ..maj(true) };
        let (fixed, report) =
            super::prefer_majority(&c, &g, &s, &walked(), &q("ACGTAG")).unwrap();
        assert!(fixed.is_some(), "5% walked against a 95% majority must be corrected");
        assert_eq!((report.swaps, report.minority_rescues), (1, 1));
    }

    /// No minority rescue where the fork is clipped.
    ///
    /// The rescue path's whole safety rail is "the other side holds an outright majority of
    /// the fork", and that is a ratio against the fork total. Where a branch sits at Logan's
    /// 25,000 ceiling the total is a lower bound, so the ratio is overstated and the rail is
    /// not where it claims to be. Measured on 5,000 walker20k accessions, 62% of them
    /// saturated: without this guard corrections went from 470 sites to 97,068, the median
    /// corrected accession hit the loop's own backstop, and 62 gained an internal stop codon
    /// while none lost one.
    #[test]
    fn a_minority_rescue_is_refused_at_a_saturated_fork() {
        // 5% walked against 95%, the rescue case, but the majority is clipped.
        let g = diamond(1_250.0, crate::types::KA_CEILING);
        let s = subgraph_of(&g);
        let c = Config { rescue_minority: true, ..maj(true) };
        let (fixed, report) =
            super::prefer_majority(&c, &g, &s, &walked(), &q("ACGTAG")).unwrap();
        assert!(fixed.is_none(), "the majority test rests on a clipped total");
        assert_eq!(report.left.saturated, 1, "and the refusal is reported: {:?}", report.left);

        // One below the ceiling, same ratio, and the rescue goes ahead.
        let g = diamond(1_250.0, crate::types::KA_CEILING - 1.0);
        let s = subgraph_of(&g);
        let (fixed, report) =
            super::prefer_majority(&c, &g, &s, &walked(), &q("ACGTAG")).unwrap();
        assert!(fixed.is_some(), "24,999 is a real number; the rescue is sound");
        assert_eq!(report.minority_rescues, 1);
        assert_eq!(report.left.saturated, 0);
    }

    /// An ORDINARY two-allele correction is not refused at a saturated fork, because it needs
    /// the alternative to beat the walked branch, which clipping cannot manufacture: two
    /// clipped branches compare equal and nothing moves.
    #[test]
    fn an_ordinary_correction_is_unaffected_by_saturation() {
        // Both branches well over --min-vaf, the alternative clipped and clearly larger.
        let g = diamond(12_000.0, crate::types::KA_CEILING);
        let s = subgraph_of(&g);
        let (fixed, report) =
            super::prefer_majority(&maj(true), &g, &s, &walked(), &q("ACGTAG")).unwrap();
        assert!(fixed.is_some(), "not a minority site, so the guard does not apply");
        assert_eq!(report.minority_rescues, 0, "an ordinary correction, not a rescue");
        assert_eq!(report.left.saturated, 0);
    }

    /// The bar that replaces the old exemption. A plurality is not enough: with the route on
    /// 5% and two alternatives at 47.5% each, there is no majority to move onto, and moving
    /// on a plurality is the greedy rewrite this must not become.
    #[test]
    fn a_walked_minority_is_left_alone_without_an_outright_majority() {
        let fa = ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+  L:+:3:+\nACGT\n\
                  >T_1 ka:f:5.0    L:-:0:-  L:+:4:+\nCGTA\n\
                  >T_2 ka:f:47.0   L:-:0:-  L:+:4:+\nCGTC\n\
                  >T_3 ka:f:47.0   L:-:0:-  L:+:4:+\nCGTG\n\
                  >T_4 ka:f:100.0  L:-:1:-  L:-:2:-  L:-:3:-\nGTAG\n";
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX)
            .unwrap();
        let s = subgraph_of(&g);
        let walk = GafPath {
            steps: vec![
                PathStep { orient: Orient::Forward, node: 0 },
                PathStep { orient: Orient::Forward, node: 1 },
                PathStep { orient: Orient::Forward, node: 4 },
            ],
            ..walked()
        };
        let c = Config { rescue_minority: true, ..maj(true) };
        let (fixed, report) =
            super::prefer_majority(&c, &g, &s, &walk, &q("ACGTAG")).unwrap();
        assert!(fixed.is_none(), "47.5% is a plurality, not a majority");
        assert_eq!(report.swaps, 0);
        // Left, and said so, with the reason.
        assert_eq!(report.left.not_outright_majority, 1, "{:?}", report.left);
        assert_eq!(report.left.total(), 1);
    }



    /// `prefer_majority` and `find_alts` must agree on what an allele is. An in-frame indel
    /// is correctable, so a majority allele that inserts a codon replaces the minority one
    /// in the consensus rather than being emitted with a VAF above 0.5.
    #[test]
    fn an_in_frame_indel_majority_is_swapped_in() {
        let fa = ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n\
                  >T_1 ka:f:40.0   L:-:0:-  L:+:3:+\nCGTA\n\
                  >T_2 ka:f:60.0   L:-:0:-  L:+:3:+\nCGTCAAA\n\
                  >T_3 ka:f:100.0  L:-:1:-  L:-:2:-\nGTAG\n";
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX).unwrap();
        let s = subgraph_of(&g);
        let fixed = prefer_majority(&maj(true), &g, &s, &walked(), &q("ACGTAG"))
            .unwrap()
            .expect("the 60% allele inserts a codon and should still win");
        assert!(fixed.steps.iter().any(|p| p.node == 2));
    }

    /// A frameshifting majority is not swapped in, for the same reason it is not emitted:
    /// it would reframe the consensus that feeds the codon alignment.
    #[test]
    fn a_frameshifting_majority_is_not_swapped_in() {
        let fa = ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n\
                  >T_1 ka:f:40.0   L:-:0:-  L:+:3:+\nCGTA\n\
                  >T_2 ka:f:60.0   L:-:0:-  L:+:3:+\nCGTCA\n\
                  >T_3 ka:f:100.0  L:-:1:-  L:-:2:-\nGTAG\n";
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX).unwrap();
        let s = subgraph_of(&g);
        assert!(prefer_majority(&maj(true), &g, &s, &walked(), &q("ACGTAG")).unwrap().is_none());
        let c = Config { allow_indel_alts: false, ..maj(true) };
        assert!(prefer_majority(&c, &g, &s, &walked(), &q("ACGTAG")).unwrap().is_none(),
                "and not with indel support switched off either");
        // `--allow-frameshift-alts` admits a frameshifting allele to the side file; it must
        // not reframe `_primary`, so the rule here is deliberately not shared with
        // `find_alts`.
        let c = Config { allow_frameshift_alts: true, ..maj(true) };
        assert!(prefer_majority(&c, &g, &s, &walked(), &q("ACGTAG")).unwrap().is_none(),
                "--allow-frameshift-alts must not let the consensus be reframed");
        // And the same on a flanked window, now that such a window reaches the rule at all.
        assert!(prefer_majority(&c, &g, &s, &walked_window(1, 5), &q("ACGTAG")).unwrap().is_none(),
                "nor on a flanked window");
    }

    #[test]
    fn emits_nothing_when_switched_off() {
        assert!(Config::default().emit_alts, "on by default; the records go to their own file");
        let g = diamond(60.0, 40.0);
        let s = subgraph_of(&g);
        let alts = find_alts(&cfg(false), &g, &s, &walked(), &q("ACGTAG")).unwrap();
        assert!(alts.is_empty());
    }

    /// 60/40: both alleles clear the gates, so the minor one is a real second haplotype.
    #[test]
    fn recovers_the_second_allele_at_a_biallelic_site() {
        let g = diamond(60.0, 40.0);
        let s = subgraph_of(&g);
        let primary = stitch::stitch(&cfg(true), &g, &walked(), &s, &q("ACGTAG")).unwrap().unwrap();
        let alts = find_alts(&cfg(true), &g, &s, &walked(), &q("ACGTAG")).unwrap();

        assert_eq!(alts.len(), 1, "one site, one record");
        let a = &alts[0];
        assert_eq!(a.haplotype.id, "ACC_alt1");
        assert_eq!(a.fork, 0);
        assert_eq!(a.alt_node, 2);
        assert_eq!(a.detour_nodes, 1);
        assert!((a.alt_vaf - 0.4).abs() < 1e-9, "vaf was {}", a.alt_vaf);
        assert!((a.taken_vaf - 0.6).abs() < 1e-9);
        assert_ne!(a.haplotype.seq, primary.seq, "an alt that equals the primary is not a variant");
        assert_eq!(a.haplotype.seq.len(), primary.seq.len(), "substitution, not indel");
    }

    /// THE filter. At 99.7/0.3 the walk took a rare allele, which is the reference-pull
    /// defect, not two populations. Recovering the 99.7% branch here would be fixing path
    /// selection while calling the result a quasispecies, so nothing is emitted.
    #[test]
    fn ignores_a_site_where_the_walk_took_a_rare_allele() {
        let g = diamond(3.0, 12000.0);
        let s = subgraph_of(&g);
        let alts = find_alts(&cfg(true), &g, &s, &walked(), &q("ACGTAG")).unwrap();
        assert!(
            alts.is_empty(),
            "the walked branch is 0.02% here; that is a path-selection bug, not a variant"
        );
    }

    /// The mirror case: the untaken branch is below the floors, so it is sequencing error.
    #[test]
    fn ignores_an_untaken_branch_below_the_floors() {
        let g = diamond(1000.0, 4.0); // 4x and 0.4% both fail
        let s = subgraph_of(&g);
        assert!(find_alts(&cfg(true), &g, &s, &walked(), &q("ACGTAG")).unwrap().is_empty());
    }

    /// An indel allele whose length change is a multiple of three keeps the reading frame,
    /// so it is emitted and the window is shifted to match. Branch 2 is 3 bp longer.
    #[test]
    fn accepts_an_in_frame_indel_allele() {
        let fa = ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n\
                  >T_1 ka:f:60.0   L:-:0:-  L:+:3:+\nCGTA\n\
                  >T_2 ka:f:40.0   L:-:0:-  L:+:3:+\nCGTCAAA\n\
                  >T_3 ka:f:100.0  L:-:1:-  L:-:2:-\nGTAG\n";
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX).unwrap();
        let s = subgraph_of(&g);
        let query = q("ACGTAG");
        let primary = stitch::stitch(&cfg(true), &g, &walked(), &s, &query).unwrap().unwrap();
        let alts = find_alts(&cfg(true), &g, &s, &walked(), &query).unwrap();
        assert_eq!(alts.len(), 1, "an in-frame indel allele is recoverable");
        assert_eq!(alts[0].length_delta, 3);
        assert_eq!(
            alts[0].haplotype.seq.len(),
            primary.seq.len() + 3,
            "the window must shift by exactly the delta, or the record is reframed"
        );
    }

    /// The same allele with indel support switched off, which is the behaviour that held
    /// before `shifted_window` existed.
    #[test]
    fn refuses_length_changes_when_indel_support_is_off() {
        let fa = ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n\
                  >T_1 ka:f:60.0   L:-:0:-  L:+:3:+\nCGTA\n\
                  >T_2 ka:f:40.0   L:-:0:-  L:+:3:+\nCGTCAAA\n\
                  >T_3 ka:f:100.0  L:-:1:-  L:-:2:-\nGTAG\n";
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX).unwrap();
        let s = subgraph_of(&g);
        let c = Config { allow_indel_alts: false, ..cfg(true) };
        assert!(find_alts(&c, &g, &s, &walked(), &q("ACGTAG")).unwrap().is_empty());
    }

    /// A single-base indel reframes everything downstream. 190 of the 203 length-changing
    /// calls measured across the panel are this shape, with a median magnitude of one base,
    /// and they cannot be told apart from a homopolymer artifact, so they are opt-in.
    #[test]
    fn a_frameshifting_allele_is_refused_unless_asked_for() {
        let fa = ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n\
                  >T_1 ka:f:60.0   L:-:0:-  L:+:3:+\nCGTA\n\
                  >T_2 ka:f:40.0   L:-:0:-  L:+:3:+\nCGTCA\n\
                  >T_3 ka:f:100.0  L:-:1:-  L:-:2:-\nGTAG\n";
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX).unwrap();
        let s = subgraph_of(&g);
        assert!(
            find_alts(&cfg(true), &g, &s, &walked(), &q("ACGTAG")).unwrap().is_empty(),
            "+1 nt shifts the frame and is off by default"
        );
        let c = Config { allow_frameshift_alts: true, ..cfg(true) };
        let got = find_alts(&c, &g, &s, &walked(), &q("ACGTAG")).unwrap();
        assert_eq!(got.len(), 1, "and recoverable when explicitly enabled");
        assert_eq!(got[0].length_delta, 1);
    }

    /// A large delta is the route being sent through a different stretch of graph, not an
    /// allele. The panel's tail sits at +100, +142 and -282.
    #[test]
    fn a_large_length_change_is_a_rerouting_not_an_allele() {
        let long: String = std::iter::repeat_n("ACGT", 20).collect(); // 80 nt, delta +76
        let fa = format!(
            ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n\
             >T_1 ka:f:60.0   L:-:0:-  L:+:3:+\nCGTA\n\
             >T_2 ka:f:40.0   L:-:0:-  L:+:3:+\n{long}\n\
             >T_3 ka:f:100.0  L:-:1:-  L:-:2:-\nGTAG\n"
        );
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.into_bytes())), u64::MAX).unwrap();
        let s = subgraph_of(&g);
        let c = Config { allow_frameshift_alts: true, ..cfg(true) };
        assert!(
            find_alts(&c, &g, &s, &walked(), &q("ACGTAG")).unwrap().is_empty(),
            "past max_alt_indel it is refused even with every other gate open"
        );
    }

    /// The window arithmetic itself.
    #[test]
    fn the_window_shifts_only_when_the_splice_sits_inside_it() {
        // splice at 100, replacing 50 bases, inside a window of [26, 2764)
        assert_eq!(shifted_window(26, 2764, 100, 50, 3), Some((26, 2767)));
        assert_eq!(shifted_window(26, 2764, 100, 50, -3), Some((26, 2761)));
        // a substitution never moves the window
        assert_eq!(shifted_window(26, 2764, 100, 50, 0), Some((26, 2764)));
        // the splice begins before the window starts: the prefix is not intact
        assert_eq!(shifted_window(26, 2764, 10, 50, 3), None);
        // the replaced stretch runs past the window end
        assert_eq!(shifted_window(26, 2764, 2700, 100, 3), None);
        // a deletion that would collapse the window
        assert_eq!(shifted_window(26, 100, 30, 20, -200), None);
    }

    /// Off by default, so the floor is exactly `min_abund` and nothing changes.
    #[test]
    fn depth_relative_floor_is_off_by_default() {
        let c = Config::default();
        assert_eq!(c.min_abund_frac, 0.0);
        assert_eq!(abundance_floor(&c, 15000.0), c.min_abund, "must not bite unless asked");
    }

    /// The case this exists for. At 15,000x, five reads is 0.03% and admits error branches.
    #[test]
    fn scales_the_floor_with_sequencing_depth() {
        let c = Config { min_abund_frac: 0.01, ..Config::default() };
        assert_eq!(abundance_floor(&c, 15000.0), 150.0, "1% of a deep sample");
        assert_eq!(abundance_floor(&c, 180.0), 5.0, "1% of 180x is below min_abund, so it holds");
    }

    /// It may only tighten. A shallow accession must never end up with a floor under the
    /// one the rest of the pipeline already assumes.
    #[test]
    fn never_loosens_the_existing_floor() {
        let c = Config { min_abund_frac: 0.5, min_abund: 5.0, ..Config::default() };
        for depth in [0.0, 1.0, 9.0, 10.0, 1000.0] {
            assert!(abundance_floor(&c, depth) >= 5.0, "floor dipped at depth {depth}");
        }
    }

    /// Depth is the median over the walked path, not the mean, so one repeat-collapsed
    /// unitig cannot drag it.
    #[test]
    fn depth_is_the_median_of_the_walked_path() {
        let g = diamond(60.0, 40.0);
        let s = subgraph_of(&g);
        let walk = vec![1u32, 3, 7]; // nodes 0, 1, 3 forward: ka:f 100, 60, 100
        assert_eq!(path_depth(&g, &s, &walk), 100.0);
    }

    /// End to end: the same 60/40 site is recovered at the default floor and suppressed once
    /// the floor is raised above the alt, without touching min_vaf.
    #[test]
    fn a_marginal_allele_is_dropped_once_the_floor_scales() {
        let g = diamond(60.0, 40.0);
        let s = subgraph_of(&g);
        let query = q("ACGTAG");
        assert_eq!(find_alts(&cfg(true), &g, &s, &walked(), &query).unwrap().len(), 1);

        // path depth is 100x; 0.5 puts the floor at 50, above the 40x alt
        let c = Config { min_abund_frac: 0.5, ..cfg(true) };
        assert!(
            find_alts(&c, &g, &s, &walked(), &query).unwrap().is_empty(),
            "40x is 40% of local coverage but only 0.4x the gene depth"
        );
    }

    #[test]
    fn respects_the_cap() {
        let g = diamond(60.0, 40.0);
        let s = subgraph_of(&g);
        let c = Config { max_alts: 0, ..cfg(true) };
        assert!(find_alts(&c, &g, &s, &walked(), &q("ACGTAG")).unwrap().is_empty());
    }

    /// A well-covered exit that never returns to the route is not an allele of this gene and
    /// must not sit in the denominator.
    ///
    /// The shape is a repeat, a paralog, or a contaminant sharing the fork's 30-mer. Counting
    /// it deflates every frequency at the site: here the allele is 20 against the route's 80,
    /// which is 0.2, but with a 400x dead-end exit in the total it reads 0.04 and clears no
    /// gate, so there is no alt, no VCF row and no snvs.tsv row. `phase` already excluded
    /// such branches from its own total, so `bubbles.json` and the VCF disagreed about the
    /// same site.
    #[test]
    fn an_exit_that_never_rejoins_is_left_out_of_the_fork_total() {
        // 0 forks to 1 (walked, 80x), 2 (allele, 20x) and 3 (a 400x dead end). 1 and 2
        // rejoin at 4; 3 goes nowhere.
        let fa = ">T_0 ka:f:80.0   L:+:1:+  L:+:2:+  L:+:3:+\nACGT\n\
                  >T_1 ka:f:80.0   L:-:0:-  L:+:4:+\nCGTA\n\
                  >T_2 ka:f:20.0   L:-:0:-  L:+:4:+\nCGTC\n\
                  >T_3 ka:f:400.0  L:-:0:-\nCGTG\n\
                  >T_4 ka:f:80.0   L:-:1:-  L:-:2:-\nGTAG\n";
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX)
            .unwrap();
        let s = subgraph_of(&g);
        let walk = GafPath {
            steps: vec![
                PathStep { orient: Orient::Forward, node: 0 },
                PathStep { orient: Orient::Forward, node: 1 },
                PathStep { orient: Orient::Forward, node: 4 },
            ],
            ..walked()
        };
        let (alts, _) = super::find_alts(&cfg(true), &g, &s, &walk, &q("ACGTAG")).unwrap();
        assert_eq!(alts.len(), 1, "the 20x allele must be recovered");
        // 20 / (80 + 20), not 20 / (80 + 20 + 400) = 0.04.
        assert!(
            (alts[0].alt_vaf - 0.2).abs() < 1e-9,
            "VAF must exclude the dead-end exit, got {}",
            alts[0].alt_vaf
        );
        assert!((alts[0].taken_vaf - 0.8).abs() < 1e-9, "{}", alts[0].taken_vaf);
    }

    /// The other half: with the route on the minority, a dead-end exit in the denominator
    /// holds the majority correction back and leaves the 20% base on `_primary`.
    #[test]
    fn a_dead_end_exit_no_longer_blocks_the_majority_correction() {
        let fa = ">T_0 ka:f:80.0   L:+:1:+  L:+:2:+  L:+:3:+\nACGT\n\
                  >T_1 ka:f:20.0   L:-:0:-  L:+:4:+\nCGTA\n\
                  >T_2 ka:f:80.0   L:-:0:-  L:+:4:+\nCGTC\n\
                  >T_3 ka:f:400.0  L:-:0:-\nCGTG\n\
                  >T_4 ka:f:80.0   L:-:1:-  L:-:2:-\nGTAG\n";
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX)
            .unwrap();
        let s = subgraph_of(&g);
        let walk = GafPath {
            steps: vec![
                PathStep { orient: Orient::Forward, node: 0 },
                PathStep { orient: Orient::Forward, node: 1 },
                PathStep { orient: Orient::Forward, node: 4 },
            ],
            ..walked()
        };
        let (fixed, report) =
            super::prefer_majority(&maj(true), &g, &s, &walk, &q("ACGTAG")).unwrap();
        let fixed = fixed.expect("80x against 20x is a majority once the dead end is excluded");
        assert!(fixed.steps.iter().any(|p| p.node == 2));
        assert_eq!(report.swaps, 1);
    }

    /// What the cap leaves out is counted and reported, not dropped in silence.
    ///
    /// Breaking at the cap would give output for the 5' end of the gene with the rest looking
    /// invariant, and the stderr counters would stop with the loop and agree with it.
    #[test]
    fn alleles_past_the_cap_are_counted_and_the_scan_continues() {
        let g = diamond(60.0, 40.0);
        let s = subgraph_of(&g);
        let c = Config { max_alts: 0, ..cfg(true) };
        let (alts, stats) = super::find_alts(&c, &g, &s, &walked(), &q("ACGTAG")).unwrap();
        assert!(alts.is_empty(), "the cap is 0, so nothing is emitted");
        assert_eq!(stats.beyond_cap, 1, "but the credible allele is counted");
        // The site was still examined, which is what the old `break` prevented.
        assert_eq!(stats.forks, 1);
        assert_eq!(stats.biallelic, 1);
    }

    /// A path that leaves the subgraph cannot be spliced against it. `stitch` returns None
    /// rather than erroring in the same situation, and this must not diverge.
    #[test]
    fn a_path_outside_the_subgraph_yields_nothing() {
        let g = diamond(60.0, 40.0);
        let s = subgraph_of(&g);
        let mut p = walked();
        p.steps[1].node = 99;
        assert!(find_alts(&cfg(true), &g, &s, &p, &q("ACGTAG")).unwrap().is_empty());
    }
    /// Regression set: the alignment window sits INSIDE the path.
    ///
    /// Every other test in this module uses `walked()`, whose window spans the whole path.
    /// Real GAF rows almost never look like that; the one the repository records is `7..1281`
    /// of 1288. Comparing a restitched sequence against the whole stitched path rather than
    /// against the window silently drops every correction on a flanked alignment, and only a
    /// flanked fixture catches it.
    ///
    /// The route takes the 40% branch; after correction the window must carry the 60% one.
    #[test]
    fn the_majority_is_swapped_in_when_the_window_is_inside_the_path() {
        let g = diamond(40.0, 60.0);
        let s = subgraph_of(&g);
        let query = q("ACGTAG");
        let c = maj(true);
        for (start, end, before_seq, after_seq) in [
            // Flanked at both ends, the shape of a real row.
            (1usize, 5usize, "CGTA", "CGTC"),
            // Truncated only at the end: still unequal to the path, so still affected.
            (0, 5, "ACGTA", "ACGTC"),
        ] {
            let w = walked_window(start, end);
            let before = stitch::stitch(&c, &g, &w, &s, &query).unwrap().unwrap();
            assert_eq!(
                before.seq,
                before_seq.as_bytes(),
                "window [{start},{end}) of the 40% route"
            );
            let fixed = prefer_majority(&c, &g, &s, &w, &query)
                .unwrap()
                .unwrap_or_else(|| panic!("window [{start},{end}) should be corrected"));
            let after = stitch::stitch(&c, &g, &fixed, &s, &query).unwrap().unwrap();
            assert_eq!(
                after.seq,
                after_seq.as_bytes(),
                "window [{start},{end}) should publish the 60% allele"
            );
            assert_eq!(after.seq.len(), before.seq.len(), "substitution, not indel");
            assert!(fixed.steps.iter().any(|p| p.node == 2), "must traverse the 60% branch");
            assert!(!fixed.steps.iter().any(|p| p.node == 1), "and not the 40% one");
        }
    }

    /// The same, for a length-changing allele: the window must move by exactly the delta.
    /// This is the case the old whole-path comparison could not express at all, because the
    /// path length and the window length changed by the same amount but from different bases.
    #[test]
    fn an_in_frame_indel_majority_is_swapped_in_with_a_flanked_window() {
        let fa = ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n\
                  >T_1 ka:f:40.0   L:-:0:-  L:+:3:+\nCGTA\n\
                  >T_2 ka:f:60.0   L:-:0:-  L:+:3:+\nCGTCAAA\n\
                  >T_3 ka:f:100.0  L:-:1:-  L:-:2:-\nGTAG\n";
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX)
            .unwrap();
        let s = subgraph_of(&g);
        let c = maj(true);
        let fixed = prefer_majority(&c, &g, &s, &walked_window(1, 5), &q("ACGTAG"))
            .unwrap()
            .expect("the 60% allele inserts a codon and should still win on a flanked window");
        assert!(fixed.steps.iter().any(|p| p.node == 2));
        let after = stitch::stitch(&c, &g, &fixed, &s, &q("ACGTAG")).unwrap().unwrap();
        assert_eq!(after.seq.len(), 4 + 3, "the window grew by the inserted codon");
    }

    /// Pins `find_alts` on the same shape of input as `prefer_majority`: a regression that
    /// made it compare whole paths rather than windows would otherwise be caught by nothing.
    #[test]
    fn recovers_the_second_allele_when_the_window_is_inside_the_path() {
        let g = diamond(60.0, 40.0);
        let s = subgraph_of(&g);
        let w = walked_window(1, 5);
        let primary = stitch::stitch(&cfg(true), &g, &w, &s, &q("ACGTAG")).unwrap().unwrap();
        let alts = find_alts(&cfg(true), &g, &s, &w, &q("ACGTAG")).unwrap();
        assert_eq!(alts.len(), 1, "one site, one record");
        assert_eq!(alts[0].alt_node, 2);
        assert_eq!(alts[0].haplotype.seq, b"CGTC", "the alt over the same window");
        assert_ne!(alts[0].haplotype.seq, primary.seq);
        assert_eq!(alts[0].haplotype.seq.len(), primary.seq.len());
    }
    /// A triallelic fork where the most abundant alternative is refused: a +1 nt frameshift
    /// at 45%, which never reaches `_primary`. The next one down, a clean substitution at
    /// 30%, still beats the walked 25% branch and must be taken. Considering only the top
    /// alternative would abandon the fork and leave the 25% allele published.
    #[test]
    fn a_refused_top_alternative_falls_through_to_the_next() {
        let fa = ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+  L:+:4:+\nACGT\n\
                  >T_1 ka:f:25.0   L:-:0:-  L:+:3:+\nCGTA\n\
                  >T_2 ka:f:45.0   L:-:0:-  L:+:3:+\nCGTCA\n\
                  >T_3 ka:f:100.0  L:-:1:-  L:-:2:-  L:-:4:-\nGTAG\n\
                  >T_4 ka:f:30.0   L:-:0:-  L:+:3:+\nCGTC\n";
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX)
            .unwrap();
        let s = subgraph_of(&g);
        let c = maj(true);
        let fixed = prefer_majority(&c, &g, &s, &walked(), &q("ACGTAG"))
            .unwrap()
            .expect("the 30% substitution should be taken after the 45% frameshift is refused");
        assert!(fixed.steps.iter().any(|p| p.node == 4), "must take the 30% substitution");
        assert!(!fixed.steps.iter().any(|p| p.node == 2), "never the frameshifting 45% branch");
        assert!(!fixed.steps.iter().any(|p| p.node == 1), "and not the walked 25% one");
        let after = stitch::stitch(&c, &g, &fixed, &s, &q("ACGTAG")).unwrap().unwrap();
        assert_eq!(after.seq, b"ACGTCG", "in frame, one base changed");

        // The same on a flanked window, which is where the window comparison matters.
        let w = walked_window(1, 5);
        let fixed = prefer_majority(&c, &g, &s, &w, &q("ACGTAG"))
            .unwrap()
            .expect("the fall-through must work on a flanked window too");
        assert!(fixed.steps.iter().any(|p| p.node == 4));
        let after = stitch::stitch(&c, &g, &fixed, &s, &q("ACGTAG")).unwrap().unwrap();
        assert_eq!(after.seq, b"CGTC", "the 30% allele over the flanked window");
    }
}
