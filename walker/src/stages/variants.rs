//! Step 8: recover the intra-accession haplotypes the primary walk leaves behind.
//!
//! `stitch` emits one sequence per accession. Where the sample genuinely carries two alleles
//! at a site, the second one is discarded, and nothing downstream ever sees it. This stage
//! puts it back, as one extra record per site.
//!
//! IN   UnitigGraph, Subgraph, GafPath, Query
//! OUT  Vec<AltHaplotype>, empty unless `cfg.emit_alts`
//!
//! # What counts as a site worth recovering
//!
//! A fork the walked path passes through, where **both** the branch taken and the branch not
//! taken clear `min_abund` and `min_vaf`. That "both" is the whole filter, and it is what
//! separates real quasispecies from the far more common failure next door.
//!
//! Measured on `SRR8209082` (HIV-1 *pol*), over the 58 on-path forks that have a credible
//! untaken branch:
//!
//! ```text
//! walked branch is itself under 10%   54   the aligner took a rare allele; not a variant
//! both branches between 10% and 90%    4   genuinely two populations
//! ```
//!
//! The 54 are the reference-pull defect: GraphAligner scores paths by residue matches to the
//! reference, so where the sample's dominant allele differs from the reference and an error
//! branch matches it, the walk takes the error branch. Emitting those as "variants" would
//! publish sequencing error as biology. They are a bug to fix in path selection, not signal
//! to recover here, and `work/issue_intra_accession_variance_lost.md` covers them.
//!
//! # Why one variant per record
//!
//! With *n* independent sites there are up to 2^n haplotypes, and a de Bruijn graph cannot
//! say which alternatives sit on the same molecule. Emitting combinations would invent
//! sequences that may not exist in any organism, which is the chimera problem the walker
//! exists to avoid. Every record here is the primary path with **exactly one** substitution,
//! so it is a real connected path through the graph and asserts only the minimal thing.
//!
//! In practice the cap is not the binding constraint: the measured site count is 4, 4 and 0
//! across the accessions looked at, not 50.
//!
//! # Why substitutions only
//!
//! An accepted detour must leave the stitched sequence exactly as long as the primary. That
//! keeps `path_start`/`path_end`, which are offsets into the stitched sequence, valid without
//! adjustment, and it rejects the dangerous shape: a detour that rejoins far downstream is a
//! path *shortcut* that silently deletes gene. Measured while prototyping, a coverage-greedy
//! re-walk without this constraint cut a 1,279-node path to 334 and the sequence from 2,738
//! to 653 nt. Every genuine allele measured so far is an equal-length substitution, so this
//! costs nothing real.

use crate::config::Config;
use crate::error::Result;
use crate::stages::{parse, stitch};
use crate::types::{
    AltHaplotype, GafPath, Orient, PathStep, Query, Strand, Subgraph, UnitigGraph, UnitigId,
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
    local: HashMap<UnitigId, usize>,
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
    Some(View { out, walk, at, local })
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
/// only where **both** branches clear `min_abund` and `min_vaf`, which is a handful of sites
/// per accession rather than every fork, and each one must still rejoin the route and leave
/// the sequence exactly as long.
///
/// Returns the corrected route, or `None` if nothing needed changing.
pub fn prefer_majority(
    cfg: &Config,
    graph: &UnitigGraph,
    subgraph: &Subgraph,
    path: &GafPath,
    query: &Query,
) -> Result<Option<GafPath>> {
    if !cfg.prefer_majority {
        return Ok(None);
    }
    let mut baseline_len = stitched_len(cfg, &path.steps, subgraph);
    let mut current = GafPath { steps: path.steps.clone(), ..copy_coords(path) };
    let mut swaps = 0usize;

    // Each swap fixes one site and makes that site's taken branch the majority, so it cannot
    // be chosen again. The cap is a backstop against a shape that defeats that reasoning.
    for _ in 0..cfg.max_alts.max(1) * 2 {
        let Some(v) = view(graph, subgraph, &current) else { break };
        let abund = |o: Oriented| graph.abundance[subgraph.nodes[(o as usize) / 2] as usize];
        let depth = path_depth(graph, subgraph, &v.walk);
        let floor = abundance_floor(cfg, depth);

        let mut applied = false;
        for i in 0..v.walk.len() - 1 {
            let (fork, taken) = (v.walk[i], v.walk[i + 1]);
            let branches = &v.out[fork as usize];
            if branches.len() < 2 || !branches.contains(&taken) {
                continue;
            }
            let total: f64 = branches.iter().map(|&b| abund(b) as f64).sum();
            if total <= 0.0 {
                continue;
            }
            let taken_ab = abund(taken);
            if taken_ab < floor || (taken_ab as f64 / total) < cfg.min_vaf {
                continue; // not a two-allele site; this is the error-branch case
            }
            // The best-supported credible alternative, and only if it beats what we took.
            let mut best: Option<(Oriented, f32)> = None;
            for &alt in branches {
                if alt == taken {
                    continue;
                }
                let ab = abund(alt);
                if ab < floor || (ab as f64 / total) < cfg.min_vaf || ab <= taken_ab {
                    continue;
                }
                if best.map_or(true, |(_, b)| ab > b) {
                    best = Some((alt, ab));
                }
            }
            let Some((alt, _)) = best else { continue };
            let Some((rejoin_at, detour)) = rejoin(&v.out, &v.at, alt, i, cfg.alt_max_detour)
            else {
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

            // The same length rules as `find_alts`, so the two stages agree on what counts
            // as an allele. Without this, an indel site could not be corrected and its
            // alternative would be emitted with a VAF above 0.5, meaning the primary had
            // kept the rarer base there. Measured before the rules were shared: 1 of 48
            // alternatives on `ebola_gp`, every one of them an indel.
            let spliced_len = stitched_len(cfg, &steps, subgraph);
            let delta = spliced_len as i64 - baseline_len as i64;
            if delta != 0 {
                if !cfg.allow_indel_alts
                    || delta.unsigned_abs() as usize > cfg.max_alt_indel
                    || (delta % 3 != 0 && !cfg.allow_frameshift_alts)
                {
                    continue;
                }
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
                continue;
            };
            let candidate = GafPath {
                steps,
                path_start: new_start,
                path_end: new_end,
                ..copy_coords(path)
            };
            if stitch::stitch(cfg, &candidate, subgraph, query)?.is_none() {
                continue; // would fail the coverage gate; leave the site alone
            }
            // A length-changing swap moves the baseline every later swap is measured
            // against, so it has to be updated rather than left at the original route's.
            baseline_len = spliced_len;
            current = candidate;
            swaps += 1;
            applied = true;
            break;
        }
        if !applied {
            break;
        }
    }

    if swaps == 0 {
        return Ok(None);
    }
    eprintln!(
        "[{}] consensus: {} site(s) where the route carried the minority allele, corrected",
        subgraph.accession, swaps
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
    if swaps > 12 {
        eprintln!(
            "[{}] WARNING: {} consensus corrections is far more than the 1-11 seen across \
             the rest of the panel. The route may have been rerouted through a divergent \
             segment rather than having single alleles swapped; treat this accession as \
             provisional and compare it against the unswapped result.",
            subgraph.accession, swaps
        );
    }
    Ok(Some(current))
}

/// Recover one alternative haplotype per genuinely biallelic site on the walked path.
pub fn find_alts(
    cfg: &Config,
    graph: &UnitigGraph,
    subgraph: &Subgraph,
    path: &GafPath,
    query: &Query,
) -> Result<Vec<AltHaplotype>> {
    // The sites feed both alt_haplotypes.fasta and the VCF, so compute them when either is
    // wanted. Skipping the work only when neither output will use it.
    if (!cfg.emit_alts && !cfg.emit_vcf) || path.steps.len() < 2 || subgraph.nodes.is_empty() {
        return Ok(Vec::new());
    }

    let n = subgraph.nodes.len();
    let mut local: HashMap<UnitigId, usize> = HashMap::with_capacity(n * 2);
    for (i, &id) in subgraph.nodes.iter().enumerate() {
        local.insert(id, i);
    }

    // Same oriented adjacency as `phase`: an edge in u's header is the traversal
    // (u, s1) -> (t, s2), and the reverse arrives independently from t's own header.
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

    // The walked path as oriented ids, and where each one sits.
    let mut walk: Vec<Oriented> = Vec::with_capacity(path.steps.len());
    for step in &path.steps {
        let Some(&i) = local.get(&step.node) else {
            return Ok(Vec::new()); // path leaves the subgraph; stitch would bail too
        };
        let fwd = matches!(step.orient, Orient::Forward);
        walk.push((i * 2 + usize::from(fwd)) as Oriented);
    }
    let mut at: HashMap<Oriented, usize> = HashMap::with_capacity(walk.len());
    for (i, &o) in walk.iter().enumerate() {
        // First occurrence wins; a repeated node makes "where does it rejoin" ambiguous and
        // the rejoin search below only ever moves forward.
        at.entry(o).or_insert(i);
    }

    let Some(primary) = stitch::stitch(cfg, path, subgraph, query)? else {
        return Ok(Vec::new());
    };
    let baseline_len = stitched_len(cfg, &path.steps, subgraph);

    let abund = |o: Oriented| graph.abundance[subgraph.nodes[(o as usize) / 2] as usize];
    let depth = path_depth(graph, subgraph, &walk);
    let floor = abundance_floor(cfg, depth);

    // Counters, so a zero result can be told apart from a filtered one. Calibration needs
    // to know which gate is doing the work, not just the total.
    let (mut forks, mut biallelic, mut no_rejoin, mut not_substitution, mut identical) =
        (0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut rejected_large, mut rejected_frameshift, mut rejected_window, mut indels) =
        (0usize, 0usize, 0usize, 0usize);
    let mut deltas: Vec<i64> = Vec::new();
    let mut alts: Vec<AltHaplotype> = Vec::new();
    for i in 0..walk.len() - 1 {
        if alts.len() >= cfg.max_alts {
            break;
        }
        let fork = walk[i];
        let taken = walk[i + 1];
        let branches = &out[fork as usize];
        if branches.len() < 2 || !branches.contains(&taken) {
            continue;
        }
        forks += 1;
        let total: f64 = branches.iter().map(|&b| abund(b) as f64).sum();
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

        for &alt in branches {
            if alt == taken || alts.len() >= cfg.max_alts {
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

            // Substitutions only. See the module note: this keeps path_start/path_end valid
            // and rejects shortcuts that delete gene.
            let spliced_len = stitched_len(cfg, &steps, subgraph);
            let delta = spliced_len as i64 - baseline_len as i64;

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
            let Some(mut hap) = stitch::stitch(cfg, &spliced, subgraph, query)? else {
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
    Ok(alts)
}

/// The gene's own sequencing depth: median `ka:f` over the nodes the walk uses.
///
/// The walked path is the gene, so this measures the quantity the floor is relative to.
/// Median rather than mean because a single repeat-collapsed unitig at 25,000x should not
/// drag the estimate for an otherwise ordinary accession.
fn path_depth(graph: &UnitigGraph, subgraph: &Subgraph, walk: &[Oriented]) -> f32 {
    let mut ab: Vec<f32> = walk
        .iter()
        .map(|&o| graph.abundance[subgraph.nodes[(o as usize) / 2] as usize])
        .collect();
    if ab.is_empty() {
        return 0.0;
    }
    ab.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    ab[ab.len() / 2]
}

/// The abundance a branch must clear to count as an allele.
///
/// Never below `min_abund`: this can only tighten the existing gate, never loosen it.
fn abundance_floor(cfg: &Config, depth: f32) -> f32 {
    if cfg.min_abund_frac > 0.0 {
        (cfg.min_abund_frac as f32 * depth).max(cfg.min_abund)
    } else {
        cfg.min_abund
    }
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
    let overlap = cfg.graph_k.saturating_sub(1);
    let mut total = 0usize;
    for (i, step) in steps.iter().enumerate() {
        let Ok(idx) = subgraph.nodes.binary_search(&step.node) else {
            return usize::MAX;
        };
        let len = subgraph.seqs[idx].len();
        total += if i == 0 {
            len
        } else {
            len.saturating_sub(overlap)
        };
    }
    total
}

/// Shortest way from `alt` back onto the walked path, strictly after the fork.
///
/// Returns the path index it rejoins at and the detour, excluding the rejoin node itself
/// because the path supplies that. Breadth-first, so the detour found is the shortest, and
/// bounded so a tangled neighbourhood cannot make this expensive.
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
        parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.into_bytes()))).unwrap()
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
        }
    }
    fn walked() -> GafPath {
        GafPath {
            steps: vec![
                PathStep { orient: Orient::Forward, node: 0 },
                PathStep { orient: Orient::Forward, node: 1 },
                PathStep { orient: Orient::Forward, node: 3 },
            ],
            path_len: 7,
            path_start: 0,
            path_end: 7,
            num_matches: 99,
            aln_len: 100,
        }
    }
    fn q(seq: &str) -> Query {
        Query {
            name: "q".into(),
            record_id: "q".into(),
            seq: seq.as_bytes().to_vec(),
            path: PathBuf::from("q.fa"),
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
        let before = stitch::stitch(&maj(true), &walked(), &s, &query).unwrap().unwrap();
        let fixed = prefer_majority(&maj(true), &g, &s, &walked(), &query)
            .unwrap()
            .expect("the site should be corrected");
        let after = stitch::stitch(&maj(true), &fixed, &s, &query).unwrap().unwrap();
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

    /// The precondition that separates this from the whole-route repair that failed. A 0.02%
    /// branch is sequencing error, not an allele. Swapping there is `--prune-frac`'s job,
    /// done by removing the branch before the aligner sees it; doing it here would be the
    /// greedy rewrite that cut a 1,279-unitig route to 334.
    #[test]
    fn an_error_branch_site_is_not_touched() {
        let g = diamond(3.0, 12000.0);
        let s = subgraph_of(&g);
        assert!(
            prefer_majority(&maj(true), &g, &s, &walked(), &q("ACGTAG")).unwrap().is_none(),
            "the branch we took is 0.02%; that is not a two-allele site"
        );
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
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec()))).unwrap();
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
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec()))).unwrap();
        let s = subgraph_of(&g);
        assert!(prefer_majority(&maj(true), &g, &s, &walked(), &q("ACGTAG")).unwrap().is_none());
        let c = Config { allow_indel_alts: false, ..maj(true) };
        assert!(prefer_majority(&c, &g, &s, &walked(), &q("ACGTAG")).unwrap().is_none(),
                "and not with indel support switched off either");
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
        let primary = stitch::stitch(&cfg(true), &walked(), &s, &q("ACGTAG")).unwrap().unwrap();
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
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec()))).unwrap();
        let s = subgraph_of(&g);
        let query = q("ACGTAG");
        let primary = stitch::stitch(&cfg(true), &walked(), &s, &query).unwrap().unwrap();
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
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec()))).unwrap();
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
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec()))).unwrap();
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
        let long: String = std::iter::repeat("ACGT").take(20).collect(); // 80 nt, delta +76
        let fa = format!(
            ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n\
             >T_1 ka:f:60.0   L:-:0:-  L:+:3:+\nCGTA\n\
             >T_2 ka:f:40.0   L:-:0:-  L:+:3:+\n{long}\n\
             >T_3 ka:f:100.0  L:-:1:-  L:-:2:-\nGTAG\n"
        );
        let g = parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.into_bytes()))).unwrap();
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
}
