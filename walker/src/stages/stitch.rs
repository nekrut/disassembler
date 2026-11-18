//! Step 6: turn a path through the graph into a sequence.
//!
//! The path is a list of node traversals; the sequence is their concatenation with the k-1
//! overlap removed from every node after the first, and each `<` node reverse-complemented
//! before it is appended.
//!
//! The arithmetic, checked against a real `ERR2028820` path:
//!
//! ```text
//! 115 nodes, lengths summing to 4,708 bp
//! 4,708 - (30 x 114) = 1,288 bp    == GAF column 7
//! ```
//!
//! Requiring a connected path is what keeps chimeras out. Concatenating BLAST hits in
//! coordinate order instead produces sequences that never existed in any organism, which
//! then read as strong positive selection downstream.
//!
//! IN   GafPath, Subgraph, Query
//! OUT  Haplotype, if coverage clears `cfg.min_cov`

use crate::config::{Config, PrimaryMetric};
use crate::error::{Error, Result};
use crate::types::{GafPath, Haplotype, Orient, PathStep, Query, Subgraph, UnitigGraph};

/// Reverse complement. Sequences come from the 2-bit store so they are pure ACGT; anything
/// else would already have been rejected by `parse::build`.
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// Walk the path, trimming `k-1` from each node after the first.
///
/// # The overlap
///
/// Consecutive nodes on a path share exactly `k-1` bases, because that is what an edge in a
/// de Bruijn graph of order k means. Concatenating naively would duplicate those bases, so
/// every node after the first contributes only its tail. Checked against a real path:
///
/// ```text
/// 115 nodes, lengths summing to 4,708 bp
/// 4,708 - (30 x 114) = 1,288 bp   == GafPath::path_len
/// ```
///
/// # The slice
///
/// `path_start` and `path_end` are offsets into this stitched sequence, not into the query.
/// GraphAligner reports where along the path the alignment actually sits, and the flanks
/// outside it are graph context that is not part of the gene.
///
/// # Why a path at all
///
/// Concatenating BLAST hits in coordinate order instead produces sequences assembled from
/// fragments of different strains, which never existed in any organism. README section 4
/// measures that at 87% chimeric artifacts, and a chimera looks exactly like a rapidly
/// evolving gene to a selection test. Requiring a connected path makes them impossible.
/// Concatenate unitigs along path steps with k-1 overlap removed.
pub fn stitch_steps(
    k_overlap: usize,
    steps: &[PathStep],
    subgraph: &Subgraph,
) -> Option<Vec<u8>> {
    if steps.is_empty() {
        return Some(Vec::new());
    }
    let mut stitched: Vec<u8> = Vec::new();

    for (i, step) in steps.iter().enumerate() {
        let Ok(idx) = subgraph.nodes.binary_search(&step.node) else {
            return None;
        };
        let seq = &subgraph.seqs[idx];
        let oriented;
        let s: &[u8] = match step.orient {
            Orient::Forward => seq,
            Orient::Reverse => {
                oriented = revcomp(seq);
                &oriented
            }
        };
        if i == 0 {
            stitched.extend_from_slice(s);
        } else if s.len() > k_overlap {
            stitched.extend_from_slice(&s[k_overlap..]);
        }
    }
    Some(stitched)
}

/// One reference's candidate haplotype and the GAF path it was stitched from, tagged with
/// the reference id and the reference's position in the anchor pool. `pick_primary` uses
/// the pool index for tie-breaking so the choice is stable across runs.
///
/// The `path` field lets the caller carry both halves through the fan-out without having
/// to re-align to recover the winning path (which would double the align cost).
#[derive(Clone)]
pub struct Candidate {
    pub record_id: String,
    pub pool_index: usize,
    /// Length of this reference, the denominator of the `Matched` metric.
    pub ref_len: usize,
    pub path: Option<GafPath>,
    pub haplotype: Option<Haplotype>,
}

/// Rank per-reference candidates and pick the winning primary.
///
/// `metric` picks the primary sort key ([`PrimaryMetric`]):
/// * `Matched` (default): `num_matches / ref_len` desc, then pident desc. The fraction of the
///   reference's own bases the alignment matched; never above 1.0.
/// * `Cov`: walked length / reference length desc, then pident desc. Uncapped, so a shorter
///   reference can win on a ratio above 1.0.
/// * `Pident`: pident desc, then coverage desc.
/// * `CovXPident`: `cov * pident/100`, then cov, then pident.
///
/// `stitch` already returns `None` for a candidate whose coverage fell below `min_cov`,
/// so an under-covered reference is filtered out before it can win regardless of metric.
/// Candidates within `tie_eps` of the best score are treated as tied and the **deepest** of
/// them wins, by `Haplotype::mean_ka`. Without that, a sample holding two strains that each
/// equal a panel record gives two candidates scoring 1.0 and pident 100, and the record
/// listed first in the panel wins: a 5% strain can be published as `_primary`, taking its
/// `anchor=` tag and every VCF coordinate with it, while the 95% strain leaves no trace.
/// Nothing downstream repairs that, because `prefer_majority` skips a fork whose taken branch
/// is below `--min-vaf`, so a 5% route is never moved to the 95% one.
///
/// The band is measured on the metric's own value, normalised to a fraction (pident is
/// divided by 100) so one `tie_eps` means the same thing whichever metric is in use.
/// `tie_eps = 0` disables the depth tie-break entirely and restores the pre-0.2.0 ordering,
/// including for exact ties.
///
/// Every metric falls back to `pool_index` as the final tie-break so the first reference
/// in the anchor pool wins on a true tie, which keeps primary selection stable. Depth is
/// compared with a stable sort, so candidates of equal depth keep that order too.
///
/// Returns `None` when no reference produced a stitchable haplotype.
pub fn pick_primary(
    candidates: Vec<Candidate>,
    metric: PrimaryMetric,
    tie_eps: f64,
) -> Option<(String, GafPath, Haplotype)> {
    let mut with_haps: Vec<(usize, String, GafPath, Haplotype, f64)> = candidates
        .into_iter()
        .filter_map(|c| match (c.path, c.haplotype) {
            (Some(p), Some(h)) => {
                let matched = if c.ref_len > 0 {
                    p.num_matches as f64 / c.ref_len as f64
                } else {
                    0.0
                };
                Some((c.pool_index, c.record_id, p, h, matched))
            }
            _ => None,
        })
        .collect();
    // `partial_cmp` on the (cov, pident) f64s never sees NaN in practice: both are
    // computed from bounded integer ratios in `stitch`. `Equal` on NaN is a defensive
    // fallback that keeps the sort total.
    let cmpf = |a: f64, b: f64| b.partial_cmp(&a).unwrap_or(std::cmp::Ordering::Equal);
    with_haps.sort_by(|a, b| {
        let (ac, ap) = (a.3.cov, a.3.pident);
        let (bc, bp) = (b.3.cov, b.3.pident);
        let primary = match metric {
            PrimaryMetric::Matched => cmpf(a.4, b.4).then_with(|| cmpf(ap, bp)),
            PrimaryMetric::Cov => cmpf(ac, bc).then_with(|| cmpf(ap, bp)),
            PrimaryMetric::Pident => cmpf(ap, bp).then_with(|| cmpf(ac, bc)),
            PrimaryMetric::CovXPident => {
                cmpf(ac * ap / 100.0, bc * bp / 100.0)
                    .then_with(|| cmpf(ac, bc))
                    .then_with(|| cmpf(ap, bp))
            }
        };
        primary.then(a.0.cmp(&b.0))
    });
    // The metric's own value, as a fraction, so one `tie_eps` means the same for each.
    let key = |c: &(usize, String, GafPath, Haplotype, f64)| match metric {
        PrimaryMetric::Matched => c.4,
        PrimaryMetric::Cov => c.3.cov,
        PrimaryMetric::Pident => c.3.pident / 100.0,
        PrimaryMetric::CovXPident => c.3.cov * c.3.pident / 100.0,
    };
    if tie_eps > 0.0 {
        if let Some(first) = with_haps.first() {
            let best = key(first);
            // The list is already sorted by the metric, so the band is a prefix.
            let n = with_haps.iter().take_while(|c| best - key(c) <= tie_eps).count();
            // Stable, so equal depths keep the metric-then-pool order decided above.
            with_haps[..n].sort_by(|a, b| cmpf(a.3.mean_ka, b.3.mean_ka));
        }
    }
    with_haps.into_iter().next().map(|(_, name, p, h, _)| (name, p, h))
}

pub fn stitch(
    cfg: &Config,
    graph: &UnitigGraph,
    path: &GafPath,
    subgraph: &Subgraph,
    query: &Query,
) -> Result<Option<Haplotype>> {
    if path.steps.is_empty() {
        return Ok(None);
    }
    if query.seq.is_empty() {
        return Err(Error::Usage("query has no sequence".into()));
    }

    let k_overlap = cfg.graph_k.saturating_sub(1);
    let Some(stitched) = stitch_steps(k_overlap, &path.steps, subgraph) else {
        return Ok(None);
    };

    // Clamp rather than fail: a path_end past the stitched length takes the tail.
    let start = path.path_start.min(stitched.len());
    let end = if path.path_end <= stitched.len() {
        path.path_end
    } else {
        stitched.len()
    };
    let target: Vec<u8> = stitched[start..end.max(start)].to_vec();

    let cov = target.len() as f64 / query.seq.len() as f64;
    // What `--min-cov` actually gates on (WA-05). `cov` above is the walked length, which an
    // indel in the sample moves in the wrong direction: a deletion makes a gene that is fully
    // there look absent, and an insertion makes a partial alignment look complete. The
    // reference span is the fraction of the reference the alignment found, which is what the
    // flag has always been documented as meaning.
    let ref_cov = path.query_span() as f64 / query.seq.len() as f64;

    // Tripwire for the multi-record-query bug fixed in `align::run`: that published another
    // gene's span under this gene's length and reached cov 2.097. A genuine alignment can
    // exceed 1.0 slightly, because insertions in the path push `path_end` past
    // `path_start + qlen`, so this warns on a real excess rather than asserting.
    const COV_EXCESS: f64 = 1.02;
    if cov > COV_EXCESS {
        eprintln!(
            "[warn] {}: reconstructed span is {:.3}x the query ({} nt from a {} nt query {:?}); \
             a coverage this far above 1.0 means the wrong alignment row was used",
            subgraph.accession, cov, target.len(), query.seq.len(), query.record_id
        );
    }

    let pident = path.pident();

    if ref_cov < cfg.min_cov {
        return Ok(None);
    }

    Ok(Some(Haplotype {
        id: format!("{}_primary", subgraph.accession),
        seq: target,
        cov,
        ref_cov,
        ref_span: (path.query_start, path.query_end),
        pident,
        // Set by the caller when it splices a path; `stitch` itself is always handed the
        // path the GAF describes.
        pident_inherited: false,
        nodes: path.steps.len(),
        // Was hardcoded 0.0 while being published in every FASTA header and in summary.json.
        mean_ka: crate::stages::extract::path_mean_ka(graph, &path.steps),
        vaf: 1.0,
        // stitch does not know how big the anchor pool is; the pipeline sets this field
        // to the winning reference id only when a pool of >1 reference is in play.
        anchor_used: None,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{PathStep, UnitigId};
    use std::path::PathBuf;

    fn sub(nodes: Vec<UnitigId>, seqs: Vec<&str>) -> Subgraph {
        let s = Subgraph {
            accession: "ACC".into(),
            nodes,
            seqs: seqs.iter().map(|x| x.as_bytes().to_vec()).collect(),
            // Points at nothing, and `owns_gfa: false` means Drop leaves it alone. Before
            // TW-01 this relied on the unlink failing, which it does not do as root.
            gfa_path: PathBuf::from("/tmp/logan-stitch-test-nonexistent.gfa"),
            n_edges: 0,
            owns_gfa: false,
        };
        s
    }
    fn q(seq: &str) -> Query {
        Query { name: "q".into(), record_id: "q".into(), seq: seq.as_bytes().to_vec(), path: PathBuf::from("q.fa"), n_records: 1 }
    }
    fn path(steps: &[(Orient, u32)], len: usize, start: usize, end: usize) -> GafPath {
        GafPath {
            steps: steps.iter().map(|&(orient, node)| PathStep { orient, node }).collect(),
            query_start: 0, query_end: end.saturating_sub(start),
            path_len: len, path_start: start, path_end: end,
            num_matches: 99, aln_len: 100,
        }
    }
    fn cfg(k: usize, min_cov: f64) -> Config {
        Config { graph_k: k, min_cov, ..Config::default() }
    }
    /// A graph whose unitig `i` carries `ka[i]`. `stitch` reads only `abundance` from it;
    /// the sequences it stitches come from the `Subgraph`.
    fn gph(ka: &[f32]) -> UnitigGraph {
        let fa: String = ka
            .iter()
            .enumerate()
            .map(|(i, a)| format!(">T_{i} ka:f:{a}\nACGT\n"))
            .collect();
        crate::stages::parse::build(&"T".into(), Box::new(std::io::Cursor::new(fa.into_bytes())), u64::MAX)
            .unwrap()
    }
    /// For the cases where abundance is irrelevant; wide enough for every node id used here.
    fn anyg() -> UnitigGraph {
        gph(&[10.0; 78])
    }

    #[test]
    fn revcomp_is_correct_and_involutive() {
        assert_eq!(revcomp(b"ACGT"), b"ACGT".to_vec(), "ACGT is its own reverse complement");
        assert_eq!(revcomp(b"AAAC"), b"GTTT".to_vec());
        assert_eq!(revcomp(&revcomp(b"ACGGTCA")), b"ACGGTCA".to_vec());
    }

    // ---- pick_primary (multi-anchor) ----

    fn cand(name: &str, pool_index: usize, cov: f64, pident: f64) -> Candidate {
        cand_m(name, pool_index, cov, pident, 1, 1000)
    }

    /// `cand` with the alignment's matched-base count and the reference length, which is what
    /// the `Matched` metric ranks on.
    fn cand_m(name: &str, pool_index: usize, cov: f64, pident: f64,
              num_matches: usize, ref_len: usize) -> Candidate {
        // Real GafPath and Haplotype shells: pick_primary reads path.num_matches and
        // hap.cov / hap.pident for ranking. Everything else is placeholder.
        let path = GafPath {
            steps: vec![PathStep { orient: Orient::Forward, node: 1 }],
            query_start: 0, query_end: 1,
            path_len: 1, path_start: 0, path_end: 1,
            num_matches, aln_len: 1,
        };
        let hap = Haplotype {
            id: name.into(),
            seq: vec![],
            cov,
            ref_cov: cov,
            ref_span: (0, 0),
            pident,
            pident_inherited: false,
            nodes: 1,
            mean_ka: 0.0,
            vaf: 1.0,
            anchor_used: None,
        };
        Candidate {
            record_id: name.into(),
            pool_index,
            ref_len,
            path: Some(path),
            haplotype: Some(hap),
        }
    }

    #[test]
    fn pick_primary_returns_none_on_empty_pool() {
        assert!(pick_primary(vec![], PrimaryMetric::Cov, 0.0).is_none());
    }

    #[test]
    fn pick_primary_returns_none_when_all_haplotypes_are_none() {
        let none_cand = Candidate {
            record_id: "a".into(),
            pool_index: 0,
            ref_len: 1000,
            path: None,
            haplotype: None,
        };
        assert!(pick_primary(vec![none_cand], PrimaryMetric::Cov, 0.0).is_none());
    }

    #[test]
    fn pick_primary_picks_highest_coverage() {
        let cs = vec![
            cand("low", 0, 0.30, 99.5),
            cand("high", 1, 0.80, 99.5),
            cand("mid", 2, 0.55, 99.5),
        ];
        let (name, _, hap) = pick_primary(cs, PrimaryMetric::Cov, 0.0).unwrap();
        assert_eq!(name, "high");
        assert!((hap.cov - 0.80).abs() < 1e-9);
    }

    #[test]
    fn pick_primary_breaks_ties_on_pident() {
        let cs = vec![
            cand("weaker", 0, 0.80, 98.5),
            cand("stronger", 1, 0.80, 99.9),
        ];
        let (name, _, _) = pick_primary(cs, PrimaryMetric::Cov, 0.0).unwrap();
        assert_eq!(name, "stronger");
    }

    #[test]
    fn pick_primary_breaks_full_ties_on_pool_order() {
        // Same cov and pident: the earlier pool_index wins, so the anchor set's ordering
        // is what a downstream sieve or reader can rely on.
        let cs = vec![
            cand("second", 5, 0.60, 99.0),
            cand("first",  1, 0.60, 99.0),
        ];
        let (name, _, _) = pick_primary(cs, PrimaryMetric::Cov, 0.0).unwrap();
        assert_eq!(name, "first");
    }

    #[test]
    fn pick_primary_ignores_candidates_without_a_haplotype() {
        let mut cs = vec![
            cand("valid", 3, 0.50, 99.0),
        ];
        cs.push(Candidate {
            record_id: "no_stitch".into(),
            pool_index: 0,
            ref_len: 1000,
            path: None,
            haplotype: None,
        });
        let (name, _, _) = pick_primary(cs, PrimaryMetric::Cov, 0.0).unwrap();
        assert_eq!(name, "valid", "None candidates must not win by having pool_index 0");
    }

    /// The measured SARS-CoV-2 case that motivated the flag: a longer walk on a mid-cluster
    /// reference beats a shorter walk on the phylogenetically closest one under `Cov`, and
    /// the flag flips the pick without either candidate having to change.
    #[test]
    fn primary_metric_flips_the_winner_between_cov_and_pident() {
        // XBB.1.5-like: fuller walk, moderate identity.
        // JN.1-like: shorter walk, higher identity.
        let cs = || vec![
            cand("xbb15", 0, 0.99, 99.6),
            cand("jn1",   1, 0.71, 99.9),
        ];
        assert_eq!(pick_primary(cs(), PrimaryMetric::Cov, 0.0).unwrap().0, "xbb15");
        assert_eq!(pick_primary(cs(), PrimaryMetric::Pident, 0.0).unwrap().0, "jn1");
    }

    /// Cov-x-pident is a joint score, so a small identity gap can be overridden by a
    /// large coverage advantage and vice versa. The two candidates below score:
    ///   xbb15: 0.99 * 99.6 / 100 = 0.9860
    ///   jn1:   0.71 * 99.9 / 100 = 0.7093
    /// so xbb15 still wins. Flip the values around and the winner flips too.
    #[test]
    fn primary_metric_cov_x_pident_uses_the_joint_score() {
        let cs = vec![
            cand("xbb15", 0, 0.99, 99.6),
            cand("jn1",   1, 0.71, 99.9),
        ];
        assert_eq!(pick_primary(cs, PrimaryMetric::CovXPident, 0.0).unwrap().0, "xbb15");

        // Now a case where jn1 has enough of an identity edge that the joint score flips.
        let cs = vec![
            cand("xbb15", 0, 0.71, 90.0),
            cand("jn1",   1, 0.70, 99.9),
        ];
        assert_eq!(pick_primary(cs, PrimaryMetric::CovXPident, 0.0).unwrap().0, "jn1");
    }

    #[test]
    fn primary_metric_still_breaks_ties_on_pool_order() {
        // Identical numbers under every metric; the pool_index tie-break must survive.
        let cs = vec![
            cand("second", 5, 0.90, 99.0),
            cand("first",  1, 0.90, 99.0),
        ];
        for m in [PrimaryMetric::Matched, PrimaryMetric::Cov, PrimaryMetric::Pident,
                  PrimaryMetric::CovXPident] {
            assert_eq!(pick_primary(cs.clone(), m, 0.002).unwrap().0, "first",
                "metric {m:?} must still fall back to pool order");
        }
    }

    /// `cand_m` with a route depth, which is what breaks a tie from 0.2.0 on.
    fn cand_ka(name: &str, pool_index: usize, num_matches: usize, ref_len: usize,
               ka: f64) -> Candidate {
        let mut c = cand_m(name, pool_index, 1.0, 100.0, num_matches, ref_len);
        if let Some(h) = c.haplotype.as_mut() {
            h.mean_ka = ka;
        }
        c
    }

    /// Two strains that each match a panel record exactly: the deeper route wins, and the
    /// answer does not depend on which record the panel happens to list first.
    ///
    /// This is the co-infection and wastewater case. Both candidates score 1.0 matched and
    /// pident 100, so before 0.2.0 panel order decided, and a 5% strain could be published as
    /// `_primary` with its `anchor=` tag and every VCF coordinate. Nothing downstream
    /// repaired it: `prefer_majority` skips a fork whose taken branch is under `--min-vaf`,
    /// so a 5% route was never moved to the 95% one.
    #[test]
    fn matched_tie_goes_to_the_deeper_route_not_panel_order() {
        for (a, b) in [((0, 500.0), (1, 9500.0)), ((0, 9500.0), (1, 500.0))] {
            let (shallow, deep) = if a.1 < b.1 { (a, b) } else { (b, a) };
            let cs = vec![
                cand_ka("shallow", shallow.0, 900, 900, shallow.1),
                cand_ka("deep", deep.0, 900, 900, deep.1),
            ];
            let got = pick_primary(cs, PrimaryMetric::Matched, 0.002).unwrap().0;
            assert_eq!(got, "deep", "panel order must not decide between tied references");
        }
        // Off, the pre-0.2.0 behaviour: the first in the panel wins, whatever its depth.
        let cs = vec![
            cand_ka("shallow", 0, 900, 900, 500.0),
            cand_ka("deep", 1, 900, 900, 9500.0),
        ];
        assert_eq!(pick_primary(cs, PrimaryMetric::Matched, 0.0).unwrap().0, "shallow");
    }

    /// Depth breaks ties and nothing else. A deeper but less similar reference must not
    /// displace one the sample actually matches, or the tie-break would become the metric.
    #[test]
    fn a_deeper_but_worse_matching_reference_still_loses() {
        let cs = vec![
            // 0.98 against 1.00, two orders of magnitude outside the 0.002 band.
            cand_ka("deep_but_98pc", 0, 882, 900, 9500.0),
            cand_ka("exact", 1, 900, 900, 100.0),
        ];
        assert_eq!(pick_primary(cs, PrimaryMetric::Matched, 0.002).unwrap().0, "exact");
    }

    /// The case that made `Matched` the default. The sample carries the full 1,200 nt gene.
    /// `exact` is that gene; `short` is the same gene with a 9 nt deletion (1,191 nt), listed
    /// first in the pool. The walk is 1,200 nt against both, so `Cov` gives `short`
    /// 1200/1191 = 1.0076 and picks it over the exact match. Both references have every one of
    /// their own bases matched (1.0), and the 9 extra sample bases lower `short`'s identity to
    /// 1191/1200, so `Matched` picks `exact`.
    #[test]
    fn matched_prefers_the_exact_reference_over_one_with_a_deletion() {
        let cs = || vec![
            cand_m("short", 0, 1200.0 / 1191.0, 1191.0 / 1200.0 * 100.0, 1191, 1191),
            cand_m("exact", 1, 1.0, 100.0, 1200, 1200),
        ];
        assert_eq!(pick_primary(cs(), PrimaryMetric::Cov, 0.0).unwrap().0, "short");
        assert_eq!(pick_primary(cs(), PrimaryMetric::Matched, 0.0).unwrap().0, "exact");
    }

    /// `Pident` can prefer a truncated walk that matches its reference cleanly (measured:
    /// SRR17435695, 0.999-coverage walk at 99.1% vs a 0.919 walk at 99.4%). `Matched` charges
    /// the truncated walk for the reference bases it never reached.
    #[test]
    fn matched_does_not_prefer_a_truncated_walk() {
        let cs = || vec![
            cand_m("full", 0, 0.999, 99.1, 3784, 3822),
            cand_m("truncated", 1, 0.919, 99.4, 3492, 3822),
        ];
        assert_eq!(pick_primary(cs(), PrimaryMetric::Pident, 0.0).unwrap().0, "truncated");
        assert_eq!(pick_primary(cs(), PrimaryMetric::Matched, 0.0).unwrap().0, "full");
    }

    #[test]
    fn matched_is_the_default_metric() {
        assert_eq!(Config::default().primary_metric, PrimaryMetric::Matched);
    }

    /// `mean_ka` was hardcoded to 0.0 while being published as `ka=` in every FASTA header
    /// and in summary.json.
    #[test]
    fn mean_ka_is_the_mean_over_the_path() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1, 2], vec!["ACGTA", "GTACC"]);
        let p = path(&[(Orient::Forward, 1), (Orient::Forward, 2)], 7, 0, 7);

        let h = stitch(&c, &gph(&[0.0, 10.0, 20.0]), &p, &s, &q("ACGTACC")).unwrap().unwrap();
        assert!((h.mean_ka - 15.0).abs() < 1e-9, "mean of 10 and 20, got {}", h.mean_ka);

        // It tracks the graph rather than being any fixed value.
        let h = stitch(&c, &gph(&[0.0, 4.0, 100.0]), &p, &s, &q("ACGTACC")).unwrap().unwrap();
        assert!((h.mean_ka - 52.0).abs() < 1e-9, "got {}", h.mean_ka);

        // A repeated step counts twice: unweighted per step, as documented.
        let p2 = path(&[(Orient::Forward, 1), (Orient::Forward, 1), (Orient::Forward, 2)], 7, 0, 7);
        let h = stitch(&c, &gph(&[0.0, 10.0, 40.0]), &p2, &s, &q("ACGTACC")).unwrap().unwrap();
        assert!((h.mean_ka - 20.0).abs() < 1e-9, "(10+10+40)/3, got {}", h.mean_ka);
    }

    /// Mean, not `path_depth`'s median: one deep repeat-collapsed unitig must move `ka=`,
    /// or the two statistics would be indistinguishable on the data that separates them.
    #[test]
    fn mean_ka_is_not_the_median() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1, 2, 3], vec!["ACGTA", "GTACC", "ACCGT"]);
        let p = path(
            &[(Orient::Forward, 1), (Orient::Forward, 2), (Orient::Forward, 3)],
            9, 0, 9,
        );
        let g = gph(&[0.0, 1.0, 2.0, 300.0]);
        let h = stitch(&c, &g, &p, &s, &q("ACGTACCGT")).unwrap().unwrap();
        assert!((h.mean_ka - 101.0).abs() < 1e-9, "mean of 1,2,300, got {}", h.mean_ka);
        assert_eq!(
            crate::stages::extract::path_depth(&g, &p.steps),
            2.0,
            "path_depth stays the median and must not move with the outlier"
        );
    }

    /// Consecutive nodes share k-1 bases, so only the first contributes in full.
    #[test]
    fn trims_the_overlap_from_every_node_after_the_first() {
        let c = cfg(4, 0.0); // k-1 = 3
        let s = sub(vec![1, 2], vec!["ACGTA", "GTACC"]);
        // node1 "ACGTA" whole, node2 contributes "CC" (drops first 3)
        let p = path(&[(Orient::Forward, 1), (Orient::Forward, 2)], 7, 0, 7);
        let h = stitch(&c, &anyg(), &p, &s, &q("ACGTACC")).unwrap().unwrap();
        assert_eq!(h.seq, b"ACGTACC");
    }

    #[test]
    fn reverse_steps_are_reverse_complemented() {
        let c = cfg(4, 0.0);
        // revcomp("GGTAC") = "GTACC"; after trimming 3 it contributes "CC"
        let s = sub(vec![1, 2], vec!["ACGTA", "GGTAC"]);
        let p = path(&[(Orient::Forward, 1), (Orient::Reverse, 2)], 7, 0, 7);
        let h = stitch(&c, &anyg(), &p, &s, &q("ACGTACC")).unwrap().unwrap();
        assert_eq!(h.seq, b"ACGTACC");
    }

    /// The arithmetic that must hold for any path: sum of lengths minus (k-1) per junction.
    #[test]
    fn stitched_length_matches_the_overlap_formula() {
        let c = cfg(31, 0.0);
        let a = "A".repeat(40);
        let b = "C".repeat(45);
        let d = "G".repeat(50);
        let s = sub(vec![1, 2, 3], vec![&a, &b, &d]);
        let expected = 40 + 45 + 50 - 30 * 2;
        let p = path(&[(Orient::Forward, 1), (Orient::Forward, 2), (Orient::Forward, 3)],
                     expected, 0, expected);
        let h = stitch(&c, &anyg(), &p, &s, &q(&"N".repeat(expected))).unwrap().unwrap();
        assert_eq!(h.seq.len(), expected);
    }

    /// path_start and path_end index the stitched path, not the query.
    #[test]
    fn slices_with_path_coordinates() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1], vec!["ACGTACGTAC"]);
        let p = path(&[(Orient::Forward, 1)], 10, 2, 7);
        let h = stitch(&c, &anyg(), &p, &s, &q("ACGTA")).unwrap().unwrap();
        assert_eq!(h.seq, b"GTACG", "expected stitched[2..7]");
    }

    /// Clamps rather than failing when path_end runs past the stitched sequence.
    #[test]
    fn path_end_past_the_end_takes_the_tail() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1], vec!["ACGTACGTAC"]);
        let p = path(&[(Orient::Forward, 1)], 10, 6, 9999);
        let h = stitch(&c, &anyg(), &p, &s, &q("ACGT")).unwrap().unwrap();
        assert_eq!(h.seq, b"GTAC");
    }

    #[test]
    fn coverage_is_span_over_query_and_identity_is_from_the_gaf() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1], vec!["ACGTACGTAC"]);
        let mut p = path(&[(Orient::Forward, 1)], 10, 0, 10);
        p.num_matches = 95;
        p.aln_len = 100;
        let h = stitch(&c, &anyg(), &p, &s, &q(&"N".repeat(20))).unwrap().unwrap();
        assert!((h.cov - 0.5).abs() < 1e-12, "10 bases over a 20 bp query");
        assert!((h.pident - 95.0).abs() < 1e-12);
        assert_eq!(h.nodes, 1);
        assert_eq!(h.id, "ACC_primary");
    }

    /// `path`, with the reference span given separately from the walked span.
    ///
    /// The ordinary helper ties the two together, which is exactly the assumption WA-05 was
    /// about: an indel in the sample makes them differ.
    fn path_spanning(
        steps: &[(Orient, u32)],
        len: usize, start: usize, end: usize,
        qstart: usize, qend: usize,
    ) -> GafPath {
        GafPath {
            steps: steps.iter().map(|&(orient, node)| PathStep { orient, node }).collect(),
            query_start: qstart, query_end: qend,
            path_len: len, path_start: start, path_end: end,
            num_matches: 99, aln_len: 100,
        }
    }

    /// WA-05. A sample that carries the whole gene with an in-frame deletion is present, and
    /// `--min-cov` is documented as the fraction of the reference a primary must cover.
    ///
    /// The gate used to measure the walked length instead: a 300 nt gene with a 60 nt
    /// deletion walks 240 bases, scores 0.80 against the default 0.85 and was published as
    /// `no_path`, "the gene is not here". That biases presence calls against exactly the
    /// alleles a deletion study is looking for. The alignment spans the whole reference, so
    /// the reference coverage is 1.00 and it is kept.
    #[test]
    fn a_deletion_does_not_make_a_present_gene_look_absent() {
        let c = cfg(4, 0.85);
        let s = sub(vec![1], vec![&"A".repeat(240)]);
        // Walks 240 bases; the alignment covers reference bases 0..300 of 300.
        let p = path_spanning(&[(Orient::Forward, 1)], 240, 0, 240, 0, 300);
        let h = stitch(&c, &anyg(), &p, &s, &q(&"N".repeat(300))).unwrap().expect("kept");
        assert_eq!(h.seq.len(), 240);
        assert!((h.cov - 0.8).abs() < 1e-12, "cov still reports the walked length");
        assert!((h.ref_cov - 1.0).abs() < 1e-12, "the whole reference was found");
        assert_eq!(h.ref_span, (0, 300));
    }

    /// The mirror image, and the reason this is a correctness fix rather than a loosening:
    /// an alignment that reaches only reference base 240 of 300 is 0.80 of the reference
    /// however much sequence it walked. Walking a 40 nt insertion used to lift it to 0.933
    /// and publish it as a primary.
    #[test]
    fn an_insertion_does_not_make_a_partial_alignment_look_complete() {
        let c = cfg(4, 0.85);
        let s = sub(vec![1], vec![&"A".repeat(280)]);
        let p = path_spanning(&[(Orient::Forward, 1)], 280, 0, 280, 0, 240);
        assert!(
            stitch(&c, &anyg(), &p, &s, &q(&"N".repeat(300))).unwrap().is_none(),
            "0.80 of the reference is below --min-cov 0.85, whatever was walked"
        );
    }

    /// The span is reported as the interval, not only as a fraction: a walk that stops 60
    /// bases from the 3' end is a different finding from one missing 60 in the middle, and
    /// `ref_cov` alone cannot tell them apart.
    #[test]
    fn the_reference_interval_is_published_not_just_the_fraction() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1], vec![&"A".repeat(100)]);
        let p = path_spanning(&[(Orient::Forward, 1)], 100, 0, 100, 40, 160);
        let h = stitch(&c, &anyg(), &p, &s, &q(&"N".repeat(200))).unwrap().expect("kept");
        assert_eq!(h.ref_span, (40, 160));
        assert!((h.ref_cov - 0.6).abs() < 1e-12);
    }

    #[test]
    fn below_min_cov_yields_no_haplotype() {
        let c = cfg(4, 0.85);
        let s = sub(vec![1], vec!["ACGTACGTAC"]);
        let p = path(&[(Orient::Forward, 1)], 10, 0, 10);
        // 10 bases against a 100 bp query is 10% coverage.
        assert!(stitch(&c, &anyg(), &p, &s, &q(&"N".repeat(100))).unwrap().is_none());
    }

    /// GraphAligner naming a node absent from the GFA means the alignment cannot be
    /// realised. That yields no haplotype rather than an error.
    #[test]
    fn unknown_node_in_path_yields_none() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1], vec!["ACGTACGTAC"]);
        let p = path(&[(Orient::Forward, 1), (Orient::Forward, 77)], 10, 0, 10);
        assert!(stitch(&c, &anyg(), &p, &s, &q("ACGT")).unwrap().is_none());
    }

    #[test]
    fn empty_path_yields_none() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1], vec!["ACGT"]);
        assert!(stitch(&c, &anyg(), &path(&[], 0, 0, 0), &s, &q("ACGT")).unwrap().is_none());
    }

    #[test]
    fn empty_query_is_an_error_not_a_division_by_zero() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1], vec!["ACGT"]);
        let p = path(&[(Orient::Forward, 1)], 4, 0, 4);
        assert!(stitch(&c, &anyg(), &p, &s, &q("")).is_err());
    }
}

/// Full-chain comparison against the committed golden fixture. This is the end-to-end
/// concordance check: fetch through stitch, in Rust. The fixture was first captured under
/// GraphAligner 1.0.20 and recaptured under 1.0.17b on 2026-10-18, which
/// added the final base that 1.0.19+ clips (docs/issue_graphaligner_end_clipping.md). Run it
/// with GraphAligner 1.0.17b on PATH.
///
/// ```text
/// PATH=$REPO/.condaenv-walker/bin:$PATH LOGAN_TEST_DIR=/dir/with/zst \
///     cargo test --release -- --ignored --nocapture golden_haplotype
/// ```
#[cfg(test)]
mod golden_haplotype {
    use super::*;
    use crate::io::input::read_query;
    use crate::stages::{align, anchor, extract, parse, stream};
    use crate::types::FetchedUnitigs;
    use sha2::{Digest, Sha256};
    use std::path::Path;

    #[test]
    #[ignore]
    fn rust_reproduces_the_golden_haplotype() {
        let Ok(dir) = std::env::var("LOGAN_TEST_DIR") else {
            eprintln!("SKIP: set LOGAN_TEST_DIR");
            return;
        };
        let Ok(compressed) = std::fs::read(format!("{dir}/ERR2028820.unitigs.fa.zst")) else {
            eprintln!("SKIP: no data");
            return;
        };
        let mut cfg = Config::default();
        // `write_gfa` no longer creates a missing scratch directory (WF-17),
        // so a test that names its own creates it, as a run does.
        cfg.scratch_dir = std::env::temp_dir().join(format!("logan_golden_{}", std::process::id()));
        std::fs::create_dir_all(&cfg.scratch_dir).unwrap();
        if align::resolve_binary(&cfg).is_err() {
            eprintln!("SKIP: GraphAligner not on PATH");
            return;
        }

        let f = FetchedUnitigs {
            index: 0,
            accession: "ERR2028820".into(),
            compressed_len: compressed.len(),
            compressed: crate::types::Compressed::Memory(compressed),
            fetch_secs: 0.0,
            payload: crate::types::Payload::default(),
        };
        let graph = { let rd = stream::decoder(&f).unwrap(); parse::build(&f.accession, rd, u64::MAX).unwrap() };
        let query = read_query(Path::new("../queries/pfcrt_only.fasta")).unwrap();
        let idx = anchor::index_query(&cfg, &query).unwrap();
        let a = anchor::find(&cfg, &graph, &idx).unwrap();
        let nodes = extract::expand(&cfg, &graph, &a).unwrap();
        let seqs = extract::collect_seqs(&graph, &nodes).unwrap();
        let subgraph = extract::write_gfa(&cfg, &graph, nodes, seqs).unwrap();
        let gaf = subgraph.gfa_path.with_extension("gaf");
        let mut gaf_map = align::run(&cfg, &subgraph, std::slice::from_ref(&query), &gaf).unwrap();
        let gp = gaf_map.remove(query.record_id.as_str()).expect("a path");
        let h = stitch(&cfg, &graph, &gp, &subgraph, &query).unwrap().expect("a haplotype");

        let sha: String = Sha256::digest(&h.seq)
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        eprintln!(
            "  {} len {} cov {:.6} pident {:.3} nodes {}",
            h.id, h.seq.len(), h.cov, h.pident, h.nodes
        );
        eprintln!("  sha256 {sha}");

        // Values of bench/golden/ERR2028820.pfcrt.primary.fasta (2026-10-18, GraphAligner
        // 1.0.17b). The original 1.0.20 capture in ERR2028820.pfcrt.golden.json is the same
        // sequence without its final base: 1,274 nt, cov 0.999216, the same identity.
        assert_eq!(h.id, "ERR2028820_primary");
        assert_eq!(h.seq.len(), 1275, "haplotype length");
        assert!((h.cov - 1.0).abs() < 1e-6, "coverage {}", h.cov);
        assert!((h.pident - 99.529).abs() < 1e-3, "identity {}", h.pident);
        assert_eq!(
            sha, "a1a952438ec40865c944ccb3dadf0095860f13f4e4a082dbb4416d93c4b1a2e3",
            "haplotype does not match the golden fixture"
        );
    }
}
