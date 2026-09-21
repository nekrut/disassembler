//! Step 7: phase forks into co-existing variants.
//!
//! Replaces `logan_path_walker.py:238-275` (`phase_superbubbles`). Runs off the graph and
//! subgraph, NOT the GAF, so it is independent of GraphAligner.
//!
//! A fork whose branches both clear `ka:f >= min_abund` and `vaf >= min_vaf` represents two
//! populations genuinely present in one sample. That is the pipeline's whole reason for
//! working on unitigs rather than contigs, which would have discarded the minor branch.
//!
//! # Bug to fix here, not to port
//!
//! The Python version builds its adjacency while discarding the strand fields (`:244-247`):
//!
//! ```python
//! for s1, tgt, s2 in edge_list:   # s1 and s2 thrown away
//!     outgoing[u].add(tgt)
//! ```
//!
//! Because a bidirected edge is recorded in both endpoints' headers, `A->B` and `B->A` are
//! the same edge seen twice, and the reconvergence test then reports the fork as its own
//! merge node. Measured: `fork == merge` on 93.1% of `ERR2028820`'s 202 bubbles and 58.5%
//! of `SRR8209082`'s 11,825.
//!
//! `UnitigGraph::targets` keeps the strand bits precisely so this can be done correctly.
//! Sequence the work: port to byte-identical behaviour first, verify against
//! `bench/golden/`, then fix the strand handling as a separate reviewable change
//! with a deliberate fixture regeneration.
//!
//! Also unlike the Python version, keep the per-bubble detail. There, only `len(bubbles)`
//! survives into the summary (`:317`, `:385`), so the quasispecies claim cannot be audited
//! from the output at all.
//!
//! IN   UnitigGraph, Subgraph
//! OUT  Vec<Bubble>

use crate::config::Config;
use crate::error::Result;
use crate::stages::parse;
use crate::types::{Branch, Bubble, Strand, Subgraph, UnitigGraph, UnitigId};
use std::collections::HashMap;

/// Detect bubbles and phase their branches.
///
/// # What a bubble is
///
/// One end of a unitig has two or more ways forward, and those alternatives meet again at a
/// common next step. That diamond means the sample contains two versions of the same region.
/// Comparing their read coverage says whether both are real.
///
/// # The orientation fix
///
/// Logan edges are bidirected: `L:s1:t:s2` leaves one specific *end* of the source and
/// enters one specific *end* of the target. A unitig has two ends, so a neighbour at each is
/// an ordinary chain, not a choice.
///
/// `logan_path_walker.py:244-247` builds its adjacency while discarding `s1` and `s2`:
///
/// ```python
/// for s1, tgt, s2 in edge_list:   # s1 and s2 thrown away
///     outgoing[u].add(tgt)
/// ```
///
/// Two consequences. A plain chain `x <- u -> y` reads as a two-way fork. And because Logan
/// writes every edge into both endpoints' headers, `u -> v` and `v -> u` are the same edge
/// seen twice, so both "branches" appear to lead back to `u`, which is then reported as its
/// own merge node. Measured: `fork == merge` on 93.1% of `ERR2028820`'s 202 bubbles and
/// 58.5% of `SRR8209082`'s 11,825.
///
/// Here the graph is over *oriented* nodes, `(unitig, strand)`. An edge in `u`'s header is
/// the traversal `(u, s1) -> (t, s2)`; the reverse traversal `(t, !s2) -> (u, !s1)` arrives
/// independently from `t`'s own header, so the oriented graph is closed under reversal with
/// each direction recorded exactly once.
///
/// Verified on `ERR2028820`: unitig 2933 carries `L:-:443893:+` and `L:+:1239391:-`, which
/// the Python version reads as a fork with two branches. Oriented, those are one exit from
/// each end, so `(2933,-)` and `(2933,+)` each have out-degree 1 and there is no fork.
///
/// # Output differs from Python deliberately
///
/// This is the behaviour change the plan sequences after a faithful port. The golden
/// fixtures pin the old numbers and must be regenerated, with the diff reviewed.
pub fn phase(cfg: &Config, graph: &UnitigGraph, subgraph: &Subgraph) -> Result<Vec<Bubble>> {
    let n = subgraph.nodes.len();
    if n == 0 {
        return Ok(Vec::new());
    }

    // Local index for the subgraph, so the oriented adjacency is a dense Vec rather than a
    // map. Oriented node (local i, strand s) lives at 2*i + s.
    let mut local: HashMap<UnitigId, usize> = HashMap::with_capacity(n * 2);
    for (i, &id) in subgraph.nodes.iter().enumerate() {
        local.insert(id, i);
    }

    let mut out: Vec<Vec<u32>> = vec![Vec::new(); n * 2];
    for (i, &id) in subgraph.nodes.iter().enumerate() {
        for &packed in graph.edges(id) {
            let (target, src_fwd, dst_fwd) = parse::unpack_edge(packed);
            let Some(&j) = local.get(&target) else { continue };
            let from = i * 2 + usize::from(src_fwd);
            let to = (j * 2 + usize::from(dst_fwd)) as u32;
            out[from].push(to);
        }
    }
    // Sorting and deduplicating makes the whole stage order-independent, and collapses the
    // parallel edges Logan emits when two unitigs meet at more than one end.
    for v in out.iter_mut() {
        v.sort_unstable();
        v.dedup();
    }

    let unlocal = |oriented: u32| -> (UnitigId, Strand) {
        let i = (oriented as usize) / 2;
        let s = if oriented % 2 == 1 { Strand::Forward } else { Strand::Reverse };
        (subgraph.nodes[i], s)
    };

    let mut bubbles = Vec::new();
    for fork in 0..(n * 2) {
        let branches = &out[fork];
        if branches.len() < 2 {
            continue;
        }

        // Which oriented nodes are reachable from at least two of the branches.
        let mut hits: HashMap<u32, u32> = HashMap::new();
        for &b in branches {
            for &c in &out[b as usize] {
                *hits.entry(c).or_insert(0) += 1;
            }
        }
        let mut reconverged: Vec<u32> = hits
            .into_iter()
            .filter(|&(c, count)| count >= 2 && c != fork as u32)
            .map(|(c, _)| c)
            .collect();
        if reconverged.is_empty() {
            continue;
        }
        // Sorted for a deterministic emission order, unlike Python's hash-ordered Counter.
        // NOT to pick a single winner: this used to take `reconverged[0]` and discard the
        // rest, so a fork with several reconvergence points reported one bubble and silently
        // lost the others, in every orientation, before the dedup below was even reached.
        reconverged.sort_unstable();

        for &merge in &reconverged {
            // Every bubble exists twice in a bidirected graph, once per traversal direction:
            // `(0,+) -> {1,2} -> (3,+)` is the same diamond as `(3,-) -> {1,2} -> (0,-)`.
            // Flipping an oriented node's strand is `^ 1`, so the mirror of `(fork, merge)`
            // is `(merge ^ 1, fork ^ 1)`, and the graph being closed under reversal means
            // that pair is enumerated too.
            //
            // Keep the lexicographically smaller of the two *pairs*: exactly one of
            // `(f, m)` and `(m^1, f^1)` is smaller, and a self-mirrored bubble (`m == f^1`)
            // compares equal and is kept once.
            //
            // This is EQUIVALENT to the older `fork > merge ^ 1`, which compared only the
            // forks: given both directions are enumerated with the same merge, `f <= m^1`
            // and `m^1 <= f` cannot both fail. Checked exhaustively over 300k random
            // reversal-closed graphs, zero disagreements. The older form was therefore never
            // the defect here; it only looked wrong because taking `reconverged[0]` meant the
            // two directions were comparing *different* merges. The pair form is kept because
            // it does not depend on that invariant, so a future change that again reports one
            // merge per fork cannot silently turn the dedup into a bubble-loser.
            let f = fork as u32;
            if (f, merge) > (merge ^ 1, f ^ 1) {
                continue;
            }

            // Only the exits that actually reach THIS merge are alleles of this bubble.
            // `total` and `kept` used to iterate every exit from the fork, so a third branch
            // heading somewhere else was published as an allele of a bubble it is not part
            // of, and its coverage diluted the real branches' VAFs, which can push a true
            // minority allele under `min_vaf` and suppress the bubble via `kept.len() < 2`.
            // `out[b]` is sorted and deduped above, so this is a binary search.
            let in_bubble: Vec<u32> = branches
                .iter()
                .copied()
                .filter(|&b| out[b as usize].binary_search(&merge).is_ok())
                .collect();
            if in_bubble.len() < 2 {
                continue;
            }

            // Coverage share is computed over every branch of this bubble, then the gates are
            // applied, which matches the Python ordering: a branch filtered out by a gate
            // still contributed to the total.
            let total: f64 = in_bubble
                .iter()
                .map(|&b| graph.abundance[unlocal(b).0 as usize] as f64)
                .sum();
            if total <= 0.0 {
                continue;
            }

            let mut kept: Vec<Branch> = Vec::new();
            for &b in &in_bubble {
                let (node, strand) = unlocal(b);
                let ab = graph.abundance[node as usize];
                let vaf = ab as f64 / total;
                if ab >= cfg.min_abund && vaf >= cfg.min_vaf {
                    kept.push(Branch { node, strand, abundance: ab, vaf });
                }
            }
            if kept.len() < 2 {
                continue;
            }
            kept.sort_by_key(|a| (a.node, a.strand));

            let (fork_id, fork_strand) = unlocal(f);
            let (merge_id, merge_strand) = unlocal(merge);
            bubbles.push(Bubble {
                fork: fork_id,
                fork_strand,
                merge: merge_id,
                merge_strand,
                branches: kept,
            });
        }
    }

    // A fork can now contribute more than one bubble, so the merge is part of the key or the
    // order would depend on discovery.
    bubbles.sort_by_key(|a| (a.fork, a.fork_strand, a.merge, a.merge_strand));
    Ok(bubbles)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stages::extract;
    
    use std::path::PathBuf;

    fn build(fa: &str) -> UnitigGraph {
        parse::build(&"T".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec()))).unwrap()
    }
    fn cfg() -> Config {
        let mut d = std::env::temp_dir();
        d.push(format!("logan_phase_test_{}", std::process::id()));
        Config { k: 31, min_abund: 5.0, min_vaf: 0.10, scratch_dir: d, ..Config::default() }
    }
    /// Build the whole graph as one subgraph, so phasing sees every node.
    fn whole(cfg: &Config, g: &UnitigGraph) -> Subgraph {
        let nodes: Vec<UnitigId> = (0..g.n_unitigs as u32).collect();
        let seqs = extract::collect_seqs(g, &nodes).unwrap();
        
        extract::write_gfa(cfg, g, nodes, seqs).unwrap()
    }
    fn cleanup(s: &Subgraph) { std::fs::remove_file(&s.gfa_path).ok(); }

    /// The exact shape the Python implementation miscounts. Unitig 0 has one neighbour off
    /// each end, which is a chain, not a choice.
    #[test]
    fn a_neighbour_at_each_end_is_not_a_fork() {
        let c = cfg();
        let g = build(concat!(
            ">T_0 ka:f:50.0   L:-:1:+  L:+:2:-\nACGT\n",
            ">T_1 ka:f:50.0   L:-:0:+\nTTTT\n",
            ">T_2 ka:f:50.0   L:+:0:-\nGGGG\n",
        ));
        let sub = whole(&c, &g);
        let b = phase(&c, &g, &sub).unwrap();
        cleanup(&sub);
        assert!(b.is_empty(), "python reports a bubble here; oriented, there is no fork");
    }

    /// A real bubble: one end of unitig 0 leads to either 1 or 2, and both reach 3.
    #[test]
    fn detects_a_genuine_two_branch_bubble() {
        let c = cfg();
        let g = build(concat!(
            ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n",
            ">T_1 ka:f:60.0   L:-:0:-  L:+:3:+\nTTTT\n",
            ">T_2 ka:f:40.0   L:-:0:-  L:+:3:+\nGGGG\n",
            ">T_3 ka:f:100.0  L:-:1:-  L:-:2:-\nCCCC\n",
        ));
        let sub = whole(&c, &g);
        let b = phase(&c, &g, &sub).unwrap();
        cleanup(&sub);
        assert_eq!(b.len(), 1, "expected exactly one bubble, got {}", b.len());
        assert_eq!(b[0].fork, 0);
        assert_eq!(b[0].merge, 3, "merge must be a real node, never the fork itself");
        assert_ne!(b[0].merge, b[0].fork);
        let nodes: Vec<UnitigId> = b[0].branches.iter().map(|x| x.node).collect();
        assert_eq!(nodes, vec![1, 2]);
        assert!((b[0].branches[0].vaf - 0.6).abs() < 1e-9);
        assert!((b[0].branches[1].vaf - 0.4).abs() < 1e-9);
    }

    /// Every oriented traversal Logan records implies its reverse, because the edge appears
    /// in both endpoints' headers. A fixture that is not closed under reversal is not a graph
    /// this stage can be handed, so assert it rather than assume it.
    fn assert_reversal_closed(g: &UnitigGraph) {
        use std::collections::HashSet;
        let mut edges: HashSet<(u32, u32)> = HashSet::new();
        for i in 0..g.n_unitigs as u32 {
            for &packed in g.edges(i) {
                let (target, src_fwd, dst_fwd) = parse::unpack_edge(packed);
                edges.insert((i * 2 + u32::from(src_fwd), target * 2 + u32::from(dst_fwd)));
            }
        }
        for &(f, to) in &edges {
            assert!(
                edges.contains(&(to ^ 1, f ^ 1)),
                "fixture is not reversal-closed: {f}->{to} has no mirror {}->{}",
                to ^ 1,
                f ^ 1
            );
        }
    }

    /// A fork with two reconvergence points must report BOTH bubbles.
    ///
    /// The loop used to sort `reconverged` and take element 0, so every reconvergence point
    /// but the smallest was discarded before the dedup was even reached.
    #[test]
    fn a_fork_with_two_merges_reports_both() {
        let c = cfg();
        // 0+ forks to 1+ and 2+; both reach 3+ AND both reach 4+. Two genuine diamonds
        // sharing one fork.
        let g = build(concat!(
            ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n",
            ">T_1 ka:f:60.0   L:-:0:-  L:+:3:+  L:+:4:+\nTTTT\n",
            ">T_2 ka:f:40.0   L:-:0:-  L:+:3:+  L:+:4:+\nGGGG\n",
            ">T_3 ka:f:100.0  L:-:1:-  L:-:2:-\nCCCC\n",
            ">T_4 ka:f:100.0  L:-:1:-  L:-:2:-\nAAAA\n",
        ));
        assert_reversal_closed(&g);
        let sub = whole(&c, &g);
        let b = phase(&c, &g, &sub).unwrap();
        cleanup(&sub);

        let merges: Vec<UnitigId> = b
            .iter()
            .filter(|x| x.fork == 0 && x.fork_strand == Strand::Forward)
            .map(|x| x.merge)
            .collect();
        assert_eq!(merges, vec![3, 4], "both reconvergence points, got {merges:?}");
    }

    /// A bubble whose two traversal directions choose different merges survives.
    ///
    /// Fork `(2,-)` reconverges at both `(2,+)` and `(3,-)`; taking only the smallest meant
    /// `(2,-) -> (3,-)` was never considered, and the mirror direction `(3,+)` picked a
    /// different merge and was then skipped by the dedup, so the diamond appeared in NEITHER
    /// direction of `bubbles.json`. Found by brute-force search over reversal-closed graphs;
    /// see review finding 2.1.
    #[test]
    fn a_bubble_whose_directions_disagree_on_the_merge_survives() {
        let c = cfg();
        let g = build(concat!(
            ">T_0 ka:f:50.0  L:-:1:-  L:-:1:+  L:-:2:+  L:-:3:-  L:+:2:-  L:+:2:+\nACGT\n",
            ">T_1 ka:f:50.0  L:-:0:+  L:-:3:-  L:+:0:+  L:+:2:+\nTTTT\n",
            ">T_2 ka:f:50.0  L:-:0:-  L:-:0:+  L:-:1:-  L:+:0:-\nGGGG\n",
            ">T_3 ka:f:50.0  L:+:0:+  L:+:1:+\nCCCC\n",
        ));
        assert_reversal_closed(&g);
        let sub = whole(&c, &g);
        let b = phase(&c, &g, &sub).unwrap();
        cleanup(&sub);

        // fork (2,-) reconverges at both (2,+) and (3,-); the mirror direction (3,+) picks
        // (2,+). The (2,-) -> (3,-) diamond appeared in neither direction before.
        let found = b.iter().any(|x| {
            x.fork == 2
                && x.fork_strand == Strand::Reverse
                && x.merge == 3
                && x.merge_strand == Strand::Reverse
        });
        let seen: Vec<(UnitigId, Strand, UnitigId, Strand)> =
            b.iter().map(|x| (x.fork, x.fork_strand, x.merge, x.merge_strand)).collect();
        assert!(found, "the (2,-) -> (3,-) diamond must be reported; got {seen:?}");
    }

    /// A branch that leaves the fork but never reaches the chosen merge is not an allele of
    /// that bubble, and must not dilute the real alleles' VAFs.
    #[test]
    fn a_branch_that_misses_the_merge_is_not_an_allele() {
        let c = cfg();
        // 0+ exits to 1+, 2+ and 3+. Only 1+ and 2+ reach 4+; 3+ goes to 5+ and carries
        // enough coverage to push the minority allele under min_vaf if it were counted.
        let g = build(concat!(
            ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+  L:+:3:+\nACGT\n",
            ">T_1 ka:f:80.0   L:-:0:-  L:+:4:+\nTTTT\n",
            ">T_2 ka:f:20.0   L:-:0:-  L:+:4:+\nGGGG\n",
            ">T_3 ka:f:400.0  L:-:0:-  L:+:5:+\nCCCC\n",
            ">T_4 ka:f:100.0  L:-:1:-  L:-:2:-\nAAAA\n",
            ">T_5 ka:f:400.0  L:-:3:-\nACAC\n",
        ));
        assert_reversal_closed(&g);
        let sub = whole(&c, &g);
        let b = phase(&c, &g, &sub).unwrap();
        cleanup(&sub);

        let bubble = b
            .iter()
            .find(|x| x.fork == 0 && x.fork_strand == Strand::Forward && x.merge == 4)
            .expect("the 1/2 diamond must be reported");
        let nodes: Vec<UnitigId> = bubble.branches.iter().map(|x| x.node).collect();
        assert_eq!(nodes, vec![1, 2], "unitig 3 does not reach merge 4");

        // VAFs over 80 and 20, not over 80/20/400: with 3 counted, the minority allele is
        // 20/500 = 0.04, under the 0.10 floor, and `kept.len() < 2` would drop the bubble.
        assert!((bubble.branches[0].vaf - 0.8).abs() < 1e-9, "got {}", bubble.branches[0].vaf);
        assert!((bubble.branches[1].vaf - 0.2).abs() < 1e-9, "got {}", bubble.branches[1].vaf);
    }

    /// A branch below the coverage floor is a sequencing error, not a variant.
    #[test]
    fn low_coverage_branch_is_dropped() {
        let c = cfg();
        let g = build(concat!(
            ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n",
            ">T_1 ka:f:95.0   L:-:0:-  L:+:3:+\nTTTT\n",
            ">T_2 ka:f:1.8    L:-:0:-  L:+:3:+\nGGGG\n",
            ">T_3 ka:f:100.0  L:-:1:-  L:-:2:-\nCCCC\n",
        ));
        let sub = whole(&c, &g);
        let b = phase(&c, &g, &sub).unwrap();
        cleanup(&sub);
        assert!(b.is_empty(), "1.8x coverage is below min_abund, leaving one branch");
    }

    /// Real but rare: above the coverage floor yet below the frequency floor.
    #[test]
    fn branch_below_vaf_floor_is_dropped() {
        let c = cfg();
        let g = build(concat!(
            ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n",
            ">T_1 ka:f:1000.0 L:-:0:-  L:+:3:+\nTTTT\n",
            ">T_2 ka:f:9.0    L:-:0:-  L:+:3:+\nGGGG\n",
            ">T_3 ka:f:100.0  L:-:1:-  L:-:2:-\nCCCC\n",
        ));
        let sub = whole(&c, &g);
        let b = phase(&c, &g, &sub).unwrap();
        cleanup(&sub);
        assert!(b.is_empty(), "9/1009 is under 10% VAF");
    }

    /// Diverging without meeting again is not a bubble.
    #[test]
    fn fork_without_reconvergence_is_not_a_bubble() {
        let c = cfg();
        let g = build(concat!(
            ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n",
            ">T_1 ka:f:60.0   L:-:0:-\nTTTT\n",
            ">T_2 ka:f:40.0   L:-:0:-\nGGGG\n",
        ));
        let sub = whole(&c, &g);
        let b = phase(&c, &g, &sub).unwrap();
        cleanup(&sub);
        assert!(b.is_empty());
    }

    #[test]
    fn empty_subgraph_yields_nothing() {
        let c = cfg();
        let g = build(">T_0 ka:f:9.0\nACGT\n");
        let sub = Subgraph {
            accession: "T".into(), nodes: vec![], seqs: vec![],
            gfa_path: PathBuf::from("/dev/null"), n_edges: 0,
        };
        assert!(phase(&c, &g, &sub).unwrap().is_empty());
    }

    #[test]
    fn output_is_deterministic() {
        let c = cfg();
        let g = build(concat!(
            ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+\nACGT\n",
            ">T_1 ka:f:60.0   L:-:0:-  L:+:3:+\nTTTT\n",
            ">T_2 ka:f:40.0   L:-:0:-  L:+:3:+\nGGGG\n",
            ">T_3 ka:f:100.0  L:-:1:-  L:-:2:-\nCCCC\n",
        ));
        let sub = whole(&c, &g);
        let a = phase(&c, &g, &sub).unwrap();
        let b = phase(&c, &g, &sub).unwrap();
        cleanup(&sub);
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!((x.fork, x.merge), (y.fork, y.merge));
            assert_eq!(
                x.branches.iter().map(|v| v.node).collect::<Vec<_>>(),
                y.branches.iter().map(|v| v.node).collect::<Vec<_>>()
            );
        }
    }
}

/// What the orientation fix changes on real data, against the Python counts.
///
/// `LOGAN_TEST_DIR=/dir/with/zst cargo test --release -- --ignored --nocapture phase_real`
#[cfg(test)]
mod phase_real {
    use super::*;
    use crate::io::input::read_query;
    use crate::stages::{anchor, extract, stream};
    use crate::types::FetchedUnitigs;
    use std::path::Path;

    #[test]
    #[ignore]
    fn compares_against_python_bubble_counts() {
        let Ok(dir) = std::env::var("LOGAN_TEST_DIR") else {
            eprintln!("SKIP: set LOGAN_TEST_DIR");
            return;
        };
        // (accession, query, python bubble count, python fork==merge count)
        let cases = [
            ("ERR2028820", "../queries/pfcrt_only.fasta", 202usize, 188usize),
            ("SRR8209082", "../queries/query_hiv1_pol.fasta", 11825, 6917),
        ];
        for (acc, qpath, py_total, py_degenerate) in cases {
            let Ok(compressed) = std::fs::read(format!("{dir}/{acc}.unitigs.fa.zst")) else {
                eprintln!("SKIP: no data for {acc}");
                continue;
            };
            let f = FetchedUnitigs { index: 0, accession: acc.into(), compressed_len: compressed.len(), compressed, fetch_secs: 0.0 };
            let graph = { let rd = stream::decoder(&f).unwrap(); parse::build(&f.accession, rd).unwrap() };
            let query = read_query(Path::new(qpath)).unwrap();
            let mut cfg = Config::default();
            cfg.scratch_dir = std::env::temp_dir().join(format!("logan_phase_real_{}", std::process::id()));

            let idx = anchor::index_query(&cfg, &query).unwrap();
            let a = anchor::find(&cfg, &graph, &idx).unwrap();
            let nodes = extract::expand(&cfg, &graph, &a).unwrap();
            let seqs = extract::collect_seqs(&graph, &nodes).unwrap();
            let sub = extract::write_gfa(&cfg, &graph, nodes, seqs).unwrap();

            let t0 = std::time::Instant::now();
            let bubbles = phase(&cfg, &graph, &sub).unwrap();
            let elapsed = t0.elapsed();
            std::fs::remove_file(&sub.gfa_path).ok();

            let degenerate = bubbles.iter().filter(|b| b.fork == b.merge).count();
            let multi = bubbles.iter().filter(|b| b.branches.len() > 2).count();
            eprintln!(
                "  {acc}: subgraph {} nodes -> {} bubbles in {:.4}s  \
                 (python: {} bubbles, {} of them fork==merge)",
                sub.nodes.len(), bubbles.len(), elapsed.as_secs_f64(), py_total, py_degenerate
            );
            eprintln!("     fork==merge: {degenerate}   branches>2: {multi}");
            if let Some(b) = bubbles.first() {
                eprintln!(
                    "     example: fork {}{:?} -> merge {}{:?}, branches {:?}",
                    b.fork, b.fork_strand, b.merge, b.merge_strand,
                    b.branches.iter().map(|x| (x.node, x.abundance, (x.vaf * 1000.0).round() / 1000.0)).collect::<Vec<_>>()
                );
            }
            assert_eq!(degenerate, 0, "a fork must never be its own merge once orientation is honoured");
            for b in &bubbles {
                assert!(b.branches.len() >= 2);
                assert!(b.branches.iter().all(|x| x.abundance >= cfg.min_abund));
                assert!(b.branches.iter().all(|x| x.vaf >= cfg.min_vaf));
            }
        }
    }
}

#[cfg(test)]
mod phase_diagnostic {
    use super::*;
    use crate::io::input::read_query;
    use crate::stages::{anchor, extract, stream};
    use crate::types::FetchedUnitigs;
    use std::path::Path;

    #[test]
    #[ignore]
    fn shows_why_the_count_drops() {
        let Ok(dir) = std::env::var("LOGAN_TEST_DIR") else { return };
        for (acc, qpath) in [
            ("ERR2028820", "../queries/pfcrt_only.fasta"),
            ("SRR8209082", "../queries/query_hiv1_pol.fasta"),
        ] {
            let Ok(compressed) = std::fs::read(format!("{dir}/{acc}.unitigs.fa.zst")) else { continue };
            let f = FetchedUnitigs { index: 0, accession: acc.into(), compressed_len: compressed.len(), compressed, fetch_secs: 0.0 };
            let graph = { let rd = stream::decoder(&f).unwrap(); parse::build(&f.accession, rd).unwrap() };
            let query = read_query(Path::new(qpath)).unwrap();
            let mut cfg = Config::default();
            cfg.scratch_dir = std::env::temp_dir().join(format!("logan_diag_{}", std::process::id()));
            let idx = anchor::index_query(&cfg, &query).unwrap();
            let a = anchor::find(&cfg, &graph, &idx).unwrap();
            let nodes = extract::expand(&cfg, &graph, &a).unwrap();
            let seqs = extract::collect_seqs(&graph, &nodes).unwrap();
            let sub = extract::write_gfa(&cfg, &graph, nodes, seqs).unwrap();

            // Count, over the same subgraph, how the candidate pool shrinks at each condition.
            let n = sub.nodes.len();
            let mut local = std::collections::HashMap::new();
            for (i, &id) in sub.nodes.iter().enumerate() { local.insert(id, i); }
            let mut out: Vec<Vec<u32>> = vec![Vec::new(); n * 2];
            let mut strandless: Vec<std::collections::HashSet<u32>> = vec![Default::default(); n];
            for (i, &id) in sub.nodes.iter().enumerate() {
                for &p in graph.edges(id) {
                    let (t, sf, df) = parse::unpack_edge(p);
                    if let Some(&j) = local.get(&t) {
                        out[i * 2 + usize::from(sf)].push((j * 2 + usize::from(df)) as u32);
                        strandless[i].insert(j as u32);
                    }
                }
            }
            for v in out.iter_mut() { v.sort_unstable(); v.dedup(); }
            let strandless_forks = strandless.iter().filter(|s| s.len() >= 2).count();
            let oriented_forks = out.iter().filter(|v| v.len() >= 2).count();

            let bubbles = phase(&cfg, &graph, &sub).unwrap();
            std::fs::remove_file(&sub.gfa_path).ok();
            eprintln!("{acc}: {n} subgraph nodes");
            eprintln!("  strand-less 'forks' (>=2 neighbours, what python counts): {strandless_forks}");
            eprintln!("  oriented forks (one end with >=2 exits):                  {oriented_forks}");
            eprintln!("  of those, reconverging + passing coverage/VAF gates:      {}", bubbles.len());
            for b in bubbles.iter().take(12) {
                eprintln!(
                    "    fork {:>7}{:?} -> merge {:>7}{:?}  branches {:?}",
                    b.fork, b.fork_strand, b.merge, b.merge_strand,
                    b.branches.iter().map(|x| (x.node, x.abundance)).collect::<Vec<_>>()
                );
            }
        }
    }
}
