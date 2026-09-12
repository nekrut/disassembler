//! Step 6b: Extract all phased haplotypes via closed-loop detour splicing.
//!
//! Reconstructs genuine intra-host quasispecies / diploid alleles and overcomes reference bias.
//! When GraphAligner chooses a reference-matching minor branch over a divergent majority branch,
//! this stage explores alternative branches at each fork along the walked path. Any branch that
//! clears abundance (`ka:f >= min_abund`) and frequency (`VAF >= min_vaf`) and rejoins the path
//! within `max_detour_hops` is spliced into the path to emit a full-length alternative haplotype.

use crate::config::Config;
use crate::error::Result;
use crate::stages::parse;
use crate::stages::stitch;
use crate::types::{GafPath, Haplotype, HaplotypeMode, Orient, PathStep, Query, Subgraph, UnitigGraph, UnitigId};
use std::collections::{HashMap, HashSet, VecDeque};

/// Compute the average k-mer abundance (ka:f) across a sequence of path steps.
pub fn calculate_mean_ka(graph: &UnitigGraph, steps: &[PathStep]) -> f64 {
    if steps.is_empty() {
        return 0.0;
    }
    let total: f64 = steps
        .iter()
        .map(|s| graph.abundance[s.node as usize] as f64)
        .sum();
    total / steps.len() as f64
}

/// Extract all valid haplotypes for an accession.
///
/// In `HaplotypeMode::One`, returns only the primary path.
/// In `HaplotypeMode::Many`, discovers all closed-loop detours, splices them, deduplicates
/// sequences, ranks them by `mean_ka` descending, and assigns clean ids:
/// - `{accession}_major` for the dominant allele
/// - `{accession}_minor1_fork{id}{strand}`, `{accession}_minor2_...` for alternative alleles.
pub fn extract_haplotypes(
    cfg: &Config,
    graph: &UnitigGraph,
    subgraph: &Subgraph,
    query: &Query,
    gaf_path: &GafPath,
) -> Result<Vec<Haplotype>> {
    // 1. Primary path reconstruction
    let Some(mut primary_hap) = stitch::stitch(cfg, gaf_path, subgraph, query)? else {
        return Ok(Vec::new());
    };

    let steps = &gaf_path.steps;
    if steps.is_empty() {
        return Ok(Vec::new());
    }

    let primary_mean_ka = calculate_mean_ka(graph, steps);
    primary_hap.mean_ka = primary_mean_ka;
    primary_hap.vaf = 1.0;

    // In single-haplotype mode or if path is too short to fork, return primary
    if cfg.haplotypes == HaplotypeMode::One || steps.len() < 2 {
        return Ok(vec![primary_hap]);
    }

    let n = subgraph.nodes.len();
    if n == 0 {
        return Ok(vec![primary_hap]);
    }

    // Local index for subgraph nodes
    let mut local: HashMap<UnitigId, usize> = HashMap::with_capacity(n * 2);
    for (i, &id) in subgraph.nodes.iter().enumerate() {
        local.insert(id, i);
    }

    // Build oriented adjacency graph within subgraph
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
    for v in out.iter_mut() {
        v.sort_unstable();
        v.dedup();
    }

    // Map each step in path to oriented node id (2 * local_idx + strand_bit)
    let mut path_oriented: Vec<Option<u32>> = Vec::with_capacity(steps.len());
    for step in steps {
        if let Some(&j) = local.get(&step.node) {
            let strand_bit = match step.orient {
                Orient::Forward => 1,
                Orient::Reverse => 0,
            };
            path_oriented.push(Some((j * 2 + strand_bit) as u32));
        } else {
            path_oriented.push(None);
        }
    }

    // Index all occurrences of oriented nodes on the path
    let mut path_pos: HashMap<u32, Vec<usize>> = HashMap::new();
    for (idx, &opt) in path_oriented.iter().enumerate() {
        if let Some(oid) = opt {
            path_pos.entry(oid).or_default().push(idx);
        }
    }

    let k_overlap = cfg.graph_k.saturating_sub(1);
    let mut candidates: Vec<(String, Vec<PathStep>, Vec<u8>, f64, f64, f64, usize)> = Vec::new();
    // (tag, steps, seq, cov, pident, mean_ka, local_vaf_percent)
    candidates.push((
        "primary".to_string(),
        steps.clone(),
        primary_hap.seq.clone(),
        primary_hap.cov,
        primary_hap.pident,
        primary_mean_ka,
        100,
    ));

    let mut seen_seqs: HashSet<Vec<u8>> = HashSet::new();
    seen_seqs.insert(primary_hap.seq.clone());

    // Traverse every fork along the path
    for i in 0..(steps.len() - 1) {
        let Some(fork_oid) = path_oriented[i] else { continue };
        let Some(taken_oid) = path_oriented[i + 1] else { continue };

        let branches = &out[fork_oid as usize];
        if branches.len() < 2 || !branches.contains(&taken_oid) {
            continue;
        }

        // Compute total abundance across all outgoing branches from this fork
        let mut total_fork_abund = 0.0f32;
        for &b in branches {
            let b_node = subgraph.nodes[(b as usize) / 2];
            total_fork_abund += graph.abundance[b_node as usize];
        }
        if total_fork_abund <= 0.0 {
            continue;
        }

        let fork_unitig = subgraph.nodes[(fork_oid as usize) / 2];
        let fork_strand_char = if fork_oid % 2 == 1 { '+' } else { '-' };

        // Test EVERY alternative branch (supports N >= 2 branches)
        for &alt in branches {
            if alt == taken_oid {
                continue;
            }
            let alt_node = subgraph.nodes[(alt as usize) / 2];
            let alt_abund = graph.abundance[alt_node as usize];
            let alt_vaf = (alt_abund as f64) / (total_fork_abund as f64);

            if alt_abund < cfg.min_abund || alt_vaf < cfg.min_vaf {
                continue;
            }

            // Forward BFS from alt to find where it rejoins the path at j > i + 1
            let mut queue = VecDeque::new();
            queue.push_back((alt, vec![alt]));
            let mut visited = HashSet::new();
            visited.insert(alt);

            let mut best_detour: Option<(usize, Vec<u32>)> = None;

            while let Some((curr, trail)) = queue.pop_front() {
                if trail.len() > cfg.max_detour_hops {
                    continue;
                }
                if let Some(indices) = path_pos.get(&curr) {
                    if let Some(&j) = indices.iter().find(|&&idx| idx > i + 1) {
                        best_detour = Some((j, trail));
                        break;
                    }
                }
                for &nxt in &out[curr as usize] {
                    if !visited.contains(&nxt) {
                        visited.insert(nxt);
                        let mut next_trail = trail.clone();
                        next_trail.push(nxt);
                        queue.push_back((nxt, next_trail));
                    }
                }
            }

            if let Some((j, detour)) = best_detour {
                // Construct spliced path: steps[..=i] + detour[..detour.len()-1] + steps[j..]
                let mut spliced_steps = Vec::new();
                spliced_steps.extend_from_slice(&steps[..=i]);
                for &d in &detour[..detour.len() - 1] {
                    let d_local = (d as usize) / 2;
                    let d_node = subgraph.nodes[d_local];
                    let d_orient = if d % 2 == 1 { Orient::Forward } else { Orient::Reverse };
                    spliced_steps.push(PathStep { orient: d_orient, node: d_node });
                }
                spliced_steps.extend_from_slice(&steps[j..]);

                let Some(stitched) = stitch::stitch_steps(k_overlap, &spliced_steps, subgraph) else {
                    continue;
                };

                let start = gaf_path.path_start.min(stitched.len());
                let end = if gaf_path.path_end <= stitched.len() {
                    gaf_path.path_end
                } else {
                    stitched.len()
                };
                let target: Vec<u8> = stitched[start..end.max(start)].to_vec();

                let cov = target.len() as f64 / query.seq.len() as f64;
                if cov < cfg.min_cov {
                    continue;
                }

                if seen_seqs.contains(&target) {
                    continue;
                }
                seen_seqs.insert(target.clone());

                let mean_ka = calculate_mean_ka(graph, &spliced_steps);
                let tag = format!("fork{fork_unitig}{fork_strand_char}");
                candidates.push((tag, spliced_steps, target, cov, primary_hap.pident, mean_ka, (alt_vaf * 100.0).round() as usize));
            }
        }
    }

    // Sort all candidates descending by mean_ka (highest abundance is major)
    candidates.sort_by(|a, b| b.5.partial_cmp(&a.5).unwrap_or(std::cmp::Ordering::Equal));

    let total_candidates = candidates.len();
    let mut final_haplotypes = Vec::with_capacity(total_candidates);

    for (rank, (tag, spliced_steps, target, cov, pident, mean_ka, vaf_int)) in candidates.into_iter().enumerate() {
        let id = if total_candidates == 1 {
            format!("{}_primary", subgraph.accession)
        } else if rank == 0 {
            format!("{}_major", subgraph.accession)
        } else {
            format!("{}_minor{}_{}", subgraph.accession, rank, tag)
        };

        final_haplotypes.push(Haplotype {
            id,
            seq: target,
            cov,
            pident,
            nodes: spliced_steps.len(),
            mean_ka,
            vaf: (vaf_int as f64) / 100.0,
        });
    }

    Ok(final_haplotypes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stages::extract;
    use std::path::PathBuf;

    fn build(text: &str) -> UnitigGraph {
        parse::build(&"T".into(), Box::new(std::io::Cursor::new(text.as_bytes().to_vec()))).unwrap()
    }

    fn whole(cfg: &Config, g: &UnitigGraph) -> Subgraph {
        let nodes: Vec<UnitigId> = (0..g.n_unitigs as u32).collect();
        let seqs = extract::collect_seqs(g, &nodes).unwrap();
        extract::write_gfa(cfg, g, nodes, seqs).unwrap()
    }

    fn cleanup(s: &Subgraph) {
        let _ = std::fs::remove_file(&s.gfa_path);
    }

    #[test]
    fn multi_branch_fork_extracts_all_haplotypes_ranked_by_abundance() {
        let mut c = Config::default();
        c.graph_k = 4; // k-1 = 3 bp overlap
        c.min_cov = 0.50;
        c.min_abund = 5.0;
        c.min_vaf = 0.05;
        c.haplotypes = HaplotypeMode::Many;

        // Tri-way fork from node 0:
        // Branch 1: node 1 (abundance 10.0, walked by GAF)
        // Branch 2: node 2 (abundance 80.0, major allele)
        // Branch 3: node 3 (abundance 40.0, secondary allele)
        // All branches merge into node 4.
        let g = build(concat!(
            ">T_0 ka:f:100.0  L:+:1:+  L:+:2:+  L:+:3:+\nACGT\n",
            ">T_1 ka:f:10.0   L:-:0:-  L:+:4:+\nCGTA\n",
            ">T_2 ka:f:80.0   L:-:0:-  L:+:4:+\nCGTC\n",
            ">T_3 ka:f:40.0   L:-:0:-  L:+:4:+\nCGTG\n",
            ">T_4 ka:f:100.0  L:-:1:-  L:-:2:-  L:-:3:-\nGTAC\n",
        ));
        let sub = whole(&c, &g);

        let query = Query {
            name: "q".into(),
            record_id: "q".into(),
            seq: b"ACGTATAC".to_vec(),
            path: PathBuf::from("q.fa"),
        };

        // GAF path walked branch 1: 0 -> 1 -> 4
        let gaf = GafPath {
            steps: vec![
                PathStep { orient: Orient::Forward, node: 0 },
                PathStep { orient: Orient::Forward, node: 1 },
                PathStep { orient: Orient::Forward, node: 4 },
            ],
            path_len: 6,
            path_start: 0,
            path_end: 6,
            num_matches: 6,
            aln_len: 6,
        };

        // 1. Test Many mode
        let haps = extract_haplotypes(&c, &g, &sub, &query, &gaf).unwrap();
        assert_eq!(haps.len(), 3, "Expected 3 haplotypes from tri-way fork, got {}", haps.len());
        // Highest abundance (node 2, ka=80) should be major!
        assert_eq!(haps[0].id, "T_major");
        assert!(haps[0].mean_ka > haps[1].mean_ka);
        assert!(haps[1].mean_ka > haps[2].mean_ka);
        assert_eq!(haps[1].id, "T_minor1_fork0+");
        assert_eq!(haps[2].id, "T_minor2_primary");

        // 2. Test One mode
        c.haplotypes = HaplotypeMode::One;
        let single_hap = extract_haplotypes(&c, &g, &sub, &query, &gaf).unwrap();
        assert_eq!(single_hap.len(), 1, "Expected exactly 1 haplotype in One mode");
        assert_eq!(single_hap[0].id, "T_primary");

        cleanup(&sub);
    }
}
