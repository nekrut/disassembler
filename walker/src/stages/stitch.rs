//! Step 6: turn a path through the graph into a sequence.
//!
//! Replaces `logan_path_walker.py:195-236`. The path is a list of node traversals; the
//! sequence is their concatenation with the k-1 overlap removed from every node after the
//! first, and each `<` node reverse-complemented before it is appended.
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

use crate::config::Config;
use crate::error::{Error, Result};
use crate::types::{GafPath, Haplotype, Orient, PathStep, Query, Subgraph};

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
/// Mirrors `logan_path_walker.py:202-234`.
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

pub fn stitch(
    cfg: &Config,
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

    // Clamp exactly as Python does: a path_end past the stitched length takes the tail.
    let start = path.path_start.min(stitched.len());
    let end = if path.path_end <= stitched.len() {
        path.path_end
    } else {
        stitched.len()
    };
    let target: Vec<u8> = stitched[start..end.max(start)].to_vec();

    let cov = target.len() as f64 / query.seq.len() as f64;
    let pident = if path.aln_len > 0 {
        (path.num_matches as f64 / path.aln_len as f64) * 100.0
    } else {
        0.0
    };

    if cov < cfg.min_cov {
        return Ok(None);
    }

    Ok(Some(Haplotype {
        id: format!("{}_primary", subgraph.accession),
        seq: target,
        cov,
        pident,
        nodes: path.steps.len(),
        mean_ka: 0.0,
        vaf: 1.0,
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
            // Points at nothing; Drop's remove_file is allowed to fail.
            gfa_path: PathBuf::from("/tmp/logan-stitch-test-nonexistent.gfa"),
            n_edges: 0,
        };
        s
    }
    fn q(seq: &str) -> Query {
        Query { name: "q".into(), record_id: "q".into(), seq: seq.as_bytes().to_vec(), path: PathBuf::from("q.fa") }
    }
    fn path(steps: &[(Orient, u32)], len: usize, start: usize, end: usize) -> GafPath {
        GafPath {
            steps: steps.iter().map(|&(orient, node)| PathStep { orient, node }).collect(),
            path_len: len, path_start: start, path_end: end,
            num_matches: 99, aln_len: 100,
        }
    }
    fn cfg(k: usize, min_cov: f64) -> Config {
        Config { graph_k: k, min_cov, ..Config::default() }
    }

    #[test]
    fn revcomp_is_correct_and_involutive() {
        assert_eq!(revcomp(b"ACGT"), b"ACGT".to_vec(), "ACGT is its own reverse complement");
        assert_eq!(revcomp(b"AAAC"), b"GTTT".to_vec());
        assert_eq!(revcomp(&revcomp(b"ACGGTCA")), b"ACGGTCA".to_vec());
    }

    /// Consecutive nodes share k-1 bases, so only the first contributes in full.
    #[test]
    fn trims_the_overlap_from_every_node_after_the_first() {
        let c = cfg(4, 0.0); // k-1 = 3
        let s = sub(vec![1, 2], vec!["ACGTA", "GTACC"]);
        // node1 "ACGTA" whole, node2 contributes "CC" (drops first 3)
        let p = path(&[(Orient::Forward, 1), (Orient::Forward, 2)], 7, 0, 7);
        let h = stitch(&c, &p, &s, &q("ACGTACC")).unwrap().unwrap();
        assert_eq!(h.seq, b"ACGTACC");
    }

    #[test]
    fn reverse_steps_are_reverse_complemented() {
        let c = cfg(4, 0.0);
        // revcomp("GGTAC") = "GTACC"; after trimming 3 it contributes "CC"
        let s = sub(vec![1, 2], vec!["ACGTA", "GGTAC"]);
        let p = path(&[(Orient::Forward, 1), (Orient::Reverse, 2)], 7, 0, 7);
        let h = stitch(&c, &p, &s, &q("ACGTACC")).unwrap().unwrap();
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
        let h = stitch(&c, &p, &s, &q(&"N".repeat(expected))).unwrap().unwrap();
        assert_eq!(h.seq.len(), expected);
    }

    /// path_start and path_end index the stitched path, not the query.
    #[test]
    fn slices_with_path_coordinates() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1], vec!["ACGTACGTAC"]);
        let p = path(&[(Orient::Forward, 1)], 10, 2, 7);
        let h = stitch(&c, &p, &s, &q("ACGTA")).unwrap().unwrap();
        assert_eq!(h.seq, b"GTACG", "expected stitched[2..7]");
    }

    /// Python clamps rather than failing when path_end runs past the stitched sequence.
    #[test]
    fn path_end_past_the_end_takes_the_tail() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1], vec!["ACGTACGTAC"]);
        let p = path(&[(Orient::Forward, 1)], 10, 6, 9999);
        let h = stitch(&c, &p, &s, &q("ACGT")).unwrap().unwrap();
        assert_eq!(h.seq, b"GTAC");
    }

    #[test]
    fn coverage_and_identity_are_computed_as_python_does() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1], vec!["ACGTACGTAC"]);
        let mut p = path(&[(Orient::Forward, 1)], 10, 0, 10);
        p.num_matches = 95;
        p.aln_len = 100;
        let h = stitch(&c, &p, &s, &q(&"N".repeat(20))).unwrap().unwrap();
        assert!((h.cov - 0.5).abs() < 1e-12, "10 bases over a 20 bp query");
        assert!((h.pident - 95.0).abs() < 1e-12);
        assert_eq!(h.nodes, 1);
        assert_eq!(h.id, "ACC_primary");
    }

    #[test]
    fn below_min_cov_yields_no_haplotype() {
        let c = cfg(4, 0.85);
        let s = sub(vec![1], vec!["ACGTACGTAC"]);
        let p = path(&[(Orient::Forward, 1)], 10, 0, 10);
        // 10 bases against a 100 bp query is 10% coverage.
        assert!(stitch(&c, &p, &s, &q(&"N".repeat(100))).unwrap().is_none());
    }

    /// GraphAligner naming a node absent from the GFA means the alignment cannot be
    /// realised. Python returns no haplotype rather than erroring; same here.
    #[test]
    fn unknown_node_in_path_yields_none() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1], vec!["ACGTACGTAC"]);
        let p = path(&[(Orient::Forward, 1), (Orient::Forward, 77)], 10, 0, 10);
        assert!(stitch(&c, &p, &s, &q("ACGT")).unwrap().is_none());
    }

    #[test]
    fn empty_path_yields_none() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1], vec!["ACGT"]);
        assert!(stitch(&c, &path(&[], 0, 0, 0), &s, &q("ACGT")).unwrap().is_none());
    }

    #[test]
    fn empty_query_is_an_error_not_a_division_by_zero() {
        let c = cfg(4, 0.0);
        let s = sub(vec![1], vec!["ACGT"]);
        let p = path(&[(Orient::Forward, 1)], 4, 0, 4);
        assert!(stitch(&c, &p, &s, &q("")).is_err());
    }
}

/// Full-chain comparison against the committed golden fixture, which was produced by the
/// Python walker. This is the end-to-end concordance check: fetch through stitch, in Rust,
/// reproducing a haplotype Python generated.
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
    fn rust_reproduces_the_python_golden_haplotype() {
        let Ok(dir) = std::env::var("LOGAN_TEST_DIR") else {
            eprintln!("SKIP: set LOGAN_TEST_DIR");
            return;
        };
        let Ok(compressed) = std::fs::read(format!("{dir}/ERR2028820.unitigs.fa.zst")) else {
            eprintln!("SKIP: no data");
            return;
        };
        let mut cfg = Config::default();
        cfg.scratch_dir = std::env::temp_dir().join(format!("logan_golden_{}", std::process::id()));
        if align::resolve_binary(&cfg).is_err() {
            eprintln!("SKIP: GraphAligner not on PATH");
            return;
        }

        let f = FetchedUnitigs {
            accession: "ERR2028820".into(),
            compressed_len: compressed.len(),
            compressed,
            fetch_secs: 0.0,
        };
        let graph = { let rd = stream::decoder(&f).unwrap(); parse::build(&f.accession, rd).unwrap() };
        let query = read_query(Path::new("../queries/pfcrt_only.fasta")).unwrap();
        let idx = anchor::index_query(&cfg, &query).unwrap();
        let a = anchor::find(&cfg, &graph, &idx).unwrap();
        let nodes = extract::expand(&cfg, &graph, &a).unwrap();
        let seqs = extract::collect_seqs(&graph, &nodes).unwrap();
        let subgraph = extract::write_gfa(&cfg, &graph, nodes, seqs).unwrap();
        let gaf = subgraph.gfa_path.with_extension("gaf");
        let gp = align::run(&cfg, &subgraph, &query, &gaf).unwrap().expect("a path");
        let h = stitch(&cfg, &gp, &subgraph, &query).unwrap().expect("a haplotype");

        let sha: String = Sha256::digest(&h.seq)
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect();
        eprintln!(
            "  {} len {} cov {:.6} pident {:.3} nodes {}",
            h.id, h.seq.len(), h.cov, h.pident, h.nodes
        );
        eprintln!("  sha256 {sha}");

        // Values recorded in bench/golden/ERR2028820.pfcrt.golden.json,
        // produced by the unmodified Python walker.
        assert_eq!(h.id, "ERR2028820_primary");
        assert_eq!(h.seq.len(), 1274, "haplotype length");
        assert!((h.cov - 0.999216).abs() < 1e-6, "coverage {}", h.cov);
        assert!((h.pident - 99.529).abs() < 1e-3, "identity {}", h.pident);
        assert_eq!(
            sha, "a83dd93100279b42b5a1812982b268ef9c2cb0a5970a61b23679f3d1e761c8e7",
            "haplotype does not match the golden fixture"
        );
    }
}
