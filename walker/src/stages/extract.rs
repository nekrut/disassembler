//! Step 4: 5-hop BFS around the anchors, emitted as GFA.
//!
//! Replaces `logan_path_walker.py:101-165` (`extract_local_subgraph`). The Python version
//! indexes the whole graph before consulting the anchors, then re-reads the entire file for
//! sequences. Here the CSR from `parse` is already in memory, so the BFS is a pointer chase
//! and only the sequence pass needs the bytes again.
//!
//! Scale, measured: 76 anchors expand to 356 nodes for `ERR2028820`; median 330 across the
//! ERG11 panel; 26,323 for `SRR8209082`, which is about a third of that accession's entire
//! graph. The issue's "20-100 local unitigs" is wrong by one to two orders of magnitude,
//! so do not size buffers from it.
//!
//! GFA shape:
//!
//! ```text
//! H	VN:Z:1.0
//! S	2933	CTGTATACACCC…
//! L	2933	-	443893	+	30M      <- overlap is always k-1
//! ```
//!
//! Write it to tmpfs. It is ~38 KB, and that is the whole of what Frontier 4 called
//! "intermediate disk I/O".
//!
//! IN   UnitigGraph, Anchors, decompressed unitig bytes
//! OUT  Subgraph (+ a GFA file on disk for GraphAligner)

use crate::config::Config;
use crate::error::{Error, Result};
use crate::stages::parse;
use std::io::{BufWriter, Write};
use crate::types::{Anchors, PathStep, Subgraph, UnitigGraph, UnitigId};

/// Breadth-first expansion from the anchors, `cfg.hops` times.
///
/// Replaces the BFS at `logan_path_walker.py:128-133`. There the adjacency is a dict of
/// sets built over the whole accession before the anchors are even consulted; here the CSR
/// from `parse::build` is already in memory, so this is a pointer chase with no allocation
/// beyond the frontier.
///
/// # Direction
///
/// Logan writes every edge into *both* endpoints' headers, so the CSR already contains each
/// connection from both sides and following `targets` forward is equivalent to Python's
/// explicitly symmetrised `adj`. Verified on `ERR2028820`: 2,578,430 distinct link tokens,
/// **zero** with a missing reciprocal. Also present are 1,156 self-loops, which the visited
/// set absorbs.
///
/// If a future accession did carry a one-sided edge this would quietly explore slightly
/// less than Python does. That has not been observed, and the concordance fixtures would
/// catch it on the two golden accessions.
///
/// # Scale
///
/// Measured: 77 anchors expand to 422 nodes on `ERR2028820`, out of 1,277,362. The issue's
/// "20-100 local unitigs" understates this by one to two orders of magnitude, and for
/// `SRR8209082` the neighbourhood reaches 8,466 nodes, about a tenth of that accession's
/// entire graph. Do not size buffers from the issue's figure.
pub fn expand(cfg: &Config, graph: &UnitigGraph, anchors: &Anchors) -> Result<Vec<UnitigId>> {
    if anchors.ids.is_empty() {
        return Ok(Vec::new());
    }

    // A bitset over unitig ids is the cheap membership test here: one bit each, so 156 KB
    // for a 1.28M-unitig accession, versus a hash set that would allocate per insert.
    let mut visited = vec![0u64; graph.n_unitigs.div_ceil(64)];
    let mark = |v: &mut Vec<u64>, id: u32| -> bool {
        let (w, b) = ((id >> 6) as usize, id & 63);
        let bit = 1u64 << b;
        if v[w] & bit != 0 {
            false
        } else {
            v[w] |= bit;
            true
        }
    };

    let mut all: Vec<UnitigId> = Vec::with_capacity(anchors.ids.len() * 8);
    let mut frontier: Vec<UnitigId> = Vec::with_capacity(anchors.ids.len());
    for &a in &anchors.ids {
        if (a as usize) >= graph.n_unitigs {
            return Err(Error::Parse {
                accession: graph.accession.clone(),
                detail: format!("anchor {a} is outside the graph ({} unitigs)", graph.n_unitigs),
            });
        }
        if mark(&mut visited, a) {
            all.push(a);
            frontier.push(a);
        }
    }

    // Only the newly reached nodes need expanding each round. Python re-scans every node it
    // has ever seen, which is the same closure but quadratic in the neighbourhood size.
    let mut next: Vec<UnitigId> = Vec::new();
    for _ in 0..cfg.hops {
        next.clear();
        for &node in &frontier {
            for &packed in graph.edges(node) {
                let (target, _, _) = parse::unpack_edge(packed);
                if (target as usize) >= graph.n_unitigs {
                    return Err(Error::Parse {
                        accession: graph.accession.clone(),
                        detail: format!("edge target {target} is outside the graph"),
                    });
                }
                if mark(&mut visited, target) {
                    all.push(target);
                    next.push(target);
                }
            }
        }
        if next.is_empty() {
            break; // the component is closed; further hops cannot add anything
        }
        std::mem::swap(&mut frontier, &mut next);
    }

    // Ascending order keeps GFA output deterministic regardless of discovery order.
    all.sort_unstable();
    Ok(all)
}

/// The gene's own sequencing depth, measured from an alignment that already succeeded.
///
/// Median `ka:f` over the unitigs the alignment uses. The alignment *is* the gene, so this
/// measures the right quantity, and it is precisely the estimate a prune running *before*
/// alignment cannot have: the anchor set's median on `SRR8209082` is 70x against a true
/// 14,500x, because `min_shared_kmers = 3` admits thousands of shallow incidental matches.
/// Every floor derived from the anchor set collapsed to `min_abund` and could not be swept.
///
/// Median rather than mean, so one repeat-collapsed unitig cannot drag it.
pub fn path_depth(graph: &UnitigGraph, steps: &[PathStep]) -> f32 {
    let mut ab: Vec<f32> = steps.iter().map(|s| graph.abundance[s.node as usize]).collect();
    if ab.is_empty() {
        return 0.0;
    }
    ab.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    ab[ab.len() / 2]
}

/// Mean `ka:f` over the unitigs the path traverses, for reporting.
///
/// Deliberately the arithmetic mean, not `path_depth`'s median: the field it fills is named
/// `mean_ka` and is published as `ka=` in the FASTA header, so it must be what it says it is.
/// The median is the right statistic for *deciding* a prune floor, where one repeat-collapsed
/// unitig must not drag the estimate, and `path_depth` stays that. Use this only to describe
/// a haplotype, never to gate one.
///
/// Unweighted per step, so a 50 bp unitig counts as much as a 50 kb one, and a unitig visited
/// twice counts twice.
pub fn path_mean_ka(graph: &UnitigGraph, steps: &[PathStep]) -> f64 {
    if steps.is_empty() {
        return 0.0;
    }
    let sum: f64 = steps.iter().map(|s| graph.abundance[s.node as usize] as f64).sum();
    sum / steps.len() as f64
}

/// The coverage a unitig must carry to survive a prune, given the gene's depth.
///
/// Never below `min_abund`: this can only tighten the existing floor, never loosen it.
pub fn depth_floor(cfg: &Config, depth: f32) -> f32 {
    if cfg.prune_frac > 0.0 {
        (cfg.prune_frac as f32 * depth).max(cfg.min_abund)
    } else {
        cfg.min_abund
    }
}

/// Keep only the unitigs carrying at least `floor` coverage.
///
/// Node-level, not edge-level: an error branch is a *unitig* that should not be in the
/// search space at all, and removing the unitig removes every way into it at once.
pub fn prune_nodes(
    graph: &UnitigGraph,
    nodes: Vec<UnitigId>,
    seqs: Vec<Vec<u8>>,
    floor: f32,
) -> (Vec<UnitigId>, Vec<Vec<u8>>) {
    let mut kn = Vec::with_capacity(nodes.len());
    let mut ks = Vec::with_capacity(seqs.len());
    for (id, seq) in nodes.into_iter().zip(seqs) {
        if graph.abundance[id as usize] >= floor {
            kn.push(id);
            ks.push(seq);
        }
    }
    (kn, ks)
}

/// Share of the whole graph pulled into the local neighbourhood.
///
/// This is the guard signal, and it replaces the anchor-fraction check that used to live in
/// `anchor::find`. That one was calibrated on BLAST anchor counts and then applied to k-mer
/// anchors, so at the default `min_shared_kmers = 3` it sat 9x below its threshold and
/// stayed silent on the very accession it existed to catch.
///
/// Subgraph fraction measures the thing that actually causes trouble: a neighbourhood that
/// swallows a large share of the graph offers many near-equal paths, leaving the haplotype
/// under-determined. Measured at the default threshold:
///
/// | accession | anchors | subgraph | fraction |
/// |---|---|---|---|
/// | `ERR2028820` | 77 | 422 | 0.033% |
/// | `SRR8209082` | 92 | 8,466 | 10.4% |
///
/// A separation of about 300x, against 1.2x for the anchor-fraction version.
pub fn subgraph_fraction(graph: &UnitigGraph, nodes: &[UnitigId]) -> f64 {
    if graph.n_unitigs == 0 {
        0.0
    } else {
        nodes.len() as f64 / graph.n_unitigs as f64
    }
}

/// Decode sequences for the selected nodes out of the graph's 2-bit store.
///
/// `O(nodes kept)`, not `O(all unitigs)`: `seq_offset` gives each unitig's slice directly, so
/// this touches 422 entries rather than rescanning 1,277,362. The Python version has to
/// re-read the entire unitig file here (`logan_path_walker.py:138`) because it kept no
/// sequence index and had already discarded the text.
///
/// Sequences come back expanded to one byte per base because that is what GFA output and
/// the stitching step both need. The packed form stays in the graph; this is a copy of the
/// small subset, typically a few hundred KB.
pub fn collect_seqs(graph: &UnitigGraph, nodes: &[UnitigId]) -> Result<Vec<Vec<u8>>> {
    let mut out = Vec::with_capacity(nodes.len());
    for &id in nodes {
        if (id as usize) >= graph.n_unitigs {
            return Err(Error::Parse {
                accession: graph.accession.clone(),
                detail: format!(
                    "node {id} is outside the graph ({} unitigs)",
                    graph.n_unitigs
                ),
            });
        }
        let mut seq = Vec::new();
        graph.decode_seq(id, &mut seq);
        if seq.is_empty() {
            return Err(Error::Parse {
                accession: graph.accession.clone(),
                detail: format!("node {id} has an empty sequence"),
            });
        }
        out.push(seq);
    }
    Ok(out)
}

/// A scratch basename no other worker or process can collide with.
///
/// The GFA, the GAF derived from it, and the single-record query FASTA derived from that all
/// hang off this stem. Keying them on the accession alone meant two `logan-walker` processes,
/// or one process with the same accession twice in its list, truncated and unlinked each
/// other's files: measured as GraphAligner SIGABRTs, a "non supported overlap" on a
/// half-written GFA, and a run that produced zero records at exit 0.
///
/// `--scratch-dir` defaults to `/dev/shm/logan-walker` for *every* process, so two runs over
/// different genes collide by default. The PID separates processes and the counter separates
/// duplicate accessions inside one process.
fn scratch_stem(accession: &str) -> String {
    use std::sync::atomic::{AtomicU64, Ordering};
    static SEQ: AtomicU64 = AtomicU64::new(0);
    let n = SEQ.fetch_add(1, Ordering::Relaxed);
    format!("{accession}.{}.{n}", std::process::id())
}

/// Serialise nodes and edges to GFA and return the populated `Subgraph`.
///
/// Byte-for-byte compatible with `logan_path_walker.py:153-165`, because GraphAligner's
/// output depends on its input and the golden fixtures pin that output. Specifically:
///
/// * `H	VN:Z:1.0` header
/// * one `S` line per node, ascending id, sequence expanded
/// * `L` lines emitted per source node in ascending id order, and within a node in the order
///   the links appeared in that unitig's header
/// * only edges whose *both* endpoints are in the subgraph
/// * overlap always `{k-1}M`, which is exact for a de Bruijn graph of order k
///
/// Each edge therefore appears twice, once from each endpoint, because Logan records it in
/// both headers. Python does the same, so the duplication is preserved deliberately.
///
/// The file goes to `cfg.scratch_dir`, which should be tmpfs: it is ~38 KB for a typical
/// accession and GraphAligner is the only consumer.
pub fn write_gfa(
    cfg: &Config,
    graph: &UnitigGraph,
    nodes: Vec<UnitigId>,
    seqs: Vec<Vec<u8>>,
) -> Result<Subgraph> {
    if nodes.len() != seqs.len() {
        return Err(Error::Parse {
            accession: graph.accession.clone(),
            detail: format!("{} nodes but {} sequences", nodes.len(), seqs.len()),
        });
    }

    std::fs::create_dir_all(&cfg.scratch_dir)?;
    let gfa_path = cfg.scratch_dir.join(format!("{}.sub.gfa", scratch_stem(&graph.accession)));
    let file = std::fs::File::create(&gfa_path)?;
    let mut w = BufWriter::with_capacity(1 << 16, file);

    // Membership test for "is this edge target inside the subgraph". `nodes` is sorted, so a
    // bitset keyed by unitig id is both smaller and faster than a binary search per edge.
    let mut inside = vec![0u64; graph.n_unitigs.div_ceil(64)];
    for &id in &nodes {
        inside[(id >> 6) as usize] |= 1u64 << (id & 63);
    }
    // An edge target may point past the graph in a malformed (non-reciprocal) input. Valid
    // Logan output never does, but every other function here turns that into a clean skip or
    // error rather than an out-of-bounds panic, so guard the range: an out-of-range target
    // cannot be one of `nodes` anyway, so it is simply not inside the subgraph.
    let is_inside =
        |id: u32| (id as usize) < graph.n_unitigs && inside[(id >> 6) as usize] & (1u64 << (id & 63)) != 0;

    w.write_all(b"H\tVN:Z:1.0\n")?;
    for (i, &id) in nodes.iter().enumerate() {
        w.write_all(b"S\t")?;
        write!(w, "{id}")?;
        w.write_all(b"\t")?;
        w.write_all(&seqs[i])?;
        w.write_all(b"\n")?;
    }

    let overlap = cfg.graph_k.saturating_sub(1);
    let mut n_edges = 0usize;
    for &id in &nodes {
        for &packed in graph.edges(id) {
            let (target, src_fwd, dst_fwd) = parse::unpack_edge(packed);
            if !is_inside(target) {
                continue;
            }
            writeln!(
                w,
                "L\t{id}\t{}\t{target}\t{}\t{overlap}M",
                if src_fwd { '+' } else { '-' },
                if dst_fwd { '+' } else { '-' },
            )?;
            n_edges += 1;
        }
    }
    w.flush()?;

    Ok(Subgraph {
        accession: graph.accession.clone(),
        nodes,
        seqs,
        gfa_path,
        n_edges,
    })
}

#[cfg(test)]
mod collect_tests {
    use super::*;
    use crate::stages::parse;

    fn graph(fa: &str) -> UnitigGraph {
        parse::build(&"T".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec()))).unwrap()
    }

    #[test]
    fn decodes_the_requested_nodes_in_order() {
        let g = graph(">T_0\nACGT\n>T_1\nTTGCA\n>T_2\nG\n>T_3\nCCCCCCCC\n");
        let got = collect_seqs(&g, &[0, 2, 3]).unwrap();
        assert_eq!(got.len(), 3);
        assert_eq!(got[0], b"ACGT");
        assert_eq!(got[1], b"G");
        assert_eq!(got[2], b"CCCCCCCC");
    }

    /// The 2-bit store packs across unitig boundaries, so most of these begin mid-byte.
    /// This is the case a byte-aligned implementation would pass by accident.
    #[test]
    fn handles_unaligned_starts() {
        let g = graph(">T_0\nA\n>T_1\nCC\n>T_2\nGGG\n>T_3\nTTTT\n>T_4\nACGTA\n");
        let got = collect_seqs(&g, &[0, 1, 2, 3, 4]).unwrap();
        assert_eq!(got, vec![
            b"A".to_vec(), b"CC".to_vec(), b"GGG".to_vec(),
            b"TTTT".to_vec(), b"ACGTA".to_vec(),
        ]);
    }

    #[test]
    fn empty_node_list_is_fine() {
        let g = graph(">T_0\nACGT\n");
        assert!(collect_seqs(&g, &[]).unwrap().is_empty());
    }

    #[test]
    fn out_of_range_node_is_rejected() {
        let g = graph(">T_0\nACGT\n");
        assert!(collect_seqs(&g, &[0, 5]).is_err());
    }

    #[test]
    fn repeated_nodes_are_decoded_each_time() {
        // collect_seqs does not deduplicate; expand already guarantees uniqueness.
        let g = graph(">T_0\nACGT\n");
        let got = collect_seqs(&g, &[0, 0]).unwrap();
        assert_eq!(got.len(), 2);
        assert_eq!(got[0], got[1]);
    }

    #[test]
    fn agrees_with_base_at_for_every_position() {
        let g = graph(">T_0\nACGTACGTA\n>T_1\nTTGCAGGCTA\n");
        const BASES: [u8; 4] = *b"ACGT";
        for id in 0..2u32 {
            let seq = &collect_seqs(&g, &[id]).unwrap()[0];
            for (j, &b) in seq.iter().enumerate() {
                assert_eq!(b, BASES[g.base_at(id, j) as usize], "unitig {id} base {j}");
            }
        }
    }
}

#[cfg(test)]
mod expand_tests {
    use super::*;
    use crate::stages::parse;
    use crate::types::Anchors;

    fn graph(fa: &str) -> UnitigGraph {
        parse::build(&"T".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec()))).unwrap()
    }
    fn anchors(ids: &[UnitigId], n: usize) -> Anchors {
        Anchors {
            ids: ids.to_vec(),
            shared_kmers: vec![1; ids.len()],
            fraction_of_graph: ids.len() as f64 / n as f64,
        }
    }
    fn cfg_hops(h: usize) -> Config {
        Config { hops: h, ..Config::default() }
    }

    /// A path 0-1-2-3-4-5. Starting at 0, each hop should reach exactly one further node.
    fn chain() -> UnitigGraph {
        graph(concat!(
            ">T_0 L:+:1:+\nAAAA\n",
            ">T_1 L:+:0:+ L:+:2:+\nCCCC\n",
            ">T_2 L:+:1:+ L:+:3:+\nGGGG\n",
            ">T_3 L:+:2:+ L:+:4:+\nTTTT\n",
            ">T_4 L:+:3:+ L:+:5:+\nAACC\n",
            ">T_5 L:+:4:+\nGGTT\n",
        ))
    }

    #[test]
    fn each_hop_advances_exactly_one_step() {
        let g = chain();
        for (hops, want) in [(0usize, 1usize), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)] {
            let got = expand(&cfg_hops(hops), &g, &anchors(&[0], 6)).unwrap();
            assert_eq!(got.len(), want, "hops={hops}");
        }
    }

    #[test]
    fn stops_early_when_the_component_is_closed() {
        // 20 hops cannot exceed the 6 nodes that exist.
        let g = chain();
        assert_eq!(expand(&cfg_hops(20), &g, &anchors(&[0], 6)).unwrap().len(), 6);
    }

    #[test]
    fn does_not_cross_into_a_disconnected_component() {
        let g = graph(concat!(
            ">T_0 L:+:1:+\nAAAA\n",
            ">T_1 L:+:0:+\nCCCC\n",
            ">T_2 L:+:3:+\nGGGG\n",
            ">T_3 L:+:2:+\nTTTT\n",
        ));
        assert_eq!(expand(&cfg_hops(5), &g, &anchors(&[0], 4)).unwrap(), vec![0, 1]);
    }

    /// Logan emits self-loops (1,156 of them in ERR2028820). They must not loop forever.
    #[test]
    fn self_loops_terminate() {
        let g = graph(">T_0 L:+:0:+ L:+:1:+\nAAAA\n>T_1 L:+:0:+\nCCCC\n");
        assert_eq!(expand(&cfg_hops(5), &g, &anchors(&[0], 2)).unwrap(), vec![0, 1]);
    }

    #[test]
    fn output_is_sorted_and_deduplicated() {
        let g = chain();
        // Overlapping anchors, given out of order.
        let got = expand(&cfg_hops(1), &g, &anchors(&[3, 0, 3], 6)).unwrap();
        assert_eq!(got, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn no_anchors_gives_an_empty_neighbourhood() {
        assert!(expand(&cfg_hops(5), &chain(), &anchors(&[], 6)).unwrap().is_empty());
    }

    #[test]
    fn anchor_outside_the_graph_is_rejected() {
        assert!(expand(&cfg_hops(1), &chain(), &anchors(&[99], 6)).is_err());
    }

    #[test]
    fn fraction_is_computed_over_the_whole_graph() {
        let g = chain();
        let nodes = expand(&cfg_hops(1), &g, &anchors(&[0], 6)).unwrap();
        assert!((subgraph_fraction(&g, &nodes) - 2.0 / 6.0).abs() < 1e-9);
    }
}

/// Cross-validation against the Python anchor-recall study's recorded subgraph sizes.
///
/// `LOGAN_TEST_DIR=/dir/with/*.unitigs.fa.zst cargo test --release -- --ignored --nocapture`
#[cfg(test)]
mod python_concordance_expand {
    use super::*;
    use crate::io::input::read_query;
    use crate::stages::{anchor, parse, stream};
    use crate::types::FetchedUnitigs;
    use std::path::Path;

    fn load(dir: &str, acc: &str) -> Option<UnitigGraph> {
        let p = format!("{dir}/{acc}.unitigs.fa.zst");
        let compressed = std::fs::read(&p).ok()?;
        let f = FetchedUnitigs { index: 0, accession: acc.into(), compressed_len: compressed.len(), compressed, fetch_secs: 0.0 };
        let g = {
            let rd = stream::decoder(&f).ok()?;
            parse::build(&f.accession, rd).ok()
        };
        g
    }

    #[test]
    #[ignore]
    fn matches_python_subgraph_sizes() {
        let Ok(dir) = std::env::var("LOGAN_TEST_DIR") else {
            eprintln!("SKIP: set LOGAN_TEST_DIR");
            return;
        };
        // (accession, query, [(min_shared_kmers, expected anchors, expected subgraph nodes)])
        let cases: [(&str, &str, &[(u32, usize, usize)]); 2] = [
            ("ERR2028820", "../queries/pfcrt_only.fasta",
             &[(1, 104, 452), (2, 83, 422), (3, 77, 422), (5, 59, 386), (10, 36, 348)]),
            ("SRR8209082", "../queries/query_hiv1_pol.fasta",
             &[(1, 1285, 24331), (3, 92, 8466)]),
        ];

        for (acc, qpath, expected) in cases {
            let Some(graph) = load(&dir, acc) else {
                eprintln!("SKIP: no data for {acc}");
                continue;
            };
            let query = read_query(Path::new(qpath)).unwrap();
            eprintln!("{acc} ({} unitigs)", graph.n_unitigs);
            for &(m, want_anchors, want_nodes) in expected {
                let cfg = Config { min_shared_kmers: m, ..Config::default() };
                let idx = anchor::index_query(&cfg, &query).unwrap();
                let a = anchor::find(&cfg, &graph, &idx).unwrap();
                let t0 = std::time::Instant::now();
                let nodes = expand(&cfg, &graph, &a).unwrap();
                eprintln!(
                    "  M={:<3} anchors {:>5} (py {:>5})  subgraph {:>6} (py {:>6})  {:>7.4}s  frac {:>7.3}%",
                    m, a.ids.len(), want_anchors, nodes.len(), want_nodes,
                    t0.elapsed().as_secs_f64(), subgraph_fraction(&graph, &nodes) * 100.0
                );
                assert_eq!(a.ids.len(), want_anchors, "{acc} anchors at M={m}");
                assert_eq!(nodes.len(), want_nodes, "{acc} subgraph at M={m}");
                assert!(nodes.windows(2).all(|w| w[0] < w[1]), "must be sorted and unique");
            }
        }
    }
}

/// True when the local neighbourhood covers so much of the graph that the eventual haplotype
/// is likely under-determined. See `Config::subgraph_fraction_guard` for the calibration.
pub fn exceeds_guard(cfg: &Config, graph: &UnitigGraph, nodes: &[UnitigId]) -> bool {
    subgraph_fraction(graph, nodes) > cfg.subgraph_fraction_guard
}

/// Validation of the 2-bit sequence store against the original text.
///
/// ```text
/// LOGAN_TEST_ZST=.../ERR2028820.unitigs.fa.zst \
/// LOGAN_TEST_PLAIN=.../ERR2028820.unitigs.fa \
///     cargo test --release -- --ignored --nocapture
/// ```
#[cfg(test)]
mod real_data_collect {
    use super::*;
    use crate::io::input::read_query;
    use crate::stages::{anchor, parse, stream};
    use crate::types::FetchedUnitigs;
    use std::io::{BufRead, BufReader};
    use std::path::Path;

    /// Every unitig, compared against the plaintext record by record. Ids are sequential, so
    /// this streams both sides and holds only one record at a time.
    #[test]
    #[ignore]
    fn every_sequence_round_trips_against_the_original_text() {
        let (Ok(zst), Ok(plain)) = (
            std::env::var("LOGAN_TEST_ZST"),
            std::env::var("LOGAN_TEST_PLAIN"),
        ) else {
            eprintln!("SKIP: set LOGAN_TEST_ZST and LOGAN_TEST_PLAIN");
            return;
        };
        let Ok(compressed) = std::fs::read(&zst) else {
            eprintln!("SKIP: cannot read {zst}");
            return;
        };
        let f = FetchedUnitigs {
            index: 0,
            accession: "ERR2028820".into(),
            compressed_len: compressed.len(),
            compressed,
            fetch_secs: 0.0,
        };
        let graph = {
            let rd = stream::decoder(&f).unwrap();
            parse::build(&f.accession, rd).unwrap()
        };

        let Ok(fh) = std::fs::File::open(&plain) else {
            eprintln!("SKIP: cannot read {plain}");
            return;
        };
        let mut rd = BufReader::with_capacity(1 << 20, fh);
        let mut line = String::new();
        let mut id: u32 = 0;
        let mut decoded = Vec::new();
        let mut checked = 0usize;

        loop {
            line.clear();
            if rd.read_line(&mut line).unwrap() == 0 {
                break;
            }
            let t = line.trim_end();
            if t.starts_with('>') || t.is_empty() {
                continue;
            }
            graph.decode_seq(id, &mut decoded);
            assert_eq!(
                decoded.as_slice(),
                t.as_bytes(),
                "unitig {id} does not match the original text"
            );
            id += 1;
            checked += 1;
        }
        eprintln!("  round-tripped {checked} sequences byte-for-byte");
        assert_eq!(checked, graph.n_unitigs, "record count mismatch");

        // And the subgraph path specifically.
        let query = read_query(Path::new("../queries/pfcrt_only.fasta")).unwrap();
        let cfg = Config::default();
        let idx = anchor::index_query(&cfg, &query).unwrap();
        let a = anchor::find(&cfg, &graph, &idx).unwrap();
        let nodes = expand(&cfg, &graph, &a).unwrap();
        let seqs = collect_seqs(&graph, &nodes).unwrap();
        assert_eq!(seqs.len(), nodes.len());
        assert_eq!(nodes.len(), 422);
        let total: usize = seqs.iter().map(|s| s.len()).sum();
        eprintln!(
            "  subgraph {} nodes, {} bp, mean {:.1} bp",
            nodes.len(),
            total,
            total as f64 / nodes.len() as f64
        );
        for (i, s) in seqs.iter().enumerate() {
            assert_eq!(s.len(), graph.seq_len(nodes[i]) as usize);
            assert!(s.iter().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')));
        }
    }
}

#[cfg(test)]
mod gfa_tests {
    use super::*;
    use crate::stages::parse;
    use std::io::Read;

    fn graph(fa: &str) -> UnitigGraph {
        parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec()))).unwrap()
    }
    fn tmp_cfg(tag: &str) -> Config {
        let mut d = std::env::temp_dir();
        d.push(format!("logan_gfa_test_{}_{}", std::process::id(), tag));
        Config { scratch_dir: d, ..Config::default() }
    }
    fn render(cfg: &Config, g: &UnitigGraph, nodes: Vec<UnitigId>) -> String {
        let seqs = collect_seqs(g, &nodes).unwrap();
        let sub = write_gfa(cfg, g, nodes, seqs).unwrap();
        let mut s = String::new();
        std::fs::File::open(&sub.gfa_path).unwrap().read_to_string(&mut s).unwrap();
        std::fs::remove_file(&sub.gfa_path).ok();
        s
    }

    /// The overlap in an L line comes from `graph_k`, not `k`. Those were one field until
    /// upstream's `3aabf8f` split them, and this test kept setting `k`, which the GFA writer
    /// does not read. It went on passing only because `graph_k` also defaults to 31, so it
    /// asserted `30M` without controlling it: had the default moved, the test would have
    /// failed while pointing at the wrong knob.
    #[test]
    fn matches_the_python_line_format_exactly() {
        let g = graph(">ACC_0 ka:f:5.3   L:-:1:+  L:+:1:-\nACGT\n>ACC_1 ka:f:9.0   L:+:0:+\nTTGCA\n");
        let cfg = Config { graph_k: 31, ..tmp_cfg("fmt") };
        let got = render(&cfg, &g, vec![0, 1]);
        assert_eq!(
            got,
            "H\tVN:Z:1.0\n\
             S\t0\tACGT\n\
             S\t1\tTTGCA\n\
             L\t0\t-\t1\t+\t30M\n\
             L\t0\t+\t1\t-\t30M\n\
             L\t1\t+\t0\t+\t30M\n"
        );
    }

    /// Python writes an L line only when both endpoints survived into the subgraph
    /// (`if tgt in seqs`). Keeping a dangling edge would make GraphAligner reference a node
    /// the file never declares.
    #[test]
    fn drops_edges_leaving_the_subgraph() {
        let g = graph(
            ">ACC_0 L:+:1:+ L:+:2:+\nACGT\n>ACC_1 L:+:0:+\nTTTT\n>ACC_2 L:+:0:+\nGGGG\n",
        );
        let cfg = tmp_cfg("dangling");
        let got = render(&cfg, &g, vec![0, 1]);
        assert!(got.contains("L\t0\t+\t1\t+\t30M\n"), "kept edge missing:\n{got}");
        assert!(!got.contains("\t2\t"), "edge to excluded node 2 leaked:\n{got}");
        assert_eq!(got.matches("S\t").count(), 2);
    }

    #[test]
    fn overlap_follows_k() {
        let g = graph(">ACC_0 L:+:1:+\nACGT\n>ACC_1 L:+:0:+\nTTTT\n");
        for (k, want) in [(31usize, "30M"), (21, "20M")] {
            let cfg = Config { graph_k: k, ..tmp_cfg(&format!("k{k}")) };
            let got = render(&cfg, &g, vec![0, 1]);
            assert!(got.contains(want), "k={k} expected {want} in:\n{got}");
        }
    }

    #[test]
    fn reports_edge_count_and_path() {
        let g = graph(">ACC_0 L:+:1:+\nACGT\n>ACC_1 L:+:0:+\nTTTT\n");
        let cfg = tmp_cfg("count");
        let seqs = collect_seqs(&g, &[0, 1]).unwrap();
        let sub = write_gfa(&cfg, &g, vec![0, 1], seqs).unwrap();
        assert_eq!(sub.n_edges, 2, "one connection, recorded from both ends");
        assert_eq!(sub.nodes.len(), 2);
        // The name carries the accession plus a PID and a per-task counter, so two runs
        // cannot write the same file; see `scratch_stem`.
        let name = sub.gfa_path.file_name().unwrap().to_string_lossy().into_owned();
        assert!(name.starts_with("ACC."), "accession stays legible: {name}");
        assert!(name.ends_with(".sub.gfa"), "suffix is unchanged: {name}");
        assert!(sub.gfa_path.exists());
        std::fs::remove_file(&sub.gfa_path).ok();
    }

    #[test]
    fn empty_subgraph_writes_a_header_only_file() {
        // GraphAligner aborts with SIGABRT on a nodeless graph, so align must refuse to run
        // rather than be handed this. Writing it is still valid GFA.
        let g = graph(">ACC_0\nACGT\n");
        let cfg = tmp_cfg("empty");
        assert_eq!(render(&cfg, &g, vec![]), "H\tVN:Z:1.0\n");
    }

    #[test]
    fn mismatched_nodes_and_seqs_is_rejected() {
        let g = graph(">ACC_0\nACGT\n>ACC_1\nTTTT\n");
        let cfg = tmp_cfg("mismatch");
        assert!(write_gfa(&cfg, &g, vec![0, 1], vec![b"ACGT".to_vec()]).is_err());
    }
}

/// Emit a real GFA plus its anchor list so the Python implementation can be run on exactly
/// the same input and the two files diffed.
///
/// ```text
/// LOGAN_TEST_ZST=.../ERR2028820.unitigs.fa.zst LOGAN_GFA_OUT=/tmp/cmp \
///     cargo test --release -- --ignored --nocapture emits_gfa_for_cross_check
/// ```
#[cfg(test)]
mod cross_check {
    use super::*;
    use crate::io::input::read_query;
    use crate::stages::{anchor, parse, stream};
    use crate::types::FetchedUnitigs;
    use std::path::{Path, PathBuf};

    #[test]
    #[ignore]
    fn emits_gfa_for_cross_check() {
        let (Ok(zst), Ok(outdir)) = (
            std::env::var("LOGAN_TEST_ZST"),
            std::env::var("LOGAN_GFA_OUT"),
        ) else {
            eprintln!("SKIP: set LOGAN_TEST_ZST and LOGAN_GFA_OUT");
            return;
        };
        let Ok(compressed) = std::fs::read(&zst) else {
            eprintln!("SKIP: cannot read {zst}");
            return;
        };
        let f = FetchedUnitigs {
            index: 0,
            accession: "ERR2028820".into(),
            compressed_len: compressed.len(),
            compressed,
            fetch_secs: 0.0,
        };
        let graph = {
            let rd = stream::decoder(&f).unwrap();
            parse::build(&f.accession, rd).unwrap()
        };
        let query = read_query(Path::new("../queries/pfcrt_only.fasta")).unwrap();
        let cfg = Config {
            scratch_dir: PathBuf::from(&outdir),
            ..Config::default()
        };
        let idx = anchor::index_query(&cfg, &query).unwrap();
        let a = anchor::find(&cfg, &graph, &idx).unwrap();
        let nodes = expand(&cfg, &graph, &a).unwrap();
        let seqs = collect_seqs(&graph, &nodes).unwrap();
        let sub = write_gfa(&cfg, &graph, nodes, seqs).unwrap();

        // The exact anchor ids, so Python can be driven with the same input.
        let list: Vec<String> = a.ids.iter().map(|i| i.to_string()).collect();
        std::fs::write(format!("{outdir}/anchors.txt"), list.join("\n")).unwrap();
        eprintln!(
            "  wrote {} ({} nodes, {} edges) and {} anchors",
            sub.gfa_path.display(),
            sub.nodes.len(),
            sub.n_edges,
            a.ids.len()
        );
    }
}

/// The depth estimate and floor that the second-pass prune is built on.
#[cfg(test)]
mod depth_tests {
    use super::*;
    use crate::stages::parse;
    use crate::types::{Orient, PathStep};

    fn graph(fa: &str) -> UnitigGraph {
        parse::build(&"T".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec()))).unwrap()
    }
    fn steps(ids: &[UnitigId]) -> Vec<PathStep> {
        ids.iter().map(|&node| PathStep { node, orient: Orient::Forward }).collect()
    }

    #[test]
    fn depth_is_the_median_along_the_alignment() {
        let g = graph(concat!(
            ">T_0 ka:f:100.0\nACGT\n",
            ">T_1 ka:f:200.0\nACGT\n",
            ">T_2 ka:f:300.0\nACGT\n",
        ));
        assert_eq!(path_depth(&g, &steps(&[0, 1, 2])), 200.0);
    }

    /// One repeat-collapsed unitig must not drag the estimate, which is why this is a
    /// median and not a mean.
    #[test]
    fn one_huge_unitig_does_not_move_the_estimate() {
        let g = graph(concat!(
            ">T_0 ka:f:100.0\nACGT\n",
            ">T_1 ka:f:100.0\nACGT\n",
            ">T_2 ka:f:25000.0\nACGT\n",
        ));
        assert_eq!(path_depth(&g, &steps(&[0, 1, 2])), 100.0, "mean would be 8400");
    }

    #[test]
    fn an_empty_alignment_has_no_depth() {
        let g = graph(">T_0 ka:f:100.0\nACGT\n");
        assert_eq!(path_depth(&g, &[]), 0.0);
    }

    /// Pruning is on by default at 1% of depth; `0.0` is how it is switched off, and that
    /// path must stay an exact no-op.
    #[test]
    fn the_floor_is_off_when_prune_frac_is_zero() {
        assert_eq!(Config::default().prune_frac, 0.01, "the panel-calibrated default");
        let off = Config { prune_frac: 0.0, ..Config::default() };
        assert_eq!(depth_floor(&off, 14500.0), off.min_abund);
        assert_eq!(depth_floor(&Config::default(), 14500.0), 145.0);
    }

    /// The whole point: five reads means something different at 180x and at 14,500x.
    #[test]
    fn the_floor_scales_with_depth() {
        let c = Config { prune_frac: 0.01, ..Config::default() };
        assert_eq!(depth_floor(&c, 14500.0), 145.0);
        assert_eq!(depth_floor(&c, 180.0), 5.0, "1% of 180x is under min_abund, which holds");
    }

    #[test]
    fn the_floor_never_drops_below_min_abund() {
        let c = Config { prune_frac: 0.9, min_abund: 5.0, ..Config::default() };
        for d in [0.0, 1.0, 5.0, 5.5, 1000.0] {
            assert!(depth_floor(&c, d) >= 5.0, "floor dipped at depth {d}");
        }
    }

    /// A clonal sample sits above its own floor everywhere, so nothing is removed.
    #[test]
    fn a_uniformly_covered_gene_loses_nothing() {
        let g = graph(concat!(
            ">T_0 ka:f:180.0\nACGT\n",
            ">T_1 ka:f:170.0\nACGT\n",
            ">T_2 ka:f:190.0\nACGT\n",
        ));
        let c = Config { prune_frac: 0.01, ..Config::default() };
        let floor = depth_floor(&c, path_depth(&g, &steps(&[0, 1, 2])));
        let (kept, seqs) = prune_nodes(&g, vec![0, 1, 2], vec![vec![b'A']; 3], floor);
        assert_eq!(kept, vec![0, 1, 2]);
        assert_eq!(seqs.len(), 3);
    }

    /// The shape the prune exists for: a deep gene with an error branch hanging off it.
    #[test]
    fn an_error_branch_beside_a_deep_gene_is_removed() {
        let g = graph(concat!(
            ">T_0 ka:f:14000.0\nACGT\n",
            ">T_1 ka:f:15000.0\nACGT\n",
            ">T_2 ka:f:2.0\nACGT\n",
        ));
        let c = Config { prune_frac: 0.01, ..Config::default() };
        // median of an even-length list takes the upper middle, so depth is 15000 here
        let floor = depth_floor(&c, path_depth(&g, &steps(&[0, 1])));
        assert_eq!(floor, 150.0);
        let (kept, _) = prune_nodes(&g, vec![0, 1, 2], vec![vec![b'A']; 3], floor);
        assert_eq!(kept, vec![0, 1], "the 2x branch is 0.014% of depth and goes");
    }

    /// Sequences must stay aligned with the ids they belong to.
    #[test]
    fn sequences_follow_their_unitigs() {
        let g = graph(concat!(
            ">T_0 ka:f:1000.0\nAAAA\n",
            ">T_1 ka:f:1.0\nCCCC\n",
            ">T_2 ka:f:1000.0\nGGGG\n",
        ));
        let (kept, seqs) = prune_nodes(
            &g, vec![0, 1, 2],
            vec![b"AAAA".to_vec(), b"CCCC".to_vec(), b"GGGG".to_vec()], 10.0,
        );
        assert_eq!(kept, vec![0, 2]);
        assert_eq!(seqs, vec![b"AAAA".to_vec(), b"GGGG".to_vec()]);
    }
}

#[cfg(test)]
mod scratch_stem_tests {
    use super::scratch_stem;

    /// The bug: two `logan-walker` processes, or one list with a duplicate accession, wrote
    /// the same `<acc>.sub.gfa` and unlinked each other's files mid-alignment.
    #[test]
    fn two_stems_for_one_accession_differ() {
        let a = scratch_stem("SRR8209082");
        let b = scratch_stem("SRR8209082");
        assert_ne!(a, b, "duplicate accessions must not share a scratch path");
        for s in [&a, &b] {
            assert!(s.starts_with("SRR8209082."), "accession stays legible in {s:?}");
            assert!(
                s.contains(&std::process::id().to_string()),
                "PID separates concurrent processes in {s:?}"
            );
        }
    }

    /// The GAF and the single-record query FASTA are derived from the GFA path, so
    /// uniqueness has to survive `with_extension`.
    #[test]
    fn derived_scratch_paths_stay_unique() {
        let p = |acc: &str| std::path::PathBuf::from(format!("/tmp/{}.sub.gfa", scratch_stem(acc)));
        let (g1, g2) = (p("ERR2028820"), p("ERR2028820"));
        assert_ne!(g1.with_extension("gaf"), g2.with_extension("gaf"));
        assert_ne!(g1.with_extension("query.fa"), g2.with_extension("query.fa"));
        assert_eq!(
            g1.with_extension("gaf").extension().unwrap(),
            "gaf",
            "the stem's dots must not swallow the extension"
        );
    }
}
