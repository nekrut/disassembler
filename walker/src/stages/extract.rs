//! Step 4: 5-hop BFS around the anchors, emitted as GFA.
//!
//! The CSR from `parse` is already in memory, so the BFS is a pointer chase and only the
//! sequence pass needs the bytes again.
//!
//! Scale, measured: 76 anchors expand to 356 nodes for `ERR2028820`; median 330 across the
//! ERG11 panel; 26,323 for `SRR8209082`, which is about a third of that accession's entire
//! graph. A "20-100 local unitigs" estimate is wrong by one to two orders of magnitude,
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
//! Write it to tmpfs. It is ~38 KB, and that is the whole of the intermediate disk I/O.
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
/// The CSR from `parse::build` is already in memory, so this is a pointer chase with no
/// allocation beyond the frontier.
///
/// # Direction
///
/// Logan writes every edge into *both* endpoints' headers, so the CSR already contains each
/// connection from both sides and following `targets` forward is equivalent to an
/// explicitly symmetrised adjacency. Verified on `ERR2028820`: 2,578,430 distinct link tokens,
/// **zero** with a missing reciprocal. Also present are 1,156 self-loops, which the visited
/// set absorbs.
///
/// If a future accession did carry a one-sided edge this would quietly explore slightly
/// less than a symmetrised adjacency would. That has not been observed, and the concordance
/// fixtures would catch it on the two golden accessions.
///
/// # Scale
///
/// Measured: 77 anchors expand to 422 nodes on `ERR2028820`, out of 1,277,362. An estimate
/// of "20-100 local unitigs" understates this by one to two orders of magnitude, and for
/// `SRR8209082` the neighbourhood reaches 8,466 nodes, about a tenth of that accession's
/// entire graph. Do not size buffers from that figure.
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

    // Only the newly reached nodes need expanding each round. Re-scanning every node ever
    // seen gives the same closure but is quadratic in the neighbourhood size.
    let mut next: Vec<UnitigId> = Vec::new();
    // One buffer for the whole walk: the rows may be on disk, and a fresh Vec per node would
    // allocate once per visited unitig.
    let mut row: Vec<u32> = Vec::new();
    for _ in 0..cfg.hops {
        next.clear();
        for &node in &frontier {
            graph.edges_by_id_into(node, &mut row).map_err(|e| Error::Parse {
                accession: graph.accession.clone(),
                detail: format!("reading adjacency for unitig {node}: {e}"),
            })?;
            for &packed in &row {
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
///
/// **Unweighted, deliberately.** WP-01 asks for a k-mer-weighted median, on the reasoning that
/// a route is not a list of equal things and short error unitigs can drag the estimate down.
/// That was implemented and **measured, and it is wrong on this walker's data**: on
/// `SRR8209082`, whose true depth is about 14,500x (see above), the step median gives 15,000x
/// and the k-mer-weighted median gives **91.7x**, which collapses the floor to `--min-abund`
/// and stops the prune removing anything. Forks went 136 -> 1,303, the primary changed, and
/// its identity to the reference rose from 97.48% to 98.50%, which on a less-pruned graph is
/// the reference pull the prune exists to remove rather than a better answer.
///
/// The reason is that a tangled neighbourhood's route threads long, low-coverage unitigs that
/// are not the gene, while the gene itself is cut into short high-coverage pieces by its own
/// variation. Weighting by sequence therefore weights the noise. Reconsider only with a
/// measurement on an accession whose depth is independently known.
pub fn path_depth(cfg: &Config, graph: &UnitigGraph, steps: &[PathStep]) -> f32 {
    let _ = cfg;
    median_abundance(steps.iter().map(|s| graph.abundance_of(s.node)).collect())
}

/// Median of a set of abundances, or 0.0 when empty.
///
/// WV-23. `variants` measures the same quantity over the oriented walk rather than over
/// `PathStep`s, and had its own copy of this. Identical today, but a change to one (a NaN-safe
/// sort, length weighting) would have made the prune floor and the allele floor measure gene
/// depth differently, which is the kind of divergence nothing downstream would reveal.
pub fn median_abundance(mut ab: Vec<f32>) -> f32 {
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
    let sum: f64 = steps.iter().map(|s| graph.abundance_of(s.node) as f64).sum();
    sum / steps.len() as f64
}

/// The coverage a unitig must carry to survive a prune, given the gene's depth.
///
/// `prune_frac` of the gene's depth, never below `prune_min_abund`, which is why a shallow
/// gene is not pruned at a fraction of a read.
///
/// WG-05: the lower bound used to be `min_abund`, the *allele* floor, so `--min-abund 20` on
/// a 60x gene silently moved this floor from 5.0 to 20.0 and could change the published
/// primary, while on a gene under 500x every `--prune-frac` below `5/depth` meant the same
/// thing. The two are separate settings now, with the same default each, so a default run
/// prunes exactly as it did.
pub fn depth_floor(cfg: &Config, depth: f32) -> f32 {
    relative_floor(cfg.prune_frac, depth, cfg.prune_min_abund)
}

/// A floor at `frac` of the gene's depth, never below `min_abund`.
///
/// WV-23. The prune uses it with `--prune-frac` and the allele gate with `--min-abund-frac`;
/// the two differ only in which fraction they pass, so they share the arithmetic.
pub fn relative_floor(frac: f64, depth: f32, min_abund: f32) -> f32 {
    if frac > 0.0 {
        (frac as f32 * depth).max(min_abund)
    } else {
        min_abund
    }
}

/// Whether a unitig at `abund` clears a coverage floor. At the floor exactly, it does.
///
/// TW-06. Three sites pruned on this comparison: `prune_nodes`, the surviving-node count in
/// `pipeline::run_one` (which carried the comment "Predicate must match `prune_nodes`
/// exactly") and the branch gate in `phase`. A comment is not a guarantee: if one of them
/// had drifted to `>`, the count reported in `PruneReport` would have disagreed with the
/// pruning actually done, and a unitig sitting exactly on the floor would have been both
/// kept and counted as dropped. Stated once, it cannot drift.
pub fn survives(abund: f32, floor: f32) -> bool {
    abund >= floor
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
    protected: &std::collections::HashSet<UnitigId>,
) -> (Vec<UnitigId>, Vec<Vec<u8>>) {
    let mut kn = Vec::with_capacity(nodes.len());
    let mut ks = Vec::with_capacity(seqs.len());
    for (id, seq) in nodes.into_iter().zip(seqs) {
        if kept_by_prune(graph, id, floor, protected) {
            kn.push(id);
            ks.push(seq);
        }
    }
    (kn, ks)
}

/// Whether a prune at `floor` keeps unitig `id`. **The one predicate**: `prune_nodes` and the
/// surviving-node count in `pipeline::run_one` must agree, or the count the summary reports
/// describes a different prune from the one performed (TW-06).
pub fn kept_by_prune(
    graph: &UnitigGraph,
    id: UnitigId,
    floor: f32,
    protected: &std::collections::HashSet<UnitigId>,
) -> bool {
    survives(graph.abundance_of(id), floor) || protected.contains(&id)
}

/// Route unitigs that are below `floor` but are the only way through at that point (WP-01).
///
/// The prune floor is one number for the whole gene, and `prune_nodes` drops every unitig
/// below it. That is right for a reference-matching error branch sitting beside a deep true
/// allele, which is what the prune exists to remove. It is wrong for a stretch of the gene
/// that is genuinely covered below the floor: amplicon dropout, a 3' fall-off, any uneven
/// coverage. Removing it cuts the route in two, the realignment cannot span the reference,
/// the prune is refused, and the reference-pulled first-pass route stands with its reference
/// alleles intact. `prefer_majority` then skips those sites too, so nothing corrects them and
/// the only trace is `prune.outcome`.
///
/// The test is whether the route had a **choice**. For each sub-floor route unitig, look at
/// the branch points on either side of it: the edges leaving the predecessor from the end the
/// route leaves by, and the edges entering the successor from the end the route enters by. If
/// any sibling there clears the floor, the aligner could have gone that way and this unitig is
/// a droppable alternative. If no sibling does, this unitig is the gene and it is protected.
///
/// Route ends count as "no sibling": there is nothing on that side to choose between.
pub fn protected_route_nodes(
    graph: &UnitigGraph,
    steps: &[PathStep],
    floor: f32,
) -> Result<std::collections::HashSet<UnitigId>> {
    let mut protected = std::collections::HashSet::new();
    let mut row: Vec<u32> = Vec::new();
    // A sibling the extracted neighbourhood does not hold is not an alternative the aligner
    // could have taken: it is not in the GFA that was handed to it. `edges` deliberately
    // returns whole rows, targets outside the neighbourhood included, so every sibling has to
    // be checked for membership before its abundance is read. Reading it unchecked panicked on
    // exactly the accessions this protection exists for, the ones tangled enough that a 5-hop
    // neighbourhood stops short of some route unitig's siblings.
    let sibling_survives = |graph: &UnitigGraph, t: UnitigId| -> bool {
        graph.has(t) && survives(graph.abundance_of(t), floor)
    };
    let read = |graph: &UnitigGraph, id: UnitigId, row: &mut Vec<u32>| -> Result<()> {
        graph.edges_by_id_into(id, row).map_err(|e| Error::Parse {
            accession: graph.accession.clone(),
            detail: format!("reading adjacency for unitig {id}: {e}"),
        })
    };
    for (i, step) in steps.iter().enumerate() {
        if survives(graph.abundance_of(step.node), floor) {
            continue; // clears the floor on its own; nothing to protect it from
        }
        let mut has_alternative = false;
        // Leaving the predecessor: siblings are the other targets of the edges that leave the
        // predecessor by the same end this step was reached through.
        if let Some(prev) = steps.get(i.wrapping_sub(1)) {
            read(graph, prev.node, &mut row)?;
            let by = row
                .iter()
                .map(|&p| crate::stages::parse::unpack_edge(p))
                .find(|&(t, _, _)| t == step.node)
                .map(|(_, src_fwd, _)| src_fwd);
            if let Some(src_fwd) = by {
                has_alternative |= row
                    .iter()
                    .map(|&p| crate::stages::parse::unpack_edge(p))
                    .any(|(t, s, _)| {
                        s == src_fwd && t != step.node && sibling_survives(graph, t)
                    });
            }
        }
        // Entering the successor: the mirror image, read from the successor's own row because
        // Logan records every edge from both endpoints.
        if let Some(next) = steps.get(i + 1) {
            read(graph, next.node, &mut row)?;
            let by = row
                .iter()
                .map(|&p| crate::stages::parse::unpack_edge(p))
                .find(|&(t, _, _)| t == step.node)
                .map(|(_, src_fwd, _)| src_fwd);
            if let Some(src_fwd) = by {
                has_alternative |= row
                    .iter()
                    .map(|&p| crate::stages::parse::unpack_edge(p))
                    .any(|(t, s, _)| {
                        s == src_fwd && t != step.node && sibling_survives(graph, t)
                    });
            }
        }
        if !has_alternative {
            protected.insert(step.node);
        }
    }
    Ok(protected)
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
/// this touches 422 entries rather than rescanning 1,277,362. Without a sequence index this
/// step would have to re-read the entire unitig file.
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
/// The byte layout is fixed, because GraphAligner's output depends on its input and the
/// golden fixtures pin that output. Specifically:
///
/// * `H	VN:Z:1.0` header
/// * one `S` line per node, ascending id, sequence expanded
/// * `L` lines emitted per source node in ascending id order, and within a node in the order
///   the links appeared in that unitig's header
/// * only edges whose *both* endpoints are in the subgraph
/// * overlap always `{k-1}M`, which is exact for a de Bruijn graph of order k
///
/// Each edge therefore appears twice, once from each endpoint, because Logan records it in
/// both headers. The duplication is preserved deliberately.
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

    // WF-01, then WF-17. This used to be a bare `create_dir_all`, which silently put back a
    // scratch directory another run had deleted mid-flight, and put it back with the default
    // mode 0755 instead of 0700. That silence is half of why the deletion surfaced as the
    // biological answer `no_path` rather than as a fault. WF-01 made it warn and recreate;
    // it now refuses instead, because recreating is wrong in every case that is left.
    //
    // `Scratch::create` makes this directory and holds its `.lock` for the life of the run,
    // so another live run cannot sweep it any more. What can still remove it mid-run is this
    // run's own signal handler, on its way out, and a worker that recreates the directory
    // then puts it straight back after the cleanup that was supposed to remove it, with a
    // GFA and a GAF inside. The remaining case, somebody deleting the directory by hand, is
    // a fault too. So: fail the accession and say what happened.
    if !cfg.scratch_dir.exists() {
        return Err(crate::error::file_err("write into", &cfg.scratch_dir)(
            std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "the scratch directory vanished during the run. Either this run is shutting \
                 down on a signal, or another process removed it; if several runs share a \
                 --scratch-dir, give each its own",
            ),
        ));
    }
    let gfa_path = cfg.scratch_dir.join(format!("{}.sub.gfa", scratch_stem(&graph.accession)));
    // A failed write (a full tmpfs) removes the partial file instead of leaving it behind.
    let n_edges = write_gfa_file(cfg, graph, &nodes, &seqs, &gfa_path).map_err(|e| {
        let _ = std::fs::remove_file(&gfa_path);
        crate::error::file_err("write", &gfa_path)(e)
    })?;

    Ok(Subgraph {
        accession: graph.accession.clone(),
        nodes,
        seqs,
        gfa_path,
        n_edges,
        // This function just wrote the file, so it is the one place that may unlink it.
        owns_gfa: true,
    })
}

/// The GFA bytes for `write_gfa`. Returns the number of `L` lines.
/// `create_dir_all` that leaves the final component mode 0700, not the default 0755.
#[cfg(unix)]
fn private_dir_all(path: &std::path::Path) -> std::io::Result<()> {
    use std::os::unix::fs::DirBuilderExt;
    std::fs::DirBuilder::new().mode(0o700).recursive(true).create(path)
}

#[cfg(not(unix))]
fn private_dir_all(path: &std::path::Path) -> std::io::Result<()> {
    std::fs::create_dir_all(path)
}

fn write_gfa_file(
    cfg: &Config,
    graph: &UnitigGraph,
    nodes: &[UnitigId],
    seqs: &[Vec<u8>],
    gfa_path: &std::path::Path,
) -> std::io::Result<usize> {
    let file = std::fs::File::create(gfa_path)?;
    let mut w = BufWriter::with_capacity(1 << 16, file);

    // Membership test for "is this edge target inside the subgraph". `nodes` is sorted, so a
    // bitset keyed by unitig id is both smaller and faster than a binary search per edge.
    let mut inside = vec![0u64; graph.n_unitigs.div_ceil(64)];
    for &id in nodes {
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
    for &id in nodes {
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
    Ok(n_edges)
}

#[cfg(test)]
mod collect_tests {
    use super::*;
    use crate::stages::parse;

    fn graph(fa: &str) -> UnitigGraph {
        parse::build(&"T".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX).unwrap()
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
        parse::build(&"T".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX).unwrap()
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

    /// The guard measures the neighbourhood against the reference, not against the sample's
    /// graph, so a targeted library whose graph is almost all target does not trip it.
    #[test]
    fn the_guard_is_per_reference_base_not_per_graph() {
        let g = chain();
        let nodes = expand(&cfg_hops(1), &g, &anchors(&[0], 6)).unwrap();
        assert_eq!(nodes.len(), 2);
        let mut cfg = cfg_hops(1);
        cfg.subgraph_nodes_per_base_guard = 3.0;

        // 2 unitigs over a 1,000 base reference: a third of a percent of the guard, even
        // though those 2 unitigs are a third of the whole graph.
        assert!(!exceeds_guard(&cfg, &nodes, 1_000));
        // The same 2 unitigs against a 1-base reference are 2.0 per base, still under.
        assert!(!exceeds_guard(&cfg, &nodes, 1));
        // Strictly above the guard trips it; exactly at it does not.
        cfg.subgraph_nodes_per_base_guard = 2.0;
        assert!(!exceeds_guard(&cfg, &nodes, 1));
        cfg.subgraph_nodes_per_base_guard = 1.9;
        assert!(exceeds_guard(&cfg, &nodes, 1));
    }

    /// A query length of zero would otherwise divide by zero and flag every accession.
    #[test]
    fn an_empty_reference_does_not_trip_the_guard() {
        let g = chain();
        let nodes = expand(&cfg_hops(1), &g, &anchors(&[0], 6)).unwrap();
        assert_eq!(nodes_per_ref_base(&nodes, 0), 0.0);
        assert!(!exceeds_guard(&cfg_hops(1), &nodes, 0));
    }
}

/// Cross-validation against the anchor-recall study's recorded subgraph sizes.
///
/// `LOGAN_TEST_DIR=/dir/with/*.unitigs.fa.zst cargo test --release -- --ignored --nocapture`
#[cfg(test)]
mod recall_study_subgraphs {
    use super::*;
    use crate::io::input::read_query;
    use crate::stages::{anchor, parse, stream};
    use crate::types::FetchedUnitigs;
    use std::path::Path;

    /// Load a study fixture. `None` means the file is absent; anything else panics.
    ///
    /// TW-12: read, decode and parse failures all used to become `None`, and the caller
    /// printed `SKIP: no data for {acc}` and carried on. A truncated or non-zstd fixture
    /// therefore reported `ok` locally with a false reason, and in CI the `SKIP:` guard
    /// failed the job with a message pointing at missing data. The sibling tests panic.
    fn load(dir: &str, acc: &str) -> Option<UnitigGraph> {
        let p = format!("{dir}/{acc}.unitigs.fa.zst");
        let compressed = match std::fs::read(&p) {
            Ok(b) => b,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return None,
            Err(e) => panic!("{p}: {e}"),
        };
        let f = FetchedUnitigs { index: 0, accession: acc.into(), compressed_len: compressed.len(), compressed: crate::types::Compressed::Memory(compressed), fetch_secs: 0.0, payload: crate::types::Payload::default() };
        let rd = stream::decoder(&f).unwrap_or_else(|e| panic!("{p}: not decodable: {e}"));
        Some(parse::build(&f.accession, rd, u64::MAX).unwrap_or_else(|e| panic!("{p}: not parseable: {e}")))
    }

    /// TW-12: a corrupt fixture must fail loudly, not be reported as "no data".
    #[test]
    // The zstd reader is lazy, so a bad frame surfaces at parse rather than at decode. Match
    // the filename, which both panics carry, rather than pinning which stage noticed.
    #[should_panic(expected = "ERRCORRUPT.unitigs.fa.zst")]
    fn a_corrupt_fixture_panics_rather_than_reading_as_absent() {
        let dir = std::env::temp_dir().join(format!("logan-recall-corrupt-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("ERRCORRUPT.unitigs.fa.zst");
        std::fs::write(&p, b"this is plainly not a zstd frame").unwrap();
        let out = std::panic::catch_unwind(|| load(dir.to_str().unwrap(), "ERRCORRUPT"));
        let _ = std::fs::remove_dir_all(&dir);
        match out {
            Err(e) => std::panic::resume_unwind(e),
            Ok(_) => panic!("a corrupt fixture must not load"),
        }
    }

    #[test]
    #[ignore]
    fn matches_the_recall_study_subgraph_sizes() {
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
                // The study counted canonical 31-mers; the default -k has since moved to 21.
                let cfg = Config { k: 31, min_shared_kmers: m, ..Config::default() };
                let idx = anchor::index_query(&cfg, &query).unwrap();
                let a = anchor::find(&cfg, &graph, &idx).unwrap();
                let t0 = std::time::Instant::now();
                let nodes = expand(&cfg, &graph, &a).unwrap();
                eprintln!(
                    "  M={:<3} anchors {:>5} (study {:>5})  subgraph {:>6} (study {:>6})  {:>7.4}s  \
                     frac {:>7.3}%",
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

/// Unitigs in the neighbourhood per base of the reference. The measure the guard is calibrated
/// on; see `Config::subgraph_nodes_per_base_guard`.
pub fn nodes_per_ref_base(nodes: &[UnitigId], ref_len: usize) -> f64 {
    if ref_len == 0 {
        0.0
    } else {
        nodes.len() as f64 / ref_len as f64
    }
}

/// True when the local neighbourhood is tangled enough around the gene that the eventual
/// haplotype is likely under-determined. See `Config::subgraph_nodes_per_base_guard`.
pub fn exceeds_guard(cfg: &Config, nodes: &[UnitigId], ref_len: usize) -> bool {
    nodes_per_ref_base(nodes, ref_len) > cfg.subgraph_nodes_per_base_guard
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
            compressed: crate::types::Compressed::Memory(compressed),
            fetch_secs: 0.0,
            payload: crate::types::Payload::default(),
        };
        let graph = {
            let rd = stream::decoder(&f).unwrap();
            parse::build(&f.accession, rd, u64::MAX).unwrap()
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

        // And the subgraph path specifically. 422 nodes is the k=31 subgraph (the recall
        // study's, above); the default -k has since moved to 21.
        let query = read_query(Path::new("../queries/pfcrt_only.fasta")).unwrap();
        let cfg = Config { k: 31, ..Config::default() };
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
        parse::build(&"ACC".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX).unwrap()
    }
    /// A config and the guard that removes its scratch directory (TW-10). Bind the guard,
    /// usually as `_d`, for as long as the config is in use.
    fn tmp_cfg(tag: &str) -> (Config, crate::scratch::TestDir) {
        let d = crate::scratch::TestDir::new(tag);
        (Config { scratch_dir: d.path(), ..Config::default() }, d)
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
    fn writes_the_exact_gfa_line_format() {
        let g = graph(">ACC_0 ka:f:5.3   L:-:1:+  L:+:1:-\nACGT\n>ACC_1 ka:f:9.0   L:+:0:+\nTTGCA\n");
        let (base, _d) = tmp_cfg("fmt");
        let cfg = Config { graph_k: 31, ..base };
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

    /// An L line is written only when both endpoints survived into the subgraph. Keeping a
    /// dangling edge would make GraphAligner reference a node the file never declares.
    #[test]
    fn drops_edges_leaving_the_subgraph() {
        let g = graph(
            ">ACC_0 L:+:1:+ L:+:2:+\nACGT\n>ACC_1 L:+:0:+\nTTTT\n>ACC_2 L:+:0:+\nGGGG\n",
        );
        let (cfg, _d) = tmp_cfg("dangling");
        let got = render(&cfg, &g, vec![0, 1]);
        assert!(got.contains("L\t0\t+\t1\t+\t30M\n"), "kept edge missing:\n{got}");
        assert!(!got.contains("\t2\t"), "edge to excluded node 2 leaked:\n{got}");
        assert_eq!(got.matches("S\t").count(), 2);
    }

    #[test]
    fn overlap_follows_k() {
        let g = graph(">ACC_0 L:+:1:+\nACGT\n>ACC_1 L:+:0:+\nTTTT\n");
        for (k, want) in [(31usize, "30M"), (21, "20M")] {
            let (base, _d) = tmp_cfg(&format!("k{k}"));
            let cfg = Config { graph_k: k, ..base };
            let got = render(&cfg, &g, vec![0, 1]);
            assert!(got.contains(want), "k={k} expected {want} in:\n{got}");
        }
    }

    #[test]
    fn reports_edge_count_and_path() {
        let g = graph(">ACC_0 L:+:1:+\nACGT\n>ACC_1 L:+:0:+\nTTTT\n");
        let (cfg, _d) = tmp_cfg("count");
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

    /// An edge pointing past the end of the graph: refused by `expand`, skipped by the GFA
    /// writer, a panic in neither (TW-06).
    ///
    /// Logan's headers are reciprocal and in range, so this is malformed input rather than
    /// something seen in the wild; the point is that the two readers of `graph.edges` treat
    /// it the way they treat every other structural fault, which is to say not with an
    /// out-of-bounds index. They differ in kind deliberately: `expand` is choosing the
    /// search space and cannot quietly drop part of it, while the writer is emitting the
    /// space already chosen, and a target outside the graph cannot be one of `nodes`, so
    /// there is nothing to emit and nothing to report.
    ///
    /// The fixture is assembled by hand because **the parser no longer hands one over**:
    /// WG-10 made a link past the last record a parse error, since the dense id rule means
    /// such a link can only mean records are missing. These two guards stay because they are
    /// what keeps a structural fault from becoming an out-of-bounds index, and the parser
    /// being strict is not a reason to stop testing them.
    #[test]
    fn an_edge_past_the_end_of_the_graph_is_refused_then_skipped() {
        use crate::types::Anchors;
        // WG-10: the same text through the parser is now refused outright.
        let text = concat!(
            ">ACC_0 ka:f:10.0  L:+:1:+  L:+:99:+\nACGT\n",
            ">ACC_1 ka:f:10.0  L:-:0:-\nTTTT\n",
        );
        match parse::build(
            &"ACC".into(),
            Box::new(std::io::Cursor::new(text.as_bytes().to_vec())),
            u64::MAX,
        ) {
            Err(e) => assert!(e.to_string().contains("truncated or incomplete"), "{e}"),
            Ok(g) => panic!("a link past the last record was parsed into {} unitigs", g.n_unitigs),
        }

        // So build the valid graph and splice the dangling edge into node 0's CSR row.
        let mut g = graph(concat!(
            ">ACC_0 ka:f:10.0  L:+:1:+\nACGT\n",
            ">ACC_1 ka:f:10.0  L:-:0:-\nTTTT\n",
        ));
        let row_end = g.offsets[1] as usize;
        g.targets.insert(row_end, parse::pack_edge(99, true, true));
        for o in g.offsets[1..].iter_mut() {
            *o += 1;
        }
        assert_eq!(g.n_unitigs, 2);

        let (cfg, _d) = tmp_cfg("oob");
        let anchors = Anchors { ids: vec![0], shared_kmers: vec![1], fraction_of_graph: 0.5 };
        let err = expand(&cfg, &g, &anchors).expect_err("an edge out of range must be reported");
        assert!(
            matches!(err, Error::Parse { .. }) && err.to_string().contains("99"),
            "expected a parse error naming the target, got {err}"
        );

        // The writer is handed both real nodes and must drop only the dangling L line.
        let seqs = collect_seqs(&g, &[0, 1]).unwrap();
        let sub = write_gfa(&cfg, &g, vec![0, 1], seqs).unwrap();
        let mut got = String::new();
        std::fs::File::open(&sub.gfa_path).unwrap().read_to_string(&mut got).unwrap();
        assert!(got.contains("L\t0\t+\t1\t+\t30M\n"), "the in-range edge is kept:\n{got}");
        assert!(got.contains("L\t1\t-\t0\t-\t30M\n"), "its mirror is kept:\n{got}");
        assert!(!got.contains("99"), "the dangling edge leaked:\n{got}");
        assert_eq!(sub.n_edges, 2, "two L lines, not three");
    }

    #[test]
    fn empty_subgraph_writes_a_header_only_file() {
        // GraphAligner aborts with SIGABRT on a nodeless graph, so align must refuse to run
        // rather than be handed this. Writing it is still valid GFA.
        let g = graph(">ACC_0\nACGT\n");
        let (cfg, _d) = tmp_cfg("empty");
        assert_eq!(render(&cfg, &g, vec![]), "H\tVN:Z:1.0\n");
    }

    #[test]
    fn mismatched_nodes_and_seqs_is_rejected() {
        let g = graph(">ACC_0\nACGT\n>ACC_1\nTTTT\n");
        let (cfg, _d) = tmp_cfg("mismatch");
        assert!(write_gfa(&cfg, &g, vec![0, 1], vec![b"ACGT".to_vec()]).is_err());
    }
}

/// Emit a real GFA plus its anchor list so another implementation can be run on exactly
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
            compressed: crate::types::Compressed::Memory(compressed),
            fetch_secs: 0.0,
            payload: crate::types::Payload::default(),
        };
        let graph = {
            let rd = stream::decoder(&f).unwrap();
            parse::build(&f.accession, rd, u64::MAX).unwrap()
        };
        let query = read_query(Path::new("../queries/pfcrt_only.fasta")).unwrap();
        let cfg = Config {
            scratch_dir: PathBuf::from(&outdir),
            ..Config::default()
        };
        // `write_gfa` no longer creates a missing scratch directory (WF-17), so a test that
        // names its own creates it, as a run does.
        std::fs::create_dir_all(&cfg.scratch_dir).unwrap();
        let idx = anchor::index_query(&cfg, &query).unwrap();
        let a = anchor::find(&cfg, &graph, &idx).unwrap();
        let nodes = expand(&cfg, &graph, &a).unwrap();
        let seqs = collect_seqs(&graph, &nodes).unwrap();
        let sub = write_gfa(&cfg, &graph, nodes, seqs).unwrap();

        // TW-13: copy the GFA to a stable name BEFORE `sub` drops. `Drop for Subgraph`
        // unlinks the file it created, so this helper used to report "wrote <path>" and then
        // delete it, leaving `LOGAN_GFA_OUT` holding only anchors.txt. The name it printed
        // also embedded the PID, so it was not a name another tool could be pointed at.
        //
        // This is not TW-01: here the walker really did create the file, so the ownership
        // flag is correct and the copy is the fix.
        let kept = format!("{outdir}/{}.sub.gfa", f.accession);
        std::fs::copy(&sub.gfa_path, &kept)
            .unwrap_or_else(|e| panic!("could not keep {}: {e}", sub.gfa_path.display()));

        // The exact anchor ids, so another implementation can be driven with the same input.
        let list: Vec<String> = a.ids.iter().map(|i| i.to_string()).collect();
        std::fs::write(format!("{outdir}/anchors.txt"), list.join("\n")).unwrap();
        eprintln!(
            "  wrote {} ({} nodes, {} edges) and {} anchors",
            kept,
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
        parse::build(&"T".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX).unwrap()
    }
    fn steps(ids: &[UnitigId]) -> Vec<PathStep> {
        ids.iter().map(|&node| PathStep { node, orient: Orient::Forward }).collect()
    }

    /// WP-01 asked for a k-mer-weighted median here. It was implemented and then removed:
    /// see `path_depth`'s own note for the measurement that refutes it. This pins the
    /// statistic that survived, so a future attempt starts from the number it has to beat.
    #[test]
    fn depth_is_the_step_median_not_the_kmer_weighted_one() {
        let c = Config { graph_k: 4, ..Config::default() };
        let g = graph(concat!(
            ">T_0 ka:f:10.0\nACGT\n",
            ">T_1 ka:f:10.0\nACGT\n",
            ">T_2 ka:f:1000.0\n",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n",
        ));
        assert_eq!(
            path_depth(&c, &g, &steps(&[0, 1, 2])),
            10.0,
            "two steps at 10 outvote one at 1000, however long that one is"
        );
    }

    /// WP-01. A sub-floor route unitig with a deeper sibling is an error branch and is
    /// pruned; one with no deeper sibling is the only way through and is kept.
    ///
    /// Removing the second cuts the gene's route in two. The realignment then cannot span the
    /// reference, the prune is refused, and the reference-pulled first-pass route stands with
    /// its reference alleles intact, with nothing in the output but `prune.outcome` to say so.
    #[test]
    fn a_thin_stretch_of_the_gene_is_protected_but_an_error_branch_is_not() {
        // 0 -> 1 -> 2 is the route; 3 is a sibling of 1, reachable from 0 by the same end.
        let with_sibling_at = |ka: &str| {
            graph(&format!(
                ">T_0 ka:f:100.0   L:+:1:+  L:+:3:+\nACGTACGTACGT\n\
                 >T_1 ka:f:1.0   L:-:0:-  L:+:2:+\nACGTACGTACGT\n\
                 >T_2 ka:f:100.0   L:-:1:-\nACGTACGTACGT\n\
                 >T_3 ka:f:{ka}   L:-:0:-\nACGTACGTACGT\n"
            ))
        };
        let route = steps(&[0, 1, 2]);

        // The sibling clears the floor, so the aligner could have gone that way: unitig 1 is
        // a droppable alternative, which is exactly what the prune is for.
        let deep = with_sibling_at("100.0");
        assert!(
            protected_route_nodes(&deep, &route, 10.0).unwrap().is_empty(),
            "an error branch beside a deep sibling must still be pruned"
        );

        // The sibling is below the floor too, so there is no deeper way round: unitig 1 is
        // the gene, thinly covered.
        let thin = with_sibling_at("1.0");
        let protected = protected_route_nodes(&thin, &route, 10.0).unwrap();
        assert_eq!(protected.len(), 1);
        assert!(protected.contains(&1), "{protected:?}");

        // And the predicate the prune and its reported count share honours it.
        assert!(!kept_by_prune(&thin, 1, 10.0, &Default::default()));
        assert!(kept_by_prune(&thin, 1, 10.0, &protected));
    }

    /// A sibling outside the extracted neighbourhood is not an alternative, and asking the
    /// graph for its abundance used to panic the accession.
    ///
    /// `edges` returns whole rows on purpose, so a route unitig at the edge of the
    /// neighbourhood has siblings the graph does not hold. The aligner never saw them, since
    /// they are not in the GFA, so they cannot be the choice this protection is looking for.
    /// Treating them as absent protects the unitig, which is the conservative direction: the
    /// prune keeps a unitig it might have dropped rather than cutting the route in two.
    #[test]
    fn a_sibling_outside_the_neighbourhood_does_not_panic_and_is_not_an_alternative() {
        // 0 -> 1 -> 2 is the route; 3 is a deep sibling of 1, reachable from 0 by the same
        // end, and is left out of the collected neighbourhood.
        const FA: &str = ">T_0 ka:f:100.0   L:+:1:+  L:+:3:+\nACGTACGTACGT\n\
                          >T_1 ka:f:1.0   L:-:0:-  L:+:2:+\nACGTACGTACGT\n\
                          >T_2 ka:f:100.0   L:-:1:-\nACGTACGTACGT\n\
                          >T_3 ka:f:100.0   L:-:0:-\nACGTACGTACGT\n";
        let cursor = || Box::new(std::io::Cursor::new(FA.as_bytes().to_vec())) as Box<dyn std::io::Read>;
        let skeleton = parse::scan(&"T".into(), cursor(), u64::MAX, &mut ()).expect("scan");
        let sub = parse::collect(skeleton, cursor(), &[0, 1, 2], u64::MAX).expect("collect");
        assert!(!sub.has(3), "unitig 3 must be outside the neighbourhood for this test");

        let protected = protected_route_nodes(&sub, &steps(&[0, 1, 2]), 10.0).unwrap();
        assert!(
            protected.contains(&1),
            "a sibling the neighbourhood does not hold is not a choice the aligner had: {protected:?}"
        );

        // The same graph with unitig 3 retained does see the alternative, so the rule above
        // is about membership and not about having quietly stopped looking at siblings.
        let whole = parse::build(&"T".into(), cursor(), u64::MAX).unwrap();
        assert!(
            protected_route_nodes(&whole, &steps(&[0, 1, 2]), 10.0).unwrap().is_empty(),
            "with the sibling retained, unitig 1 is a droppable error branch"
        );
    }

    /// A route whose every unitig clears the floor protects nothing: there is nothing at risk.
    #[test]
    fn nothing_is_protected_when_the_route_clears_the_floor() {
        let g = graph(concat!(
            ">T_0 ka:f:100.0   L:+:1:+\nACGTACGTACGT\n",
            ">T_1 ka:f:100.0   L:-:0:-\nACGTACGTACGT\n",
        ));
        assert!(protected_route_nodes(&g, &steps(&[0, 1]), 10.0).unwrap().is_empty());
    }

    #[test]
    fn depth_is_the_median_along_the_alignment() {
        let g = graph(concat!(
            ">T_0 ka:f:100.0\nACGT\n",
            ">T_1 ka:f:200.0\nACGT\n",
            ">T_2 ka:f:300.0\nACGT\n",
        ));
        assert_eq!(path_depth(&Config::default(), &g, &steps(&[0, 1, 2])), 200.0);
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
        assert_eq!(path_depth(&Config::default(), &g, &steps(&[0, 1, 2])), 100.0, "mean would be 8400");
    }

    #[test]
    fn an_empty_alignment_has_no_depth() {
        let g = graph(">T_0 ka:f:100.0\nACGT\n");
        assert_eq!(path_depth(&Config::default(), &g, &[]), 0.0);
    }

    /// At the floor exactly, a unitig is kept (TW-06).
    ///
    /// `depth_floor` computes the floor as an f32 product, so "exactly at the floor" is a
    /// value the pipeline really does produce, not a contrived one: 1% of a 14,500x gene is
    /// 145.0, and a unitig at 145.0 is the commonest kind of boundary case here because the
    /// floor is derived from the median of these same abundances. The rule is `>=`, and it
    /// is `survives` that says so for all three sites that prune.
    #[test]
    fn a_unitig_exactly_at_the_floor_survives() {
        let c = Config { prune_frac: 0.01, min_abund: 5.0, ..Config::default() };
        let floor = depth_floor(&c, 14_500.0);
        assert_eq!(floor, 145.0, "the floor this fixture is built around");

        assert!(survives(floor, floor), "at the floor is in");
        assert!(survives(floor + 0.1, floor));
        assert!(!survives(144.9, floor), "below the floor is out");

        let g = graph(concat!(
            ">T_0 ka:f:145.0\nACGT\n",
            ">T_1 ka:f:144.9\nTTTT\n",
            ">T_2 ka:f:145.1\nGGGG\n",
        ));
        let (kept, seqs) =
            prune_nodes(&g, vec![0, 1, 2], vec![vec![b'A'], vec![b'C'], vec![b'G']], floor, &Default::default());
        assert_eq!(kept, vec![0, 2], "only the unitig under the floor is removed");
        assert_eq!(seqs, vec![vec![b'A'], vec![b'G']], "sequences stay paired with their nodes");
    }

    /// Pruning is on by default at 1% of depth; `0.0` is how it is switched off, and that
    /// path must stay an exact no-op.
    #[test]
    fn the_floor_is_off_when_prune_frac_is_zero() {
        assert_eq!(Config::default().prune_frac, 0.01, "the panel-calibrated default");
        let off = Config { prune_frac: 0.0, ..Config::default() };
        assert_eq!(depth_floor(&off, 14500.0), off.prune_min_abund);
        assert_eq!(depth_floor(&Config::default(), 14500.0), 145.0);
    }

    /// The whole point: five reads means something different at 180x and at 14,500x.
    #[test]
    fn the_floor_scales_with_depth() {
        let c = Config { prune_frac: 0.01, ..Config::default() };
        assert_eq!(depth_floor(&c, 14500.0), 145.0);
        assert_eq!(
            depth_floor(&c, 180.0),
            5.0,
            "1% of 180x is under prune_min_abund, which holds"
        );
    }

    #[test]
    fn the_floor_never_drops_below_prune_min_abund() {
        let c = Config { prune_frac: 0.9, prune_min_abund: 5.0, ..Config::default() };
        for d in [0.0, 1.0, 5.0, 5.5, 1000.0] {
            assert!(depth_floor(&c, d) >= 5.0, "floor dipped at depth {d}");
        }
    }

    /// WG-05. `--min-abund` is the allele floor. It used to be the prune's lower bound as
    /// well, so a user raising it to be stricter about which forks count as variants also
    /// raised the coverage at which unitigs are deleted from the graph the primary is walked
    /// on. On a 60x gene with a 42x/18x fork that moved the floor from 5.0 to 20.0, removed
    /// the reference-matching branch and published a different sequence, with nothing in the
    /// help text, the docs or the summary saying the two settings were one.
    #[test]
    fn min_abund_does_not_move_the_prune_floor() {
        let strict = Config { min_abund: 20.0, ..Config::default() };
        assert_eq!(
            depth_floor(&strict, 60.0),
            depth_floor(&Config::default(), 60.0),
            "the allele floor must not decide which unitigs the prune deletes"
        );
        assert_eq!(depth_floor(&strict, 60.0), 5.0);

        // The prune's own setting still does move it, which is the point of splitting them.
        let pruned_hard = Config { prune_min_abund: 20.0, ..Config::default() };
        assert_eq!(depth_floor(&pruned_hard, 60.0), 20.0);

        // And the difference is visible where it matters: the 18x branch of a 60x gene.
        let g = graph(concat!(
            ">T_0 ka:f:60.0\nAAAA\n",
            ">T_1 ka:f:18.0\nCCCC\n",
            ">T_2 ka:f:42.0\nGGGG\n",
        ));
        let seqs = || vec![vec![b'A'], vec![b'C'], vec![b'G']];
        let floor = depth_floor(&strict, 60.0);
        let (kept, _) = prune_nodes(&g, vec![0, 1, 2], seqs(), floor, &Default::default());
        assert_eq!(kept, vec![0, 1, 2], "--min-abund 20 must leave the 18x unitig in the graph");
        let floor = depth_floor(&pruned_hard, 60.0);
        let (kept, _) = prune_nodes(&g, vec![0, 1, 2], seqs(), floor, &Default::default());
        assert_eq!(kept, vec![0, 2], "--prune-min-abund 20 is how that unitig is removed");
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
        let floor = depth_floor(&c, path_depth(&c, &g, &steps(&[0, 1, 2])));
        let (kept, seqs) = prune_nodes(&g, vec![0, 1, 2], vec![vec![b'A']; 3], floor, &Default::default());
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
        let floor = depth_floor(&c, path_depth(&c, &g, &steps(&[0, 1])));
        assert_eq!(floor, 150.0);
        let (kept, _) = prune_nodes(&g, vec![0, 1, 2], vec![vec![b'A']; 3], floor, &Default::default());
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
            &Default::default(),
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
