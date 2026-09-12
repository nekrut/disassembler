//! Threading structure and stage wiring. This is the file to read to understand the flow.
//!
//! ```text
//!  io::input                                                          io::output
//!  ─────────                                                          ──────────
//!  accession list                                              reconstructed.fasta
//!  reference CDS                                               bubbles.json
//!       │                                                      summary.json
//!       │                                                            ▲
//!       ▼                                                            │
//!  ┌─────────────────────────── pipeline::run ────────────────────────────┐
//!  │                                                                      │
//!  │  work queue (accessions)                                             │
//!  │       │                                                              │
//!  │       ▼                                                              │
//!  │  FETCH POOL  cfg.fetch_workers (default 8, = bandwidth saturation)   │
//!  │   ┌────────┐ ┌────────┐        ┌────────┐                            │
//!  │   │fetch 0 │ │fetch 1 │  ...   │fetch F │   stages::stream           │
//!  │   └───┬────┘ └───┬────┘        └───┬────┘                            │
//!  │       └──────────┴─────────────────┘                                 │
//!  │                  │                                                   │
//!  │        sync_channel(cfg.channel_capacity)   ◄── BOUNDED. capacity x  │
//!  │                  │                              buffer size is the   │
//!  │                  │                              hard memory ceiling. │
//!  │                  ▼                                                   │
//!  │  COMPUTE POOL  cfg.compute_workers (= cores)                         │
//!  │   ┌────────┐ ┌────────┐        ┌────────┐                            │
//!  │   │ work 0 │ │ work 1 │  ...   │ work W │   run_one() below          │
//!  │   └───┬────┘ └───┬────┘        └───┬────┘                            │
//!  │       │          │                 │                                 │
//!  │       │   per accession, in order:                                   │
//!  │       │     stages::parse    -> UnitigGraph (CSR, ~20 MB)            │
//!  │       │     stages::anchor   -> Anchors                              │
//!  │       │     stages::extract  -> Subgraph + GFA on tmpfs              │
//!  │       │     stages::align    -> GraphAligner subprocess -> GafPath   │
//!  │       │     stages::stitch   -> Haplotype                            │
//!  │       │     stages::phase    -> Vec<Bubble>                          │
//!  │       │          │                 │                                 │
//!  │       └──────────┴─────────────────┘                                 │
//!  │                  │  results channel (unbounded; results are small)   │
//!  │                  ▼                                                   │
//!  │            COLLECTOR (calling thread) -> Summary                     │
//!  └──────────────────────────────────────────────────────────────────────┘
//! ```
//!
//! # Why two pools and not one
//!
//! Fetch is bandwidth-bound and saturates around 8 concurrent streams at 17.4 MB/s.
//! Compute is CPU-bound and wants one worker per core. `logan_path_walker.py` conflates
//! them behind a single `--workers`, which is part of why `--workers 16` was so
//! destructive: it multiplied memory 16x while buying no extra bandwidth.
//!
//! # Why the channel is bounded
//!
//! Capacity times the largest retained buffer is a provable ceiling on resident memory.
//! Fetch threads block once the queue is full, so a fast network cannot outrun slow compute
//! and pile up buffers. The Python `ThreadPoolExecutor` had no such limit, which is how a
//! single run reached 21.3 GB and took the machine down with it.
//!
//! # Ordering constraint
//!
//! Do not raise concurrency before `stages::parse` produces the bounded CSR. Threads sharing
//! one heap is what capped the original damage; independent workers at the old ~700 bytes
//! per unitig would be worse, not better.

use crate::config::Config;
use crate::error::Result;
use crate::stages;
use crate::types::{Accession, AccessionResult, FetchedUnitigs, Outcome, Query, Timing};
use std::time::Instant;

use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::sync::{Arc, Mutex};

/// `std::sync::mpsc` is single-consumer, so the compute pool shares one receiver behind a
/// mutex. Swap for `crossbeam_channel::Receiver`, which is natively multi-consumer, when
/// dependencies are added.
type SharedReceiver<T> = Arc<Mutex<Receiver<T>>>;

/// What the whole run produced.
pub struct Summary {
    pub total: usize,
    pub success: usize,
    pub results: Vec<AccessionResult>,
}

/// Fan accessions through the two pools and collect results.
pub fn run(cfg: &Config, accessions: Vec<Accession>, query: &Query) -> Result<Summary> {
    let total = accessions.len();
    let (fetch_tx, fetch_rx) = sync_channel::<FetchedUnitigs>(cfg.channel_capacity);
    let (result_tx, result_rx) = std::sync::mpsc::channel::<AccessionResult>();
    let fetch_rx: SharedReceiver<FetchedUnitigs> = Arc::new(Mutex::new(fetch_rx));
    let work: SharedReceiver<Accession> = {
        let (tx, rx) = std::sync::mpsc::channel();
        for a in accessions {
            let _ = tx.send(a);
        }
        Arc::new(Mutex::new(rx))
    };

    // Built once per run, not per accession: it depends only on cfg and query.
    let query_index = stages::anchor::index_query(cfg, query)?;
    let query_index = &query_index;

    // One pooled HTTP client shared by every fetch thread. The bucket is a single host,
    // so per-thread agents would just multiply TLS handshakes.
    let fetcher = stages::stream::Fetcher::new(cfg);
    let fetcher = &fetcher;

    let results = std::thread::scope(|scope| {
        // ---- fetch pool -------------------------------------------------
        for _ in 0..cfg.fetch_workers {
            let work = Arc::clone(&work);
            let tx: SyncSender<FetchedUnitigs> = fetch_tx.clone();
            let rtx = result_tx.clone();
            scope.spawn(move || fetch_worker(fetcher, work, tx, rtx));
        }
        drop(fetch_tx);

        // ---- compute pool -----------------------------------------------
        for _ in 0..cfg.compute_workers {
            let rx = Arc::clone(&fetch_rx);
            let tx = result_tx.clone();
            scope.spawn(move || compute_worker(cfg, query, query_index, rx, tx));
        }
        drop(result_tx);

        // ---- collector (this thread) ------------------------------------
        // Report each accession as it lands. The Python walker prints a per-accession line
        // and so must this: a run that quietly reconstructs nothing looks identical to a
        // successful one otherwise, which is exactly how the int/str regression hid.
        let mut done = 0usize;
        result_rx
            .into_iter()
            .inspect(|r| {
                done += 1;
                report(done, total, r);
            })
            .collect::<Vec<_>>()
    });

    let success = results
        .iter()
        .filter(|r| r.outcome == Outcome::Success)
        .count();
    Ok(Summary {
        total,
        success,
        results,
    })
}

/// One line per finished accession, mirroring `logan_path_walker.py:386-390`.
fn report(done: usize, total: usize, r: &AccessionResult) {
    match r.outcome {
        Outcome::Success => {
            let (cov, pid) = r
                .haplotype
                .as_ref()
                .map(|h| (h.cov, h.pident))
                .unwrap_or((0.0, 0.0));
            let haps_info = if r.haplotypes.len() > 1 {
                format!(", haps={}", r.haplotypes.len())
            } else {
                String::new()
            };
            println!(
                "[{done}/{total}] {}: DONE in {:.1}s (nodes={}, cov={:.2}, pid={:.1}%, bubbles={}{})",
                r.accession,
                r.timing.total,
                r.subgraph_nodes,
                cov,
                pid,
                r.bubbles.len(),
                haps_info
            );
        }
        Outcome::NoAnchors => println!("[{done}/{total}] {}: no_anchors", r.accession),
        Outcome::NoPath => println!("[{done}/{total}] {}: no_path", r.accession),
        Outcome::Failed => println!(
            "[{done}/{total}] {}: ERROR: {}",
            r.accession,
            r.error.as_deref().unwrap_or("unknown")
        ),
    }
}

/// Pull accessions, stream them, hand the bytes downstream. Blocks when the queue is full,
/// which is the backpressure that bounds memory.
fn fetch_worker(
    fetcher: &stages::stream::Fetcher,
    work: SharedReceiver<Accession>,
    out: SyncSender<FetchedUnitigs>,
    result_tx: std::sync::mpsc::Sender<AccessionResult>,
) {
    loop {
        let accession = match work.lock().unwrap().recv() {
            Ok(a) => a,
            Err(_) => break, // queue drained
        };
        let t = Instant::now();
        match fetcher.fetch(&accession) {
            Ok(mut fetched) => {
                fetched.fetch_secs = t.elapsed().as_secs_f64();
                if out.send(fetched).is_err() {
                    break; // compute pool gone
                }
            }
            Err(e) => {
                eprintln!("[{accession}] fetch failed: {e}");
                let _ = result_tx.send(AccessionResult::failed(accession, e.to_string()));
            }
        }
    }
}

/// Pull fetched buffers and run every compute stage for that accession.
fn compute_worker(
    cfg: &Config,
    query: &Query,
    query_index: &stages::anchor::QueryIndex,
    input: SharedReceiver<FetchedUnitigs>,
    out: std::sync::mpsc::Sender<AccessionResult>,
) {
    loop {
        let fetched = match input.lock().unwrap().recv() {
            Ok(f) => f,
            Err(_) => break,
        };
        let accession = fetched.accession.clone();

        // A stage returning `Err` is already handled below. A stage that *panics* (e.g. an
        // index out of bounds on a malformed graph) must not take the worker down with it: the
        // accession's buffer has already been consumed from the queue, so an uncaught panic
        // would send no result. That silently drops the accession, and if a systematic bug
        // kills every compute worker, the fetch pool then blocks forever on a full bounded
        // channel and the whole run hangs. Catching here turns a panic into one `Failed`
        // result and keeps the worker alive.
        let panic_acc = accession.clone();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            run_one(cfg, query, query_index, fetched)
                .unwrap_or_else(|e| AccessionResult::failed(accession, e.to_string()))
        }))
        .unwrap_or_else(|_| {
            AccessionResult::failed(panic_acc, "internal error: a stage panicked".to_string())
        });
        if out.send(result).is_err() {
            break;
        }
    }
}

/// Every compute stage for a single accession, in order. This is the whole per-accession
/// pipeline in one readable place.
///
/// # Memory shape
///
/// `parse::build` consumes the compressed bytes as a stream and returns everything the rest
/// of the pipeline needs, so the decompressed text (145 MB for `ERR2028820`, ~1.44 GB for the
/// largest ERG11 accessions) is never held. After the scope below closes, resident state per
/// worker is one `UnitigGraph` at ~43 MB.
///
/// Peak is therefore roughly `compressed + graph`, about 76 MB for `ERR2028820`, against
/// ~199 MB for a design that materialises the decompressed buffer. Multiply by
/// `cfg.compute_workers` for the real ceiling.
fn run_one(
    cfg: &Config,
    query: &Query,
    query_index: &stages::anchor::QueryIndex,
    fetched: FetchedUnitigs,
) -> Result<AccessionResult> {
    let t_start = Instant::now();
    let accession = fetched.accession.clone();
    let fetch_secs = fetched.fetch_secs;
    let mut timing = Timing::default();

    // Single streaming decode. `fetched` is dropped at the end of this scope; nothing
    // downstream needs the raw file again.
    let t = Instant::now();
    let graph = {
        let reader = stages::stream::decoder(&fetched)?;
        stages::parse::build(&fetched.accession, reader)?
    };
    drop(fetched);
    timing.stream = fetch_secs;
    timing.parse = t.elapsed().as_secs_f64();
    let n_unitigs = graph.n_unitigs;

    let t = Instant::now();
    let anchors = stages::anchor::find(cfg, &graph, query_index)?;
    timing.anchor = t.elapsed().as_secs_f64();

    // No anchors means the gene is not present in this accession. A legitimate outcome,
    // not a failure, and there is nothing to expand from.
    if anchors.ids.is_empty() {
        timing.total = t_start.elapsed().as_secs_f64();
        return Ok(AccessionResult {
            accession,
            outcome: Outcome::NoAnchors,
            n_unitigs,
            n_anchors: 0,
            subgraph_nodes: 0,
            haplotype: None,
            haplotypes: Vec::new(),
            alts: Vec::new(),
            bubbles: Vec::new(),
            timing,
            error: None,
        });
    }

    let t = Instant::now();
    let nodes = stages::extract::expand(cfg, &graph, &anchors)?;
    if stages::extract::exceeds_guard(cfg, &graph, &nodes) {
        eprintln!(
            "[{}] WARNING: 5-hop neighbourhood is {:.2}% of the graph ({} of {} unitigs), \
             above the {:.2}% guard. The haplotype is likely under-determined; treat this \
             accession's result as provisional.",
            graph.accession,
            stages::extract::subgraph_fraction(&graph, &nodes) * 100.0,
            nodes.len(),
            graph.n_unitigs,
            cfg.subgraph_fraction_guard * 100.0
        );
    }
    let seqs = stages::extract::collect_seqs(&graph, &nodes)?;
    let mut subgraph = stages::extract::write_gfa(cfg, &graph, nodes, seqs)?;
    timing.extract = t.elapsed().as_secs_f64();
    let subgraph_nodes = subgraph.nodes.len();

    // phase does not depend on the GAF, so it is computed regardless of whether a path was
    // found. This differs from Python, where a no_path accession reports no bubbles.
    let t = Instant::now();
    let bubbles = stages::phase::phase(cfg, &graph, &subgraph)?;
    timing.phase = t.elapsed().as_secs_f64();

    // align returns None when GraphAligner found no alignment; that is Outcome::NoPath,
    // not an error. An empty GFA making GraphAligner abort IS an error, and must not be
    // silently folded into this branch.
    let gaf_path = subgraph.gfa_path.with_extension("gaf");
    let t = Instant::now();
    let gaf = stages::align::run(cfg, &subgraph, query, &gaf_path)?;
    timing.align = t.elapsed().as_secs_f64();

    let t = Instant::now();
    // Alternative-haplotype recovery is `stages::variants`: it corrects the consensus to the
    // majority allele and emits second alleles at biallelic sites. See that module for why it
    // replaced the upstream `detour` engine (gap rate and VAF labelling).
    let mut final_path = gaf;
    let mut haplotype = match &final_path {
        Some(path) => stages::stitch::stitch(cfg, path, &subgraph, query)?,
        None => None,
    };
    timing.stitch = t.elapsed().as_secs_f64();

    // Second pass, off unless --prune-frac is set.
    //
    // The first alignment exists only to measure the gene's depth. With that number the
    // noise can be removed from the search space and the alignment redone, which is the
    // only ordering that works: the estimate is not available before an alignment, and
    // every attempt to derive it from the anchor set produced 70x against a true 14,500x.
    //
    // # What this fixes
    //
    // GraphAligner scores a route by how well it matches the *reference*. Where the sample's
    // real sequence differs from the reference and some branch carried by a handful of
    // erroneous reads happens to match it, the aligner prefers the error. Measured on
    // `SRR8209082`, whose gene sits at 14,500x, the first-pass alignment spends 118 of its
    // 1,319 unitigs below 1% of that depth, some as low as 2x. Removing sub-threshold
    // unitigs and realigning leaves a route whose *minimum* coverage is 147x:
    //
    // ```text
    //            nodes  alignments   cov     min ka:f  under 1% of depth
    // unpruned   34582      3        0.9996      2.0        118
    // 1% floor    6930      1        0.9993    147.5          2
    // 5% floor    5343      1        0.9993   7000.0          0
    // ```
    //
    // The three alignments collapsing to one is the point, not a loss: they were three ways
    // of wandering through the same noise, not three haplotypes.
    //
    // # Why it is safe on clonal samples
    //
    // The floor is a fraction of the sample's *own* depth, so a 180x isolate gets a 1.8x
    // floor and nothing is removed. Measured on `ERR2028820` at every setting from 0.01 to
    // 0.05: byte-identical haplotype, same alignment count, zero nodes dropped from the path.
    //
    // # Never a regression
    //
    // If the second pass finds no haplotype, the first pass's result stands. Pruning can
    // improve an accession and cannot break one. The cost is one extra alignment.
    if cfg.prune_frac > 0.0 {
        if let (Some(path), Some(_)) = (&final_path, &haplotype) {
            let depth = stages::extract::path_depth(&graph, &path.steps);
            let floor = stages::extract::depth_floor(cfg, depth);
            let before = subgraph.nodes.len();

            // Both subgraphs live at the same scratch path, so the first must be dropped
            // (which deletes its GFA) before the second is written.
            let nodes = std::mem::take(&mut subgraph.nodes);
            let seqs = std::mem::take(&mut subgraph.seqs);
            let (kept_nodes, kept_seqs) =
                stages::extract::prune_nodes(&graph, nodes, seqs, floor);
            let dropped = before - kept_nodes.len();

            if dropped == 0 || kept_nodes.is_empty() {
                eprintln!(
                    "[{}] prune: depth {:.0}x, floor {:.1}x, nothing to remove",
                    graph.accession, depth, floor
                );
                subgraph.nodes = kept_nodes;
                subgraph.seqs = kept_seqs;
            } else {
                drop(subgraph);
                let pruned =
                    stages::extract::write_gfa(cfg, &graph, kept_nodes, kept_seqs)?;
                let gaf2 = pruned.gfa_path.with_extension("pass2.gaf");
                let t = Instant::now();
                let second = match stages::align::run(cfg, &pruned, query, &gaf2)? {
                    Some(p) => stages::stitch::stitch(cfg, &p, &pruned, query)?.map(|h| (p, h)),
                    None => None,
                };
                timing.align += t.elapsed().as_secs_f64();
                match second {
                    Some((p, h)) => {
                        eprintln!(
                            "[{}] prune: depth {:.0}x, floor {:.1}x, {} of {} unitigs kept; \
                             realigned, min coverage on the route now {:.1}x",
                            graph.accession, depth, floor, before - dropped, before,
                            stages::extract::path_depth(&graph, &p.steps).min(
                                p.steps.iter()
                                    .map(|s| graph.abundance[s.node as usize])
                                    .fold(f32::INFINITY, f32::min)
                            )
                        );
                        haplotype = Some(h);
                        final_path = Some(p);
                        subgraph = pruned;
                    }
                    None => {
                        eprintln!(
                            "[{}] prune: depth {:.0}x, floor {:.1}x cost the alignment; \
                             keeping the unpruned result",
                            graph.accession, depth, floor
                        );
                        // Redo pass 1's subgraph so later stages have it. Cheap: the
                        // sequences are already decoded, only the GFA is rewritten.
                        drop(pruned);
                        let nodes = stages::extract::expand(cfg, &graph, &anchors)?;
                        let seqs = stages::extract::collect_seqs(&graph, &nodes)?;
                        subgraph = stages::extract::write_gfa(cfg, &graph, nodes, seqs)?;
                    }
                }
            }
        }
    }

    // Upstream's explicit scratch cleanup, applied to whichever subgraph ended up active:
    // the fallback replaces `subgraph` above when pruning cost us the haplotype.
    if std::env::var_os("LOGAN_KEEP_SCRATCH").is_none() {
        let _ = std::fs::remove_file(&subgraph.gfa_path);
    }
    // Intra-host variant recovery. `prune_frac` is not part of this: it cleans the graph
    // before this runs. `prefer_majority` and the alt/vcf/snv outputs each gate themselves, so
    // with all of them off this reduces to emitting the single stitched primary.
    let t = Instant::now();
    let mut alts: Vec<crate::types::AltHaplotype> = Vec::new();
    // Correct the consensus first, so the alternatives that follow are the minority alleles
    // relative to a corrected route rather than to a reference-biased one.
    if let (Some(path), Some(_)) = (&final_path, &haplotype) {
        if let Some(fixed) = stages::variants::prefer_majority(cfg, &graph, &subgraph, path, query)?
        {
            if let Some(h) = stages::stitch::stitch(cfg, &fixed, &subgraph, query)? {
                haplotype = Some(h);
                final_path = Some(fixed);
            }
        }
    }
    if let (Some(path), Some(_)) = (&final_path, &haplotype) {
        alts = stages::variants::find_alts(cfg, &graph, &subgraph, path, query)?;
    }
    let haplotypes: Vec<crate::types::Haplotype> = haplotype.clone().into_iter().collect();
    timing.stitch += t.elapsed().as_secs_f64();
    timing.total = t_start.elapsed().as_secs_f64();

    // Matches `logan_path_walker.py:312`: a path that failed the coverage gate is no_path.
    let haplotype = haplotypes.first().cloned();
    let outcome = if !haplotypes.is_empty() {
        Outcome::Success
    } else {
        Outcome::NoPath
    };

    Ok(AccessionResult {
        accession,
        outcome,
        alts,
        n_unitigs,
        n_anchors: anchors.ids.len(),
        subgraph_nodes,
        haplotype,
        haplotypes,
        bubbles,
        timing,
        error: None,
    })
}
