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
//!  │        sync_channel(cfg.channel_capacity)   ◄── BOUNDED, but the     │
//!  │                  │                              ceiling is NOT just  │
//!  │                  │                              capacity: see below. │
//!  │                  ▼                                                   │
//!  │  COMPUTE POOL  cfg.compute_workers (default 8, NOT auto-sized)       │
//!  │   ┌────────┐ ┌────────┐        ┌────────┐                            │
//!  │   │ work 0 │ │ work 1 │  ...   │ work W │   run_one() below          │
//!  │   └───┬────┘ └───┬────┘        └───┬────┘                            │
//!  │       │          │                 │                                 │
//!  │       │   per accession, in order:                                   │
//!  │       │     stages::parse    -> UnitigGraph (CSR, ~43 MB)            │
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
//! # Why the channel is bounded, and what the bound actually is
//!
//! Fetch threads block once the queue is full, so a fast network cannot outrun slow compute
//! and pile up buffers without limit. The Python `ThreadPoolExecutor` had no such limit,
//! which is how a single run reached 21.3 GB and took the machine down with it.
//!
//! This file used to claim the ceiling was "capacity times the largest retained buffer".
//! **That is wrong, by a factor of five at the defaults.** A compressed buffer is resident in
//! three places at once, not one:
//!
//! ```text
//!   cfg.channel_capacity     queued in the bounded channel             4
//! + cfg.fetch_workers        held by a worker blocked in send(),       8
//!                            after read_to_vec has already returned
//! + cfg.compute_workers      held by run_one for the whole parse,      8
//!                            dropped only once the CSR graph exists
//!                            -----------------------------------------
//!                            peak concurrent buffers                  20
//! ```
//!
//! `Config::peak_compressed_buffers` computes that and
//! `Config::max_compressed_resident_bytes` multiplies it by `--max-size-mb`. At the 1 GB
//! default cap the worst case is therefore **about 20 GB**, not 4 GB, and the binary prints
//! the figure at startup so it is visible before the run rather than after the OOM.
//!
//! On top of that sit one `UnitigGraph` per compute worker (~43 MB each) and whatever
//! `read_to_vec` transiently allocates while its `Vec` grows.
//!
//! The bound is worst case, not typical: it assumes every accession sits at the cap, while
//! `ERR2028820` is 33 MB and the largest of the ERG11 panel is 343 MB. But the cap is what
//! bounds the run, and raising `--max-size-mb`, `-w` or `--fetch-workers` moves all three
//! terms at once. That 343 MB panel is already ~6.9 GB in compressed buffers alone.
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

/// What to do with a depth-based prune, given how many unitigs would survive it.
///
/// Three distinct cases. The one that used to be folded in with "nothing to remove" is
/// `WouldEmpty`: a floor above every unitig in the subgraph removes the entire search space,
/// which is not a prune but a silent deletion, and it left the accession with no majority
/// correction, no alts, no VCF rows and no SNV rows under a log line saying nothing happened.
#[derive(Debug, PartialEq, Eq)]
enum PruneOutcome {
    /// Every unitig is below the floor. Skip the prune; an unpruned result beats none.
    WouldEmpty,
    /// Every unitig clears the floor. Nothing to do.
    NothingToRemove,
    /// A proper subset survives: rewrite the GFA and realign.
    Prune,
}

fn prune_outcome(surviving: usize, before: usize) -> PruneOutcome {
    if surviving == 0 {
        PruneOutcome::WouldEmpty
    } else if surviving == before {
        PruneOutcome::NothingToRemove
    } else {
        PruneOutcome::Prune
    }
}

/// --fetch-only variant: download every accession into the unitigs cache, no compute stage.
///
/// Returns `(ok, fail)`. Uses only the fetch pool; the compute pool is not spawned. Each
/// fetch's success side is a no-op beyond writing the cache entry (which `Fetcher::fetch`
/// already does), and each failure is reported inline. Bounded parallelism is
/// `cfg.fetch_workers`.
///
/// Requires `cfg.unitigs_cache` — validated by the caller (main.rs).
pub fn run_fetch_only(
    cfg: &Config,
    accessions: Vec<Accession>,
) -> Result<(usize, usize)> {
    let total = accessions.len();
    let work: SharedReceiver<(usize, Accession)> = {
        let (tx, rx) = std::sync::mpsc::channel();
        for (i, a) in accessions.into_iter().enumerate() {
            let _ = tx.send((i, a));
        }
        Arc::new(Mutex::new(rx))
    };

    let fetcher = stages::stream::Fetcher::new(cfg);
    let fetcher = &fetcher;
    let done = Arc::new(Mutex::new((0usize, 0usize))); // (ok, fail)

    std::thread::scope(|scope| {
        for _ in 0..cfg.fetch_workers {
            let work = Arc::clone(&work);
            let done = Arc::clone(&done);
            scope.spawn(move || loop {
                let (index, accession) = match work.lock().unwrap().recv() {
                    Ok(a) => a,
                    Err(_) => break,
                };
                let t = Instant::now();
                match fetcher.fetch(index, &accession) {
                    Ok(f) => {
                        let mut d = done.lock().unwrap();
                        d.0 += 1;
                        let progress = d.0 + d.1;
                        println!(
                            "[{progress}/{total}] {accession}: cached ({} bytes in {:.1}s)",
                            f.compressed_len,
                            t.elapsed().as_secs_f64()
                        );
                    }
                    Err(e) => {
                        let mut d = done.lock().unwrap();
                        d.1 += 1;
                        let progress = d.0 + d.1;
                        eprintln!("[{progress}/{total}] {accession}: FETCH FAILED: {e}");
                    }
                }
            });
        }
    });

    let (ok, fail) = *done.lock().unwrap();
    Ok((ok, fail))
}

/// Fan accessions through the two pools and collect results.
pub fn run(cfg: &Config, accessions: Vec<Accession>, query: &Query) -> Result<Summary> {
    let total = accessions.len();
    let (fetch_tx, fetch_rx) = sync_channel::<FetchedUnitigs>(cfg.channel_capacity);
    let (result_tx, result_rx) = std::sync::mpsc::channel::<AccessionResult>();
    let fetch_rx: SharedReceiver<FetchedUnitigs> = Arc::new(Mutex::new(fetch_rx));
    let work: SharedReceiver<(usize, Accession)> = {
        let (tx, rx) = std::sync::mpsc::channel();
        for (i, a) in accessions.into_iter().enumerate() {
            let _ = tx.send((i, a));
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

    // Progress is reported as each accession lands, which is what the Python walker did and
    // what makes a stalled run visible. The *stored* order is the input order, so every
    // output file is byte-reproducible across runs and `--sieve`'s representative is chosen
    // by the list the user supplied rather than by which S3 fetch happened to win.
    let mut results = results;
    results.sort_by_key(|r| r.index);

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
    work: SharedReceiver<(usize, Accession)>,
    out: SyncSender<FetchedUnitigs>,
    result_tx: std::sync::mpsc::Sender<AccessionResult>,
) {
    loop {
        let (index, accession) = match work.lock().unwrap().recv() {
            Ok(a) => a,
            Err(_) => break, // queue drained
        };
        let t = Instant::now();
        match fetcher.fetch(index, &accession) {
            Ok(mut fetched) => {
                fetched.fetch_secs = t.elapsed().as_secs_f64();
                if out.send(fetched).is_err() {
                    break; // compute pool gone
                }
            }
            Err(e) => {
                eprintln!("[{accession}] fetch failed: {e}");
                let _ = result_tx.send(AccessionResult::failed(index, accession, e.to_string()));
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
        let index = fetched.index;

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
                .unwrap_or_else(|e| AccessionResult::failed(index, accession, e.to_string()))
        }))
        .unwrap_or_else(|_| {
            AccessionResult::failed(
                index,
                panic_acc,
                "internal error: a stage panicked".to_string(),
            )
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
/// ~199 MB for a design that materialises the decompressed buffer.
///
/// For the whole-process ceiling, multiplying this by `cfg.compute_workers` is **not**
/// enough: compressed buffers are also resident in the channel and in blocked fetch workers.
/// See the module doc and `Config::peak_compressed_buffers`.
fn run_one(
    cfg: &Config,
    query: &Query,
    query_index: &stages::anchor::QueryIndex,
    fetched: FetchedUnitigs,
) -> Result<AccessionResult> {
    let t_start = Instant::now();
    let accession = fetched.accession.clone();
    // Captured before `fetched` is dropped below; every result this function returns carries
    // it so the collector can restore input order.
    let index = fetched.index;
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
            index,
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
        Some(path) => stages::stitch::stitch(cfg, &graph, path, &subgraph, query)?,
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

            // Count survivors before consuming anything. `prune_nodes` takes ownership of the
            // node and sequence vectors, so a prune that turns out to be unusable cannot be
            // undone afterwards.
            //
            // Three outcomes, and the old code collapsed two of them. It tested
            // `dropped == 0 || kept_nodes.is_empty()` and, in the second case, assigned the
            // EMPTY vectors back while logging "nothing to remove". Splitting the empty case
            // out keeps the unpruned subgraph in that branch, which is the same "never a
            // regression" rule the failed-realign path applies.
            //
            // Predicate must match `prune_nodes` exactly.
            let surviving = subgraph
                .nodes
                .iter()
                .filter(|&&id| graph.abundance[id as usize] >= floor)
                .count();

            match prune_outcome(surviving, before) {
                PruneOutcome::WouldEmpty => {
                // Pruning to nothing is not a prune. Keep the unpruned subgraph, which is the
                // same "never a regression" rule the failed-realign branch below applies.
                eprintln!(
                    "[{}] prune: depth {:.0}x, floor {:.1}x would remove all {} unitigs; \
                     skipping the prune and keeping the unpruned result",
                    graph.accession, depth, floor, before
                );
                }
                PruneOutcome::NothingToRemove => {
                eprintln!(
                    "[{}] prune: depth {:.0}x, floor {:.1}x, nothing to remove",
                    graph.accession, depth, floor
                );
                }
                PruneOutcome::Prune => {
                // Both subgraphs live at the same scratch path, so the first must be dropped
                // (which deletes its GFA) before the second is written.
                let nodes = std::mem::take(&mut subgraph.nodes);
                let seqs = std::mem::take(&mut subgraph.seqs);
                let (kept_nodes, kept_seqs) =
                    stages::extract::prune_nodes(&graph, nodes, seqs, floor);
                let dropped = before - kept_nodes.len();
                drop(subgraph);
                let pruned =
                    stages::extract::write_gfa(cfg, &graph, kept_nodes, kept_seqs)?;
                let gaf2 = pruned.gfa_path.with_extension("pass2.gaf");
                let t = Instant::now();
                let second = match stages::align::run(cfg, &pruned, query, &gaf2)? {
                    Some(p) => stages::stitch::stitch(cfg, &graph, &p, &pruned, query)?.map(|h| (p, h)),
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
            if let Some(mut h) = stages::stitch::stitch(cfg, &graph, &fixed, &subgraph, query)? {
                // The route was spliced onto the majority allele, so the GAF's identity
                // describes the sequence we just replaced, not this one.
                h.pident_inherited = true;
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
        index,
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

#[cfg(test)]
mod ordering_tests {
    use crate::types::{AccessionResult, Outcome};

    /// Results arrive in completion order, decided by S3 and by how long each graph takes.
    /// `run` sorts on `index` before returning so every output file is byte-reproducible and
    /// `--sieve`'s representative is the first record of the *input*, not of a race.
    #[test]
    fn sorting_on_index_restores_input_order() {
        // The measured completion order from the nine-accession Zika run in the review.
        let input = [
            "SRR6251963", "SRR6251966", "SRR6251967", "SRR6251969", "SRR6251973",
            "SRR6251979", "SRR6251984", "SRR6251992", "SRR6251993",
        ];
        let completion = [5usize, 7, 4, 1, 6, 8, 3, 0, 2];

        let mut results: Vec<AccessionResult> = completion
            .iter()
            .map(|&i| AccessionResult::failed(i, input[i].to_string(), "x".into()))
            .collect();
        assert_ne!(
            results.iter().map(|r| r.index).collect::<Vec<_>>(),
            (0..9).collect::<Vec<_>>(),
            "fixture must start scrambled or it proves nothing"
        );

        results.sort_by_key(|r| r.index);

        let got: Vec<&str> = results.iter().map(|r| r.accession.as_str()).collect();
        assert_eq!(got, input, "stored order must be the input order");
    }

    /// Every failure path constructs its own result, so the index has to survive all of them.
    #[test]
    fn a_failed_result_keeps_its_index() {
        let r = AccessionResult::failed(41, "SRR8209082".into(), "fetch failed".into());
        assert_eq!(r.index, 41);
        assert_eq!(r.outcome, Outcome::Failed);
    }
}

#[cfg(test)]
mod prune_outcome_tests {
    use super::{prune_outcome, PruneOutcome};
    use crate::config::Config;
    use crate::stages::extract::{depth_floor, prune_nodes};
    use crate::stages::parse;
    use crate::types::UnitigGraph;

    fn graph(fa: &str) -> UnitigGraph {
        parse::build(&"T".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec()))).unwrap()
    }

    /// The three cases must stay distinct. Folding `WouldEmpty` back into
    /// `NothingToRemove` is the defect this replaced.
    #[test]
    fn the_three_cases_are_distinct() {
        assert_eq!(prune_outcome(0, 10), PruneOutcome::WouldEmpty);
        assert_eq!(prune_outcome(10, 10), PruneOutcome::NothingToRemove);
        assert_eq!(prune_outcome(3, 10), PruneOutcome::Prune);
        // A one-node subgraph that survives is "nothing to remove", not a prune.
        assert_eq!(prune_outcome(1, 1), PruneOutcome::NothingToRemove);
        // An empty subgraph never reaches a prune, but must not claim one.
        assert_eq!(prune_outcome(0, 0), PruneOutcome::WouldEmpty);
    }

    /// The case is reachable on real data: a gene sitting entirely under `min_abund`.
    /// `prune_nodes` returns nothing, which is exactly why the caller must not apply it.
    #[test]
    fn a_gene_below_min_abund_would_empty_the_subgraph() {
        let g = graph(">T_0 ka:f:2.0\nACGT\n>T_1 ka:f:3.0\nACGT\n>T_2 ka:f:1.0\nACGT\n");
        let cfg = Config::default(); // min_abund 5.0
        let floor = depth_floor(&cfg, 2.0);
        assert!(floor >= 5.0, "floor is clamped to min_abund, got {floor}");

        let nodes = vec![0u32, 1, 2];
        let surviving = nodes
            .iter()
            .filter(|&&id| g.abundance[id as usize] >= floor)
            .count();
        assert_eq!(surviving, 0);
        assert_eq!(prune_outcome(surviving, nodes.len()), PruneOutcome::WouldEmpty);

        // And the predicate above matches what prune_nodes would actually do.
        let seqs = vec![b"ACGT".to_vec(); 3];
        let (kept, _) = prune_nodes(&g, nodes, seqs, floor);
        assert!(kept.is_empty(), "prune_nodes agrees the subgraph would be emptied");
    }
}
