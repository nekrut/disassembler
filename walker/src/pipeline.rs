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
//!  │       │     stages::parse    -> UnitigGraph (CSR, up to 1.7x compr.) │
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
//! Compute is CPU-bound and wants one worker per core. Conflating them behind a single
//! `--workers` makes `--workers 16` destructive: it multiplies memory 16x while buying no
//! extra bandwidth.
//!
//! # Why the channel is bounded, and what the bound actually is
//!
//! Fetch threads block once the queue is full, so a fast network cannot outrun slow compute
//! and pile up buffers without limit. An unbounded thread pool has no such limit, which is
//! how a single run reached 21.3 GB and took the machine down with it.
//!
//! This file used to claim the ceiling was "capacity times the largest retained buffer".
//! **That is wrong, by a factor of five at the defaults.** A compressed buffer is resident in
//! three places at once, not one:
//!
//! ```text
//!   cfg.channel_capacity     queued in the bounded channel             0
//! + cfg.fetch_workers        held by a worker blocked in send(),       8
//!                            after read_to_vec has already returned
//! + cfg.compute_workers      held by run_one for the whole parse,      8
//!                            dropped only once the CSR graph exists
//!                            -----------------------------------------
//!                            peak concurrent buffers                  16
//!
//!   capped by the accession count: a run of n accessions cannot hold
//!   more than n buffers or more than min(n, W) graphs, and --fetch-only
//!   holds F buffers and no graphs at all
//! ```
//!
//! `Config::peak_compressed_buffers_for` computes that and
//! `Config::max_compressed_resident_bytes_for` multiplies it by `--max-size-mb`: at the 1 GB
//! default cap, **about 16 GB** in compressed buffers. It was 20 GB until the queue depth
//! went to 0, which is the whole of what a positive depth was costing.
//!
//! On top of that sits one `UnitigGraph` per compute worker, which peaks at up to 1.7x the
//! compressed size of its accession (`Config::GRAPH_PEAK_PER_COMPRESSED_BYTE`; it was taken
//! to be a fixed ~43 MB), and whatever `read_to_vec` transiently allocates while its `Vec`
//! grows. `Config::max_resident_bytes_for` adds the graphs: **about 29.6 GB** at the
//! defaults. The binary prints it at startup for the run it was actually given, so it is
//! visible before the run rather than after the OOM.
//!
//! That 1.7x is a calibration, not an enforced bound, and it was measured one accession per
//! run. Across many accessions glibc's dynamic mmap threshold rises as freed graph arrays are
//! returned, so later arrays grow inside the heap and the steady state drifts well above it.
//! `main::pin_mmap_threshold` pins the threshold at 1 MiB to stop that.
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
use crate::types::{
    Accession, AccessionResult, FetchedUnitigs, GafPath, Outcome, PruneReport, PruneResult, Query,
    Timing,
};
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

/// Whether a pruned realignment may replace the unpruned one: only if it spans at least as
/// much of the reference.
///
/// Removing a low-depth stretch can cut the gene's route in two. GraphAligner then reports the
/// pieces as separate alignments, the one with the most matches wins, and the haplotype loses
/// the other piece while still passing `--min-cov`. Measured on `ERR6343483` (SARS-CoV-2
/// spike, 1,700x): the unpruned alignment spans reference bases 0-3821; with the 17x floor it
/// splits into 0-166 and 424-3821, and the primary went from 3,816 nt to 3,391.
///
/// The span is in reference coordinates, so the indels pruning exists to correct do not count
/// against it: removing a spurious insertion or recovering a real deletion changes the
/// haplotype's length, not the stretch of reference it covers.
///
/// `PRUNE_SPAN_TOLERANCE` keeps the rule from firing on a single terminal codon.
fn pruned_covers_as_much(unpruned: &GafPath, pruned: &GafPath) -> bool {
    pruned.query_span() + PRUNE_SPAN_TOLERANCE >= unpruned.query_span()
}

/// Whether the pruned haplotype is the unpruned one with bases trimmed off an end (WA-09).
///
/// `pruned_covers_as_much` allows a pruned realignment to give up `PRUNE_SPAN_TOLERANCE`
/// reference bases, because a strict rule kept a reference-biased primary over one terminal
/// codon. But that tolerance was the whole test: a realignment that changed *nothing else*
/// was accepted and logged as `applied`, which is pure loss. The maintainers' own count says
/// it happens: of 511 refusals in a 20,000-accession run, 33 gave up three reference bases or
/// fewer, 21 of them exactly one, and a lost terminal base can be a truncated stop codon.
///
/// So the tolerance now buys a shorter result only when the prune actually changed the
/// haplotype. A strict substring test, not an end-trim test: a prune that removed an interior
/// unitig cannot leave the rest contiguous and identical, so anything that is still a
/// substring of the first pass is trimming and nothing else.
fn trimmed_copy(pruned: &[u8], unpruned: &[u8]) -> bool {
    pruned.len() < unpruned.len()
        && !pruned.is_empty()
        && unpruned.windows(pruned.len()).any(|w| w == pruned)
}

/// Reference bases a pruned realignment may lose and still replace the unpruned one.
///
/// The rule above was strict through 0.1.0, which meant one terminal base was enough to keep a
/// reference-biased primary. The 2026-11-01 reruns showed that is not hypothetical: of the
/// refusals, 33 of 511 in the 20,000-accession SARS-CoV-2 run and 21 of 56 across the obelisk
/// cohorts gave up three reference bases or fewer, 21 of them exactly one. The losses this rule
/// exists to catch are nothing like that size: 421 of those 511 lost more than 30 bases, and the
/// case it was written for lost 424.
///
/// Bases at the very end are also the ones a sample is most likely to genuinely disagree with
/// the reference on, so a strict rule trades a real correction for a codon of reference bias. A
/// tolerance this small cannot admit a split alignment.
const PRUNE_SPAN_TOLERANCE: usize = 3;

/// --fetch-only variant: download every accession into the unitigs cache, no compute stage.
///
/// Returns `(ok, negatives, errors)`: cached, not cached because Logan has no data or the
/// accession is over the size cap, and not cached because of an error. Uses only the fetch
/// pool; the compute pool is not spawned. Each
/// fetch's success side is a no-op beyond writing the cache entry (which `Fetcher::fetch`
/// already does), and each failure is reported inline. Bounded parallelism is
/// `cfg.fetch_workers`.
///
/// Requires `cfg.unitigs_cache` — validated by the caller (main.rs).
pub fn run_fetch_only(
    cfg: &Config,
    accessions: Vec<Accession>,
) -> Result<Summary> {
    let total = accessions.len();
    let work: SharedReceiver<(usize, Accession)> = {
        let (tx, rx) = std::sync::mpsc::channel();
        for (i, a) in accessions.into_iter().enumerate() {
            let _ = tx.send((i, a));
        }
        Arc::new(Mutex::new(rx))
    };

    let fetcher = stages::stream::Fetcher::new_for(cfg, total);
    let fetcher = &fetcher;
    // Results, not counts. A fetch-only run used to report three integers to stderr and exit,
    // so the frozen-input route (`disassembler fetch` into a cache, then `--offline`) left no
    // record of what it had frozen: the one step whose whole purpose is reproducibility was
    // the only one that wrote no summary.
    let collected: Arc<Mutex<Vec<AccessionResult>>> = Arc::new(Mutex::new(Vec::new()));

    std::thread::scope(|scope| {
        for _ in 0..cfg.fetch_workers {
            let work = Arc::clone(&work);
            let collected = Arc::clone(&collected);
            scope.spawn(move || loop {
                let (index, accession) = match work.lock().unwrap().recv() {
                    Ok(a) => a,
                    Err(_) => break,
                };
                let t = Instant::now();
                // WP-18, as in `fetch_worker`: a panic here would end the run rather than
                // the accession, and `--fetch-only` is the step whose whole purpose is to
                // leave a usable cache behind.
                let fetched = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    fetcher.fetch(index, &accession)
                })) {
                    Ok(r) => r,
                    Err(p) => {
                        let msg =
                            format!("internal error: the fetch panicked: {}", panic_message(&*p));
                        eprintln!("[{accession}] {msg}");
                        let mut out = collected.lock().unwrap();
                        out.push(AccessionResult::failed(index, accession, msg));
                        continue;
                    }
                };
                // The lock is taken once and held across the report, which also stops two
                // workers interleaving a line.
                let mut out = collected.lock().unwrap();
                let progress = out.len() + 1;
                match fetched {
                    Ok(f) => {
                        crate::outln!(
                            "[{progress}/{total}] {accession}: cached ({} bytes in {:.1}s)",
                            f.compressed_len,
                            t.elapsed().as_secs_f64()
                        );
                        out.push(cached(index, accession, f.payload));
                    }
                    Err(e) => {
                        let r = result_for_error(index, accession.clone(), e.to_string(), &e);
                        if r.outcome.is_negative() {
                            eprintln!("[{progress}/{total}] {accession}: not cached: {e}");
                        } else {
                            eprintln!("[{progress}/{total}] {accession}: FETCH FAILED: {e}");
                        }
                        out.push(r);
                    }
                }
            });
        }
    });

    let mut results = std::mem::take(&mut *collected.lock().unwrap());
    restore_input_order(&mut results);
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

/// The result for an accession whose bytes are now in the cache.
///
/// `Outcome::Success` means the unit of work succeeded, which under `--fetch-only` is the
/// caching, not a reconstruction. The summary says which kind of run it was: `--fetch-only`
/// writes `inputs.references: null` and an empty `tool.external`, and `cohort.reconstructed`
/// is 0 because nothing was attempted.
fn cached(index: usize, accession: Accession, payload: crate::types::Payload) -> AccessionResult {
    let mut r = AccessionResult::empty(index, accession, Outcome::Success, None);
    r.payload = payload;
    r
}

/// What a caught panic actually said.
///
/// WP-18. `catch_unwind` hands back `Box<dyn Any>`, and the compute path used to throw it
/// away and report "a stage panicked", which says nothing about where. A panic payload is a
/// `&str` or a `String` in every case the standard macros produce.
fn panic_message(p: &(dyn std::any::Any + Send)) -> &str {
    p.downcast_ref::<&str>()
        .copied()
        .or_else(|| p.downcast_ref::<String>().map(String::as_str))
        .unwrap_or("non-string panic payload")
}

/// Turn a per-accession error into its result, with the right outcome.
///
/// Shared by both pools, because both can produce the same negatives. The fetch pool raises
/// `TooLarge` from a `Content-Length` or a cached file's size; the compute pool raises it from
/// `stages::parse`, when a graph exceeds the per-accession budget the startup line promised.
/// Routing the second through `AccessionResult::failed` would have reported "this run broke"
/// for what is really "this accession does not fit the cap you set", and `--allow-empty`
/// treats the two differently.
fn result_for_error(
    index: usize,
    accession: Accession,
    detail: String,
    e: &crate::error::Error,
) -> AccessionResult {
    use crate::error::Error;
    match e {
        Error::NoLoganData { .. } => AccessionResult::empty(index, accession, Outcome::NoLoganData, None),
        // The bare detail, which now has to carry the reason as well as the size: `too_large`
        // is no longer only `--max-size-mb`. It is also a graph over the parse budget, an
        // accession whose spool will not fit the spool directory, and a spool write that ran
        // out of space. Every site writes a self-describing detail, and the report line below
        // prints it rather than naming a flag that may not be the one at fault.
        //
        // The payload is what makes the skip actionable: it carries the size that tripped the
        // cap as a number, where before it survived only as prose inside `detail`.
        Error::TooLarge { detail: d, payload, .. } => {
            let mut r = AccessionResult::empty(index, accession, Outcome::TooLarge, Some(d.clone()));
            r.payload = (**payload).clone();
            r
        }
        _ => AccessionResult::failed(index, accession, detail),
    }
}

/// The variant stages' config when the prune asked for a floor it never got to apply.
///
/// `Some(raised)` only when the prune was refused (no haplotype, a shorter span, or it would
/// have emptied the subgraph) AND its floor is above `--min-abund`. In every other case the
/// variant stages run under the caller's config unchanged, so applied and pruned-off runs are
/// byte-identical to before. See the call site for what this costs when it is missing.
fn inherited_floor_config(cfg: &Config, prune: Option<&PruneReport>) -> Option<Config> {
    let r = prune?;
    let refused = !matches!(
        r.result,
        crate::types::PruneResult::Applied | crate::types::PruneResult::NothingToRemove
    );
    (refused && r.floor > cfg.min_abund)
        .then(|| Config { min_abund: r.floor, ..cfg.clone() })
}

/// Put results back into the order the user listed the accessions in.
///
/// Results arrive in completion order, which S3 and graph size decide, so without this every
/// output file is ordered differently from run to run: `summary.json`, the FASTA, the VCFs
/// and `snvs.tsv` all follow it, and `--sieve` keeps whichever duplicate finished first.
///
/// TW-05: this is a named function so that a unit test can exercise the real thing. The test
/// that claimed to cover it sorted a vector of its own and never called into the pipeline, so
/// it passed with the production sort deleted, as did every other test.
fn restore_input_order(results: &mut [AccessionResult]) {
    results.sort_by_key(|r| r.index);
}

/// Fan accessions through the two pools and collect results.
///
/// `queries` is the reference pool. In today's single-anchor path it holds one [`Query`];
/// multi-anchor mode passes several and the walker unions their canonical k-mers into one
/// anchor index. The downstream stages are handed the first [`Query`] for now — later phases
/// of the multi-anchor feature (align/stitch/pick_primary) fan out over the whole pool.
pub fn run(
    cfg: &Config,
    accessions: Vec<Accession>,
    queries: &[Query],
    journal: Option<&crate::io::output::Journal>,
) -> Result<Summary> {
    if queries.is_empty() {
        return Err(crate::error::Error::Usage("no queries supplied".into()));
    }

    let total = accessions.len();
    let (fetch_tx, fetch_rx) = sync_channel::<FetchedUnitigs>(cfg.channel_capacity);
    let (result_tx, result_rx) = std::sync::mpsc::channel::<AccessionResult>();
    let fetch_rx: SharedReceiver<FetchedUnitigs> = Arc::new(Mutex::new(fetch_rx));
    // Kept so an accession that comes back without a result can still be named below.
    let names = accessions.clone();
    let work: SharedReceiver<(usize, Accession)> = {
        let (tx, rx) = std::sync::mpsc::channel();
        for (i, a) in accessions.into_iter().enumerate() {
            let _ = tx.send((i, a));
        }
        Arc::new(Mutex::new(rx))
    };

    // Built once per run, not per accession: it depends only on cfg and the anchor pool.
    // Multi-anchor mode's whole point lives here — the anchor index is the union of every
    // reference's canonical k-mers, so a graph carrying an Omicron-like spike is anchor-hit
    // even under a Wuhan-only invocation of downstream stages (until Phase 4 fans them out).
    let query_index = stages::anchor::index_query_set(cfg, queries)?;
    let query_index = &query_index;

    // One pooled HTTP client shared by every fetch thread. The bucket is a single host,
    // so per-thread agents would just multiply TLS handshakes.
    let fetcher = stages::stream::Fetcher::new_for(cfg, total);
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
            scope.spawn(move || compute_worker(cfg, queries, query_index, rx, tx));
        }
        drop(result_tx);
        // Only the compute workers may hold the receiving end. While this thread kept a
        // handle, a compute pool that stopped draining (zero workers, or every worker gone
        // after the collector below panicked on a closed stdout) left the fetch pool blocked
        // on a full channel forever, and the scope never joined: the run hung instead of
        // ending. With the last receiver gone, a blocked `send` fails and the fetcher stops.
        drop(fetch_rx);

        // ---- collector (this thread) ------------------------------------
        // Report each accession as it lands, one line each: a run that quietly reconstructs
        // nothing looks identical to a successful one otherwise, which is exactly how the
        // int/str regression hid.
        let mut done = 0usize;
        result_rx
            .into_iter()
            .inspect(|r| {
                done += 1;
                report(done, total, r);
                // WP-06: and to disk, as it lands. Everything else this run writes waits for
                // the last accession, so without this a kill at hour eight of a nine-hour
                // shard leaves nothing at all.
                if let Some(j) = journal {
                    j.record(r);
                }
            })
            .collect::<Vec<_>>()
    });

    // Progress is reported as each accession lands, which is what makes a stalled run
    // visible. The *stored* order is the input order, so every
    // output file is byte-reproducible across runs and `--sieve`'s representative is chosen
    // by the list the user supplied rather than by which S3 fetch happened to win.
    let mut results = results;
    restore_input_order(&mut results);

    // Every accession must come back with a result, or it vanishes from every output and
    // from the success count. The pools are built so this cannot happen; if it ever does,
    // the accession is reported as an error rather than dropped.
    if results.len() < total {
        let seen: std::collections::HashSet<usize> = results.iter().map(|r| r.index).collect();
        for (i, name) in names.into_iter().enumerate() {
            if !seen.contains(&i) {
                eprintln!("[{name}] internal error: no result was produced");
                results.push(AccessionResult::failed(
                    i,
                    name,
                    "internal error: no result was produced".to_string(),
                ));
            }
        }
        restore_input_order(&mut results);
    }

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

/// One line per finished accession.
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
            crate::outln!(
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
        Outcome::NoAnchors => crate::outln!("[{done}/{total}] {}: no_anchors", r.accession),
        Outcome::NoPath => crate::outln!("[{done}/{total}] {}: no_path", r.accession),
        Outcome::NoLinks => crate::outln!(
            "[{done}/{total}] {}: no_links (the unitig file carries no L: adjacency)",
            r.accession
        ),
        Outcome::NoLoganData => crate::outln!(
            "[{done}/{total}] {}: no_logan_data (no Logan unitigs for this accession)",
            r.accession
        ),
        Outcome::SubgraphTooLarge => crate::outln!(
            "[{done}/{total}] {}: subgraph_too_large{}",
            r.accession,
            r.error.as_deref().map(|e| format!(" ({e})")).unwrap_or_default()
        ),
        Outcome::TooLarge => crate::outln!(
            "[{done}/{total}] {}: too_large, skipped ({})",
            r.accession,
            r.error.as_deref().unwrap_or("over the cap")
        ),
        Outcome::Failed => crate::outln!(
            "[{done}/{total}] {}: ERROR: {}",
            r.accession,
            r.error.as_deref().unwrap_or("unknown")
        ),
    }
}

/// Pull accessions, stream them, hand the bytes downstream. Blocks when the queue is full,
/// which is the backpressure that bounds memory.
/// What `fetch_worker` needs of a fetcher.
///
/// WP-18 put a `catch_unwind` around the call; this is how that is testable. The production
/// implementation is [`stages::stream::Fetcher`] and there is no other, so the trait buys
/// nothing at run time: it exists so a test can hand the worker a fetcher that panics and
/// watch the run survive it.
pub trait Fetch {
    fn fetch(&self, index: usize, accession: &Accession) -> Result<FetchedUnitigs>;
}

impl Fetch for stages::stream::Fetcher {
    fn fetch(&self, index: usize, accession: &Accession) -> Result<FetchedUnitigs> {
        stages::stream::Fetcher::fetch(self, index, accession)
    }
}

fn fetch_worker<F: Fetch + ?Sized>(
    fetcher: &F,
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
        // WP-18. The compute pool has caught its panics since TW-01; this one did not, and a
        // panic here does not fail one accession, it unwinds the scope and takes the whole
        // run with it: exit 101, no outputs, and the collector's "no result was produced"
        // fallback never runs. No reachable panic was found in the shipped fetch path, so
        // this is defence in depth for the pool that touches the network, zstd and the
        // filesystem.
        let fetched = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            fetcher.fetch(index, &accession)
        })) {
            Ok(r) => r,
            Err(p) => {
                let msg = format!("internal error: the fetch panicked: {}", panic_message(&*p));
                eprintln!("[{accession}] {msg}");
                let _ = result_tx.send(AccessionResult::failed(index, accession, msg));
                continue;
            }
        };
        match fetched {
            Ok(mut fetched) => {
                fetched.fetch_secs = t.elapsed().as_secs_f64();
                if out.send(fetched).is_err() {
                    break; // compute pool gone
                }
            }
            Err(e) => {
                let result = result_for_error(index, accession, e.to_string(), &e);
                if result.outcome == Outcome::Failed {
                    eprintln!("{e}");
                }
                let _ = result_tx.send(result);
            }
        }
    }
}

/// Pull fetched buffers and run every compute stage for that accession.
fn compute_worker(
    cfg: &Config,
    queries: &[Query],
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
        // Cloned off before `run_one` consumes `fetched`, and reattached below whatever the
        // stages do with it. Provenance has to survive a failure: a run that errored on an
        // accession still has to be able to say which bytes it errored on.
        let payload = fetched.payload.clone();

        // A stage returning `Err` is already handled below. A stage that *panics* (e.g. an
        // index out of bounds on a malformed graph) must not take the worker down with it: the
        // accession's buffer has already been consumed from the queue, so an uncaught panic
        // would send no result. That silently drops the accession, and if a systematic bug
        // kills every compute worker, the fetch pool then blocks forever on a full bounded
        // channel and the whole run hangs. Catching here turns a panic into one `Failed`
        // result and keeps the worker alive.
        let panic_acc = accession.clone();
        let mut result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            run_one(cfg, queries, query_index, fetched)
                .unwrap_or_else(|e| result_for_error(index, accession, e.to_string(), &e))
        }))
        .unwrap_or_else(|p| {
            // WP-18: say what it was. The payload is the only description of a panic that
            // reaches the summary, and it was being dropped.
            AccessionResult::failed(
                index,
                panic_acc,
                format!("internal error: a stage panicked: {}", panic_message(&*p)),
            )
        });
        result.payload = payload;
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
/// worker is one `UnitigGraph`: ~43 MB of arrays for `ERR2028820`, up to 1.7x the
/// compressed size at peak in general (`Config::GRAPH_PEAK_PER_COMPRESSED_BYTE`).
///
/// Peak is therefore roughly `compressed + graph`: measured at 80 MB for the whole process
/// on `ERR2028820` (33 MB), and 853 MB on `ERR9831160` (348 MB), against ~199 MB for
/// `ERR2028820` in a design that materialises the decompressed buffer.
///
/// For the whole-process ceiling, multiplying this by `cfg.compute_workers` is **not**
/// enough: compressed buffers are also resident in the channel and in blocked fetch workers.
/// See the module doc and `Config::peak_compressed_buffers`.
fn run_one(
    cfg: &Config,
    queries: &[Query],
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

    // Pass 1 of 2. Adjacency only, with anchoring fused into the decode: `Scanner` sees each
    // base as it streams past, so neither the sequences nor a second walk over them exists.
    // `fetched` stays alive until pass 2, which is why it is not dropped here.
    let t = Instant::now();
    let mut scanner = stages::anchor::Scanner::new(cfg, query_index)?;
    let mut graph = stages::parse::scan_in(
        &fetched.accession,
        stages::stream::decoder(&fetched)?,
        // The budget the startup line promised: one graph per compute worker at
        // GRAPH_PEAK_PER_COMPRESSED_BYTE times the size cap. Without it that line is a
        // calibration from real accessions rather than a limit anything enforces.
        cfg.max_graph_bytes_per_accession(),
        &mut scanner,
        cfg.spool_dir.as_deref(),
        crate::types::SPILL_ADJACENCY_AT,
    )?;
    timing.stream = fetch_secs;
    // `parse` now covers the anchoring too, and `anchor` is no longer a pass of its own. The
    // key stays in the contract rather than being dropped, because removing it would make two
    // schema-1 summaries disagree about which keys exist.
    timing.parse = t.elapsed().as_secs_f64();
    timing.anchor = 0.0;
    let n_unitigs = graph.n_unitigs;
    // WG-02 and WG-07. Read off the graph now: `collect` rewrites the graph down to the
    // neighbourhood, and both of these are statements about the file.
    let abundance_missing = graph.abundance_missing;
    let graph_has_links = graph.n_links > 0;
    let anchors = scanner.finish();

    // No anchors means the gene is not present in this accession. A legitimate outcome,
    // not a failure, and there is nothing to expand from.
    if anchors.ids.is_empty() {
        timing.total = t_start.elapsed().as_secs_f64();
        return Ok(AccessionResult {
            index,
            accession,
            // Overwritten by `compute_worker`, which is the last holder of the fetched
            // bytes; `run_one` never sees them as anything but a stream.
            payload: crate::types::Payload::default(),
            outcome: Outcome::NoAnchors,
            n_unitigs,
            abundance_missing,
            n_anchors: 0,
            subgraph_nodes: 0,
            haplotype: None,
            haplotypes: Vec::new(),
            alts: Vec::new(),
            bubbles: Vec::new(),
            timing,
            error: None,
            anchor_used: None,
            candidates: Vec::new(),
            subgraph_fraction: None,
            snv_sites: None,
            ka_max: 0.0,
            depth_saturated: false,
            provisional: false,
            prune: None,
            majority_corrections: 0,
            majority: crate::types::MajorityReport::default(),
            variants: crate::types::VariantStats::default(),
        });
    }

    let t = Instant::now();
    let nodes = stages::extract::expand(cfg, &graph, &anchors)?;
    let subgraph_fraction = stages::extract::subgraph_fraction(&graph, &nodes);
    // The memory term the two-pass parser cannot remove. Everything from here on is
    // proportional to the neighbourhood, so this is where an accession whose gene turns out
    // to be most of its graph has to be refused rather than attempted. Checked before
    // `collect`, which is the first thing that would allocate against it.
    if cfg.max_subgraph_nodes > 0 && nodes.len() > cfg.max_subgraph_nodes {
        timing.extract = t.elapsed().as_secs_f64();
        timing.total = t_start.elapsed().as_secs_f64();
        let mut r = AccessionResult::empty(
            index,
            accession,
            Outcome::SubgraphTooLarge,
            Some(format!(
                "neighbourhood is {} unitigs of {n_unitigs} ({:.1}% of the graph), over \
                 --max-subgraph-nodes {}",
                nodes.len(),
                subgraph_fraction * 100.0,
                cfg.max_subgraph_nodes
            )),
        );
        r.n_unitigs = n_unitigs;
        r.n_anchors = anchors.ids.len();
        r.subgraph_nodes = nodes.len();
        r.subgraph_fraction = Some(subgraph_fraction);
        r.timing = timing;
        return Ok(r);
    }
    // Pass 2 of 2: sequence and abundance for the neighbourhood alone. This is where the
    // graph stops being O(accession) and becomes O(subgraph), and it is the last thing the
    // compressed bytes are needed for.
    graph = stages::parse::collect(
        graph,
        stages::stream::decoder(&fetched)?,
        &nodes,
        cfg.max_graph_bytes_per_accession(),
    )?;
    drop(fetched);
    // In multi-anchor mode the pool's references are the same gene at slightly different
    // lengths, and the guard runs before one of them has won, so the longest stands for the
    // gene. It is also the most forgiving choice: the largest denominator flags the fewest.
    let ref_len = queries.iter().map(|q| q.seq.len()).max().unwrap_or(0);
    let provisional = stages::extract::exceeds_guard(cfg, &nodes, ref_len);
    if provisional {
        eprintln!(
            "[{}] WARNING: the 5-hop neighbourhood holds {:.1} unitigs per reference base \
             ({} unitigs over {} bases, {:.2}% of the graph), above the {:.1} guard. The \
             haplotype is likely under-determined; treat this accession's result as \
             provisional.",
            graph.accession,
            stages::extract::nodes_per_ref_base(&nodes, ref_len),
            nodes.len(),
            ref_len,
            subgraph_fraction * 100.0,
            cfg.subgraph_nodes_per_base_guard
        );
    }
    let seqs = stages::extract::collect_seqs(&graph, &nodes)?;
    let mut subgraph = stages::extract::write_gfa(cfg, &graph, nodes, seqs)?;
    timing.extract = t.elapsed().as_secs_f64();
    let subgraph_nodes = subgraph.nodes.len();

    // align now returns one best GafPath per reference in the anchor pool. Multi-anchor
    // mode fans out over the pool: stitch each reference against its own length and let
    // `pick_primary` choose the winning haplotype. Under the single-anchor path (the
    // legacy shape used by every golden fixture) this reduces to one candidate and picks
    // it iff it clears `min_cov`, which is exactly what the old code did.
    //
    // An empty map == no alignment == Outcome::NoPath (not an error). An empty GFA making
    // GraphAligner abort IS an error, and must not be silently folded into this branch.
    let gaf_path = subgraph.gfa_path.with_extension("gaf");
    let t = Instant::now();
    let mut gaf_by_ref = stages::align::run(cfg, &subgraph, queries, &gaf_path)?;
    timing.align = t.elapsed().as_secs_f64();

    let t = Instant::now();
    let candidates: Vec<stages::stitch::Candidate> = queries
        .iter()
        .enumerate()
        .map(|(pool_index, q)| {
            let path = gaf_by_ref.remove(q.record_id.as_str());
            let hap = match &path {
                Some(p) => {
                    // WA-12. The first point a row from the aligner is believed. A row that
                    // is not a walk of the GFA we handed out stitches into a sequence anyway,
                    // and carries the identity the row claims while describing a different
                    // molecule, so it is refused here rather than published.
                    stages::stitch::check_gaf_walk(cfg, &graph, &subgraph, p).map_err(
                        |detail| crate::error::Error::GraphAligner {
                            accession: graph.accession.clone(),
                            detail: format!("GAF row for {:?}: {detail}", q.record_id),
                        },
                    )?;
                    stages::stitch::stitch(cfg, &graph, p, &subgraph, q)?
                }
                None => None,
            };
            Ok::<_, crate::error::Error>(stages::stitch::Candidate {
                record_id: q.record_id.clone(),
                pool_index,
                ref_len: q.seq.len(),
                path,
                haplotype: hap,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    // Snapshot per-reference stats before `pick_primary` consumes the candidates. Every
    // number here was computed by `stitch`/`align` for this accession; without the snapshot
    // the non-winning references' cov/pident are dropped on the floor. Only populated in
    // multi-anchor mode so single-anchor `summary.json` stays byte-identical.
    let multi_anchor = queries.len() > 1;
    let candidate_reports: Vec<crate::types::CandidateReport> = if multi_anchor {
        candidates
            .iter()
            .map(|c| crate::types::CandidateReport {
                record_id: c.record_id.clone(),
                pool_index: c.pool_index,
                pident: c.haplotype.as_ref().map(|h| h.pident).or_else(|| {
                    c.path.as_ref().map(|p| p.pident())
                }),
                cov: c.haplotype.as_ref().map(|h| h.cov),
                // From the path, not from the haplotype, so a candidate the `--min-cov` gate
                // rejected still reports the number that rejected it (WA-05). `stitch`
                // returns `None` for those, so `haplotype` is `None` and `cov` is null.
                ref_cov: c
                    .path
                    .as_ref()
                    .map(|p| p.query_span() as f64 / c.ref_len.max(1) as f64),
                ref_span: c.path.as_ref().map(|p| (p.query_start, p.query_end)),
                length: c.haplotype.as_ref().map(|h| h.seq.len()),
                nodes: c.haplotype.as_ref().map(|h| h.nodes),
                path_len: c.path.as_ref().map(|p| p.path_len),
                num_matches: c.path.as_ref().map(|p| p.num_matches),
                aln_len: c.path.as_ref().map(|p| p.aln_len),
                selected: false, // set below once pick_primary names a winner
            })
            .collect()
    } else {
        Vec::new()
    };
    let winner = stages::stitch::pick_primary(candidates, cfg.primary_metric, cfg.primary_tie_eps);
    // Alternative-haplotype recovery is `stages::variants`: it corrects the consensus to the
    // majority allele and emits second alleles at biallelic sites. See that module for why it
    // replaced the upstream `detour` engine (gap rate and VAF labelling).
    let (mut final_path, mut haplotype, winning_ref): (
        Option<GafPath>,
        Option<crate::types::Haplotype>,
        Option<String>,
    ) = match winner {
        Some((name, path, hap)) => (Some(path), Some(hap), Some(name)),
        None => (None, None, None),
    };
    // The winning reference's Query is what every subsequent stage aligns against — prune
    // pass, prefer_majority, find_alts. `queries[0]` is a defensible fallback when nothing
    // won: the downstream code paths that follow all bail on `final_path == None` before
    // touching `query`.
    let query: &Query = winning_ref
        .as_deref()
        .and_then(|rid| queries.iter().find(|q| q.record_id == rid))
        .unwrap_or(&queries[0]);
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
    // The floor is a fraction of the sample's *own* depth, never below `--prune-min-abund`,
    // so on a shallow isolate the lower bound is the whole rule: `ERR2028820` at 180x gets
    // 1.8x from the fraction and prunes at 5.0. That is not "nothing is removed", which this
    // comment used to claim (WG-05): it removes 268 of the neighbourhood's 747 unitigs. What
    // is measured, and what makes it safe, is that the *haplotype* does not move. At every
    // setting from 0.01 to 0.05 `ERR2028820` comes back byte-identical, with the same
    // alignment count and no node dropped from the walked path; the 268 are unitigs hanging
    // off the route, which is what the prune is for.
    //
    // # Never a regression
    //
    // If the second pass finds no haplotype, or one spanning less of the reference (see
    // `pruned_covers_as_much`), the first pass's result stands. Pruning can improve an
    // accession and cannot break one. The cost is one extra alignment.
    let mut prune: Option<PruneReport> = None;
    if cfg.prune_frac > 0.0 {
        if let (Some(path), Some(_)) = (&final_path, &haplotype) {
            let depth = stages::extract::path_depth(cfg, &graph, &path.steps);
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
            // Route unitigs the floor would remove although the route had no deeper way
            // round them: a genuinely thin stretch of the gene rather than an error branch
            // (WP-01). Computed before the count below, because they are part of the
            // predicate it must agree with.
            //
            // A failure to read adjacency here is not worth the accession: it means the prune
            // cannot be sized, and an unpruned result is the "never a regression" answer.
            let protected = match stages::extract::protected_route_nodes(&graph, &path.steps, floor) {
                Ok(p) => p,
                Err(e) => {
                    eprintln!(
                        "[{}] prune: could not check the route against the floor ({e}); \
                         keeping the unpruned result",
                        graph.accession
                    );
                    Default::default()
                }
            };

            // The predicate is `extract::kept_by_prune`, the same function `prune_nodes`
            // calls, so this count cannot drift from the pruning it reports on (TW-06).
            let surviving = subgraph
                .nodes
                .iter()
                .filter(|&&id| stages::extract::kept_by_prune(&graph, id, floor, &protected))
                .count();

            let unpruned_span = (path.query_start, path.query_end);
            let report = |result, pruned_span| PruneReport {
                result,
                depth,
                floor,
                unitigs: before,
                unitigs_kept: surviving,
                unpruned_span,
                pruned_span,
            };

            match prune_outcome(surviving, before) {
                PruneOutcome::WouldEmpty => {
                prune = Some(report(PruneResult::WouldEmpty, None));
                // Pruning to nothing is not a prune. Keep the unpruned subgraph, which is the
                // same "never a regression" rule the failed-realign branch below applies.
                eprintln!(
                    "[{}] prune: depth {:.0}x, floor {:.1}x would remove all {} unitigs; \
                     skipping the prune and keeping the unpruned result",
                    graph.accession, depth, floor, before
                );
                }
                PruneOutcome::NothingToRemove => {
                prune = Some(report(PruneResult::NothingToRemove, None));
                eprintln!(
                    "[{}] prune: depth {:.0}x, floor {:.1}x, nothing to remove",
                    graph.accession, depth, floor
                );
                }
                PruneOutcome::Prune => {
                // The first subgraph is consumed here: its node and sequence vectors move
                // into the prune, and its `Drop` deletes its GFA. The refusal paths below
                // rebuild it from `graph`, which still holds every decoded sequence.
                let nodes = std::mem::take(&mut subgraph.nodes);
                let seqs = std::mem::take(&mut subgraph.seqs);
                let (kept_nodes, kept_seqs) =
                    stages::extract::prune_nodes(&graph, nodes, seqs, floor, &protected);
                let dropped = before - kept_nodes.len();
                drop(subgraph);
                let t = Instant::now();
                // WP-04. Every step of the second pass used to be `?`, so GraphAligner
                // aborting or being OOM-killed on the pruned graph, or a full tmpfs during
                // the second GFA write, turned an accession that was already reconstructed
                // into `ERROR`: missing from the FASTA, `prune: null` in the summary, and
                // exit 3 for a one-accession run. That breaks this stage's own contract, that
                // pruning can improve an accession and cannot break one. The second pass is
                // now fallible as a whole, and failing it is a refusal like any other.
                type PassTwo = (
                    crate::types::Subgraph,
                    Option<(GafPath, crate::types::Haplotype)>,
                );
                let attempt = (|| -> Result<PassTwo> {
                    let pruned = stages::extract::write_gfa(cfg, &graph, kept_nodes, kept_seqs)?;
                    let gaf2 = pruned.gfa_path.with_extension("pass2.gaf");
                    let mut gaf2_map =
                        stages::align::run(cfg, &pruned, std::slice::from_ref(query), &gaf2)?;
                    let second = match gaf2_map.remove(query.record_id.as_str()) {
                        Some(p) => {
                            // WA-12, on the pruned graph this time. An inconsistent second row
                            // is an `Err` out of this closure, which WP-04 turns into a
                            // refused prune: the first pass's haplotype stands and the
                            // accession is not lost.
                            stages::stitch::check_gaf_walk(cfg, &graph, &pruned, &p).map_err(
                                |detail| crate::error::Error::GraphAligner {
                                    accession: graph.accession.clone(),
                                    detail: format!("pass 2 GAF row: {detail}"),
                                },
                            )?;
                            stages::stitch::stitch(cfg, &graph, &p, &pruned, query)?
                                .map(|h| (p, h))
                        }
                        None => None,
                    };
                    Ok((pruned, second))
                })();
                timing.align += t.elapsed().as_secs_f64();
                if let Err(e) = &attempt {
                    eprintln!(
                        "[{}] prune: realignment failed ({e}); keeping the unpruned result",
                        graph.accession
                    );
                    prune = Some(report(PruneResult::RefusedError, None));
                }
                let (pruned, second) = match attempt {
                    Ok(v) => (Some(v.0), v.1),
                    Err(_) => (None, None),
                };
                let failed = pruned.is_none();
                let pruned_span = second.as_ref().map(|(p, _)| (p.query_start, p.query_end));
                match second {
                    Some((p, h))
                        if pruned_covers_as_much(path, &p)
                            && !trimmed_copy(
                                &h.seq,
                                haplotype.as_ref().map(|u| u.seq.as_slice()).unwrap_or(&[]),
                            ) =>
                    {
                        prune = Some(report(PruneResult::Applied, pruned_span));
                        eprintln!(
                            "[{}] prune: depth {:.0}x, floor {:.1}x, {} of {} unitigs kept; \
                             realigned, min coverage on the route now {:.1}x",
                            graph.accession, depth, floor, before - dropped, before,
                            stages::extract::path_depth(cfg, &graph, &p.steps).min(
                                p.steps.iter()
                                    .map(|s| graph.abundance_of(s.node))
                                    .fold(f32::INFINITY, f32::min)
                            )
                        );
                        haplotype = Some(h);
                        final_path = Some(p);
                        subgraph = pruned.expect("an applied prune has a pruned subgraph");
                    }
                    rejected => {
                        // `failed` has already logged and recorded `RefusedError`; the rest
                        // of this arm is the rebuild, which every refusal needs.
                        if !failed {
                            let trimmed = rejected.as_ref().is_some_and(|(p, h)| {
                                pruned_covers_as_much(path, p)
                                    && trimmed_copy(
                                        &h.seq,
                                        haplotype.as_ref().map(|u| u.seq.as_slice()).unwrap_or(&[]),
                                    )
                            });
                            let result = match (&rejected, trimmed) {
                                (Some(_), true) => PruneResult::RefusedTrimmed,
                                (Some(_), false) => PruneResult::RefusedShorterSpan,
                                (None, _) => PruneResult::RefusedNoHaplotype,
                            };
                            prune = Some(report(result, pruned_span));
                            match (&rejected, trimmed) {
                                (Some((_, h)), true) => eprintln!(
                                    "[{}] prune: depth {:.0}x, floor {:.1}x; the realigned \
                                     haplotype is the unpruned one {} base(s) shorter and \
                                     otherwise identical, so the prune changed nothing worth \
                                     the loss; keeping the unpruned result",
                                    graph.accession, depth, floor,
                                    haplotype.as_ref().map_or(0, |u| u.seq.len())
                                        .saturating_sub(h.seq.len())
                                ),
                                (Some((p, _)), false) => eprintln!(
                                    "[{}] prune: depth {:.0}x, floor {:.1}x; the realigned \
                                     route spans reference bases {}-{} of {}, {} short of the \
                                     unpruned {}-{} and past the {}-base tolerance; keeping \
                                     the unpruned result",
                                    graph.accession, depth, floor, p.query_start, p.query_end,
                                    query.seq.len(),
                                    path.query_span().saturating_sub(p.query_span()),
                                    path.query_start, path.query_end, PRUNE_SPAN_TOLERANCE
                                ),
                                (None, _) => eprintln!(
                                    "[{}] prune: depth {:.0}x, floor {:.1}x cost the \
                                     alignment; keeping the unpruned result",
                                    graph.accession, depth, floor
                                ),
                            }
                        }
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
    if !crate::scratch::keep_scratch() {
        let _ = std::fs::remove_file(&subgraph.gfa_path);
    }
    // Intra-host variant recovery. `prune_frac` is not part of this: it cleans the graph
    // before this runs. `prefer_majority` and the alt/vcf/snv outputs each gate themselves, so
    // with all of them off this reduces to emitting the single stitched primary.
    let t = Instant::now();
    let mut alts: Vec<crate::types::AltHaplotype> = Vec::new();
    // Correct the consensus first, so the alternatives that follow are the minority alleles
    // relative to a corrected route rather than to a reference-biased one.
    // Largest abundance anywhere in the neighbourhood, for the saturation flag.
    let ka_max = subgraph
        .nodes
        .iter()
        .map(|&id| graph.abundance_of(id))
        .fold(0.0f32, f32::max);
    // WV-05. The allele floor is a flat `--min-abund` (5x) unless `--min-abund-frac` is set,
    // which it is not by default. That is safe only because the prune is supposed to have
    // removed sub-floor unitigs first. When the prune is REFUSED, the unpruned subgraph is
    // what the variant stages see, and the 5x gate lets through exactly the branches the
    // prune would have deleted: a gene at 1000x with a fork at 30x against 8x (below the 10x
    // prune floor) published `alt1 vaf=0.211`, a PASS VCF row and an snvs.tsv row, with only
    // `prune.outcome` in the summary to hint at it.
    //
    // So where the prune wanted a floor and did not get to apply it, the variant stages
    // inherit that floor. Applied and pruned-off runs are untouched.
    let raised_cfg = inherited_floor_config(cfg, prune.as_ref());
    let vcfg: &Config = raised_cfg.as_ref().unwrap_or(cfg);
    let mut majority_corrections = 0usize;
    let mut majority = crate::types::MajorityReport::default();
    let mut variant_stats = crate::types::VariantStats::default();
    if let (Some(path), Some(_)) = (&final_path, &haplotype) {
        let (fixed, report) = stages::variants::prefer_majority(vcfg, &graph, &subgraph, path, query)?;
        // The sites left on the minority allele are reported whether or not anything was
        // swapped: a route with no corrections and four forks it declined to move is not the
        // same result as a route with nothing to correct (WV-02).
        majority = report;
        if let Some(fixed) = fixed {
            if let Some(mut h) = stages::stitch::stitch(vcfg, &graph, &fixed, &subgraph, query)? {
                // The route was spliced onto the majority allele, so the GAF's identity
                // describes the sequence we just replaced, not this one.
                h.pident_inherited = true;
                haplotype = Some(h);
                final_path = Some(fixed);
                // Only count corrections that reached the output. A swap whose restitch
                // failed left the published primary unchanged, so reporting it would
                // overstate what happened to this record.
                majority_corrections = report.swaps;
            }
        }
    }
    if let (Some(path), Some(_)) = (&final_path, &haplotype) {
        let (found, stats) = stages::variants::find_alts(vcfg, &graph, &subgraph, path, query)?;
        alts = found;
        variant_stats = stats;
        // WA-04. The mixture numbers belong to the sequence that is published, so they are
        // read off the route AFTER `prefer_majority` has decided what that route is, and
        // copied onto the haplotype the FASTA writer will see. They describe the primary;
        // an alt is one allele changed at one site and carries its own frequencies already.
        if let Some(h) = haplotype.as_mut() {
            h.two_allele_sites = stats.two_allele_sites;
            h.min_site_vaf = stats.min_site_vaf;
        }
    }
    // Bubbles are computed here, after the route is settled, for two reasons (WG-03): the
    // floor they are gated at comes from the gene's depth, and there is no depth without a
    // route; and `subgraph` is the pruned graph by now when the prune applied. The stage
    // itself still never looks at the GAF, it reports the neighbourhood's diamonds, which is
    // why a `no_path` accession still reports bubbles, at the absolute floor, since nothing
    // measured its depth.
    //
    // The floor is the strictest of the three the run already uses: `--prune-frac` of depth
    // (a branch the prune would delete is not an allele), `--min-abund-frac` of depth (the
    // variant stages' own relative floor), and the absolute `--min-abund`.
    //
    // Two things decide what is left to report, and it is worth being clear about which does
    // the work. The ORDER does most of it: an applied prune replaces `subgraph` with the
    // pruned one above, so by here the sub-floor branches are gone and the diamonds built
    // from them cannot be found. Running before the alignment, as this stage used to, meant
    // reporting the neighbourhood as extracted, including everything the prune then deleted.
    //
    // The FLOOR covers the cases where no prune cleaned the graph: `--prune-frac 0`, where it
    // reduces to the variant stages' own floor, and a refused prune, where WV-05 raises their
    // floor to this same number because the sub-floor branches are still there. Taking the
    // prune's fraction is that same argument: a branch the prune judged to be error is not an
    // allele, whether or not the prune got to act on it.
    let t = Instant::now();
    let bubble_floor = final_path.as_ref().map_or(cfg.min_abund, |p| {
        let depth = stages::extract::path_depth(cfg, &graph, &p.steps);
        stages::extract::relative_floor(cfg.prune_frac.max(cfg.min_abund_frac), depth, cfg.min_abund)
    });
    let bubbles = stages::phase::phase_with_floor(cfg, &graph, &subgraph, bubble_floor)?;
    timing.phase = t.elapsed().as_secs_f64();

    // Tag the primary with the winning reference when the anchor pool held more than one.
    // Single-anchor mode leaves anchor_used = None on Haplotype, alts, and AccessionResult
    // so the FASTA header, alt header, snvs.tsv and summary.json stay byte-identical to
    // today's shape. Multi-anchor mode tags them all with the winning reference id.
    let anchor_used_for_output: Option<String> = if multi_anchor { winning_ref.clone() } else { None };
    // Mark the winning candidate in the snapshot so the summary writer can flag it.
    let mut candidate_reports = candidate_reports;
    if let Some(rid) = winning_ref.as_deref() {
        for c in candidate_reports.iter_mut() {
            if c.record_id == rid {
                c.selected = true;
            }
        }
    }
    let mut haplotypes: Vec<crate::types::Haplotype> = haplotype.clone().into_iter().collect();
    if multi_anchor {
        for h in haplotypes.iter_mut() {
            h.anchor_used = winning_ref.clone();
        }
        for a in alts.iter_mut() {
            a.haplotype.anchor_used = winning_ref.clone();
        }
    }
    timing.stitch += t.elapsed().as_secs_f64();
    timing.total = t_start.elapsed().as_secs_f64();

    // A path that failed the coverage gate is no_path.
    let haplotype = haplotypes.first().cloned();
    let outcome = if !haplotypes.is_empty() {
        Outcome::Success
    } else if !graph_has_links {
        // WG-07. The file carried no `L:` token anywhere, so there was no adjacency to walk
        // and nothing could have spanned more than one unitig. Reporting that as `no_path`
        // says "this accession does not carry the gene", which is a claim about the sample
        // this run is in no position to make. Logan publishes such objects: 0 links in 15,864
        // headers for `ERR3988483`.
        //
        // Checked only when nothing was reconstructed. A gene that fits inside a single
        // unitig needs no edge and is a `success` like any other.
        eprintln!(
            "[{accession}] no_links: the unitig file carries no L: adjacency ({n_unitigs} \
             records, 0 links), so no path spanning more than one unitig exists. This is a \
             property of the Logan object, not evidence that the gene is absent."
        );
        Outcome::NoLinks
    } else {
        Outcome::NoPath
    };

    Ok(AccessionResult {
        index,
        accession,
        payload: crate::types::Payload::default(),
        outcome,
        alts,
        n_unitigs,
        abundance_missing,
        n_anchors: anchors.ids.len(),
        subgraph_nodes,
        haplotype,
        haplotypes,
        bubbles,
        timing,
        error: None,
        anchor_used: anchor_used_for_output,
        candidates: candidate_reports,
        subgraph_fraction: Some(subgraph_fraction),
        // WV-03. Measured over the whole extracted neighbourhood rather than the walked route,
        // because the prune floor and every fork total are drawn from it, not just the path.
        snv_sites: None,
        ka_max,
        depth_saturated: ka_max >= crate::types::KA_CEILING,
        provisional,
        prune,
        majority_corrections,
        majority,
        variants: variant_stats,
    })
}

#[cfg(test)]
mod floor_tests {
    use super::inherited_floor_config;
    use crate::config::Config;
    use crate::types::{PruneReport, PruneResult};

    fn rep(result: PruneResult, floor: f32) -> PruneReport {
        PruneReport {
            result,
            depth: 1000.0,
            floor,
            unitigs: 100,
            unitigs_kept: 90,
            unpruned_span: (0, 100),
            pruned_span: None,
        }
    }

    /// WV-05. The allele floor is a flat 5x by default, which is only safe because the prune
    /// is supposed to have removed sub-floor unitigs first. When the prune is refused, the
    /// variant stages see the unpruned subgraph, and the 5x gate admits exactly the branches
    /// the prune meant to delete: a gene at 1000x with a fork at 30x against 8x published
    /// `vaf=0.211`, a PASS VCF row and an snvs.tsv row, with only `prune.outcome` to hint at
    /// it. Where the prune wanted a floor and could not apply it, the variant stages inherit
    /// it.
    #[test]
    fn a_refused_prune_hands_its_floor_to_the_variant_stages() {
        let cfg = Config { min_abund: 5.0, ..Config::default() };
        for refused in [
            PruneResult::RefusedShorterSpan,
            PruneResult::RefusedNoHaplotype,
            PruneResult::WouldEmpty,
        ] {
            let raised = inherited_floor_config(&cfg, Some(&rep(refused, 10.0)))
                .unwrap_or_else(|| panic!("{refused:?} must raise the floor"));
            assert_eq!(raised.min_abund, 10.0);
            // Only the floor moves; nothing else about the run changes.
            assert_eq!(raised.min_vaf, cfg.min_vaf);
            assert_eq!(raised.max_alts, cfg.max_alts);
        }
    }

    /// A prune that ran, or had nothing to do, already removed what it wanted to. Those runs
    /// must stay byte-identical, so the config is left alone.
    #[test]
    fn an_applied_or_unnecessary_prune_changes_nothing() {
        let cfg = Config { min_abund: 5.0, ..Config::default() };
        for ok in [PruneResult::Applied, PruneResult::NothingToRemove] {
            assert!(inherited_floor_config(&cfg, Some(&rep(ok, 10.0))).is_none(), "{ok:?}");
        }
        // No prune at all (`--prune-frac 0`), and a floor that is already the minimum.
        assert!(inherited_floor_config(&cfg, None).is_none());
        assert!(
            inherited_floor_config(&cfg, Some(&rep(PruneResult::WouldEmpty, 5.0))).is_none(),
            "a floor equal to --min-abund is not a raise"
        );
    }
}

#[cfg(test)]
mod ordering_tests {
    use crate::types::{AccessionResult, Outcome};

    /// Results arrive in completion order, decided by S3 and by how long each graph takes.
    /// `run` sorts on `index` before returning so every output file is byte-reproducible and
    /// `--sieve`'s representative is the first record of the *input*, not of a race.
    #[test]
    fn sorting_on_index_restores_input_order() {
        // The measured completion order from a nine-accession Zika run.
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

        super::restore_input_order(&mut results);

        let got: Vec<&str> = results.iter().map(|r| r.accession.as_str()).collect();
        assert_eq!(got, input, "stored order must be the input order");
    }

    /// No Logan data and the size cap are negatives, so `--allow-empty` can accept them; any
    /// other fetch failure is an error it never hides.
    #[test]
    fn a_fetch_failure_is_classified_by_its_cause() {
        use crate::error::Error;
        let acc = || "SRR1".to_string();
        let cases = [
            (Error::NoLoganData { accession: acc() }, Outcome::NoLoganData),
            (Error::TooLarge { accession: acc(), detail: "big".into(),
              payload: Box::new(crate::types::Payload::default()) }, Outcome::TooLarge),
            (Error::Fetch { accession: acc(), detail: "HTTP 503".into() }, Outcome::Failed),
        ];
        for (e, want) in cases {
            let r = super::result_for_error(7, acc(), e.to_string(), &e);
            assert_eq!((r.index, r.outcome), (7, want), "{e}");
        }
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
    use super::{prune_outcome, pruned_covers_as_much, trimmed_copy, PruneOutcome};
    use crate::config::Config;
    use crate::stages::extract::{depth_floor, prune_nodes};
    use crate::stages::parse;
    use crate::types::UnitigGraph;

    fn graph(fa: &str) -> UnitigGraph {
        parse::build(&"T".into(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX).unwrap()
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

    fn spanning(query_start: usize, query_end: usize, path_len: usize) -> crate::types::GafPath {
        crate::types::GafPath {
            steps: Vec::new(),
            query_start,
            query_end,
            path_len,
            path_start: 0,
            path_end: path_len,
            num_matches: 0,
            aln_len: 0,
        }
    }

    /// `ERR6343483`'s numbers: the pruned route lost the first 424 reference bases, and the
    /// haplotype was still long enough to pass `--min-cov`.
    #[test]
    fn a_pruned_route_that_stops_short_is_refused() {
        let unpruned = spanning(0, 3821, 3816);
        assert!(!pruned_covers_as_much(&unpruned, &spanning(424, 3821, 3391)));
    }

    /// Up to `PRUNE_SPAN_TOLERANCE` reference bases at an end may go; the next one may not.
    /// A strict rule refused a prune over a single terminal base on real cohorts.
    #[test]
    fn a_pruned_route_may_give_up_a_terminal_codon_but_no_more() {
        let unpruned = spanning(0, 3821, 3816);
        for lost in 1..=super::PRUNE_SPAN_TOLERANCE {
            assert!(pruned_covers_as_much(&unpruned, &spanning(0, 3821 - lost, 3816 - lost)),
                    "{lost} base(s) short at the 3' end should be kept");
            assert!(pruned_covers_as_much(&unpruned, &spanning(lost, 3821, 3816 - lost)),
                    "{lost} base(s) short at the 5' end should be kept");
        }
        let over = super::PRUNE_SPAN_TOLERANCE + 1;
        assert!(!pruned_covers_as_much(&unpruned, &spanning(0, 3821 - over, 3816 - over)),
                "one base past the tolerance must be refused");
    }

    /// An indel changes the haplotype's length but not the reference span, so pruning away a
    /// spurious insertion, or onto a real deletion, is kept.
    #[test]
    fn an_indel_is_not_lost_coverage() {
        let unpruned = spanning(0, 3821, 3822);
        assert!(pruned_covers_as_much(&unpruned, &spanning(0, 3821, 3816)), "6-nt deletion");
        assert!(pruned_covers_as_much(&unpruned, &spanning(0, 3821, 3823)), "insertion");
        assert!(pruned_covers_as_much(&spanning(4, 3820, 3816), &spanning(0, 3821, 3821)), "longer");
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
            .filter(|&&id| g.abundance_of(id) >= floor)
            .count();
        assert_eq!(surviving, 0);
        assert_eq!(prune_outcome(surviving, nodes.len()), PruneOutcome::WouldEmpty);

        // And the predicate above matches what prune_nodes would actually do.
        let seqs = vec![b"ACGT".to_vec(); 3];
        let (kept, _) = prune_nodes(&g, nodes, seqs, floor, &Default::default());
        assert!(kept.is_empty(), "prune_nodes agrees the subgraph would be emptied");
    }

    /// WA-09. The three-base tolerance buys a shorter result only when the prune changed the
    /// haplotype; a shorter copy of the same sequence is pure loss.
    #[test]
    fn a_shorter_copy_of_the_same_sequence_is_trimming_and_nothing_else() {
        let unpruned = b"ACGTACGTAA";
        // Trimmed at either end, or both.
        assert!(trimmed_copy(b"ACGTACGTA", unpruned), "one base off the 3' end");
        assert!(trimmed_copy(b"CGTACGTAA", unpruned), "one base off the 5' end");
        assert!(trimmed_copy(b"CGTACGTA", unpruned), "both ends");
        // Same length, so nothing was given up: not this rule's business.
        assert!(!trimmed_copy(unpruned, unpruned));
        // Shorter AND different: the prune did something, which is what the tolerance is for.
        assert!(!trimmed_copy(b"ACGTAGGTA", unpruned), "an interior substitution");
        assert!(!trimmed_copy(b"ACGTCGTAA", unpruned), "an interior deletion");
        // Degenerate inputs must not claim a match.
        assert!(!trimmed_copy(b"", unpruned));
        assert!(!trimmed_copy(b"ACGT", b""));
    }
}

#[cfg(test)]
mod fetch_panic_tests {
    use super::*;

    /// A fetcher that panics for one accession and reports the rest as absent.
    struct Panicky(&'static str);

    impl Fetch for Panicky {
        fn fetch(&self, _index: usize, accession: &Accession) -> Result<FetchedUnitigs> {
            if accession == self.0 {
                panic!("synthetic fetch panic");
            }
            Err(crate::error::Error::NoLoganData {
                accession: accession.clone(),
            })
        }
    }

    /// WP-18. The compute pool has caught its panics since TW-01, the fetch pool did not,
    /// and the difference is the whole run: an unwinding fetch thread re-panics in
    /// `thread::scope`, so `main` exits 101 with no FASTA, no summary and no `skipped.tsv`,
    /// and the collector never runs. At the scale this release is for, that is hours of
    /// fetching and alignment for one panic in the pool that touches the network, zstd and
    /// the filesystem.
    #[test]
    fn a_panicking_fetch_fails_its_own_accession_and_the_others_still_land() {
        let accs = ["SRR1".to_string(), "SRRBOOM".to_string(), "SRR3".to_string()];
        let work: SharedReceiver<(usize, Accession)> = {
            let (tx, rx) = std::sync::mpsc::channel();
            for (i, a) in accs.iter().cloned().enumerate() {
                tx.send((i, a)).unwrap();
            }
            Arc::new(Mutex::new(rx))
        };
        let (fetch_tx, fetch_rx) = sync_channel::<FetchedUnitigs>(4);
        let (result_tx, result_rx) = std::sync::mpsc::channel::<AccessionResult>();

        // The panic is caught inside the worker, so this returns normally.
        fetch_worker(&Panicky("SRRBOOM"), work, fetch_tx, result_tx);
        drop(fetch_rx);

        let mut results: Vec<AccessionResult> = result_rx.into_iter().collect();
        results.sort_by_key(|r| r.index);
        assert_eq!(results.len(), 3, "every accession must come back with a result");

        assert_eq!(results[1].accession, "SRRBOOM");
        assert_eq!(results[1].outcome, Outcome::Failed);
        let detail = results[1].error.as_deref().unwrap_or("");
        assert!(detail.contains("panicked"), "{detail}");
        // The payload is the only description of a panic that reaches the summary, so it is
        // carried rather than dropped.
        assert!(detail.contains("synthetic fetch panic"), "{detail}");

        // And the panic ended one accession, not the pool.
        for i in [0, 2] {
            assert_eq!(results[i].outcome, Outcome::NoLoganData, "{:?}", results[i].accession);
        }
    }
}
