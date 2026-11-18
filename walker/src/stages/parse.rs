//! Step 2: two streaming passes over the unitig file.
//!
//! Builds a compact CSR rather than dicts keyed by every unitig in the accession, which cost
//! about 700 bytes per unitig against ~114 bytes on disk.
//!
//! A record looks like:
//!
//! ```text
//! >ERR2028820_2933 ka:f:5.3   L:-:443893:+  L:+:1239391:-
//! CTGTATACACCCTTCCATTTTTAAAACAACGTAAGAATTAATTAGGAAAG…
//! ```
//!
//! `L:` tokens reference the bare integer id, not the accession-prefixed one. Edges are
//! bidirected and each edge appears in *both* endpoints' headers, which is what the strand
//! bits are for. Dropping them makes bubble detection mistake one edge read twice for a
//! reconvergence.
//!
//! # Two passes, not one
//!
//! `build` below is the single-pass form, kept for tests. The production path is `scan` then
//! `collect`, because almost all of a full CSR is thrown away: `extract::expand` keeps the
//! anchor neighbourhood, 747 of 1,277,362 unitigs for `ERR2028820`.
//!
//! * `scan` reads adjacency and nothing else, feeding each record's bases to a [`RecordSink`]
//!   as they stream past. `stages::anchor::Scanner` is that sink, so anchoring costs no
//!   storage and no second walk over 68 million bases.
//! * `extract::expand` then finds the neighbourhood from adjacency alone.
//! * `collect` decodes a second time and keeps sequence and abundance for those ids only,
//!   stopping once the last one is past.
//!
//! The second decode is cheap, about 0.25 s per 146 MB of text, because pass 1 packs no 2-bit
//! sequence and stores no abundance. What it buys is 21.5 of the 33.5 bytes per unitig a
//! single-pass graph holds.
//!
//! No regex. The cost of parsing this format is the data structure, not the tokenizer.
//!
//! IN   accession, streaming zstd decoder (twice)
//! OUT  UnitigGraph holding full adjacency and the neighbourhood's sequences

use crate::error::{Error, Result};
use crate::types::{Accession, UnitigGraph, UnitigId};
use std::io::{BufReader, Read, BufRead};

/// Target id occupies the low 30 bits of a packed edge; the two strand flags sit above it.
pub const ID_BITS: u32 = 30;
pub const ID_MASK: u32 = (1 << ID_BITS) - 1;
const SRC_FWD: u32 = 1 << 30;
const DST_FWD: u32 = 1 << 31;

/// Pack a target id with its source and target strand bits into one `u32`.
#[inline]
pub fn pack_edge(target: u32, src_fwd: bool, dst_fwd: bool) -> u32 {
    debug_assert!(target <= ID_MASK);
    target | if src_fwd { SRC_FWD } else { 0 } | if dst_fwd { DST_FWD } else { 0 }
}

/// Inverse of `pack_edge`, returning `(target, src_fwd, dst_fwd)`.
#[inline]
pub fn unpack_edge(packed: u32) -> (u32, bool, bool) {
    (
        packed & ID_MASK,
        packed & SRC_FWD != 0,
        packed & DST_FWD != 0,
    )
}

/// Receives every record as `scan` streams past it.
///
/// Monomorphised rather than a trait object: `base` is called once per base, about 68 million
/// times for a mid-sized accession, and a virtual call there would cost more than the work it
/// dispatches. `stages::anchor::Scanner` is the production implementation; `()` is the no-op
/// for callers that only want adjacency.
pub trait RecordSink {
    /// A new record begins. Any per-record state resets here.
    fn begin(&mut self, id: UnitigId);
    /// One base, as its 2-bit code. Already validated as ACGT by the caller.
    fn base(&mut self, code: u8);
    /// The record is complete.
    fn end(&mut self);
}

impl RecordSink for () {
    #[inline]
    fn begin(&mut self, _id: UnitigId) {}
    #[inline]
    fn base(&mut self, _code: u8) {}
    #[inline]
    fn end(&mut self) {}
}

/// 2-bit code for one base, or an error naming it.
///
/// Shared by all three passes so that what counts as a valid base is defined once. A stray N
/// packed as some other nucleotide would propagate silently into the haplotype, which is why
/// this is an error and not a substitution.
#[inline]
fn base_code(b: u8, record: u64) -> std::result::Result<u8, String> {
    match b {
        b'A' => Ok(0),
        b'C' => Ok(1),
        b'G' => Ok(2),
        b'T' => Ok(3),
        _ => Err(format!(
            "record {record}: non-ACGT base {:?}; 2-bit packing would corrupt it",
            b as char
        )),
    }
}

/// A decode that produced no unitig records: absent data, or a file the walker cannot read.
///
/// An **empty** stream is a property of the data, not a failure of the run. Logan publishes a
/// 55-byte zstd frame, valid and holding nothing, for accessions it assembled nothing from,
/// about 0.4% of the bucket. That is `NoLoganData`, the same negative the 404 path returns,
/// because it is the same fact arriving by a different route.
///
/// A stream that had **content** but no `>` record stays an error. That is a malformed file
/// rather than an absent one, and collapsing the two would hide a real defect: the first is
/// Logan having nothing to say, the second is Logan saying something the walker cannot read.
fn no_records_error(accession: &Accession, bytes_read: u64) -> Error {
    if bytes_read == 0 {
        Error::NoLoganData { accession: accession.clone() }
    } else {
        Error::Parse {
            accession: accession.clone(),
            detail: format!("no unitig records found in {bytes_read} bytes of decompressed data"),
        }
    }
}

/// Classify a failed spill write: out of space is a negative, anything else is an error.
///
/// A spool directory that fills up is the run exceeding a budget the operator set, exactly
/// like a `Content-Length` over `--max-size-mb`, so it gets the same `too_large` outcome
/// rather than one that makes the run look broken. The distinction is not cosmetic: a
/// `too_large` accession is a negative that `--allow-empty` tolerates and that leaves the
/// exit code alone, while a `Parse` error counts towards `AllErrored` and exit code 3. A
/// cohort where one 11 GB accession fills the disk should not report that the walker failed.
///
/// Everything else keeps the error. A short write to a healthy filesystem, an EIO, or a
/// handle that went away are not budget decisions and must stay loud.
///
/// Note that `AdjacencySpool::create` failing is handled separately and more gently, in
/// `scan_in`: a spool that cannot be *created* falls back to RAM with a warning, because the
/// accession can still be parsed the old way. Once spilling has begun that door is shut, as
/// the rows already written are the only copy.
fn spill_error(accession: &Accession, what: &str, after: u64, e: std::io::Error) -> Error {
    match e.kind() {
        std::io::ErrorKind::StorageFull | std::io::ErrorKind::QuotaExceeded => Error::TooLarge {
            accession: accession.clone(),
            detail: format!(
                "the spool directory ran out of space writing adjacency {what} after {after} \
                 unitigs ({e}); lower -w, lower --max-size-mb, or point --spool-dir at a \
                 larger filesystem"
            ),
            // Nothing is known here about the object's provenance: `scan_in` sees a reader,
            // not the fetch that produced it. The size that matters is in `detail`, and it is
            // the spool's, not the payload's.
            payload: Box::new(crate::types::Payload::default()),
        },
        _ => Error::Parse {
            accession: accession.clone(),
            detail: format!("spill {what} after {after} unitigs: {e}"),
        },
    }
}

/// Pass 1: adjacency only, with every base handed to `sink` as it goes.
///
/// Adjacency stays in RAM. Use `scan_in` to let it spill to disk past
/// [`crate::types::SPILL_ADJACENCY_AT`], which is what the pipeline does.
///
/// Returns a `UnitigGraph` whose `nodes`, `abundance`, `seq_2bit` and `seq_offset` are empty.
/// That graph is enough for `extract::expand`, which needs `n_unitigs` and `edges` and nothing
/// else; `collect` fills in the rest for the ids that survive.
pub fn scan<S: RecordSink>(
    accession: &Accession,
    reader: Box<dyn Read + '_>,
    budget_bytes: u64,
    sink: &mut S,
) -> Result<UnitigGraph> {
    scan_in(accession, reader, budget_bytes, sink, None, crate::types::SPILL_ADJACENCY_AT)
}

/// As `scan`, with `spool_dir` naming where adjacency may spill once it holds more than
/// `spill_at` bytes. `None` keeps it in RAM however large it gets, which is what `build` and
/// most tests want. `spill_at` is a parameter rather than the constant so a test can prove
/// the threshold is what decides, without a fixture large enough to cross the real one.
pub fn scan_in<S: RecordSink>(
    accession: &Accession,
    reader: Box<dyn Read + '_>,
    budget_bytes: u64,
    sink: &mut S,
    spool_dir: Option<&std::path::Path>,
    spill_at: usize,
) -> Result<UnitigGraph> {
    let bad = |detail: String| Error::Parse {
        accession: accession.clone(),
        detail,
    };

    let mut rd = BufReader::with_capacity(1 << 18, reader);
    let mut line: Vec<u8> = Vec::with_capacity(512);
    let mut offsets: Vec<u32> = Vec::new();
    let mut targets: Vec<u32> = Vec::new();

    let mut expected_id: u64 = 0;
    // Decompressed bytes seen. What tells an absent file from an unreadable one when no
    // records turn up; see `no_records_error`.
    let mut bytes_read: u64 = 0;
    let mut in_record = false;
    let mut iters: u64 = 0;
    // The ceiling warning used to scan the stored abundances. Nothing stores them here, so the
    // maximum is tracked as it streams past and the warning is unchanged.
    let mut max_abundance = 0.0f32;
    // Records carrying neither `ka:f:` nor `km:f:`. Pass 1 is where this is countable,
    // because pass 2 only reads the neighbourhood's records.
    let mut abundance_missing: u64 = 0;
    // The largest link target seen. Compared against the record count at the end,
    // where the count is finally known. Tracked as rows are appended rather than scanned at
    // the end, because `targets` is a flush buffer once adjacency spills.
    let mut max_target: u32 = 0;
    // Adjacency may end up on disk, so `offsets` cannot be indexed by `targets.len()`: the
    // Vecs are flush buffers and get cleared. This is the only running count of edges.
    let mut n_targets: u64 = 0;
    let mut spool: Option<crate::types::AdjacencySpool> = None;
    let mut spill_failed = false;

    loop {
        // Adjacency is the only thing growing, so the budget is measured over it alone. See
        // `build` for why this exists at all.
        iters += 1;
        if iters.is_multiple_of(BUDGET_CHECK_EVERY) && budget_bytes != u64::MAX {
            let used = (offsets.capacity() + targets.capacity()) as u64 * size_of::<u32>() as u64;
            if used > budget_bytes {
                return Err(Error::TooLarge {
                    accession: accession.clone(),
                    detail: format!(
                        "graph adjacency reached {used} bytes after {expected_id} unitigs, over \
                         the {budget_bytes}-byte per-accession graph budget, which is \
                         the memory bound printed at startup"
                    ),
                    payload: Box::new(crate::types::Payload::default()),
                });
            }
        }

        line.clear();
        let n = rd
            .read_until(b'\n', &mut line)
            .map_err(|e| bad(format!("read failed after {expected_id} records: {e}")))?;
        if n == 0 {
            break;
        }
        bytes_read += n as u64;
        while matches!(line.last(), Some(b'\n') | Some(b'\r')) {
            line.pop();
        }
        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            if in_record {
                sink.end();
            }
            // As in `build`: the offset is the edge count on entry to the record, counted
            // globally because `targets` is a buffer that gets flushed and cleared.
            let edges_before = u32::try_from(n_targets)
                .map_err(|_| bad("edge count exceeds u32 addressing".into()))?;
            let before = targets.len();
            let h = header(&line, expected_id, Some(&mut targets), &bad)?;
            for &packed in &targets[before..] {
                max_target = max_target.max(unpack_edge(packed).0);
            }
            n_targets += (targets.len() - before) as u64;
            offsets.push(edges_before);

            // Spill once adjacency outgrows the threshold, then keep the Vecs as small write
            // buffers. A spill that cannot be created is a warning, not a failure: the
            // accession simply stays in RAM, which is what it would have done before.
            let held = (offsets.len() + targets.len()) * size_of::<u32>();
            if spool.is_none() && !spill_failed && held > spill_at {
                match spool_dir.map(crate::types::AdjacencySpool::create) {
                    Some(Ok(s)) => spool = Some(s),
                    Some(Err(e)) => {
                        eprintln!(
                            "[{accession}] warning: could not spill adjacency to disk ({e}); \
                             keeping it in memory"
                        );
                        spill_failed = true;
                    }
                    None => spill_failed = true,
                }
            }
            if let Some(s) = spool.as_mut() {
                if held >= SPOOL_FLUSH_AT {
                    s.push_offsets(&offsets)
                        .map_err(|e| spill_error(accession, "offsets", expected_id, e))?;
                    offsets.clear();
                    s.push_targets(&targets)
                        .map_err(|e| spill_error(accession, "targets", expected_id, e))?;
                    targets.clear();
                }
            }
            if h.abundance > max_abundance {
                max_abundance = h.abundance;
            }
            if !h.has_abundance {
                abundance_missing += 1;
            }
            sink.begin(h.id);
            expected_id += 1;
            in_record = true;
        } else {
            if !in_record {
                return Err(bad("sequence data before the first header".into()));
            }
            for &b in &line {
                sink.base(base_code(b, expected_id - 1).map_err(&bad)?);
            }
        }
    }
    if in_record {
        sink.end();
    }

    if expected_id == 0 {
        return Err(no_records_error(accession, bytes_read));
    }
    offsets.push(
        u32::try_from(n_targets).map_err(|_| bad("edge count exceeds u32 addressing".into()))?,
    );
    if n_targets > 0 && max_target as u64 >= expected_id {
        return Err(truncated_error(accession, max_target, expected_id));
    }
    warn_ka_ceiling(accession, max_abundance);
    warn_missing_abundance(accession, abundance_missing, expected_id);

    if let Some(s) = spool.as_mut() {
        s.push_offsets(&offsets).map_err(|e| spill_error(accession, "offsets", expected_id, e))?;
        s.push_targets(&targets).map_err(|e| spill_error(accession, "targets", expected_id, e))?;
        offsets = Vec::new();
        targets = Vec::new();
    }

    Ok(UnitigGraph {
        accession: accession.clone(),
        n_unitigs: expected_id as usize,
        offsets,
        targets,
        nodes: Vec::new(),
        abundance: Vec::new(),
        seq_2bit: Vec::new(),
        seq_offset: Vec::new(),
        spool,
        abundance_missing,
        // The running count, not `targets.len()`: once adjacency spills, the Vec is a
        // flush buffer that gets cleared, and on a spilled accession it is empty here.
        n_links: n_targets,
    })
}

/// How much adjacency accumulates in the write buffers between flushes once spilled.
const SPOOL_FLUSH_AT: usize = 1 << 20;

/// Pass 2: fill in sequence and abundance for `wanted`, which must be sorted and unique.
///
/// Everything else in the file is parsed past without being stored, and the pass stops once
/// the last wanted id has been read, so a neighbourhood near the start of a large accession
/// costs a fraction of a decode.
pub fn collect(
    mut graph: UnitigGraph,
    reader: Box<dyn Read + '_>,
    wanted: &[UnitigId],
    budget_bytes: u64,
) -> Result<UnitigGraph> {
    let accession = graph.accession.clone();
    let bad = |detail: String| Error::Parse {
        accession: accession.clone(),
        detail,
    };
    debug_assert!(wanted.windows(2).all(|w| w[0] < w[1]), "wanted must be sorted and unique");

    graph.nodes = wanted.to_vec();
    graph.abundance = Vec::with_capacity(wanted.len());
    graph.seq_offset = Vec::with_capacity(wanted.len() + 1);
    graph.seq_2bit.clear();
    if wanted.is_empty() {
        graph.seq_offset.push(0);
        return Ok(graph);
    }
    let last_wanted = *wanted.last().expect("non-empty");

    let mut rd = BufReader::with_capacity(1 << 18, reader);
    let mut line: Vec<u8> = Vec::with_capacity(512);
    let mut expected_id: u64 = 0;
    let mut total_bases: u64 = 0;
    // Index into `wanted`, which advances monotonically because ids are dense and ascending.
    let mut next = 0usize;
    let mut keeping = false;

    loop {
        line.clear();
        let n = rd
            .read_until(b'\n', &mut line)
            .map_err(|e| bad(format!("read failed after {expected_id} records: {e}")))?;
        if n == 0 {
            break;
        }
        while matches!(line.last(), Some(b'\n') | Some(b'\r')) {
            line.pop();
        }
        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            // Links were taken in pass 1 and the adjacency is already built, and a record
            // this pass is skipping needs nothing out of its header but the id. `wanted` is
            // ascending and ids are dense, so `wanted[next] == expected_id` decides it before
            // anything has been parsed; `header_light` then re-establishes `id == expected_id`
            // and the two agree by construction.
            let want = next < wanted.len() && wanted[next] as u64 == expected_id;
            let h = header_light(&line, expected_id, want, &bad)?;
            expected_id += 1;
            if want {
                keeping = true;
                next += 1;
                graph.abundance.push(h.abundance);
                graph.seq_offset.push(
                    u32::try_from(total_bases)
                        .map_err(|_| bad("retained sequence exceeds u32 addressing".into()))?,
                );
                // The same budget as pass 1, now over everything the graph holds. Pass 2 is
                // where the neighbourhood's sequences land, and a neighbourhood can be most of
                // the graph: on amplicon-like data `extract::expand` keeps 42% of the unitigs.
                //
                // Checked once per retained record rather than on a line counter. A 747-node
                // neighbourhood is 747 records, so a counter with any useful period would
                // never have fired and the budget would have been inert exactly where it is
                // cheapest to enforce.
                if budget_bytes != u64::MAX {
                    let used = graph.resident_bytes();
                    if used > budget_bytes {
                        return Err(Error::TooLarge {
                            accession: accession.clone(),
                            detail: format!(
                                "graph reached {used} bytes while collecting {next} of {} \
                                 neighbourhood unitigs, over the {budget_bytes}-byte \
                                 per-accession graph budget, which is the memory bound \
                                 printed at startup",
                                wanted.len()
                            ),
                            payload: Box::new(crate::types::Payload::default()),
                        });
                    }
                }
            } else {
                keeping = false;
                if h.id > last_wanted {
                    break; // everything wanted has been read
                }
            }
        } else if keeping {
            for &b in &line {
                let code = base_code(b, expected_id - 1).map_err(&bad)?;
                if total_bases & 3 == 0 {
                    graph.seq_2bit.push(0);
                }
                let at = graph.seq_2bit.len() - 1;
                graph.seq_2bit[at] |= code << (2 * (total_bases & 3));
                total_bases += 1;
            }
        }
    }

    if next != wanted.len() {
        return Err(bad(format!(
            "pass 2 found {next} of {} requested unitigs; the two decodes disagree",
            wanted.len()
        )));
    }
    graph.seq_offset.push(
        u32::try_from(total_bases)
            .map_err(|_| bad("retained sequence exceeds u32 addressing".into()))?,
    );
    shrink_adjacency(&mut graph).map_err(|e| bad(format!("reading spooled adjacency: {e}")))?;
    Ok(graph)
}

/// Rewrite the full graph's adjacency down to the neighbourhood's rows, and free the rest.
///
/// `extract::expand` is the last thing that needs an arbitrary unitig's neighbours. After it,
/// `write_gfa`, `phase` and `variants` only ever ask about nodes in the neighbourhood, so
/// holding a row per unitig is holding the one part of the graph that is still O(accession).
/// For `ERR2028820` that is 1,277,363 offsets and 2,582,206 targets, 15.4 MB, against 747
/// rows kept: with eight compute workers it is 15.4 MB each that no longer has to coexist.
///
/// Rows are kept whole, including edges that leave the neighbourhood, because `write_gfa`
/// decides which links to drop and cannot do that from a truncated row.
fn shrink_adjacency(graph: &mut UnitigGraph) -> std::io::Result<()> {
    // A graph that kept every unitig and never spilled is already its own neighbourhood;
    // rewriting it would copy 15 MB to produce the same thing.
    if graph.spool.is_none() && graph.nodes.len() == graph.n_unitigs {
        return Ok(());
    }
    let mut offsets: Vec<u32> = Vec::with_capacity(graph.nodes.len() + 1);
    let mut targets: Vec<u32> = Vec::new();
    let mut row: Vec<u32> = Vec::new();
    let nodes = std::mem::take(&mut graph.nodes);
    for &id in &nodes {
        offsets.push(targets.len() as u32);
        graph.edges_by_id_into(id, &mut row)?;
        targets.extend_from_slice(&row);
    }
    offsets.push(targets.len() as u32);
    graph.nodes = nodes;
    // Assigning drops the old allocations here rather than at the end of the accession, and
    // dropping the spool closes the two unlinked handles, which is what frees the disk.
    graph.offsets = offsets;
    graph.targets = targets;
    graph.spool = None;
    Ok(())
}

/// `KA_CEILING` is empirical, read off three cached accessions, not something Logan
/// publishes in a schema. If a file ever carries a higher value the constant is wrong and
/// every saturation flag drawn from it is wrong with it, so say so loudly rather than
/// silently under-flagging.
fn warn_ka_ceiling(accession: &Accession, max_abundance: f32) {
    if max_abundance > crate::types::KA_CEILING {
        eprintln!(
            "[{accession}] WARNING: ka:f {max_abundance} exceeds the assumed ceiling of {}. \
             `types::KA_CEILING` is measured, not declared by Logan; update it, because the \
             `depth_saturated` and `KaSat` flags are computed against it.",
            crate::types::KA_CEILING
        );
    }
}

/// A link whose target is past the last record the file holds.
///
/// Logan objects are several concatenated zstd frames, and a copy cut exactly at a frame
/// boundary decodes **without error**: zstd has no idea more frames were meant to follow. The
/// graph that comes out is smaller and internally plausible, so the walk runs on it and the
/// accession comes back `no_anchors`, or with a partial gene, with nothing saying the input
/// was half a file. `write_gfa` drops the dangling links in silence, so even the GFA looks
/// well formed.
///
/// The links are what give it away: ids are dense and sequential, so a link to a unitig the
/// file does not hold means records are missing. The walker cannot write such a file itself
/// (the cache entry is renamed into place only once the body is complete, and the hash is
/// checked against the object's ETag), so this is about a cache directory filled by something
/// else: a Galaxy-staged dataset, a partial copy, a shared cache.
fn truncated_error(accession: &Accession, target: u32, n_unitigs: u64) -> Error {
    Error::Parse {
        accession: accession.clone(),
        detail: format!(
            "a link names unitig {target} but the file holds only {n_unitigs} records, so              records are missing: this copy of the object is truncated or incomplete. A zstd              stream cut at a frame boundary decodes cleanly, which is why this is caught here              rather than by the decoder"
        ),
    }
}

/// A file whose records carry no abundance tag parses cleanly and reconstructs, but
/// with every depth-based stage inert: the prune floor is 0 so nothing is pruned, no fork has
/// a majority to correct towards, no branch clears `--min-abund`, and the VCF is empty. None
/// of that is an error, and before 0.2.0 the only trace was `mean_ka: 0.0` in a summary nobody
/// reads per accession. Counted over the whole file in pass 1, which is the only pass that
/// sees every record.
fn warn_missing_abundance(accession: &Accession, missing: u64, total: u64) {
    if missing == 0 {
        return;
    }
    eprintln!(
        "[{accession}] WARNING: {missing} of {total} unitigs carry no ka:f/km:f abundance, so          the prune, the majority correction, bubbles and the variant stages see 0x there.          `abundance_missing` in summary.json carries this count."
    );
}

/// What one header line carries.
struct Header {
    id: u32,
    abundance: f32,
    /// Whether the record carried an abundance tag at all, `ka:f:` or its `km:f:` alias.
    ///
    /// A record with neither keeps abundance 0, which is indistinguishable from a real
    /// zero and silently disables every depth-based stage. `scan` and `build` count these;
    /// `header_light` leaves the flag false when it was not asked for an abundance, because
    /// pass 2 does not count anything.
    has_abundance: bool,
}

/// Whether this header field is an abundance tag.
///
/// Logan's `Unitigs.md` says that "in some large accessions, the `ka:f:[xx]` field is
/// replaced by `km:f:[xx]` ... consider that they are the same information". Reading only
/// `ka:f:` left every unitig of such an accession at abundance 0, which is an error nowhere:
/// the prune, the majority correction, bubble detection and the variant stages all simply
/// have nothing to work with, and the run reports `status: success` with `mean_ka: 0.0` and
/// empty alts. Both tags are five bytes, so the value starts at the same offset.
#[inline]
fn is_abundance_tag(f: &[u8]) -> bool {
    f.starts_with(b"ka:f:") || f.starts_with(b"km:f:")
}

/// The value of an abundance tag, rejecting what cannot be a coverage.
///
/// `parse::<f32>()` accepts `inf` and `nan`, and both then flow into every VAF, the
/// prune floor and `summary.json`, where they are written as bare tokens no JSON parser
/// accepts. A negative abundance is equally meaningless. One definition, so the three decodes
/// of a file cannot disagree about which values it holds.
fn abundance_value(f: &[u8], id: u32, bad: &impl Fn(String) -> Error) -> Result<f32> {
    std::str::from_utf8(&f[5..])
        .ok()
        .and_then(|s| s.parse::<f32>().ok())
        .filter(|v| v.is_finite() && *v >= 0.0)
        .ok_or_else(|| {
            bad(format!(
                "record {id}: bad abundance {:?}",
                String::from_utf8_lossy(f)
            ))
        })
}

/// `collect`'s header reader: the id always, the abundance only when asked, the links never.
///
/// `header` parses every field because pass 1 needs every field. Pass 2 needs almost none:
/// the adjacency is already built, and a record pass 2 will skip needs no abundance. Dropping
/// them takes `collect` down by about 34%.
///
/// **What pass 2 therefore stops validating**, and why that is safe: `header` rejects a
/// malformed link and a non-finite, negative or unparseable `ka:f:`. Pass 1 ran those checks
/// over the same bytes, from the same `Compressed` handle, which neither decode can change,
/// so a file that would fail them has already failed in `scan`. The id rules stay here,
/// because they are what makes a disagreement between the two decodes detectable.
fn header_light(
    line: &[u8],
    expected_id: u64,
    want_abundance: bool,
    bad: &impl Fn(String) -> Error,
) -> Result<Header> {
    // Leading whitespace skipped first, so this splits the name exactly where `header`'s
    // `split(is_ascii_whitespace).filter(non-empty)` would.
    let rest = &line[1..];
    let start = rest.iter().position(|b| !b.is_ascii_whitespace()).unwrap_or(rest.len());
    let name_end = rest[start..]
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .map(|p| p + start)
        .unwrap_or(rest.len());
    let name = &rest[start..name_end];
    let id_bytes = match name.iter().rposition(|&b| b == b'_') {
        Some(p) => &name[p + 1..],
        None => name,
    };
    let id = parse_u32(id_bytes).ok_or_else(|| {
        bad(format!(
            "record {expected_id}: cannot read a unitig id from {:?}",
            String::from_utf8_lossy(name)
        ))
    })?;
    if id as u64 != expected_id {
        return Err(bad(format!(
            "unitig ids are not dense and sequential: expected {expected_id}, found {id}. \
             CSR construction relies on record order; this file needs the fallback path."
        )));
    }
    if !want_abundance {
        return Ok(Header { id, abundance: 0.0, has_abundance: false });
    }
    // Same default as `header`: a record with no `ka:f:`/`km:f:` field has abundance 0.
    let mut abundance = 0.0f32;
    let mut has_abundance = false;
    for f in rest[name_end..].split(|b| b.is_ascii_whitespace()).filter(|f| !f.is_empty()) {
        if is_abundance_tag(f) {
            abundance = abundance_value(f, id, bad)?;
            has_abundance = true;
            break;
        }
    }
    Ok(Header { id, abundance, has_abundance })
}

/// Parse one `>` line, pushing its links into `targets` when given one.
///
/// Shared by `scan` and `build`, so the id rules, the abundance validation and the link
/// packing have one definition. Two decodes of the same file that disagreed about any of them
/// would be a silent corruption rather than an error. `collect` uses `header_light`, which
/// shares the id rules and skips the fields pass 2 has no use for.
fn header(
    line: &[u8],
    expected_id: u64,
    targets: Option<&mut Vec<u32>>,
    bad: &impl Fn(String) -> Error,
) -> Result<Header> {
    let mut fields = line[1..].split(|b| b.is_ascii_whitespace()).filter(|f| !f.is_empty());
    let name = fields
        .next()
        .ok_or_else(|| bad(format!("record {expected_id}: empty header")))?;
    let id_bytes = match name.iter().rposition(|&b| b == b'_') {
        Some(p) => &name[p + 1..],
        None => name,
    };
    let id = parse_u32(id_bytes).ok_or_else(|| {
        bad(format!(
            "record {expected_id}: cannot read a unitig id from {:?}",
            String::from_utf8_lossy(name)
        ))
    })?;
    if id as u64 != expected_id {
        return Err(bad(format!(
            "unitig ids are not dense and sequential: expected {expected_id}, found {id}. \
             CSR construction relies on record order; this file needs the fallback path."
        )));
    }
    if id > ID_MASK {
        return Err(bad(format!(
            "unitig id {id} exceeds the {ID_BITS}-bit packing limit"
        )));
    }

    let mut abundance = 0.0f32;
    let mut has_abundance = false;
    let mut sink = targets;
    for f in fields {
        if is_abundance_tag(f) {
            abundance = abundance_value(f, id, bad)?;
            has_abundance = true;
        } else if f.starts_with(b"L:") {
            let mut parts = f[2..].split(|&b| b == b':');
            let s1 = parts.next().unwrap_or(b"");
            let tgt = parts.next().unwrap_or(b"");
            let s2 = parts.next().unwrap_or(b"");
            let (Some(src_fwd), Some(tid), Some(dst_fwd)) = (strand(s1), parse_u32(tgt), strand(s2))
            else {
                return Err(bad(format!(
                    "record {id}: malformed link {:?}",
                    String::from_utf8_lossy(f)
                )));
            };
            if tid > ID_MASK {
                return Err(bad(format!("record {id}: link target {tid} too large")));
            }
            if let Some(t) = sink.as_deref_mut() {
                t.push(pack_edge(tid, src_fwd, dst_fwd));
            }
        }
    }
    Ok(Header { id, abundance, has_abundance })
}

/// Heap held by the five CSR arrays, including the slack doubling leaves behind.
///
/// Capacity, not length: the slack is resident and is most of what the measured
/// `GRAPH_PEAK_PER_COMPRESSED_BYTE` of 1.7 is made of.
fn csr_bytes(
    offsets: &Vec<u32>,
    targets: &Vec<u32>,
    abundance: &Vec<f32>,
    seq_2bit: &Vec<u8>,
    seq_offset: &Vec<u32>,
) -> u64 {
    let w = size_of::<u32>() as u64;
    offsets.capacity() as u64 * w
        + targets.capacity() as u64 * w
        + abundance.capacity() as u64 * size_of::<f32>() as u64
        + seq_2bit.capacity() as u64
        + seq_offset.capacity() as u64 * w
}

/// Consume the decoder and build the complete `UnitigGraph`.
///
/// # Structural assumptions, all enforced
///
/// Verified against real Logan output (`ERR2028820`, 1,277,362 records): ids run strictly
/// 0, 1, 2, ... with no gaps, and every sequence is on a single line. Strict ordering is
/// what makes CSR construction trivial here: a record's `L:` tokens are all the edges
/// incident to that unitig, so appending them in record order already produces the CSR
/// layout, with no edge sort and no temporary edge list.
///
/// These are assumptions about someone else's output format, so each is checked and
/// produces a clear error rather than a silently wrong graph. Multi-line sequences are
/// handled regardless, since that costs nothing.
/// How often the CSR arrays are measured against the budget, in loop iterations.
///
/// Every line, header or sequence, is one iteration, so this catches both an accession with
/// very many records and one with very few enormous ones. The arrays grow by doubling, so the
/// footprint can at most double between two checks; the refusal therefore happens within one
/// doubling of the budget rather than exactly at it, which is why the printed bound says
/// "up to".
///
/// Short enough that the check actually runs: at 65,536 a graph of fewer than about 32,000
/// records would never be measured at all. The cost is six `capacity()` reads per 4096 lines.
const BUDGET_CHECK_EVERY: u64 = 1 << 12;

pub fn build(
    accession: &Accession,
    reader: Box<dyn Read + '_>,
    budget_bytes: u64,
) -> Result<UnitigGraph> {
    let bad = |detail: String| Error::Parse {
        accession: accession.clone(),
        detail,
    };

    let mut rd = BufReader::with_capacity(1 << 18, reader);
    let mut line: Vec<u8> = Vec::with_capacity(512);

    let mut offsets: Vec<u32> = Vec::new();
    let mut targets: Vec<u32> = Vec::new();
    let mut abundance: Vec<f32> = Vec::new();
    let mut seq_2bit: Vec<u8> = Vec::new();
    let mut seq_offset: Vec<u32> = Vec::new();

    let mut expected_id: u64 = 0;
    let mut total_bases: u64 = 0;
    // As in `scan_in`: absent data and unreadable data are different outcomes.
    let mut bytes_read: u64 = 0;
    let mut abundance_missing: u64 = 0;
    let mut in_record = false;
    let mut iters: u64 = 0;

    loop {
        // Enforce the memory the run told the user it would use. `GRAPH_PEAK_PER_COMPRESSED_BYTE`
        // is a calibration from real accessions, not a property of the format: a dense enough
        // file can build a graph many times its compressed size, and before this check the only
        // thing standing between such a file and the OOM killer was that nobody had uploaded
        // one. Refusing it is `too_large`, the same answer the size cap gives, because it is
        // the same judgement: this accession does not fit the budget the user set.
        iters += 1;
        if iters.is_multiple_of(BUDGET_CHECK_EVERY) && budget_bytes != u64::MAX {
            let used = csr_bytes(&offsets, &targets, &abundance, &seq_2bit, &seq_offset);
            if used > budget_bytes {
                return Err(Error::TooLarge {
                    accession: accession.clone(),
                    detail: format!(
                        "parsed graph reached {used} bytes after {expected_id} unitigs, over \
                         the {budget_bytes}-byte per-accession graph budget, which is \
                         the memory bound printed at startup"
                    ),
                    // Filled in by the caller, which is the only place that knows where the
                    // bytes came from.
                    payload: Box::new(crate::types::Payload::default()),
                });
            }
        }
        line.clear();
        let n = rd
            .read_until(b'\n', &mut line)
            .map_err(|e| bad(format!("read failed after {expected_id} records: {e}")))?;
        if n == 0 {
            break;
        }
        bytes_read += n as u64;
        while matches!(line.last(), Some(b'\n') | Some(b'\r')) {
            line.pop();
        }
        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            // ---- header: id, abundance, edges ----
            // Shared with `scan` and `collect`, so three decodes of the same file cannot
            // disagree about ids, abundance or link packing.
            // Taken before the header's links are appended: `offsets[i]` is where record
            // `i`'s edges start, so it is the count as it stands on entry.
            let edges_before = u32::try_from(targets.len())
                .map_err(|_| bad("edge count exceeds u32 addressing".into()))?;
            let h = header(&line, expected_id, Some(&mut targets), &bad)?;
            offsets.push(edges_before);
            seq_offset.push(
                u32::try_from(total_bases)
                    .map_err(|_| bad("total sequence exceeds u32 addressing".into()))?,
            );
            abundance.push(h.abundance);
            if !h.has_abundance {
                abundance_missing += 1;
            }
            expected_id += 1;
            in_record = true;
        } else {
            // ---- sequence, possibly wrapped over several lines ----
            if !in_record {
                return Err(bad("sequence data before the first header".into()));
            }
            for &b in &line {
                let code = base_code(b, expected_id - 1).map_err(&bad)?;
                if total_bases & 3 == 0 {
                    seq_2bit.push(0);
                }
                let last = seq_2bit.len() - 1;
                seq_2bit[last] |= code << (2 * (total_bases & 3));
                total_bases += 1;
            }
        }
    }

    if expected_id == 0 {
        return Err(no_records_error(accession, bytes_read));
    }
    let n_unitigs = expected_id as usize;
    offsets.push(
        u32::try_from(targets.len()).map_err(|_| bad("edge count exceeds u32 addressing".into()))?,
    );
    seq_offset.push(
        u32::try_from(total_bases).map_err(|_| bad("total sequence exceeds u32 addressing".into()))?,
    );

    if let Some(t) = targets.iter().map(|&p| unpack_edge(p).0).max() {
        if t as u64 >= expected_id {
            return Err(truncated_error(accession, t, expected_id));
        }
    }
    warn_ka_ceiling(
        accession,
        abundance.iter().copied().fold(0.0f32, f32::max),
    );
    warn_missing_abundance(accession, abundance_missing, expected_id);
    let n_links = targets.len() as u64;

    Ok(UnitigGraph {
        accession: accession.clone(),
        n_unitigs,
        offsets,
        targets,
        // `build` keeps everything, so the retained set is the whole graph and `slot` is the
        // identity. The production path is `scan` + `collect`, which keeps the neighbourhood.
        nodes: (0..n_unitigs as UnitigId).collect(),
        abundance,
        seq_2bit,
        seq_offset,
        // `build` never spills: it is the test and fixture path, and it keeps everything.
        spool: None,
        abundance_missing,
        n_links,
    })
}

#[inline]
fn strand(b: &[u8]) -> Option<bool> {
    match b {
        b"+" => Some(true),
        b"-" => Some(false),
        _ => None,
    }
}

#[inline]
fn parse_u32(b: &[u8]) -> Option<u32> {
    if b.is_empty() {
        return None;
    }
    let mut v: u32 = 0;
    for &c in b {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v.checked_mul(10)?.checked_add((c - b'0') as u32)?;
    }
    Some(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_str(s: &str) -> Result<UnitigGraph> {
        build_str_budget(s, u64::MAX)
    }

    fn build_str_budget(s: &str, budget: u64) -> Result<UnitigGraph> {
        build(&"TEST".to_string(), Box::new(std::io::Cursor::new(s.as_bytes().to_vec())), budget)
    }

    /// A graph wide enough to need more than the budget is refused, and refused as
    /// `too_large` rather than as a parse failure.
    ///
    /// `GRAPH_PEAK_PER_COMPRESSED_BYTE` is a calibration from seven real accessions, not a
    /// property of the format. Nothing enforced it, so a file dense enough to build a graph
    /// many times its compressed size would have been parsed until the allocator or the OOM
    /// killer stopped it, whatever the startup line had promised.
    #[test]
    fn a_graph_over_the_budget_is_too_large_not_a_parse_error() {
        // Enough records that the budget check (every 1 << 16 iterations, two per record)
        // runs several times.
        let mut fa = String::new();
        for i in 0..200_000u32 {
            fa.push_str(&format!(">T_{i} ka:f:5.0\nACGTACGTAC\n"));
        }
        match build_str_budget(&fa, 1 << 16) {
            Err(Error::TooLarge { detail, .. }) => {
                assert!(detail.contains("over the"), "{detail}");
                assert!(detail.contains("budget"), "{detail}");
            }
            Err(e) => panic!("wrong error kind, this must be a negative not a failure: {e}"),
            Ok(g) => panic!("parsed {} unitigs with a 64 KiB budget", g.n_unitigs),
        }
        // The same input is fine when nothing is budgeted, so the refusal is the budget and
        // not the data.
        assert_eq!(build_str_budget(&fa, u64::MAX).expect("unbudgeted").n_unitigs, 200_000);
    }

    /// A graph with everything the format can throw at the parser: links on both strands,
    /// sequence wrapped over several lines, a unitig shorter than any k, varied abundance,
    /// and a record with no links at all.
    const MIXED: &str = ">T_0 ka:f:5.5 L:+:1:+ L:-:2:-\n\
                         ACGTACGTAC\n\
                         GTTTACGGCA\n\
                         >T_1 ka:f:17.25 L:-:0:+\n\
                         TTTT\n\
                         >T_2 ka:f:0\n\
                         GGGGCCCCAAAATTTTGGGGCCCC\n\
                         >T_3 ka:f:3.5 L:+:0:- L:+:1:- L:-:2:+\n\
                         ACGT\n";

    fn cursor(s: &str) -> Box<dyn Read> {
        Box::new(std::io::Cursor::new(s.as_bytes().to_vec()))
    }

    /// A Logan object is several concatenated zstd frames, and a copy cut exactly at a
    /// frame boundary decodes **without error**: zstd sees a complete frame and stops. The
    /// graph that comes out is smaller and internally plausible, so without a check the
    /// accession comes back `no_anchors` or with a partial gene and nothing says half the
    /// file was missing.
    ///
    /// The links give it away, because ids are dense and sequential: a link to a record the
    /// file does not hold can only mean records are missing.
    #[test]
    fn a_file_cut_at_a_frame_boundary_is_refused_as_truncated() {
        let first = ">T_0 ka:f:9.0 L:+:1:+\nACGT\n>T_1 ka:f:9.0 L:-:0:- L:+:2:+\nACGT\n";
        let second = ">T_2 ka:f:9.0 L:-:1:- L:+:3:+\nACGT\n>T_3 ka:f:9.0 L:-:2:-\nACGT\n";
        let mut whole = zstd::encode_all(first.as_bytes(), 3).unwrap();
        whole.extend(zstd::encode_all(second.as_bytes(), 3).unwrap());

        // The whole object is four unitigs and parses.
        let all = zstd::decode_all(&whole[..]).expect("two frames decode");
        assert_eq!(
            build(&"T".into(), Box::new(std::io::Cursor::new(all)), u64::MAX)
                .expect("the whole file")
                .n_unitigs,
            4
        );

        // Cut after the first frame. This is the part that makes the defect invisible: zstd
        // reports no error at all, so nothing upstream of the parser can catch it.
        let cut = zstd::encode_all(first.as_bytes(), 3).unwrap();
        let half = zstd::decode_all(&cut[..]).expect("a frame boundary decodes cleanly");
        assert_eq!(half, first.as_bytes(), "the cut decodes to exactly the first frame");

        for pass in ["build", "scan"] {
            let rd: Box<dyn Read> = Box::new(std::io::Cursor::new(half.clone()));
            let got = if pass == "build" {
                build(&"T".into(), rd, u64::MAX).map(|g| g.n_unitigs)
            } else {
                scan(&"T".to_string(), rd, u64::MAX, &mut ()).map(|g| g.n_unitigs)
            };
            match got {
                Err(Error::Parse { detail, .. }) => {
                    assert!(detail.contains("truncated or incomplete"), "{pass}: {detail}");
                    assert!(detail.contains("unitig 2"), "{pass}: name the target: {detail}");
                }
                Err(e) => panic!("{pass}: wrong error kind: {e}"),
                Ok(n) => panic!("{pass}: parsed {n} unitigs from half a file"),
            }
        }
    }

    /// A file with no `L:` token anywhere parses into a graph with no adjacency. That
    /// is legal input and the parser keeps accepting it, because a gene inside a single
    /// unitig needs no edge; what changed is that the count survives to the pipeline, which
    /// reports `no_links` rather than the biological negative `no_path`.
    #[test]
    fn a_file_with_no_links_parses_and_says_it_has_none() {
        let linkless = ">T_0 ka:f:9.0\nACGT\n>T_1 ka:f:9.0\nACGT\n>T_2 ka:f:9.0\nACGT\n";
        let g = build_str(linkless).expect("a link-less file is still parseable");
        assert_eq!(g.n_unitigs, 3);
        assert_eq!(g.n_links, 0, "nothing in the file links anything");
        // Pass 1 is where the production path learns it, and it must agree with `build`.
        let skeleton = scan(&"T".to_string(), cursor(linkless), u64::MAX, &mut ()).expect("scan");
        assert_eq!(skeleton.n_links, 0);

        // The same records with links report them, counted in both endpoints' headers as
        // Logan writes them.
        let linked = ">T_0 ka:f:9.0 L:+:1:+\nACGT\n>T_1 ka:f:9.0 L:-:0:-\nACGT\n";
        assert_eq!(build_str(linked).expect("build").n_links, 2);
        assert_eq!(
            scan(&"T".to_string(), cursor(linked), u64::MAX, &mut ()).expect("scan").n_links,
            2
        );
    }

    /// Logan documents `km:f:` as an alias of `ka:f:` "in some large accessions", and
    /// the parser read only the first. A file using the alias parsed cleanly with every
    /// abundance at 0, which is not an error anywhere downstream: it is a 0x gene, so the
    /// prune floor collapses, no fork has a majority, nothing clears `--min-abund`, and the
    /// run reports success with an empty VCF. The two spellings must be the same file.
    #[test]
    fn km_f_is_read_as_abundance_exactly_like_ka_f() {
        let ka = MIXED;
        let km = MIXED.replace("ka:f:", "km:f:");
        let a = build_str(ka).expect("ka:f");
        let b = build_str(&km).expect("km:f");
        assert_eq!(a.abundance, b.abundance, "the alias must carry the same abundances");
        assert!(a.abundance.iter().any(|v| *v > 0.0), "the fixture must have real abundances");
        assert_eq!(a.offsets, b.offsets, "offsets");
        assert_eq!(a.targets, b.targets, "targets");
        assert_eq!(a.seq_2bit, b.seq_2bit, "sequence");
        // Pass 1 and pass 2 read the alias through different functions (`header` and
        // `header_light`), so a fix to one and not the other would be a silent disagreement
        // between the two decodes of the same file.
        let skeleton = scan(&"T".to_string(), cursor(&km), u64::MAX, &mut ()).expect("scan");
        let all: Vec<UnitigId> = (0..b.n_unitigs as UnitigId).collect();
        let two = collect(skeleton, cursor(&km), &all, u64::MAX).expect("collect");
        assert_eq!(b.abundance, two.abundance, "pass 2 must read the alias too");
    }

    /// The alias is validated exactly as `ka:f:` is: non-finite and negative values are
    /// rejected rather than flowing into the floor, the VAFs and `summary.json`.
    #[test]
    fn a_bad_km_f_value_is_refused_like_a_bad_ka_f() {
        for bad in ["nan", "inf", "-7.5", "x"] {
            let fa = format!(">T_0 km:f:{bad}\nACGT\n");
            match build_str(&fa) {
                Err(Error::Parse { detail, .. }) => {
                    assert!(detail.contains("bad abundance"), "{detail}");
                }
                other => panic!("km:f:{bad} must be refused, got {:?}", other.map(|g| g.abundance)),
            }
        }
    }

    /// A file carrying neither tag is counted, so the summary can say
    /// that the depth-based stages ran against zeros rather than against a 0x gene.
    #[test]
    fn records_with_no_abundance_tag_are_counted_not_guessed() {
        let fa = ">T_0 ka:f:9.0 L:+:1:+\nACGT\n>T_1\nACGT\n>T_2 km:f:4.0\nACGT\n";
        let one = build_str(fa).expect("build");
        assert_eq!(one.abundance_missing, 1, "one of three records has no tag");
        assert_eq!(one.abundance, vec![9.0, 0.0, 4.0]);
        // Counted in pass 1 as well, which is the production path and the only pass that
        // sees every record.
        let skeleton = scan(&"T".to_string(), cursor(fa), u64::MAX, &mut ()).expect("scan");
        assert_eq!(skeleton.abundance_missing, 1);
        // An ordinary object leaves it at zero, which is what keeps the summary field absent.
        assert_eq!(build_str(MIXED).expect("build").abundance_missing, 0);
    }

    /// Two decodes must produce exactly what one did.
    ///
    /// `build` and `scan` + `collect` are separate walks over the same bytes. A disagreement
    /// between them about ids, links, abundance or packing would be a silent corruption, not
    /// an error, and the production path uses the second while almost every test fixture uses
    /// the first. Collecting the whole id range makes them directly comparable.
    #[test]
    fn scan_then_collect_everything_equals_build() {
        let one = build_str(MIXED).expect("build");
        let skeleton = scan(&"T".to_string(), cursor(MIXED), u64::MAX, &mut ()).expect("scan");
        // Pass 1 holds adjacency and nothing else. That is the saving.
        assert!(skeleton.nodes.is_empty(), "pass 1 must not retain per-unitig data");
        assert!(skeleton.seq_2bit.is_empty());
        assert!(skeleton.abundance.is_empty());
        assert_eq!(skeleton.n_unitigs, one.n_unitigs);

        let all: Vec<UnitigId> = (0..one.n_unitigs as UnitigId).collect();
        let two = collect(skeleton, cursor(MIXED), &all, u64::MAX).expect("collect");
        assert_eq!(one.n_unitigs, two.n_unitigs, "n_unitigs");
        assert_eq!(one.offsets, two.offsets, "offsets");
        assert_eq!(one.targets, two.targets, "targets");
        assert_eq!(one.nodes, two.nodes, "nodes");
        assert_eq!(one.abundance, two.abundance, "abundance");
        assert_eq!(one.seq_offset, two.seq_offset, "seq_offset");
        assert_eq!(one.seq_2bit, two.seq_2bit, "seq_2bit");
    }

    /// Collecting a subset yields exactly the sequences and abundances `build` has for those
    /// ids, and nothing for the rest. This is what every downstream stage reads.
    #[test]
    fn collecting_a_subset_matches_the_full_graph_for_those_ids() {
        let one = build_str(MIXED).expect("build");
        let wanted: Vec<UnitigId> = vec![1, 3];
        let skeleton = scan(&"T".to_string(), cursor(MIXED), u64::MAX, &mut ()).expect("scan");
        let sub = collect(skeleton, cursor(MIXED), &wanted, u64::MAX).expect("collect");

        assert_eq!(sub.n_unitigs, one.n_unitigs, "the graph is still that big");
        // Adjacency has been narrowed to the neighbourhood's rows, which is the saving, but
        // each row is whole: `write_gfa` decides which links leave the neighbourhood and
        // cannot do that from a truncated row.
        assert!(sub.targets.len() < one.targets.len(), "adjacency must have been narrowed");
        assert_eq!(sub.offsets.len(), wanted.len() + 1);
        for &id in &wanted {
            assert_eq!(sub.edges(id), one.edges(id), "unitig {id} edges");
        }

        let (mut a, mut b) = (Vec::new(), Vec::new());
        for &id in &wanted {
            one.decode_seq(id, &mut a);
            sub.decode_seq(id, &mut b);
            assert_eq!(a, b, "unitig {id} sequence");
            assert_eq!(one.seq_len(id), sub.seq_len(id), "unitig {id} length");
            assert_eq!(one.abundance_of(id), sub.abundance_of(id), "unitig {id} abundance");
        }
        for id in [0u32, 2] {
            assert!(!sub.has(id), "unitig {id} was not asked for and must not be retained");
        }
        assert!(sub.seq_2bit.len() < one.seq_2bit.len(), "the subset must be smaller");
    }

    /// The two decodes must agree about what is in the file.
    ///
    /// `collect` walks the same bytes `scan` did and takes it on trust that the ids it is
    /// asked for are there. If they are not, the graph comes back with fewer sequences than
    /// nodes and every downstream index is off by an unknown amount, which is a silent
    /// corruption rather than an error. The check that catches it is the only thing standing
    /// between a disagreeing pair of decodes and a wrong haplotype.
    #[test]
    fn collect_refuses_when_the_two_decodes_disagree() {
        let skeleton = scan(&"T".to_string(), cursor(MIXED), u64::MAX, &mut ()).expect("scan");
        // Id 9 is not in this file at all, which is what a disagreement looks like from here.
        let msg = match collect(skeleton, cursor(MIXED), &[1, 9], u64::MAX) {
            Err(e) => e.to_string(),
            Ok(g) => panic!("accepted a graph with {} of 2 nodes", g.nodes.len()),
        };
        assert!(msg.contains("found 1 of 2"), "{msg}");
        assert!(msg.contains("decodes disagree"), "{msg}");
    }

    /// The reader pass 2 uses must parse exactly what the shared one does.
    ///
    /// `header_light` exists only to skip work. The moment it disagrees with `header` about an
    /// id or an abundance, the two decodes are reading different graphs out of the same bytes,
    /// and nothing downstream could tell. These are the shapes where two implementations of
    /// the same rules drift apart: the name split, a name with no `_` and one with two,
    /// whitespace after `>`, `ka:f:` before and after the links, and no `ka:f:` at all.
    #[test]
    fn header_light_agrees_with_header_on_every_shape() {
        let bad = |detail: String| Error::Parse {
            accession: "T".to_string(),
            detail,
        };
        let lines: &[&str] = &[
            ">T_7 ka:f:5.5 L:+:1:+ L:-:2:-",
            ">7 ka:f:0.25",
            ">  T_7  ka:f:12  L:+:3:+",
            ">T_7 L:+:1:+ ka:f:9.75 L:-:2:-",
            ">T_7 L:+:1:+",
            ">SRR000_sub_7 ka:f:1e2",
        ];
        for line in lines {
            let full = header(line.as_bytes(), 7, None, &bad).expect(line);
            let light = header_light(line.as_bytes(), 7, true, &bad).expect(line);
            assert_eq!(full.id, light.id, "id for {line}");
            assert_eq!(full.abundance, light.abundance, "abundance for {line}");
            // Not asking for the abundance must not change which record this is.
            let id_only = header_light(line.as_bytes(), 7, false, &bad).expect(line);
            assert_eq!(full.id, id_only.id, "id without the abundance for {line}");
        }

        // They must also refuse the same lines. Pass 1 would have rejected any of these
        // first, so without this the fast reader's own validation is a guard no test can
        // fail, and a guard like that may as well not be there.
        for line in [
            ">T_7 ka:f:nan",
            ">T_7 ka:f:-1",
            ">T_7 ka:f:x",
            ">notanid ka:f:1",
            // A header with no name at all: `header` calls that an empty header and
            // `header_light` cannot read an id out of it, so they refuse it for different
            // stated reasons. What matters is that neither invents a record 0.
            ">",
            ">   ",
        ] {
            assert!(header(line.as_bytes(), 7, None, &bad).is_err(), "header took {line}");
            assert!(
                header_light(line.as_bytes(), 7, true, &bad).is_err(),
                "header_light took {line}"
            );
        }
    }

    /// Pass 2 kept the id rules, and they have to fire on a record it is skipping.
    ///
    /// Everything else `header_light` checks was dropped on purpose; this is what is left, and
    /// it is what makes a second decode that disagrees with the first an error rather than a
    /// wrong answer. Without it the collection below returns the two nodes it was asked for
    /// and never mentions that the record between them was not the one `scan` saw.
    #[test]
    fn pass_2_refuses_a_file_whose_ids_moved_under_it() {
        let skeleton = scan(&"T".to_string(), cursor(MIXED), u64::MAX, &mut ()).expect("scan");
        // The same four records with the third renumbered: what a decode of changed bytes
        // looks like from here. It is a record `collect` skips, so only `header_light`'s own
        // check can catch it.
        let moved = MIXED.replace(">T_2 ", ">T_9 ");
        assert_ne!(moved, MIXED, "the fixture must actually have been renumbered");
        let msg = match collect(skeleton, cursor(&moved), &[1, 3], u64::MAX) {
            Err(e) => e.to_string(),
            Ok(g) => panic!("accepted {} nodes from a renumbered file", g.nodes.len()),
        };
        assert!(msg.contains("not dense and sequential"), "{msg}");
    }

    /// What pass 2 stopped checking, pass 1 still checks.
    ///
    /// `header_light` does not look at a skipped record's `ka:f:`, nor at anybody's links. The
    /// safety argument for that is entirely that `scan` ran first over the same bytes, so this
    /// test holds up both halves of it: pass 1 refuses each file, and pass 2 walks past the
    /// same record without noticing. If the first half ever stops being true, the second half
    /// is a hole.
    #[test]
    fn the_fields_pass_2_skips_are_still_validated_by_pass_1() {
        for (what, bad_record) in [
            ("abundance", ">T_1 ka:f:nan"),
            ("link", ">T_1 ka:f:1.0 L:?:0:+"),
        ] {
            let good = ">T_0 ka:f:1.0\nACGT\n>T_1 ka:f:1.0\nACGT\n".to_string();
            let bad_fa = format!(">T_0 ka:f:1.0\nACGT\n{bad_record}\nACGT\n");
            let msg = match scan(&"T".to_string(), cursor(&bad_fa), u64::MAX, &mut ()) {
                Err(e) => e.to_string(),
                Ok(_) => panic!("pass 1 must refuse a bad {what}"),
            };
            assert!(msg.contains("bad abundance") || msg.contains("malformed link"), "{msg}");

            // Pass 2 reaches that record only after pass 1 has passed the file, so it is
            // entitled to skip it. The skeleton comes from the file's valid twin, which is the
            // only way to get one at all.
            let skeleton = scan(&"T".to_string(), cursor(&good), u64::MAX, &mut ()).expect("scan");
            let g = collect(skeleton, cursor(&bad_fa), &[0], u64::MAX)
                .unwrap_or_else(|e| panic!("pass 2 must not re-read the {what}: {e}"));
            assert_eq!(g.nodes, vec![0], "only the wanted record is collected");
            assert_eq!(g.abundance_of(0), 1.0, "the kept record's abundance");
        }
    }

    /// Pass 2 is bounded too.
    ///
    /// The neighbourhood is where the sequences land, and it can be most of the graph: on
    /// amplicon-like data `extract::expand` keeps 42% of the unitigs. A budget that bounded
    /// only pass 1's adjacency would hold the printed bound right up until the point the
    /// memory was actually spent.
    #[test]
    fn collect_refuses_a_neighbourhood_over_the_budget() {
        let mut fa = String::new();
        for i in 0..4000u32 {
            fa.push_str(&format!(">T_{i} ka:f:5.0\nACGTACGTACGTACGTACGTACGTACGTACGT\n"));
        }
        let all: Vec<UnitigId> = (0..4000).collect();
        let skeleton = scan(&"T".to_string(), cursor(&fa), u64::MAX, &mut ()).expect("scan");
        match collect(skeleton, cursor(&fa), &all, 1 << 14) {
            Err(Error::TooLarge { detail, .. }) => {
                assert!(detail.contains("neighbourhood unitigs"), "{detail}");
                assert!(detail.contains("budget"), "{detail}");
            }
            Err(e) => panic!("must be the too_large negative, not {e}"),
            Ok(g) => panic!("collected {} nodes inside a 16 KiB budget", g.nodes.len()),
        }
        // The same collection is fine with room, so the refusal is the budget not the data.
        let skeleton = scan(&"T".to_string(), cursor(&fa), u64::MAX, &mut ()).expect("scan");
        assert_eq!(
            collect(skeleton, cursor(&fa), &all, u64::MAX).expect("unbudgeted").nodes.len(),
            4000
        );
    }

    /// The budget check has to run often enough to fire on the graphs it exists for.
    ///
    /// This is the defect `BUDGET_CHECK_EVERY = 1 << 16` had: every line is one iteration, so
    /// a graph of fewer than about 32,000 records was never measured once and the budget did
    /// not apply to it at all. 10,000 link-dense records are about 200 KB of adjacency, which
    /// a 64 KiB budget must refuse; at the old period the check would not have run.
    #[test]
    fn the_budget_fires_on_a_graph_too_small_for_a_long_period() {
        let mut fa = String::new();
        for i in 0..10_000u32 {
            fa.push_str(&format!(">T_{i} ka:f:5.0 L:+:1:+ L:+:2:+ L:-:3:- L:+:4:-\nACGTACGT\n"));
        }
        match scan(&"T".to_string(), cursor(&fa), 1 << 16, &mut ()) {
            Err(Error::TooLarge { detail, .. }) => assert!(detail.contains("budget"), "{detail}"),
            Err(e) => panic!("must be the too_large negative, not {e}"),
            Ok(g) => panic!(
                "scanned {} unitigs ({} bytes of adjacency) inside a 64 KiB budget",
                g.n_unitigs,
                g.resident_bytes()
            ),
        }
        // Unbudgeted, the same file parses, so the refusal is the period and not the data.
        assert_eq!(
            scan(&"T".to_string(), cursor(&fa), u64::MAX, &mut ()).expect("unbudgeted").n_unitigs,
            10_000
        );
    }

    /// Spilling adjacency to disk must change nothing but where the bytes are.
    ///
    /// The spool is read back with `pread` into caller buffers rather than sliced out of a
    /// Vec, so an off-by-one in the row bounds or a mix-up between the two files would give
    /// plausible-looking edges for the wrong unitig. Comparing against the in-memory path over
    /// the same input is what rules that out.
    #[test]
    fn a_spooled_graph_has_the_same_adjacency_as_an_in_memory_one() {
        // Rows must DIFFER from each other, in both their targets and their length. With
        // every record carrying the same two links, reading the wrong row, or reading targets
        // from byte 0, returns identical bytes and the comparison below proves nothing: two
        // mutants of the row addressing survived this test before the fixture was varied.
        let mut fa = String::new();
        for i in 0..40_000u32 {
            fa.push_str(&format!(">T_{i} ka:f:5.0"));
            for j in 0..(i % 4) + 1 {
                fa.push_str(&format!(" L:+:{}:-", (i * 7 + j * 13) % 40_000));
            }
            fa.push_str("\nACGTACGTACGT\n");
        }
        let dir = std::env::temp_dir().join(format!("logan-spool-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let plain = scan(&"T".to_string(), cursor(&fa), u64::MAX, &mut ()).expect("in memory");
        let spooled = scan_in(&"T".to_string(), cursor(&fa), u64::MAX, &mut (), Some(&dir), 1 << 12)
            .expect("spooled");
        assert!(plain.spool.is_none(), "no spool dir was given");
        assert!(spooled.spool.is_some(), "200,000 records must outgrow the spill threshold");
        assert_eq!(plain.n_unitigs, spooled.n_unitigs);
        // The spooled graph holds its adjacency on disk, so it costs almost nothing resident.
        assert!(
            spooled.resident_bytes() * 4 < plain.resident_bytes(),
            "spooled {} bytes against {} in memory",
            spooled.resident_bytes(),
            plain.resident_bytes()
        );

        // Every row must agree, including the empty ones and the last.
        // Every row, not a sample: an addressing bug that only bites some rows is still a
        // wrong haplotype for those accessions.
        let (mut a, mut b) = (Vec::new(), Vec::new());
        let mut lengths = std::collections::BTreeSet::new();
        for id in 0..40_000u32 {
            plain.edges_by_id_into(id, &mut a).expect("memory row");
            spooled.edges_by_id_into(id, &mut b).expect("spooled row");
            assert_eq!(a, b, "unitig {id} edges differ between memory and spool");
            assert!(!a.is_empty(), "fixture should give unitig {id} edges");
            lengths.insert(a.len());
        }
        // Collecting a neighbourhood must read the rows back off disk and then release the
        // handles. Asserting this on a graph that never spilled, which is what an earlier
        // version of this did, proves nothing: `spool` is already `None` there.
        let wanted: Vec<UnitigId> = (0..40_000u32).step_by(997).collect();
        let spooled2 = scan_in(&"T".to_string(), cursor(&fa), u64::MAX, &mut (), Some(&dir), 1 << 12)
            .expect("spooled");
        assert!(spooled2.spool.is_some(), "must have spilled before collect");
        let collected = collect(spooled2, cursor(&fa), &wanted, u64::MAX).expect("collect");
        assert!(
            collected.spool.is_none(),
            "collect must release the spool handles, or every accession leaks two fds"
        );
        assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 0);
        // And the rows it kept are the ones the in-memory graph has.
        for &id in &wanted {
            plain.edges_by_id_into(id, &mut a).expect("memory row");
            assert_eq!(collected.edges(id), a.as_slice(), "unitig {id} after collect");
        }

        // The fixture is only a test of row addressing if the rows actually vary.
        assert!(lengths.len() >= 4, "rows must differ in length, saw {lengths:?}");
        assert_ne!(
            {
                plain.edges_by_id_into(0, &mut a).unwrap();
                a.clone()
            },
            {
                plain.edges_by_id_into(1, &mut b).unwrap();
                b.clone()
            },
            "neighbouring rows must differ, or reading the wrong one is undetectable"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A row read by id must still be that row after `collect` re-indexes by slot.
    ///
    /// `collect` keeps only the neighbourhood's rows and indexes them by slot, but the
    /// pipeline expands the neighbourhood a second time when pruning is rejected, and that
    /// runs after `collect`. Indexing `offsets` by id there walks off the end of an array
    /// that now holds one row per retained node rather than one per unitig: on a
    /// 5,000-accession run all 8 accessions that reached that path panicked, and an id small
    /// enough to land inside the shorter array would have silently returned another unitig's
    /// edges. The neighbourhood here is deliberately sparse so slot and id cannot coincide.
    #[test]
    fn a_row_read_by_id_still_reads_that_row_after_collect() {
        let mut fa = String::new();
        for i in 0..5_000u32 {
            fa.push_str(&format!(">T_{i} ka:f:5.0"));
            for j in 0..(i % 4) + 1 {
                fa.push_str(&format!(" L:+:{}:-", (i * 7 + j * 13) % 5_000));
            }
            fa.push_str("\nACGTACGTACGT\n");
        }

        let before = scan(&"T".to_string(), cursor(&fa), u64::MAX, &mut ()).expect("in memory");
        assert!(before.spool.is_none(), "the in-memory rows are the ones that broke");

        let wanted: Vec<UnitigId> = (0..5_000u32).step_by(97).collect();
        let g = scan(&"T".to_string(), cursor(&fa), u64::MAX, &mut ()).expect("in memory");
        let after = collect(g, cursor(&fa), &wanted, u64::MAX).expect("collect");
        assert!(
            after.nodes.len() * 10 < after.n_unitigs,
            "the fixture must leave slot and id far apart, or indexing by id passes by luck"
        );

        let (mut a, mut b) = (Vec::new(), Vec::new());
        for &id in &wanted {
            before.edges_by_id_into(id, &mut a).expect("row before collect");
            after.edges_by_id_into(id, &mut b).expect("row after collect");
            assert_eq!(a, b, "unitig {id} reads differently after collect");
            assert_eq!(after.edges(id), b.as_slice(), "and must agree with edges()");
        }
    }

    /// The spool files are unlinked at creation, so a run leaves nothing to clean up.
    #[test]
    fn spilling_leaves_no_files_behind() {
        let mut fa = String::new();
        for i in 0..200_000u32 {
            fa.push_str(&format!(">T_{i} ka:f:5.0 L:+:1:+\nACGT\n"));
        }
        let dir = std::env::temp_dir().join(format!("logan-spool-unlink-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        {
            let g = scan_in(&"T".to_string(), cursor(&fa), u64::MAX, &mut (), Some(&dir), 1 << 12)
                .expect("spooled");
            assert!(g.spool.is_some());
            // Already gone from the directory while still readable through the handle.
            assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 0, "spool must be unlinked");
            let mut row = Vec::new();
            g.edges_by_id_into(5, &mut row).expect("still readable");
            assert_eq!(row.len(), 1);
        }
        assert_eq!(std::fs::read_dir(&dir).unwrap().count(), 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A spool directory that cannot be written is a warning, not a failed accession.
    ///
    /// Spilling is an optimisation, so an unwritable temp directory degrades to memory rather
    /// than refusing the work.
    #[test]
    fn an_unusable_spool_directory_falls_back_to_memory() {
        let mut fa = String::new();
        for i in 0..200_000u32 {
            fa.push_str(&format!(">T_{i} ka:f:5.0 L:+:1:+\nACGT\n"));
        }
        // A path whose parent is a file, so `create_dir_all` cannot make it.
        let file = std::env::temp_dir().join(format!("logan-spool-blocked-{}", std::process::id()));
        std::fs::write(&file, b"x").unwrap();
        let blocked = file.join("nope");
        let g = scan_in(&"T".to_string(), cursor(&fa), u64::MAX, &mut (), Some(&blocked), 1 << 12)
            .expect("must still parse");
        assert!(g.spool.is_none(), "it must have fallen back to memory");
        assert_eq!(g.n_unitigs, 200_000);
        let mut row = Vec::new();
        g.edges_by_id_into(5, &mut row).expect("readable from memory");
        assert_eq!(row.len(), 1);
        std::fs::remove_file(&file).ok();
    }

    /// An empty Logan object is absent data, not a broken run; an unreadable one is not.
    ///
    /// Logan serves a 55-byte zstd frame, valid and holding nothing, for accessions it
    /// assembled nothing from, about 0.4% of the bucket. Reporting that as `error` is what a
    /// scheduler reads as "this run broke", and what `--allow-empty` refuses to forgive.
    ///
    /// The second half keeps that from being a blanket excuse: a stream that carried bytes
    /// but no `>` record is a file the walker could not read, so it stays an `error`.
    #[test]
    fn an_empty_unitig_file_is_absent_data_and_a_contentful_one_without_records_is_not() {
        // "   \n" is deliberately not here: a line of spaces survives the newline trim and
        // is read as sequence data, so it fails earlier with "sequence data before the first
        // header". That is also an error, and also right, but it is a different check.
        for (name, input, want_absent) in [
            ("an empty object", "", true),
            ("blank lines only", "\n\n\n", false),
        ] {
            // Both parsers: `scan` is pass 1 of the two-pass path, `build` the single-pass one.
            for pass in ["scan", "build"] {
                let got = if pass == "scan" {
                    scan(&"T".to_string(), cursor(input), u64::MAX, &mut ())
                } else {
                    build(&"T".to_string(), cursor(input), u64::MAX)
                };
                match got {
                    Err(Error::NoLoganData { accession }) => {
                        assert!(want_absent, "{pass}/{name}: absent, but it had content");
                        assert_eq!(accession, "T");
                    }
                    Err(Error::Parse { detail, .. }) => {
                        assert!(!want_absent, "{pass}/{name}: an empty file must be a negative");
                        assert!(detail.contains("no unitig records"), "{detail}");
                        assert!(detail.contains("bytes"), "it must say how much it saw: {detail}");
                    }
                    other => panic!("{pass}/{name}: unexpected {:?}", other.map(|_| ())),
                }
            }
        }
    }

    /// The fallback to RAM is still policed by the budget the run announced.
    ///
    /// Before the budget was re-based, this was the hole in the memory promise: a spool
    /// directory that could not be created meant adjacency stayed in RAM and grew to
    /// `--max-size-mb` x 1.7, while the startup line had announced a figure about fifty times
    /// smaller. Now the budget is the announced one, so the same accession is refused instead.
    /// The warning still names the directory, so the cause is visible either way.
    #[test]
    fn a_fallback_to_memory_over_the_budget_is_refused_rather_than_grown_into() {
        let mut fa = String::new();
        for i in 0..200_000u32 {
            fa.push_str(&format!(">T_{i} ka:f:5.0 L:+:1:+\nACGT\n"));
        }
        // Same unusable path as the test above: the parent is a file.
        let file = std::env::temp_dir().join(format!("logan-fallback-budget-{}", std::process::id()));
        std::fs::write(&file, b"x").unwrap();
        let blocked = file.join("nope");

        // With room, the fallback parses, which is what the sibling test asserts.
        let ok = scan_in(&"T".to_string(), cursor(&fa), u64::MAX, &mut (), Some(&blocked), 1 << 12);
        assert!(ok.is_ok(), "the fallback must still work when the budget allows it");

        // With the budget the startup line would have promised, it does not.
        let err = scan_in(&"T".to_string(), cursor(&fa), 64 * 1024, &mut (), Some(&blocked), 1 << 12);
        match err {
            Err(Error::TooLarge { detail, .. }) => {
                assert!(detail.contains("graph adjacency reached"), "{detail}");
                assert!(detail.contains("per-accession graph budget"), "{detail}");
            }
            Err(e) => panic!("a budget refusal must be a negative, got {e}"),
            Ok(_) => panic!("the fallback grew past the budget the run announced"),
        }
        std::fs::remove_file(&file).ok();
    }

    /// A spool that fills the disk is a budget the run exceeded, not a broken file.
    ///
    /// The two outcomes are not interchangeable. `TooLarge` is a negative: the accession is
    /// skipped, `--allow-empty` tolerates it, and a cohort where one giant accession filled
    /// the disk still reports the others as found. `Parse` is an error: it counts towards
    /// `AllErrored` and exit code 3, so a disk that filled would have told a scheduler the
    /// walker itself was broken.
    ///
    /// Exercised through the classifier rather than through a real `ENOSPC`, which would need
    /// a filesystem sized for the test. Both directions are asserted, because mapping every
    /// write failure to a negative would be just as wrong in the other direction: an EIO is
    /// not a budget decision and must stay loud.
    #[test]
    fn a_spool_write_that_runs_out_of_space_is_a_negative_not_an_error() {
        let acc = "ACC".to_string();
        let full = |raw: i32| spill_error(&acc, "offsets", 42, std::io::Error::from_raw_os_error(raw));

        // ENOSPC, and EDQUOT, which is the same decision made by a quota rather than by the
        // device: both arrive at `write_all` on a spool that cannot grow.
        for raw in [libc::ENOSPC, libc::EDQUOT] {
            match full(raw) {
                Error::TooLarge { accession, detail, .. } => {
                    assert_eq!(accession, acc);
                    assert!(detail.contains("ran out of space"), "{detail}");
                    assert!(detail.contains("after 42 unitigs"), "{detail}");
                    // The message has to say what to do, because the remedy is not obvious:
                    // the accession is within every size cap the run set.
                    assert!(detail.contains("--spool-dir"), "{detail}");
                }
                other => panic!("errno {raw} must be a negative, got {other:?}"),
            }
        }

        // Anything else keeps the error.
        for raw in [libc::EIO, libc::EBADF, libc::EFBIG] {
            match full(raw) {
                Error::Parse { accession, detail } => {
                    assert_eq!(accession, acc);
                    assert!(detail.contains("spill offsets"), "{detail}");
                }
                other => panic!("errno {raw} must stay an error, got {other:?}"),
            }
        }
    }

    /// Asking for an id the graph does not hold is a bug in the caller, and says so.
    #[test]
    fn reading_an_unretained_unitig_panics_with_a_useful_message() {
        let skeleton = scan(&"T".to_string(), cursor(MIXED), u64::MAX, &mut ()).expect("scan");
        let sub = collect(skeleton, cursor(MIXED), &[1], u64::MAX).expect("collect");
        let err = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| sub.abundance_of(2)))
            .expect_err("must panic");
        let msg = err
            .downcast_ref::<String>()
            .map(String::as_str)
            .unwrap_or("");
        assert!(msg.contains("not in this graph"), "unhelpful panic: {msg}");
    }

    /// Capacity, not length.
    ///
    /// The slack doubling leaves behind is resident memory, and it is most of what the
    /// measured `GRAPH_PEAK_PER_COMPRESSED_BYTE` of 1.7 is made of. A budget computed from
    /// `len` would under-count by up to half, so a graph could reach twice the budget the
    /// startup line promised before anything refused it. The over-budget test above cannot
    /// tell the two apart, because its budget is small enough that either measure trips.
    #[test]
    fn the_budget_counts_the_slack_not_only_the_data() {
        let mut offsets: Vec<u32> = Vec::with_capacity(1024);
        offsets.push(1);
        let (targets, abundance, seq_offset): (Vec<u32>, Vec<f32>, Vec<u32>) =
            (Vec::new(), Vec::new(), Vec::new());
        let seq_2bit: Vec<u8> = Vec::new();
        let got = csr_bytes(&offsets, &targets, &abundance, &seq_2bit, &seq_offset);
        assert!(
            got >= 1024 * size_of::<u32>() as u64,
            "counted {got} bytes for a Vec holding one element in a 1024-slot allocation"
        );
    }

    /// What `collect` measures must fit inside what `--max-subgraph-nodes` budgets for it.
    ///
    /// The graph budget is now `max_subgraph_nodes` x `BYTES_PER_SUBGRAPH_NODE`, and `collect`
    /// polices it with `csr_bytes`. If a node costs more in `csr_bytes` than the constant
    /// allows, a neighbourhood inside the node cap would be refused for being over the byte
    /// budget, which is a false negative on legitimate data. The two are measured differently
    /// (`BYTES_PER_SUBGRAPH_NODE` is a slope of peak RSS, `csr_bytes` is a sum of capacities),
    /// so this checks that the margin between them is real rather than assumed.
    #[test]
    fn a_neighbourhood_inside_the_node_cap_fits_the_byte_budget() {
        // A neighbourhood with the shape real ones have: ~2.6 edges per node and unitigs long
        // enough to matter, since the 2-bit sequence is the largest per-node term.
        const N: usize = 50_000;
        const BASES: usize = 200;
        let mut offsets: Vec<u32> = Vec::with_capacity(N + 1);
        let mut targets: Vec<u32> = Vec::with_capacity(N * 3);
        let mut seq_offset: Vec<u32> = Vec::with_capacity(N + 1);
        for i in 0..N {
            offsets.push(targets.len() as u32);
            for e in 0..3 {
                targets.push(((i + e + 1) % N) as u32);
            }
            seq_offset.push((i * BASES / 4) as u32);
        }
        offsets.push(targets.len() as u32);
        seq_offset.push((N * BASES / 4) as u32);
        let abundance: Vec<f32> = vec![5.0; N];
        let seq_2bit: Vec<u8> = vec![0; N * BASES / 4];

        let got = csr_bytes(&offsets, &targets, &abundance, &seq_2bit, &seq_offset);
        let budget = N as u64 * crate::config::BYTES_PER_SUBGRAPH_NODE;
        assert!(
            got < budget,
            "a {N}-node neighbourhood measures {got} bytes against a {budget}-byte budget; \
             BYTES_PER_SUBGRAPH_NODE no longer covers what csr_bytes counts"
        );
        // And the margin is not a rounding accident.
        assert!(got * 2 < budget * 3, "the margin is under 1.5x: {got} against {budget}");
    }

    /// The budget must not refuse what fits. An ordinary small graph parses under a budget
    /// far below the 1.7 GB the defaults imply, so the check costs normal data nothing.
    #[test]
    fn a_graph_inside_the_budget_is_untouched() {
        let fa = ">T_0 ka:f:5.0\nACGTACGT\n>T_1 ka:f:6.0\nTTTTGGGG\n";
        assert_eq!(build_str_budget(fa, 1 << 20).expect("in budget").n_unitigs, 2);
        // And byte for byte the same graph as with no budget at all.
        let unbudgeted = build_str(fa).expect("unbudgeted");
        let budgeted = build_str_budget(fa, 1 << 20).expect("budgeted");
        assert_eq!(unbudgeted.n_unitigs, budgeted.n_unitigs);
        for id in 0..2u32 {
            let (mut a, mut b) = (Vec::new(), Vec::new());
            unbudgeted.decode_seq(id, &mut a);
            budgeted.decode_seq(id, &mut b);
            assert_eq!(a, b, "unitig {id} differs under a budget");
        }
    }

    #[test]
    fn edge_packing_round_trips() {
        for &(t, a, b) in &[(0u32, false, false), (1, true, false), (7, false, true), (ID_MASK, true, true)] {
            assert_eq!(unpack_edge(pack_edge(t, a, b)), (t, a, b));
        }
    }

    #[test]
    fn parses_a_minimal_graph() {
        let g = build_str(">T_0 ka:f:5.3   L:-:1:+  L:+:1:-\nACGT\n>T_1 ka:f:9.0   L:+:0:+\nTTTTAA\n").unwrap();
        assert_eq!(g.n_unitigs, 2);
        assert_eq!(g.abundance, vec![5.3, 9.0]);
        // CSR: node 0 has two edges, node 1 has one.
        assert_eq!(g.offsets, vec![0, 2, 3]);
        assert_eq!(g.edges(0).len(), 2);
        assert_eq!(g.edges(1).len(), 1);
        assert_eq!(unpack_edge(g.edges(0)[0]), (1, false, true));
        assert_eq!(unpack_edge(g.edges(0)[1]), (1, true, false));
        assert_eq!(unpack_edge(g.edges(1)[0]), (0, true, true));
    }

    #[test]
    fn sequences_round_trip_through_2bit() {
        let g = build_str(">T_0\nACGT\n>T_1\nTTTTAAAACCCCGGGGA\n").unwrap();
        let mut out = Vec::new();
        g.decode_seq(0, &mut out);
        assert_eq!(out, b"ACGT");
        g.decode_seq(1, &mut out);
        assert_eq!(out, b"TTTTAAAACCCCGGGGA");
        assert_eq!(g.seq_len(0), 4);
        assert_eq!(g.seq_len(1), 17);
    }

    /// Packing is contiguous, so unitig 1 starts mid-byte. This is the case that a
    /// byte-aligned implementation would get right by accident and this one must get right
    /// on purpose.
    #[test]
    fn unaligned_sequence_boundaries_decode_correctly() {
        let g = build_str(">T_0\nACG\n>T_1\nTTGCA\n>T_2\nG\n>T_3\nCCCCCCCC\n").unwrap();
        let mut out = Vec::new();
        for (i, want) in [(0, "ACG"), (1, "TTGCA"), (2, "G"), (3, "CCCCCCCC")] {
            g.decode_seq(i, &mut out);
            assert_eq!(out, want.as_bytes(), "unitig {i}");
        }
        assert_eq!(g.base_at(1, 2), 2, "base_at must agree with decode_seq (G=2)");
    }

    #[test]
    fn handles_multi_line_sequences() {
        let g = build_str(">T_0\nACGT\nACGT\nAC\n").unwrap();
        let mut out = Vec::new();
        g.decode_seq(0, &mut out);
        assert_eq!(out, b"ACGTACGTAC");
    }

    #[test]
    fn handles_crlf() {
        let g = build_str(">T_0 ka:f:2.0   L:+:0:+\r\nACGT\r\n").unwrap();
        assert_eq!(g.n_unitigs, 1);
        assert_eq!(g.edges(0).len(), 1);
    }

    #[test]
    fn record_without_links_or_abundance_is_fine() {
        let g = build_str(">T_0\nACGT\n").unwrap();
        assert_eq!(g.abundance, vec![0.0]);
        assert!(g.edges(0).is_empty());
    }

    /// CSR construction relies on record order. A gap must be a loud error, never a
    /// silently misaligned graph.
    #[test]
    fn non_sequential_ids_are_rejected() {
        let err = build_str(">T_0\nACGT\n>T_2\nACGT\n").err().expect("should fail");
        assert!(err.to_string().contains("not dense and sequential"), "{err}");
    }

    #[test]
    fn non_acgt_base_is_rejected_not_coerced() {
        let err = build_str(">T_0\nACGN\n").err().expect("should fail");
        assert!(err.to_string().contains("non-ACGT"), "{err}");
    }

    #[test]
    fn malformed_inputs_are_errors() {
        assert!(build_str("").is_err(), "empty input");
        assert!(build_str("ACGT\n").is_err(), "sequence before header");
        assert!(build_str(">T_0 L:x:1:+\nACGT\n").is_err(), "bad strand");
        assert!(build_str(">T_0 ka:f:notanumber\nACGT\n").is_err(), "bad abundance");
        assert!(build_str(">T_zero\nACGT\n").is_err(), "non-numeric id");
    }
}

/// Validation against a real Logan file. Ground truth was measured independently of this parser:
/// `ERR2028820` has 1,277,362 unitigs, 2,582,206 `L:` tokens and 68,745,081 bp.
///
/// `LOGAN_TEST_ZST=/path/to/ERR2028820.unitigs.fa.zst cargo test --release -- --ignored --nocapture`
#[cfg(test)]
mod real_data {
    use super::*;
    use crate::stages::stream;
    use crate::types::FetchedUnitigs;

    #[test]
    #[ignore]
    fn builds_err2028820_and_matches_measured_ground_truth() {
        let Ok(path) = std::env::var("LOGAN_TEST_ZST") else {
            eprintln!("SKIP: set LOGAN_TEST_ZST");
            return;
        };
        let Ok(compressed) = std::fs::read(&path) else {
            eprintln!("SKIP: cannot read {path}");
            return;
        };
        let fetched = FetchedUnitigs {
            index: 0,
            accession: "ERR2028820".to_string(),
            compressed_len: compressed.len(),
            compressed: crate::types::Compressed::Memory(compressed),
            fetch_secs: 0.0,
            payload: crate::types::Payload::default(),
        };

        let t0 = std::time::Instant::now();
        let g = build(
            &fetched.accession,
            stream::decoder(&fetched).expect("decoder"),
            u64::MAX,
        )
        .expect("build");
        let elapsed = t0.elapsed();

        let mem = g.offsets.len() * 4 + g.targets.len() * 4 + g.abundance.len() * 4
            + g.seq_2bit.len() + g.seq_offset.len() * 4;
        eprintln!(
            "  unitigs {:>10}  edges {:>10}  bases {:>11}  {:.2}s  graph {:.1} MB",
            g.n_unitigs,
            g.targets.len(),
            g.seq_offset[g.n_unitigs],
            elapsed.as_secs_f64(),
            mem as f64 / 1_048_576.0
        );

        assert_eq!(g.n_unitigs, 1_277_362, "unitig count");
        assert_eq!(g.targets.len(), 2_582_206, "edge count");
        assert_eq!(g.seq_offset[g.n_unitigs], 68_745_081, "total bases");
        assert_eq!(g.offsets.len(), g.n_unitigs + 1);
        assert_eq!(g.abundance.len(), g.n_unitigs);

        // Known record: >ERR2028820_2933 ka:f:5.3   L:-:443893:+  L:+:1239391:-
        assert!((g.abundance[2933] - 5.3).abs() < 1e-6, "abundance {}", g.abundance[2933]);
        let e: Vec<_> = g.edges(2933).iter().map(|&p| unpack_edge(p)).collect();
        assert_eq!(e, vec![(443893, false, true), (1239391, true, false)]);

        // Sequence content, decoded back out of the 2-bit store.
        let mut out = Vec::new();
        g.decode_seq(2933, &mut out);
        assert!(
            out.starts_with(b"CTGTATACACCCTTCCATTTTTAAAACAACGTAAGAATTAATTAGGAAAG"),
            "unitig 2933 decoded as {:?}",
            String::from_utf8_lossy(&out[..50.min(out.len())])
        );

        // Every unitig is at least k long and the largest matches the measured maximum.
        let mut min_len = usize::MAX;
        let mut max_len = 0usize;
        for i in 0..g.n_unitigs as u32 {
            let l = g.seq_len(i);
            min_len = min_len.min(l);
            max_len = max_len.max(l);
        }
        eprintln!("  unitig length min {min_len} max {max_len}");
        assert_eq!(min_len, 31, "shortest unitig should be exactly k");
        assert_eq!(max_len, 1588, "longest unitig");
    }
}
