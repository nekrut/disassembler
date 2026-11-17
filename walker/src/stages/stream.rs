//! Step 1: fetch one accession's unitigs from the Logan S3 bucket.
//!
//! No `curl … | zstd -dc` to a decompressed file on disk that is then read back: only the
//! compressed bytes are held, and a streaming decoder runs over them.
//!
//! Measured on this hardware: 6.9 MB/s single stream, 17.4 MB/s aggregate across 8 parallel
//! range requests, which is where the link saturates. A median accession is ~71 MB compressed,
//! so this stage has a hard floor of roughly 4.1 s that no amount of code can remove.
//!
//! IN   accession
//! OUT  FetchedUnitigs (compressed bytes only), plus a streaming decoder over them

use crate::config::Config;
use crate::error::{Error, Result};
use crate::types::{Accession, Compressed, FetchedUnitigs, Payload, PayloadSource};
use std::io::Read;
use std::path::{Path, PathBuf};
use std::time::Duration;

/// URL for an accession's unitig file.
pub fn unitigs_url(accession: &str) -> String {
    format!("https://logan-pub.s3.amazonaws.com/u/{accession}/{accession}.unitigs.fa.zst")
}

/// Cache path for an accession's compressed payload.
///
/// Flat layout, one file per accession, byte-identical to the S3 object. Callers use the
/// same function whether they are reading from the cache or writing to it.
pub fn cache_path(dir: &Path, accession: &str) -> PathBuf {
    dir.join(format!("{accession}.unitigs.fa.zst"))
}

/// A pooled HTTP client, shared by every fetch thread.
///
/// Built once in `pipeline::run`. Connection reuse matters: TLS setup is a meaningful share
/// of the small-accession case, and the S3 bucket is a single host, so one pool serves all
/// fetch threads.
pub struct Fetcher {
    agent: ureq::Agent,
    max_compressed_bytes: u64,
    attempts: u32,
    cache_dir: Option<PathBuf>,
    offline: bool,
    /// Where adjacency will spill, so the fetch can refuse an accession whose spool will not
    /// fit. `None` disables the check: no spool directory means adjacency stays in RAM, where
    /// `max_graph_bytes_per_accession` already bounds it.
    spool_dir: Option<PathBuf>,
    /// Spools that can exist at once, which is `peak_graphs_for` the run's own work. Zero
    /// under `--fetch-only`, which parses nothing and so spools nothing.
    spool_slots: usize,
}

impl Fetcher {
    /// A fetcher sized for a pipeline that could fill every slot. Tests and callers with no
    /// accession count in hand.
    pub fn new(cfg: &Config) -> Self {
        Self::new_for(cfg, usize::MAX)
    }

    /// As `new`, told how much work the run actually has.
    ///
    /// Only the spool gate uses it, and it matters there: quoting `compute_workers` concurrent
    /// spools at a two-accession run would refuse on four times the disk the run can use. The
    /// same correction `max_resident_bytes_for` makes to the printed memory bound.
    pub fn new_for(cfg: &Config, n_accessions: usize) -> Self {
        let agent_config = ureq::Agent::config_builder()
            .timeout_global(Some(Duration::from_secs(cfg.fetch_timeout_secs)))
            .timeout_connect(Some(Duration::from_secs(cfg.fetch_connect_timeout_secs)))
            // Inspect the status ourselves so 404 (no Logan data) is distinguishable from
            // 5xx (transient, worth retrying).
            .http_status_as_error(false)
            .build();
        Self {
            agent: ureq::Agent::new_with_config(agent_config),
            max_compressed_bytes: cfg.max_compressed_bytes,
            attempts: cfg.fetch_attempts.max(1),
            cache_dir: cfg.unitigs_cache.clone(),
            offline: cfg.offline,
            spool_dir: cfg.spool_dir.clone(),
            spool_slots: cfg.peak_graphs_for(n_accessions),
        }
    }

    /// Refuse an accession whose spooled adjacency will not fit the spool directory.
    ///
    /// Decided on the accession's real size, never on `--max-size-mb`. The cap is a worst case
    /// that assumes every accession sits at it, and real ones do not: at the 1 GB default it
    /// predicts 4.4 GiB of spool for a cohort that writes a few hundred MB. That number is
    /// honest as a printed bound and dishonest as a gate, because it would refuse runs needing
    /// a fiftieth of it. Here the size is known exactly, from `Content-Length` or from the
    /// cached file, so the prediction is `SPOOL_BYTES_PER_COMPRESSED_BYTE` applied to a real
    /// number.
    ///
    /// Multiplied by `spool_slots` because that many accessions can be parsing at once.
    /// `available_bytes` already sees every spool that exists (they are unlinked, not hidden),
    /// so the only unknown is forward growth, and the assumption is that the other workers
    /// hold accessions about this size. For the homogeneous cohorts Logan runs are made of
    /// that is close; for a mixed one it errs in both directions, which is why the write path
    /// in `parse::spill_error` still has to handle a disk that fills anyway. This is the cheap
    /// layer that catches the obvious case before 11 GB is transferred, not a guarantee.
    fn spool_room(&self, len: u64) -> std::result::Result<(), String> {
        let Some(dir) = self.spool_dir.as_ref() else {
            return Ok(());
        };
        if self.spool_slots == 0 {
            return Ok(());
        }
        let per = crate::config::spool_bytes_for_compressed(len);
        // Below the spill threshold adjacency never reaches disk at all, so there is nothing
        // to refuse and no reason to stat a filesystem.
        if per <= crate::types::SPILL_ADJACENCY_AT as u64 {
            return Ok(());
        }
        // A filesystem that will not answer disables the gate rather than failing the
        // accession: refusing on an unknown is worse than the behaviour this replaces.
        let Some(avail) = crate::types::available_bytes(dir) else {
            return Ok(());
        };
        let need = per.saturating_mul(self.spool_slots as u64);
        if need <= avail {
            return Ok(());
        }
        Err(format!(
            "spooled adjacency for {len} compressed bytes needs about {per} bytes, {need} \
             across {} concurrent worker(s), and {} has {avail} free",
            self.spool_slots,
            dir.display()
        ))
    }

    /// Get one accession's compressed unitig file, as a handle wherever possible.
    ///
    /// **Nothing here reads an object into memory.** A cache hit is opened and streamed; a
    /// download is written to disk as it arrives, into the cache entry when there is one and
    /// into an unlinked file under `--spool-dir` otherwise. Only a caller that configured
    /// neither gets a `Vec`, which no CLI run can be.
    ///
    /// With `--unitigs-cache`: read from `<dir>/<ACC>.unitigs.fa.zst` if present, otherwise
    /// (under `--offline`) fail, otherwise download straight into a tempfile that is renamed
    /// into place once the body is complete. A cached file's contents are trusted, because a
    /// truncated file could only have arrived from a killed writer and the rename prevents
    /// that. **Its size is not**, and that check outlived the reason it was added: see
    /// `the_size_cap_decides_the_same_way_wherever_the_bytes_came_from`.
    ///
    /// Returns the compressed bytes only. Decompression is `decoder`'s job, and happens in
    /// the compute pool so the fetch thread can immediately start the next transfer.
    pub fn fetch(&self, index: usize, accession: &Accession) -> Result<FetchedUnitigs> {
        if let Some(dir) = &self.cache_dir {
            let path = cache_path(dir, accession);
            if path.exists() {
                // Opened, not read. The handle is what crosses the channel and what `decoder`
                // streams from, so a cached accession never costs its own size in RAM. The
                // size is taken from the handle rather than from a separate `stat`, which also
                // closes the window between checking the size and reading the file.
                let mut file = std::fs::File::open(&path).map_err(|e| Error::Fetch {
                    accession: accession.clone(),
                    detail: format!("cache read {}: {e}", path.display()),
                })?;
                let len = file.metadata().map(|m| m.len()).map_err(|e| Error::Fetch {
                    accession: accession.clone(),
                    detail: format!("cache read {}: {e}", path.display()),
                })?;
                // `--max-size-mb` applies here exactly as it applies to a download, and it
                // is no longer about memory: a cache hit costs a file handle, and since 0.2.0
                // nothing in the walk scales with accession size. What it is about now is that
                // the flag names a property of the accession, not of the transport. It is
                // threaded from the Python cohort pipelines into every walk, so if the cache
                // ignored it, a warm run and a cold run of the same command would include
                // different accessions and the summary, which records `max_size_mb` and not
                // the cache's history, could not explain the difference. `disassembler fetch`
                // has its own `--max-size-mb`, so filling a cache at a high cap and walking at
                // a low one is a supported thing to do; this is what makes the walk's cap the
                // one that decides.
                if len > self.max_compressed_bytes {
                    return Err(Error::TooLarge {
                        accession: accession.clone(),
                        detail: format!(
                            "cached file {} is {len} bytes, over --max-size-mb ({} bytes); the \
                             cap applies wherever the bytes came from, so raise it to include \
                             this accession",
                            path.display(),
                            self.max_compressed_bytes
                        ),
                        payload: Box::new(Payload {
                            source: PayloadSource::Cache,
                            bytes: Some(len),
                            ..Payload::default()
                        }),
                    });
                }
                if len == 0 {
                    return Err(Error::Fetch {
                        accession: accession.clone(),
                        detail: format!("cache file is empty: {}", path.display()),
                    });
                }
                // Checked before the file is hashed, which is a full pass over it: an
                // accession that cannot be spooled is not worth reading 11 GB to fingerprint.
                // A cache hit is exactly where this bites, because the bytes are already
                // here and the size cap has already let them through.
                if let Err(detail) = self.spool_room(len) {
                    return Err(Error::TooLarge {
                        accession: accession.clone(),
                        detail,
                        payload: Box::new(Payload {
                            source: PayloadSource::Cache,
                            bytes: Some(len),
                            ..Payload::default()
                        }),
                    });
                }
                // Hashed here rather than at the point of download, so the digest covers
                // every path into the pipeline. A cache hit and an `--offline` run have no
                // HTTP response to take an ETag from, and they are exactly the runs whose
                // inputs are least self-evident. Streamed in chunks and rewound, so this
                // costs a pass over the file, not the file.
                let payload = Payload::read_file(PayloadSource::Cache, &mut file).map_err(|e| {
                    Error::Fetch {
                        accession: accession.clone(),
                        detail: format!("cache read {}: {e}", path.display()),
                    }
                })?;
                return Ok(FetchedUnitigs {
                    index,
                    accession: accession.clone(),
                    compressed_len: len as usize,
                    compressed: Compressed::File(file),
                    fetch_secs: 0.0,
                    payload,
                });
            }
            if self.offline {
                return Err(Error::Fetch {
                    accession: accession.clone(),
                    detail: format!(
                        "cache miss under --offline; expected {}",
                        path.display()
                    ),
                });
            }
        }

        let url = unitigs_url(accession);
        let mut last: String = String::from("no attempt made");

        for attempt in 1..=self.attempts {
            match self.try_once(index, accession, &url) {
                // No cache write here any more. `try_once` streams the body straight into the
                // cache entry, so by the time it returns the bytes are already on disk and
                // already renamed into place. The old shape read the object into a `Vec`,
                // wrote that `Vec` out, and then reopened the file to drop the `Vec`: three
                // copies of the same decision and one of them the size of the accession.
                Ok(f) => return Ok(f),
                Err(FetchAttempt::Permanent(detail)) => {
                    return Err(Error::Fetch {
                        accession: accession.clone(),
                        detail,
                    })
                }
                Err(FetchAttempt::Missing) => {
                    return Err(Error::NoLoganData {
                        accession: accession.clone(),
                    })
                }
                Err(FetchAttempt::TooLarge(detail, payload)) => {
                    return Err(Error::TooLarge {
                        accession: accession.clone(),
                        detail,
                        payload,
                    })
                }
                Err(FetchAttempt::Transient(detail)) => {
                    last = detail;
                    if attempt < self.attempts {
                        // Linear backoff is enough here; S3 throttling is rare and the
                        // dominant transient is a dropped connection.
                        std::thread::sleep(Duration::from_millis(500 * attempt as u64));
                    }
                }
            }
        }
        Err(Error::Fetch {
            accession: accession.clone(),
            detail: format!("giving up after {} attempts: {last}", self.attempts),
        })
    }

    fn try_once(&self, index: usize, accession: &Accession, url: &str) -> std::result::Result<FetchedUnitigs, FetchAttempt> {
        let mut resp = self
            .agent
            .get(url)
            .call()
            .map_err(|e| FetchAttempt::Transient(format!("request failed: {e}")))?;

        let status = resp.status().as_u16();
        match status {
            200 => {}
            404 => return Err(FetchAttempt::Missing),
            // 4xx other than 404 is our fault and will not improve on retry.
            400..=499 => return Err(FetchAttempt::Permanent(format!("HTTP {status}"))),
            _ => return Err(FetchAttempt::Transient(format!("HTTP {status}"))),
        }

        // Provenance off the response, read before the body so it survives a refusal: a
        // `too_large` skip can still say which object was too big, and how big.
        //
        // The ETag is the object's MD5 for everything below S3's 5 GB single-PUT limit, which
        // `--max-size-mb` (default 1024) keeps the walker well inside. Above that limit S3
        // appends `-<parts>` and the value is opaque; it is recorded either way, and the
        // suffix is self-identifying.
        let mut payload = Payload {
            source: PayloadSource::S3,
            etag: header(&resp, ureq::http::header::ETAG)
                .map(|v| v.trim_matches('"').to_string()),
            last_modified: header(&resp, ureq::http::header::LAST_MODIFIED)
                .as_deref()
                .and_then(http_date_to_iso8601),
            fetched_at: Some(crate::io::output::iso8601(std::time::SystemTime::now())),
            ..Payload::default()
        };

        // Refuse oversized accessions before transferring them, when the server says how big.
        if let Some(len) = resp.body().content_length() {
            payload.bytes = Some(len);
            if len > self.max_compressed_bytes {
                return Err(FetchAttempt::TooLarge(
                    format!(
                        "compressed size {len} exceeds max_compressed_bytes {}",
                        self.max_compressed_bytes
                    ),
                    Box::new(payload),
                ));
            }
            // The same pre-transfer refusal, for the disk the parse will need rather than
            // the memory the download will need. Both decided from one `Content-Length`, in
            // the ~50 ms before any body is read.
            if let Err(detail) = self.spool_room(len) {
                return Err(FetchAttempt::TooLarge(detail, Box::new(payload)));
            }
        }

        // Decided before a byte is read, so the bytes have somewhere to go the moment they
        // arrive and are never assembled in RAM first.
        let mut sink = self.open_sink(accession);

        // One pass over the body: counted, hashed and written as it streams. The count is
        // ours rather than ureq's `limit`, so the cap is enforced by code that can say which
        // limit was reached and stop the transfer at exactly that point. The hash is free
        // here, where the bytes are already in a register; computing it afterwards would be
        // a second full pass over a file that can be 11 GB (the nn1 probe measured that pass
        // at 91 s).
        use sha2::{Digest, Sha256};
        let mut hasher = Sha256::new();
        let mut buf = vec![0u8; crate::config::DOWNLOAD_BUFFER_BYTES as usize];
        let mut total: u64 = 0;
        let mut reader = resp.body_mut().as_reader();
        loop {
            let got = match reader.read(&mut buf) {
                Ok(0) => break,
                Ok(n) => n,
                // A dropped connection mid-body. Worth retrying, and the partial sink is
                // discarded so the next attempt starts from nothing.
                Err(e) => {
                    sink.abandon();
                    return Err(FetchAttempt::Transient(format!(
                        "body read failed after {total} bytes: {e}"
                    )));
                }
            };
            total += got as u64;
            // The guard for a missing or lying `Content-Length`. A body of exactly the cap
            // passes, as it does the pre-transfer check above; the transfer stops here rather
            // than running to completion first.
            if total > self.max_compressed_bytes {
                sink.abandon();
                return Err(FetchAttempt::TooLarge(
                    format!(
                        "compressed size exceeds max_compressed_bytes {} (found by \
                         downloading: the server sent no Content-Length, or an untrue one)",
                        self.max_compressed_bytes
                    ),
                    // All that is known here is "over the cap". A Content-Length cannot have
                    // been both present and honest, or the check above would have fired, so
                    // its value is dropped rather than recorded as a size.
                    Box::new(Payload { bytes: None, ..payload.clone() }),
                ));
            }
            hasher.update(&buf[..got]);
            if let Err(e) = sink.write_all(&buf[..got]) {
                let where_ = sink.describe();
                sink.abandon();
                return Err(download_write_error(where_, total, e, &payload));
            }
        }

        if total == 0 {
            sink.abandon();
            return Err(FetchAttempt::Transient("empty body".to_string()));
        }

        // Content-Length is the server's claim; this is the bytes that were actually read,
        // and it is what every downstream stage consumed.
        payload.bytes = Some(total);

        // The spool gate again, for a server that sent no `Content-Length`: then the size is
        // only knowable by downloading, and refusing here still beats failing mid-parse with
        // the transfer already spent.
        if let Err(detail) = self.spool_room(total) {
            sink.abandon();
            return Err(FetchAttempt::TooLarge(detail, Box::new(payload)));
        }

        let mut hex = String::with_capacity(64);
        for b in hasher.finalize() {
            hex.push_str(&format!("{b:02x}"));
        }
        payload.sha256 = Some(hex);

        // The rename into the cache happens here, after the last byte and after every
        // refusal, so a skipped accession never leaves a cache entry behind.
        let compressed = sink.finish().map_err(|e| {
            FetchAttempt::Transient(format!("could not commit the downloaded body: {e}"))
        })?;

        Ok(FetchedUnitigs {
            index,
            accession: accession.clone(),
            compressed_len: total as usize,
            compressed,
            fetch_secs: 0.0,
            payload,
        })
    }

    /// Choose where this download's bytes will land.
    ///
    /// Preference order, and both fallbacks are warnings rather than failures, for the same
    /// reason `scan_in` falls back when a spool cannot be created: a storage problem should
    /// cost performance, not the accession.
    ///
    /// 1. `--unitigs-cache`, where the bytes are wanted anyway, so the download and the cache
    ///    write are the same write.
    /// 2. `--spool-dir`, as an unlinked file. The run still never holds the object.
    /// 3. RAM, which is what every download used to do.
    fn open_sink(&self, accession: &Accession) -> BodySink {
        if let Some(dir) = &self.cache_dir {
            match open_cache_tmp(dir, accession) {
                Ok(sink) => return sink,
                Err(e) => eprintln!(
                    "[{accession}] warning: could not open a cache entry in {} for writing \
                     ({e}); this accession will not be cached",
                    dir.display()
                ),
            }
        }
        if let Some(dir) = &self.spool_dir {
            match crate::types::transient_file(dir, "dl") {
                Ok(file) => return BodySink::Transient(file),
                Err(e) => eprintln!(
                    "[{accession}] warning: could not open a download file in {} ({e}); \
                     holding the compressed object in memory instead",
                    dir.display()
                ),
            }
        }
        BodySink::Memory(Vec::new())
    }
}

/// Classify a failed write of a downloaded body.
///
/// The same decision `parse::spill_error` makes, for the same reason: a disk that fills up is
/// a budget the run exceeded, so it is a `too_large` negative that leaves the other accessions
/// and the exit code alone. Anything else is a real failure, and it is permanent rather than
/// transient because retrying a read-only directory three times helps nobody.
fn download_write_error(
    where_: &str,
    after: u64,
    e: std::io::Error,
    payload: &Payload,
) -> FetchAttempt {
    match e.kind() {
        std::io::ErrorKind::StorageFull | std::io::ErrorKind::QuotaExceeded => {
            FetchAttempt::TooLarge(
                format!(
                    "ran out of space writing the downloaded object to {where_} after {after} \
                     bytes ({e}); lower --fetch-workers or -w, lower --max-size-mb, or point \
                     --unitigs-cache or --spool-dir at a larger filesystem"
                ),
                Box::new(payload.clone()),
            )
        }
        _ => FetchAttempt::Permanent(format!(
            "could not write the downloaded object to {where_} after {after} bytes: {e}"
        )),
    }
}

/// One response header as a `String`, or `None` when it is absent or not valid ASCII.
fn header<B>(resp: &ureq::http::Response<B>, name: ureq::http::header::HeaderName) -> Option<String> {
    resp.headers().get(name)?.to_str().ok().map(str::to_string)
}

/// `Sun, 31 Mar 2024 17:38:36 GMT` -> `2024-03-31T17:38:36Z`.
///
/// RFC 9110 IMF-fixdate, which is the form S3 sends, reformatted rather than parsed into an
/// instant: the components are already the ones `summary.json` wants, so there is no epoch
/// arithmetic and no date crate. Obsolete RFC 850 and asctime formats are not accepted, and
/// neither is anything else: a value that is not exactly this shape returns `None` and the
/// field is omitted, so `last_modified` has one format or is absent. A consumer never has to
/// sniff which it got.
fn http_date_to_iso8601(s: &str) -> Option<String> {
    let rest = s.trim().split_once(", ")?.1;
    let mut it = rest.split(' ');
    let (day, mon, year, time) = (it.next()?, it.next()?, it.next()?, it.next()?);
    if it.next()? != "GMT" || it.next().is_some() {
        return None;
    }
    if day.len() != 2 || year.len() != 4 || time.len() != 8 {
        return None;
    }
    if !year.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let month = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ]
    .iter()
    .position(|m| *m == mon)?
        + 1;
    if !time
        .bytes()
        .enumerate()
        .all(|(i, b)| if i == 2 || i == 5 { b == b':' } else { b.is_ascii_digit() })
    {
        return None;
    }
    // Range-checked so the result is a timestamp rather than a plausible-looking string.
    // 60 seconds is allowed: it is a leap second, not a malformed value.
    let d: u8 = day.parse().ok()?;
    let (h, mi, sec): (u8, u8, u8) = (
        time[0..2].parse().ok()?,
        time[3..5].parse().ok()?,
        time[6..8].parse().ok()?,
    );
    if d == 0 || d > 31 || h > 23 || mi > 59 || sec > 60 {
        return None;
    }
    Some(format!("{year}-{month:02}-{day}T{time}Z"))
}

/// Whether a failed attempt is worth repeating.
enum FetchAttempt {
    /// Transport error or 5xx. Retry.
    Transient(String),
    /// A 4xx other than 404. Retrying cannot help.
    Permanent(String),
    /// HTTP 404: Logan has no unitigs for this accession. Not an error of the run.
    Missing,
    /// Over `--max-size-mb`, by Content-Length or, without one, by the bytes read. Skipped by
    /// the user's cap, not an error. Carries what the response said about the object, so the
    /// skip is recorded with a size and an identity rather than only a sentence.
    TooLarge(String, Box<Payload>),
}

/// Where a download's bytes go as they arrive.
///
/// The object is never assembled in memory first. This is what takes `--max-size-mb` out of
/// the memory bound entirely: before it, a download was `read_to_vec`'d whole, which the nn1
/// probe measured at **10.5 GB of peak RSS on an 11 GB accession**, and at the defaults
/// sixteen of those can be in flight at once.
enum BodySink {
    /// Into `--unitigs-cache`, through a tempfile renamed into place only once the whole
    /// body has arrived. A partial file that survived a crash would look valid to a later
    /// `--offline` run and quietly corrupt the analysis, so nothing is ever renamed early.
    Cache { file: std::fs::File, tmp: PathBuf, final_path: PathBuf },
    /// Into an unlinked file under `--spool-dir`. No cache was configured, so the bytes are
    /// wanted only for this run; they die with the process and need no sweeping.
    Transient(std::fs::File),
    /// Held in RAM, as every download used to be. Only reachable when neither a cache nor a
    /// spool directory exists, which no CLI run can produce.
    Memory(Vec<u8>),
}

impl BodySink {
    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        use std::io::Write;
        match self {
            BodySink::Cache { file, .. } | BodySink::Transient(file) => file.write_all(buf),
            BodySink::Memory(v) => {
                v.extend_from_slice(buf);
                Ok(())
            }
        }
    }

    /// Commit the body and hand back the handle the pipeline will decode from.
    ///
    /// The rename happens here and nowhere else, so the cache only ever gains whole files.
    /// A renamed file keeps its inode, so the handle opened on the tempfile is already a
    /// handle on the final path and no reopen is needed; it is rewound because the caller
    /// just finished writing to it.
    fn finish(self) -> std::io::Result<Compressed> {
        match self {
            BodySink::Cache { mut file, tmp, final_path } => {
                file.flush_and_rewind()?;
                if let Err(e) = std::fs::rename(&tmp, &final_path) {
                    let _ = std::fs::remove_file(&tmp);
                    return Err(e);
                }
                Ok(Compressed::File(file))
            }
            BodySink::Transient(mut file) => {
                file.flush_and_rewind()?;
                Ok(Compressed::File(file))
            }
            BodySink::Memory(v) => Ok(Compressed::Memory(v)),
        }
    }

    /// Give up without committing. The tempfile must not outlive a failed download.
    fn abandon(self) {
        if let BodySink::Cache { tmp, .. } = self {
            let _ = std::fs::remove_file(tmp);
        }
    }

    /// What `payload.source` should say. A download is `S3` wherever its bytes landed.
    fn describe(&self) -> &'static str {
        match self {
            BodySink::Cache { .. } => "cache",
            BodySink::Transient(_) => "spool",
            BodySink::Memory(_) => "memory",
        }
    }
}

/// `flush` then `seek(0)`, which every file-backed sink wants and neither trait provides.
trait FlushAndRewind {
    fn flush_and_rewind(&mut self) -> std::io::Result<()>;
}

impl FlushAndRewind for std::fs::File {
    fn flush_and_rewind(&mut self) -> std::io::Result<()> {
        use std::io::{Seek, SeekFrom, Write};
        self.flush()?;
        self.seek(SeekFrom::Start(0))?;
        Ok(())
    }
}

/// A readable-writable tempfile beside the cache entry it will become.
///
/// Same directory as the final path so the `rename` is same-filesystem and therefore atomic.
/// PID plus nanoseconds is enough uniqueness for one host writing its own cache; two
/// processes racing on the same accession may both download, and one rename wins, which is
/// harmless because the bytes are the same.
fn open_cache_tmp(dir: &Path, accession: &str) -> std::io::Result<BodySink> {
    std::fs::create_dir_all(dir)?;
    let final_path = cache_path(dir, accession);
    let ns = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let tmp = dir.join(format!(".{accession}.{}.{ns}.tmp", std::process::id()));
    let file = std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create_new(true)
        .open(&tmp)?;
    Ok(BodySink::Cache { file, tmp, final_path })
}

/// A streaming zstd decoder over the compressed buffer.
///
/// Deliberately a `Read`, not a `Vec<u8>`. `ERR2028820` decompresses to 145 MB and the
/// largest ERG11 accessions to about 1.44 GB; materialising that per worker is what makes
/// concurrency expensive. `stages::parse` consumes this once and keeps only what it needs.
///
/// # Frames
///
/// A zstd file may hold several concatenated frames. This decoder reads through all of
/// them, which matters because stopping at the first would silently truncate the graph
/// rather than fail. `single_frame()` is deliberately NOT called. Covered by
/// `decodes_concatenated_frames`.
///
/// # Where errors appear
///
/// Only the frame header is validated here, so this returns `Ok` for a stream that is
/// corrupt further in. Truncated or damaged data surfaces as an `io::Error` from `read()`
/// during `parse::build`. Callers must not treat a short read as a clean end of file.
///
/// # Reading
///
/// The returned reader is unbuffered on the output side. Read in large chunks
/// (`zstd::Decoder::recommended_output_size()`, ~128 KB) or wrap in a `BufReader`;
/// byte-at-a-time reads through `dyn Read` are slow.
pub fn decoder<'a>(fetched: &'a FetchedUnitigs) -> Result<Box<dyn Read + 'a>> {
    let fail = |e: std::io::Error| Error::Decompress {
        accession: fetched.accession.clone(),
        detail: format!("could not start zstd decoder: {e}"),
    };
    match &fetched.compressed {
        // Cursor<&[u8]> is already BufRead, so `with_buffer` avoids the redundant BufReader
        // that `Decoder::new` would insert over an in-memory slice.
        Compressed::Memory(v) => {
            let cursor = std::io::Cursor::new(v.as_slice());
            Ok(Box::new(
                zstd::stream::read::Decoder::with_buffer(cursor).map_err(fail)?,
            ))
        }
        // `&File` is `Read`, so the file is streamed rather than loaded: this is the whole
        // point of `Compressed::File`. Rewound first, because the handle is shared with
        // whoever hashed it and a caller may legitimately decode twice.
        Compressed::File(f) => {
            let mut handle = f;
            use std::io::{Seek, SeekFrom};
            handle.seek(SeekFrom::Start(0)).map_err(fail)?;
            Ok(Box::new(
                zstd::stream::read::Decoder::new(handle).map_err(fail)?,
            ))
        }
    }
}

/// Suggested read size for the decoder returned by [`decoder`].
pub fn recommended_chunk() -> usize {
    zstd::stream::read::Decoder::<'_, std::io::Cursor<&[u8]>>::recommended_output_size()
}

#[cfg(test)]
mod date_tests {
    use super::http_date_to_iso8601 as parse;

    /// The form S3 actually sends, taken from a live `HEAD` of
    /// `u/SRR8209082/SRR8209082.unitigs.fa.zst`.
    #[test]
    fn reformats_the_imf_fixdate_s3_sends() {
        assert_eq!(
            parse("Sun, 31 Mar 2024 17:38:36 GMT").as_deref(),
            Some("2024-03-31T17:38:36Z")
        );
        // Single-digit month, and the leap second HTTP permits.
        assert_eq!(parse("Wed, 01 Jan 2025 00:00:00 GMT").as_deref(), Some("2025-01-01T00:00:00Z"));
        assert_eq!(parse("Tue, 30 Jun 2026 23:59:60 GMT").as_deref(), Some("2026-06-30T23:59:60Z"));
    }

    /// Anything else is dropped rather than stored raw, so `last_modified` has exactly one
    /// format or is absent. A consumer never has to sniff which it got.
    #[test]
    fn refuses_everything_that_is_not_that_shape() {
        for bad in [
            "",
            "not a date",
            "2024-03-31T17:38:36Z",                   // already ISO: not what S3 sends
            "Sunday, 31-Mar-24 17:38:36 GMT",         // RFC 850
            "Sun Mar 31 17:38:36 2024",               // asctime
            "Sun, 31 Mar 2024 17:38:36 UTC",          // wrong zone token
            "Sun, 31 Mar 2024 17:38:36 GMT extra",    // trailing junk
            "Sun, 31 Zzz 2024 17:38:36 GMT",          // bad month
            "Sun, 3 Mar 2024 17:38:36 GMT",           // one-digit day
            "Sun, 31 Mar 24 17:38:36 GMT",            // two-digit year
            "Sun, 00 Mar 2024 17:38:36 GMT",          // day 0
            "Sun, 32 Mar 2024 17:38:36 GMT",          // day 32
            "Sun, 31 Mar 2024 24:00:00 GMT",          // hour 24
            "Sun, 31 Mar 2024 17:60:00 GMT",          // minute 60
            "Sun, 31 Mar 2024 17:38:61 GMT",          // second 61
            "Sun, 31 Mar 2024 17-38-36 GMT",          // wrong separators
        ] {
            assert_eq!(parse(bad), None, "accepted {bad:?}");
        }
    }
}

#[cfg(test)]
mod url_tests {
    use super::*;

    #[test]
    fn builds_the_logan_url() {
        assert_eq!(
            unitigs_url("ERR2028820"),
            "https://logan-pub.s3.amazonaws.com/u/ERR2028820/ERR2028820.unitigs.fa.zst"
        );
    }
}

/// Live network checks. Ignored by default so `cargo test` stays offline and deterministic.
/// Run with `cargo test --release -- --ignored --nocapture`.
#[cfg(test)]
mod network_tests {
    use super::*;
    use crate::config::Config;

    #[test]
    #[ignore]
    fn fetches_a_real_accession_and_it_decodes() {
        let f = Fetcher::new(&Config::default());
        let got = f.fetch(0, &"SRR8209082".to_string()).expect("fetch");
        eprintln!("fetched {} compressed bytes", got.compressed_len);
        assert_eq!(got.compressed_len, 1_861_058, "unexpected compressed size");
        assert_eq!(
            got.compressed_len,
            got.compressed
                .as_slice()
                .expect("neither a cache nor a spool dir was set, so this is in memory")
                .len()
        );

        // The bytes must be usable by the next stage, not merely non-empty.
        let mut out = Vec::new();
        decoder(&got).expect("decoder").read_to_end(&mut out).expect("decode");
        assert_eq!(out.len(), 9_873_585);
        assert!(out.starts_with(b">SRR8209082_"));
    }

    #[test]
    #[ignore]
    fn missing_accession_is_a_permanent_404_not_a_retry_storm() {
        let f = Fetcher::new(&Config::default());
        let started = std::time::Instant::now();
        let err = f.fetch(0, &"NOTAREALACCESSION999".to_string()).err().expect("should fail");
        eprintln!("{err}");
        // Matched on the variant, not on the message. This asserted `msg.contains("404")`
        // until 0.2.0, when the bucket's 404 and its 55-byte empty object became one fact
        // and the Display text stopped naming the HTTP status: the test then failed on every
        // golden run, and because it is `#[ignore]` nothing offline noticed. The variant is
        // what the pipeline branches on, so it is also what this should pin.
        assert!(
            matches!(err, Error::NoLoganData { .. }),
            "a missing accession must be absent Logan data, got: {err}"
        );
        // A permanent failure must not sleep through the retry schedule.
        assert!(started.elapsed() < Duration::from_secs(5), "404 should fail fast");
    }

    #[test]
    #[ignore]
    fn oversized_accession_is_refused_before_transfer() {
        // SRR12016917 is ~3.9 GB compressed. With a small cap this must be rejected from
        // Content-Length alone, without downloading the body.
        let cfg = Config {
            max_compressed_bytes: 64 * 1024 * 1024,
            ..Config::default()
        };
        let f = Fetcher::new(&cfg);
        let started = std::time::Instant::now();
        let err = f.fetch(0, &"SRR12016917".to_string()).err().expect("should fail");
        let msg = err.to_string();
        eprintln!("{msg}");
        assert!(msg.contains("exceeds max_compressed_bytes"), "got: {msg}");
        assert!(
            started.elapsed() < Duration::from_secs(20),
            "must reject on Content-Length, not after downloading 3.9 GB"
        );
    }
}

/// The size cap against a local server, for the case the live test cannot reach: a response
/// without a Content-Length.
#[cfg(test)]
mod size_cap_tests {
    use super::*;
    use crate::config::Config;
    use std::io::{BufRead, BufReader, Write};

    /// Serve one request with `body` and no Content-Length (the body ends when the
    /// connection closes, as with a chunked or HTTP/1.0 response). Returns the URL.
    fn serve_without_length(body: Vec<u8>) -> String {
        let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
        let port = listener.local_addr().unwrap().port();
        std::thread::spawn(move || {
            let (mut sock, _) = listener.accept().unwrap();
            let mut reader = BufReader::new(sock.try_clone().unwrap());
            let mut line = String::new();
            while reader.read_line(&mut line).unwrap_or(0) > 0 && line != "\r\n" {
                line.clear();
            }
            let _ = sock.write_all(b"HTTP/1.1 200 OK\r\nConnection: close\r\n\r\n");
            let _ = sock.write_all(&body);
        });
        format!("http://127.0.0.1:{port}/u/X/X.unitigs.fa.zst")
    }

    /// A fetcher that talks to the loopback server directly.
    ///
    /// TW-09. `Fetcher::new` leaves ureq's `Proxy::try_from_env()` in place, which is correct
    /// in production. In these tests it means that with `HTTPS_PROXY` set and no `NO_PROXY`
    /// entry for loopback, the request goes to the proxy and both tests fail with "refused at
    /// the cap", which reads like a size-cap regression rather than an environment problem.
    fn fetcher(cap: u64) -> Fetcher {
        let mut f = Fetcher::new(&Config { max_compressed_bytes: cap, ..Config::default() });
        f.agent = ureq::Agent::new_with_config(
            ureq::Agent::config_builder()
                .proxy(None)
                .http_status_as_error(false)
                .build(),
        );
        f
    }

    #[test]
    fn a_body_over_the_cap_without_content_length_is_too_large_not_transient() {
        let url = serve_without_length(vec![b'A'; 100]);
        match fetcher(10).try_once(0, &"X".to_string(), &url) {
            Err(FetchAttempt::TooLarge(d, _)) => assert!(d.contains("no Content-Length"), "{d}"),
            Err(FetchAttempt::Transient(d)) => panic!("retried as transient: {d}"),
            Err(_) => panic!("wrong failure"),
            Ok(f) => panic!("accepted {} bytes over a 10-byte cap", f.compressed_len),
        }
    }

    /// Serve one request with `extra` headers and a correct Content-Length, so the size
    /// check runs on the headers rather than on the transfer.
    fn serve_once(extra: &str, body: Vec<u8>) -> String {
        let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
        let port = listener.local_addr().unwrap().port();
        let head = format!(
            "HTTP/1.1 200 OK\r\nContent-Length: {}\r\n{extra}Connection: close\r\n\r\n",
            body.len()
        );
        std::thread::spawn(move || {
            let (mut sock, _) = listener.accept().unwrap();
            let mut reader = BufReader::new(sock.try_clone().unwrap());
            let mut line = String::new();
            while reader.read_line(&mut line).unwrap_or(0) > 0 && line != "\r\n" {
                line.clear();
            }
            let _ = sock.write_all(head.as_bytes());
            let _ = sock.write_all(&body);
        });
        format!("http://127.0.0.1:{port}/u/X/X.unitigs.fa.zst")
    }

    /// The Logan objects are unversioned and mutable, so a run has to record which bytes it
    /// read. Everything available from the response is captured, and the digest is taken
    /// from the bytes themselves rather than trusted from the header.
    #[test]
    fn a_live_fetch_records_the_object_it_read() {
        let body = vec![b'Z'; 32];
        let url = serve_once(
            "ETag: \"d2dd288b238c3c88b46acdcd78860dce\"\r\nLast-Modified: Sun, 31 Mar 2024 17:38:36 GMT\r\n",
            body.clone(),
        );
        let got = fetcher(1024).try_once(0, &"X".to_string(), &url).ok().expect("fetch");
        let p = &got.payload;
        assert_eq!(p.source, PayloadSource::S3);
        assert_eq!(p.bytes, Some(32));
        assert_eq!(p.sha256.as_deref(), Some(crate::types::sha256_hex(&body).as_str()));
        // Quotes stripped: the value is compared with other ETags, not with a header line.
        assert_eq!(p.etag.as_deref(), Some("d2dd288b238c3c88b46acdcd78860dce"));
        // Normalised, so two summaries sort and diff against each other.
        assert_eq!(p.last_modified.as_deref(), Some("2024-03-31T17:38:36Z"));
        assert!(p.fetched_at.is_some(), "a live fetch records when it read the object");
    }

    /// A response with neither header still yields a usable record: the digest is the field
    /// that does not depend on the transport.
    #[test]
    fn a_response_without_those_headers_still_hashes_what_it_read() {
        let body = vec![b'Q'; 8];
        let got = fetcher(1024)
            .try_once(0, &"X".to_string(), &serve_once("", body.clone()))
            .ok()
            .expect("fetch");
        assert_eq!(got.payload.sha256.as_deref(), Some(crate::types::sha256_hex(&body).as_str()));
        assert_eq!(got.payload.etag, None);
        assert_eq!(got.payload.last_modified, None);
    }

    /// The one negative where the walker knows something about an object it did not read.
    /// Until 0.2.0 the size survived only as prose inside the error string.
    #[test]
    fn a_refusal_on_content_length_records_what_was_refused() {
        let url = serve_once("ETag: \"abc123\"\r\n", vec![b'Z'; 100]);
        match fetcher(10).try_once(0, &"X".to_string(), &url) {
            Err(FetchAttempt::TooLarge(_, p)) => {
                assert_eq!(p.bytes, Some(100), "the size that tripped the cap must be a number");
                assert_eq!(p.etag.as_deref(), Some("abc123"));
                assert_eq!(p.source, PayloadSource::S3);
                assert!(p.sha256.is_none(), "no bytes were read, so there is nothing to hash");
            }
            Err(FetchAttempt::Transient(d)) => panic!("retried as transient: {d}"),
            Err(_) => panic!("wrong failure"),
            Ok(_) => panic!("accepted 100 bytes over a 10-byte cap"),
        }
    }

    #[test]
    fn a_body_of_exactly_the_cap_is_accepted() {
        let url = serve_without_length(vec![b'A'; 10]);
        let got = fetcher(10).try_once(0, &"X".to_string(), &url).ok().expect("refused at the cap");
        assert_eq!(got.compressed_len, 10);
    }
}

#[cfg(test)]
mod cache_tests {
    use super::*;
    use crate::config::Config;

    /// What the next stage would actually get. A cache hit now hands on a file handle rather
    /// than a buffer, so comparing raw bytes no longer proves anything about the handle;
    /// decoding through `decoder` proves it is positioned at the start of the right file.
    fn decoded(f: &FetchedUnitigs) -> Vec<u8> {
        let mut out = Vec::new();
        decoder(f).expect("decoder").read_to_end(&mut out).expect("decode");
        out
    }

    /// Random-ish tempdir under the system temp root; no external `tempfile` dep.
    fn tempdir(tag: &str) -> PathBuf {
        let ns = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        let dir = std::env::temp_dir().join(format!(
            "logan-walker-{tag}-{}-{ns}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).expect("create tempdir");
        dir
    }

    /// The size cap decides the same way wherever the bytes came from.
    ///
    /// This check was added because `--max-size-mb` bounded memory and a cached file was read
    /// whole. Neither is true since 0.2.0: a cache hit is a file handle, and nothing in the
    /// walk scales with accession size. The check stays for a different and stronger reason.
    /// `--max-size-mb` names a property of the accession, not of the transport, and the
    /// Python cohort pipelines thread it into every walk; if the cache ignored it, the same
    /// command over the same list would include different accessions depending on what
    /// happened to be cached, and `summary.json` records the cap but not the cache's history,
    /// so nothing in the output could explain the difference.
    ///
    /// Three cases, because the property is about agreement, not about one refusal: the
    /// boundary is exact, `--offline` does not relax it, and neither does the origin being
    /// reachable (nothing here has a network, so a fall-through to S3 would hang or fail
    /// rather than quietly pass).
    #[test]
    fn the_size_cap_decides_the_same_way_wherever_the_bytes_came_from() {
        let dir = tempdir("cap");
        let bytes = zstd::encode_all(&b">CACHED0002_0\nACGTACGTACGT\n"[..], 3).unwrap();
        std::fs::write(cache_path(&dir, "CACHED0002"), &bytes).unwrap();
        let fetcher = |cap: u64, offline: bool| Fetcher::new(&Config {
            unitigs_cache: Some(dir.clone()),
            offline,
            max_compressed_bytes: cap,
            ..Config::default()
        });
        let acc = "CACHED0002".to_string();
        let len = bytes.len() as u64;

        for offline in [true, false] {
            match fetcher(len - 1, offline).fetch(0, &acc) {
                Err(Error::TooLarge { detail, payload, .. }) => {
                    assert!(detail.contains("cached file"), "{detail}");
                    assert!(detail.contains("--max-size-mb"), "the skip must name the flag: {detail}");
                    // The summary has to be able to explain the skip on its own.
                    assert_eq!(payload.bytes, Some(len));
                    assert_eq!(payload.source, PayloadSource::Cache);
                }
                other => panic!(
                    "offline={offline}: expected TooLarge, got {:?}",
                    other.map(|f| f.compressed_len)
                ),
            }
        }

        // Exact at the boundary, the same as the `Content-Length` check on a download: a file
        // of exactly the cap passes.
        let got = fetcher(len, true).fetch(0, &acc).unwrap();
        assert_eq!(decoded(&got), b">CACHED0002_0\nACGTACGTACGT\n");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// An accession whose spool will not fit is refused as a negative, before the bytes are
    /// read, hashed or parsed.
    ///
    /// The cache hit is where this matters most: the size cap has already let the object
    /// through, the bytes are already local, and nothing else between here and `scan_in`
    /// looks at disk at all. Without the gate the run would download or open an 11 GB
    /// accession, parse for ten minutes, and fail on a write.
    ///
    /// The demand is made unsatisfiable with the worker count rather than with a small
    /// filesystem, which a unit test cannot create: `usize::MAX` concurrent spools saturates
    /// the multiplication, so the refusal does not depend on how much space the machine
    /// running the suite happens to have. What the test still proves is that the gate is what
    /// refused, because the identical file and cap pass once `spool_dir` is `None`.
    #[test]
    fn an_accession_whose_spool_will_not_fit_is_too_large() {
        let dir = tempdir("spool-gate");
        // Over the spill threshold once multiplied by SPOOL_BYTES_PER_COMPRESSED_BYTE: a
        // smaller file never reaches disk and must not be gated. Contents are irrelevant,
        // because the refusal happens before anything decodes them.
        let big = vec![0u8; 16 * 1024 * 1024];
        std::fs::write(cache_path(&dir, "BIGACC0001"), &big).unwrap();
        let acc = "BIGACC0001".to_string();

        let base = Config {
            unitigs_cache: Some(dir.clone()),
            offline: true,
            spool_dir: Some(dir.clone()),
            compute_workers: usize::MAX,
            ..Config::default()
        };

        match Fetcher::new(&base).fetch(0, &acc) {
            Err(Error::TooLarge { detail, payload, .. }) => {
                assert!(detail.contains("spooled adjacency"), "{detail}");
                assert!(detail.contains("free"), "{detail}");
                // The payload is what makes the skip actionable, the same as for the size
                // cap: the size that tripped it, as a number, and where the bytes came from.
                assert_eq!(payload.bytes, Some(big.len() as u64));
                assert_eq!(payload.source, PayloadSource::Cache);
            }
            other => panic!("expected TooLarge, got {:?}", other.map(|f| f.compressed_len)),
        }

        // Same file, same cap, no spool directory: adjacency stays in RAM, where the graph
        // budget bounds it, so there is nothing for this gate to decide.
        let no_spool = Config { spool_dir: None, ..base.clone() };
        let got = Fetcher::new(&no_spool).fetch(0, &acc).expect("no spool dir, no gate");
        assert_eq!(got.compressed_len, big.len());

        // And `--fetch-only` parses nothing, so it spools nothing: the gate must not refuse
        // the one mode whose whole job is to put large accessions in the cache.
        let fetch_only = Config { fetch_only: true, ..base.clone() };
        let got = Fetcher::new(&fetch_only).fetch(0, &acc).expect("--fetch-only never spools");
        assert_eq!(got.compressed_len, big.len());

        // A spool directory nothing can be learned about disables the gate rather than
        // refusing on it. An empty path has no existing ancestor to walk up to, so
        // `available_bytes` cannot answer; treating that as zero free would refuse every
        // accession on a filesystem the walker simply failed to stat.
        let unknowable = Config { spool_dir: Some(PathBuf::from("")), ..base.clone() };
        let got = Fetcher::new(&unknowable).fetch(0, &acc).expect("an unanswerable fs must not refuse");
        assert_eq!(got.compressed_len, big.len());

        std::fs::remove_dir_all(&dir).ok();
    }

    /// Below the spill threshold adjacency never reaches disk, so the gate must not fire
    /// however little space there is. Without this the check would refuse every accession on
    /// a nearly full filesystem, including the ones that would have parsed entirely in RAM.
    #[test]
    fn an_accession_too_small_to_spill_is_never_gated() {
        let dir = tempdir("spool-gate-small");
        let small = vec![0u8; 1024 * 1024];
        std::fs::write(cache_path(&dir, "SMALLACC01"), &small).unwrap();
        // The same unsatisfiable demand as above; only the file size differs.
        let cfg = Config {
            unitigs_cache: Some(dir.clone()),
            offline: true,
            spool_dir: Some(dir.clone()),
            compute_workers: usize::MAX,
            ..Config::default()
        };
        let got = Fetcher::new(&cfg).fetch(0, &"SMALLACC01".to_string()).expect("must not gate");
        assert_eq!(got.compressed_len, small.len());
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A short list cannot fill the pipeline, so it cannot hold `compute_workers` spools at
    /// once, and quoting that many would refuse on disk the run has no way to use.
    #[test]
    fn the_gate_counts_the_spools_the_run_can_actually_hold() {
        let dir = tempdir("spool-gate-slots");
        let cfg = Config {
            spool_dir: Some(dir.clone()),
            compute_workers: 8,
            ..Config::default()
        };
        assert_eq!(Fetcher::new_for(&cfg, 2).spool_slots, 2);
        assert_eq!(Fetcher::new_for(&cfg, 64).spool_slots, 8);
        assert_eq!(Fetcher::new(&cfg).spool_slots, 8);
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A cache hit returns the exact bytes on disk. Verified against a preloaded file so
    /// nothing about this test needs the network.
    #[test]
    fn a_cache_hit_returns_the_bytes_on_disk() {
        let dir = tempdir("hit");
        let bytes = zstd::encode_all(&b">CACHED0001_0\nACGT\n"[..], 3).unwrap();
        std::fs::write(cache_path(&dir, "CACHED0001"), &bytes).unwrap();

        let cfg = Config {
            unitigs_cache: Some(dir.clone()),
            offline: true,
            ..Config::default()
        };
        let f = Fetcher::new(&cfg);
        let got = f.fetch(0, &"CACHED0001".to_string()).expect("cache hit");
        assert_eq!(decoded(&got), b">CACHED0001_0\nACGT\n");
        assert_eq!(got.accession, "CACHED0001");
        assert_eq!(got.compressed_len, bytes.len());
        // The point of the change: a cached accession costs a file handle, not its own size
        // in RAM, and it is the compute worker holding one of these for the whole parse that
        // used to dominate the bound.
        assert!(
            matches!(got.compressed, Compressed::File(_)),
            "a cache hit must not be read into memory"
        );
        assert!(got.compressed.as_slice().is_none());
        // Decoding twice must work: the handle is rewound, not consumed.
        assert_eq!(decoded(&got), decoded(&got));

        std::fs::remove_dir_all(&dir).ok();
    }

    /// A cache hit has no HTTP response to take an ETag from, which is exactly why the
    /// digest is computed from the bytes: `--offline` runs are the ones whose inputs are
    /// least self-evident, and they used to record nothing at all.
    #[test]
    fn a_cache_hit_records_the_bytes_it_read() {
        let dir = tempdir("prov");
        let bytes = zstd::encode_all(&b">CACHED0003_0\nACGT\n"[..], 3).unwrap();
        std::fs::write(cache_path(&dir, "CACHED0003"), &bytes).unwrap();
        let f = Fetcher::new(&Config {
            unitigs_cache: Some(dir.clone()),
            offline: true,
            ..Config::default()
        });
        let p = f.fetch(0, &"CACHED0003".to_string()).expect("cache hit").payload;
        assert_eq!(p.source, PayloadSource::Cache);
        assert_eq!(p.bytes, Some(bytes.len() as u64));
        assert_eq!(p.sha256.as_deref(), Some(crate::types::sha256_hex(&bytes).as_str()));
        assert_eq!(p.etag, None, "there was no response to take one from");
        assert_eq!(
            p.fetched_at, None,
            "cached bytes may predate the run by any amount; run.started is the honest clock"
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A cached file over the cap is skipped, and the skip carries its size.
    #[test]
    fn a_cached_file_over_the_cap_records_its_size() {
        let dir = tempdir("capprov");
        let bytes = zstd::encode_all(&b">CACHED0004_0\nACGTACGT\n"[..], 3).unwrap();
        std::fs::write(cache_path(&dir, "CACHED0004"), &bytes).unwrap();
        let f = Fetcher::new(&Config {
            unitigs_cache: Some(dir.clone()),
            offline: true,
            max_compressed_bytes: bytes.len() as u64 - 1,
            ..Config::default()
        });
        match f.fetch(0, &"CACHED0004".to_string()) {
            Err(Error::TooLarge { payload, .. }) => {
                assert_eq!(payload.bytes, Some(bytes.len() as u64));
                assert_eq!(payload.source, PayloadSource::Cache);
            }
            other => panic!("expected TooLarge, got {:?}", other.map(|f| f.compressed_len)),
        }
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A cache miss under --offline fails fast, without touching the network. The error
    /// message names the accession and the expected cache path so a user can act on it.
    #[test]
    fn a_cache_miss_under_offline_errors_without_network() {
        let dir = tempdir("miss");
        let cfg = Config {
            unitigs_cache: Some(dir.clone()),
            offline: true,
            fetch_attempts: 1,
            ..Config::default()
        };
        let f = Fetcher::new(&cfg);
        let started = std::time::Instant::now();
        let err = f
            .fetch(0, &"NOTINCACHE".to_string())
            .err()
            .expect("must fail");
        assert!(started.elapsed() < Duration::from_secs(2), "must not hit the network");
        let msg = err.to_string();
        assert!(msg.contains("NOTINCACHE"), "error names accession: {msg}");
        assert!(msg.contains("cache miss under --offline"), "error names cause: {msg}");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// An empty cache file is rejected explicitly rather than being handed to zstd as an
    /// empty payload. A killed writer normally cannot produce this because writes are
    /// atomic (tempfile + rename), but a manually placed empty file must still surface.
    #[test]
    fn an_empty_cache_file_is_rejected_not_returned_as_success() {
        let dir = tempdir("empty");
        std::fs::write(cache_path(&dir, "EMPTY"), b"").unwrap();
        let cfg = Config {
            unitigs_cache: Some(dir.clone()),
            offline: true,
            ..Config::default()
        };
        let f = Fetcher::new(&cfg);
        let err = f.fetch(0, &"EMPTY".to_string()).err().expect("must fail");
        assert!(err.to_string().contains("empty"), "{err}");
        std::fs::remove_dir_all(&dir).ok();
    }

    fn tmpfiles_in(dir: &Path) -> Vec<String> {
        std::fs::read_dir(dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .filter(|n| n.ends_with(".tmp"))
            .collect()
    }

    /// The cache gains whole files or nothing. The rename is the commit, and it happens only
    /// in `finish`.
    ///
    /// This matters more than it did. The old code assembled the object in memory and wrote
    /// it out in one `fs::write`, so a partial cache entry was close to impossible; now the
    /// bytes land in the tempfile as they arrive, and the only thing standing between a
    /// half-downloaded 11 GB object and a cache entry that a later `--offline` run would
    /// trust is that the rename has not happened yet.
    #[test]
    fn a_committed_download_lands_whole_and_an_abandoned_one_lands_not_at_all() {
        let dir = tempdir("atomic");

        let mut sink = open_cache_tmp(&dir, "ATOMIC1").unwrap();
        sink.write_all(b"pay").unwrap();
        sink.write_all(b"load").unwrap();
        // Mid-download: the final path must not exist yet, and the bytes must be in a
        // tempfile rather than nowhere.
        assert!(!cache_path(&dir, "ATOMIC1").exists(), "renamed before the body finished");
        assert_eq!(tmpfiles_in(&dir).len(), 1);
        sink.finish().unwrap();
        assert_eq!(std::fs::read(cache_path(&dir, "ATOMIC1")).unwrap(), b"payload");
        assert!(tmpfiles_in(&dir).is_empty(), "leftover tempfiles: {:?}", tmpfiles_in(&dir));

        // An abandoned download leaves neither a cache entry nor a tempfile. Every refusal
        // and every transport failure in `try_once` goes through this path.
        let mut sink = open_cache_tmp(&dir, "ATOMIC2").unwrap();
        sink.write_all(b"half of an object").unwrap();
        sink.abandon();
        assert!(!cache_path(&dir, "ATOMIC2").exists(), "an abandoned download was cached");
        assert!(tmpfiles_in(&dir).is_empty(), "leftover tempfiles: {:?}", tmpfiles_in(&dir));

        std::fs::remove_dir_all(&dir).ok();
    }

    /// Where a download lands, decided before a byte of it arrives.
    ///
    /// The ordering is the whole of the memory fix: an object only reaches RAM when there is
    /// nowhere on disk to put it, and `main` always resolves a spool directory, so no CLI run
    /// can reach that case. Checked here rather than through a real download because the
    /// decision is made from configuration alone and the suite has no network.
    #[test]
    fn a_download_lands_in_the_cache_then_the_spool_then_memory() {
        let cache = tempdir("sink-cache");
        let spool = tempdir("sink-spool");
        let acc = "SINKACC001".to_string();

        // 1. A cache is where the bytes are wanted anyway, so the download writes it directly.
        let cfg = Config {
            unitigs_cache: Some(cache.clone()),
            spool_dir: Some(spool.clone()),
            ..Config::default()
        };
        let sink = Fetcher::new(&cfg).open_sink(&acc);
        assert!(matches!(sink, BodySink::Cache { .. }), "a cache must win");
        assert_eq!(sink.describe(), "cache");
        sink.abandon();

        // 2. No cache, but a spool directory: an unlinked file, so still never in RAM.
        let cfg = Config { spool_dir: Some(spool.clone()), ..Config::default() };
        let sink = Fetcher::new(&cfg).open_sink(&acc);
        assert!(matches!(sink, BodySink::Transient(_)), "a spool dir must be used");
        assert_eq!(sink.describe(), "spool");

        // 3. Neither: the old behaviour, and the only path that costs the object's size.
        let cfg = Config::default();
        assert!(!cfg.downloads_stream_to_disk());
        let sink = Fetcher::new(&cfg).open_sink(&acc);
        assert!(matches!(sink, BodySink::Memory(_)), "nowhere on disk, so RAM");

        // 4. An unusable cache directory falls through rather than failing the accession,
        //    the same way an unusable spool directory does in `scan_in`.
        let blocked = cache.join("not-a-dir");
        std::fs::write(&blocked, b"x").unwrap();
        let cfg = Config {
            unitigs_cache: Some(blocked.join("under-a-file")),
            spool_dir: Some(spool.clone()),
            ..Config::default()
        };
        let sink = Fetcher::new(&cfg).open_sink(&acc);
        assert!(matches!(sink, BodySink::Transient(_)), "an unusable cache must fall through");

        std::fs::remove_dir_all(&cache).ok();
        std::fs::remove_dir_all(&spool).ok();
    }

    /// Name a `FetchAttempt` variant for a panic message; it has no `Debug` and does not
    /// need one outside tests.
    fn kind_of(a: &FetchAttempt) -> &'static str {
        match a {
            FetchAttempt::Transient(_) => "Transient",
            FetchAttempt::Permanent(_) => "Permanent",
            FetchAttempt::Missing => "Missing",
            FetchAttempt::TooLarge(..) => "TooLarge",
        }
    }

    /// A download that fills the disk is the user's budget, not a broken run, and it is
    /// classified exactly as `parse::spill_error` classifies the same failure one stage later.
    #[test]
    fn a_download_that_fills_the_disk_is_a_negative_not_an_error() {
        let p = Payload::default();
        for raw in [libc::ENOSPC, libc::EDQUOT] {
            match download_write_error("cache", 99, std::io::Error::from_raw_os_error(raw), &p) {
                FetchAttempt::TooLarge(detail, _) => {
                    assert!(detail.contains("ran out of space"), "{detail}");
                    assert!(detail.contains("after 99 bytes"), "{detail}");
                    assert!(detail.contains("--unitigs-cache"), "{detail}");
                }
                other => panic!("errno {raw} must be a negative, got {}", kind_of(&other)),
            }
        }
        // Anything else is a real failure, and permanent: retrying a read-only directory
        // three times would just download the object three times to fail the same way.
        for raw in [libc::EACCES, libc::EIO] {
            match download_write_error("spool", 99, std::io::Error::from_raw_os_error(raw), &p) {
                FetchAttempt::Permanent(detail) => {
                    assert!(detail.contains("could not write"), "{detail}")
                }
                other => panic!("errno {raw} must be permanent, got {}", kind_of(&other)),
            }
        }
    }

    /// A download with no cache still never holds the object: it goes to an unlinked file
    /// under `--spool-dir`, which is readable, rewound and leaves nothing on the filesystem.
    #[test]
    fn a_transient_sink_is_readable_and_leaves_nothing_behind() {
        let dir = tempdir("transient");
        let before = std::fs::read_dir(&dir).unwrap().count();

        let mut sink = BodySink::Transient(crate::types::transient_file(&dir, "dl").unwrap());
        sink.write_all(b"compressed bytes").unwrap();
        let compressed = sink.finish().unwrap();

        // File-backed, so nothing may read it into a Vec by accident...
        assert!(compressed.as_slice().is_none());
        // ...but it is rewound and readable, which is what `decoder` needs.
        match &compressed {
            Compressed::File(f) => {
                let mut got = Vec::new();
                let mut h = f;
                h.read_to_end(&mut got).unwrap();
                assert_eq!(got, b"compressed bytes");
            }
            _ => panic!("a transient sink must be file-backed"),
        }
        // Unlinked at creation, so the directory never saw it.
        assert_eq!(std::fs::read_dir(&dir).unwrap().count(), before);
        drop(compressed);
        std::fs::remove_dir_all(&dir).ok();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fetched(compressed: Vec<u8>) -> FetchedUnitigs {
        FetchedUnitigs {
            index: 0,
            accession: "TEST0001".to_string(),
            compressed_len: compressed.len(),
            compressed: crate::types::Compressed::Memory(compressed),
            fetch_secs: 0.0,
            payload: crate::types::Payload::default(),
        }
    }

    fn read_all(f: &FetchedUnitigs) -> std::io::Result<Vec<u8>> {
        let mut out = Vec::new();
        decoder(f).expect("decoder").read_to_end(&mut out)?;
        Ok(out)
    }

    #[test]
    fn round_trips_a_unitig_record() {
        let plain = b">TEST0001_0 ka:f:5.1   L:+:1:-  L:-:2:-\nACGTACGTACGT\n".to_vec();
        let f = fetched(zstd::encode_all(&plain[..], 3).unwrap());
        assert_eq!(read_all(&f).unwrap(), plain);
    }

    /// A zstd file may hold several concatenated frames. Stopping after the first would
    /// silently truncate the graph, so this asserts we read through all of them.
    #[test]
    fn decodes_concatenated_frames() {
        let a = zstd::encode_all(&b">TEST0001_0\nAAAA\n"[..], 3).unwrap();
        let b = zstd::encode_all(&b">TEST0001_1\nCCCC\n"[..], 3).unwrap();
        let mut both = a.clone();
        both.extend_from_slice(&b);
        let f = fetched(both);
        assert_eq!(read_all(&f).unwrap(), b">TEST0001_0\nAAAA\n>TEST0001_1\nCCCC\n");
    }

    #[test]
    fn handles_large_payload_across_many_reads() {
        // Larger than one internal buffer, so this exercises repeated read() calls.
        let mut plain = Vec::new();
        for i in 0..50_000 {
            plain.extend_from_slice(format!(">TEST0001_{i} ka:f:3.0\nACGTACGTACGTACGTACGT\n").as_bytes());
        }
        let f = fetched(zstd::encode_all(&plain[..], 3).unwrap());
        assert_eq!(read_all(&f).unwrap(), plain);
    }

    #[test]
    fn empty_input_is_rejected_or_empty_not_silently_wrong() {
        let f = fetched(Vec::new());
        let produced = match decoder(&f) {
            Err(_) => Vec::new(),
            Ok(mut r) => {
                let mut out = Vec::new();
                let _ = r.read_to_end(&mut out);
                out
            }
        };
        assert!(produced.is_empty(), "empty input must not yield data");
    }

    #[test]
    fn garbage_is_an_error_not_data() {
        let f = fetched(b"this is definitely not a zstd frame".to_vec());
        let produced = decoder(&f).ok().and_then(|mut r| {
            let mut out = Vec::new();
            r.read_to_end(&mut out).ok().map(|_| out)
        });
        assert!(produced.is_none(), "garbage must not decode to data");
    }

    /// The decoder validates only the frame header, so corruption further in shows up as an
    /// io::Error during read, not at construction. parse::build must not mistake this for EOF.
    #[test]
    fn truncated_frame_errors_during_read_not_at_construction() {
        let plain: Vec<u8> = (0..40_000).flat_map(|i| format!(">T_{i}\nACGTACGT\n").into_bytes()).collect();
        let full = zstd::encode_all(&plain[..], 3).unwrap();
        let f = fetched(full[..full.len() / 2].to_vec());
        let mut r = decoder(&f).expect("header alone should still open");
        let mut out = Vec::new();
        assert!(r.read_to_end(&mut out).is_err(), "truncation must surface as an error");
    }
}

/// Integration checks against real Logan files. Ignored by default because they need data
/// on disk. Point `LOGAN_TEST_ZST` at a real `.unitigs.fa.zst`:
///
/// ```text
/// LOGAN_TEST_ZST=/path/to/ERR2028820.unitigs.fa.zst \
///     cargo test --release -- --ignored --nocapture
/// ```
///
/// The expected plaintext is produced by shelling out to the `zstd` CLI, so this compares
/// against an independent implementation rather than against ourselves. Override with
/// `LOGAN_TEST_PLAIN` if you already have the decompressed file.
#[cfg(test)]
mod real_data_tests {
    use super::*;

    /// Decompress with the zstd CLI, so the comparison is against a separate implementation.
    fn reference_plaintext(zst: &str) -> Option<Vec<u8>> {
        if let Ok(p) = std::env::var("LOGAN_TEST_PLAIN") {
            return std::fs::read(p).ok();
        }
        let out = std::process::Command::new("zstd")
            .arg("-dc")
            .arg(zst)
            .output()
            .ok()?;
        out.status.success().then_some(out.stdout)
    }

    #[test]
    #[ignore]
    fn matches_reference_zstd_on_a_real_logan_file() {
        let Ok(path) = std::env::var("LOGAN_TEST_ZST") else {
            eprintln!("SKIP: set LOGAN_TEST_ZST to a real .unitigs.fa.zst");
            return;
        };
        let compressed = match std::fs::read(&path) {
            Ok(c) => c,
            Err(e) => {
                eprintln!("SKIP: cannot read LOGAN_TEST_ZST={path}: {e}");
                return;
            }
        };
        let Some(expected) = reference_plaintext(&path) else {
            eprintln!("SKIP: could not produce reference plaintext (is `zstd` on PATH?)");
            return;
        };

        let f = FetchedUnitigs {
            index: 0,
            accession: "REAL".to_string(),
            compressed_len: compressed.len(),
            compressed: crate::types::Compressed::Memory(compressed),
            fetch_secs: 0.0,
            payload: crate::types::Payload::default(),
        };
        let mut out = Vec::new();
        decoder(&f)
            .expect("decoder")
            .read_to_end(&mut out)
            .expect("decode");

        eprintln!(
            "decoded {} bytes from {} compressed ({:.1}x) from {}",
            out.len(),
            f.compressed_len,
            out.len() as f64 / f.compressed_len as f64,
            path
        );
        assert_eq!(out.len(), expected.len(), "length mismatch vs zstd CLI");
        assert!(out == expected, "content mismatch vs zstd CLI");
    }
}
