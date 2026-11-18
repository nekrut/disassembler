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

/// The provenance sidecar beside a cached object. See `types::Sidecar`.
pub fn sidecar_path(object: &Path) -> PathBuf {
    let mut p = object.as_os_str().to_os_string();
    p.push(".provenance");
    PathBuf::from(p)
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
    /// `--revalidate`: confirm each cache hit against the bucket before using it.
    revalidate: bool,
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
            revalidate: cfg.revalidate,
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
    /// into place once the body is complete. A cached file's contents are trusted as far as
    /// this layer goes: the walker's own writer cannot leave a short file behind, because the
    /// rename happens only after the body is complete. It is not the only writer, though, and
    /// a directory filled by something else (a staged dataset, a partial copy) can hold one.
    /// The magic-number check catches an entry that is not zstd at all, and one cut at a
    /// frame boundary decodes cleanly and is caught in `parse` from the links. **The size is
    /// not trusted**: see
    /// `the_size_cap_decides_the_same_way_wherever_the_bytes_came_from`.
    ///
    /// Returns the compressed bytes only. Decompression is `decoder`'s job, and happens in
    /// the compute pool so the fetch thread can immediately start the next transfer.
    pub fn fetch(&self, index: usize, accession: &Accession) -> Result<FetchedUnitigs> {
        if let Some(dir) = &self.cache_dir {
            let path = cache_path(dir, accession);
            // What `--revalidate` decided about this entry, if it ran. `None` here means the
            // cached bytes are used as they always were.
            let mut revalidated: Option<&'static str> = None;
            if path.exists() && self.revalidate && !self.offline {
                match self.check_cached(index, accession, &path, &unitigs_url(accession)) {
                    Revalidation::Current => revalidated = Some("current"),
                    Revalidation::Unavailable => revalidated = Some("unavailable"),
                    // The bytes changed, and the conditional GET already brought the new
                    // ones down and rewrote the entry. This is the event the flag exists to
                    // surface: an earlier run of the same command saw different data.
                    Revalidation::Replaced(fresh) => {
                        let mut fresh = *fresh;
                        eprintln!(
                            "[logan-walker] [{accession}] REPLACED: the bucket now serves \
                             different bytes than the cache held; the cache has been \
                             refreshed and this run uses the new object."
                        );
                        fresh.payload.revalidated = Some("replaced");
                        return Ok(fresh);
                    }
                    // Nothing recorded what the cached bytes were a copy of, so there is
                    // nothing to revalidate against. Falling through to the ordinary
                    // download both answers the question and leaves a sidecar, so the cache
                    // heals on first use rather than staying unverifiable forever.
                    Revalidation::NoOrigin => {
                        eprintln!(
                            "[logan-walker] [{accession}] no recorded origin for the cached \
                             object, so it cannot be revalidated; downloading it again."
                        );
                        return self.download(index, accession, Some("refetched"));
                    }
                }
            }
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
                // "Not empty" is not enough: a proxy's error page cached under the
                // accession's name fails that accession with `parse failed: Unknown frame
                // descriptor` on every later run, online ones included, until someone finds
                // and deletes it. Four bytes say whether this is a zstd stream at all.
                if !entry_is_zstd(&mut file) {
                    drop(file);
                    self.discard_bad_entry(accession, &path)?;
                    return self.download_recording(index, accession, Some("refetched"));
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
                let mut payload =
                    Payload::read_file(PayloadSource::Cache, &mut file).map_err(|e| {
                        Error::Fetch {
                            accession: accession.clone(),
                            detail: format!("cache read {}: {e}", path.display()),
                        }
                    })?;
                // What the object was downloaded from, when a sidecar beside it says so.
                // Believed only when its digest matches the bytes just hashed: a sidecar that
                // disagrees is describing an object that has since been replaced underneath
                // it, which is exactly the substitution this provenance exists to expose.
                // An entry written before 0.2.0 has no sidecar and reads as it always did.
                if let Some(side) = read_sidecar(&path) {
                    if side.describes(payload.sha256.as_deref()) {
                        payload.adopt(&side);
                    } else {
                        eprintln!(
                            "[logan-walker] warn: [{accession}] {} describes different bytes \
                             than the file beside it; ignoring it. The cached object is \
                             intact, but its recorded origin is not trustworthy.",
                            sidecar_path(&path).display()
                        );
                    }
                }
                payload.revalidated = revalidated;
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
                // Not every miss is a failure. `disassembler fetch` writes nothing for
                // an accession Logan has no unitigs for, or one over the cap, because there
                // is nothing to write; the documented `--offline` run that follows then saw
                // an ordinary miss and called it an error. The same list and the same flags
                // gave `no_logan_data` and exit 0 online, `error` and exit 1 or 3 offline,
                // and since 0.2.0 any error makes the run non-zero, so a shard of 404s now
                // fails outright. The fetch leaves a one-line note saying which negative it
                // was, and this reads it back.
                if let Some(negative) =
                    read_absent(dir, accession).and_then(|a| a.into_error(accession, self.max_compressed_bytes))
                {
                    return Err(negative);
                }
                return Err(Error::Fetch {
                    accession: accession.clone(),
                    detail: format!(
                        "cache miss under --offline; expected {}",
                        path.display()
                    ),
                });
            }
        }

        self.download_recording(index, accession, None)
    }

    /// What to do with a cache entry that is not a zstd stream.
    ///
    /// Offline there is nothing to replace it with, so it is refused **by name**: the user
    /// has to be told which file to delete, because the walker will not delete a file it was
    /// told to treat as the only copy. Online it is removed, with its provenance sidecar,
    /// and the caller downloads the object again.
    fn discard_bad_entry(&self, accession: &Accession, path: &Path) -> Result<()> {
        if self.offline {
            return Err(Error::Fetch {
                accession: accession.clone(),
                detail: format!(
                    "the cached file {} is not zstd data, so it cannot be the Logan object. \
                     Delete it, or rerun without --offline to replace it",
                    path.display()
                ),
            });
        }
        eprintln!(
            "[logan-walker] warn: [{accession}] the cached file {} is not zstd data; \
             removing it and downloading the object again.",
            path.display()
        );
        let _ = std::fs::remove_file(path);
        let _ = std::fs::remove_file(sidecar_path(path));
        Ok(())
    }

    /// `download`, recording a negative in the cache so an `--offline` run can tell a 404
    /// from a failure.
    ///
    /// Only when there is a cache to record it in, and only for the two outcomes that are
    /// properties of the accession rather than of this run: no Logan data, and over the cap.
    /// A success removes any note, because Logan does add data for accessions it had none
    /// for, and because a cap raised between runs must not keep refusing.
    fn download_recording(
        &self,
        index: usize,
        accession: &Accession,
        mark: Option<&'static str>,
    ) -> Result<FetchedUnitigs> {
        let out = self.download(index, accession, mark);
        if let Some(dir) = &self.cache_dir {
            note_outcome(dir, accession, &out, self.max_compressed_bytes);
        }
        out
    }

    /// The cold path: fetch the object, with retries, and mark how it came to be fetched.
    ///
    /// `mark` lands on `payload.revalidated`. It is `Some("refetched")` only when
    /// `--revalidate` sent us here because the cache held an object with no recorded origin,
    /// and `None` for an ordinary download.
    fn download(
        &self,
        index: usize,
        accession: &Accession,
        mark: Option<&'static str>,
    ) -> Result<FetchedUnitigs> {
        let url = unitigs_url(accession);
        let mut last: String = String::from("no attempt made");
        // One transfer across all the attempts, so an attempt that was interrupted after
        // real progress resumes where it stopped rather than starting again.
        let mut xfer = Transfer::new(self.open_sink(accession));
        // Failures that made no progress. `attempt` counts every try, including the free
        // ones, and only guards against a server that dribbles forever.
        let mut failures: u32 = 0;

        for _attempt in 1..=MAX_FETCH_ATTEMPTS_EVER {
            let before = xfer.got;
            match self.try_once(index, accession, &url, None, &mut xfer) {
                // No cache write here any more. `try_once` streams the body straight into the
                // cache entry, so by the time it returns the bytes are already on disk and
                // already renamed into place. The old shape read the object into a `Vec`,
                // wrote that `Vec` out, and then reopened the file to drop the `Vec`: three
                // copies of the same decision and one of them the size of the accession.
                Ok(mut f) => {
                    f.payload.revalidated = mark;
                    return Ok(f);
                }
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
                // Unreachable: this loop never sends `If-None-Match`, so the server has
                // nothing to answer 304 against. Treated as transport noise rather than
                // panicking, because a proxy inventing one must not take down a run.
                Err(FetchAttempt::NotModified) => {
                    last = "unconditional request answered 304".to_string();
                    failures += 1;
                }
                Err(FetchAttempt::TooLarge(detail, payload)) => {
                    return Err(Error::TooLarge {
                        accession: accession.clone(),
                        detail,
                        payload,
                    })
                }
                Err(FetchAttempt::Transient(detail, after)) => {
                    let progressed = xfer.got.saturating_sub(before) >= RESUME_PROGRESS_BYTES;
                    last = detail;
                    if !progressed {
                        failures += 1;
                    }
                    if failures >= self.attempts {
                        break;
                    }
                    // A server that asked for a particular wait gets it; otherwise back off.
                    // An interrupted transfer that is making progress retries at once, since
                    // the delay is for a server in trouble, not for a dropped connection.
                    if let Some(wait) = after {
                        std::thread::sleep(wait);
                    } else if !progressed {
                        std::thread::sleep(backoff(failures));
                    }
                    continue;
                }
            }
            if failures >= self.attempts {
                break;
            }
        }
        Err(Error::Fetch {
            accession: accession.clone(),
            detail: format!(
                "giving up after {failures} failed attempts ({} bytes transferred): {last}",
                xfer.got
            ),
        })
    }

    /// Ask the bucket whether the cached object is still what it serves.
    ///
    /// A conditional GET, so the three answers arrive in one round trip: 304 (still current,
    /// no body), 200 (replaced, and the new bytes come with it and have already been written
    /// into the cache by `try_once`), or a failure. A failure is deliberately not fatal: the
    /// cached bytes are still a usable analysis, and refusing to run because a revalidation
    /// could not be made would turn a network blip into a lost cohort. It is recorded as
    /// `unavailable` so the summary does not claim a check that did not happen.
    fn check_cached(
        &self,
        index: usize,
        accession: &Accession,
        path: &Path,
        url: &str,
    ) -> Revalidation {
        let Some(tag) = read_sidecar(path).and_then(|s| s.etag) else {
            return Revalidation::NoOrigin;
        };
        // A 200 here carries the replacement object, so this needs somewhere to put it. One
        // attempt only: a revalidation that cannot be completed falls back to the cached
        // bytes, which is the whole point of `Unavailable`.
        let mut xfer = Transfer::new(self.open_sink(accession));
        match self.try_once(index, accession, url, Some(&tag), &mut xfer) {
            Err(FetchAttempt::NotModified) => Revalidation::Current,
            Ok(fresh) => Revalidation::Replaced(Box::new(fresh)),
            // The object is gone from the bucket. The cached copy is still readable and is
            // still the best answer available, so the run continues on it and says so.
            Err(FetchAttempt::Missing) => {
                eprintln!(
                    "[logan-walker] warn: [{accession}] the bucket no longer holds this \
                     accession; using the cached copy."
                );
                Revalidation::Unavailable
            }
            Err(FetchAttempt::TooLarge(detail, _)) => {
                eprintln!(
                    "[logan-walker] warn: [{accession}] could not revalidate ({detail}); \
                     using the cached copy."
                );
                Revalidation::Unavailable
            }
            Err(FetchAttempt::Permanent(detail)) | Err(FetchAttempt::Transient(detail, _)) => {
                eprintln!(
                    "[logan-walker] warn: [{accession}] could not revalidate ({detail}); \
                     using the cached copy."
                );
                Revalidation::Unavailable
            }
        }
    }

    /// One attempt at the object. `if_none_match` makes it conditional: the server answers
    /// 304 with no body when the cached `ETag` is still the current one, which is what makes
    /// `--revalidate` cost a round trip rather than a download.
    fn try_once(
        &self,
        index: usize,
        accession: &Accession,
        url: &str,
        if_none_match: Option<&str>,
        xfer: &mut Transfer,
    ) -> std::result::Result<FetchedUnitigs, FetchAttempt> {
        let req = self.agent.get(url);
        // Quoted on the way out: the stored value has them stripped for comparison, but the
        // header is an entity-tag and S3 will not match a bare one.
        let req = match if_none_match {
            Some(tag) => req.header("If-None-Match", &format!("\"{tag}\"")),
            None => req,
        };
        // Resume. `If-Match` is what makes it safe: without it a resume could staple
        // the tail of a replaced object onto the head of the one we started with, and the
        // Logan objects are mutable (3.1% of a 5,000-key sample were rewritten after the
        // 2024-03 load). A server that does not know the entity answers 412 and we start over.
        let resuming = xfer.got > 0;
        let req = if resuming {
            let r = req.header("Range", &format!("bytes={}-", xfer.got));
            match &xfer.etag {
                Some(tag) => r.header("If-Match", &format!("\"{tag}\"")),
                None => r,
            }
        } else {
            req
        };
        let mut resp = req
            .call()
            .map_err(|e| FetchAttempt::Transient(format!("request failed: {e}"), None))?;

        let status = resp.status().as_u16();
        // How much of the object this response carries, and from where.
        let mut range_total: Option<u64> = None;
        match status {
            200 => {
                // The whole object, whatever we asked for. Anything already transferred is
                // the head of a different response and cannot be kept.
                if resuming {
                    xfer.restart(self.open_sink(accession));
                }
            }
            206 => {
                let start = header(&resp, ureq::http::header::CONTENT_RANGE)
                    .as_deref()
                    .and_then(parse_content_range);
                match start {
                    Some((first, total)) if resuming && first == xfer.got => range_total = total,
                    // A partial response to a request we did not make one for, or one that
                    // starts somewhere else. Neither can be appended to what we hold.
                    other => {
                        xfer.restart(self.open_sink(accession));
                        return Err(FetchAttempt::Transient(
                            format!("unusable 206: Content-Range {other:?} for a resume at {}", xfer.got),
                            None,
                        ));
                    }
                }
            }
            // Only reachable when `if_none_match` was sent, and it is the answer that run
            // wanted: the cached bytes are the object the bucket still serves.
            304 => return Err(FetchAttempt::NotModified),
            404 => return Err(FetchAttempt::Missing),
            // The object changed under us (412), or the range is no longer inside it (416).
            // Both mean the bytes in hand are stale; drop them and fetch the object afresh.
            412 | 416 => {
                xfer.restart(self.open_sink(accession));
                return Err(FetchAttempt::Transient(
                    format!("HTTP {status}: the object changed while it was being fetched"),
                    None,
                ));
            }
            // Request timeout and "slow down". The server is asking to be asked again, and
            // may have said when; before 0.2.0 both were permanent on the first try.
            408 | 429 => {
                return Err(FetchAttempt::Transient(
                    format!("HTTP {status}"),
                    retry_after(&resp),
                ))
            }
            // 4xx other than those is our fault and will not improve on retry.
            400..=499 => return Err(FetchAttempt::Permanent(format!("HTTP {status}"))),
            _ => return Err(FetchAttempt::Transient(format!("HTTP {status}"), None)),
        }

        // Provenance off the response, read before the body so it survives a refusal: a
        // `too_large` skip can still say which object was too big, and how big.
        //
        // The ETag is the object's MD5 for everything below S3's 5 GB single-PUT limit, which
        // `--max-size-mb` (default 1024) keeps the walker well inside. Above that limit S3
        // appends `-<parts>` and the value is opaque; it is recorded either way, and the
        // suffix is self-identifying.
        // Captured on the response that opened the object and kept across resumes, so a
        // transfer that took four connections still reports the identity it was started on.
        if xfer.payload.etag.is_none() && xfer.payload.last_modified.is_none() {
            xfer.payload = Payload {
                source: PayloadSource::S3,
                etag: header(&resp, ureq::http::header::ETAG)
                    .map(|v| v.trim_matches('"').to_string()),
                last_modified: header(&resp, ureq::http::header::LAST_MODIFIED)
                    .as_deref()
                    .and_then(http_date_to_iso8601),
                fetched_at: Some(crate::io::output::iso8601(std::time::SystemTime::now())),
                ..Payload::default()
            };
            xfer.etag = xfer.payload.etag.clone();
        }
        let mut payload = xfer.payload.clone();

        // Refuse oversized accessions before transferring them, when the server says how big.
        // On a 206 the content length is what is left, so the object's size is the total from
        // `Content-Range` when the server gave one, and what we hold plus what is coming
        // otherwise.
        let object_len = resp
            .body()
            .content_length()
            .map(|len| range_total.unwrap_or(xfer.got + len));
        if let Some(len) = object_len {
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

        // The sink was opened before the first attempt and survives them, because that is
        // what a resume appends to. Only a restart replaces it.
        let sink = xfer.sink.as_mut().expect("the transfer's sink is taken only by finish");

        // One pass over the body: counted, hashed and written as it streams. The count is
        // ours rather than ureq's `limit`, so the cap is enforced by code that can say which
        // limit was reached and stop the transfer at exactly that point. The hash is free
        // here, where the bytes are already in a register; computing it afterwards would be
        // a second full pass over a file that can be 11 GB (the nn1 probe measured that pass
        // at 91 s).
        use sha2::Digest;
        let mut buf = vec![0u8; crate::config::DOWNLOAD_BUFFER_BYTES as usize];
        // The object's first four bytes, collected across however many reads it takes, so
        // the zstd check below is about the object and not about one buffer. Only
        // for a transfer that starts at byte 0; a resume already passed it.
        let mut head: Vec<u8> = Vec::with_capacity(4);
        let mut checked_magic = xfer.got > 0;
        let mut reader = resp.body_mut().as_reader();
        loop {
            let got = match reader.read(&mut buf) {
                Ok(0) => break,
                Ok(n) => n,
                // A dropped connection mid-body. Worth retrying, and what has been written
                // is kept: the next attempt asks for the rest with a `Range` header.
                Err(e) => {
                    return Err(FetchAttempt::Transient(
                        format!("body read failed after {} bytes: {e}", xfer.got),
                        None,
                    ))
                }
            };
            xfer.got += got as u64;
            // The guard for a missing or lying `Content-Length`. A body of exactly the cap
            // passes, as it does the pre-transfer check above; the transfer stops here rather
            // than running to completion first.
            if xfer.got > self.max_compressed_bytes {
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
            // Refuse a body that is not a zstd stream, before it can be renamed into
            // the cache. A proxy's error page used to be cached under the accession's name
            // and then fail that accession on every run, forever, even online.
            if !checked_magic {
                head.extend_from_slice(&buf[..got.min(4 - head.len())]);
                if head.len() == 4 {
                    checked_magic = true;
                    if !looks_like_zstd(&head) {
                        return Err(FetchAttempt::Transient(
                            format!(
                                "the body is not zstd ({:02x?}); a proxy or error page, not \
                                 a Logan object",
                                &head[..]
                            ),
                            None,
                        ));
                    }
                }
            }
            xfer.hasher.update(&buf[..got]);
            if let Err(e) = xfer.sink.as_mut().expect("sink").write_all(&buf[..got]) {
                let where_ = xfer.sink.as_ref().expect("sink").describe();
                return Err(download_write_error(where_, xfer.got, e, &payload));
            }
        }

        if xfer.got == 0 {
            return Err(FetchAttempt::Transient("empty body".to_string(), None));
        }
        // A short read that the connection closed cleanly on: the server said how big the
        // object is and sent less. Resumable like any other interruption.
        if let Some(len) = object_len {
            if xfer.got < len {
                return Err(FetchAttempt::Transient(
                    format!("body ended at {} of {len} bytes", xfer.got),
                    None,
                ));
            }
        }
        // A one-byte object cannot be zstd, and the loop above only checks once it has four.
        if !checked_magic && !looks_like_zstd(&head) {
            return Err(FetchAttempt::Transient(
                format!("the body is {} bytes and is not zstd", xfer.got),
                None,
            ));
        }

        // Content-Length is the server's claim; this is the bytes that were actually read,
        // and it is what every downstream stage consumed.
        let total = xfer.got;
        payload.bytes = Some(total);

        // The spool gate again, for a server that sent no `Content-Length`: then the size is
        // only knowable by downloading, and refusing here still beats failing mid-parse with
        // the transfer already spent.
        if let Err(detail) = self.spool_room(total) {
            return Err(FetchAttempt::TooLarge(detail, Box::new(payload)));
        }

        let mut hex = String::with_capacity(64);
        let hasher = std::mem::replace(&mut xfer.hasher, <sha2::Sha256 as sha2::Digest>::new());
        for b in hasher.finalize() {
            hex.push_str(&format!("{b:02x}"));
        }
        payload.sha256 = Some(hex);

        // The rename into the cache happens here, after the last byte and after every
        // refusal, so a skipped accession never leaves a cache entry behind. Taking the sink
        // is what marks the transfer finished: `Drop` abandons anything still in it.
        let sink = xfer.sink.take().expect("sink");
        let compressed = sink.finish(&payload).map_err(|e| {
            FetchAttempt::Transient(format!("could not commit the downloaded body: {e}"), None)
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
    /// 3. RAM, the fallback when neither is available.
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

/// Whether a cache entry opens with the zstd magic.
///
/// Reads four bytes and rewinds, so the caller's handle is still at the start for the hash
/// pass that follows. An unreadable handle is treated as "not zstd": the read that follows
/// would fail anyway, and this way it fails with a message about the file.
fn entry_is_zstd(file: &mut std::fs::File) -> bool {
    use std::io::{Read, Seek, SeekFrom};
    let mut head = [0u8; 4];
    let ok = file.read_exact(&mut head).is_ok() && looks_like_zstd(&head);
    file.seek(SeekFrom::Start(0)).is_ok() && ok
}

/// The note recording why an accession is absent from the cache.
///
/// `<DIR>/<ACC>.unitigs.absent`, one line. Deliberately not the provenance sidecar's format:
/// that one describes an object that exists, and this one exists because no object does.
/// Deleting it forces the next online run to ask the bucket again.
pub fn absent_path(dir: &Path, accession: &str) -> PathBuf {
    dir.join(format!("{accession}.unitigs.absent"))
}

/// What a `.absent` note says.
pub enum Absent {
    /// Logan has no unitigs for this accession.
    NoLoganData,
    /// The object was over the cap that was in force when it was recorded, at this size.
    TooLarge(u64),
}

impl Absent {
    /// The error this note stands for, or `None` when it no longer applies.
    ///
    /// A `too_large` note stops applying when the run's cap is at least the recorded size,
    /// which is what happens when the user raises `--max-size-mb` between the fetch and the
    /// walk. The caller then falls through to its ordinary cache-miss error, which says the
    /// object is not here, because it is not.
    pub fn into_error(self, accession: &str, cap: u64) -> Option<Error> {
        match self {
            Absent::NoLoganData => Some(Error::NoLoganData { accession: accession.to_string() }),
            Absent::TooLarge(bytes) if bytes > cap => Some(Error::TooLarge {
                accession: accession.to_string(),
                detail: format!(
                    "compressed size {bytes} exceeds max_compressed_bytes {cap}; recorded \
                     when the cache was filled, and not fetched"
                ),
                payload: Box::new(Payload { bytes: Some(bytes), ..Payload::default() }),
            }),
            Absent::TooLarge(_) => None,
        }
    }
}

/// Read a `.absent` note. Anything unparseable is no note at all, so a damaged one degrades
/// to the ordinary cache-miss error rather than inventing an outcome.
fn read_absent(dir: &Path, accession: &str) -> Option<Absent> {
    let text = std::fs::read_to_string(absent_path(dir, accession)).ok()?;
    let line = text.lines().next()?.trim();
    if line == "no_logan_data" {
        return Some(Absent::NoLoganData);
    }
    let bytes = line.strip_prefix("too_large ")?.trim().parse().ok()?;
    Some(Absent::TooLarge(bytes))
}

/// Update an accession's `.absent` note from what the fetch returned.
///
/// Only the two outcomes that are properties of the accession rather than of this run are
/// recorded: no Logan data, and over the cap. A success removes any note, because Logan does
/// add data for accessions it had none for, and because a cap raised between runs must not
/// keep refusing. Every other failure leaves the note alone: a transport error says nothing
/// about the accession.
fn note_outcome(dir: &Path, accession: &str, out: &Result<FetchedUnitigs>, cap: u64) {
    match out {
        Ok(_) => {
            let _ = std::fs::remove_file(absent_path(dir, accession));
        }
        Err(Error::NoLoganData { .. }) => write_absent(dir, accession, "no_logan_data"),
        Err(Error::TooLarge { payload, .. }) => {
            // The size is a lower bound when the server sent no Content-Length: all that was
            // learned is "more than the cap", so record the cap plus one, which is the
            // weakest true statement and still refuses at this cap.
            let bytes = payload.bytes.unwrap_or(cap + 1);
            write_absent(dir, accession, &format!("too_large {bytes}"));
        }
        Err(_) => {}
    }
}

/// Write a `.absent` note, through a tempfile so a reader never sees half of one.
///
/// Failure is ignored: the note is an optimisation of the offline run's diagnosis, and a
/// cache directory that cannot take it has already produced a louder error elsewhere.
fn write_absent(dir: &Path, accession: &str, line: &str) {
    let ns = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let tmp = dir.join(format!(".{accession}.absent.{}.{ns}.tmp", std::process::id()));
    if std::fs::write(&tmp, format!("{line}\n")).is_ok()
        && std::fs::rename(&tmp, absent_path(dir, accession)).is_err()
    {
        let _ = std::fs::remove_file(&tmp);
    }
}

/// Bytes an attempt must add for the resume to cost nothing against `--fetch-attempts`.
///
/// A connection that dropped after transferring real data has not failed in the way
/// the attempt budget is counting; it has made progress and been interrupted. Without this, a
/// 4 GB object on a link that drops every few minutes exhausts six attempts having downloaded
/// most of it six times.
const RESUME_PROGRESS_BYTES: u64 = 8 * 1024 * 1024;

/// Attempts at one object, however much progress each makes.
///
/// The free-resume rule above has to stop somewhere, or a server dribbling exactly
/// `RESUME_PROGRESS_BYTES` and hanging up would loop until the size cap. At the 4 GB default
/// cap a transfer that keeps making progress needs at most 512 resumes, so this cannot refuse
/// a download that was going to finish.
const MAX_FETCH_ATTEMPTS_EVER: u32 = 1000;

/// Whether `b` opens a zstd stream.
///
/// Accepting a cache entry on "not empty", and caching any 200 body, turns a proxy's HTML
/// error page into an entry that fails its accession with `parse failed: Unknown frame
/// descriptor` on every later run, online ones included. Four bytes settle it: `28 b5 2f fd`
/// is the frame magic, and `184d2a5?` is a skippable frame, which a legitimate file may
/// begin with.
fn looks_like_zstd(b: &[u8]) -> bool {
    if b.len() < 4 {
        return false;
    }
    let magic = u32::from_le_bytes([b[0], b[1], b[2], b[3]]);
    magic == 0xFD2F_B528 || (0x184D_2A50..=0x184D_2A5F).contains(&magic)
}

/// `Retry-After`, in seconds, capped at five minutes.
///
/// Only the delta-seconds form is read. The HTTP-date form would need epoch arithmetic
/// against a clock we do not otherwise trust, and falling back to the backoff is a safe
/// answer for it: the backoff is never longer than the date would have asked for by more than
/// a minute. S3 sends delta-seconds.
fn retry_after<B>(resp: &ureq::http::Response<B>) -> Option<Duration> {
    header(resp, ureq::http::header::RETRY_AFTER)?
        .trim()
        .parse::<u64>()
        .ok()
        .map(|secs| Duration::from_secs(secs.min(300)))
}

/// How long to wait before attempt `n` (1-based), with jitter.
///
/// This was `500 ms * attempt` over three attempts: 1.5 s in total, no jitter. A 60 s
/// network, proxy or DNS interruption therefore discarded one accession per worker every
/// 1.5 s, which at the default eight workers marks about 320 accessions `error`, and since
/// 0.2.0 any error makes the run non-zero. Exponential from 1 s, capped at a minute, so the
/// default six attempts cover about half a minute of outage and `--fetch-attempts 10` covers
/// four.
fn backoff(attempt: u32) -> Duration {
    let secs = 1u64 << attempt.clamp(1, 6).saturating_sub(1);
    let millis = (secs.min(60) * 1000) as f64 * jitter_fraction();
    Duration::from_millis(millis as u64)
}

/// A number in `0.5..1.5`, to spread retries that would otherwise arrive together.
///
/// Every fetch worker tends to fail at the same instant, because they fail for the same
/// reason, and a fixed backoff then sends them all back at the same instant. The clock alone
/// is not enough to separate them (they are inside the same microsecond), so a per-call
/// counter is mixed in and the pair is run through xorshift64*. No `rand` dependency: the
/// walker's dependency list is deliberately short, and this needs spread, not randomness.
fn jitter_fraction() -> f64 {
    use std::sync::atomic::{AtomicU64, Ordering};
    static SEQ: AtomicU64 = AtomicU64::new(0x9E37_79B9_7F4A_7C15);
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.subsec_nanos() as u64)
        .unwrap_or(0);
    let mut x = SEQ.fetch_add(0x9E37_79B9_7F4A_7C15, Ordering::Relaxed) ^ nanos;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    x = x.wrapping_mul(0x2545_F491_4F6C_DD1D);
    0.5 + (x >> 11) as f64 / (1u64 << 53) as f64
}

/// The first byte and the total size out of a `Content-Range: bytes 100-999/1000`.
///
/// The total is `None` for the `*` form, which says the server does not know it.
fn parse_content_range(v: &str) -> Option<(u64, Option<u64>)> {
    let rest = v.trim().strip_prefix("bytes ")?;
    let (range, total) = rest.split_once('/')?;
    let (start, _end) = range.split_once('-')?;
    Some((start.trim().parse().ok()?, total.trim().parse().ok()))
}

/// One object's download, carried across attempts so a dropped connection resumes.
///
/// Before this, a failed attempt threw away everything it had transferred and the next one
/// issued a bare GET from byte 0. With the default 900 s request timeout that meant any
/// object whose per-stream rate fell below size/900s could never be fetched at all, however
/// many attempts it was given: each one timed out at the same place.
struct Transfer {
    /// `None` only while a reset is swapping it, and after `finish` has taken it.
    sink: Option<BodySink>,
    hasher: sha2::Sha256,
    /// Bytes committed to the sink so far, which is where a resume asks to continue from.
    got: u64,
    /// The object's `ETag` from the first response, sent back as `If-Match` so a resume
    /// cannot staple the tail of a replaced object onto the head of the old one.
    etag: Option<String>,
    /// Provenance from the first response that carried a body.
    payload: Payload,
}

impl Transfer {
    fn new(sink: BodySink) -> Self {
        Self {
            sink: Some(sink),
            hasher: <sha2::Sha256 as sha2::Digest>::new(),
            got: 0,
            etag: None,
            payload: Payload::default(),
        }
    }

    /// Throw away what has been transferred and start again on `fresh`.
    ///
    /// Used when the server answers a resume with the whole object, with a range we did not
    /// ask for, or with "that object has changed".
    fn restart(&mut self, fresh: BodySink) {
        if let Some(old) = self.sink.replace(fresh) {
            old.abandon();
        }
        self.hasher = <sha2::Sha256 as sha2::Digest>::new();
        self.got = 0;
    }
}

impl Drop for Transfer {
    fn drop(&mut self) {
        // A download that never completed must not leave a tempfile behind, and must never
        // leave one renamed into the cache: `abandon` is the only way a sink ends unfinished.
        if let Some(sink) = self.sink.take() {
            sink.abandon();
        }
    }
}

/// Whether a failed attempt is worth repeating.
enum FetchAttempt {
    /// Transport error, 5xx, 408 or 429. Retry, after the server's `Retry-After` when it sent
    /// one and the backoff otherwise.
    Transient(String, Option<Duration>),
    /// A 4xx other than 404. Retrying cannot help.
    Permanent(String),
    /// HTTP 404: Logan has no unitigs for this accession. Not an error of the run.
    Missing,
    /// Over `--max-size-mb`, by Content-Length or, without one, by the bytes read. Skipped by
    /// the user's cap, not an error. Carries what the response said about the object, so the
    /// skip is recorded with a size and an identity rather than only a sentence.
    TooLarge(String, Box<Payload>),
    /// HTTP 304 to a conditional request: the cached object is still the current one.
    ///
    /// Not a failure. It is only ever produced when `--revalidate` sent an `If-None-Match`,
    /// and it is the answer that run was hoping for.
    NotModified,
}

/// What `--revalidate` learned about one cached object.
enum Revalidation {
    /// HTTP 304: the cache holds what the bucket still serves.
    Current,
    /// HTTP 200: different bytes now, already downloaded and written over the cache entry.
    /// Boxed because `FetchedUnitigs` is much larger than the other variants.
    Replaced(Box<FetchedUnitigs>),
    /// Nothing records what the cached bytes were a copy of, so there is nothing to ask about.
    NoOrigin,
    /// The question could not be answered. The cached bytes are used.
    Unavailable,
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
    fn finish(self, payload: &Payload) -> std::io::Result<Compressed> {
        match self {
            BodySink::Cache { mut file, tmp, final_path } => {
                file.flush_and_rewind()?;
                if let Err(e) = std::fs::rename(&tmp, &final_path) {
                    let _ = std::fs::remove_file(&tmp);
                    return Err(e);
                }
                // After the object, never before. A crash between the two leaves an object
                // with no sidecar, which is what every cache entry written before 0.2.0 looks
                // like and is handled; the other order would leave a sidecar describing a
                // file that is not there. A sidecar that cannot be written is not worth
                // failing a download over: it costs provenance, not bytes.
                let side = sidecar_path(&final_path);
                if let Err(e) = write_sidecar(&side, payload) {
                    eprintln!(
                        "[logan-walker] warn: cached {} but could not write {}: {e}. \
                         The object is usable; --revalidate will refetch it.",
                        final_path.display(),
                        side.display()
                    );
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
/// The sidecar beside a cached object, when there is a legible one.
///
/// Absent is the normal answer for any entry written before 0.2.0, and unreadable is treated
/// the same way: provenance is worth reporting, never worth failing an accession over.
fn read_sidecar(object: &Path) -> Option<crate::types::Sidecar> {
    std::fs::read_to_string(sidecar_path(object))
        .ok()
        .map(|t| crate::types::Sidecar::from_text(&t))
}

/// Write a sidecar through a tempfile, so a reader never sees a half-written one.
fn write_sidecar(path: &Path, payload: &Payload) -> std::io::Result<()> {
    use std::io::Write;
    let dir = path.parent().unwrap_or(Path::new("."));
    let ns = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let tmp = dir.join(format!(".provenance.{}.{ns}.tmp", std::process::id()));
    let mut f = std::fs::File::create(&tmp)?;
    let r = f
        .write_all(payload.sidecar().to_text().as_bytes())
        .and_then(|()| f.sync_all())
        .and_then(|()| std::fs::rename(&tmp, path));
    if r.is_err() {
        let _ = std::fs::remove_file(&tmp);
    }
    r
}

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
/// **One truncation does not surface here at all**. A file cut exactly at a frame
/// boundary is a complete, valid zstd stream of the frames it still has, so every read
/// succeeds and the decode ends cleanly on what is half an object. Nothing in this layer can
/// tell that apart from a small accession. `parse` catches it instead, from the links: unitig
/// ids are dense and sequential, so a link naming a record the file does not hold means
/// records are missing.
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

    /// One fresh attempt at `url`, with the `Transfer` that `download` would own.
    fn attempt(
        f: &Fetcher,
        acc: &str,
        url: &str,
    ) -> std::result::Result<FetchedUnitigs, FetchAttempt> {
        let acc = acc.to_string();
        let mut xfer = Transfer::new(f.open_sink(&acc));
        f.try_once(0, &acc, url, None, &mut xfer)
    }

    /// A body that is a plausible Logan object: zstd frame magic, then filler.
    fn zbody(n: usize) -> Vec<u8> {
        let mut v = vec![0x28, 0xb5, 0x2f, 0xfd];
        v.resize(n.max(4), b'Z');
        v
    }

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
    /// `Fetcher::new` leaves ureq's `Proxy::try_from_env()` in place, which is correct
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
        match attempt(&fetcher(10), "X", &url) {
            Err(FetchAttempt::TooLarge(d, _)) => assert!(d.contains("no Content-Length"), "{d}"),
            Err(FetchAttempt::Transient(d, _)) => panic!("retried as transient: {d}"),
            Err(_) => panic!("wrong failure"),
            Ok(f) => panic!("accepted {} bytes over a 10-byte cap", f.compressed_len),
        }
    }

    /// A download into the cache leaves a sidecar, and a later warm run reads it back.
    ///
    /// The `ETag` is on the wire exactly once. Without writing it down beside the object, a
    /// warm run can never say which bucket object its bytes are a copy of, and `--revalidate`
    /// has nothing to send an `If-None-Match` for.
    #[test]
    fn a_download_into_the_cache_records_where_the_bytes_came_from() {
        let ns = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        let dir = std::env::temp_dir()
            .join(format!("logan-walker-sidecar-w-{}-{ns}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let body = zstd::encode_all(&b">CACHED0011_0\nACGT\n"[..], 3).unwrap();
        let url = serve_once(
            "ETag: \"d2dd288b238c3c88b46acdcd78860dce\"\r\n\
             Last-Modified: Sun, 31 Mar 2024 17:38:36 GMT\r\n",
            body.clone(),
        );
        let mut f = fetcher(u64::MAX);
        f.cache_dir = Some(dir.clone());
        let got = match attempt(&f, "CACHED0011", &url) {
            Ok(got) => got,
            Err(FetchAttempt::Missing) => panic!("download failed: missing"),
            Err(FetchAttempt::TooLarge(d, _)) => panic!("download failed, too large: {d}"),
            Err(FetchAttempt::Permanent(d)) => panic!("download failed, permanent: {d}"),
            Err(FetchAttempt::Transient(d, _)) => panic!("download failed, transient: {d}"),
            Err(FetchAttempt::NotModified) => panic!("304 to an unconditional request"),
        };
        assert_eq!(got.payload.etag.as_deref(), Some("d2dd288b238c3c88b46acdcd78860dce"));

        let obj = cache_path(&dir, "CACHED0011");
        let side = read_sidecar(&obj).expect("a cached object must get a sidecar");
        assert!(
            side.describes(Some(&crate::types::sha256_hex(&body))),
            "the sidecar must name the digest of the object it sits beside"
        );
        assert_eq!(side.etag.as_deref(), Some("d2dd288b238c3c88b46acdcd78860dce"));
        assert_eq!(side.last_modified.as_deref(), Some("2024-03-31T17:38:36Z"));
        assert!(side.fetched_at.is_some(), "and when it was fetched");

        // The run that matters: offline, no response to learn anything from, and it still
        // knows what these bytes are a copy of.
        let warm = Fetcher::new(&Config {
            unitigs_cache: Some(dir.clone()),
            offline: true,
            ..Config::default()
        });
        let p = warm.fetch(0, &"CACHED0011".to_string()).expect("cache hit").payload;
        assert_eq!(p.source, PayloadSource::Cache);
        assert_eq!(p.etag.as_deref(), Some("d2dd288b238c3c88b46acdcd78860dce"));
        assert_eq!(p.last_modified.as_deref(), Some("2024-03-31T17:38:36Z"));
        std::fs::remove_dir_all(&dir).ok();
    }

    /// Serve one request with an arbitrary status line, for the 304 a conditional GET gets.
    fn serve_status(status: &str, extra: &str, body: Vec<u8>) -> String {
        let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
        let port = listener.local_addr().unwrap().port();
        let head = format!(
            "HTTP/1.1 {status}\r\nContent-Length: {}\r\n{extra}Connection: close\r\n\r\n",
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
        format!("http://127.0.0.1:{port}/x")
    }

    /// A cache dir holding one object, with or without a sidecar naming `etag`.
    fn cached_with(tag: &str, etag: Option<&str>, body: &[u8]) -> (PathBuf, Fetcher) {
        let ns = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        let dir = std::env::temp_dir()
            .join(format!("logan-walker-reval-{tag}-{}-{ns}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let obj = cache_path(&dir, "REVAL1");
        std::fs::write(&obj, body).unwrap();
        if let Some(e) = etag {
            let side = crate::types::Sidecar {
                sha256: Some(crate::types::sha256_hex(body)),
                etag: Some(e.to_string()),
                ..Default::default()
            };
            std::fs::write(sidecar_path(&obj), side.to_text()).unwrap();
        }
        let mut f = fetcher(u64::MAX);
        f.cache_dir = Some(dir.clone());
        f.revalidate = true;
        (dir, f)
    }

    /// 304 is the answer a warm run wants: the cache holds what the bucket still serves, and
    /// no body crossed the wire to establish it.
    #[test]
    fn a_cache_entry_the_bucket_still_serves_is_current() {
        let body = zstd::encode_all(&b">REVAL1_0\nACGT\n"[..], 3).unwrap();
        let (dir, f) = cached_with("current", Some("tag-abc"), &body);
        let url = serve_status("304 Not Modified", "", Vec::new());
        let got = f.check_cached(0, &"REVAL1".to_string(), &cache_path(&dir, "REVAL1"), &url);
        assert!(matches!(got, Revalidation::Current), "a 304 must read as current");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// 200 to a conditional request means the object changed. The new bytes arrive in the
    /// same round trip and replace the cache entry, which is the point: the next run, with
    /// or without `--revalidate`, sees the current object rather than the stale one.
    #[test]
    fn a_cache_entry_the_bucket_replaced_comes_back_with_the_new_bytes() {
        let old = zstd::encode_all(&b">REVAL1_0\nACGT\n"[..], 3).unwrap();
        let new = zstd::encode_all(&b">REVAL1_0\nTTTTTTTT\n"[..], 3).unwrap();
        assert_ne!(old, new);
        let (dir, f) = cached_with("replaced", Some("tag-old"), &old);
        let url = serve_status("200 OK", "ETag: \"tag-new\"\r\n", new.clone());

        let got = f.check_cached(0, &"REVAL1".to_string(), &cache_path(&dir, "REVAL1"), &url);
        match got {
            Revalidation::Replaced(fresh) => {
                assert_eq!(fresh.compressed_len, new.len());
                assert_eq!(fresh.payload.etag.as_deref(), Some("tag-new"));
            }
            _ => panic!("a 200 to a conditional request must read as replaced"),
        }
        // The entry on disk is the new object, and its sidecar names the new ETag: a later
        // run must not find the bytes this one just proved stale.
        let obj = cache_path(&dir, "REVAL1");
        assert_eq!(std::fs::read(&obj).unwrap(), new, "the cache must hold the new bytes");
        let side = read_sidecar(&obj).expect("sidecar rewritten");
        assert_eq!(side.etag.as_deref(), Some("tag-new"));
        assert!(side.describes(Some(&crate::types::sha256_hex(&new))));
        std::fs::remove_dir_all(&dir).ok();
    }

    /// With nothing recording what the cached bytes are a copy of, there is no question to
    /// ask. Every entry written before sidecars existed is in this state.
    #[test]
    fn a_cache_entry_with_no_recorded_origin_cannot_be_revalidated() {
        let body = zstd::encode_all(&b">REVAL1_0\nACGT\n"[..], 3).unwrap();
        let (dir, f) = cached_with("noorigin", None, &body);
        let got = f.check_cached(
            0,
            &"REVAL1".to_string(),
            &cache_path(&dir, "REVAL1"),
            "http://127.0.0.1:1/never-used",
        );
        assert!(matches!(got, Revalidation::NoOrigin));
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A revalidation that cannot be made is not a failed accession.
    ///
    /// The cached bytes are still a usable analysis, and turning a network blip into a lost
    /// cohort would be a worse answer than running on them and saying the check did not
    /// happen. Recorded as `unavailable` so the summary never claims a check it did not make.
    #[test]
    fn a_revalidation_that_cannot_be_made_falls_back_to_the_cached_bytes() {
        let body = zstd::encode_all(&b">REVAL1_0\nACGT\n"[..], 3).unwrap();
        let (dir, f) = cached_with("unavail", Some("tag-abc"), &body);
        // Nothing listening: a transport error, not an HTTP answer.
        let got = f.check_cached(
            0,
            &"REVAL1".to_string(),
            &cache_path(&dir, "REVAL1"),
            "http://127.0.0.1:1/dead",
        );
        assert!(matches!(got, Revalidation::Unavailable), "must not fail the accession");
        std::fs::remove_dir_all(&dir).ok();
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
        let body = zbody(32);
        let url = serve_once(
            "ETag: \"d2dd288b238c3c88b46acdcd78860dce\"\r\nLast-Modified: Sun, 31 Mar 2024 17:38:36 GMT\r\n",
            body.clone(),
        );
        let got = attempt(&fetcher(1024), "X", &url).ok().expect("fetch");
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
        let body = zbody(8);
        let got = attempt(&fetcher(1024), "X", &serve_once("", body.clone()))
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
        match attempt(&fetcher(10), "X", &url) {
            Err(FetchAttempt::TooLarge(_, p)) => {
                assert_eq!(p.bytes, Some(100), "the size that tripped the cap must be a number");
                assert_eq!(p.etag.as_deref(), Some("abc123"));
                assert_eq!(p.source, PayloadSource::S3);
                assert!(p.sha256.is_none(), "no bytes were read, so there is nothing to hash");
            }
            Err(FetchAttempt::Transient(d, _)) => panic!("retried as transient: {d}"),
            Err(_) => panic!("wrong failure"),
            Ok(_) => panic!("accepted 100 bytes over a 10-byte cap"),
        }
    }

    #[test]
    fn a_body_of_exactly_the_cap_is_accepted() {
        // Exactly the cap, and still a zstd stream: the magic is four of the ten bytes.
        let url = serve_without_length(zbody(10));
        let got = attempt(&fetcher(10), "X", &url).ok().expect("refused at the cap");
        assert_eq!(got.compressed_len, 10);
    }
}

#[cfg(test)]
mod retry_tests {
    use super::*;
    use crate::config::Config;
    use std::io::{BufRead, BufReader, Write};
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;

    fn zbody(n: usize) -> Vec<u8> {
        let mut v = vec![0x28, 0xb5, 0x2f, 0xfd];
        v.resize(n.max(4), b'Z');
        v
    }

    /// A fetcher that talks to the loopback server directly, with no proxy from the env.
    fn fetcher(cfg: Config) -> Fetcher {
        let mut f = Fetcher::new(&cfg);
        f.agent = ureq::Agent::new_with_config(
            ureq::Agent::config_builder()
                .proxy(None)
                .http_status_as_error(false)
                .timeout_global(Some(Duration::from_secs(cfg.fetch_timeout_secs)))
                .build(),
        );
        f
    }

    /// The `Range: bytes=N-` of a request, if it has one.
    fn read_request(sock: &std::net::TcpStream) -> Option<u64> {
        let mut reader = BufReader::new(sock.try_clone().unwrap());
        let mut line = String::new();
        let mut from = None;
        while reader.read_line(&mut line).unwrap_or(0) > 0 && line != "\r\n" {
            if let Some(v) = line.to_ascii_lowercase().strip_prefix("range: bytes=") {
                from = v.split('-').next().and_then(|n| n.trim().parse().ok());
            }
            line.clear();
        }
        from
    }

    /// Serve `body`, closing the connection after `chunk` bytes of each response, and
    /// honouring `Range`. The returned counter is the number of connections made.
    fn serve_resumable(body: Vec<u8>, chunk: usize) -> (String, Arc<AtomicUsize>) {
        let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
        let port = listener.local_addr().unwrap().port();
        let seen = Arc::new(AtomicUsize::new(0));
        let counter = Arc::clone(&seen);
        std::thread::spawn(move || {
            for sock in listener.incoming() {
                let Ok(mut sock) = sock else { break };
                counter.fetch_add(1, Ordering::Relaxed);
                let from = read_request(&sock).unwrap_or(0) as usize;
                if from >= body.len() {
                    let _ = sock.write_all(b"HTTP/1.1 416 Range Not Satisfiable\r\nContent-Length: 0\r\nConnection: close\r\n\r\n");
                    continue;
                }
                let rest = &body[from..];
                let head = if from == 0 {
                    format!(
                        "HTTP/1.1 200 OK\r\nETag: \"abc\"\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                        rest.len()
                    )
                } else {
                    format!(
                        "HTTP/1.1 206 Partial Content\r\nETag: \"abc\"\r\nContent-Length: {}\r\n\
                         Content-Range: bytes {from}-{}/{}\r\nConnection: close\r\n\r\n",
                        rest.len(),
                        body.len() - 1,
                        body.len()
                    )
                };
                let _ = sock.write_all(head.as_bytes());
                // Short by design when there is more to come: this is the dropped connection.
                let _ = sock.write_all(&rest[..rest.len().min(chunk)]);
            }
        });
        (format!("http://127.0.0.1:{port}/u/X/X.unitigs.fa.zst"), seen)
    }

    /// Serve one response with a fixed head, for the status-code cases.
    fn serve_head(head: String, body: Vec<u8>) -> String {
        let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
        let port = listener.local_addr().unwrap().port();
        std::thread::spawn(move || {
            let (mut sock, _) = listener.accept().unwrap();
            let _ = read_request(&sock);
            let _ = sock.write_all(head.as_bytes());
            let _ = sock.write_all(&body);
        });
        format!("http://127.0.0.1:{port}/u/X/X.unitigs.fa.zst")
    }

    /// A transfer interrupted after real progress continues where it stopped.
    ///
    /// Before this, every attempt reissued a bare GET from byte 0, so an object whose
    /// per-stream rate fell below size/timeout could never be fetched however many attempts
    /// it was given: each one timed out in the same place. The chunk here is over
    /// `RESUME_PROGRESS_BYTES`, so the resumes are free and the run does not back off.
    #[test]
    fn an_interrupted_download_resumes_instead_of_starting_again() {
        let chunk = RESUME_PROGRESS_BYTES as usize + 1;
        let body = zbody(chunk * 3);
        let (url, seen) = serve_resumable(body.clone(), chunk);
        // `download` builds the URL from the accession, so the retry shape is driven here
        // against the loopback server instead.
        let f = fetcher(Config::default());
        let mut xfer = Transfer::new(f.open_sink(&"X".to_string()));
        let mut attempts = 0;
        let out = loop {
            attempts += 1;
            match f.try_once(0, &"X".to_string(), &url, None, &mut xfer) {
                Ok(v) => break v,
                Err(FetchAttempt::Transient(d, _)) if attempts < 10 => {
                    assert!(xfer.got > 0, "nothing was kept across the attempt: {d}");
                }
                Err(FetchAttempt::Transient(d, _)) => panic!("gave up: {d}"),
                Err(_) => panic!("wrong failure"),
            }
        };
        assert_eq!(out.compressed_len, body.len(), "the whole object must arrive");
        assert_eq!(
            out.payload.sha256.as_deref(),
            Some(crate::types::sha256_hex(&body).as_str()),
            "a resumed object must hash as the original, not as a concatenation"
        );
        assert_eq!(seen.load(Ordering::Relaxed), 3, "one connection per chunk");
    }

    /// A server that answers a resume with the whole object starts the transfer over, rather
    /// than appending the object to the part already held.
    #[test]
    fn a_200_answer_to_a_resume_starts_over() {
        let body = zbody(300);
        let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
        let port = listener.local_addr().unwrap().port();
        let whole = body.clone();
        std::thread::spawn(move || {
            let mut n = 0;
            for sock in listener.incoming() {
                let Ok(mut sock) = sock else { break };
                let _ = read_request(&sock);
                n += 1;
                let head = format!(
                    "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                    whole.len()
                );
                let _ = sock.write_all(head.as_bytes());
                // First connection: half the object, then hang up. Second: all of it.
                let send = if n == 1 { 100 } else { whole.len() };
                let _ = sock.write_all(&whole[..send]);
            }
        });
        let url = format!("http://127.0.0.1:{port}/u/X/X.unitigs.fa.zst");
        let f = fetcher(Config::default());
        let mut xfer = Transfer::new(f.open_sink(&"X".to_string()));
        let first = f.try_once(0, &"X".to_string(), &url, None, &mut xfer);
        assert!(matches!(first, Err(FetchAttempt::Transient(..))), "expected a short body");
        assert_eq!(xfer.got, 100, "the first attempt kept what it read");
        let out = f.try_once(0, &"X".to_string(), &url, None, &mut xfer).ok().expect("second");
        assert_eq!(out.compressed_len, body.len(), "400 bytes would mean the parts were glued");
        assert_eq!(
            out.payload.sha256.as_deref(),
            Some(crate::types::sha256_hex(&body).as_str())
        );
    }

    /// 408 and 429 are transient, and a `Retry-After` is obeyed.
    ///
    /// They used to fall into the `400..=499` arm and fail the accession on the first try,
    /// so a throttling bucket cost one accession per 429 rather than one pause.
    #[test]
    fn throttling_is_transient_and_honours_retry_after() {
        let f = fetcher(Config::default());
        for (status, extra, want) in [
            ("429 Too Many Requests", "Retry-After: 7\r\n", Some(7)),
            ("429 Too Many Requests", "", None),
            ("408 Request Timeout", "", None),
            // Not delta-seconds, so the backoff decides instead of a bad parse.
            ("429 Too Many Requests", "Retry-After: Wed, 21 Oct 2026 07:28:00 GMT\r\n", None),
            // Absurd values are capped rather than trusted.
            ("429 Too Many Requests", "Retry-After: 99999\r\n", Some(300)),
        ] {
            let url = serve_head(
                format!("HTTP/1.1 {status}\r\n{extra}Content-Length: 0\r\nConnection: close\r\n\r\n"),
                Vec::new(),
            );
            let mut xfer = Transfer::new(f.open_sink(&"X".to_string()));
            match f.try_once(0, &"X".to_string(), &url, None, &mut xfer) {
                Err(FetchAttempt::Transient(_, after)) => assert_eq!(
                    after.map(|d| d.as_secs()),
                    want,
                    "{status} {extra:?}"
                ),
                Err(FetchAttempt::Permanent(d)) => panic!("{status}: permanent: {d}"),
                _ => panic!("{status}: wrong outcome"),
            }
        }
    }

    /// The backoff grows, is capped, and is not the same number twice.
    #[test]
    fn backoff_grows_and_is_capped() {
        // 1, 2, 4, 8, 16, 32 seconds before jitter, then flat.
        for (attempt, base) in [(1u32, 1u64), (2, 2), (3, 4), (4, 8), (5, 16), (6, 32), (9, 32)] {
            let d = backoff(attempt).as_millis() as u64;
            assert!(
                d >= base * 500 && d <= base * 1500,
                "attempt {attempt}: {d} ms is outside +-50% of {base} s"
            );
        }
        let a: Vec<u64> = (0..8).map(|_| backoff(3).as_millis() as u64).collect();
        assert!(a.windows(2).any(|w| w[0] != w[1]), "no jitter at all: {a:?}");
    }

    /// A body that is not zstd is refused, and never cached.
    ///
    /// A proxy's HTML error page must not be written into the cache under the accession's
    /// name, where it would fail that accession with `parse failed: Unknown frame
    /// descriptor` on every later run until someone found and deleted the file.
    #[test]
    fn a_body_that_is_not_zstd_is_refused_and_not_cached() {
        let dir = std::env::temp_dir().join(format!("logan-wf13-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let page = b"<html><head><title>502 Bad Gateway</title></head></html>".to_vec();
        let url = serve_head(
            format!(
                "HTTP/1.1 200 OK\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                page.len()
            ),
            page,
        );
        let f = fetcher(Config { unitigs_cache: Some(dir.clone()), ..Config::default() });
        let mut xfer = Transfer::new(f.open_sink(&"X".to_string()));
        match f.try_once(0, &"X".to_string(), &url, None, &mut xfer) {
            Err(FetchAttempt::Transient(d, _)) => assert!(d.contains("not zstd"), "{d}"),
            Err(FetchAttempt::Permanent(d)) => panic!("not retryable: {d}"),
            _ => panic!("an HTML error page was accepted as a Logan object"),
        }
        drop(xfer);
        assert!(!cache_path(&dir, "X").exists(), "the error page was cached");
        let left: Vec<_> = std::fs::read_dir(&dir).unwrap().flatten().map(|e| e.file_name()).collect();
        assert!(left.is_empty(), "something was left in the cache: {left:?}");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A cache entry that is not zstd is refused by name under `--offline`, and replaced
    /// when the bucket can be asked.
    #[test]
    fn a_cache_entry_that_is_not_zstd_is_named_offline_and_replaced_online() {
        let dir = std::env::temp_dir().join(format!("logan-wf13b-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = cache_path(&dir, "BADENTRY01");
        std::fs::write(&path, b"<html>502</html>").unwrap();

        let offline = fetcher(Config {
            unitigs_cache: Some(dir.clone()),
            offline: true,
            ..Config::default()
        });
        match offline.fetch(0, &"BADENTRY01".to_string()) {
            Err(Error::Fetch { detail, .. }) => {
                assert!(detail.contains("not zstd data"), "{detail}");
                assert!(detail.contains(&path.display().to_string()), "must name the file: {detail}");
            }
            other => panic!("expected a named refusal, got {:?}", other.map(|f| f.compressed_len)),
        }
        assert!(path.exists(), "--offline must not delete the user's file");

        // Online the entry is removed so the next download can replace it. Driven through
        // the decision itself rather than through `fetch`, which would go to the real bucket
        // and this suite stays offline.
        let side = sidecar_path(&path);
        std::fs::write(&side, b"schema 1\n").unwrap();
        let online = fetcher(Config { unitigs_cache: Some(dir.clone()), ..Config::default() });
        online.discard_bad_entry(&"BADENTRY01".to_string(), &path).expect("online discards");
        assert!(!path.exists(), "an online run must remove the entry it cannot use");
        assert!(!side.exists(), "its provenance describes a file that is gone");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A negative recorded at fetch time keeps its status in the `--offline` run.
    ///
    /// `disassembler fetch` writes nothing for an accession Logan has no unitigs for, so the
    /// documented `--offline` run that follows saw an ordinary cache miss and called it an
    /// error: the same list and the same flags gave `no_logan_data` and exit 0 online, and
    /// `error` plus exit 1 or 3 offline. Since 0.2.0 any error makes the run non-zero, so a
    /// shard whose accessions are all 404s fails outright.
    #[test]
    fn an_offline_miss_with_a_negative_note_keeps_its_status() {
        let dir = std::env::temp_dir().join(format!("logan-wf03-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let cfg = Config {
            unitigs_cache: Some(dir.clone()),
            offline: true,
            max_compressed_bytes: 1000,
            ..Config::default()
        };

        write_absent(&dir, "NODATA0001", "no_logan_data");
        match fetcher(cfg.clone()).fetch(0, &"NODATA0001".to_string()) {
            Err(Error::NoLoganData { accession }) => assert_eq!(accession, "NODATA0001"),
            other => panic!("expected NoLoganData, got {:?}", other.map(|f| f.compressed_len)),
        }

        write_absent(&dir, "TOOBIG0001", "too_large 5000");
        match fetcher(cfg.clone()).fetch(0, &"TOOBIG0001".to_string()) {
            Err(Error::TooLarge { detail, payload, .. }) => {
                assert_eq!(payload.bytes, Some(5000));
                assert!(detail.contains("5000"), "{detail}");
            }
            other => panic!("expected TooLarge, got {:?}", other.map(|f| f.compressed_len)),
        }

        // The cap was raised since the cache was filled, so the note no longer applies and
        // this is an ordinary miss again: the object really is not here.
        let raised = Config { max_compressed_bytes: 1_000_000, ..cfg.clone() };
        match fetcher(raised).fetch(0, &"TOOBIG0001".to_string()) {
            Err(Error::Fetch { detail, .. }) => assert!(detail.contains("cache miss"), "{detail}"),
            other => panic!("expected a plain miss, got {:?}", other.map(|f| f.compressed_len)),
        }

        // A damaged note says nothing, rather than inventing an outcome.
        std::fs::write(absent_path(&dir, "GARBLED001"), b"\x00\x01 not a note at all\n").unwrap();
        match fetcher(cfg).fetch(0, &"GARBLED001".to_string()) {
            Err(Error::Fetch { detail, .. }) => assert!(detail.contains("cache miss"), "{detail}"),
            other => panic!("expected a plain miss, got {:?}", other.map(|f| f.compressed_len)),
        }
        std::fs::remove_dir_all(&dir).ok();
    }

    /// What each fetch outcome does to the note.
    ///
    /// A success clears it, because Logan does add data for accessions it had none for, and
    /// because a cap raised between runs must not keep refusing. A transport failure leaves
    /// it alone: it says nothing about the accession.
    #[test]
    fn the_note_follows_the_outcome_and_a_failure_does_not_invent_one() {
        let dir = std::env::temp_dir().join(format!("logan-wf03b-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let acc = "ONCEGONE01";
        let note = absent_path(&dir, acc);

        note_outcome(&dir, acc, &Err(Error::NoLoganData { accession: acc.into() }), 1000);
        assert_eq!(std::fs::read_to_string(&note).unwrap().trim(), "no_logan_data");

        note_outcome(
            &dir,
            acc,
            &Err(Error::TooLarge {
                accession: acc.into(),
                detail: String::new(),
                payload: Box::new(Payload { bytes: Some(4096), ..Payload::default() }),
            }),
            1000,
        );
        assert_eq!(std::fs::read_to_string(&note).unwrap().trim(), "too_large 4096");

        // Nothing is known about the accession from a dropped connection, so the note stands.
        note_outcome(
            &dir,
            acc,
            &Err(Error::Fetch { accession: acc.into(), detail: "connection reset".into() }),
            1000,
        );
        assert_eq!(std::fs::read_to_string(&note).unwrap().trim(), "too_large 4096");

        // A size with no Content-Length behind it is recorded as "more than the cap".
        note_outcome(
            &dir,
            acc,
            &Err(Error::TooLarge {
                accession: acc.into(),
                detail: String::new(),
                payload: Box::new(Payload::default()),
            }),
            1000,
        );
        assert_eq!(std::fs::read_to_string(&note).unwrap().trim(), "too_large 1001");

        note_outcome(
            &dir,
            acc,
            &Ok(FetchedUnitigs {
                index: 0,
                accession: acc.into(),
                compressed_len: 4,
                compressed: Compressed::Memory(zbody(4)),
                fetch_secs: 0.0,
                payload: Payload::default(),
            }),
            1000,
        );
        assert!(!note.exists(), "a successful fetch must clear the note");
        std::fs::remove_dir_all(&dir).ok();
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
        // A cache entry has to look like a zstd object before anything else is asked of it.
        let mut big = vec![0x28u8, 0xb5, 0x2f, 0xfd];
        big.resize(16 * 1024 * 1024, 0);
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
        let mut small = vec![0x28u8, 0xb5, 0x2f, 0xfd];
        small.resize(1024 * 1024, 0);
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

    /// A cache entry with no sidecar records its own bytes and claims nothing more.
    ///
    /// This is every entry written before 0.2.0, and the reason the digest is computed from
    /// the bytes at all: a cache hit has no HTTP response to take an `ETag` from, and
    /// `--offline` runs are the ones whose inputs are least self-evident. A sidecar, when
    /// there is one, fills the rest in; see the tests below.
    #[test]
    fn a_cache_hit_with_no_sidecar_records_the_bytes_it_read() {
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
        assert_eq!(p.etag, None, "no sidecar, so nothing says what this is a copy of");
        assert_eq!(p.fetched_at, None, "and nothing says when it was fetched");
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A sidecar beside a cached object tells a warm run what the bytes are a copy of.
    ///
    /// Without it a cache hit knows its own digest and nothing else, so it cannot say which
    /// bucket object it came from and `--revalidate` has no `ETag` to send. Logan's `u/`
    /// objects are mutable, which is what makes that gap matter.
    #[test]
    fn a_sidecar_gives_a_cache_hit_the_origin_its_bytes_came_from() {
        let dir = tempdir("sidecar-read");
        let bytes = zstd::encode_all(&b">CACHED0009_0\nACGT\n"[..], 3).unwrap();
        let obj = cache_path(&dir, "CACHED0009");
        std::fs::write(&obj, &bytes).unwrap();
        let side = crate::types::Sidecar {
            bytes: Some(bytes.len() as u64),
            sha256: Some(crate::types::sha256_hex(&bytes)),
            etag: Some("d2dd288b238c3c88b46acdcd78860dce".into()),
            last_modified: Some("2024-03-31T17:38:36Z".into()),
            fetched_at: Some("2026-11-17T22:08:24Z".into()),
        };
        std::fs::write(sidecar_path(&obj), side.to_text()).unwrap();

        let f = Fetcher::new(&Config {
            unitigs_cache: Some(dir.clone()),
            offline: true,
            ..Config::default()
        });
        let p = f.fetch(0, &"CACHED0009".to_string()).expect("cache hit").payload;
        assert_eq!(p.source, PayloadSource::Cache, "the bytes still came from the cache");
        assert_eq!(p.etag.as_deref(), Some("d2dd288b238c3c88b46acdcd78860dce"));
        assert_eq!(p.last_modified.as_deref(), Some("2024-03-31T17:38:36Z"));
        assert_eq!(p.fetched_at.as_deref(), Some("2026-11-17T22:08:24Z"));
        // Measured from the file, not adopted: the sidecar is believed about provenance
        // precisely because it is checked against these.
        assert_eq!(p.bytes, Some(bytes.len() as u64));
        assert_eq!(p.sha256.as_deref(), Some(crate::types::sha256_hex(&bytes).as_str()));
        std::fs::remove_dir_all(&dir).ok();
    }

    /// A sidecar that describes other bytes is ignored, not believed.
    ///
    /// This is the case the digest exists for: someone replaced the object and left the old
    /// sidecar, so the recorded origin now names something these bytes are not a copy of.
    /// Believing it would attach a false `ETag` to a real accession, which is worse than
    /// having none, and would make `--revalidate` compare against the wrong object.
    #[test]
    fn a_sidecar_that_does_not_match_the_object_is_ignored() {
        let dir = tempdir("sidecar-stale");
        let bytes = zstd::encode_all(&b">CACHED0010_0\nACGT\n"[..], 3).unwrap();
        let obj = cache_path(&dir, "CACHED0010");
        std::fs::write(&obj, &bytes).unwrap();
        let stale = crate::types::Sidecar {
            sha256: Some("0".repeat(64)),
            etag: Some("stale-etag".into()),
            ..Default::default()
        };
        std::fs::write(sidecar_path(&obj), stale.to_text()).unwrap();

        let f = Fetcher::new(&Config {
            unitigs_cache: Some(dir.clone()),
            offline: true,
            ..Config::default()
        });
        let p = f.fetch(0, &"CACHED0010".to_string()).expect("cache hit").payload;
        assert_eq!(p.etag, None, "a sidecar for other bytes must not be adopted");
        assert_eq!(p.sha256.as_deref(), Some(crate::types::sha256_hex(&bytes).as_str()));
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
        sink.finish(&Payload::default()).unwrap();
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
            FetchAttempt::Transient(..) => "Transient",
            FetchAttempt::Permanent(_) => "Permanent",
            FetchAttempt::Missing => "Missing",
            FetchAttempt::TooLarge(..) => "TooLarge",
            FetchAttempt::NotModified => "NotModified",
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
        let compressed = sink.finish(&Payload::default()).unwrap();

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
