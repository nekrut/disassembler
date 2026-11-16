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
}

impl Fetcher {
    pub fn new(cfg: &Config) -> Self {
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
        }
    }

    /// Download one accession's compressed unitig file into memory.
    ///
    /// With `--unitigs-cache`: read from `<dir>/<ACC>.unitigs.fa.zst` if present, otherwise
    /// (under `--offline`) fail, otherwise download and write back atomically. A cached file's
    /// contents are trusted, because a truncated file could only have arrived from a killed
    /// writer, and the atomic rename below prevents that. Its size is not: `--max-size-mb`
    /// bounds memory, and a cache filled under a higher cap used to be read whole regardless.
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
                if len > self.max_compressed_bytes {
                    return Err(Error::TooLarge {
                        accession: accession.clone(),
                        detail: format!(
                            "cached file {} is {len} bytes, over max_compressed_bytes {}",
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
                Ok(mut f) => {
                    if let Some(dir) = &self.cache_dir {
                        match f.compressed.as_slice() {
                            Some(bytes) => match write_cache_atomic(dir, accession, bytes) {
                                // Written, so the bytes exist on disk and there is no reason
                                // to keep carrying them. Swapping the buffer for the handle
                                // drops this accession's compressed size before it reaches
                                // the compute worker, which is where a buffer is held
                                // longest: for the whole parse. A reopen that fails is not
                                // worth failing the fetch over, so the buffer simply stays.
                                Ok(path) => match std::fs::File::open(&path) {
                                    Ok(file) => f.compressed = Compressed::File(file),
                                    Err(e) => eprintln!(
                                        "[{accession}] warning: cached but could not reopen \
                                         {}: {e}",
                                        path.display()
                                    ),
                                },
                                // Cache-write failure must not lose a successful fetch. Warn
                                // and hand the bytes back so the compute stage still runs;
                                // the next run will just have to fetch again.
                                Err(e) => eprintln!(
                                    "[{accession}] warning: could not write cache entry: {e}"
                                ),
                            },
                            // Already file-backed: nothing downloaded it into memory.
                            None => {}
                        }
                    }
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
        }

        // ureq caps read_to_vec at 10 MB by default, which every real accession exceeds.
        // The limit here doubles as the guard for a missing or lying Content-Length. ureq
        // refuses a body that reaches its limit, so the limit is one byte over the cap: a file
        // of exactly the cap passes, as it does the Content-Length check above.
        //
        // Running into the limit is the size cap, not a transient failure. It used to be
        // retried, downloading up to the cap three times before reporting a fetch error.
        let compressed = resp
            .body_mut()
            .with_config()
            .limit(self.max_compressed_bytes.saturating_add(1))
            .read_to_vec()
            .map_err(|e| match e {
                ureq::Error::BodyExceedsLimit(_) => FetchAttempt::TooLarge(
                    format!(
                        "compressed size exceeds max_compressed_bytes {} (the server sent no \
                         Content-Length, so this was found by downloading)",
                        self.max_compressed_bytes
                    ),
                    // All that is known here is "over the cap". A Content-Length cannot
                    // have been both present and honest, or the check above would have
                    // fired, so its value is dropped rather than recorded as a size.
                    Box::new(Payload { bytes: None, ..payload.clone() }),
                ),
                e => FetchAttempt::Transient(format!("body read failed: {e}")),
            })?;

        if compressed.is_empty() {
            return Err(FetchAttempt::Transient("empty body".to_string()));
        }

        // Content-Length is the server's claim; this is the bytes that were actually read,
        // and it is what every downstream stage consumed.
        payload.bytes = Some(compressed.len() as u64);
        payload.sha256 = Some(crate::types::sha256_hex(&compressed));

        Ok(FetchedUnitigs {
            index,
            accession: accession.clone(),
            compressed_len: compressed.len(),
            compressed: Compressed::Memory(compressed),
            fetch_secs: 0.0,
            payload,
        })
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

/// Write `bytes` to `<dir>/<ACC>.unitigs.fa.zst` via a tempfile-and-rename so a killed
/// process cannot leave a partial file behind. A partial file that survived a crash would
/// otherwise look valid to a subsequent `--offline` run and corrupt the analysis.
fn write_cache_atomic(dir: &Path, accession: &str, bytes: &[u8]) -> std::io::Result<PathBuf> {
    std::fs::create_dir_all(dir)?;
    let final_path = cache_path(dir, accession);
    // Tempfile in the same directory so `rename` is same-filesystem and atomic. PID + nanos
    // is enough uniqueness for one host writing its own cache; two processes racing on the
    // same accession may both fetch and one rename wins, which is fine (same bytes).
    let ns = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let tmp = dir.join(format!(".{accession}.{}.{ns}.tmp", std::process::id()));
    // Written-then-renamed. Errors during write must not leave the tempfile around either;
    // the finally-style branch below handles that.
    let write_result = std::fs::write(&tmp, bytes);
    match write_result {
        Ok(()) => std::fs::rename(&tmp, &final_path)
            .inspect_err(|_| {
                let _ = std::fs::remove_file(&tmp);
            })
            .map(|()| final_path),
        Err(e) => {
            let _ = std::fs::remove_file(&tmp);
            Err(e)
        }
    }
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
            got.compressed.as_slice().expect("no cache was set, so this is in memory").len()
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
        let msg = err.to_string();
        eprintln!("{msg}");
        assert!(msg.contains("404"), "expected a 404, got: {msg}");
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

    /// `--max-size-mb` bounds memory, so it applies to a cached file as it does to a
    /// download. A cache filled under a higher cap used to be read whole whatever the cap.
    #[test]
    fn a_cached_file_over_the_size_cap_is_too_large() {
        let dir = tempdir("cap");
        let bytes = zstd::encode_all(&b">CACHED0002_0\nACGTACGTACGT\n"[..], 3).unwrap();
        std::fs::write(cache_path(&dir, "CACHED0002"), &bytes).unwrap();
        let fetcher = |cap: u64| Fetcher::new(&Config {
            unitigs_cache: Some(dir.clone()),
            offline: true,
            max_compressed_bytes: cap,
            ..Config::default()
        });
        let acc = "CACHED0002".to_string();
        match fetcher(bytes.len() as u64 - 1).fetch(0, &acc) {
            Err(Error::TooLarge { detail, .. }) => assert!(detail.contains("cached file"), "{detail}"),
            other => panic!("expected TooLarge, got {:?}", other.map(|f| f.compressed_len)),
        }
        let got = fetcher(bytes.len() as u64).fetch(0, &acc).unwrap();
        assert_eq!(decoded(&got), b">CACHED0002_0\nACGTACGTACGT\n");
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

    /// `write_cache_atomic` uses a tempfile + rename, so a mid-write crash cannot leave a
    /// partial file at the final path.
    #[test]
    fn atomic_write_lands_the_final_file_and_leaves_no_tempfile() {
        let dir = tempdir("atomic");
        write_cache_atomic(&dir, "ATOMIC1", b"payload").unwrap();
        let final_path = cache_path(&dir, "ATOMIC1");
        assert!(final_path.exists());
        assert_eq!(std::fs::read(&final_path).unwrap(), b"payload");
        // No stray .tmp files linger.
        let leftover: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_name().to_string_lossy().ends_with(".tmp"))
            .collect();
        assert!(leftover.is_empty(), "leftover tempfiles: {leftover:?}");
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
