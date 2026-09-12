//! Step 1: fetch one accession's unitigs from the Logan S3 bucket.
//!
//! Replaces `logan_path_walker.py:45-67` (`stream_unitigs`), which shells out to
//! `curl … | zstd -dc`, writes the decompressed file to disk, and then reads the whole
//! file back a third time purely to count `>` characters for a reporting field.
//!
//! Measured on this hardware: 6.9 MB/s single stream, 17.4 MB/s aggregate across 8 parallel
//! range requests, which is where the link saturates. A median accession is ~71 MB compressed,
//! so this stage has a hard floor of roughly 4.1 s that no amount of code can remove.
//!
//! IN   accession
//! OUT  FetchedUnitigs (compressed bytes only), plus a streaming decoder over them

use crate::config::Config;
use crate::error::{Error, Result};
use crate::types::{Accession, FetchedUnitigs};
use std::io::Read;
use std::time::Duration;

/// URL for an accession's unitig file.
pub fn unitigs_url(accession: &str) -> String {
    format!("https://logan-pub.s3.amazonaws.com/u/{accession}/{accession}.unitigs.fa.zst")
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
        }
    }

    /// Download one accession's compressed unitig file into memory.
    ///
    /// Returns the compressed bytes only. Decompression is `decoder`'s job, and happens in
    /// the compute pool so the fetch thread can immediately start the next transfer.
    pub fn fetch(&self, accession: &Accession) -> Result<FetchedUnitigs> {
        let url = unitigs_url(accession);
        let mut last: String = String::from("no attempt made");

        for attempt in 1..=self.attempts {
            match self.try_once(accession, &url) {
                Ok(f) => return Ok(f),
                Err(FetchAttempt::Permanent(detail)) => {
                    return Err(Error::Fetch {
                        accession: accession.clone(),
                        detail,
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

    fn try_once(&self, accession: &Accession, url: &str) -> std::result::Result<FetchedUnitigs, FetchAttempt> {
        let mut resp = self
            .agent
            .get(url)
            .call()
            .map_err(|e| FetchAttempt::Transient(format!("request failed: {e}")))?;

        let status = resp.status().as_u16();
        match status {
            200 => {}
            404 => {
                return Err(FetchAttempt::Permanent(
                    "no Logan unitigs for this accession (HTTP 404)".to_string(),
                ))
            }
            // 4xx other than 404 is our fault and will not improve on retry.
            400..=499 => return Err(FetchAttempt::Permanent(format!("HTTP {status}"))),
            _ => return Err(FetchAttempt::Transient(format!("HTTP {status}"))),
        }

        // Refuse oversized accessions before transferring them, when the server says how big.
        if let Some(len) = resp.body().content_length() {
            if len > self.max_compressed_bytes {
                return Err(FetchAttempt::Permanent(format!(
                    "compressed size {len} exceeds max_compressed_bytes {}",
                    self.max_compressed_bytes
                )));
            }
        }

        // ureq caps read_to_vec at 10 MB by default, which every real accession exceeds.
        // The limit here doubles as the guard for a missing or lying Content-Length.
        let compressed = resp
            .body_mut()
            .with_config()
            .limit(self.max_compressed_bytes)
            .read_to_vec()
            .map_err(|e| FetchAttempt::Transient(format!("body read failed: {e}")))?;

        if compressed.is_empty() {
            return Err(FetchAttempt::Transient("empty body".to_string()));
        }

        Ok(FetchedUnitigs {
            accession: accession.clone(),
            compressed_len: compressed.len(),
            compressed,
            fetch_secs: 0.0,
        })
    }
}

/// Whether a failed attempt is worth repeating.
enum FetchAttempt {
    /// Transport error or 5xx. Retry.
    Transient(String),
    /// 404, oversized, or another 4xx. Retrying cannot help.
    Permanent(String),
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
    // Cursor<&[u8]> is already BufRead, so `with_buffer` avoids the redundant BufReader
    // that `Decoder::new` would insert over an in-memory slice.
    let cursor = std::io::Cursor::new(fetched.compressed.as_slice());
    let decoder = zstd::stream::read::Decoder::with_buffer(cursor).map_err(|e| {
        Error::Decompress {
            accession: fetched.accession.clone(),
            detail: format!("could not start zstd decoder: {e}"),
        }
    })?;
    Ok(Box::new(decoder))
}

/// Suggested read size for the decoder returned by [`decoder`].
pub fn recommended_chunk() -> usize {
    zstd::stream::read::Decoder::<'_, std::io::Cursor<&[u8]>>::recommended_output_size()
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
        let got = f.fetch(&"SRR8209082".to_string()).expect("fetch");
        eprintln!("fetched {} compressed bytes", got.compressed_len);
        assert_eq!(got.compressed_len, 1_861_058, "unexpected compressed size");
        assert_eq!(got.compressed_len, got.compressed.len());

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
        let err = f.fetch(&"NOTAREALACCESSION999".to_string()).err().expect("should fail");
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
        let err = f.fetch(&"SRR12016917".to_string()).err().expect("should fail");
        let msg = err.to_string();
        eprintln!("{msg}");
        assert!(msg.contains("exceeds max_compressed_bytes"), "got: {msg}");
        assert!(
            started.elapsed() < Duration::from_secs(20),
            "must reject on Content-Length, not after downloading 3.9 GB"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fetched(compressed: Vec<u8>) -> FetchedUnitigs {
        FetchedUnitigs {
            accession: "TEST0001".to_string(),
            compressed_len: compressed.len(),
            compressed,
            fetch_secs: 0.0,
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
            accession: "REAL".to_string(),
            compressed_len: compressed.len(),
            compressed,
            fetch_secs: 0.0,
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
