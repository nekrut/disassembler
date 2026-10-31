//! The walker binary end to end, on accessions served from a fake `--unitigs-cache` under
//! `--offline`: no network and no GraphAligner, because nothing in these graphs shares a
//! k-mer with the query, so every accession stops at `no_anchors`.

use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

const BIN: &str = env!("CARGO_BIN_EXE_logan-walker");
const QUERY: &str = "ATGGCTGCTAAAGAAGGCGTTCTGGAAATCGGTCAGTGGATCAAAGAGCTGCTGAAAGGCTAA";

struct Fixture {
    dir: PathBuf,
}

impl Fixture {
    /// `n` accessions cached in `<dir>/cache`, listed in `<dir>/accs.txt`.
    fn new(name: &str, n: usize) -> Self {
        let dir = std::env::temp_dir().join(format!("logan-walker-cli-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("cache")).unwrap();
        std::fs::create_dir_all(dir.join("empty_cache")).unwrap();
        std::fs::write(dir.join("q.fasta"), format!(">gene\n{QUERY}\n")).unwrap();
        let mut list = String::new();
        for i in 0..n {
            let acc = format!("SRRFAKE{i}");
            let unitigs = format!(
                ">{acc}_0 ka:f:5.0   L:+:1:+\nTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATT\n\
                 >{acc}_1 ka:f:5.0   L:-:0:-\nCCCCTTTTAAAAGGGGCCCCTTTTAAAAGGGGCC\n"
            );
            let zst = zstd::encode_all(unitigs.as_bytes(), 3).unwrap();
            std::fs::write(dir.join("cache").join(format!("{acc}.unitigs.fa.zst")), zst).unwrap();
            list.push_str(&acc);
            list.push('\n');
        }
        std::fs::write(dir.join("accs.txt"), list).unwrap();
        Fixture { dir }
    }

    fn cmd(&self, cache: &str, extra: &[&str]) -> Command {
        let mut c = Command::new(BIN);
        c.arg("-l").arg(self.dir.join("accs.txt"))
            .arg("-q").arg(self.dir.join("q.fasta"))
            .arg("--unitigs-cache").arg(self.dir.join(cache))
            .arg("--offline")
            .arg("-o").arg(self.dir.join("out"))
            .arg("-f").arg(self.dir.join("out").join("r.fasta"))
            .args(extra);
        c
    }

    fn code(&self, cache: &str, extra: &[&str]) -> (i32, String) {
        let out = self.cmd(cache, extra).output().unwrap();
        (out.status.code().unwrap_or(-1), String::from_utf8_lossy(&out.stderr).into_owned())
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

fn wait_or_kill(child: &mut std::process::Child, limit: Duration) -> Option<std::process::ExitStatus> {
    let start = Instant::now();
    while start.elapsed() < limit {
        if let Some(status) = child.try_wait().unwrap() {
            return Some(status);
        }
        std::thread::sleep(Duration::from_millis(50));
    }
    let _ = child.kill();
    let _ = child.wait();
    None
}

/// `logan-walker ... | head -1` used to panic on the closed pipe and then hang: the fetch
/// pool blocked on a full channel that nothing would drain again.
#[test]
fn a_closed_stdout_ends_the_run_instead_of_hanging() {
    use std::os::unix::process::ExitStatusExt;
    let fx = Fixture::new("pipe", 300);
    let mut child = fx.cmd("cache", &[])
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    drop(child.stdout.take()); // the reader goes away, as `head` does
    let status = wait_or_kill(&mut child, Duration::from_secs(60)).expect("walker hung on a closed stdout");
    assert_eq!(status.signal(), Some(13), "expected to end on SIGPIPE, got {status:?}");
}

#[test]
fn nothing_found_fails_unless_allow_empty() {
    let fx = Fixture::new("negatives", 3);
    let (code, err) = fx.code("cache", &[]);
    assert_eq!(code, 1, "{err}");
    assert!(err.contains("none errored"), "{err}");
    let (code, err) = fx.code("cache", &["--allow-empty"]);
    assert_eq!(code, 0, "{err}");
}

#[test]
fn every_accession_erroring_exits_3_even_with_allow_empty() {
    // An empty cache under --offline: every accession is a cache miss, an error.
    let fx = Fixture::new("errors", 3);
    for extra in [&[][..], &["--allow-empty"][..]] {
        let (code, err) = fx.code("empty_cache", extra);
        assert_eq!(code, 3, "{extra:?}: {err}");
        assert!(err.contains("all 3 accessions failed with an error"), "{err}");
    }
}

#[test]
fn zero_workers_is_a_usage_error_not_a_hang() {
    let fx = Fixture::new("workers", 6);
    let mut child = fx.cmd("cache", &["-w", "0"])
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let status = wait_or_kill(&mut child, Duration::from_secs(60)).expect("-w 0 hung");
    assert_eq!(status.code(), Some(1));
    let mut err = String::new();
    std::io::Read::read_to_string(child.stderr.as_mut().unwrap(), &mut err).unwrap();
    assert!(err.contains("--workers must be at least 1"), "{err}");
}

