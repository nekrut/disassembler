//! The walker binary end to end, on accessions served from a fake `--unitigs-cache` under
//! `--offline`, so no network. GraphAligner is a stand-in script passed as `--graphaligner`:
//! in most fixtures nothing shares a k-mer with the query, so every accession stops at
//! `no_anchors` before it is called; the `matching` ones reach it.

use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

const BIN: &str = env!("CARGO_BIN_EXE_logan-walker");
const QUERY: &str = "ATGGCTGCTAAAGAAGGCGTTCTGGAAATCGGTCAGTGGATCAAAGAGCTGCTGAAAGGCTAA";

struct Fixture {
    dir: PathBuf,
}

impl Fixture {
    /// `n` accessions cached in `<dir>/cache`, listed in `<dir>/accs.txt`, whose graphs share
    /// nothing with the query.
    fn new(name: &str, n: usize) -> Self {
        Self::build(name, n, false)
    }

    /// As `new`, but each graph contains the query, so every accession reaches GraphAligner.
    fn matching(name: &str, n: usize) -> Self {
        Self::build(name, n, true)
    }

    fn build(name: &str, n: usize, matching: bool) -> Self {
        let dir = std::env::temp_dir().join(format!("logan-walker-cli-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("cache")).unwrap();
        std::fs::create_dir_all(dir.join("empty_cache")).unwrap();
        std::fs::write(dir.join("q.fasta"), format!(">gene\n{QUERY}\n")).unwrap();
        // Stand-in GraphAligners: both answer --version; `ga_fast` then finds no alignment,
        // `ga_slow` sleeps, so a test can catch the walker mid-alignment.
        for (name, body) in [("ga_fast", "exit 0"), ("ga_slow", "sleep 5")] {
            let script = dir.join(name);
            std::fs::write(&script, format!(
                "#!/bin/sh\n[ \"$1\" = --version ] && {{ echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }}\n{body}\n"
            )).unwrap();
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).unwrap();
        }
        let mut list = String::new();
        for i in 0..n {
            let acc = format!("SRRFAKE{i}");
            let first = if matching { QUERY } else { "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATT" };
            let unitigs = format!(
                ">{acc}_0 ka:f:5.0   L:+:1:+\n{first}\n\
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
        self.cmd_ga("ga_fast", cache, extra)
    }

    fn cmd_ga(&self, ga: &str, cache: &str, extra: &[&str]) -> Command {
        let mut c = Command::new(BIN);
        c.arg("--graphaligner").arg(self.dir.join(ga))
            .arg("-l").arg(self.dir.join("accs.txt"))
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

/// A gzipped query and accession list, not named `.gz`, give the same run as plain ones.
#[test]
fn gzipped_query_and_list_are_read() {
    use std::io::Write;
    let fx = Fixture::new("gzip", 3);
    for name in ["q.fasta", "accs.txt"] {
        let path = fx.dir.join(name);
        let plain = std::fs::read(&path).unwrap();
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        e.write_all(&plain).unwrap();
        std::fs::write(&path, e.finish().unwrap()).unwrap();
    }
    let (code, err) = fx.code("cache", &["--allow-empty"]);
    assert_eq!(code, 0, "{err}");
    let summary = std::fs::read_to_string(fx.dir.join("out").join("q_walker_summary.json")).unwrap();
    assert_eq!(summary.matches("\"status\": \"no_anchors\"").count(), 3, "{summary}");
}

/// The run's scratch directories under `base`.
fn run_dirs(base: &std::path::Path) -> Vec<PathBuf> {
    std::fs::read_dir(base)
        .map(|it| {
            it.flatten()
                .map(|e| e.path())
                .filter(|p| p.file_name().unwrap().to_string_lossy().starts_with("logan-walker."))
                .collect()
        })
        .unwrap_or_default()
}

/// Wait until some accession's GFA or query file exists in the run's scratch directory.
fn wait_for_scratch_file(base: &std::path::Path) {
    let start = Instant::now();
    while start.elapsed() < Duration::from_secs(30) {
        for d in run_dirs(base) {
            if std::fs::read_dir(&d).map(|mut it| it.next().is_some()).unwrap_or(false) {
                return;
            }
        }
        std::thread::sleep(Duration::from_millis(20));
    }
    panic!("no scratch file appeared under {}", base.display());
}

/// A missing GraphAligner used to surface per accession, after each had been fetched and
/// parsed. It is now checked once, before anything is read.
#[test]
fn a_missing_graphaligner_stops_the_run_before_any_accession() {
    let fx = Fixture::matching("noga", 3);
    let out = fx.cmd("cache", &["--graphaligner", "/nonexistent/GraphAligner"]).output().unwrap();
    let (stdout, stderr) = (String::from_utf8_lossy(&out.stdout), String::from_utf8_lossy(&out.stderr));
    assert_eq!(out.status.code(), Some(1), "{stderr}");
    assert!(stderr.contains("GraphAligner not found, or not executable, at /nonexistent/GraphAligner"), "{stderr}");
    assert!(!stdout.contains("SRRFAKE"), "no accession may be processed: {stdout}");
}

/// SIGTERM mid-alignment used to leave the run's GFA and query files in /dev/shm, holding
/// RAM until reboot. The run's directory is now removed, and the run still ends on SIGTERM.
#[test]
fn sigterm_removes_the_scratch_directory() {
    use std::os::unix::process::ExitStatusExt;
    let fx = Fixture::matching("sigterm", 4);
    let base = fx.dir.join("scratch");
    let base_arg = base.to_string_lossy().into_owned();
    let mut child = fx.cmd_ga("ga_slow", "cache", &["--scratch-dir", &base_arg])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    wait_for_scratch_file(&base);
    let ok = Command::new("kill").args(["-TERM", &child.id().to_string()]).status().unwrap();
    assert!(ok.success());
    let status = wait_or_kill(&mut child, Duration::from_secs(60)).expect("walker ignored SIGTERM");
    assert_eq!(status.signal(), Some(15), "{status:?}");
    assert_eq!(run_dirs(&base), Vec::<PathBuf>::new(), "the run's scratch directory was left behind");
}

/// A run killed outright (SIGKILL, the OOM killer) cannot clean up after itself; the next run
/// over the same --scratch-dir removes what it left, and its own directory when it ends.
#[test]
fn a_killed_runs_scratch_is_swept_by_the_next_run() {
    let fx = Fixture::matching("sigkill", 4);
    let base = fx.dir.join("scratch");
    let base_arg = base.to_string_lossy().into_owned();
    let mut child = fx.cmd_ga("ga_slow", "cache", &["--scratch-dir", &base_arg])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    wait_for_scratch_file(&base);
    child.kill().unwrap(); // SIGKILL
    child.wait().unwrap();
    assert_eq!(run_dirs(&base).len(), 1, "SIGKILL leaves the directory, as expected");

    let out = fx.cmd("cache", &["--scratch-dir", &base_arg, "--allow-empty"]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "{stderr}");
    assert!(stderr.contains("removed 1 scratch directory"), "{stderr}");
    assert_eq!(run_dirs(&base), Vec::<PathBuf>::new());
}

/// `-a` alongside `-l` used to win silently, so a run meant for a list did one accession.
#[test]
fn an_accession_and_a_list_together_are_refused() {
    let fx = Fixture::new("both", 3);
    let (code, err) = fx.code("cache", &["-a", "SRRFAKE0"]);
    assert_eq!(code, 1, "{err}");
    assert!(err.contains("not both"), "{err}");
}

/// A repeated accession is processed once, so it cannot produce two `>ACC_primary` records.
#[test]
fn a_duplicated_accession_is_processed_once() {
    let fx = Fixture::new("dupes", 2);
    std::fs::write(fx.dir.join("accs.txt"), "SRRFAKE0\nSRRFAKE1\nSRRFAKE0\n").unwrap();
    let (code, err) = fx.code("cache", &["--allow-empty"]);
    assert_eq!(code, 0, "{err}");
    assert!(err.contains("1 duplicate accession(s)"), "{err}");
    let summary = std::fs::read_to_string(fx.dir.join("out").join("q_walker_summary.json")).unwrap();
    assert_eq!(summary.matches("\"acc\": \"SRRFAKE0\"").count(), 1, "{summary}");
    assert!(summary.contains("\"total\": 2"), "{summary}");
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

