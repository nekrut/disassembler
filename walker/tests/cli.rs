//! The walker binary end to end, on accessions served from a fake `--unitigs-cache` under
//! `--offline`, so no network. GraphAligner is a stand-in script passed as `--graphaligner`:
//! in most fixtures nothing shares a k-mer with the query, so every accession stops at
//! `no_anchors` before it is called; the `matching` ones reach it.

use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

const BIN: &str = env!("CARGO_BIN_EXE_logan-walker");
const QUERY: &str = "ATGGCTGCTAAAGAAGGCGTTCTGGAAATCGGTCAGTGGATCAAAGAGCTGCTGAAAGGCTAA";

/// Body of the `ga_align` stand-in. Reads the walker's own `-g`/`-f`/`-a` arguments, so it stays
/// correct if `align::command` reorders them.
const GA_ALIGN: &str = r#"#!/bin/sh
[ "$1" = --version ] && { echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }
gfa=; reads=; gaf=
while [ $# -gt 0 ]; do
  case "$1" in
    -g) gfa=$2; shift 2 ;;
    -f) reads=$2; shift 2 ;;
    -a) gaf=$2; shift 2 ;;
    *) shift ;;
  esac
done
node=$(awk -F'\t' '$1=="S" && length($3) > m { m = length($3); n = $2 } END { print n }' "$gfa")
nlen=$(awk -F'\t' '$1=="S" && length($3) > m { m = length($3) } END { print m }' "$gfa")
qname=$(awk 'sub(/^>/, "") { print $1; exit }' "$reads")
qlen=$(awk '!/^>/ { n += length($0) } END { print n + 0 }' "$reads")
printf '%s\t%s\t0\t%s\t+\t>%s\t%s\t0\t%s\t%s\t%s\t60\n' \
  "$qname" "$qlen" "$qlen" "$node" "$nlen" "$nlen" "$qlen" "$qlen" > "$gaf"
"#;

// ---- bubble fixture ----
//
// A four-unitig graph with one SNP, built to real de Bruijn geometry for the walker's default
// `graph_k: 31`: a SNP's k-mers form two parallel branch unitigs of length 2k-1 = 61 that share
// a 30 bp prefix (the source's suffix) and a 30 bp suffix (the merge's prefix).
//
//        N0 (60 nt, ka 100)
//        /              \
//   N1 (61, ka 40, ref)  N2 (61, ka 60, alt)      <- the SNP, at query offset 60
//        \              /
//        N3 (60 nt, ka 100)
//
// Stitching N0 -> N1 -> N3 with the k-1 overlap rebuilds BUBBLE_QUERY exactly; N0 -> N2 -> N3
// gives BUBBLE_MAJORITY, which differs from the query at that one base. The query has no
// repeated 21-mer (the anchor k) and no repeated 30-mer, so anchoring and the overlap
// arithmetic are unambiguous. Generated with a fixed seed; the invariants above are re-checked
// by `the_bubble_fixture_has_the_geometry_it_claims`.
//
// The minority branch is the reference-matching one on purpose: that is the reference-bias case
// `variants::prefer_majority` exists to correct (docs/issue_intra_accession_variance_lost.md).
const BUBBLE_QUERY: &str = "CACGATACAATAAGAGTAGTGATTGGCCAGTTGTGATTGTAGGAACCCATGCCTGCCTGAACCTACAACTCCTGAGACTGGAGATCCCGCACAGCCTCGATCAGACGGTGCCTAGTCCAAG";
const BUBBLE_MAJORITY: &str = "CACGATACAATAAGAGTAGTGATTGGCCAGTTGTGATTGTAGGAACCCATGCCTGCCTGAGCCTACAACTCCTGAGACTGGAGATCCCGCACAGCCTCGATCAGACGGTGCCTAGTCCAAG";
const BUBBLE_N0: &str = "CACGATACAATAAGAGTAGTGATTGGCCAGTTGTGATTGTAGGAACCCATGCCTGCCTGA";
const BUBBLE_N1: &str = "TTGTGATTGTAGGAACCCATGCCTGCCTGAACCTACAACTCCTGAGACTGGAGATCCCGCA";
const BUBBLE_N2: &str = "TTGTGATTGTAGGAACCCATGCCTGCCTGAGCCTACAACTCCTGAGACTGGAGATCCCGCA";
const BUBBLE_N3: &str = "CCTACAACTCCTGAGACTGGAGATCCCGCACAGCCTCGATCAGACGGTGCCTAGTCCAAG";

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
        // A third stand-in that actually aligns: it reads the subgraph GFA the walker hands it
        // on `-g`, picks the longest segment (in a `matching` fixture that is the one carrying
        // the query), and writes a GAF row on `-a` spanning that segment end to end. That is
        // enough for the walker to stitch, phase and write its outputs, which is what lets a
        // test reach a successful reconstruction with no GraphAligner and no network.
        //
        // The row is deliberately full-span (`path_start` 0, `path_end` == `path_len`). A row
        // whose window sits inside the path is the WV-01 reproducer and belongs with that fix,
        // not here: asserting today's output for such a row would pin the defect in place.
        {
            let script = dir.join("ga_align");
            std::fs::write(&script, GA_ALIGN).unwrap();
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

    /// The bubble fixture described above, served from `cache` under `--offline`, with a
    /// stand-in GraphAligner that replays a fixed GAF row walking the REFERENCE branch.
    fn bubble(name: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("logan-walker-cli-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("cache")).unwrap();
        std::fs::write(dir.join("q.fasta"), format!(">gene\n{BUBBLE_QUERY}\n")).unwrap();
        std::fs::write(dir.join("accs.txt"), "SRRFAKE0\n").unwrap();

        // Edges are bidirected and Logan records each one in BOTH endpoints' headers, so this
        // fixture is closed under reversal; a graph that is not is not one this pipeline can be
        // handed. Same shape as `phase::tests::detects_a_genuine_two_branch_bubble`.
        let unitigs = format!(
            ">SRRFAKE0_0 ka:f:100.0  L:+:1:+  L:+:2:+\n{BUBBLE_N0}\n\
             >SRRFAKE0_1 ka:f:40.0   L:-:0:-  L:+:3:+\n{BUBBLE_N1}\n\
             >SRRFAKE0_2 ka:f:60.0   L:-:0:-  L:+:3:+\n{BUBBLE_N2}\n\
             >SRRFAKE0_3 ka:f:100.0  L:-:1:-  L:-:2:-\n{BUBBLE_N3}\n"
        );
        let zst = zstd::encode_all(unitigs.as_bytes(), 3).unwrap();
        std::fs::write(dir.join("cache").join("SRRFAKE0.unitigs.fa.zst"), zst).unwrap();

        // The GAF row the stand-in replays. `>0>1>3` is the reference-matching branch, which is
        // the route GraphAligner takes on this graph because the query matches it base for base.
        // The window spans the whole path (`path_start` 0, `path_end` == `path_len`), so
        // `prefer_majority`'s length guard is satisfied and the swap is reached; see WV-01 for
        // what a flanked window does instead.
        let n = BUBBLE_QUERY.len();
        std::fs::write(
            dir.join("gaf.txt"),
            format!("gene\t{n}\t0\t{n}\t+\t>0>1>3\t{n}\t0\t{n}\t{n}\t{n}\t60\n"),
        )
        .unwrap();

        let script = dir.join("ga_replay");
        std::fs::write(
            &script,
            format!(
                "#!/bin/sh\n\
                 [ \"$1\" = --version ] && {{ echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }}\n\
                 gaf=\n\
                 while [ $# -gt 0 ]; do\n\
                 \x20 case \"$1\" in\n\
                 \x20   -a) gaf=$2; shift 2 ;;\n\
                 \x20   *) shift ;;\n\
                 \x20 esac\n\
                 done\n\
                 cat '{}' > \"$gaf\"\n",
                dir.join("gaf.txt").display()
            ),
        )
        .unwrap();
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).unwrap();

        Fixture { dir }
    }

    /// Overwrite the GAF row the `ga_replay` stand-in serves, so a test can choose the
    /// alignment window. Columns 3/4 are the query window, 8/9 the path window.
    fn set_gaf(&self, q_start: usize, q_end: usize, p_start: usize, p_end: usize) {
        let n = BUBBLE_QUERY.len();
        let matches = q_end - q_start;
        std::fs::write(
            self.dir.join("gaf.txt"),
            format!(
                "gene\t{n}\t{q_start}\t{q_end}\t+\t>0>1>3\t{n}\t{p_start}\t{p_end}\t{matches}\t{matches}\t60\n"
            ),
        )
        .unwrap();
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



/// A default-flags run, end to end, through to the files on disk.
///
/// Until this test there was none: every other case in this file stops at `no_anchors` or looks
/// only at an exit code, `bench/check_golden.sh` runs the walker under
/// `--prune-frac 0 --no-prefer-majority --no-emit-alts` and compares only the bare primary, and
/// the data-backed `#[ignore]` tests call stage functions directly and are skipped in CI. So
/// nothing covered the defaults every user actually runs: the prune/realign acceptance, the
/// stitch, and the writers for the reconstruction, bubbles, alts, SNV table and summary
/// (review TC-02, TW-03).
///
/// This fixture has no fork, so it gates the pipeline and the writers rather than the variant
/// logic; `prefer_majority` and the alt/VCF content need a bubble fixture, which is the rest of
/// TC-02.
#[test]
fn the_default_walk_reconstructs_and_writes_every_output() {
    let f = Fixture::matching("default-walk", 1);
    let out = f
        .cmd_ga("ga_align", "cache", &["--emit-vcf", "--emit-snv-table"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "default run should succeed; stderr:\n{stderr}");

    // The prune stage ran under its default fraction rather than being switched off.
    assert!(stderr.contains("prune: depth 5x"), "prune did not report; stderr:\n{stderr}");

    let outdir = f.dir.join("out");
    let recon = std::fs::read_to_string(outdir.join("r.fasta")).unwrap();
    assert_eq!(
        recon,
        format!(">SRRFAKE0_primary cov=1.000 pident=100.0 nodes=1 ka=5.0 vaf=1.00\n{QUERY}\n"),
        "the primary record changed"
    );

    // Each writer ran. bubbles.json and the SNV table are written even with nothing to report,
    // which is what makes their absence a real signal rather than an empty-input one.
    let bubbles = std::fs::read_to_string(outdir.join("bubbles.json")).unwrap();
    assert!(bubbles.contains("\"n_bubbles\": 0"), "bubbles.json: {bubbles}");
    let snvs = std::fs::read_to_string(outdir.join("snvs.tsv")).unwrap();
    assert_eq!(snvs, "gene\tpos\tref\talt\taccession\taf\tka_alt\tka_ref\tfork\n");
    assert!(outdir.join("alt_haplotypes.fasta").exists(), "alt_haplotypes.fasta not written");

    let summary = std::fs::read_to_string(outdir.join("q_walker_summary.json")).unwrap();
    for needle in [
        "\"logan_walker\"",
        "\"graphaligner\"",
        "\"SRRFAKE0_primary\"",
        "\"outcome\": \"nothing_to_remove\"",
        "\"subgraph_nodes\": 2",
    ] {
        assert!(summary.contains(needle), "summary missing {needle}:\n{summary}");
    }
}

/// The bubble fixture is only meaningful if its de Bruijn geometry is right, and getting the
/// `graph_k - 1` overlaps wrong would make every assertion below vacuous rather than failing.
#[test]
fn the_bubble_fixture_has_the_geometry_it_claims() {
    // The branches share the source's 30 bp suffix and the merge's 30 bp prefix.
    assert_eq!(&BUBBLE_N1[..30], &BUBBLE_N0[BUBBLE_N0.len() - 30..]);
    assert_eq!(&BUBBLE_N2[..30], &BUBBLE_N0[BUBBLE_N0.len() - 30..]);
    assert_eq!(&BUBBLE_N1[BUBBLE_N1.len() - 30..], &BUBBLE_N3[..30]);
    assert_eq!(&BUBBLE_N2[BUBBLE_N2.len() - 30..], &BUBBLE_N3[..30]);
    // 2k-1 branches differing at exactly one base: one SNP, not a longer detour.
    assert_eq!(BUBBLE_N1.len(), 61);
    assert_eq!(BUBBLE_N2.len(), 61);
    let diffs = BUBBLE_N1.bytes().zip(BUBBLE_N2.bytes()).filter(|(a, b)| a != b).count();
    assert_eq!(diffs, 1, "the branches must differ at exactly one base");
    // Stitching each route with the k-1 overlap rebuilds the two expected haplotypes.
    assert_eq!(format!("{}{}{}", BUBBLE_N0, &BUBBLE_N1[30..], &BUBBLE_N3[30..]), BUBBLE_QUERY);
    assert_eq!(format!("{}{}{}", BUBBLE_N0, &BUBBLE_N2[30..], &BUBBLE_N3[30..]), BUBBLE_MAJORITY);
}

/// The default walk over a fork, end to end: phasing, the majority correction, and the alt, VCF
/// and SNV-table writers, none of which `bench/check_golden.sh` covers (it runs with
/// `--prune-frac 0 --no-prefer-majority --no-emit-alts` and compares only the bare primary).
///
/// The route the aligner takes is the reference-matching branch at 40%, which is the
/// reference-bias case `prefer_majority` exists to correct: the primary must come back as the
/// 60% allele, with the reference allele demoted to `alt1`.
#[test]
fn a_fork_is_phased_and_the_majority_allele_becomes_the_primary() {
    let f = Fixture::bubble("fork-default");
    let out = f
        .cmd_ga("ga_replay", "cache", &["--emit-vcf", "--emit-snv-table"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "stderr:\n{stderr}");
    assert!(
        stderr.contains("consensus: 1 site(s) where the route carried the minority allele, corrected"),
        "the majority correction did not fire; stderr:\n{stderr}"
    );
    assert!(stderr.contains("1 on-path fork(s), 1 biallelic candidate(s) -> 1 recovered"), "{stderr}");

    let outdir = f.dir.join("out");
    let recon = std::fs::read_to_string(outdir.join("r.fasta")).unwrap();
    assert_eq!(
        recon.lines().nth(1).unwrap(),
        BUBBLE_MAJORITY,
        "the primary must carry the 60% allele, not the reference-matching 40% one"
    );
    assert!(recon.starts_with(">SRRFAKE0_primary cov=1.000 "), "{recon}");

    // The 40% reference allele survives as the alternative, with both frequencies recorded.
    let alts = std::fs::read_to_string(outdir.join("alt_haplotypes.fasta")).unwrap();
    assert_eq!(alts.lines().nth(1).unwrap(), BUBBLE_QUERY);
    for needle in ["vaf=0.400", "primary_vaf=0.600", "alt_kaf=40.0", "fork=0+", "alt_node=1+"] {
        assert!(alts.contains(needle), "alt record missing {needle}:\n{alts}");
    }

    let bubbles = std::fs::read_to_string(outdir.join("bubbles.json")).unwrap();
    assert!(bubbles.contains("\"n_bubbles\": 1"), "{bubbles}");
    assert!(bubbles.contains("\"fork\": 0, \"fork_strand\": \"+\""), "{bubbles}");
    assert!(bubbles.contains("\"merge\": 3"), "{bubbles}");
    assert!(bubbles.contains(r#"{"node": 1, "strand": "+", "ka": 40.0, "vaf": 0.4000}"#), "{bubbles}");
    assert!(bubbles.contains(r#"{"node": 2, "strand": "+", "ka": 60.0, "vaf": 0.6000}"#), "{bubbles}");

    // VCF and SNV table agree with each other and are in 1-based query coordinates.
    let vcf = std::fs::read_to_string(outdir.join("SRRFAKE0.vcf")).unwrap();
    assert!(vcf.contains("##contig=<ID=gene,length=121>"), "{vcf}");
    assert!(
        vcf.contains("gene\t61\t.\tA\tG\t.\tPASS\tAF=0.6000;AO=60.0;RO=40.0;FORK=0"),
        "VCF row wrong:\n{vcf}"
    );
    let snvs = std::fs::read_to_string(outdir.join("snvs.tsv")).unwrap();
    assert_eq!(
        snvs.lines().nth(1).unwrap(),
        "gene\t61\tA\tG\tSRRFAKE0\t0.6000\t60.0\t40.0\t0"
    );
}

/// WV-01, pinned. **This test asserts a known defect, not intended behaviour.**
///
/// `variants::prefer_majority` rejects its swap unless the restitched sequence is as long as
/// the whole stitched path (`variants.rs:259`), but `stitch` returns only the aligned window
/// (`stitch.rs:196`). The two agree only when `path_start == 0` and `path_end >= path_len`, so
/// the correction is silently skipped for any alignment that starts or ends inside the path.
/// The comment above that guard says it "cannot fire today"; the only real GAF row the
/// repository records (`align.rs:426`, ERR2028820) is `7..1281` of 1288, so it fires there.
///
/// The window here is a 6-base flank of the same shape. Everything else is identical to
/// `a_fork_is_phased_and_the_majority_allele_becomes_the_primary`: same graph, same
/// frequencies, same default flags. The only difference is the window, and it changes which
/// allele is published.
///
/// **When WV-01 is fixed this test must fail, and the fix is to invert it**: the primary
/// becomes the 60% allele with a `consensus:` line, matching the full-span case above. It is
/// written this way, rather than left out, so the defect cannot be fixed by accident or
/// reintroduced unnoticed.
#[test]
fn wv_01_a_flanked_alignment_window_loses_the_majority_correction() {
    let f = Fixture::bubble("fork-flanked");
    f.set_gaf(6, 115, 6, 115);
    let out = f.cmd_ga("ga_replay", "cache", &["--emit-vcf"]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "stderr:\n{stderr}");

    // The fork is still found and still reported: nothing in any output says the correction
    // was skipped, which is what makes this silent.
    assert!(stderr.contains("1 on-path fork(s), 1 biallelic candidate(s) -> 1 recovered"), "{stderr}");
    assert!(
        !stderr.contains("consensus:"),
        "WV-01 appears to be FIXED: the correction fired on a flanked window. \
         Invert this test to assert the majority allele, as the full-span case does.\n{stderr}"
    );

    let outdir = f.dir.join("out");
    let primary = std::fs::read_to_string(outdir.join("r.fasta")).unwrap();
    let seq = primary.lines().nth(1).unwrap().to_string();
    assert!(
        BUBBLE_QUERY.contains(&seq),
        "WV-01 may be FIXED: the primary is no longer the reference-matching allele.\n{primary}"
    );

    // The majority allele is published as the alternative, at a VAF above 0.5, which is the
    // outward sign of the defect: `primary_vaf` is the minority.
    let alts = std::fs::read_to_string(outdir.join("alt_haplotypes.fasta")).unwrap();
    assert!(alts.contains("vaf=0.600"), "{alts}");
    assert!(alts.contains("primary_vaf=0.400"), "{alts}");
}

/// `--fetch-only` had no test at all: nothing called `pipeline::run_fetch_only`, and no case
/// passed the flag (review TW-03). It is the path a Galaxy wrapper and the CI data step both
/// use to warm a cache before the real run, so its counts and its exit code are a contract.
#[test]
fn fetch_only_reports_its_counts_and_exit_code() {
    let f = Fixture::new("fetch-only", 2);

    // Every accession already in the cache: nothing to download, even under --offline.
    let (code, stderr) = f.code("cache", &["--fetch-only"]);
    assert_eq!(code, 0, "a full cache should satisfy --fetch-only; stderr:\n{stderr}");
    assert!(
        stderr.contains(
            "[logan-walker --fetch-only] 2 fetched, 0 not available \
             (no Logan data or over --max-size-mb), 0 failed"
        ),
        "counts line wrong:\n{stderr}"
    );

    // An offline miss is an ERROR, not a negative. The distinction matters: "not available"
    // means Logan has no data for the accession (or it is over --max-size-mb) and is a fact
    // about the accession, while this is a fact about the cache. Only the former is what
    // --allow-empty is for, which is why the next assertion holds.
    let (code, stderr) = f.code("empty_cache", &["--fetch-only"]);
    assert_eq!(code, 3, "an offline miss should exit 3; stderr:\n{stderr}");
    assert!(
        stderr.contains(
            "[logan-walker --fetch-only] 0 fetched, 0 not available \
             (no Logan data or over --max-size-mb), 2 failed"
        ),
        "counts line wrong:\n{stderr}"
    );

    // --allow-empty must not rescue an all-errored fetch. It covers an empty RESULT, not a
    // run that could not read its inputs.
    let (code, stderr) = f.code("empty_cache", &["--fetch-only", "--allow-empty"]);
    assert_eq!(code, 3, "--allow-empty must not mask an all-errored fetch; stderr:\n{stderr}");
    assert!(stderr.contains("--allow-empty does not apply to errors"), "{stderr}");
}
