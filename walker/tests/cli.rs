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
/// WA-04: the same reference with 61 more bases, so the graph can carry **two** bubbles with
/// a shared unitig between them. 182 bases; the two allele sites are at offsets 60 and 121,
/// which is further apart than k - 1 = 30, so no k-mer links them.
/// Four unitigs for the WG-03 fixture's off-route diamond: arbitrary sequence, since nothing
/// stitches through them. Only their abundances and the edges between them matter.
const OFF_N4: &str = "TGGGCGAACTTGGTCACCCCGAAGTATCTGATGAGATGATCACCGAGAGCCGGGGCGAGGA";
const OFF_N5: &str = "AGATGTACGGATACTTTCCGCACAGGGACTAGGTTAACCGCGATTTCTTATCCTGCGATAG";
const OFF_N6: &str = "CCGGCCGTGTAAACCTTTCTTAGGCATGGCAGAAAATGCAATCATATAACGGGGTTAGAAG";
const OFF_N7: &str = "GGAGCCTGTAGCATGCTGCCCGATTTCCCGTGTACCCCTGTCGCTGCGAAGTATATCCAG";
const MIX_QUERY: &str = "CACGATACAATAAGAGTAGTGATTGGCCAGTTGTGATTGTAGGAACCCATGCCTGCCTGAACCTACAACTCCTGAGACTGGAGATCCCGCACAGCCTCGATCAGACGGTGCCTAGTCCAAGTTACGGATCAGCTTAGCCATGGACTTCAGATCCGTTAACGGCATTAGCCTGAATCGGTACA";

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
        //
        // `ga_fast` writes an EMPTY GAF rather than exiting 0 with no file at all. Those are
        // different answers since WA-10: an empty GAF is "the tool ran and found nothing",
        // which is what this stand-in is for, while no file is "that was not GraphAligner"
        // and is now an error. It used to be a bare `exit 0`, which is the second.
        for (name, body) in [
            ("ga_fast", "while [ $# -gt 0 ]; do [ \"$1\" = -a ] && { : > \"$2\"; exit 0; }; shift; done\nexit 0"),
            ("ga_slow", "sleep 5"),
        ] {
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

    /// One accession whose unitigs compress to very little but build a graph many times that
    /// size: 120,000 near-identical short records, each carrying four links.
    ///
    /// Dense in *links* specifically, because adjacency is what the first parse pass holds:
    /// 480,000 packed edges plus the offsets are about 2.4 MB against a 156 KB file. This is
    /// the shape `GRAPH_PEAK_PER_COMPRESSED_BYTE` does not describe. The constant was
    /// calibrated on seven real accessions at 1.24x to 1.62x; nothing in the file format holds
    /// it there, and before the parse budget nothing stopped a file like this from being
    /// parsed until the allocator gave up.
    fn dense(name: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("logan-walker-cli-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("cache")).unwrap();
        std::fs::write(dir.join("q.fasta"), format!(">gene\n{QUERY}\n")).unwrap();
        let script = dir.join("ga_fast");
        std::fs::write(
            &script,
            "#!/bin/sh\n[ \"$1\" = --version ] && { echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }\nexit 0\n",
        )
        .unwrap();
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).unwrap();

        let mut unitigs = String::with_capacity(120_000 * 40);
        for i in 0..120_000u32 {
            // Four links each, all to ids that exist, so nothing but the budget can refuse it.
            unitigs.push_str(&format!(
                ">SRRDENSE_{i} ka:f:5.0 L:+:1:+ L:+:2:+ L:-:3:- L:+:4:-\n\
                 ACGTACGTACGTACGTACGTACGTACGTACG\n"
            ));
        }
        let zst = zstd::encode_all(unitigs.as_bytes(), 3).unwrap();
        std::fs::write(dir.join("cache").join("SRRDENSE.unitigs.fa.zst"), zst).unwrap();
        std::fs::write(dir.join("accs.txt"), "SRRDENSE\n").unwrap();
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

        write_replay_script(&dir);
        Fixture { dir }
    }

    /// WG-03: a deep gene with a shallow diamond hanging off its neighbourhood.
    ///
    /// The route is two unitigs at 15,000x with no fork of its own, so nothing about the
    /// primary, the alts or the VCF depends on the diamond. The diamond is 6x against 8x,
    /// five hops away: on `SRR8209082` that is the shape of all 3,632 bubbles the absolute
    /// 5x floor reported, none above 0.5% of the gene's depth.
    fn offpath_diamond(name: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("logan-walker-cli-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("cache")).unwrap();
        let q = &BUBBLE_QUERY[0..90];
        std::fs::write(dir.join("q.fasta"), format!(">gene\n{q}\n")).unwrap();
        std::fs::write(dir.join("accs.txt"), "SRRDEEP\n").unwrap();
        let unitigs = format!(
            ">SRRDEEP_0 ka:f:15000.0  L:+:3:+\n{}\n\
             >SRRDEEP_1 ka:f:15000.0\n{OFF_N4}\n\
             >SRRDEEP_2 ka:f:15000.0\n{OFF_N7}\n\
             >SRRDEEP_3 ka:f:15000.0  L:-:0:-  L:+:4:+\n{}\n\
             >SRRDEEP_4 ka:f:15000.0  L:-:3:-  L:+:5:+  L:+:6:+\n{OFF_N4}\n\
             >SRRDEEP_5 ka:f:6.0      L:-:4:-  L:+:7:+\n{OFF_N5}\n\
             >SRRDEEP_6 ka:f:8.0      L:-:4:-  L:+:7:+\n{OFF_N6}\n\
             >SRRDEEP_7 ka:f:15000.0  L:-:5:-  L:-:6:-\n{OFF_N7}\n",
            &BUBBLE_QUERY[0..60],
            &BUBBLE_QUERY[30..90],
        );
        let zst = zstd::encode_all(unitigs.as_bytes(), 3).unwrap();
        std::fs::write(dir.join("cache").join("SRRDEEP.unitigs.fa.zst"), zst).unwrap();
        std::fs::write(
            dir.join("gaf.txt"),
            "gene\t90\t0\t90\t+\t>0>3\t90\t0\t90\t90\t90\t60\n",
        )
        .unwrap();
        write_replay_script(&dir);
        Fixture { dir }
    }

    /// WG-07: a Logan object with no `L:` token anywhere, so the graph has no adjacency.
    ///
    /// Logan publishes these. The review found 0 links in 15,864 headers for `ERR3988483`
    /// and 0 in 21,370 for `ERR3988478`, against 9,752 of 10,070 for `SRR8209082`.
    ///
    /// `split` chooses what the accession can still do. Split across two unitigs, no walk can
    /// span the gene and there is nothing to find; in one unitig the gene needs no edge and
    /// reconstructs normally, which is why the outcome cannot be decided at parse time.
    fn unlinked(name: &str, split: bool) -> Self {
        let dir =
            std::env::temp_dir().join(format!("logan-walker-cli-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("cache")).unwrap();
        std::fs::write(dir.join("q.fasta"), format!(">gene\n{BUBBLE_QUERY}\n")).unwrap();
        std::fs::write(dir.join("accs.txt"), "SRRNOLINK\n").unwrap();
        let unitigs = if split {
            format!(
                ">SRRNOLINK_0 ka:f:50.0\n{}\n>SRRNOLINK_1 ka:f:50.0\n{}\n",
                &BUBBLE_QUERY[0..70],
                &BUBBLE_QUERY[51..121]
            )
        } else {
            format!(">SRRNOLINK_0 ka:f:50.0\n{BUBBLE_QUERY}\n")
        };
        let zst = zstd::encode_all(unitigs.as_bytes(), 3).unwrap();
        std::fs::write(dir.join("cache").join("SRRNOLINK.unitigs.fa.zst"), zst).unwrap();
        for (script, body) in [
            ("ga_fast", "while [ $# -gt 0 ]; do [ \"$1\" = -a ] && { : > \"$2\"; exit 0; }; shift; done\nexit 0".to_string()),
            ("ga_align", String::new()),
        ] {
            let path = dir.join(script);
            let text = if script == "ga_align" {
                GA_ALIGN.to_string()
            } else {
                format!(
                    "#!/bin/sh\n[ \"$1\" = --version ] && {{ echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }}\n{body}\n"
                )
            };
            std::fs::write(&path, text).unwrap();
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
        }
        Fixture { dir }
    }

    /// WA-04's shape, end to end: two bubbles joined through a unitig both strains carry.
    ///
    /// `0 -> {1,2} -> 3 -> {4,5} -> 6`, with the branch pairs at 40/60 and 45/55 of their
    /// forks. The stand-in replays `>0>1>3>4>6`, the reference-matching route, so the
    /// majority correction moves the primary onto the deeper allele at both sites and the
    /// published sequence is a combination of two sites that nothing phases.
    fn two_forks(name: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("logan-walker-cli-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("cache")).unwrap();
        std::fs::write(dir.join("q.fasta"), format!(">gene\n{MIX_QUERY}\n")).unwrap();
        std::fs::write(dir.join("accs.txt"), "SRRMIX\n").unwrap();

        let (n0, a1) = (&MIX_QUERY[0..60], &MIX_QUERY[30..91]);
        let (mid, b1, tail) = (&MIX_QUERY[61..121], &MIX_QUERY[91..152], &MIX_QUERY[122..182]);
        let a2 = flip(a1, 30);
        let b2 = flip(b1, 30);
        let unitigs = format!(
            ">SRRMIX_0 ka:f:100.0  L:+:1:+  L:+:2:+\n{n0}\n\
             >SRRMIX_1 ka:f:40.0   L:-:0:-  L:+:3:+\n{a1}\n\
             >SRRMIX_2 ka:f:60.0   L:-:0:-  L:+:3:+\n{a2}\n\
             >SRRMIX_3 ka:f:100.0  L:-:1:-  L:-:2:-  L:+:4:+  L:+:5:+\n{mid}\n\
             >SRRMIX_4 ka:f:45.0   L:-:3:-  L:+:6:+\n{b1}\n\
             >SRRMIX_5 ka:f:55.0   L:-:3:-  L:+:6:+\n{b2}\n\
             >SRRMIX_6 ka:f:100.0  L:-:4:-  L:-:5:-\n{tail}\n"
        );
        let zst = zstd::encode_all(unitigs.as_bytes(), 3).unwrap();
        std::fs::write(dir.join("cache").join("SRRMIX.unitigs.fa.zst"), zst).unwrap();

        let n = MIX_QUERY.len();
        std::fs::write(
            dir.join("gaf.txt"),
            format!("gene\t{n}\t0\t{n}\t+\t>0>1>3>4>6\t{n}\t0\t{n}\t{n}\t{n}\t60\n"),
        )
        .unwrap();
        write_replay_script(&dir);
        Fixture { dir }
    }

    /// Rewrite the cached unitig file, re-compressing it. Used to spell the same graph with
    /// Logan's `km:f:` abundance alias, or to strip a record's abundance tag entirely.
    fn retag(&self, accession: &str, edit: impl Fn(&str) -> String) {
        let path = self.dir.join("cache").join(format!("{accession}.unitigs.fa.zst"));
        let plain = zstd::decode_all(std::fs::read(&path).unwrap().as_slice()).unwrap();
        let edited = edit(std::str::from_utf8(&plain).unwrap());
        std::fs::write(&path, zstd::encode_all(edited.as_bytes(), 3).unwrap()).unwrap();
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
        // TW-08: the binary reads this from the environment, and any value turns cleanup off,
        // so a developer who exports it saw four cleanup tests fail on a correct build.
        c.env_remove("LOGAN_KEEP_SCRATCH");
        // TW-10: keep scratch inside the fixture. Without this the CLI tests use, and sweep,
        // the machine-wide /dev/shm, and leave directories there when a test kills the walker.
        c.arg("--scratch-dir").arg(self.dir.join("scratch"));
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

/// TW-05. Results are stored in the order the accessions were listed, whatever order they
/// finish in.
///
/// The unit test that claimed to cover this sorted a vector of its own and never entered the
/// pipeline, so it passed with the production sort deleted. Without that sort every output
/// follows completion order, which S3 and graph size decide: `summary.json`, the FASTA, the
/// VCFs and `snvs.tsv` all reorder run to run, and `--sieve` keeps whichever duplicate
/// finished first. Reproducibility goes silently.
///
/// This forces the first-listed accession to finish LAST by giving the stand-in aligner a
/// sleep keyed on its name, which the scratch GFA path carries.
#[test]
fn results_are_stored_in_input_order_not_completion_order() {
    // `matching` so every accession actually reaches the aligner; with the default fixture
    // they stop at `no_anchors` and the stall below never fires.
    let fx = Fixture::matching("order", 4);
    // A stand-in aligner that stalls only on SRRFAKE0, so it completes after the others.
    let slow_first = fx.dir.join("ga_order");
    std::fs::write(
        &slow_first,
        // Writes an empty GAF, as `ga_fast` does: no alignment, not a missing tool (WA-10).
        "#!/bin/sh\n\
         [ \"$1\" = --version ] && { echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }\n\
         case \"$*\" in *SRRFAKE0*) sleep 2;; esac\n\
         while [ $# -gt 0 ]; do [ \"$1\" = -a ] && { : > \"$2\"; exit 0; }; shift; done\n\
         exit 0\n",
    )
    .unwrap();
    use std::os::unix::fs::PermissionsExt;
    std::fs::set_permissions(&slow_first, std::fs::Permissions::from_mode(0o755)).unwrap();

    let out = fx
        .cmd_ga("ga_order", "cache", &["-w", "4", "--fetch-workers", "4", "--allow-empty"])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(0), "{}", String::from_utf8_lossy(&out.stderr));
    // The progress lines go to stdout, and are emitted as each accession lands, so they show
    // completion order.
    let progress = String::from_utf8_lossy(&out.stdout);
    let finished: Vec<&str> = progress
        .lines()
        .filter_map(|l| l.split(']').nth(1))
        .filter_map(|l| l.split(':').next())
        .map(str::trim)
        .filter(|s| s.starts_with("SRRFAKE"))
        .collect();
    assert_eq!(
        finished.last(),
        Some(&"SRRFAKE0"),
        "the fixture must make SRRFAKE0 finish last or it proves nothing: {finished:?}"
    );

    // ... and the stored order is the input order regardless.
    let summary = std::fs::read_to_string(fx.dir.join("out").join("q.summary.json")).unwrap();
    let order: Vec<String> = summary
        .lines()
        .filter_map(|l| l.trim().strip_prefix("\"acc\": \""))
        .filter_map(|l| l.split('"').next().map(str::to_string))
        .collect();
    assert_eq!(
        order,
        vec!["SRRFAKE0", "SRRFAKE1", "SRRFAKE2", "SRRFAKE3"],
        "summary.json must be in input order"
    );
}

/// `logan-walker ... | head -1` used to panic on the closed pipe and then hang: the fetch
/// pool blocked on a full channel that nothing would drain again.
#[test]
fn a_closed_stdout_ends_the_run_instead_of_hanging() {
    use std::os::unix::fs::PermissionsExt;
    use std::os::unix::process::ExitStatusExt;
    let fx = Fixture::matching("pipe", 4);
    let base = fx.dir.join("scratch");
    let base_arg = base.to_string_lossy().into_owned();
    let pids = fx.dir.join("aligner_pids");

    // SRRFAKE0's alignment never returns; the other three take a second and then report. It
    // is that first progress line, written to a pipe nobody is reading, that ends the run,
    // so there is always an aligner running when the signal arrives.
    let script = fx.dir.join("ga_pipe");
    std::fs::write(&script, format!(
        "#!/bin/sh\n\
         [ \"$1\" = --version ] && {{ echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }}\n\
         for a in \"$@\"; do case \"$a\" in *SRRFAKE0.*) echo $$ >> '{pids}'; exec sleep 60;; esac; done\n\
         sleep 1\n\
         exec '{ga}' \"$@\"\n",
        pids = pids.display(),
        ga = fx.dir.join("ga_align").display(),
    )).unwrap();
    std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).unwrap();

    let mut child = fx
        .cmd_ga("ga_pipe", "cache", &["--scratch-dir", &base_arg])
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    drop(child.stdout.take()); // the reader goes away, as `head` does
    let status = wait_or_kill(&mut child, Duration::from_secs(60)).expect("walker hung on a closed stdout");
    assert_eq!(status.signal(), Some(13), "expected to end on SIGPIPE, got {status:?}");

    // WF-17 and WP-26: this ending used to skip cleanup altogether. The scratch directory
    // stayed in /dev/shm holding RAM until some later run over the same base swept it, and
    // the aligners carried on writing GAFs into it.
    assert_eq!(run_dirs(&base), Vec::<PathBuf>::new(), "SIGPIPE left the scratch directory");
    let started = wait_for_aligner_pids(&pids);
    assert_all_gone(&started, "after SIGPIPE");
}

/// WF-17: the signals a scheduler uses to announce a limit clean up like SIGTERM.
///
/// SGE announces a soft wall-time limit with SIGUSR1 or SIGUSR2, LSF with SIGUSR2, and a CPU
/// ulimit or a cgroup sends SIGXCPU. The handler watched only SIGINT, SIGTERM and SIGHUP, so
/// each of these killed the walker with its directory still in `/dev/shm`. SIGQUIT is in the
/// list too but is not exercised here: its default action dumps core, and a test that writes
/// core files is its own problem.
#[test]
fn a_schedulers_limit_signal_cleans_up_like_sigterm() {
    use std::os::unix::process::ExitStatusExt;
    for (name, number) in [
        ("USR1", libc::SIGUSR1),
        ("USR2", libc::SIGUSR2),
        ("XCPU", libc::SIGXCPU),
    ] {
        let fx = Fixture::matching(&format!("sig{name}"), 4);
        let base = fx.dir.join("scratch");
        let base_arg = base.to_string_lossy().into_owned();
        let pids = fx.dir.join("aligner_pids");
        pid_recording_aligner(&fx.dir, "ga_children", &pids);

        let mut child = fx
            .cmd_ga("ga_children", "cache", &["--scratch-dir", &base_arg])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .unwrap();
        let started = wait_for_aligner_pids(&pids);

        let ok = Command::new("kill")
            .args([&format!("-{name}"), &child.id().to_string()])
            .status()
            .unwrap();
        assert!(ok.success());
        let status = wait_or_kill(&mut child, Duration::from_secs(60))
            .unwrap_or_else(|| panic!("walker ignored SIG{name}"));
        assert_eq!(status.signal(), Some(number), "SIG{name}: {status:?}");
        assert_eq!(
            run_dirs(&base),
            Vec::<PathBuf>::new(),
            "SIG{name} left the scratch directory behind"
        );
        assert_all_gone(&started, &format!("after SIG{name}"));
    }
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

/// Any accession that errors makes the run non-zero, even beside a reconstruction.
///
/// Until 0.2.0 one success was enough for exit 0, so a cohort in which some accessions
/// crashed reported a green job to whatever launched it. A 5,000-accession run on 2026-11-17
/// did exactly that: 8 accessions panicked, 156 reconstructed, exit 0. The output the run did
/// produce is still written, which is why this is its own code rather than 1 ("nothing
/// reconstructed") or 3 ("nothing worked").
#[test]
fn an_accession_that_errors_fails_the_run_that_reconstructed_the_others() {
    let f = Fixture::bubble("mixed-exit");
    // SRRFAKE0 is in the cache and reconstructs. SRRGONE is not, and under --offline a cache
    // miss is an error, not one of the negatives --allow-empty forgives.
    std::fs::write(f.dir.join("accs.txt"), "SRRFAKE0\nSRRGONE\n").unwrap();

    for extra in [&[][..], &["--allow-empty"][..]] {
        let out = f.cmd_ga("ga_replay", "cache", extra).output().unwrap();
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert_eq!(out.status.code(), Some(4), "{extra:?} stderr:\n{stderr}");
        assert!(stderr.contains("1/2 accessions reconstructed"), "{stderr}");
    }

    // Non-zero, but not empty-handed: the reconstruction is still on disk and the drop is a
    // row a workflow can branch on rather than a line in a log.
    let outdir = f.dir.join("out");
    let recon = std::fs::read_to_string(outdir.join("r.fasta")).unwrap();
    assert!(recon.contains(">SRRFAKE0_primary"), "the success must still be written:\n{recon}");
    let skipped = std::fs::read_to_string(outdir.join("q.skipped.tsv")).unwrap();
    let rows: Vec<&str> = skipped.lines().collect();
    assert_eq!(rows[0], "accession\toutcome\tbytes\tdetail");
    assert!(rows[1].starts_with("SRRGONE\terror\t"), "{skipped}");
    assert!(!skipped.contains("SRRFAKE0"), "a reconstructed accession is not a drop:\n{skipped}");
    assert_eq!(rows.len(), 2, "{skipped}");

    // The summary must state the code the process exits with. It used to derive the field
    // from `status`, which made it 0 here and 1 for a run that exits 3.
    let summary = std::fs::read_to_string(outdir.join("q.summary.json")).unwrap();
    assert!(summary.contains("\"status\": \"partial\""), "{summary}");
    assert!(summary.contains("\"exit_code\": 4"), "{summary}");
}

/// WA-01: a query that can only be read once is read once.
///
/// `README.md` says every input file may be a pipe, and `cat q.fa | logan-walker -q /dev/stdin`
/// is the form it invites. Two places reopened the path after the first read: `align::run`, on
/// every accession (fixed with WA-07 and WA-08, and covered by an `align` unit test), and
/// `main`, which read the file to count its records and then read it again to load the whole
/// set when there was more than one. On a pipe that second read returned nothing, so a
/// multi-record query failed with "no FASTA records found" naming a file that was perfectly
/// good. Both record counts go through here, since the two are one decision in `main` now.
#[test]
fn a_query_arriving_on_a_pipe_is_read_once_whatever_its_record_count() {
    let f = Fixture::matching("pipe-query", 1);

    for (what, body, multi) in [
        ("one record", format!(">gene\n{QUERY}\n"), false),
        ("two records", format!(">gene\n{QUERY}\n>other\n{QUERY}TTT\n"), true),
    ] {
        // `/dev/stdin` on a pipe: the second open sees EOF, which is the whole defect. A
        // process substitution behaves the same way and is harder to arrange from a test.
        let mut child = Command::new(BIN)
            .env_remove("LOGAN_KEEP_SCRATCH")
            .arg("--scratch-dir").arg(f.dir.join("scratch"))
            .arg("--graphaligner").arg(f.dir.join("ga_fast"))
            .arg("-l").arg(f.dir.join("accs.txt"))
            .arg("-q").arg("/dev/stdin")
            .arg("--unitigs-cache").arg(f.dir.join("cache"))
            .arg("--offline")
            .arg("-o").arg(f.dir.join("out"))
            .arg("-f").arg(f.dir.join("out").join("r.fasta"))
            // `ga_fast` writes no GAF, so the accession is `no_path` and the run is empty.
            // What is under test is the query reaching the walker at all.
            .arg("--allow-empty")
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .unwrap();
        {
            use std::io::Write;
            // Dropped at the end of this block, which is the EOF the walker reads to.
            let mut stdin = child.stdin.take().unwrap();
            stdin.write_all(body.as_bytes()).unwrap();
        }
        let out = child.wait_with_output().unwrap();
        let stderr = String::from_utf8_lossy(&out.stderr);

        assert!(!stderr.contains("no FASTA records found"),
                "{what}: the query was opened a second time:\n{stderr}");
        assert_eq!(out.status.code(), Some(0), "{what}: stderr:\n{stderr}");
        assert_eq!(stderr.contains("multi-anchor mode: 2 references"), multi,
                   "{what}: wrong pool size:\n{stderr}");
    }
}

/// WA-06: a GraphAligner call that never returns fails its own accession instead of holding a
/// compute worker for the rest of the job.
///
/// Every output is written after the run loop finishes, so before the deadline the only thing
/// that ended such a run was the scheduler's wall-time kill, and that kill discarded every
/// accession the shard had already reconstructed: no FASTA, no summary, no VCF. Here the stuck
/// accession becomes one `error` row and the other two are written, which is exit 4.
#[test]
fn a_stuck_alignment_fails_its_own_accession_and_the_rest_are_still_written() {
    let f = Fixture::matching("ga-timeout", 3);

    // Aligns like `ga_align`, except on SRRFAKE1, where it hangs. `exec` so the process the
    // walker spawned is the one that stalls: `kill` reaches our child, not its descendants.
    let stall = f.dir.join("ga_stall");
    std::fs::write(&stall, format!(
        "#!/bin/sh\n\
         [ \"$1\" = --version ] && {{ echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }}\n\
         for a in \"$@\"; do case \"$a\" in *SRRFAKE1.*) exec sleep 60;; esac; done\n\
         exec '{ga}' \"$@\"\n",
        ga = f.dir.join("ga_align").display(),
    )).unwrap();
    use std::os::unix::fs::PermissionsExt;
    std::fs::set_permissions(&stall, std::fs::Permissions::from_mode(0o755)).unwrap();

    let started = std::time::Instant::now();
    let out = f
        .cmd_ga("ga_stall", "cache", &["--graphaligner-timeout", "1"])
        .output()
        .unwrap();
    let waited = started.elapsed();
    let stderr = String::from_utf8_lossy(&out.stderr);

    assert_eq!(out.status.code(), Some(4), "stderr:\n{stderr}");
    assert!(waited < std::time::Duration::from_secs(30), "waited {waited:?} for a 1 s deadline");
    assert!(stderr.contains("2/3 accessions reconstructed"), "{stderr}");

    let outdir = f.dir.join("out");
    let recon = std::fs::read_to_string(outdir.join("r.fasta")).unwrap();
    for acc in ["SRRFAKE0", "SRRFAKE2"] {
        assert!(recon.contains(&format!(">{acc}_primary")), "{acc} was lost:\n{recon}");
    }
    assert!(!recon.contains("SRRFAKE1"), "the stuck accession has no haplotype:\n{recon}");

    let skipped = std::fs::read_to_string(outdir.join("q.skipped.tsv")).unwrap();
    assert!(skipped.contains("SRRFAKE1\terror\t"), "{skipped}");
    // The row has to say what happened and which knob changes it: the per-accession detail
    // is the only place that reaches a workflow, since the run log carries the tally only.
    assert!(skipped.contains("timed out after 1 s"), "{skipped}");
    assert!(skipped.contains("--graphaligner-timeout"), "the row must name the knob:\n{skipped}");

    let summary = std::fs::read_to_string(outdir.join("q.summary.json")).unwrap();
    assert!(summary.contains("\"status\": \"partial\""), "{summary}");
    assert!(summary.contains("\"exit_code\": 4"), "{summary}");
}

/// WF-03: an accession the cache has a recorded negative for keeps that negative offline.
///
/// `disassembler fetch` writes no object for an accession Logan has no unitigs for, because
/// there is nothing to write. The documented `--offline` run that follows then saw an
/// ordinary cache miss and reported `error`: the same list and the same flags gave
/// `no_logan_data` and exit 0 online. Since 0.2.0 any error makes the run non-zero, so a
/// shard whose accessions are all 404s fails outright rather than reporting what it found.
#[test]
fn an_offline_miss_with_a_recorded_negative_keeps_its_status() {
    let f = Fixture::matching("absent-note", 1);
    // SRRFAKE0 is cached and reconstructs; the other two are not, and the cache says why.
    std::fs::write(f.dir.join("accs.txt"), "SRRFAKE0\nSRRNODATA1\nSRRTOOBIG1\n").unwrap();
    std::fs::write(f.dir.join("cache").join("SRRNODATA1.unitigs.absent"), "no_logan_data\n").unwrap();
    std::fs::write(f.dir.join("cache").join("SRRTOOBIG1.unitigs.absent"), "too_large 9000000000\n").unwrap();

    let out = f.cmd_ga("ga_align", "cache", &[]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "a run of negatives must not fail:\n{stderr}");

    let outdir = f.dir.join("out");
    let summary = std::fs::read_to_string(outdir.join("q.summary.json")).unwrap();
    assert!(summary.contains("\"no_logan_data\": 1"), "{summary}");
    assert!(summary.contains("\"too_large\": 1"), "{summary}");
    assert!(summary.contains("\"error\": 0"), "neither is a failure of the run:\n{summary}");

    let skipped = std::fs::read_to_string(outdir.join("q.skipped.tsv")).unwrap();
    assert!(skipped.contains("SRRNODATA1\tno_logan_data\t"), "{skipped}");
    assert!(skipped.contains("SRRTOOBIG1\ttoo_large\t"), "{skipped}");

    // A note for an accession that is over a cap this run has raised past no longer applies:
    // the object really is absent, and that is a plain miss. SRRFAKE0 still reconstructs and
    // SRRNODATA1 is still a negative, so this is the one-success-one-error case, exit 4.
    let out = f
        .cmd_ga("ga_align", "cache", &["--max-size-mb", "16384"])
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(4), "a raised cap turns the note back into a miss:\n{err}");
    assert!(err.contains("cache miss under --offline"), "{err}");
}

/// WA-05: `--min-cov` is a gate on the reference span, so a deletion does not make a gene
/// that is entirely present look absent.
///
/// The stand-in reports an alignment covering the whole 63 nt reference while walking only 24
/// bases, which is what a sample with a 39 nt deletion looks like. Under the old gate that was
/// `cov` 0.38 against the default 0.85, so the accession was published as `no_path`, "the gene
/// is not here", and the run exited 1 with nothing found.
#[test]
fn a_deletion_still_reconstructs_and_the_summary_says_what_the_gate_measured() {
    use std::os::unix::fs::PermissionsExt;
    let f = Fixture::matching("ref-cov", 1);
    let ga = f.dir.join("ga_deletion");
    std::fs::write(&ga, "#!/bin/sh\n\
         [ \"$1\" = --version ] && { echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }\n\
         gfa=; reads=; gaf=\n\
         while [ $# -gt 0 ]; do case \"$1\" in\n\
           -g) gfa=$2; shift 2 ;; -f) reads=$2; shift 2 ;; -a) gaf=$2; shift 2 ;; *) shift ;;\n\
         esac; done\n\
         node=$(awk -F'\\t' '$1==\"S\" && length($3) > m { m = length($3); n = $2 } END { print n }' \"$gfa\")\n\
         nlen=$(awk -F'\\t' '$1==\"S\" && length($3) > m { m = length($3) } END { print m }' \"$gfa\")\n\
         qname=$(awk 'sub(/^>/, \"\") { print $1; exit }' \"$reads\")\n\
         qlen=$(awk '!/^>/ { n += length($0) } END { print n + 0 }' \"$reads\")\n\
         # Full reference span, short walk: the deletion case.\n\
         printf '%s\\t%s\\t0\\t%s\\t+\\t>%s\\t%s\\t0\\t24\\t%s\\t%s\\t60\\n' \\\n\
           \"$qname\" \"$qlen\" \"$qlen\" \"$node\" \"$nlen\" \"$qlen\" \"$qlen\" > \"$gaf\"\n")
        .unwrap();
    std::fs::set_permissions(&ga, std::fs::Permissions::from_mode(0o755)).unwrap();

    let out = f.cmd_ga("ga_deletion", "cache", &[]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "the gene is present; stderr:\n{stderr}");

    let outdir = f.dir.join("out");
    let recon = std::fs::read_to_string(outdir.join("r.fasta")).unwrap();
    assert!(recon.contains(">SRRFAKE0_primary"), "{recon}");
    let seq: String = recon.lines().filter(|l| !l.starts_with('>')).collect();
    assert_eq!(seq.len(), 24, "the walked haplotype is the short one");

    // Both numbers are published, because they answer different questions and the gate uses
    // only one of them. Without `ref_cov` a reader sees cov 0.71 beside --min-cov 0.85 and
    // cannot tell why the record exists.
    let summary = std::fs::read_to_string(outdir.join("q.summary.json")).unwrap();
    assert!(summary.contains("\"cov\": 0.3810"), "{summary}");
    assert!(summary.contains("\"ref_cov\": 1.0000"), "{summary}");
    assert!(summary.contains("\"ref_span\": [0, 63]"), "{summary}");
}

/// Write a stand-in that replays one GAF row on the first pass and another on the second.
///
/// `ga_replay` serves the same row to both, which is right for most tests and useless for the
/// ones about what the second pass returns.
fn two_pass_aligner(dir: &std::path::Path, name: &str, pass2: &str) -> PathBuf {
    use std::os::unix::fs::PermissionsExt;
    let script = dir.join(name);
    std::fs::write(&script, format!(
        "#!/bin/sh\n\
         [ \"$1\" = --version ] && {{ echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }}\n\
         gaf=\n\
         while [ $# -gt 0 ]; do case \"$1\" in -a) gaf=$2; shift 2 ;; *) shift ;; esac; done\n\
         case \"$gaf\" in\n\
           *.pass2.gaf) {pass2} ;;\n\
           *) cat '{first}' > \"$gaf\" ;;\n\
         esac\n",
        first = dir.join("gaf.txt").display(),
    )).unwrap();
    std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).unwrap();
    script
}

/// WP-04: a second pass that fails leaves the first pass's haplotype standing.
///
/// Every step of the prune's realignment used to be `?`, so GraphAligner aborting or being
/// OOM-killed on the pruned graph turned an accession that was already reconstructed into
/// `error`: missing from the FASTA, `prune: null` in the summary, and exit 3 for a
/// one-accession run. That contradicts the stage's own contract, that pruning can improve an
/// accession and cannot break one.
#[test]
fn a_prune_whose_realignment_fails_keeps_the_first_pass() {
    let f = Fixture::bubble("prune-error");
    two_pass_aligner(&f.dir, "ga_pass2_dies", "exit 134");

    let out = f
        .cmd_ga("ga_pass2_dies", "cache", &["--prune-frac", "0.5"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "the first pass is still good:\n{stderr}");
    assert!(stderr.contains("prune: realignment failed"), "{stderr}");

    let outdir = f.dir.join("out");
    let recon = std::fs::read_to_string(outdir.join("r.fasta")).unwrap();
    assert!(recon.contains(">SRRFAKE0_primary"), "the haplotype must survive:\n{recon}");
    let summary = std::fs::read_to_string(outdir.join("q.summary.json")).unwrap();
    assert!(summary.contains("\"outcome\": \"refused_error\""), "{summary}");
    // A refused prune means the first-pass route stands, reference pull and all (WP-01).
    assert!(summary.contains("\"bias_uncorrected\": true"), "{summary}");
}

/// WA-09: a realignment that is the first pass minus an end base is refused, not "applied".
///
/// `pruned_covers_as_much` allows three reference bases of slack, because a strict rule kept a
/// reference-biased primary over one terminal codon. But the slack was the whole test, so a
/// realignment that changed nothing else was accepted and logged as `applied`: pure loss, and
/// a lost terminal base can be a truncated stop codon.
///
/// The floor here is above the off-route branch and above the route's thin unitig, so the
/// route survives the prune only because WP-01 protects it, which is what lets the second
/// pass walk the same route and return a trimmed copy of the first.
#[test]
fn a_prune_that_only_trims_an_end_base_is_refused() {
    let f = Fixture::bubble("prune-trim");
    let n = BUBBLE_QUERY.len();
    // Same route, one base shorter.
    let row = format!(
        "printf 'gene\\t{n}\\t0\\t{m}\\t+\\t>0>1>3\\t{n}\\t0\\t{m}\\t{m}\\t{m}\\t60\\n' > \"$gaf\"",
        m = n - 1
    );
    two_pass_aligner(&f.dir, "ga_pass2_trims", &row);

    let out = f
        .cmd_ga("ga_pass2_trims", "cache", &["--prune-frac", "0.7"])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "{stderr}");
    assert!(stderr.contains("otherwise identical"), "{stderr}");

    let outdir = f.dir.join("out");
    let summary = std::fs::read_to_string(outdir.join("q.summary.json")).unwrap();
    assert!(summary.contains("\"outcome\": \"refused_trimmed\""), "{summary}");
    let recon = std::fs::read_to_string(outdir.join("r.fasta")).unwrap();
    let seq: String = recon.lines().filter(|l| !l.starts_with('>')).collect();
    assert_eq!(seq.len(), n, "the full-length first pass is what gets published");
}

/// WA-10: a program that exits 0 without writing a GAF is not a biological negative.
///
/// GraphAligner 1.0.17b creates the GAF before it aligns anything and aborts if it cannot, so
/// this state means the thing on `--graphaligner` is not GraphAligner, or a wrapper swallowed
/// its failure. The startup `--version` probe accepts any program that prints a line and exits
/// 0, which `/usr/bin/true` does; the run then reported `no_path` for every accession and, with
/// `--allow-empty`, exited 0, which reads as "the gene is in none of them".
#[test]
fn an_aligner_that_writes_no_gaf_is_an_error_not_a_negative() {
    use std::os::unix::fs::PermissionsExt;
    let f = Fixture::matching("no-gaf", 2);
    let ga = f.dir.join("ga_silent");
    std::fs::write(&ga, "#!/bin/sh\n\
         [ \"$1\" = --version ] && { echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }\n\
         exit 0\n").unwrap();
    std::fs::set_permissions(&ga, std::fs::Permissions::from_mode(0o755)).unwrap();

    // `--allow-empty` forgives "nothing was found". This is not that.
    let out = f.cmd_ga("ga_silent", "cache", &["--allow-empty"]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(3), "{stderr}");

    let skipped = std::fs::read_to_string(f.dir.join("out").join("q.skipped.tsv")).unwrap();
    assert!(skipped.contains("exited 0 but wrote no GAF"), "{skipped}");
    assert!(!skipped.contains("no_path"), "this is not a negative result:\n{skipped}");

    // An EMPTY GAF still is one: the tool ran, wrote its file and found nothing.
    let ok = f.cmd_ga("ga_fast", "cache", &["--allow-empty"]).output().unwrap();
    assert_eq!(ok.status.code(), Some(0), "{}", String::from_utf8_lossy(&ok.stderr));
    let s2 = std::fs::read_to_string(f.dir.join("out").join("q.summary.json")).unwrap();
    assert!(s2.contains("\"no_path\": 2"), "{s2}");
}

/// `--revalidate` is refused where it cannot do anything, rather than silently ignored.
///
/// A flag that quietly does nothing is worse than one that is rejected: the summary records
/// `revalidate: true` either way, so a run that could not revalidate would look like one
/// that did and found everything current.
#[test]
fn revalidate_requires_a_cache_and_refuses_offline() {
    let fx = Fixture::new("reval-usage", 1);

    // Nothing to revalidate without a cache.
    let out = Command::new(BIN)
        .arg("--revalidate")
        .arg("-l")
        .arg(fx.dir.join("accs.txt"))
        .arg("-q")
        .arg(fx.dir.join("q.fasta"))
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(1), "{err}");
    assert!(err.contains("--revalidate requires --unitigs-cache"), "{err}");

    // `Fixture::cmd` always passes --offline, which forbids the network this needs.
    let (code, err) = fx.code("cache", &["--revalidate"]);
    assert_eq!(code, 1, "{err}");
    assert!(err.contains("mutually exclusive"), "{err}");
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
    let summary = std::fs::read_to_string(fx.dir.join("out").join("q.summary.json")).unwrap();
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

/// Whether `pid` names a process that is still running.
///
/// A reaped child's PID is gone, but an orphan whose reaper has not got to it yet is a zombie,
/// and `kill(pid, 0)` succeeds for a zombie. That would make the tests below pass or fail on
/// the host's init rather than on the walker, so on Linux the state is read directly.
fn process_is_alive(pid: u32) -> bool {
    if let Ok(stat) = std::fs::read_to_string(format!("/proc/{pid}/stat")) {
        // "pid (comm) state ...", and `comm` may itself hold spaces and brackets.
        return match stat.rsplit_once(") ") {
            Some((_, rest)) => !rest.starts_with('Z'),
            None => true,
        };
    }
    if std::path::Path::new("/proc").is_dir() {
        return false; // Linux, and the process is gone
    }
    // SAFETY: signal 0 only checks for the process; `kill` has no other preconditions.
    unsafe { libc::kill(pid as i32, 0) == 0 }
}

/// Write a stand-in aligner that records its PID and then hangs, and return its path.
fn pid_recording_aligner(dir: &std::path::Path, name: &str, pids: &std::path::Path) -> PathBuf {
    use std::os::unix::fs::PermissionsExt;
    let script = dir.join(name);
    std::fs::write(&script, format!(
        "#!/bin/sh\n\
         [ \"$1\" = --version ] && {{ echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }}\n\
         echo $$ >> '{pids}'\n\
         exec sleep 60\n",
        pids = pids.display(),
    )).unwrap();
    std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).unwrap();
    script
}

/// Every PID recorded in `pids`, once at least one is there.
fn wait_for_aligner_pids(pids: &std::path::Path) -> Vec<u32> {
    let start = Instant::now();
    while start.elapsed() < Duration::from_secs(30) {
        if let Ok(text) = std::fs::read_to_string(pids) {
            let found: Vec<u32> = text.lines().filter_map(|l| l.trim().parse().ok()).collect();
            if !found.is_empty() {
                return found;
            }
        }
        std::thread::sleep(Duration::from_millis(20));
    }
    panic!("no aligner recorded a PID in {}", pids.display());
}

/// Assert every PID is gone, allowing a moment for the signal to land.
fn assert_all_gone(pids: &[u32], how: &str) {
    let start = Instant::now();
    while start.elapsed() < Duration::from_secs(10) {
        if pids.iter().all(|&p| !process_is_alive(p)) {
            return;
        }
        std::thread::sleep(Duration::from_millis(20));
    }
    let left: Vec<u32> = pids.iter().copied().filter(|&p| process_is_alive(p)).collect();
    panic!("{how}: aligners still running: {left:?}");
}

/// WF-15: the signal that ends the run ends the aligners it started.
///
/// The handler used to remove the scratch directory and exit, and nothing tracked the
/// children, so after `kill <walker pid>` or `timeout N logan-walker` up to `-w` GraphAligner
/// processes kept burning CPU and RAM, reparented to init, working on files that had just been
/// deleted. Ctrl-C hid it, because a terminal signals the whole foreground process group.
#[test]
fn sigterm_stops_the_aligners_the_run_started() {
    use std::os::unix::process::ExitStatusExt;
    let fx = Fixture::matching("ga-children", 4);
    let base = fx.dir.join("scratch");
    let base_arg = base.to_string_lossy().into_owned();
    let pids = fx.dir.join("aligner_pids");
    pid_recording_aligner(&fx.dir, "ga_children", &pids);

    let mut child = fx
        .cmd_ga("ga_children", "cache", &["--scratch-dir", &base_arg])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    let started = wait_for_aligner_pids(&pids);

    let ok = Command::new("kill").args(["-TERM", &child.id().to_string()]).status().unwrap();
    assert!(ok.success());
    let status = wait_or_kill(&mut child, Duration::from_secs(60)).expect("walker ignored SIGTERM");
    assert_eq!(status.signal(), Some(15), "{status:?}");

    assert_all_gone(&started, "after SIGTERM to the walker");
    assert_eq!(run_dirs(&base), Vec::<PathBuf>::new(), "the scratch directory went too");
}

/// SIGKILL never runs the handler, so the kernel is asked to do it: the aligner is spawned
/// with `PR_SET_PDEATHSIG`, which fires when the thread that spawned it dies. Linux only.
#[cfg(target_os = "linux")]
#[test]
fn sigkill_of_the_walker_still_takes_the_aligners_with_it() {
    let fx = Fixture::matching("ga-pdeath", 4);
    let base = fx.dir.join("scratch");
    let base_arg = base.to_string_lossy().into_owned();
    let pids = fx.dir.join("aligner_pids");
    pid_recording_aligner(&fx.dir, "ga_children", &pids);

    let mut child = fx
        .cmd_ga("ga_children", "cache", &["--scratch-dir", &base_arg])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();
    let started = wait_for_aligner_pids(&pids);

    child.kill().unwrap(); // SIGKILL: no handler runs, and the scratch directory stays
    child.wait().unwrap();

    assert_all_gone(&started, "after SIGKILL of the walker");
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
    let summary = std::fs::read_to_string(fx.dir.join("out").join("q.summary.json")).unwrap();
    assert_eq!(summary.matches("\"acc\": \"SRRFAKE0\"").count(), 1, "{summary}");
    // The duplicate was dropped with a warning and nothing else. The summary now accounts
    // for it: three entries were listed, two accessions ran, and the gap is a discard.
    assert!(summary.contains("\"n_requested\": 3"), "{summary}");
    assert!(summary.contains("\"n_records\": 2"), "{summary}");
    assert!(summary.contains("\"never_reached\": 1"), "{summary}");
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
        // WA-04 added `min_vaf=`, the measured counterpart of the constant `vaf=`. This
        // route passes no fork at all, so it is 1.000 and there is no `mixed=` tag.
        format!(
            ">SRRFAKE0_primary cov=1.000 pident=100.0 nodes=1 ka=5.0 vaf=1.00 min_vaf=1.000\n{QUERY}\n"
        ),
        "the primary record changed"
    );

    // Each writer ran. bubbles.json and the SNV table are written even with nothing to report,
    // which is what makes their absence a real signal rather than an empty-input one.
    let bubbles = std::fs::read_to_string(outdir.join("bubbles.json")).unwrap();
    assert!(bubbles.contains("\"n_bubbles\": 0"), "bubbles.json: {bubbles}");
    let snvs = std::fs::read_to_string(outdir.join("snvs.tsv")).unwrap();
    assert_eq!(
        snvs,
        "gene\tpos\tref\talt\taccession\taf\tka_alt\tka_ref\tfork\tka_saturated\tblock\n"
    );
    assert!(outdir.join("alt_haplotypes.fasta").exists(), "alt_haplotypes.fasta not written");

    let summary = std::fs::read_to_string(outdir.join("q.summary.json")).unwrap();
    for needle in [
        "\"schema\": 2",
        "\"name\": \"logan-walker\"",
        "\"GraphAligner\"",
        "\"SRRFAKE0_primary\"",
        "\"seq_sha256\"",
        "\"outcome\": \"nothing_to_remove\"",
        "\"subgraph_nodes\": 2",
        "\"status\": \"ok\"",
        "\"reconstructed\": 1",
        "\"distinct_primaries\": 1",
        // The accession's record names the bytes it was built from. These came from the
        // fixture cache under --offline, which has no HTTP response to take an ETag from,
        // so the digest is the whole of the provenance and it has to be there.
        "\"payload\": {\"source\": \"cache\", \"bytes\": ",
    ] {
        assert!(summary.contains(needle), "summary missing {needle}:\n{summary}");
    }
    // Threaded from the fetch stage through the compute pool onto the finished record, not
    // recomputed by the writer: a 64-hex digest of the cached object.
    let sha = summary
        .split("\"payload\": {")
        .nth(1)
        .and_then(|s| s.split("\"sha256\": \"").nth(1))
        .and_then(|s| s.split('"').next())
        .unwrap_or("");
    assert_eq!(sha.len(), 64, "payload sha256 missing or malformed:\n{summary}");
    assert!(sha.bytes().all(|b| b.is_ascii_hexdigit()), "not hex: {sha}");
    let cached = std::fs::read(f.dir.join("cache").join("SRRFAKE0.unitigs.fa.zst")).unwrap();
    assert!(
        summary.contains(&format!("\"bytes\": {}", cached.len())),
        "payload bytes should be the cache file's size ({}):\n{summary}",
        cached.len()
    );
    // A cache hit has no response behind it, so these must be absent rather than invented.
    assert!(!summary.contains("\"etag\""), "no ETag exists for a cache hit:\n{summary}");
    assert!(!summary.contains("\"fetched_at\""), "cached bytes predate the run:\n{summary}");
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
        // Trailing 0: this fork is at 40x and 60x, nowhere near Logan's 25,000 ceiling,
        // so its frequencies are real ratios rather than lower bounds (WV-03).
        // Trailing empty `block`: one substitution, so it is not part of a linked block.
        "gene\t61\tA\tG\tSRRFAKE0\t0.6000\t60.0\t40.0\t0\t0\t"
    );
}

/// One base changed, so a branch spells a different allele than its sibling.
fn flip(s: &str, i: usize) -> String {
    let mut v = s.as_bytes().to_vec();
    v[i] = if v[i] == b'A' { b'G' } else { b'A' };
    String::from_utf8(v).unwrap()
}

/// The stand-in GraphAligner shared by the bubble fixtures: it replays `gaf.txt` into
/// whatever `-a` names, so a test chooses the route by writing that file.
fn write_replay_script(dir: &std::path::Path) {
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
}

/// WG-03: `bubbles.json` is gated against the gene's own depth, not an absolute 5x.
///
/// The gate was `ka:f >= --min-abund`, an absolute 5x, applied to the whole unpruned
/// neighbourhood before the alignment even ran. On a 15,000x sample that admits anything two
/// erroneous reads agree on, anywhere near the gene: `SRR8209082` reported 3,632 bubbles whose
/// strongest allele was 66x, median 11.9x, while the 3,600x second population it does have was
/// not among them. Read as quasispecies evidence, which is what the file is for, that is noise.
///
/// Same graph, same alignment, same primary in all three runs here. Two things changed, and
/// the runs separate them: the stage now runs AFTER the prune, on the graph the prune left,
/// and its floor is depth-relative rather than an absolute 5x.
#[test]
fn a_shallow_diamond_beside_a_deep_gene_is_not_reported_as_a_bubble() {
    let bubbles_of = |f: &Fixture, flags: &[&str]| -> String {
        let out = f.cmd_ga("ga_replay", "cache", flags).output().unwrap();
        assert_eq!(
            out.status.code(),
            Some(0),
            "stderr:\n{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let recon = std::fs::read_to_string(f.dir.join("out").join("r.fasta")).unwrap();
        assert_eq!(recon.lines().nth(1).unwrap(), &BUBBLE_QUERY[0..90], "the primary moved");
        std::fs::read_to_string(f.dir.join("out").join("bubbles.json")).unwrap()
    };

    // A default run prunes at 1% of 15,000x, so the 6x and 8x unitigs are not in the graph
    // the stage sees at all. This is the ordering half of the fix: the stage used to run
    // before the alignment, on the neighbourhood as extracted.
    let f = Fixture::offpath_diamond("deep");
    assert!(
        bubbles_of(&f, &[]).contains("\"n_bubbles\": 0"),
        "a pruned graph has nothing to report here"
    );

    // With the prune off, the graph still holds the diamond, and the absolute 5x floor
    // admits it. This is what `bubbles.json` held for every sample, however deep.
    let g = Fixture::offpath_diamond("deep-unpruned");
    let kept = bubbles_of(&g, &["--prune-frac", "0"]);
    assert!(kept.contains("\"n_bubbles\": 1"), "{kept}");
    assert!(kept.contains("\"ka\": 6.0"), "{kept}");
    // The record says where that allele set was found (WG-04).
    assert!(kept.contains("\"flanks\": [{\"fork\": 4"), "{kept}");

    // Same unpruned graph, same alignment, a depth-relative floor: this is the floor half,
    // and it is what covers a refused prune, where the graph keeps its sub-floor branches.
    let h = Fixture::offpath_diamond("deep-relative");
    assert!(
        bubbles_of(&h, &["--prune-frac", "0", "--min-abund-frac", "0.01"])
            .contains("\"n_bubbles\": 0"),
        "1% of 15,000x is 150x, and the branches are 6x and 8x"
    );
}

/// WP-06: a run that is killed leaves the reconstructions it had already made.
///
/// Everything a run produces was written after the last accession finished, so a scheduler's
/// walltime SIGTERM, an OOM kill or a dropped connection at hour eight of a nine-hour shard
/// threw away every reconstruction and left nothing to resume from. It has happened here
/// already: a `timeout 1500` killed a 50-accession measurement run at 47 of 50 and wrote
/// nothing at all.
///
/// Four accessions, an aligner that takes a second each, one compute worker: SIGTERM after
/// the first ones land.
#[test]
fn a_killed_run_keeps_the_haplotypes_it_had_already_reconstructed() {
    use std::os::unix::fs::PermissionsExt;
    let fx = Fixture::matching("journal", 4);
    // An aligner that is slow enough that the kill lands mid-cohort.
    let script = fx.dir.join("ga_slowish");
    std::fs::write(&script, format!(
        "#!/bin/sh\n\
         [ \"$1\" = --version ] && {{ echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }}\n\
         sleep 1\n\
         exec '{ga}' \"$@\"\n",
        ga = fx.dir.join("ga_align").display(),
    )).unwrap();
    std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).unwrap();

    let outdir = fx.dir.join("out");
    let partial = outdir.join("q.partial.fasta");
    let journal = outdir.join("q.journal.jsonl");
    let mut child = fx
        .cmd_ga("ga_slowish", "cache", &["-w", "1"])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .unwrap();

    // Wait for the first accession to be journalled, then kill the run mid-cohort.
    let start = Instant::now();
    while start.elapsed() < Duration::from_secs(30) {
        if std::fs::read_to_string(&journal).map(|t| !t.is_empty()).unwrap_or(false) {
            break;
        }
        std::thread::sleep(Duration::from_millis(20));
    }
    let ok = Command::new("kill").args(["-TERM", &child.id().to_string()]).status().unwrap();
    assert!(ok.success());
    wait_or_kill(&mut child, Duration::from_secs(60)).expect("walker ignored SIGTERM");

    // The canonical outputs belong to a run that finished, and this one did not.
    assert!(!outdir.join("q.summary.json").exists(), "a summary from a killed run");
    let recon = std::fs::read_to_string(fx.dir.join("out").join("r.fasta")).unwrap_or_default();
    assert!(!recon.contains(">SRRFAKE"), "the final FASTA was written by a killed run");

    // The salvage is there, and it is the real thing: FASTA records byte-identical to the
    // ones a finished run writes, and one journal line per accession that got through.
    let saved = std::fs::read_to_string(&partial).expect("no partial FASTA");
    let lines = std::fs::read_to_string(&journal).expect("no journal");
    let n_journal = lines.lines().filter(|l| !l.trim().is_empty()).count();
    assert!(n_journal >= 1, "journal is empty");
    assert!(n_journal < 4, "the run was supposed to be killed mid-cohort, not to finish");
    assert_eq!(saved.matches(">SRRFAKE").count(), n_journal, "{saved}");
    assert!(saved.contains("_primary cov=1.000 "), "{saved}");
    for line in lines.lines() {
        assert!(line.contains("\"status\": \"success\""), "{line}");
        assert!(line.contains("\"seq_sha256\": \""), "{line}");
    }
}

/// The other half: a run that finishes takes its salvage copies away again, so they only
/// ever mean "a run is in flight, or one died here".
#[test]
fn a_finished_run_leaves_no_journal_behind() {
    let f = Fixture::matching("journal-clean", 2);
    let out = f.cmd_ga("ga_align", "cache", &[]).output().unwrap();
    assert_eq!(out.status.code(), Some(0), "{}", String::from_utf8_lossy(&out.stderr));
    let outdir = f.dir.join("out");
    assert!(outdir.join("q.summary.json").is_file());
    assert!(!outdir.join("q.partial.fasta").exists(), "partial FASTA survived a clean run");
    assert!(!outdir.join("q.journal.jsonl").exists(), "journal survived a clean run");
    // And what it held while the run was going is what the real FASTA holds now.
    let recon = std::fs::read_to_string(outdir.join("r.fasta")).unwrap();
    assert_eq!(recon.matches(">SRRFAKE").count(), 2, "{recon}");
}

/// WP-06, the cheap half: an unusable output location is found before the cohort runs, not
/// after every accession has been fetched and aligned.
#[test]
fn an_unusable_out_dir_fails_before_any_accession_is_fetched() {
    let f = Fixture::matching("badout", 2);
    // A regular file where --out-dir should be.
    let blocked = f.dir.join("afile");
    std::fs::write(&blocked, "not a directory\n").unwrap();
    let marker = f.dir.join("aligner_ran");
    let script = f.dir.join("ga_marker");
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::write(&script, format!(
            "#!/bin/sh\n\
             [ \"$1\" = --version ] && {{ echo 'GraphAligner bioconda 1.0.17b-'; exit 0; }}\n\
             touch '{}'\n\
             exec '{}' \"$@\"\n",
            marker.display(),
            f.dir.join("ga_align").display(),
        )).unwrap();
        std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).unwrap();
    }
    let mut c = Command::new(BIN);
    c.env_remove("LOGAN_KEEP_SCRATCH");
    c.arg("--scratch-dir").arg(f.dir.join("scratch"))
        .arg("--graphaligner").arg(f.dir.join("ga_marker"))
        .arg("-l").arg(f.dir.join("accs.txt"))
        .arg("-q").arg(f.dir.join("q.fasta"))
        .arg("--unitigs-cache").arg(f.dir.join("cache"))
        .arg("--offline")
        .arg("-o").arg(&blocked)
        .arg("-f").arg(f.dir.join("out").join("r.fasta"));
    let out = c.output().unwrap();
    assert_ne!(out.status.code(), Some(0));
    assert!(!marker.exists(), "the cohort ran before the output location was checked");
}

/// WF-12: a lowercase accession is the accession it spells, not a missing one.
///
/// Logan's bucket is case-sensitive, so `srr8209082` is fetched verbatim, 404s, and the 404 is
/// recorded as `no_logan_data`: "Logan has nothing for this run". For a run Logan does have
/// that is a false negative nothing in the output distinguishes from a real one, and
/// `--allow-empty` turns it into exit 0 with an empty FASTA. The fixture here is a cache hit,
/// which is the same lookup by name without needing the network.
#[test]
fn a_lowercase_accession_finds_the_run_it_names() {
    let dir = std::env::temp_dir().join(format!("logan-walker-cli-{}-lowercase", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(dir.join("cache")).unwrap();
    std::fs::write(dir.join("q.fasta"), format!(">gene\n{BUBBLE_QUERY}\n")).unwrap();
    // The list is lowercase; the cached object, as on the bucket, is not.
    std::fs::write(dir.join("accs.txt"), "srr9999001\n").unwrap();
    let unitigs = format!(">SRR9999001_0 ka:f:50.0\n{BUBBLE_QUERY}\n");
    std::fs::write(
        dir.join("cache").join("SRR9999001.unitigs.fa.zst"),
        zstd::encode_all(unitigs.as_bytes(), 3).unwrap(),
    )
    .unwrap();
    std::fs::write(dir.join("ga_align"), GA_ALIGN).unwrap();
    use std::os::unix::fs::PermissionsExt;
    std::fs::set_permissions(dir.join("ga_align"), std::fs::Permissions::from_mode(0o755)).unwrap();
    let f = Fixture { dir };

    let out = f.cmd_ga("ga_align", "cache", &[]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "stderr:\n{stderr}");
    assert!(
        stderr.contains("were not uppercase"),
        "the run must say it changed the spelling; stderr:\n{stderr}"
    );
    let recon = std::fs::read_to_string(f.dir.join("out").join("r.fasta")).unwrap();
    assert!(
        recon.starts_with(">SRR9999001_primary"),
        "the output names the accession in its canonical form:\n{recon}"
    );
    let summary = std::fs::read_to_string(f.dir.join("out").join("q.summary.json")).unwrap();
    assert!(summary.contains(r#""acc": "SRR9999001""#), "{summary}");
    assert!(
        summary.contains(r#""status": "success""#)
            && !summary.contains(r#""status": "no_logan_data""#),
        "a run Logan holds must not be reported as absent:\n{summary}"
    );
}

/// WA-12: a GAF row that is not a walk of the graph is an error, not a primary.
///
/// The walker writes the GFA, hands it to GraphAligner and then believes the row that comes
/// back. `stitch_steps` appends each step's sequence minus the k - 1 overlap whatever the
/// step is, so a row naming two branches of the same fork in series stitches cleanly into a
/// sequence that exists nowhere, and publishes it with the row's own `pident`. The review
/// reproduced a primary with 89 differences from the truth at `pident=100.0`.
///
/// This fixture's graph is `0 -> {1,2} -> 3`. The stand-in is told to emit `>0>1>2`: every
/// segment is real and the lengths work out, but 1 and 2 are siblings and nothing joins them.
#[test]
fn a_gaf_row_that_is_not_a_walk_of_the_graph_fails_the_accession() {
    let f = Fixture::bubble("bad-walk");
    let n = BUBBLE_QUERY.len();
    std::fs::write(
        f.dir.join("gaf.txt"),
        format!("gene\t{n}\t0\t{n}\t+\t>0>1>2\t{n}\t0\t{n}\t{n}\t{n}\t60\n"),
    )
    .unwrap();
    let out = f.cmd_ga("ga_replay", "cache", &[]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    // One accession, and it errored: the run did not work (exit 3), rather than reporting a
    // sequence nobody can check.
    assert_eq!(out.status.code(), Some(3), "stderr:\n{stderr}");

    let outdir = f.dir.join("out");
    let summary = std::fs::read_to_string(outdir.join("q.summary.json")).unwrap();
    assert!(summary.contains(r#""status": "error""#), "{summary}");
    assert!(
        summary.contains("not joined by an edge in the graph"),
        "the record must say what was wrong with the row:\n{summary}"
    );
    // Nothing was published from it.
    let recon = std::fs::read_to_string(outdir.join("r.fasta")).unwrap_or_default();
    assert!(!recon.contains(">SRRFAKE0_primary"), "a chimera reached the FASTA:\n{recon}");

    // The same fixture with the real walk still reconstructs, so the check is about the row
    // and not about this graph.
    let g = Fixture::bubble("good-walk");
    let out = g.cmd_ga("ga_replay", "cache", &[]).output().unwrap();
    assert_eq!(out.status.code(), Some(0), "{}", String::from_utf8_lossy(&out.stderr));
}

/// WG-07: an object with no adjacency is `no_links`, not "the gene is absent".
///
/// `parse` builds edges only from `L:` header fields and never checked that any record had
/// one. A file with none parses into an edgeless graph, the BFS adds no neighbours, the GFA
/// comes out with `S` lines and no `L` lines, and the accession ends `no_path`, which is the
/// same word the walker uses for a sample that does not carry the gene. Logan publishes such
/// objects, so a cohort counting them as absences understates gene presence.
#[test]
fn an_object_with_no_links_is_its_own_outcome_not_a_biological_negative() {
    let f = Fixture::unlinked("no-links", true);
    // Without --allow-empty a run that reconstructed nothing still exits 1: the outcome is a
    // negative, which is about this accession, not a verdict on the run.
    let out = f.cmd_ga("ga_fast", "cache", &[]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(1), "stderr:\n{stderr}");
    assert!(
        stderr.contains("no L: adjacency"),
        "the run must say what it found; stderr:\n{stderr}"
    );

    let outdir = f.dir.join("out");
    let summary = std::fs::read_to_string(outdir.join("q.summary.json")).unwrap();
    assert!(summary.contains(r#""status": "no_links""#), "{summary}");
    assert!(summary.contains(r#""no_links": 1"#), "the discard counter:\n{summary}");
    assert!(
        !summary.contains(r#""status": "no_path""#),
        "no_path would claim something about the sample:\n{summary}"
    );
    // The tabular form a workflow branches on carries it too.
    let skipped = std::fs::read_to_string(outdir.join("q.skipped.tsv")).unwrap();
    assert!(skipped.contains("SRRNOLINK\tno_links"), "{skipped}");

    // It is a negative like the others, so --allow-empty accepts a run of them.
    let g = Fixture::unlinked("no-links-allow", true);
    let (code, _) = g.code("cache", &["--allow-empty"]);
    assert_eq!(code, 0, "no_links must be tolerated like no_path and too_large");
}

/// WG-07's other half: an edgeless graph is not refused, because a gene that fits inside one
/// unitig needs no edge. Refusing at parse time, which the review proposed, would have turned
/// these into negatives.
#[test]
fn a_gene_inside_one_unitig_still_reconstructs_without_any_links() {
    let f = Fixture::unlinked("no-links-one-unitig", false);
    let out = f.cmd_ga("ga_align", "cache", &[]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "stderr:\n{stderr}");
    let recon = std::fs::read_to_string(f.dir.join("out").join("r.fasta")).unwrap();
    assert_eq!(recon.lines().nth(1).unwrap(), BUBBLE_QUERY, "{recon}");
    let summary = std::fs::read_to_string(f.dir.join("out").join("q.summary.json")).unwrap();
    assert!(summary.contains(r#""status": "success""#), "{summary}");
    assert!(summary.contains(r#""no_links": 0"#), "{summary}");
}

/// WA-04: a primary stitched through two unphased sites says so./// WA-04: a primary stitched through two unphased sites says so.
///
/// `stitch.rs` claimed that requiring a connected path "makes chimeras impossible". That holds
/// for a join no k-mer supports, which is the 87% artifact rate the README measures against
/// concatenated BLAST hits. It does not hold across sites further apart than k - 1: those are
/// linked only through unitigs both strains carry, so every combination of their alleles is a
/// connected path, the walk picks one at each site independently, and the published sequence
/// can be one no molecule carries. The header said `vaf=1.00` and nothing else.
///
/// This fixture is that shape: two bubbles at 40/60 and 45/55, joined through a shared
/// unitig. The primary ends up carrying the 60% allele and the 55% allele, a pair nothing in
/// the graph links.
#[test]
fn a_primary_crossing_two_unphased_sites_is_labelled_mixed() {
    let f = Fixture::two_forks("mixed");
    let out = f.cmd_ga("ga_replay", "cache", &[]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "stderr:\n{stderr}");

    let recon = std::fs::read_to_string(f.dir.join("out").join("r.fasta")).unwrap();
    let header = recon.lines().next().unwrap();
    assert!(
        header.contains(" min_vaf=0.550"),
        "the weakest allele the sequence carries is the 55% one; header:\n{header}"
    );
    assert!(
        header.ends_with(" mixed=1"),
        "two unphased sites is a mixture and must be flagged; header:\n{header}"
    );
    // `vaf=` keeps its old constant for one release, which is exactly why `min_vaf=` had to
    // exist: the two numbers are in the same header and only one of them was ever measured.
    assert!(header.contains(" vaf=1.00 "), "{header}");

    // The sequence is the combination: both majority alleles, carried by no observed molecule.
    let expect = flip(&flip(MIX_QUERY, 60), 121);
    assert_eq!(recon.lines().nth(1).unwrap(), expect, "the primary is not the mosaic");

    let summary = std::fs::read_to_string(f.dir.join("out").join("q.summary.json")).unwrap();
    assert!(
        summary.contains(r#""two_allele_sites": 2, "min_site_vaf": 0.5500"#),
        "the record must carry the numbers the header is drawn from:\n{summary}"
    );
    assert!(summary.contains(r#""mixed": 1,"#), "the cohort count:\n{summary}");

    // One site is a variant, not a phasing problem: the single-bubble fixture must not be
    // flagged, and must still publish the measured frequency rather than 1.00.
    let g = Fixture::bubble("mixed-one-site");
    let out = g.cmd_ga("ga_replay", "cache", &[]).output().unwrap();
    assert_eq!(out.status.code(), Some(0), "{}", String::from_utf8_lossy(&out.stderr));
    let recon = std::fs::read_to_string(g.dir.join("out").join("r.fasta")).unwrap();
    let header = recon.lines().next().unwrap();
    assert!(!header.contains("mixed="), "one site is not a mixture; header:\n{header}");
    assert!(header.contains(" min_vaf=0.600"), "header:\n{header}");
}

/// WA-04: the flag must not depend on whether alternatives were asked for.
///
/// `find_alts` used to return early when no variant output was wanted, so under
/// `--no-emit-alts` nothing scanned the route and a mixed primary would have been published
/// with the same header a clonal one gets.
#[test]
fn the_mixed_flag_survives_no_emit_alts() {
    let f = Fixture::two_forks("mixed-quiet");
    let out = f.cmd_ga("ga_replay", "cache", &["--no-emit-alts"]).output().unwrap();
    assert_eq!(out.status.code(), Some(0), "{}", String::from_utf8_lossy(&out.stderr));
    let recon = std::fs::read_to_string(f.dir.join("out").join("r.fasta")).unwrap();
    let header = recon.lines().next().unwrap();
    assert!(header.ends_with(" mixed=1"), "header:\n{header}");
    assert!(header.contains(" min_vaf=0.550"), "header:\n{header}");
    assert!(
        !f.dir.join("out").join("alt_haplotypes.fasta").exists(),
        "no allele file was asked for"
    );
}

/// WP-02: a `-q` pool that is not one gene must say so./// WP-02: a `-q` pool that is not one gene must say so.
///
/// Multi-anchor mode ranks every reference together and keeps one primary per accession. For
/// strains of a gene that is the feature; for a panel of different genes it answers a
/// different question than the user asked, and 0.1.1 did it silently (0.1.0 at least printed
/// "One gene per query"). The accession here carries the first gene only, which is the
/// ordinary case: the second gene is simply absent, and the user has no way to tell that from
/// a gene that lost the ranking.
#[test]
fn a_multi_anchor_pool_of_unrelated_genes_warns() {
    let f = Fixture::bubble("mixed-genes");
    // An unrelated coding sequence, in frame and long enough to carry 21-mers of its own.
    let other = "ATGGCTGCTAAAGAAGGCGTTCTGGAAATCGGTCAGTGGATCAAAGAGCTGCTGAAAGGCTAAATGGCTGCTAAAGAAGGC";
    std::fs::write(f.dir.join("other.fasta"), format!(">bla_NDM\n{other}\n")).unwrap();
    let mut c = f.cmd_ga("ga_replay", "cache", &[]);
    c.arg("-q").arg(f.dir.join("other.fasta"));
    let out = c.output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("[warn] multi-anchor:") && stderr.contains("2 groups"),
        "a two-gene pool must be called out before the run; stderr:\n{stderr}"
    );
    assert!(stderr.contains("[bla_NDM]"), "the warning must name the genes; stderr:\n{stderr}");

    // The same run with one reference, or with a panel of one gene, stays quiet.
    let g = Fixture::bubble("one-gene");
    let out = g.cmd_ga("ga_replay", "cache", &[]).output().unwrap();
    let quiet = String::from_utf8_lossy(&out.stderr);
    assert!(!quiet.contains("multi-anchor"), "single-anchor runs must not warn:\n{quiet}");
}

/// WG-05: `--min-abund` is the allele floor and must not move the prune floor.
///
/// `depth_floor` returned `max(prune_frac * depth, min_abund)`, so a user tightening the
/// variant gate also raised the coverage at which unitigs are deleted from the graph the
/// primary is walked on. Nothing said so: the `--min-abund` help calls it the minimum for a
/// variant, the `--prune-frac` help calls the prune a fraction of depth, and `summary.json`
/// reported one `floor` without saying which setting had set it.
///
/// This fixture's route is 100x, 40x, 100x, so the default floor is the 5.0 lower bound. Under
/// the old rule `--min-abund 50` moved it to 50.0 and deleted the 40x branch the aligner
/// walked; under the new one the floor is still 5.0, and `--prune-min-abund 50` is how a user
/// asks for the old effect on purpose.
#[test]
fn the_allele_floor_does_not_move_the_prune_floor() {
    let f = Fixture::bubble("prune-floor");
    let out = f.cmd_ga("ga_replay", "cache", &["--min-abund", "50"]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "stderr:\n{stderr}");
    let summary = std::fs::read_to_string(f.dir.join("out").join("q.summary.json")).unwrap();
    // The two floors, side by side in one record: the allele gate is where the user put it,
    // and the prune is where its own default puts it.
    assert!(
        summary.contains(r#""prune": {"outcome": "nothing_to_remove", "depth": 100.0, "floor": 5.0, "unitigs": 4, "unitigs_kept": 4"#),
        "--min-abund 50 must leave the prune floor and the graph alone:\n{summary}"
    );
    assert!(summary.contains(r#""variants": {"depth": 100.0, "floor": 50.0"#), "{summary}");
    assert!(summary.contains(r#""min_abund": 50.0"#), "{summary}");
    assert!(summary.contains(r#""prune_min_abund": 5.0"#), "{summary}");

    // The prune's own lower bound still does what the joined setting used to do.
    let g = Fixture::bubble("prune-floor-explicit");
    let out = g
        .cmd_ga("ga_replay", "cache", &["--prune-min-abund", "50"])
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(0), "{}", String::from_utf8_lossy(&out.stderr));
    let summary = std::fs::read_to_string(g.dir.join("out").join("q.summary.json")).unwrap();
    assert!(
        summary.contains(r#""floor": 50.0, "unitigs": 4, "unitigs_kept": 3"#),
        "--prune-min-abund must be the setting that moves the floor and deletes the 40x \
         branch:\n{summary}"
    );
}

/// WG-02: Logan's `km:f:` abundance alias must be read as the abundance it is.
///
/// Logan's `Unitigs.md` says the `ka:f:` field "is replaced by `km:f:`" in some large
/// accessions and that they are the same information. The parser matched only `ka:f:`, so such
/// an accession parsed cleanly with every unitig at 0x. Nothing downstream calls that an
/// error: it is simply a gene with no coverage, so the prune floor collapses, no fork has a
/// majority to correct towards, no branch clears `--min-abund`, and the VCF comes back empty,
/// while `status` stays `success`.
///
/// Same graph, same frequencies and same default flags as
/// `a_fork_is_phased_and_the_majority_allele_becomes_the_primary`. Only the spelling of the
/// tag differs, so every assertion below is that one test's, and the alias must not change a
/// single byte of the answer.
#[test]
fn the_km_f_abundance_alias_is_read_like_ka_f() {
    let f = Fixture::bubble("km-alias");
    f.retag("SRRFAKE0", |fa| fa.replace("ka:f:", "km:f:"));
    let out = f.cmd_ga("ga_replay", "cache", &["--emit-vcf"]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "stderr:\n{stderr}");
    assert!(
        stderr.contains("consensus: 1 site(s) where the route carried the minority allele, corrected"),
        "the alias left the graph at 0x, so the majority correction never ran; stderr:\n{stderr}"
    );

    let outdir = f.dir.join("out");
    let recon = std::fs::read_to_string(outdir.join("r.fasta")).unwrap();
    assert_eq!(
        recon.lines().nth(1).unwrap(),
        BUBBLE_MAJORITY,
        "with the alias read, the primary must still carry the 60% allele"
    );
    // The depth-based stages are the ones that go silent, so check one of each: phasing,
    // the alt recovery and the variant writer.
    let bubbles = std::fs::read_to_string(outdir.join("bubbles.json")).unwrap();
    assert!(bubbles.contains("\"n_bubbles\": 1"), "{bubbles}");
    assert!(bubbles.contains(r#"{"node": 2, "strand": "+", "ka": 60.0, "vaf": 0.6000}"#), "{bubbles}");
    let alts = std::fs::read_to_string(outdir.join("alt_haplotypes.fasta")).unwrap();
    assert!(alts.contains("alt_kaf=40.0"), "{alts}");
    let vcf = std::fs::read_to_string(outdir.join("SRRFAKE0.vcf")).unwrap();
    assert!(vcf.contains("AF=0.6000;AO=60.0;RO=40.0"), "{vcf}");

    // Every record carries an abundance, so the WG-02 counter stays silent and the summary
    // keeps the shape an ordinary run has.
    let summary = std::fs::read_to_string(outdir.join("q.summary.json")).unwrap();
    assert!(!summary.contains("abundance_missing"), "{summary}");
    assert!(!stderr.contains("carry no ka:f/km:f abundance"), "{stderr}");
}

/// WG-02's other half: an accession whose records carry neither tag is not a 0x gene, it is an
/// accession whose abundances are unknown, and the run must say so somewhere a machine reads.
#[test]
fn unitigs_with_no_abundance_tag_are_reported_not_read_as_zero() {
    let f = Fixture::bubble("no-abundance");
    // Strip the tag from the two branch unitigs, leaving the flanks tagged: a partial loss is
    // the case a whole-file check would miss.
    f.retag("SRRFAKE0", |fa| fa.replace("ka:f:40.0  ", "").replace("ka:f:60.0  ", ""));
    let out = f.cmd_ga("ga_replay", "cache", &[]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "stderr:\n{stderr}");
    assert!(
        stderr.contains("2 of 4 unitigs carry no ka:f/km:f abundance"),
        "the run must warn about the records it could not read a depth from; stderr:\n{stderr}"
    );
    let summary =
        std::fs::read_to_string(f.dir.join("out").join("q.summary.json")).unwrap();
    assert!(
        summary.contains("\"abundance_missing\": 2"),
        "the count must reach summary.json:\n{summary}"
    );
}

/// WV-01, fixed: the majority correction must fire on a flanked window too.
///
/// `variants::prefer_majority` used to reject its swap unless the restitched sequence was as
/// long as the whole stitched path, but `stitch` returns only the aligned window. The two
/// agree only when `path_start == 0` and `path_end >= path_len`, so for any alignment that
/// started or ended inside the path, which is every real GAF row, the correction was skipped
/// silently: `_primary` kept the minority allele, no `consensus:` line was printed, and the
/// output was byte-identical to a `--no-prefer-majority` run. Present in 0.1.0 and 0.1.1.
///
/// This is the same graph, the same frequencies and the same default flags as
/// `a_fork_is_phased_and_the_majority_allele_becomes_the_primary`. Only the window differs: a
/// 6-base flank, the shape of the one real GAF row the repository records (`align.rs`,
/// ERR2028820, `7..1281` of 1288). Both must now reach the same allele.
#[test]
fn a_flanked_alignment_window_still_gets_the_majority_correction() {
    let f = Fixture::bubble("fork-flanked");
    f.set_gaf(6, 115, 6, 115);
    let out = f.cmd_ga("ga_replay", "cache", &["--emit-vcf"]).output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(0), "stderr:\n{stderr}");
    assert!(
        stderr.contains("consensus: 1 site(s) where the route carried the minority allele, corrected"),
        "the correction did not fire on a flanked window, which is WV-01; stderr:\n{stderr}"
    );

    // The published record is the window, not the whole path: 109 of 121 nt, carrying the 60%
    // allele. Slicing the expected haplotypes by the same window is what distinguishes a
    // correct shift from one that merely has the right length.
    let outdir = f.dir.join("out");
    let primary = std::fs::read_to_string(outdir.join("r.fasta")).unwrap();
    assert_eq!(
        primary.lines().nth(1).unwrap(),
        &BUBBLE_MAJORITY[6..115],
        "the primary must be the majority allele over the aligned window"
    );
    assert!(primary.starts_with(">SRRFAKE0_primary cov=0.901 "), "{primary}");

    // The reference-matching 40% allele is demoted to the alternative, over the same window.
    let alts = std::fs::read_to_string(outdir.join("alt_haplotypes.fasta")).unwrap();
    assert_eq!(alts.lines().nth(1).unwrap(), &BUBBLE_QUERY[6..115]);
    assert!(alts.contains("vaf=0.400"), "{alts}");
    assert!(alts.contains("primary_vaf=0.600"), "{alts}");
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

    // A cache with no record of where its bytes came from is not a frozen input, it is a
    // directory. Until 0.2.0 this path wrote nothing at all.
    let summary =
        std::fs::read_to_string(f.dir.join("out").join("fetch.summary.json")).unwrap_or_else(
            |e| panic!("--fetch-only must write fetch.summary.json: {e}; stderr:\n{stderr}"),
        );
    for needle in [
        "\"schema\": 2",
        "\"name\": \"logan-walker\"",
        // Nothing was reconstructed and GraphAligner never ran, so neither is claimed.
        "\"external\": []",
        "\"references\": null",
        "\"n_requested\": 2",
        "\"acc\": \"SRRFAKE0\"",
        "\"acc\": \"SRRFAKE1\"",
        "\"payload\": {\"source\": \"cache\", \"bytes\": ",
    ] {
        assert!(summary.contains(needle), "fetch summary missing {needle}:\n{summary}");
    }
    assert!(!summary.contains("GraphAligner"), "no tool that did no work:\n{summary}");
    assert_eq!(summary.matches("\"sha256\": ").count(), 2, "one digest per accession:\n{summary}");
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

/// A graph that would exceed its memory budget is refused, and refused as the `too_large`
/// NEGATIVE rather than as an error.
///
/// The distinction is the whole point. `--allow-empty` accepts negatives and not errors, and a
/// scheduler reads the exit code: mapping this to `error` would turn "this accession does not
/// fit the cap you set" into "this run broke". The compute pool used to send every error it
/// saw through `AccessionResult::failed`, so until the pool learned the same mapping the fetch
/// pool had, the budget refusal arrived as a failure.
#[test]
fn a_graph_over_its_budget_is_a_too_large_negative_not_an_error() {
    let f = Fixture::dense("dense-budget");
    let out = Command::new(BIN)
        .env_remove("LOGAN_KEEP_SCRATCH")
        .arg("--scratch-dir").arg(f.dir.join("scratch"))
        .arg("--graphaligner").arg(f.dir.join("ga_fast"))
        .arg("-l").arg(f.dir.join("accs.txt"))
        .arg("-q").arg(f.dir.join("q.fasta"))
        .arg("--unitigs-cache").arg(f.dir.join("cache"))
        .arg("--offline")
        .arg("-o").arg(f.dir.join("out"))
        .arg("-f").arg(f.dir.join("out").join("r.fasta"))
        // `--max-subgraph-nodes 0` turns off the neighbourhood cap, which is what normally
        // sets the graph budget; without it the budget falls back to `--max-size-mb` x 1.7,
        // so 1 MB means 1.7 MB. The compressed fixture is far below the size cap itself, so
        // only the graph budget can refuse it. Driving it from the size cap is the point:
        // this is the branch where the gene does not bound the graph and the file must.
        .arg("--max-subgraph-nodes").arg("0")
        .arg("--max-size-mb").arg("1")
        .arg("--allow-empty")
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    // The per-accession outcome line goes to stdout (`pipeline::report`), not stderr. Both are
    // read here so the assertion does not depend on which stream a given outcome prints to.
    let log = format!("{stderr}{}", String::from_utf8_lossy(&out.stdout));
    assert_eq!(
        out.status.code(),
        Some(0),
        "a negative under --allow-empty must exit 0; log:\n{log}"
    );
    assert!(log.contains("too_large"), "log:\n{log}");

    let summary = std::fs::read_to_string(f.dir.join("out").join("q.summary.json")).unwrap();
    assert!(summary.contains("\"status\": \"too_large\""), "{summary}");
    assert!(
        !summary.contains("\"status\": \"failed\"\n"),
        "the record must be a negative, not a failure: {summary}"
    );
    assert!(summary.contains("\"too_large\": 1"), "the discard must be counted: {summary}");
    // Refused on the way, not after building the whole graph.
    assert!(summary.contains("budget"), "the record should say what refused it: {summary}");
    // And the accession still names the bytes it was refused on.
    assert!(summary.contains("\"payload\": {\"source\": \"cache\""), "{summary}");

    // The same input is fine with room to build it, so the refusal is the budget, not the data.
    let ok = Command::new(BIN)
        .env_remove("LOGAN_KEEP_SCRATCH")
        .arg("--scratch-dir").arg(f.dir.join("scratch2"))
        .arg("--graphaligner").arg(f.dir.join("ga_fast"))
        .arg("-l").arg(f.dir.join("accs.txt"))
        .arg("-q").arg(f.dir.join("q.fasta"))
        .arg("--unitigs-cache").arg(f.dir.join("cache"))
        .arg("--offline")
        .arg("-o").arg(f.dir.join("out2"))
        .arg("-f").arg(f.dir.join("out2").join("r.fasta"))
        .arg("--max-size-mb").arg("64")
        .arg("--allow-empty")
        .output()
        .unwrap();
    let ok_log = format!(
        "{}{}",
        String::from_utf8_lossy(&ok.stderr),
        String::from_utf8_lossy(&ok.stdout)
    );
    assert_eq!(ok.status.code(), Some(0), "log:\n{ok_log}");
    // The report line, not the bare word: the startup lines mention `too_large` when they
    // explain what an unspillable accession becomes, and a substring match on that would pass
    // or fail on wording rather than on behaviour.
    assert!(
        !ok_log.contains("SRRDENSE: too_large"),
        "a 64 MB cap gives the graph room, so nothing should be refused:\n{ok_log}"
    );
}

/// A neighbourhood over `--max-subgraph-nodes` is refused, as its own negative.
///
/// This is the memory term the two-pass parser cannot remove: everything after extraction is
/// proportional to the neighbourhood, and on amplicon-like data the neighbourhood is most of
/// the graph. It gets a status of its own rather than being folded into `too_large`, because
/// the two mean different things to whoever reads them: one is about the file, the other about
/// the gene, and the remedy differs.
#[test]
fn a_neighbourhood_over_the_cap_is_its_own_negative() {
    let f = Fixture::matching("subgraph-cap", 1);
    let out = f
        .cmd_ga("ga_align", "cache", &["--max-subgraph-nodes", "1", "--allow-empty"])
        .output()
        .unwrap();
    let log = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stderr),
        String::from_utf8_lossy(&out.stdout)
    );
    assert_eq!(out.status.code(), Some(0), "a negative under --allow-empty exits 0:\n{log}");
    assert!(log.contains("subgraph_too_large"), "{log}");

    let summary = std::fs::read_to_string(f.dir.join("out").join("q.summary.json")).unwrap();
    assert!(summary.contains("\"status\": \"subgraph_too_large\""), "{summary}");
    assert!(summary.contains("\"subgraph_too_large\": 1"), "the discard is counted: {summary}");
    // The record still says how big it was and what it was measured against, so the next run
    // can be given a number rather than a guess.
    assert!(summary.contains("--max-subgraph-nodes 1"), "{summary}");
    assert!(summary.contains("\"subgraph_nodes\": 2"), "{summary}");
    // Not folded into the file-size skip, which is a different question.
    assert!(summary.contains("\"too_large\": 0"), "{summary}");

    // Default cap: the same accession reconstructs, so the refusal is the cap not the data.
    let ok = f.cmd_ga("ga_align", "cache", &[]).output().unwrap();
    assert_eq!(ok.status.code(), Some(0), "stderr:\n{}", String::from_utf8_lossy(&ok.stderr));
    let s2 = std::fs::read_to_string(f.dir.join("out").join("q.summary.json")).unwrap();
    assert!(s2.contains("\"status\": \"success\""), "{s2}");
}
