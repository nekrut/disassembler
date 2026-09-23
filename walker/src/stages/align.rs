//! Step 5: GraphAligner.
//!
//! This is the one stage that is NOT ported. GraphAligner is an external C++ binary
//! (bioconda `graphaligner`, ~29 MB) that does seeded, banded, bit-parallel dynamic
//! programming over a bidirected sequence graph. This file is only the wrapper around it.
//!
//! Why it stays external:
//!
//! - it is 0.04 s of a 14.6 s accession, and that figure already includes process spawn,
//!   the GFA read, the GAF write and the read back. The round trip is priced in
//! - the payloads are tiny: ~38 KB out, <1 KB back. On tmpfs that is a memcpy
//! - the subgraph is not a DAG. It is bidirected with cycles, which is exactly why paths
//!   carry reverse-complement steps (`ERR2028820`: 68 forward, 47 reverse of 115)
//! - there is no mature sequence-to-graph aligner crate in Rust. `rust-bio` is pairwise;
//!   spoa/abpoa bindings do partial-order alignment over DAGs
//!
//! The real cost of keeping it is deployment, not speed: the walker cannot ship as a single
//! static binary, so GraphAligner must be on PATH.
//!
//! Batching does not help. `-g` takes exactly one graph and each accession has its own;
//! concatenating subgraphs would return one global best path rather than one per component.
//!
//! IN   Subgraph (GFA path), Query
//! OUT  GafPath

use crate::config::Config;
use crate::error::{Error, Result};
use crate::types::{GafPath, Orient, PathStep, Query, Subgraph};
use std::path::{Path, PathBuf};
use std::process::Command;

/// Build the GraphAligner command line. Mirrors `logan_path_walker.py:174-177` exactly.
///
/// `query_fa` is the FASTA to align, which is `query.path` only when that file holds a single
/// record; see `run` for why.
pub fn command(bin: &Path, cfg: &Config, subgraph: &Subgraph, query_fa: &Path, gaf_out: &Path) -> Command {
    let mut cmd = Command::new(bin);
    cmd.arg("-g")
        .arg(&subgraph.gfa_path)
        .arg("-f")
        .arg(query_fa)
        .arg("-a")
        .arg(gaf_out)
        .arg("-x")
        .arg("dbg")
        .arg("-t")
        .arg(cfg.graphaligner_threads.to_string());
    cmd
}

/// Resolve the binary, preferring an explicit path then `PATH`.
/// Mirrors `find_bin` at `logan_path_walker.py:36-43`, but fails loudly instead of returning
/// the bare name and letting the spawn produce a confusing "No such file or directory".
pub fn resolve_binary(cfg: &Config) -> Result<PathBuf> {
    let name = &cfg.graphaligner_bin;
    let candidate = Path::new(name);
    if candidate.components().count() > 1 || candidate.is_absolute() {
        return if candidate.is_file() {
            Ok(candidate.to_path_buf())
        } else {
            Err(Error::Usage(format!("GraphAligner not found at {name}")))
        };
    }
    let path = std::env::var_os("PATH").ok_or_else(|| Error::Usage("PATH is unset".into()))?;
    for dir in std::env::split_paths(&path) {
        let p = dir.join(name);
        if p.is_file() {
            return Ok(p);
        }
    }
    if name == "GraphAligner" {
        let staging_bin = PathBuf::from("/data/staging/codon_env/bin/GraphAligner");
        if staging_bin.is_file() {
            return Ok(staging_bin);
        }
    }
    Err(Error::Usage(format!(
        "{name} is not on PATH. The walker shells out to GraphAligner; install it or set \
         Config::graphaligner_bin to an absolute path."
    )))
}

/// Run GraphAligner and parse its GAF output.
///
/// Returns `Ok(None)` when the tool ran cleanly but found no alignment, which is a
/// legitimate `Outcome::NoPath`. A crash is an `Err`, never folded into that case: an empty
/// GFA makes GraphAligner abort with SIGABRT, and in the Python walker that exception is
/// caught per accession at `:392` so the process still exits 0 having reconstructed nothing.
///
/// A nodeless subgraph short-circuits before spawning, since handing GraphAligner an empty
/// graph is a known abort rather than a question worth asking.
///
/// # One record in, one record out
///
/// `read_query` selects the first record of the query FASTA and warns about the rest, but
/// handing the whole file to GraphAligner undoes that: it aligns every record, and
/// `parse_gaf` used to take the best-scoring row across all of them. `stitch` then sliced
/// with that row's coordinates while dividing by the *selected* record's length, so a
/// 160-record panel such as `queries/master_planetary_targets.fasta` produced another gene's
/// path published under the first gene's name, at a coverage above 1.0.
///
/// So a multi-record query is written out as its selected record alone before spawning, and
/// `parse_gaf` filters on the GAF's query-name column either way.
pub fn run(
    cfg: &Config,
    subgraph: &Subgraph,
    query: &Query,
    gaf_out: &Path,
) -> Result<Option<GafPath>> {
    if subgraph.nodes.is_empty() {
        return Ok(None);
    }
    let bin = resolve_binary(cfg)?;

    // Single-record queries (every golden fixture, and every well-formed run) are handed to
    // GraphAligner untouched, so their behaviour is bit-for-bit what it was.
    let single_record_fa = if query.n_records > 1 {
        let p = gaf_out.with_extension("query.fa");
        let mut body = Vec::with_capacity(query.seq.len() + query.record_id.len() + 3);
        body.push(b'>');
        body.extend_from_slice(query.record_id.as_bytes());
        body.push(b'\n');
        body.extend_from_slice(&query.seq);
        body.push(b'\n');
        std::fs::write(&p, &body)?;
        eprintln!(
            "[{}] query has {} records; aligning only {:?}",
            subgraph.accession, query.n_records, query.record_id
        );
        Some(p)
    } else {
        None
    };
    let query_fa = single_record_fa.as_deref().unwrap_or(&query.path);

    let mut cmd = command(&bin, cfg, subgraph, query_fa, gaf_out);
    // GraphAligner is chatty on stdout; only stderr is kept, for error reporting.
    cmd.stdout(std::process::Stdio::null());

    let output = cmd.output().map_err(|e| Error::GraphAligner {
        accession: subgraph.accession.clone(),
        detail: format!("could not run {}: {e}", bin.display()),
    })?;

    if !output.status.success() {
        // Remove any partial GAF so a later run cannot read a truncated file.
        let _ = std::fs::remove_file(gaf_out);
        if let Some(p) = &single_record_fa {
            let _ = std::fs::remove_file(p);
        }
        let tail = String::from_utf8_lossy(&output.stderr);
        let tail = tail.lines().rev().take(3).collect::<Vec<_>>().join(" | ");
        return Err(Error::GraphAligner {
            accession: subgraph.accession.clone(),
            detail: format!(
                "exited with {} on a {}-node graph: {tail}",
                output.status,
                subgraph.nodes.len()
            ),
        });
    }

    let bytes = match std::fs::read(gaf_out) {
        Ok(b) => b,
        // Exited cleanly but wrote nothing: no alignment found. This mapping is only sound
        // because `extract::scratch_stem` keys the GAF on PID and a per-task counter, so no
        // other worker or process can have deleted it out from under us; keyed on the
        // accession alone, a file another run unlinked arrived here as a biological negative.
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            if let Some(p) = &single_record_fa {
                let _ = std::fs::remove_file(p);
            }
            return Ok(None);
        }
        Err(e) => {
            let _ = std::fs::remove_file(gaf_out);
            if let Some(p) = &single_record_fa {
                let _ = std::fs::remove_file(p);
            }
            return Err(Error::Io(e));
        }
    };
    if std::env::var_os("LOGAN_KEEP_SCRATCH").is_none() {
        let _ = std::fs::remove_file(gaf_out);
        if let Some(p) = &single_record_fa {
            let _ = std::fs::remove_file(p);
        }
    }
    parse_gaf(&bytes, &query.record_id)
}

/// Parse the GAF GraphAligner emits and return the best alignment.
///
/// Column 6 is the path, e.g. `>860013>579233<27101`. The best row is the one with the most
/// residue matches (column 10), chosen with a stable sort so ties keep file order, matching
/// `logan_path_walker.py:193`.
pub fn parse_gaf(gaf_bytes: &[u8], expect: &str) -> Result<Option<GafPath>> {
    let text = std::str::from_utf8(gaf_bytes)
        .map_err(|e| Error::Usage(format!("GAF is not valid UTF-8: {e}")))?;

    let mut rows: Vec<Vec<&str>> = Vec::new();
    for line in text.lines() {
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        // 12 columns is the GAF minimum; anything shorter cannot carry the fields we need.
        if f.len() >= 11 {
            rows.push(f);
        }
    }
    if rows.is_empty() {
        return Ok(None);
    }

    // GAF column 1 is the query name, which GraphAligner copies from the FASTA record
    // (`pfcrt`, `hiv1_pol`). Taking the best row without checking it is what let one gene's
    // path be published under another gene's name.
    let mine: Vec<Vec<&str>> = rows.iter().filter(|r| r[0] == expect).cloned().collect();
    let mut rows = if !mine.is_empty() {
        mine
    } else {
        // Only one record was ever handed to GraphAligner, so these rows are this query
        // under a name the tool spelled differently. Say so rather than silently dropping
        // the alignment, and carry on with what is there.
        let mut names: Vec<&str> = rows.iter().map(|r| r[0]).collect();
        names.sort_unstable();
        names.dedup();
        eprintln!(
            "[warn] GAF carries no row for query {expect:?} (names present: {names:?}); \
             using the best row anyway"
        );
        rows
    };

    let matches_of = |r: &Vec<&str>| r[9].parse::<usize>().unwrap_or(0);
    rows.sort_by_key(|r| std::cmp::Reverse(matches_of(r)));
    let best = &rows[0];

    let num = |i: usize| -> Result<usize> {
        best[i].parse::<usize>().map_err(|e| {
            Error::Usage(format!("GAF column {} is not a number: {e}", i + 1))
        })
    };
    let steps = parse_path(best[5]);
    if steps.is_empty() {
        return Ok(None);
    }

    Ok(Some(GafPath {
        steps,
        path_len: num(6)?,
        path_start: num(7)?,
        path_end: num(8)?,
        num_matches: num(9)?,
        aln_len: num(10)?,
    }))
}

/// Split a GAF path string into oriented steps. `>` is forward, `<` reverse-complement.
/// No regex: the grammar is a sign followed by digits, repeated.
fn parse_path(s: &str) -> Vec<PathStep> {
    let b = s.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let orient = match b[i] {
            b'>' => Orient::Forward,
            b'<' => Orient::Reverse,
            _ => {
                i += 1;
                continue;
            }
        };
        i += 1;
        let start = i;
        let mut id: u64 = 0;
        while i < b.len() && b[i].is_ascii_digit() {
            id = id * 10 + (b[i] - b'0') as u64;
            i += 1;
        }
        if i > start && id <= u32::MAX as u64 {
            out.push(PathStep { orient, node: id as u32 });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A real GAF line captured from GraphAligner on ERR2028820.
    const REAL: &str = "pfcrt\t1275\t0\t1274\t+\t>860013>579233>390349<27101\t1288\t7\t1281\t1268\t1274\t60";

    #[test]
    fn parses_a_real_gaf_line() {
        let g = parse_gaf(REAL.as_bytes(), "pfcrt").unwrap().expect("a path");
        assert_eq!(g.path_len, 1288);
        assert_eq!(g.path_start, 7);
        assert_eq!(g.path_end, 1281);
        assert_eq!(g.num_matches, 1268);
        assert_eq!(g.aln_len, 1274);
        assert_eq!(g.steps.len(), 4);
        assert_eq!(g.steps[0].node, 860013);
        assert_eq!(g.steps[0].orient, Orient::Forward);
        assert_eq!(g.steps[3].node, 27101);
        assert_eq!(g.steps[3].orient, Orient::Reverse);
    }

    #[test]
    fn path_parsing_handles_long_mixed_orientation() {
        let steps = parse_path(">1<22>333<4444>5");
        let got: Vec<(Orient, u32)> = steps.iter().map(|s| (s.orient, s.node)).collect();
        assert_eq!(got, vec![
            (Orient::Forward, 1), (Orient::Reverse, 22), (Orient::Forward, 333),
            (Orient::Reverse, 4444), (Orient::Forward, 5),
        ]);
    }

    /// Best row is the one with the most residue matches (column 10), not the first.
    #[test]
    fn picks_the_row_with_most_matches() {
        let gaf = "q\t10\t0\t10\t+\t>1\t10\t0\t10\t50\t60\t60\n\
             q\t10\t0\t10\t+\t>2\t10\t0\t10\t900\t1000\t60\n\
             q\t10\t0\t10\t+\t>3\t10\t0\t10\t100\t120\t60\n".to_string();
        let g = parse_gaf(gaf.as_bytes(), "q").unwrap().unwrap();
        assert_eq!(g.num_matches, 900);
        assert_eq!(g.steps[0].node, 2);
    }

    /// Ties keep file order, matching Python's stable `sorted(..., reverse=True)`.
    #[test]
    fn ties_keep_file_order() {
        let gaf = "q\t10\t0\t10\t+\t>7\t10\t0\t10\t500\t600\t60\n\
                   q\t10\t0\t10\t+\t>9\t10\t0\t10\t500\t600\t60\n";
        let g = parse_gaf(gaf.as_bytes(), "q").unwrap().unwrap();
        assert_eq!(g.steps[0].node, 7);
    }

    #[test]
    fn empty_or_unusable_gaf_is_none_not_an_error() {
        assert!(parse_gaf(b"", "q").unwrap().is_none(), "empty file");
        assert!(parse_gaf(b"\n\n", "q").unwrap().is_none(), "blank lines");
        assert!(parse_gaf(b"too\tfew\tcolumns\n", "q").unwrap().is_none(), "short row");
        // Well-formed row whose path names no nodes.
        let g = "q\t10\t0\t10\t+\t*\t10\t0\t10\t5\t6\t60\n";
        assert!(parse_gaf(g.as_bytes(), "q").unwrap().is_none(), "empty path");
    }

    /// The bug this guards: a multi-record query used to let another gene's row win, and
    /// `stitch` then divided that row's span by the selected record's length.
    #[test]
    fn rows_for_another_query_are_ignored() {
        let gaf = "pfcrt\t10\t0\t10\t+\t>1\t10\t0\t10\t50\t60\t60\n\
                   hiv1_pol\t10\t0\t10\t+\t>2\t10\t0\t10\t900\t1000\t60\n";
        let g = parse_gaf(gaf.as_bytes(), "pfcrt").unwrap().unwrap();
        assert_eq!(g.num_matches, 50, "the higher-scoring hiv1_pol row must not win");
        assert_eq!(g.steps[0].node, 1);

        // Asking for the other gene selects the other row, from the same bytes.
        let g = parse_gaf(gaf.as_bytes(), "hiv1_pol").unwrap().unwrap();
        assert_eq!(g.num_matches, 900);
        assert_eq!(g.steps[0].node, 2);
    }

    /// Only one record is ever handed to GraphAligner, so an unrecognised name means the
    /// tool spelled it differently, not that a foreign gene crept in. Warn, do not drop.
    #[test]
    fn unknown_query_name_still_parses() {
        let gaf = "some_other_spelling\t10\t0\t10\t+\t>5\t10\t0\t10\t7\t8\t60\n";
        let g = parse_gaf(gaf.as_bytes(), "pfcrt").unwrap().expect("a path");
        assert_eq!(g.steps[0].node, 5);
    }

    #[test]
    fn non_numeric_column_is_an_error() {
        let g = "q\t10\t0\t10\t+\t>1\tNOTANUMBER\t0\t10\t5\t6\t60\n";
        assert!(parse_gaf(g.as_bytes(), "q").is_err());
    }

    #[test]
    fn missing_binary_is_a_clear_usage_error() {
        let cfg = Config {
            graphaligner_bin: "definitely-not-a-real-binary-xyz".into(),
            ..Config::default()
        };
        let e = resolve_binary(&cfg).expect_err("should fail").to_string();
        assert!(e.contains("not on PATH"), "{e}");
    }

    #[test]
    fn absolute_path_that_does_not_exist_is_reported_as_such() {
        let cfg = Config { graphaligner_bin: "/nope/GraphAligner".into(), ..Config::default() };
        let e = resolve_binary(&cfg).expect_err("should fail").to_string();
        assert!(e.contains("not found at"), "{e}");
    }

    #[test]
    fn command_matches_the_python_invocation() {
        let cfg = Config { graphaligner_threads: 2, ..Config::default() };
        let sub = Subgraph {
            accession: "ACC".into(), nodes: vec![1], seqs: vec![b"ACGT".to_vec()],
            gfa_path: "/tmp/nonexistent-logan-test/ACC.sub.gfa".into(), n_edges: 0,
        };
        let cmd = command(Path::new("/usr/bin/GraphAligner"), &cfg, &sub,
                          Path::new("/tmp/q.fasta"), Path::new("/tmp/out.gaf"));
        let args: Vec<String> = cmd.get_args().map(|a| a.to_string_lossy().into_owned()).collect();
        assert_eq!(args, vec![
            "-g", "/tmp/nonexistent-logan-test/ACC.sub.gfa",
            "-f", "/tmp/q.fasta",
            "-a", "/tmp/out.gaf",
            "-x", "dbg",
            "-t", "2",
        ]);
        std::mem::forget(sub); // avoid Drop trying to unlink a path that was never created
    }
}

/// End-to-end against the real GraphAligner binary. Needs it on PATH.
///
/// ```text
/// PATH=$REPO/.condaenv-walker/bin:$PATH LOGAN_TEST_DIR=/dir/with/zst \
///     cargo test --release -- --ignored --nocapture align_real
/// ```
#[cfg(test)]
mod align_real {
    use super::*;
    use crate::io::input::read_query;
    use crate::stages::{anchor, extract, parse, stream};
    use crate::types::FetchedUnitigs;

    #[test]
    #[ignore]
    fn aligns_err2028820_against_pfcrt() {
        let Ok(dir) = std::env::var("LOGAN_TEST_DIR") else {
            eprintln!("SKIP: set LOGAN_TEST_DIR");
            return;
        };
        let Ok(compressed) = std::fs::read(format!("{dir}/ERR2028820.unitigs.fa.zst")) else {
            eprintln!("SKIP: no data");
            return;
        };
        let mut cfg = Config::default();
        cfg.scratch_dir = std::env::temp_dir().join(format!("logan_align_{}", std::process::id()));
        if resolve_binary(&cfg).is_err() {
            eprintln!("SKIP: GraphAligner not on PATH");
            return;
        }

        let f = FetchedUnitigs {
            index: 0,
            accession: "ERR2028820".into(),
            compressed_len: compressed.len(),
            compressed,
            fetch_secs: 0.0,
        };
        let graph = { let rd = stream::decoder(&f).unwrap(); parse::build(&f.accession, rd).unwrap() };
        let query = read_query(Path::new("../queries/pfcrt_only.fasta")).unwrap();
        let idx = anchor::index_query(&cfg, &query).unwrap();
        let a = anchor::find(&cfg, &graph, &idx).unwrap();
        let nodes = extract::expand(&cfg, &graph, &a).unwrap();
        let seqs = extract::collect_seqs(&graph, &nodes).unwrap();
        let sub = extract::write_gfa(&cfg, &graph, nodes, seqs).unwrap();

        let gaf = sub.gfa_path.with_extension("gaf");
        let t0 = std::time::Instant::now();
        let path = run(&cfg, &sub, &query, &gaf).expect("align").expect("a path");
        let elapsed = t0.elapsed();

        let fwd = path.steps.iter().filter(|s| s.orient == Orient::Forward).count();
        eprintln!(
            "  {} nodes in path ({} forward, {} reverse), path_len {}, {}..{}, \
             {}/{} matches = {:.2}% in {:.3}s",
            path.steps.len(), fwd, path.steps.len() - fwd, path.path_len,
            path.path_start, path.path_end, path.num_matches, path.aln_len,
            100.0 * path.num_matches as f64 / path.aln_len as f64,
            elapsed.as_secs_f64()
        );

        assert!(!path.steps.is_empty());
        assert!(path.path_end > path.path_start);
        assert!(path.num_matches <= path.aln_len);
        // Every node named by the path must exist in the GFA we handed over.
        let present: std::collections::HashSet<u32> = sub.nodes.iter().copied().collect();
        for s in &path.steps {
            assert!(present.contains(&s.node), "path names node {} not in the subgraph", s.node);
        }
        // The GAF is transient and must not survive the call.
        assert!(!gaf.exists(), "GAF was left behind");
    }

    /// An empty subgraph must not reach GraphAligner: it aborts with SIGABRT on one.
    #[test]
    #[ignore]
    fn empty_subgraph_short_circuits() {
        let mut cfg = Config::default();
        cfg.scratch_dir = std::env::temp_dir().join(format!("logan_align_e_{}", std::process::id()));
        if resolve_binary(&cfg).is_err() {
            eprintln!("SKIP: GraphAligner not on PATH");
            return;
        }
        std::fs::create_dir_all(&cfg.scratch_dir).unwrap();
        let gfa = cfg.scratch_dir.join("EMPTY.sub.gfa");
        std::fs::write(&gfa, b"H\tVN:Z:1.0\n").unwrap();
        let sub = Subgraph {
            accession: "EMPTY".into(), nodes: vec![], seqs: vec![],
            gfa_path: gfa.clone(), n_edges: 0,
        };
        let q = Query {
            name: "q".into(), record_id: "q".into(), seq: b"ACGTACGT".to_vec(),
            path: "../queries/pfcrt_only.fasta".into(),
            n_records: 1,
        };
        let got = run(&cfg, &sub, &q, &cfg.scratch_dir.join("EMPTY.sub.gaf")).unwrap();
        assert!(got.is_none(), "must report no path rather than invoking the tool");
    }
}

/// GraphAligner must give the same answer twice, or nothing downstream can be pinned.
///
/// `LOGAN_TEST_DIR=/dir/with/zst cargo test --release -- --ignored --nocapture determinism`
#[cfg(test)]
mod determinism {
    use super::*;
    use crate::config::Config;
    use crate::io::input::read_query;
    use crate::stages::{anchor, extract, parse, stream};
    use crate::types::FetchedUnitigs;
    use std::path::Path;

    /// Five identical invocations must produce the same best alignment.
    ///
    /// This guards the reason `graphaligner_threads` defaults to 1. At the previous default
    /// of 2 this test fails on `SRR8209082`: the best score comes back 2691 on most runs and
    /// 2686 on some, and the pipeline then emits a different haplotype depending on which.
    #[test]
    #[ignore]
    fn the_same_graph_aligns_the_same_way_every_time() {
        let Ok(dir) = std::env::var("LOGAN_TEST_DIR") else {
            eprintln!("SKIP: set LOGAN_TEST_DIR");
            return;
        };
        for (acc, qpath) in [
            ("SRR8209082", "../queries/query_hiv1_pol.fasta"),
            ("ERR2028820", "../queries/pfcrt_only.fasta"),
        ] {
            let Ok(compressed) = std::fs::read(format!("{dir}/{acc}.unitigs.fa.zst")) else {
                eprintln!("SKIP: no data for {acc}");
                continue;
            };
            let mut cfg = Config::default();
            cfg.scratch_dir = std::env::temp_dir().join(format!("logan_det_{}", std::process::id()));
            if resolve_binary(&cfg).is_err() {
                eprintln!("SKIP: GraphAligner not on PATH");
                return;
            }
            let f = FetchedUnitigs {
                index: 0,
                accession: acc.into(),
                compressed_len: compressed.len(),
                compressed,
                fetch_secs: 0.0,
            };
            let graph = {
                let rd = stream::decoder(&f).unwrap();
                parse::build(&f.accession, rd).unwrap()
            };
            let query = read_query(Path::new(qpath)).unwrap();
            let idx = anchor::index_query(&cfg, &query).unwrap();
            let a = anchor::find(&cfg, &graph, &idx).unwrap();
            let nodes = extract::expand(&cfg, &graph, &a).unwrap();
            let seqs = extract::collect_seqs(&graph, &nodes).unwrap();
            let sub = extract::write_gfa(&cfg, &graph, nodes, seqs).unwrap();

            let mut seen: Vec<(usize, String)> = Vec::new();
            for i in 0..5 {
                let gaf = sub.gfa_path.with_extension(format!("det{i}.gaf"));
                let path = run(&cfg, &sub, &query, &gaf).unwrap().expect("an alignment");
                let steps: String = path
                    .steps
                    .iter()
                    .map(|s| format!("{:?}{}", s.orient, s.node))
                    .collect();
                seen.push((path.num_matches, steps));
            }
            let first = &seen[0];
            for (i, s) in seen.iter().enumerate() {
                assert_eq!(
                    s.0, first.0,
                    "{acc}: run {i} scored {} against {} on run 0; raise of \
                     graphaligner_threads above 1 causes exactly this",
                    s.0, first.0
                );
                assert_eq!(s.1, first.1, "{acc}: run {i} took a different path at the same score");
            }
            eprintln!("  {acc}: 5/5 identical, {} matches, {} steps", first.0, path_len(&first.1));
        }
    }

    fn path_len(s: &str) -> usize {
        s.matches("Forward").count() + s.matches("Reverse").count()
    }
}
