//! Step 5: GraphAligner.
//!
//! This is the one stage the walker does not implement itself. GraphAligner is an external
//! C++ binary (bioconda `graphaligner`, ~29 MB) that does seeded, banded, bit-parallel dynamic
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
use crate::error::{file_err, Error, Result};
use crate::types::{GafPath, Orient, PathStep, Query, Subgraph};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Build the GraphAligner command line. These are the arguments the `bench/golden/` fixtures
/// were captured with, so a change here can move the output.
///
/// `query_fa` is the FASTA to align: the temp `.fa` that `run` writes next to the GAF, never
/// the user's file itself; see `run` for why.
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

/// Resolve the binary: `--graphaligner` as given if it has a directory part, otherwise
/// searched on `PATH`. Fails loudly instead of returning the bare name and letting the spawn
/// produce a confusing "No such file or directory". There is no hardcoded fallback location:
/// the one there was (`/data/staging/codon_env/bin`) belonged to a single old machine.
pub fn resolve_binary(cfg: &Config) -> Result<PathBuf> {
    let name = &cfg.graphaligner_bin;
    let candidate = Path::new(name);
    if candidate.components().count() > 1 || candidate.is_absolute() {
        return if is_executable(candidate) {
            Ok(candidate.to_path_buf())
        } else {
            Err(Error::Usage(format!("GraphAligner not found, or not executable, at {name}")))
        };
    }
    let path = std::env::var_os("PATH").ok_or_else(|| Error::Usage("PATH is unset".into()))?;
    for dir in std::env::split_paths(&path) {
        let p = dir.join(name);
        if is_executable(&p) {
            return Ok(p);
        }
    }
    Err(Error::Usage(format!(
        "{name} is not on PATH. The walker shells out to GraphAligner (bioconda \
         graphaligner=1.0.17b); install it, or pass --graphaligner <PATH>."
    )))
}

fn is_executable(p: &Path) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::metadata(p).map(|m| m.is_file() && m.permissions().mode() & 0o111 != 0).unwrap_or(false)
    }
    #[cfg(not(unix))]
    {
        p.is_file()
    }
}

/// `GraphAligner --version`, first line. Fails if the binary does not run at all (a missing
/// shared library, the wrong architecture), which is better found at startup than on the
/// first accession.
pub fn version(bin: &Path) -> Result<String> {
    let out = Command::new(bin)
        .arg("--version")
        .output()
        .map_err(|e| Error::Usage(format!("cannot run GraphAligner at {}: {e}", bin.display())))?;
    let text = String::from_utf8_lossy(if out.stdout.is_empty() { &out.stderr } else { &out.stdout });
    match text.lines().next().map(str::trim).filter(|l| !l.is_empty()) {
        Some(line) if out.status.success() => Ok(line.to_string()),
        _ => Err(Error::Usage(format!(
            "GraphAligner at {} did not report a version ({}): {}",
            bin.display(),
            out.status,
            text.trim()
        ))),
    }
}

/// Run GraphAligner and parse its GAF output, one best alignment per reference in the pool.
///
/// Returns `Ok(HashMap::new())` when the tool ran cleanly but found no alignment against any
/// reference — the caller reads that as `Outcome::NoPath`. A crash is an `Err`, never folded
/// into that case: an empty GFA makes GraphAligner abort with SIGABRT, and folding that into
/// no-path would let the process exit 0 having reconstructed nothing.
///
/// A nodeless subgraph short-circuits before spawning, since handing GraphAligner an empty
/// graph is a known abort rather than a question worth asking.
///
/// # One record in, one record out
///
/// `read_query` selects the first record of a query FASTA and warns about the rest, but
/// handing the whole file to GraphAligner undoes that: it aligns every record, and
/// `parse_gaf` used to take the best-scoring row across all of them. `stitch` then sliced
/// with that row's coordinates while dividing by the *selected* record's length, so a
/// 160-record panel such as `queries/master_planetary_targets.fasta` produced another gene's
/// path published under the first gene's name, at a coverage above 1.0.
///
/// Every reference in `queries` is therefore written out as its own `>record_id` record in a
/// merged temp FASTA, and `parse_gaf_multi` keeps the best row per query name so the caller
/// can pick a primary against a known reference length.
///
/// # Normalised, not copied
///
/// That merged FASTA is built from the sequence `read_query` parsed at startup, never from
/// the user's file. Until 0.2.0 a single [`Query`] whose source file held exactly one record
/// took a "byte-compat" path that re-opened that file and copied it through, so GraphAligner
/// read the user's bytes exactly as it would have if handed the file itself. That cost three
/// defects, all of them silent:
///
/// - **The file was re-read on every accession** and again on every prune realignment. A pipe
///   (`-q <(zcat ref.fa.gz)`, `-q /dev/stdin`), which `README.md` invites, is empty after the
///   first read, so GraphAligner got a 0-byte FASTA and every accession came back `no_path`:
///   a run that looks like a clean biological negative. Editing the query mid-run moved the
///   reference under the accessions that had not started yet (WA-01).
/// - **The header reached GraphAligner verbatim.** 1.0.17b takes the whole header line as the
///   read name and writes it unescaped into GAF column 1, so a tab in `>gene<TAB>desc` shifted
///   every column and the row was dropped (`no_path` on every accession), and a non-UTF-8
///   header failed the GAF parse outright (exit 3). `read_query` had already reduced the same
///   header to its first token, so the two stages disagreed about the reference's name (WA-07).
/// - **Whitespace inside sequence lines reached GraphAligner.** `read_query` strips it;
///   1.0.17b strips only a trailing `\r` and trailing spaces, and its seeding restarts at
///   every non-ACGT byte while a blank in the DP reaches `ambiguousMatch`'s `default:`, which
///   calls `std::abort()`. Codon-spaced `ATG GCT AAA` gave `no_path` everywhere; GenBank
///   10-base blocks gave SIGABRT everywhere (WA-08).
///
/// Writing `>{record_id}\n{seq}\n` fixes all three at once, and makes the aligner and the
/// anchoring stage agree on the reference by construction. The sequence is the same bases in
/// the same order, uppercased and unwrapped, so GraphAligner reads the same read: the golden
/// fixtures are unchanged by it.
///
/// # Never the user's file name
///
/// GraphAligner chooses its input parser from the file name and aborts on anything that does
/// not end in `.fa`, `.fasta`, `.fq`, `.fastq` (optionally `.gz`), `.bam` or `.sam`. The walker
/// already validated the query as FASTA by content in `read_query`, so the name the user's
/// file happens to have is irrelevant: a Galaxy dataset arrives as `dataset_<uuid>.dat`. The
/// temp `.fa` is what `-f` gets, and the user's path never reaches it.
pub fn run(
    cfg: &Config,
    subgraph: &Subgraph,
    queries: &[Query],
    gaf_out: &Path,
) -> Result<HashMap<String, GafPath>> {
    if subgraph.nodes.is_empty() {
        return Ok(HashMap::new());
    }
    if queries.is_empty() {
        return Err(Error::Usage("no queries in anchor pool for align".into()));
    }
    let bin = resolve_binary(cfg)?;

    // GraphAligner only sees this file, whatever the user's query is called. The name is
    // keyed on the GAF, which `extract::scratch_stem` makes unique per worker and task.
    let query_fa = gaf_out.with_extension("queries.fa");
    // One `>record_id\n<seq>\n` block per query, in input order, built from what
    // `read_query` parsed at startup. The user's file is never re-opened here; see
    // "Normalised, not copied" above for the three defects that cost.
    let mut body: Vec<u8> = Vec::with_capacity(
        queries.iter().map(|q| q.seq.len() + q.record_id.len() + 3).sum(),
    );
    for q in queries {
        body.push(b'>');
        body.extend_from_slice(q.record_id.as_bytes());
        body.push(b'\n');
        body.extend_from_slice(&q.seq);
        body.push(b'\n');
    }
    std::fs::write(&query_fa, &body).map_err(file_err("write", &query_fa))?;
    // One record out of a panel is the multi-anchor prune pass realigning the winner, on
    // every accession; it is expected, so it is not logged.
    if queries.len() > 1 {
        let names: Vec<&str> = queries.iter().map(|q| q.record_id.as_str()).collect();
        eprintln!(
            "[{}] multi-anchor: aligning against {} references {:?}",
            subgraph.accession, queries.len(), names
        );
    }

    let mut cmd = command(&bin, cfg, subgraph, &query_fa, gaf_out);
    // GraphAligner is chatty on stdout; only stderr is kept, for error reporting.
    cmd.stdout(std::process::Stdio::null());
    cmd.stderr(std::process::Stdio::piped());
    #[cfg(unix)]
    set_up_child(&mut cmd);

    let mut child = cmd.spawn().map_err(|e| {
        let _ = std::fs::remove_file(&query_fa);
        Error::GraphAligner {
            accession: subgraph.accession.clone(),
            detail: format!("could not run {}: {e}", bin.display()),
        }
    })?;
    // Held until the child has been reaped, so a signal arriving meanwhile can stop it
    // (WF-15). Every return below happens after the wait.
    #[cfg(unix)]
    let watch = crate::scratch::watch_child(child.id());

    // Drain stderr on its own thread. `Command::output`, which this used to be, did that for
    // us; waiting with a deadline means doing it by hand, and it is not optional. A pipe
    // nobody reads fills at 64 KiB and blocks the writer there for good, so a chatty aligner
    // would hit the deadline below and be reported as stuck when what is stuck is us.
    let mut pipe = child.stderr.take();
    let drain = std::thread::spawn(move || {
        let mut buf = Vec::new();
        if let Some(p) = pipe.as_mut() {
            let _ = std::io::Read::read_to_end(p, &mut buf);
        }
        buf
    });

    let status = match wait_with_deadline(&mut child, cfg.graphaligner_timeout_secs) {
        Ok(status) => status,
        Err(detail) => {
            #[cfg(unix)]
            drop(watch);
            // The child is dead, so the pipe is closed and this returns at once.
            let _ = drain.join();
            let _ = std::fs::remove_file(gaf_out);
            let _ = std::fs::remove_file(&query_fa);
            return Err(Error::GraphAligner {
                accession: subgraph.accession.clone(),
                detail: format!("{detail} on a {}-node graph", subgraph.nodes.len()),
            });
        }
    };
    #[cfg(unix)]
    drop(watch);
    let stderr = drain.join().unwrap_or_default();

    if !status.success() {
        // Remove any partial GAF so a later run cannot read a truncated file.
        let _ = std::fs::remove_file(gaf_out);
        let _ = std::fs::remove_file(&query_fa);
        let tail = String::from_utf8_lossy(&stderr);
        let tail = tail.lines().rev().take(3).collect::<Vec<_>>().join(" | ");
        return Err(Error::GraphAligner {
            accession: subgraph.accession.clone(),
            detail: format!(
                "exited with {} on a {}-node graph: {tail}",
                status,
                subgraph.nodes.len()
            ),
        });
    }

    let bytes = match std::fs::read(gaf_out) {
        Ok(b) => b,
        // Exited cleanly but wrote nothing: no alignment found. That reading holds only
        // while the scratch directory is still ours, which `extract::scratch_stem`'s PID and
        // per-task counter do not guarantee on their own: until WF-01 another run sharing the
        // base could delete the whole directory mid-alignment, and this branch then published
        // the deletion as the biological answer `no_path`, with exit 0.
        //
        // So check what is actually missing. If the INPUTS are gone too, nothing was aligned
        // and this is not a negative result; it is an error, and the run must say so.
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            if !subgraph.gfa_path.exists() || !query_fa.exists() {
                return Err(Error::GraphAligner {
                    accession: subgraph.accession.to_string(),
                    detail: format!(
                        "scratch files vanished while GraphAligner ran ({} / {}). Another \
                         process removed this run's scratch directory; if several runs share \
                         a --scratch-dir, give each its own.",
                        subgraph.gfa_path.display(),
                        query_fa.display(),
                    ),
                });
            }
            let _ = std::fs::remove_file(&query_fa);
            return Ok(HashMap::new());
        }
        Err(e) => {
            let _ = std::fs::remove_file(gaf_out);
            let _ = std::fs::remove_file(&query_fa);
            return Err(file_err("read", gaf_out)(e));
        }
    };
    if !crate::scratch::keep_scratch() {
        let _ = std::fs::remove_file(gaf_out);
        let _ = std::fs::remove_file(&query_fa);
    }
    let expects: Vec<&str> = queries.iter().map(|q| q.record_id.as_str()).collect();
    let mut map = parse_gaf_multi(&bytes, &expects)?;
    discount_ambiguous_matches(&mut map, queries);
    Ok(map)
}

/// Put the aligner where a signal can reach it, and nothing else can reach the aligner.
///
/// **Its own process group** (WF-15). The signal handler can then stop the whole group with
/// one `kill`, which matters when `--graphaligner` names a wrapper script: killing the process
/// we spawned leaves its children running, and `wait_with_deadline` has the same limitation
/// for the same reason. A group also means the terminal no longer delivers Ctrl-C to the
/// aligner directly, which is why the handler has to do it: that is a change of mechanism, not
/// of behaviour, and it buys the `kill <walker pid>` and `timeout N logan-walker` cases, where
/// the aligners used to survive the run, reparented to init, working on deleted files.
///
/// **No stdin.** A process in its own group that reads from the controlling terminal is
/// stopped with SIGTTIN, and GraphAligner has no reason to read stdin: its reads come from
/// `-f`. Closing it also keeps it clear of the walker's own stdin, which a query may be
/// arriving on (`-q /dev/stdin`).
///
/// **`PR_SET_PDEATHSIG` on Linux.** A SIGKILL or the OOM killer never runs the handler, so the
/// kernel is asked to kill the aligner when the thread that spawned it dies. That thread is
/// the compute worker, which blocks in `wait_with_deadline` until the child is reaped, so the
/// only way it dies first is the walker dying with it. It is set after the fork, so a walker
/// that dies in the window between fork and `prctl` still leaves that one child; nothing
/// cheap closes that, and the deadline bounds it.
#[cfg(unix)]
fn set_up_child(cmd: &mut Command) {
    use std::os::unix::process::CommandExt;
    cmd.stdin(std::process::Stdio::null());
    cmd.process_group(0);
    #[cfg(target_os = "linux")]
    // SAFETY: `pre_exec` runs between fork and exec, where only async-signal-safe calls are
    // allowed. `prctl` is a bare syscall with no allocation, no locking and no libc state.
    unsafe {
        cmd.pre_exec(|| {
            // A failure here is not worth refusing the alignment over: the deadline and the
            // signal handler still cover every case except SIGKILL of the walker.
            libc::prctl(libc::PR_SET_PDEATHSIG, libc::SIGKILL);
            Ok(())
        });
    }
}

/// Wait for GraphAligner, giving up after `timeout_secs`. `0` waits forever.
///
/// Every other wait in the walker has a deadline; this one did not, and it is the only stage
/// that hands control to a program the walker did not write. One call that never returns (a
/// pathological neighbourhood, a stalled filesystem, a wrapper script waiting on a lock) holds
/// its compute worker for good, and because the walker writes the FASTA, summary and VCF after
/// `pipeline::run` returns, a scheduler's wall-time kill then throws away every accession the
/// shard had already finished. The deadline turns that into one errored accession: the rest
/// finish and every output is written.
///
/// `SIGKILL`, not `SIGTERM`: there is nothing to flush, and a process that has stopped
/// answering is not more likely to answer a polite signal. It is reaped immediately so it
/// cannot be left as a zombie.
///
/// Two limits worth knowing. `kill` reaches the process we spawned, so if `--graphaligner`
/// names a wrapper script its grandchildren survive it; killing the process group would fix
/// that and is a larger change. And the walker's own signal handler still does not stop a
/// running aligner, so a SIGTERM to the walker leaves this child reparented to init.
fn wait_with_deadline(
    child: &mut std::process::Child,
    timeout_secs: u64,
) -> std::result::Result<std::process::ExitStatus, String> {
    let deadline = (timeout_secs > 0).then(|| {
        std::time::Instant::now() + std::time::Duration::from_secs(timeout_secs)
    });
    // Poll quickly at first, then back off. A typical alignment is 0.04 s, so a flat 50 ms
    // poll would add more latency than the alignment costs; this keeps the common case under
    // a millisecond and settles at 20 wakeups a second for a long one.
    let mut nap = std::time::Duration::from_micros(200);
    let cap = std::time::Duration::from_millis(50);
    loop {
        match child.try_wait() {
            Ok(Some(status)) => return Ok(status),
            Ok(None) => {}
            Err(e) => return Err(format!("could not wait for GraphAligner: {e}")),
        }
        if deadline.is_some_and(|d| std::time::Instant::now() >= d) {
            let _ = child.kill();
            let _ = child.wait();
            return Err(format!(
                "timed out after {timeout_secs} s and was killed \
                 (raise or disable --graphaligner-timeout, 0 waits forever)"
            ));
        }
        std::thread::sleep(nap);
        nap = (nap * 2).min(cap);
    }
}

/// Take ambiguous reference bases back out of GAF column 10.
///
/// GraphAligner counts column 10 with `characterMatch`, whose `ambiguousMatch` returns true
/// for `N` against any base, and for R/Y/K/M/S/W/B/D/H/V against each of their member bases
/// (`GraphAlignerCommon.h` in 1.0.17b). Every `N` in a reference is therefore a guaranteed
/// match, so a reference padded with `N` scores as if it had matched bases it does not have.
///
/// That inflates `pident` everywhere, and in multi-anchor mode it decides the wrong thing: a
/// reference with a run of `N` beats the genuinely nearest reference wherever the sample
/// differs inside the `N` window, and the winner then drives pruning, majority correction,
/// the alternatives and every VCF coordinate. `queries/master_planetary_targets.fasta` has
/// records with 659 `N` (`zika_ns5_var3`) and 862 `N`/`Y` (`ebola_gp_var2`), so this is the
/// shipped panel, not a hypothetical one.
///
/// Corrected here, once, immediately after parsing, so the first pass, the prune pass,
/// `matched`, `pident` and the candidate report all see the same corrected number.
///
/// The count is exact for an ambiguous base aligned as a match, and conservative otherwise:
/// one placed in an insertion, or an IUPAC code that genuinely mismatched, was not counted in
/// column 10 to begin with, so subtracting it understates rather than overstates identity.
/// Understating is the safe direction here, because the number decides which reference wins.
fn discount_ambiguous_matches(
    map: &mut std::collections::HashMap<String, GafPath>,
    queries: &[Query],
) {
    for q in queries {
        let Some(p) = map.get_mut(q.record_id.as_str()) else {
            continue;
        };
        // The window comes from GAF columns written by another program. A row naming a span
        // outside the reference is malformed, and slicing on it would panic; there is nothing
        // to discount in that case, so leave the row as it is rather than take down the run.
        let Some(window) = q.seq.get(p.query_start..p.query_end) else {
            continue;
        };
        let ambiguous = window
            .iter()
            .filter(|b| !matches!(b, b'A' | b'C' | b'G' | b'T'))
            .count();
        p.num_matches = p.num_matches.saturating_sub(ambiguous);
    }
}

/// Parse a GAF and keep the best alignment for every expected reference name.
///
/// Multi-anchor's stitch stage needs one `GafPath` per reference so it can score each
/// reference against its own length and pick a primary. This is the same logic as
/// [`parse_gaf`] applied per record id: rows are grouped by column 1, the one with the most
/// residue matches (column 10) wins per group, ties keep file order.
///
/// A row is matched on the first word of its query name, which is the record id: GraphAligner
/// names a read by its whole header line, description included. Matching the whole name
/// dropped every row for a header such as `>pfcrt chloroquine resistance transporter`, so a
/// single-record query with a description (10 of the 26 in `queries/markers`) came back
/// `no_path` on every accession, and in multi-anchor mode the prune pass's realignment
/// against a described panel never matched, so pruning never took effect. Rows for other
/// names are dropped; a caller that finds a reference missing from the returned map knows
/// nothing aligned against it.
pub fn parse_gaf_multi(
    gaf_bytes: &[u8],
    expects: &[&str],
) -> Result<HashMap<String, GafPath>> {
    let text = std::str::from_utf8(gaf_bytes)
        .map_err(|e| Error::Usage(format!("GAF is not valid UTF-8: {e}")))?;

    let expected: std::collections::HashSet<&str> = expects.iter().copied().collect();
    let mut best_by_name: HashMap<String, (usize, Vec<&str>)> = HashMap::new();
    for line in text.lines() {
        let line = line.trim_end();
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 {
            continue;
        }
        let name = f[0].split_whitespace().next().unwrap_or("");
        if !expected.contains(name) {
            continue;
        }
        let matches: usize = f[9].parse().unwrap_or(0);
        let take = match best_by_name.get(name) {
            None => true,
            // Strictly greater beats existing; ties keep file order (first-seen wins).
            Some((prev, _)) => matches > *prev,
        };
        if take {
            best_by_name.insert(name.to_string(), (matches, f));
        }
    }

    let mut out: HashMap<String, GafPath> = HashMap::with_capacity(best_by_name.len());
    for (name, (_matches, best)) in best_by_name {
        let num = |i: usize| -> Result<usize> {
            best[i].parse::<usize>().map_err(|e| {
                Error::Usage(format!("GAF column {} is not a number: {e}", i + 1))
            })
        };
        let steps = parse_path(best[5]);
        if steps.is_empty() {
            continue; // path with no nodes: no usable alignment against this ref
        }
        out.insert(name, GafPath {
            steps,
            query_start: num(2)?,
            query_end: num(3)?,
            path_len: num(6)?,
            path_start: num(7)?,
            path_end: num(8)?,
            num_matches: num(9)?,
            aln_len: num(10)?,
        });
    }
    Ok(out)
}

/// Parse the GAF GraphAligner emits and return the best alignment for one reference.
///
/// Column 6 is the path, e.g. `>860013>579233<27101`. The best row is the one with the most
/// residue matches (column 10), chosen with a stable sort so ties keep file order.
///
/// TW-07: this used to be a second, independent implementation with no production caller,
/// while `align::run` has used [`parse_gaf_multi`] since multi-anchor mode landed. Seven
/// tests pinned the dead copy, so the parser that actually runs was covered by three. It is
/// now a thin wrapper, which keeps those tests' single-reference shape and points every one
/// of them at the live code.
pub fn parse_gaf(gaf_bytes: &[u8], expect: &str) -> Result<Option<GafPath>> {
    Ok(parse_gaf_multi(gaf_bytes, &[expect])?.remove(expect))
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

    /// Ties keep file order: the sort is stable.
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
        /// TW-07. A row whose query name is not the one we asked for is NOT returned.
    ///
    /// This test used to assert the opposite, and passed, because it exercised a second
    /// `parse_gaf` implementation that no production path had called since multi-anchor mode
    /// landed. The live parser attributes rows to references by name, which it must: in
    /// multi-anchor mode a row for `Delta` is not an answer about `Omicron`. The lenient
    /// behaviour the old test pinned was never available to a user.
    ///
    /// Upstream, `run` turns this into a warning naming the rows that were present, so a
    /// mangled query header is diagnosable rather than silently a `no_path`.
    #[test]
    fn a_row_for_another_query_is_not_returned() {
        let gaf = "some_other_spelling\t10\t0\t10\t+\t>5\t10\t0\t10\t7\t8\t60\n";
        assert!(
            parse_gaf(gaf.as_bytes(), "pfcrt").unwrap().is_none(),
            "a row naming another reference is not an alignment of this one"
        );
        // The same bytes, asked for by the name they carry, do parse.
        let g = parse_gaf(gaf.as_bytes(), "some_other_spelling").unwrap().expect("a path");
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
        assert!(e.contains("not found, or not executable, at /nope/GraphAligner"), "{e}");
    }

    #[cfg(unix)]
    #[test]
    fn a_file_that_is_not_executable_is_refused_and_a_version_is_read() {
        use std::os::unix::fs::PermissionsExt;
        let dir = std::env::temp_dir().join(format!("logan_align_ver_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let ga = dir.join("GraphAligner");
        std::fs::write(&ga, "#!/bin/sh\necho 'GraphAligner bioconda 1.0.17b-'\n").unwrap();
        let cfg = Config { graphaligner_bin: ga.to_string_lossy().into_owned(), ..Config::default() };
        std::fs::set_permissions(&ga, std::fs::Permissions::from_mode(0o644)).unwrap();
        assert!(resolve_binary(&cfg).is_err(), "a non-executable file is not GraphAligner");
        std::fs::set_permissions(&ga, std::fs::Permissions::from_mode(0o755)).unwrap();
        let bin = resolve_binary(&cfg).unwrap();
        assert_eq!(version(&bin).unwrap(), "GraphAligner bioconda 1.0.17b-");
        std::fs::write(&ga, "#!/bin/sh\nexit 3\n").unwrap();
        assert!(version(&bin).is_err(), "a binary that does not run is found at startup");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Ambiguous reference bases must not count as matched.
    ///
    /// GraphAligner's `ambiguousMatch` makes every `N` a guaranteed match against any base,
    /// so a reference padded with `N` scores as if it had matched bases it does not have.
    /// The repository's own panel carries such records.
    #[test]
    fn ambiguous_reference_bases_are_taken_back_out_of_the_match_count() {
        // 300 nt, 120 of them N in the middle, exactly the shape the review reproduced.
        let mut seq = vec![b'A'; 300];
        seq[90..210].fill(b'N');
        let q = Query {
            name: "panel".into(),
            record_id: "refN".into(),
            seq,
            path: std::path::PathBuf::from("/dev/null"),
            n_records: 2,
        };
        let path = |matches: usize, start: usize, end: usize| GafPath {
            steps: vec![PathStep { orient: Orient::Forward, node: 1 }],
            query_start: start,
            query_end: end,
            path_len: 1,
            path_start: 0,
            path_end: 1,
            num_matches: matches,
            aln_len: 300,
        };

        // The whole reference aligned: GraphAligner reports 300 matches, 120 of them free.
        let mut map = HashMap::from([("refN".to_string(), path(300, 0, 300))]);
        super::discount_ambiguous_matches(&mut map, std::slice::from_ref(&q));
        assert_eq!(map["refN"].num_matches, 180, "the 120 N must not count");

        // Only the clean prefix aligned: nothing to discount.
        let mut map = HashMap::from([("refN".to_string(), path(90, 0, 90))]);
        super::discount_ambiguous_matches(&mut map, std::slice::from_ref(&q));
        assert_eq!(map["refN"].num_matches, 90);

        // A row naming a span past the end of the reference is malformed. There is nothing to
        // discount and nothing to panic over; slicing it would have done the latter.
        let mut map = HashMap::from([("refN".to_string(), path(300, 0, 10_000))]);
        super::discount_ambiguous_matches(&mut map, std::slice::from_ref(&q));
        assert_eq!(map["refN"].num_matches, 300, "left alone, not panicked on");

        // Never below zero, however many ambiguous bases the window holds.
        let mut map = HashMap::from([("refN".to_string(), path(5, 0, 300))]);
        super::discount_ambiguous_matches(&mut map, std::slice::from_ref(&q));
        assert_eq!(map["refN"].num_matches, 0);
    }

    /// An N-laden reference must stop beating the genuinely nearest one.
    ///
    /// The review's reproduction: `refN` is `refA` with 120 N over positions 90-209 and the
    /// sample differs from `refA` at 6 sites inside that window. Uncorrected, `refN` scores
    /// 300 matches to `refA`'s 294 and wins the primary, taking `anchor=` and every VCF
    /// coordinate with it.
    #[test]
    fn an_n_padded_reference_no_longer_wins_the_primary() {
        use crate::stages::stitch::{pick_primary, Candidate};
        let cand = |name: &str, pool: usize, matches: usize| {
            let path = GafPath {
                steps: vec![PathStep { orient: Orient::Forward, node: 1 }],
                query_start: 0, query_end: 300,
                path_len: 1, path_start: 0, path_end: 1,
                num_matches: matches, aln_len: 300,
            };
            let hap = crate::types::Haplotype {
                id: name.into(), seq: vec![], cov: 1.0,
                ref_cov: 1.0, ref_span: (0, 300),
                pident: matches as f64 / 300.0 * 100.0,
                pident_inherited: false, nodes: 1, mean_ka: 100.0, vaf: 1.0, anchor_used: None,
            };
            Candidate {
                record_id: name.into(), pool_index: pool, ref_len: 300,
                path: Some(path), haplotype: Some(hap),
            }
        };
        // What GraphAligner reports, with the N window counted as matched.
        let uncorrected = vec![cand("refA", 0, 294), cand("refN", 1, 300)];
        assert_eq!(
            pick_primary(uncorrected, crate::config::PrimaryMetric::Matched, 0.002).unwrap().0,
            "refN",
            "this is the defect: without the discount the N-padded record wins"
        );
        // After `discount_ambiguous_matches` takes the 120 N back out of refN's count.
        let corrected = vec![cand("refA", 0, 294), cand("refN", 1, 180)];
        assert_eq!(
            pick_primary(corrected, crate::config::PrimaryMetric::Matched, 0.002).unwrap().0,
            "refA"
        );
    }

    // ---- parse_gaf_multi (multi-anchor) ----

    /// One best row per expected query name, from the same bytes.
    #[test]
    fn multi_returns_best_row_per_expected_name() {
        let gaf = "\
wuhan\t10\t0\t10\t+\t>1\t10\t0\t10\t50\t60\t60
wuhan\t10\t0\t10\t+\t>2\t10\t0\t10\t900\t1000\t60
jn1\t10\t0\t10\t+\t>7\t10\t0\t10\t200\t260\t60
jn1\t10\t0\t10\t+\t>8\t10\t0\t10\t150\t180\t60
";
        let m = parse_gaf_multi(gaf.as_bytes(), &["wuhan", "jn1"]).unwrap();
        assert_eq!(m.len(), 2);
        // Best for `wuhan` is the 900-matches row (node 2).
        assert_eq!(m["wuhan"].num_matches, 900);
        assert_eq!(m["wuhan"].steps[0].node, 2);
        // Best for `jn1` is the 200-matches row (node 7).
        assert_eq!(m["jn1"].num_matches, 200);
        assert_eq!(m["jn1"].steps[0].node, 7);
    }

    /// A reference the GAF has no row for is absent from the map. In multi-anchor mode this
    /// means the walker walked nothing against that particular reference; the caller reads
    /// the missing key as "not a candidate for primary".
    #[test]
    fn multi_omits_references_with_no_rows() {
        let gaf = "wuhan\t10\t0\t10\t+\t>1\t10\t0\t10\t50\t60\t60\n";
        let m = parse_gaf_multi(gaf.as_bytes(), &["wuhan", "jn1"]).unwrap();
        assert_eq!(m.len(), 1);
        assert!(m.contains_key("wuhan"));
        assert!(!m.contains_key("jn1"));
    }

    /// Rows whose query name is not in `expects` are ignored. The single-query codepath
    /// used to warn here; multi-anchor mode has no single "expected" to warn against, so
    /// we drop silently and trust the caller.
    #[test]
    fn multi_ignores_rows_whose_name_is_not_expected() {
        let gaf = "\
wuhan\t10\t0\t10\t+\t>1\t10\t0\t10\t50\t60\t60
some_other_gene\t10\t0\t10\t+\t>2\t10\t0\t10\t900\t1000\t60
";
        let m = parse_gaf_multi(gaf.as_bytes(), &["wuhan", "jn1"]).unwrap();
        assert_eq!(m.len(), 1);
        assert_eq!(m["wuhan"].num_matches, 50);
    }

    /// GraphAligner names a read by its whole header line. Matching that whole name dropped
    /// every row for a described header, and the accession came back no_path.
    #[test]
    fn a_header_description_does_not_hide_the_row() {
        let gaf = "pfcrt chloroquine resistance transporter\t10\t0\t10\t+\t>1\t10\t0\t10\t9\t10\t60\n";
        let m = parse_gaf_multi(gaf.as_bytes(), &["pfcrt"]).unwrap();
        assert_eq!(m["pfcrt"].num_matches, 9);
        let other = "pfcrtX extra\t10\t0\t10\t+\t>1\t10\t0\t10\t9\t10\t60\n";
        assert!(parse_gaf_multi(other.as_bytes(), &["pfcrt"]).unwrap().is_empty(),
                "the first word must equal the id, not merely start with it");
    }

    /// A row whose path names no nodes is not a usable alignment.
    #[test]
    fn multi_drops_paths_with_no_nodes() {
        let gaf = "wuhan\t10\t0\t10\t+\t*\t10\t0\t10\t50\t60\t60\n";
        let m = parse_gaf_multi(gaf.as_bytes(), &["wuhan"]).unwrap();
        assert!(m.is_empty());
    }

    /// Empty or blank input is an empty map, not an error.
    #[test]
    fn multi_empty_input_is_empty_map() {
        assert!(parse_gaf_multi(b"", &["wuhan"]).unwrap().is_empty());
        assert!(parse_gaf_multi(b"\n\n", &["wuhan"]).unwrap().is_empty());
    }

    #[test]
    fn the_graphaligner_command_line_is_exact() {
        let cfg = Config { graphaligner_threads: 2, ..Config::default() };
        let sub = Subgraph {
            accession: "ACC".into(), nodes: vec![1], seqs: vec![b"ACGT".to_vec()],
            gfa_path: "/tmp/nonexistent-logan-test/ACC.sub.gfa".into(), n_edges: 0,
            owns_gfa: false,
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
    }

    /// `run`, retrying while the stand-in binary is still "Text file busy".
    ///
    /// Spawning a script this process wrote moments ago can fail with `ETXTBSY`: another test
    /// thread forked while the write fd was still open, its child inherited the fd, and the
    /// kernel refuses to exec a file that anyone holds open for writing. The window closes as
    /// soon as that child execs, so a short retry settles it. This is a property of a parallel
    /// test suite that writes its own binaries, not of the walker.
    #[cfg(unix)]
    fn run_stand_in(
        cfg: &Config,
        sub: &Subgraph,
        queries: &[Query],
        gaf: &Path,
    ) -> Result<HashMap<String, GafPath>> {
        for _ in 0..50 {
            match run(cfg, sub, queries, gaf) {
                Err(Error::GraphAligner { ref detail, .. }) if detail.contains("Text file busy") => {
                    std::thread::sleep(std::time::Duration::from_millis(20));
                }
                other => return other,
            }
        }
        run(cfg, sub, queries, gaf)
    }

    /// GraphAligner picks its parser from the file name, so the user's name for the query
    /// must never reach `-f`: a Galaxy dataset is `dataset_<uuid>.dat`, which it refuses. A
    /// stand-in binary records what it was handed.
    ///
    /// What it is handed is the query the walker parsed, not the user's bytes (WA-07, WA-08).
    /// Each source below is a different way of writing the same reference, and all of them
    /// must arrive as the same two lines. Under the old byte-copy path they did not, and
    /// 1.0.17b failed differently on each: it takes the whole header line as the read name and
    /// writes it unescaped into GAF column 1, so the tab shifted every column and the row was
    /// dropped as `no_path`; it strips only trailing whitespace from a sequence line, so the
    /// internal blanks survived into the DP, where a blank reaches `ambiguousMatch`'s
    /// `default:` and calls `std::abort()`; and the Latin-1 byte came back in the GAF, which
    /// `parse_gaf_multi` reads as UTF-8, failing the accession outright.
    #[cfg(unix)]
    #[test]
    fn graphaligner_gets_the_parsed_query_not_the_users_bytes() {
        use std::os::unix::fs::PermissionsExt;
        let dir = std::env::temp_dir().join(format!("logan_align_name_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let fake = dir.join("FakeGraphAligner");
        std::fs::write(&fake, format!(
            "#!/bin/sh\n\
             while [ $# -gt 0 ]; do case \"$1\" in\n\
             -f) cp \"$2\" '{d}/seen.fa'; echo \"$2\" > '{d}/seen_name'; \
                 name=$(head -1 \"$2\" | cut -c2- | tr -d '\\r'); shift 2;;\n\
             -a) gaf=\"$2\"; shift 2;;\n\
             *) shift;;\n\
             esac; done\n\
             printf '%s\\t12\\t0\\t12\\t+\\t>1\\t12\\t0\\t12\\t12\\t12\\t60\\n' \"$name\" > \"$gaf\"\n",
            d = dir.display(),
        )).unwrap();
        std::fs::set_permissions(&fake, std::fs::Permissions::from_mode(0o755)).unwrap();

        let cfg = Config { graphaligner_bin: fake.to_string_lossy().into_owned(), ..Config::default() };
        let gfa = dir.join("ACC.sub.gfa");
        std::fs::write(&gfa, b"H\tVN:Z:1.0\nS\t1\tACGTACGTAAAA\n").unwrap();
        let sub = Subgraph {
            accession: "ACC".into(), nodes: vec![1], seqs: vec![b"ACGTACGTAAAA".to_vec()],
            gfa_path: gfa, n_edges: 0,
            // This test wrote that file, so Drop may remove it.
            owns_gfa: true,
        };
        let gaf = dir.join("ACC.sub.gaf");

        // The last entry is the first one gzipped: GraphAligner reads plain FASTA only, and
        // was once handed the gzip bytes under a `.fa` name, which it read as no alignment.
        let plain: &[u8] = b">pol some description\nacgtACGT\r\nAAAA\n";
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        std::io::Write::write_all(&mut gz, plain).unwrap();
        let sources: Vec<(&str, Vec<u8>)> = vec![
            ("description, mixed case and CRLF", plain.to_vec()),
            ("a tab in the header", b">pol\tplasmid, 2026\nACGTACGTAAAA\n".to_vec()),
            ("codon-spaced sequence", b">pol\nACG TAC GTA AAA\n".to_vec()),
            // GenBank-style blocks, wrapped, with a trailing tab on one line.
            ("10-base blocks", b">pol\nACGTACGTAA\t\nAA\n".to_vec()),
            // Latin-1 in the description: not UTF-8, and `read_query` reads it lossily.
            ("a Latin-1 header", b">pol isolate S\xe3o\nACGTACGTAAAA\n".to_vec()),
            ("gzipped", gz.finish().unwrap()),
        ];

        for (what, stored) in sources {
            // The user's own name for the file, which GraphAligner would refuse outright.
            let src = dir.join("dataset_1234.dat");
            std::fs::write(&src, &stored).unwrap();
            let q = crate::io::input::read_query(&src).unwrap();
            let got = run_stand_in(&cfg, &sub, std::slice::from_ref(&q), &gaf).unwrap();

            let seen_name = std::fs::read_to_string(dir.join("seen_name")).unwrap();
            assert!(seen_name.trim_end().ends_with(".fa"),
                    "{what}: GraphAligner was handed {seen_name:?}");
            assert_eq!(std::fs::read(dir.join("seen.fa")).unwrap(), b">pol\nACGTACGTAAAA\n",
                       "{what}: GraphAligner must get the parsed record, normalised");
            assert!(got.contains_key("pol"), "{what}: no row for the reference");
            assert!(!dir.join("ACC.sub.queries.fa").exists(),
                    "{what}: the temp query was left behind");
        }

        drop(sub);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// WA-06: an aligner that never returns fails its own accession rather than holding its
    /// compute worker for good. Before the deadline, the only thing that ended such a run was
    /// the scheduler's wall-time kill, and because every output is written after
    /// `pipeline::run` returns, that kill threw away every accession the shard had finished.
    ///
    /// The stand-in `exec`s `sleep` so the process the walker spawned is the one that hangs:
    /// `kill` reaches our child, not its descendants, which is the limitation recorded on
    /// `wait_with_deadline`.
    #[cfg(unix)]
    #[test]
    fn a_stuck_graphaligner_is_killed_and_fails_only_its_own_accession() {
        use std::os::unix::fs::PermissionsExt;
        let dir = std::env::temp_dir().join(format!("logan_align_stuck_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let fake = dir.join("FakeGraphAligner");
        std::fs::write(&fake, "#!/bin/sh\nexec sleep 30\n").unwrap();
        std::fs::set_permissions(&fake, std::fs::Permissions::from_mode(0o755)).unwrap();

        let cfg = Config {
            graphaligner_bin: fake.to_string_lossy().into_owned(),
            graphaligner_timeout_secs: 1,
            ..Config::default()
        };
        let gfa = dir.join("ACC.sub.gfa");
        std::fs::write(&gfa, b"H\tVN:Z:1.0\nS\t1\tACGTACGTAAAA\n").unwrap();
        let sub = Subgraph {
            accession: "ACC".into(), nodes: vec![1], seqs: vec![b"ACGTACGTAAAA".to_vec()],
            gfa_path: gfa, n_edges: 0, owns_gfa: true,
        };
        let q = Query {
            name: "q".into(), record_id: "pol".into(), seq: b"ACGTACGTAAAA".to_vec(),
            path: dir.join("q.fa"), n_records: 1,
        };
        let gaf = dir.join("ACC.sub.gaf");

        let started = std::time::Instant::now();
        // `GafPath` is not `Debug`, so the Ok side cannot be unwrapped away.
        let err = match run_stand_in(&cfg, &sub, std::slice::from_ref(&q), &gaf) {
            Err(e) => e,
            Ok(map) => panic!("expected a timeout, got {} rows", map.len()),
        };
        let waited = started.elapsed();

        let Error::GraphAligner { accession, detail } = &err else {
            panic!("expected a GraphAligner error, got {err}");
        };
        assert_eq!(accession, "ACC");
        assert!(detail.contains("timed out after 1 s"), "{detail}");
        assert!(detail.contains("--graphaligner-timeout"), "{detail}");
        // The 30 s sleep must not have been waited out, and the poll must not have slept
        // past the deadline either.
        assert!(waited < std::time::Duration::from_secs(5), "waited {waited:?}");
        assert!(!gaf.exists(), "a partial GAF was left behind");
        assert!(!dir.join("ACC.sub.queries.fa").exists(), "the temp query was left behind");

        drop(sub);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// 0 means wait, which is what every release before 0.2.0 did. The stand-in takes long
    /// enough that a deadline of 1 s would have killed it.
    #[cfg(unix)]
    #[test]
    fn a_zero_timeout_waits() {
        use std::os::unix::fs::PermissionsExt;
        let dir = std::env::temp_dir().join(format!("logan_align_nodl_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let fake = dir.join("FakeGraphAligner");
        std::fs::write(&fake,
            "#!/bin/sh\n\
             while [ $# -gt 0 ]; do case \"$1\" in -a) gaf=\"$2\"; shift 2;; *) shift;; esac; done\n\
             sleep 2\n\
             printf 'pol\\t12\\t0\\t12\\t+\\t>1\\t12\\t0\\t12\\t12\\t12\\t60\\n' > \"$gaf\"\n")
            .unwrap();
        std::fs::set_permissions(&fake, std::fs::Permissions::from_mode(0o755)).unwrap();

        let cfg = Config {
            graphaligner_bin: fake.to_string_lossy().into_owned(),
            graphaligner_timeout_secs: 0,
            ..Config::default()
        };
        let gfa = dir.join("ACC.sub.gfa");
        std::fs::write(&gfa, b"H\tVN:Z:1.0\nS\t1\tACGTACGTAAAA\n").unwrap();
        let sub = Subgraph {
            accession: "ACC".into(), nodes: vec![1], seqs: vec![b"ACGTACGTAAAA".to_vec()],
            gfa_path: gfa, n_edges: 0, owns_gfa: true,
        };
        let q = Query {
            name: "q".into(), record_id: "pol".into(), seq: b"ACGTACGTAAAA".to_vec(),
            path: dir.join("q.fa"), n_records: 1,
        };
        let got = run_stand_in(&cfg, &sub, std::slice::from_ref(&q), &dir.join("ACC.sub.gaf")).unwrap();
        assert!(got.contains_key("pol"));

        drop(sub);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Waiting with a deadline means reading stderr ourselves, which `Command::output` used to
    /// do. A pipe nobody reads fills at 64 KiB and blocks the writer there, so without the
    /// drain thread a chatty aligner would stop dead and then be reported as timed out: a
    /// stall of our own making, attributed to GraphAligner. The stand-in writes ~256 KiB.
    #[cfg(unix)]
    #[test]
    fn a_chatty_aligner_does_not_deadlock_on_a_full_stderr_pipe() {
        use std::os::unix::fs::PermissionsExt;
        let dir = std::env::temp_dir().join(format!("logan_align_chatty_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let fake = dir.join("FakeGraphAligner");
        std::fs::write(&fake,
            "#!/bin/sh\n\
             i=0\n\
             while [ $i -lt 4000 ]; do\n\
               echo 'progress progress progress progress progress progress progress' >&2\n\
               i=$((i+1))\n\
             done\n\
             echo 'the last thing it said' >&2\n\
             exit 1\n")
            .unwrap();
        std::fs::set_permissions(&fake, std::fs::Permissions::from_mode(0o755)).unwrap();

        let cfg = Config {
            graphaligner_bin: fake.to_string_lossy().into_owned(),
            // Short enough that a deadlock shows up as a timeout rather than as a hung test.
            graphaligner_timeout_secs: 20,
            ..Config::default()
        };
        let gfa = dir.join("ACC.sub.gfa");
        std::fs::write(&gfa, b"H\tVN:Z:1.0\nS\t1\tACGTACGTAAAA\n").unwrap();
        let sub = Subgraph {
            accession: "ACC".into(), nodes: vec![1], seqs: vec![b"ACGTACGTAAAA".to_vec()],
            gfa_path: gfa, n_edges: 0, owns_gfa: true,
        };
        let q = Query {
            name: "q".into(), record_id: "pol".into(), seq: b"ACGTACGTAAAA".to_vec(),
            path: dir.join("q.fa"), n_records: 1,
        };
        let err = match run_stand_in(&cfg, &sub, std::slice::from_ref(&q), &dir.join("ACC.sub.gaf")) {
            Err(e) => e,
            Ok(map) => panic!("expected the stand-in's exit 1, got {} rows", map.len()),
        };
        let Error::GraphAligner { detail, .. } = &err else {
            panic!("expected a GraphAligner error, got {err}");
        };
        assert!(!detail.contains("timed out"), "the drain thread is not keeping up: {detail}");
        // The whole pipe was read, so the tail is the end of it and not the first 64 KiB.
        assert!(detail.contains("the last thing it said"), "{detail}");

        drop(sub);
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The query file is read once, at startup, so a query that can only be read once still
    /// works (WA-01, which shares its cause with WA-07 and WA-08). A FIFO stands in for the
    /// `<(zcat ref.fa.gz)` and `/dev/stdin` forms `README.md` invites: the second read of one
    /// blocks or returns nothing, and the byte-copy path did exactly that on every accession
    /// after the first, so every accession came back `no_path` and the run reported a clean
    /// biological negative.
    #[cfg(unix)]
    #[test]
    fn a_query_that_can_only_be_read_once_still_aligns_every_accession() {
        use std::os::unix::fs::PermissionsExt;
        let dir = std::env::temp_dir().join(format!("logan_align_once_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let fake = dir.join("FakeGraphAligner");
        std::fs::write(&fake, "#!/bin/sh\n\
             while [ $# -gt 0 ]; do case \"$1\" in\n\
             -f) name=$(head -1 \"$2\" | cut -c2-); shift 2;;\n\
             -a) gaf=\"$2\"; shift 2;;\n\
             *) shift;;\n\
             esac; done\n\
             printf '%s\\t12\\t0\\t12\\t+\\t>1\\t12\\t0\\t12\\t12\\t12\\t60\\n' \"$name\" > \"$gaf\"\n")
            .unwrap();
        std::fs::set_permissions(&fake, std::fs::Permissions::from_mode(0o755)).unwrap();

        // A FIFO: `read_query` drains it, and anything that reopens it afterwards hangs.
        let fifo = dir.join("piped.fa");
        assert!(std::process::Command::new("mkfifo").arg(&fifo).status().unwrap().success());
        let writer = {
            let fifo = fifo.clone();
            std::thread::spawn(move || {
                let _ = std::fs::write(&fifo, b">pol\nACGTACGTAAAA\n");
            })
        };
        let q = crate::io::input::read_query(&fifo).unwrap();
        writer.join().unwrap();
        // Nothing may reopen it from here on; make sure of it.
        std::fs::remove_file(&fifo).unwrap();

        let cfg = Config { graphaligner_bin: fake.to_string_lossy().into_owned(), ..Config::default() };
        let gfa = dir.join("ACC.sub.gfa");
        std::fs::write(&gfa, b"H\tVN:Z:1.0\nS\t1\tACGTACGTAAAA\n").unwrap();
        let sub = Subgraph {
            accession: "ACC".into(), nodes: vec![1], seqs: vec![b"ACGTACGTAAAA".to_vec()],
            gfa_path: gfa, n_edges: 0, owns_gfa: true,
        };

        // Three accessions' worth of calls, as a run would make, plus the prune realignment.
        for i in 0..4 {
            let gaf = dir.join(format!("ACC{i}.sub.gaf"));
            let got = run_stand_in(&cfg, &sub, std::slice::from_ref(&q), &gaf).unwrap();
            assert!(got.contains_key("pol"), "call {i} found nothing");
        }

        drop(sub);
        let _ = std::fs::remove_dir_all(&dir);
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
        // `write_gfa` no longer creates a missing scratch directory (WF-17),
        // so a test that names its own creates it, as a run does.
        cfg.scratch_dir = std::env::temp_dir().join(format!("logan_align_{}", std::process::id()));
        std::fs::create_dir_all(&cfg.scratch_dir).unwrap();
        if resolve_binary(&cfg).is_err() {
            eprintln!("SKIP: GraphAligner not on PATH");
            return;
        }

        let f = FetchedUnitigs {
            index: 0,
            accession: "ERR2028820".into(),
            compressed_len: compressed.len(),
            compressed: crate::types::Compressed::Memory(compressed),
            fetch_secs: 0.0,
            payload: crate::types::Payload::default(),
        };
        let graph = { let rd = stream::decoder(&f).unwrap(); parse::build(&f.accession, rd, u64::MAX).unwrap() };
        let query = read_query(Path::new("../queries/pfcrt_only.fasta")).unwrap();
        let idx = anchor::index_query(&cfg, &query).unwrap();
        let a = anchor::find(&cfg, &graph, &idx).unwrap();
        let nodes = extract::expand(&cfg, &graph, &a).unwrap();
        let seqs = extract::collect_seqs(&graph, &nodes).unwrap();
        let sub = extract::write_gfa(&cfg, &graph, nodes, seqs).unwrap();

        let gaf = sub.gfa_path.with_extension("gaf");
        let t0 = std::time::Instant::now();
        let mut gaf_map = run(&cfg, &sub, std::slice::from_ref(&query), &gaf).expect("align");
        let path = gaf_map.remove(query.record_id.as_str()).expect("a path");
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
        // `write_gfa` no longer creates a missing scratch directory (WF-17),
        // so a test that names its own creates it, as a run does.
        cfg.scratch_dir = std::env::temp_dir().join(format!("logan_align_e_{}", std::process::id()));
        std::fs::create_dir_all(&cfg.scratch_dir).unwrap();
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
            owns_gfa: true,
        };
        let q = Query {
            name: "q".into(), record_id: "q".into(), seq: b"ACGTACGT".to_vec(),
            path: "../queries/pfcrt_only.fasta".into(),
            n_records: 1,
        };
        let got = run(&cfg, &sub, std::slice::from_ref(&q), &cfg.scratch_dir.join("EMPTY.sub.gaf")).unwrap();
        assert!(got.is_empty(), "must report no path rather than invoking the tool");
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
            // `write_gfa` no longer creates a missing scratch directory (WF-17),
            // so a test that names its own creates it, as a run does.
            cfg.scratch_dir = std::env::temp_dir().join(format!("logan_det_{}", std::process::id()));
            std::fs::create_dir_all(&cfg.scratch_dir).unwrap();
            if resolve_binary(&cfg).is_err() {
                eprintln!("SKIP: GraphAligner not on PATH");
                return;
            }
            let f = FetchedUnitigs {
                index: 0,
                accession: acc.into(),
                compressed_len: compressed.len(),
                compressed: crate::types::Compressed::Memory(compressed),
                fetch_secs: 0.0,
                payload: crate::types::Payload::default(),
            };
            let graph = {
                let rd = stream::decoder(&f).unwrap();
                parse::build(&f.accession, rd, u64::MAX).unwrap()
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
                let mut gaf_map = run(&cfg, &sub, std::slice::from_ref(&query), &gaf).unwrap();
                let path = gaf_map.remove(query.record_id.as_str()).expect("an alignment");
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
