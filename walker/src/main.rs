//! logan-walker: Logan unitigs to coding haplotypes.
//!
//! ```text
//!   main.rs        CLI, config, calls pipeline::run
//!   pipeline.rs    threading structure and stage wiring   <- read this for the flow
//!   config.rs      tunables
//!   types.rs       what moves between stages              <- read this for the shapes
//!   error.rs       error type
//!   io/input.rs    accession list, reference CDS          <- INTO the pipeline
//!   io/output.rs   fasta, vcf, snv table, bubbles, summary <- OUT of the pipeline
//!   stages/stream.rs    S3 + zstd
//!   stages/parse.rs     BCALM2 headers -> CSR
//!   stages/anchor.rs    canonical k-mers (replaces BLAST)
//!   stages/extract.rs   5-hop BFS + GFA + prune
//!   stages/align.rs     GraphAligner wrapper (external C++ binary)
//!   stages/stitch.rs    path -> sequence
//!   stages/variants.rs  prefer-majority + intra-host alt recovery
//!   stages/phase.rs     forks -> bubbles.json
//! ```
//!
//! Concordance gate: `bench/check_golden.sh`.

// A few items are retained for API symmetry / tests.
#![allow(dead_code, unused_variables)]
// Intentional style, not defects: tab-delimited GFA/VCF examples in doc comments; test setup
// that builds a Config::default() then tweaks a field; and a couple of wide tuples in tests.
// None touch the runtime path.
#![allow(
    clippy::tabs_in_doc_comments,
    clippy::field_reassign_with_default,
    clippy::type_complexity
)]

mod config;
mod error;
mod io;
mod pipeline;
mod scratch;
mod stages;
mod types;

use config::{Config, PrimaryMetric};
use std::path::PathBuf;
use std::process::ExitCode;

const USAGE: &str = "\
logan-walker

USAGE:
    logan-walker [OPTIONS]

INPUT:
    -a, --accession <ACC>         Single SRA accession
    -l, --accession-list <FILE>   File with a list of SRA accessions
    -q, --query <FASTA>           Reference CDS FASTA. Repeatable, and each file may hold
                                  multiple records; two or more references switches on
                                  multi-anchor mode (union anchor pool, primary is the
                                  reference whose alignment gives the best cov/pident)

OUTPUT:
    -o, --out-dir <DIR>           Working directory
    -f, --out-fasta <FILE>        Reconstructed haplotypes FASTA

TUNING:
    -k, --k <N>                   Anchor k-mer size, 1..=32 [default: 21]
        --hops <N>                Hop radius for subgraph expansion [default: 5]
        --min-abund <F>           Minimum ka:f for a variant [default: 5.0]
        --min-vaf <F>             Minimum VAF for an alternative branch [default: 0.10]
        --min-abund-frac <F>      Raise the variant floor to this fraction of gene depth,
                                  median ka:f on the walked path [default: 0.0 = off]
        --min-cov <F>             Minimum reference coverage [default: 0.85]
        --primary-metric <M>      Multi-anchor primary selection: matched | cov | pident |
                                  cov_x_pident [default: matched, the fraction of the
                                  reference's bases matched]. Ignored in single-anchor mode
        --min-shared-kmers <N>    Anchor threshold, validated 1..=5 [default: 3]
        --max-size-mb <N>         Maximum compressed SRA size in MB [default: 1024]
        --prune-frac <F>          Drop unitigs below this fraction of the gene's depth and
                                  realign, removing the aligner's reference bias. The
                                  realignment is kept only if it spans as much of the
                                  reference. 0.0 disables it [default: 0.01]
        --alt-max-detour <N>      Unitigs an alternative allele may traverse before it must
                                  rejoin the main route [default: 100]
    -t, --threads <N>             Threads per GraphAligner invocation, at least 1 [default: 1;
                                  above 1 is not reproducible, see Config]
        --graphaligner <PATH>     GraphAligner executable [default: GraphAligner on PATH].
                                  Checked once at startup; results are validated against
                                  1.0.17b
    -w, --workers <N>             Compute workers, at least 1 [default: 8]
        --fetch-workers <N>       Fetch workers, at least 1 [default: 8]

INTRA-HOST VARIANTS:
        --no-prefer-majority      Leave the primary on the reference-matching allele at
                                  two-allele sites, instead of the more abundant one
        --no-emit-alts            Do not write alt_haplotypes.fasta
        --emit-vcf                Also write <out-dir>/<acc>.vcf: intra-host SNPs with their
                                  ka:f frequency, in the coordinates of the reference the
                                  primary was reconstructed against [default: off]
        --emit-snv-table          Also write <out-dir>/snvs.tsv: one long table of every
                                  intra-host SNP across the run, one row per (accession, site,
                                  allele). The cohort output that scales; load into DuckDB
                                  [default: off]
        --no-indel-alts           Substitution alleles only; refuse any that change length
        --allow-frameshift-alts   Also accept indel alleles whose length change is not a
                                  multiple of 3. These reframe the CDS; off by default
        --max-alt-indel <N>       Largest length change still treated as an indel allele
                                  rather than a rerouting [default: 30]

        --allow-empty             Exit 0 when no accession reconstructed because nothing was
                                  found: no Logan data, over --max-size-mb, no anchors or no
                                  path. Errors still fail the run. Without it, a run that
                                  reconstructs nothing exits non-zero, so a wrapper cannot
                                  mistake an empty shard for a good one

SCRATCH:
        --scratch-dir <DIR>       Where the run's private scratch directory is made; it is
                                  removed when the run ends or is interrupted [default:
                                  /dev/shm where it exists, else the system temp directory]

OFFLINE / CACHE:
        --unitigs-cache <DIR>     Read-through cache for compressed unitig payloads.
                                  <DIR>/<ACC>.unitigs.fa.zst is used if present; otherwise
                                  fetched from S3 and written back atomically. Byte-identical
                                  to the S3 object
        --offline                 Refuse the S3 fallback. A cache miss becomes a hard error
                                  naming the missing accession. Requires --unitigs-cache
        --fetch-only              Fetch every accession into --unitigs-cache and exit before
                                  the compute stages. Requires --unitigs-cache; -q is not
                                  required in this mode

    -h, --help                    Print this help
    -V, --version                 Print version

EXIT STATUS:
    0   at least one accession reconstructed, or (with --allow-empty) none did because
        nothing was found
    1   a usage or I/O error, or nothing reconstructed
    3   every accession failed with an error (fetch, decode, GraphAligner): the run did
        not work. --allow-empty never turns this into 0
";

struct Args {
    accession: Option<String>,
    accession_list: Option<PathBuf>,
    /// One or more reference FASTAs. A single-file, single-record invocation stays on the
    /// legacy path (byte-identical output, multi-record warning intact). Anything wider
    /// switches the walker into multi-anchor mode.
    query: Vec<PathBuf>,
    out_dir: PathBuf,
    out_fasta: PathBuf,
    /// `--scratch-dir`: where the run's private scratch directory is made (see scratch.rs).
    scratch_base: Option<PathBuf>,
    /// Exit 0 even when nothing was reconstructed. Off by default: a shard that
    /// reconstructs nothing is a failure a scheduler needs to see.
    allow_empty: bool,
    cfg: Config,
}

fn parse_args(argv: &[String]) -> error::Result<Args> {
    let mut args = Args {
        accession: None,
        accession_list: None,
        query: Vec::new(),
        out_dir: PathBuf::from("walker_output"),
        out_fasta: PathBuf::from("reconstructed_haplotypes.fasta"),
        scratch_base: None,
        allow_empty: false,
        cfg: Config::default(),
    };
    let mut i = 0;
    while i < argv.len() {
        let flag = argv[i].as_str();
        let value = || -> error::Result<String> {
            argv.get(i + 1)
                .cloned()
                .ok_or_else(|| error::Error::Usage(format!("{flag} needs a value")))
        };
        match flag {
            "-a" | "--accession" => args.accession = Some(value()?),
            "-l" | "--accession-list" => args.accession_list = Some(PathBuf::from(value()?)),
            "-q" | "--query" => args.query.push(PathBuf::from(value()?)),
            "-o" | "--out-dir" => args.out_dir = PathBuf::from(value()?),
            "-f" | "--out-fasta" => args.out_fasta = PathBuf::from(value()?),
            "-k" | "--k" => {
                args.cfg.k = value()?.parse().map_err(|e| error::Error::Usage(format!("--k: {e}")))?;
            }
            "--hops" => {
                args.cfg.hops = value()?.parse().map_err(|e| error::Error::Usage(format!("--hops: {e}")))?;
            }
            "--min-abund" => {
                args.cfg.min_abund = value()?.parse().map_err(|e| error::Error::Usage(format!("--min-abund: {e}")))?;
            }
            "--min-vaf" => {
                args.cfg.min_vaf = value()?.parse().map_err(|e| error::Error::Usage(format!("--min-vaf: {e}")))?;
            }
            "--min-cov" => {
                args.cfg.min_cov = value()?.parse().map_err(|e| error::Error::Usage(format!("--min-cov: {e}")))?;
            }
            "--primary-metric" => {
                let raw = value()?;
                args.cfg.primary_metric = PrimaryMetric::parse(&raw).ok_or_else(|| {
                    error::Error::Usage(format!(
                        "--primary-metric: expected matched | cov | pident | cov_x_pident, got {raw:?}"
                    ))
                })?;
            }
            "--min-shared-kmers" => {
                args.cfg.min_shared_kmers = value()?.parse().map_err(|e| error::Error::Usage(format!("--min-shared-kmers: {e}")))?;
            }
            "--max-size-mb" => {
                let mb: u64 = value()?.parse().map_err(|e| error::Error::Usage(format!("--max-size-mb: {e}")))?;
                // saturate rather than wrap: an absurd value should mean "no cap", not a tiny one.
                args.cfg.max_compressed_bytes = mb.saturating_mul(1024 * 1024);
            }
            "-t" | "--threads" => args.cfg.graphaligner_threads = at_least_one("--threads", &value()?)?,
            "-w" | "--workers" => args.cfg.compute_workers = at_least_one("--workers", &value()?)?,
            "--fetch-workers" => args.cfg.fetch_workers = at_least_one("--fetch-workers", &value()?)?,
            "--scratch-dir" => args.scratch_base = Some(PathBuf::from(value()?)),
            "--graphaligner" => args.cfg.graphaligner_bin = value()?,
            // The prune sweep is the reason this flag exists, so it validates its range
            // rather than just parsing.
            // A boolean, so it consumes no value and must not fall through to `value()`.
            "--alt-max-detour" => {
                args.cfg.alt_max_detour = value()?.parse().map_err(|e| {
                    error::Error::Usage(format!("--alt-max-detour: {e}"))
                })?;
            }
            // Both are on by default, so the useful switch is the one that turns them off.
            // The positive forms stay accepted so existing command lines keep working.
            "--prefer-majority" | "--no-prefer-majority" => {
                args.cfg.prefer_majority = flag == "--prefer-majority";
                i += 1;
                continue;
            }
            "--no-indel-alts" => {
                args.cfg.allow_indel_alts = false;
                i += 1;
                continue;
            }
            "--allow-frameshift-alts" => {
                args.cfg.allow_frameshift_alts = true;
                i += 1;
                continue;
            }
            "--max-alt-indel" => {
                args.cfg.max_alt_indel = value()?.parse().map_err(|e| {
                    error::Error::Usage(format!("--max-alt-indel: {e}"))
                })?;
            }
            "--allow-empty" => {
                args.allow_empty = true;
                i += 1;
                continue;
            }
            "--unitigs-cache" => {
                args.cfg.unitigs_cache = Some(PathBuf::from(value()?));
            }
            "--offline" => {
                args.cfg.offline = true;
                i += 1;
                continue;
            }
            "--fetch-only" => {
                args.cfg.fetch_only = true;
                i += 1;
                continue;
            }
            "--emit-alts" | "--no-emit-alts" => {
                args.cfg.emit_alts = flag == "--emit-alts";
                i += 1;
                continue;
            }
            "--emit-vcf" | "--no-emit-vcf" => {
                args.cfg.emit_vcf = flag == "--emit-vcf";
                i += 1;
                continue;
            }
            "--emit-snv-table" | "--no-emit-snv-table" => {
                args.cfg.emit_snv_table = flag == "--emit-snv-table";
                i += 1;
                continue;
            }
            "--min-abund-frac" => {
                let raw = value()?;
                let v: f64 = raw.parse().map_err(|_| {
                    error::Error::Usage(format!("--min-abund-frac needs a number, got {raw:?}"))
                })?;
                if !(0.0..1.0).contains(&v) {
                    return Err(error::Error::Usage(format!(
                        "--min-abund-frac must be in [0.0, 1.0), got {v}"
                    )));
                }
                args.cfg.min_abund_frac = v;
            }
            "--prune-frac" => {
                let raw = value()?;
                let v: f64 = raw.parse().map_err(|_| {
                    error::Error::Usage(format!("--prune-frac needs a number, got {raw:?}"))
                })?;
                if !(0.0..1.0).contains(&v) {
                    return Err(error::Error::Usage(format!(
                        "--prune-frac must be in [0.0, 1.0), got {v}"
                    )));
                }
                args.cfg.prune_frac = v;
            }
            other => return Err(error::Error::Usage(format!("unknown flag {other}"))),
        }
        i += 2;
    }
    Ok(args)
}

/// A count that must be at least 1. Zero compute workers left nothing to drain the fetch
/// queue, so `-w 0` hung on six accessions and silently dropped three.
fn at_least_one(flag: &str, raw: &str) -> error::Result<usize> {
    match raw.parse::<usize>() {
        Ok(0) => Err(error::Error::Usage(format!("{flag} must be at least 1, got 0"))),
        Ok(n) => Ok(n),
        Err(e) => Err(error::Error::Usage(format!("{flag}: {e}"))),
    }
}

/// End quietly when stdout is closed, as other Unix tools do. Rust ignores SIGPIPE, so a
/// write to a closed pipe (`logan-walker ... | head -1`) returns EPIPE, and `println!`
/// panics on it.
#[cfg(unix)]
fn default_sigpipe() {
    // SAFETY: called once at startup, before any thread is spawned.
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
}

fn main() -> ExitCode {
    #[cfg(unix)]
    default_sigpipe();
    let argv: Vec<String> = std::env::args().skip(1).collect();

    if argv.is_empty() || argv.iter().any(|a| a == "-h" || a == "--help") {
        print!("{USAGE}");
        return ExitCode::SUCCESS;
    }
    if argv.iter().any(|a| a == "-V" || a == "--version") {
        println!("logan-walker {}", env!("CARGO_PKG_VERSION"));
        return ExitCode::SUCCESS;
    }

    let args = match parse_args(&argv) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };

    if let Err(e) = real_main(args) {
        eprintln!("error: {e}");
        return ExitCode::from(e.exit_code());
    }
    ExitCode::SUCCESS
}

/// Whether a finished run counts as a success. One reconstruction is enough: a bad
/// accession must not sink a 10,000-accession shard. With none, the run fails, unless
/// `allow_empty` and every accession was a negative (nothing found, not an error). A run
/// in which every accession errored fails with its own exit code, whatever `allow_empty`.
fn check_run(results: &[types::AccessionResult], allow_empty: bool) -> error::Result<()> {
    let total = results.len();
    let success = results.iter().filter(|r| r.outcome == types::Outcome::Success).count();
    let errors = results.iter().filter(|r| r.outcome == types::Outcome::Failed).count();
    if total == 0 || success > 0 {
        return Ok(());
    }
    if errors == total {
        return Err(error::Error::AllErrored { total });
    }
    if errors > 0 || !allow_empty {
        return Err(error::Error::NothingReconstructed { total, errors });
    }
    Ok(())
}

/// `--fetch-only` succeeds when every accession is in the cache. `allow_empty` accepts the
/// ones Logan has no data for or that are over the size cap, never an error.
fn check_fetch(ok: usize, negatives: usize, errors: usize, allow_empty: bool) -> error::Result<()> {
    let total = ok + negatives + errors;
    if errors > 0 && errors == total {
        return Err(error::Error::AllErrored { total });
    }
    if errors > 0 || (negatives > 0 && !allow_empty) {
        return Err(error::Error::FetchIncomplete { total, errors, negatives });
    }
    Ok(())
}

fn real_main(mut args: Args) -> error::Result<()> {
    // ---- validate offline / cache flag combinations -----------------------
    if args.cfg.offline && args.cfg.unitigs_cache.is_none() {
        return Err(error::Error::Usage(
            "--offline requires --unitigs-cache <DIR>".into(),
        ));
    }
    if args.cfg.fetch_only && args.cfg.unitigs_cache.is_none() {
        return Err(error::Error::Usage(
            "--fetch-only requires --unitigs-cache <DIR>".into(),
        ));
    }

    // ---- INTO the pipeline ------------------------------------------------
    // --fetch-only does not need a query: the compute stages never run.
    //
    // Multi-anchor is auto-enabled when the collected reference set has more than one
    // record, from any combination: several -q files, one multi-record file, or both.
    // The single-file/single-record legacy path stays byte-identical (uses `read_query`).
    // Any wider intake uses `read_query_set`, and is announced below.
    let queries: Option<Vec<crate::types::Query>> = if args.cfg.fetch_only {
        None
    } else {
        if args.query.is_empty() {
            return Err(error::Error::Usage("--query is required".into()));
        }
        // Try the legacy single-file path first for its n_records semantics.
        // If that file turns out to hold more than one record, re-load through
        // read_query_set so the whole pool goes into anchor union / align fan-out.
        Some(if args.query.len() == 1 {
            let first = io::input::read_query(&args.query[0])?;
            if first.n_records > 1 {
                io::input::read_query_set(&args.query)?
            } else {
                vec![first]
            }
        } else {
            io::input::read_query_set(&args.query)?
        })
    };

    let accessions = match (&args.accession, &args.accession_list) {
        // `-a` used to win silently, so a run meant for a whole list did one accession.
        (Some(_), Some(_)) => {
            return Err(error::Error::Usage(
                "give --accession or --accession-list, not both".into(),
            ))
        }
        (Some(a), None) => {
            // The list path has always validated; `-a` did not, though the value reaches the
            // same URL and the same two filesystem paths.
            io::input::validate_accession(a, "--accession")?;
            vec![a.clone()]
        }
        (None, Some(p)) => io::input::read_accession_list(p)?,
        (None, None) => {
            return Err(error::Error::Usage(
                "provide --accession or --accession-list".into(),
            ))
        }
    };

    // ---- the pipeline -----------------------------------------------------
    // State the memory bound up front. This used to be a comment in `pipeline`, and it was
    // wrong by 5x; a number that moves with -w, --fetch-workers and --max-size-mb is worth
    // printing rather than documenting.
    let gb = |b: u64| b as f64 / (1024.0 * 1024.0 * 1024.0);
    eprintln!(
        "[logan-walker] memory bound: {} compressed buffers ({} queued + {} fetching + {} \
         parsing) x {} MB cap = {:.1} GB, plus {} parsed graphs at up to {:.1}x their \
         compressed size = {:.1} GB; {:.1} GB worst case",
        args.cfg.peak_compressed_buffers(),
        args.cfg.channel_capacity,
        args.cfg.fetch_workers,
        args.cfg.compute_workers,
        args.cfg.max_compressed_bytes / (1024 * 1024),
        gb(args.cfg.max_compressed_resident_bytes()),
        args.cfg.compute_workers,
        config::GRAPH_PEAK_PER_COMPRESSED_BYTE,
        gb(args.cfg.max_graph_resident_bytes()),
        gb(args.cfg.max_resident_bytes()),
    );

    // --fetch-only takes the short path: fill the cache and exit. No compute pool, no
    // output files, no query. Success = every accession is in the cache after this.
    if args.cfg.fetch_only {
        let (ok, negatives, errors) = pipeline::run_fetch_only(&args.cfg, accessions)?;
        eprintln!(
            "[logan-walker --fetch-only] {ok} fetched, {negatives} not available \
             (no Logan data or over --max-size-mb), {errors} failed"
        );
        return check_fetch(ok, negatives, errors, args.allow_empty);
    }

    // Unwrap the query set for the full pipeline. The guard above ensures it is Some here.
    let queries = queries.expect("queries is Some when not --fetch-only");

    // Checked before any download. GraphAligner used to be looked up per accession, after
    // its unitigs were fetched and parsed, so a missing binary failed each one in turn.
    let bin = stages::align::resolve_binary(&args.cfg)?;
    let version = stages::align::version(&bin)?;
    eprintln!("[logan-walker] GraphAligner: {version} ({})", bin.display());
    if !version.contains("1.0.17b") {
        eprintln!(
            "[warn] the walker is validated against GraphAligner 1.0.17b; 1.0.19 and later \
             clip the query's last base (docs/issue_graphaligner_end_clipping.md)"
        );
    }
    args.cfg.graphaligner_bin = bin.to_string_lossy().into_owned();

    // Removed when `scratch` goes out of scope at the end of this function, or on a signal.
    let scratch = scratch::Scratch::create(args.scratch_base.as_deref())?;
    scratch.remove_on_signal()?;
    args.cfg.scratch_dir = scratch.path().to_path_buf();
    if queries.len() > 1 {
        eprintln!(
            "[logan-walker] multi-anchor mode: {} references from {} -q file(s); each \
             accession's primary is the reference chosen by --primary-metric {}",
            queries.len(),
            args.query.len(),
            args.cfg.primary_metric.as_str()
        );
    }

    let summary = pipeline::run(&args.cfg, accessions, &queries)?;

    // ---- OUT of the pipeline ----------------------------------------------
    std::fs::create_dir_all(&args.out_dir).map_err(error::file_err("create", &args.out_dir))?;
    io::output::write_fasta(&summary, &args.out_fasta).map_err(|e| e.at("write", &args.out_fasta))?;
    let bubbles = args.out_dir.join("bubbles.json");
    io::output::write_bubbles(&summary, &bubbles).map_err(|e| e.at("write", &bubbles))?;
    // One alignment pass per accession, shared by both writers.
    if args.cfg.emit_vcf || args.cfg.emit_snv_table {
        let sites = io::output::snv_sites_all(&summary, &queries);
        if args.cfg.emit_vcf {
            io::output::write_vcfs(&summary, &queries, &args.out_dir, &sites)
                .map_err(|e| e.at("write VCFs in", &args.out_dir))?;
        }
        if args.cfg.emit_snv_table {
            let snvs = args.out_dir.join("snvs.tsv");
            io::output::write_snv_table(&summary, &queries, &snvs, &sites)
                .map_err(|e| e.at("write", &snvs))?;
        }
    }
    if args.cfg.emit_alts {
        let alts = args.out_dir.join("alt_haplotypes.fasta");
        io::output::write_alts(&summary, &alts).map_err(|e| e.at("write", &alts))?;
    }
    let summary_path = args.out_dir.join(format!("{}_walker_summary.json", queries[0].name));
    io::output::write_summary(&summary, &summary_path).map_err(|e| e.at("write", &summary_path))?;

    println!(
        "[logan-walker] {}/{} accessions reconstructed",
        summary.success, summary.total
    );

    // Per-accession failures are collected rather than propagated, which is what lets one bad
    // accession not sink a 10,000-accession shard. The cost is that a run where *every*
    // accession failed used to write a 0-byte FASTA and exit 0, which any scheduler reads as
    // a good shard; until 0.1.1, --allow-empty still did that for a run where every
    // accession errored. Outputs are already written above, so they stay on disk for
    // inspection.
    check_run(&summary.results, args.allow_empty)
}

#[cfg(test)]
mod tests {
    use super::{at_least_one, check_fetch, check_run, parse_args};
    use crate::error::Error;
    use crate::types::{AccessionResult, Outcome};

    fn run_of(outcomes: &[Outcome]) -> Vec<AccessionResult> {
        outcomes
            .iter()
            .enumerate()
            .map(|(i, &o)| AccessionResult::empty(i, format!("SRR{i}"), o, None))
            .collect()
    }

    fn argv(s: &str) -> Vec<String> {
        s.split_whitespace().map(String::from).collect()
    }

    #[test]
    fn one_success_is_enough_whatever_else_happened() {
        let r = run_of(&[Outcome::Failed, Outcome::Success, Outcome::NoAnchors]);
        assert!(check_run(&r, false).is_ok());
    }

    #[test]
    fn allow_empty_accepts_a_run_that_found_nothing() {
        let r = run_of(&[Outcome::NoAnchors, Outcome::NoPath, Outcome::NoLoganData, Outcome::TooLarge]);
        assert!(check_run(&r, true).is_ok());
        assert!(matches!(
            check_run(&r, false),
            Err(Error::NothingReconstructed { total: 4, errors: 0 })
        ));
    }

    #[test]
    fn allow_empty_does_not_hide_a_run_where_every_accession_errored() {
        let r = run_of(&[Outcome::Failed, Outcome::Failed]);
        for allow_empty in [false, true] {
            let e = check_run(&r, allow_empty).unwrap_err();
            assert!(matches!(e, Error::AllErrored { total: 2 }));
            assert_eq!(e.exit_code(), 3);
        }
    }

    #[test]
    fn allow_empty_does_not_hide_errors_mixed_with_negatives() {
        let r = run_of(&[Outcome::NoAnchors, Outcome::Failed]);
        let e = check_run(&r, true).unwrap_err();
        assert!(matches!(e, Error::NothingReconstructed { total: 2, errors: 1 }));
        assert_eq!(e.exit_code(), 1);
    }

    #[test]
    fn fetch_only_accepts_missing_data_only_with_allow_empty() {
        assert!(check_fetch(3, 0, 0, false).is_ok());
        assert!(check_fetch(2, 1, 0, true).is_ok());
        assert!(matches!(
            check_fetch(2, 1, 0, false),
            Err(Error::FetchIncomplete { total: 3, errors: 0, negatives: 1 })
        ));
        assert!(matches!(
            check_fetch(2, 0, 1, true),
            Err(Error::FetchIncomplete { total: 3, errors: 1, negatives: 0 })
        ));
        assert!(matches!(check_fetch(0, 0, 3, true), Err(Error::AllErrored { total: 3 })));
    }

    #[test]
    fn fetch_only_partial_failure_does_not_say_reconstructed() {
        let msg = check_fetch(2, 0, 1, false).unwrap_err().to_string();
        assert!(msg.starts_with("1 of 3 accessions were not cached"), "{msg}");
        assert!(!msg.contains("reconstructed"), "{msg}");
    }

    #[test]
    fn worker_and_thread_counts_must_be_at_least_one() {
        for flag in ["-w", "--workers", "--fetch-workers", "-t", "--threads"] {
            let err = parse_args(&argv(&format!("{flag} 0"))).err().expect(flag);
            assert!(err.to_string().contains("must be at least 1"), "{flag}: {err}");
        }
        let args = parse_args(&argv("-w 2 --fetch-workers 3 -t 1")).unwrap();
        assert_eq!(
            (args.cfg.compute_workers, args.cfg.fetch_workers, args.cfg.graphaligner_threads),
            (2, 3, 1)
        );
        assert!(at_least_one("-w", "two").is_err());
    }
}
