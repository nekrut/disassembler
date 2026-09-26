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
mod stages;
mod types;

use config::Config;
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
        --min-shared-kmers <N>    Anchor threshold, validated 1..=5 [default: 3]
        --max-size-mb <N>         Maximum compressed SRA size in MB [default: 1024]
        --prune-frac <F>          Drop unitigs below this fraction of the gene's depth and
                                  realign, removing the aligner's reference bias. 0.0
                                  disables it [default: 0.01]
        --alt-max-detour <N>      Unitigs an alternative allele may traverse before it must
                                  rejoin the main route [default: 100]
    -t, --threads <N>             Threads per GraphAligner invocation [default: 1; above 1
                                  is not reproducible, see Config]
    -w, --workers <N>             Compute workers
        --fetch-workers <N>       Fetch workers [default: 8]

INTRA-HOST VARIANTS:
        --no-prefer-majority      Leave the primary on the reference-matching allele at
                                  two-allele sites, instead of the more abundant one
        --no-emit-alts            Do not write alt_haplotypes.fasta
        --emit-vcf                Also write <out-dir>/<acc>.vcf: intra-host SNPs with their
                                  ka:f frequency, in query coordinates [default: off]
        --emit-snv-table          Also write <out-dir>/snvs.tsv: one long table of every
                                  intra-host SNP across the run, one row per (accession, site,
                                  allele). The cohort output that scales; load into DuckDB
                                  [default: off]
        --no-indel-alts           Substitution alleles only; refuse any that change length
        --allow-frameshift-alts   Also accept indel alleles whose length change is not a
                                  multiple of 3. These reframe the CDS; off by default
        --max-alt-indel <N>       Largest length change still treated as an indel allele
                                  rather than a rerouting [default: 30]

        --allow-empty             Exit 0 even if no accession reconstructed. Without it, a
                                  run that reconstructs nothing exits non-zero, so a wrapper
                                  cannot mistake an all-failed shard for a good one

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
            "--min-shared-kmers" => {
                args.cfg.min_shared_kmers = value()?.parse().map_err(|e| error::Error::Usage(format!("--min-shared-kmers: {e}")))?;
            }
            "--max-size-mb" => {
                let mb: u64 = value()?.parse().map_err(|e| error::Error::Usage(format!("--max-size-mb: {e}")))?;
                // saturate rather than wrap: an absurd value should mean "no cap", not a tiny one.
                args.cfg.max_compressed_bytes = mb.saturating_mul(1024 * 1024);
            }
            "-t" | "--threads" => {
                args.cfg.graphaligner_threads = value()?.parse().map_err(|e| error::Error::Usage(format!("--threads: {e}")))?;
            }
            "-w" | "--workers" => {
                args.cfg.compute_workers = value()?.parse().map_err(|e| error::Error::Usage(format!("--workers: {e}")))?;
            }
            "--fetch-workers" => {
                args.cfg.fetch_workers = value()?.parse().map_err(|e| error::Error::Usage(format!("--fetch-workers: {e}")))?;
            }
            "--scratch-dir" => {
                args.cfg.scratch_dir = PathBuf::from(value()?);
            }
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

fn main() -> ExitCode {
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
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}

fn real_main(args: Args) -> error::Result<()> {
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
    // The single-file/single-record legacy path stays byte-identical (uses `read_query`,
    // which fires the tools-iuc #8388 warning). Any wider intake uses `read_query_set`.
    let queries: Option<Vec<crate::types::Query>> = if args.cfg.fetch_only {
        None
    } else {
        if args.query.is_empty() {
            return Err(error::Error::Usage("--query is required".into()));
        }
        // Try the legacy single-file path first for its warning + n_records semantics.
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
        (Some(a), _) => {
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
        "[logan-walker] memory bound: up to {} compressed buffers resident \
         ({} queued + {} fetching + {} parsing) x {} MB cap = {:.1} GB worst case, \
         plus ~43 MB per compute worker for the CSR graph",
        args.cfg.peak_compressed_buffers(),
        args.cfg.channel_capacity,
        args.cfg.fetch_workers,
        args.cfg.compute_workers,
        args.cfg.max_compressed_bytes / (1024 * 1024),
        gb(args.cfg.max_compressed_resident_bytes()),
    );

    // --fetch-only takes the short path: fill the cache and exit. No compute pool, no
    // output files, no query. Success = every accession is in the cache after this.
    if args.cfg.fetch_only {
        let (ok, fail) = pipeline::run_fetch_only(&args.cfg, accessions)?;
        eprintln!("[logan-walker --fetch-only] {ok} fetched, {fail} failed");
        if fail > 0 && !args.allow_empty {
            return Err(error::Error::NothingReconstructed { total: ok + fail });
        }
        return Ok(());
    }

    // Unwrap the query set for the full pipeline. The guard above ensures it is Some here.
    let queries = queries.expect("queries is Some when not --fetch-only");

    let summary = pipeline::run(&args.cfg, accessions, &queries)?;

    // Downstream writers that still take a single reference use the first query in the pool.
    // In multi-anchor mode this is only used for coordinate contexts where the collector
    // records the winning reference per haplotype separately.
    let query = &queries[0];

    // ---- OUT of the pipeline ----------------------------------------------
    std::fs::create_dir_all(&args.out_dir)?;
    io::output::write_fasta(&summary, &args.out_fasta)?;
    io::output::write_bubbles(&summary, &args.out_dir.join("bubbles.json"))?;
    // One alignment pass per accession, shared by both writers.
    if args.cfg.emit_vcf || args.cfg.emit_snv_table {
        let sites = io::output::snv_sites_all(&summary, query);
        if args.cfg.emit_vcf {
            io::output::write_vcfs(&summary, query, &args.out_dir, &sites)?;
        }
        if args.cfg.emit_snv_table {
            io::output::write_snv_table(
                &summary,
                query,
                &args.out_dir.join("snvs.tsv"),
                &sites,
            )?;
        }
    }
    if args.cfg.emit_alts {
        io::output::write_alts(&summary, &args.out_dir.join("alt_haplotypes.fasta"))?;
    }
    io::output::write_summary(
        &summary,
        &args.out_dir.join(format!("{}_walker_summary.json", query.name)),
    )?;

    println!(
        "[logan-walker] {}/{} accessions reconstructed",
        summary.success, summary.total
    );

    // Per-accession failures are collected rather than propagated, which is what lets one bad
    // accession not sink a 10,000-accession shard. The cost is that a run where *every*
    // accession failed used to write a 0-byte FASTA and exit 0, which any scheduler reads as
    // a good shard. Outputs are already written above, so they stay on disk for inspection.
    if summary.success == 0 && summary.total > 0 && !args.allow_empty {
        return Err(error::Error::NothingReconstructed { total: summary.total });
    }
    Ok(())
}
