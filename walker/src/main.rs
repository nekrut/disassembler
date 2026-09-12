//! logan-walker: Logan unitigs to coding haplotypes.
//!
//! Skeleton. No algorithms are implemented; every stage returns
//! `Error::NotImplemented` naming itself, so running the binary shows how far the flow gets.
//!
//! ```text
//!   main.rs        CLI, config, calls pipeline::run
//!   pipeline.rs    threading structure and stage wiring   <- read this for the flow
//!   config.rs      tunables
//!   types.rs       what moves between stages              <- read this for the shapes
//!   error.rs       error type
//!   io/input.rs    accession list, reference CDS          <- INTO the pipeline
//!   io/output.rs   fasta, bubbles, summary                <- OUT of the pipeline
//!   stages/stream.rs    S3 + zstd
//!   stages/parse.rs     BCALM2 headers -> CSR
//!   stages/anchor.rs    canonical k-mers (replaces BLAST)
//!   stages/extract.rs   5-hop BFS + GFA
//!   stages/align.rs     GraphAligner wrapper (external C++ binary)
//!   stages/stitch.rs    path -> sequence
//!   stages/phase.rs     forks -> variants
//! ```
//!
//! Plan: `work/plan_rust_path_walker.md`. Concordance gate: `work/benchmark/verify_golden.py`.

// Templating stage: unused items are expected until the stages are filled in.
#![allow(dead_code, unused_variables)]
// Intentional style, not defects: tab-delimited GFA/VCF examples in doc comments; test setup
// that builds a Config::default() then tweaks a field; and a couple of wide tuples in tests and
// in stages::detour (which is slated for removal). None touch the runtime path.
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
logan-walker (skeleton)

USAGE:
    logan-walker [OPTIONS]

INPUT:
    -a, --accession <ACC>         Single SRA accession
    -l, --accession-list <FILE>   File with a list of SRA accessions
    -q, --query <FASTA>           Reference CDS FASTA

OUTPUT:
    -o, --out-dir <DIR>           Working directory
    -f, --out-fasta <FILE>        Reconstructed haplotypes FASTA

TUNING:
    -k, --k <N>                   Anchor k-mer size, 1..=32 [default: 21]
    -H, --haplotypes <one|many>   Haplotype mode: 'one' (primary) or 'many' (all phased detours) [default: many]
        --hops <N>                Hop radius for subgraph expansion [default: 5]
        --max-detour-hops <N>     Maximum hops for detour reconvergence [default: 8]
        --min-abund <F>           Minimum ka:f for a variant [default: 5.0]
        --min-vaf <F>             Minimum VAF for an alternative branch [default: 0.10]
        --min-abund-frac <F>      Raise the variant floor to this fraction of gene depth,
                                  median ka:f on the walked path [default: 0.0 = off]
        --min-cov <F>             Minimum reference coverage [default: 0.85]
        --min-shared-kmers <N>    Anchor threshold, validated 1..=5 [default: 3]
        --max-size-mb <N>         Maximum compressed SRA size in MB [default: 350]
        --prune-frac <F>          Drop unitigs below this fraction of the gene's depth and
                                  realign, removing the aligner's reference bias. 0.0
                                  disables it [default: 0.01]
        --alt-max-detour <N>      Unitigs an alternative allele may traverse before it must
                                  rejoin the main route [default: 100]
    -t, --threads <N>             Threads per GraphAligner invocation [default: 1; above 1
                                  is not reproducible, see Config]
    -w, --workers <N>             Compute workers
        --fetch-workers <N>       Fetch workers [default: 8]

ALTERNATIVE HAPLOTYPES:
        --alt-engine <E>          Which implementation emits them: variants (local, records
                                  to alt_haplotypes.fasta), detour (upstream, records into
                                  the main FASTA), or off [default: variants]
    -H, --haplotypes <one|many>   Upstream detour mode; only read when --alt-engine=detour

SWITCHING OFF (applies to --alt-engine=variants):
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

    -h, --help                    Print this help
    -V, --version                 Print version
";

struct Args {
    accession: Option<String>,
    accession_list: Option<PathBuf>,
    query: Option<PathBuf>,
    out_dir: PathBuf,
    out_fasta: PathBuf,
    cfg: Config,
}

fn parse_args(argv: &[String]) -> error::Result<Args> {
    let mut args = Args {
        accession: None,
        accession_list: None,
        query: None,
        out_dir: PathBuf::from("walker_output"),
        out_fasta: PathBuf::from("reconstructed_haplotypes.fasta"),
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
            "-q" | "--query" => args.query = Some(PathBuf::from(value()?)),
            "-o" | "--out-dir" => args.out_dir = PathBuf::from(value()?),
            "-f" | "--out-fasta" => args.out_fasta = PathBuf::from(value()?),
            "-k" | "--k" => {
                args.cfg.k = value()?.parse().map_err(|e| error::Error::Usage(format!("--k: {e}")))?;
            }
            "--alt-engine" => {
                let raw = value()?;
                args.cfg.alt_engine = match raw.as_str() {
                    "variants" => types::AltEngine::Variants,
                    "detour" => types::AltEngine::Detour,
                    "off" | "none" => types::AltEngine::Off,
                    other => {
                        return Err(error::Error::Usage(format!(
                            "--alt-engine must be variants, detour or off, got {other:?}"
                        )))
                    }
                };
            }
            "-H" | "--haplotypes" => {
                let v = value()?;
                match v.to_lowercase().as_str() {
                    "one" | "single" => args.cfg.haplotypes = types::HaplotypeMode::One,
                    "many" | "all" | "multi" => args.cfg.haplotypes = types::HaplotypeMode::Many,
                    other => return Err(error::Error::Usage(format!("--haplotypes must be 'one' or 'many', got '{other}'"))),
                }
            }
            "--hops" => {
                args.cfg.hops = value()?.parse().map_err(|e| error::Error::Usage(format!("--hops: {e}")))?;
            }
            "--max-detour-hops" => {
                args.cfg.max_detour_hops = value()?.parse().map_err(|e| error::Error::Usage(format!("--max-detour-hops: {e}")))?;
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
                args.cfg.max_compressed_bytes = mb * 1024 * 1024;
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
    // ---- INTO the pipeline ------------------------------------------------
    let query_path = args
        .query
        .ok_or_else(|| error::Error::Usage("--query is required".into()))?;
    let query = io::input::read_query(&query_path)?;

    let accessions = match (&args.accession, &args.accession_list) {
        (Some(a), _) => vec![a.clone()],
        (None, Some(p)) => io::input::read_accession_list(p)?,
        (None, None) => {
            return Err(error::Error::Usage(
                "provide --accession or --accession-list".into(),
            ))
        }
    };

    // ---- the pipeline -----------------------------------------------------
    let summary = pipeline::run(&args.cfg, accessions, &query)?;

    // ---- OUT of the pipeline ----------------------------------------------
    std::fs::create_dir_all(&args.out_dir)?;
    io::output::write_fasta(&summary, &args.out_fasta)?;
    io::output::write_bubbles(&summary, &args.out_dir.join("bubbles.json"))?;
    if args.cfg.emit_vcf {
        io::output::write_vcfs(&summary, &query, &args.out_dir)?;
    }
    if args.cfg.emit_snv_table {
        io::output::write_snv_table(&summary, &query, &args.out_dir.join("snvs.tsv"))?;
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
    Ok(())
}
