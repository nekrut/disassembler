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
mod term;
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
                                  reference whose alignment gives the best cov/pident).
                                  References must be strains of ONE gene: one primary is kept
                                  per accession, so a pool of different genes reports only the
                                  winner. The walker warns when the pool splits into groups

OUTPUT:
    -o, --out-dir <DIR>           Working directory
    -f, --out-fasta <FILE>        Reconstructed haplotypes FASTA

TUNING:
    -k, --k <N>                   Anchor k-mer size, 1..=32 [default: 21]
        --hops <N>                Hop radius for subgraph expansion [default: 5]
        --min-abund <F>           Minimum ka:f for a variant [default: 5.0]. The allele
                                  floor only: it does not move the prune floor, which is
                                  --prune-min-abund
        --min-vaf <F>             Minimum VAF for an alternative branch [default: 0.10]
        --min-abund-frac <F>      Raise the variant floor to this fraction of gene depth,
                                  median ka:f on the walked path [default: 0.0 = off]
        --min-cov <F>             Minimum fraction of the reference the alignment must span
                                  [default: 0.85]. The reference span, not the length of what
                                  was walked: an indel in the sample moves those apart, and
                                  the summary reports both as ref_cov and cov
        --primary-metric <M>      Multi-anchor primary selection: matched | cov | pident |
                                  cov_x_pident [default: matched, the fraction of the
                                  reference's bases matched]. Ignored in single-anchor mode
        --primary-tie-eps <F>     Treat references scoring within this of the best as tied and
                                  give the primary to the deepest route [default: 0.002]. On a
                                  mixed sample two strains can each match a panel record
                                  exactly, and without this the one listed first in the panel
                                  wins, so a 5% strain can be published as the primary. 0
                                  restores the pre-0.2.0 panel-order behaviour
        --min-shared-kmers <N>    Anchor threshold, validated 1..=5 [default: 3]
        --max-size-mb <N>         Maximum compressed SRA size in MB [default: 4096, which
                                  refuses 0.16% of Logan accessions; 1024 refused 2.26%]. Costs
                                  disk, not memory: see the disk bound printed at startup
        --spool-dir <DIR>         Where adjacency spills once it outgrows 4 MiB. Defaults to
                                  --unitigs-cache, else the system temp directory. Not
                                  --scratch-dir, which is /dev/shm by default: spilling into a
                                  tmpfs is not a saving. An accession whose spool will not fit
                                  the free space here is skipped as too_large, from its real
                                  size rather than from --max-size-mb; the startup lines print
                                  what the run can need
        --max-subgraph-nodes <N>  Refuse an accession whose anchor neighbourhood exceeds N
                                  unitigs, as subgraph_too_large [default: 1000000; 0 = off].
                                  Everything after extraction is proportional to the
                                  neighbourhood, so this and not --max-size-mb is what bounds
                                  memory once the gene is found
        --prune-frac <F>          Drop unitigs below this fraction of the gene's depth and
                                  realign, removing the aligner's reference bias. The
                                  realignment is kept only if it spans as much of the
                                  reference, give or take a terminal codon. 0.0 disables it
                                  [default: 0.01]
        --prune-min-abund <F>     Floor under --prune-frac: a unitig at or above this ka:f
                                  survives the prune whatever the fraction works out to
                                  [default: 5.0]. On a gene shallower than
                                  prune-min-abund / prune-frac this is the whole rule
        --alt-max-detour <N>      Unitigs an alternative allele may traverse before it must
                                  rejoin the main route [default: 100]
    -t, --threads <N>             Threads per GraphAligner invocation, at least 1 [default: 1;
                                  above 1 is not reproducible, see Config]
        --graphaligner <PATH>     GraphAligner executable [default: GraphAligner on PATH].
                                  Checked once at startup; results are validated against
                                  1.0.17b
        --graphaligner-timeout <SECS>
                                  Kill a single GraphAligner call after this many seconds and
                                  fail that accession only, instead of letting it hold its
                                  worker until the job's wall-time kill discards the whole
                                  shard [default: 3600; 0 waits forever]
    -w, --workers <N>             Compute workers, at least 1 [default: 8]
        --fetch-workers <N>       Fetch workers, at least 1 [default: 8]
        --fetch-attempts <N>      Failed attempts per accession before it is given up, at
                                  least 1 [default: 6]. 408, 429, 5xx and transport errors
                                  are retried with a jittered exponential backoff, honouring
                                  Retry-After; an attempt that transferred real data before
                                  the connection dropped resumes and does not count
        --fetch-timeout <SECS>    Budget for one fetch attempt, including connect and body
                                  [default: 300]. An interrupted transfer resumes, so this
                                  bounds how long a stalled connection goes unnoticed, not
                                  how long an accession may take

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
        --summary-seqs            Also write each haplotype's literal sequence into the
                                  summary, as 0.1.x did. Off by default: the summary carries
                                  seq_sha256 and length, and the sequence itself is already in
                                  reconstructed.fasta [default: off]
        --rescue-minority         Also move the route off a branch the aligner walked that is
                                  itself a minority of its fork, onto an outright majority.
                                  OFF by default: on a 5,000-accession walker20k slice the
                                  correction loop did not converge with it on (97,068
                                  corrections against 470, 7.8 per net substitution). The
                                  sites are reported either way [default: off]
        --max-alts <N>            Most alternative alleles to emit per accession. Credible
                                  alleles beyond it are counted and warned about, not dropped
                                  silently, and the omitted ones are the 3' sites
                                  [default: 32]
        --no-indel-alts           Substitution alleles only; refuse any that change length
        --allow-frameshift-alts   Also accept indel alleles whose length change is not a
                                  multiple of 3. These reframe the CDS; off by default
        --max-alt-indel <N>       Largest length change still treated as an indel allele
                                  rather than a rerouting [default: 30]

        --allow-empty             Exit 0 when no accession reconstructed because nothing was
                                  found: no Logan data, over --max-size-mb, no anchors or no
                                  path. Errors are never forgiven: one erroring accession
                                  fails the run (exit 4) even when others reconstructed, and
                                  they are listed in <query stem>.skipped.tsv. Without this
                                  flag, a run that reconstructs nothing exits non-zero too, so
                                  a wrapper cannot mistake an empty shard for a good one

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
        --revalidate              Confirm each cache hit against the bucket before using it.
                                  Logan's objects are mutable, so a warm cache can serve
                                  bytes the bucket has since replaced and nothing in the
                                  output says so. Sends a conditional request against the
                                  ETag recorded beside the cached object: unchanged costs a
                                  round trip and no transfer, changed refreshes the entry and
                                  says so. Requires --unitigs-cache; not with --offline
        --fetch-only              Fetch every accession into --unitigs-cache and exit before
                                  the compute stages. Requires --unitigs-cache; -q is not
                                  required in this mode. Writes fetch.summary.json into
                                  --out-dir recording each object's size, sha256, ETag and
                                  Last-Modified, so the cache is a citable frozen input

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

/// A finite float in `[lo, hi]`, or a usage error naming the flag.
///
/// WV-15. These were a bare `parse()`, which accepts `inf` and `NaN`. Both then reached
/// `summary.json` through `format!("{:?}")` as the bare tokens `inf` and `NaN`, **which no
/// JSON parser accepts**, so the run exited 0 and left a summary that `json.load` rejects.
/// `--min-cov NaN` also disables the coverage gate entirely, since every comparison against
/// NaN is false, and `--min-vaf 5` was accepted although a fraction above 1 can never be met.
fn finite_in(flag: &str, raw: &str, lo: f64, hi: f64) -> error::Result<f64> {
    let v: f64 = raw
        .parse()
        .map_err(|e| error::Error::Usage(format!("{flag}: {e}")))?;
    if !v.is_finite() {
        return Err(error::Error::Usage(format!(
            "{flag}: {raw:?} is not a finite number"
        )));
    }
    if v < lo || v > hi {
        let upper = if hi == f64::MAX { "".to_string() } else { format!(" and at most {hi}") };
        return Err(error::Error::Usage(format!(
            "{flag}: {v} is out of range (must be at least {lo}{upper})"
        )));
    }
    Ok(v)
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
            "--prune-min-abund" => {
                args.cfg.prune_min_abund =
                    finite_in("--prune-min-abund", &value()?, 0.0, f64::MAX)? as f32;
            }
            "--min-abund" => {
                args.cfg.min_abund = finite_in("--min-abund", &value()?, 0.0, f64::MAX)? as f32;
            }
            "--min-vaf" => {
                args.cfg.min_vaf = finite_in("--min-vaf", &value()?, 0.0, 1.0)?;
            }
            "--min-cov" => {
                args.cfg.min_cov = finite_in("--min-cov", &value()?, 0.0, 1.0)?;
            }
            "--primary-tie-eps" => {
                args.cfg.primary_tie_eps = finite_in("--primary-tie-eps", &value()?, 0.0, 1.0)?;
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
            "--spool-dir" => args.cfg.spool_dir = Some(PathBuf::from(value()?)),
            "--max-subgraph-nodes" => {
                args.cfg.max_subgraph_nodes = value()?
                    .parse()
                    .map_err(|e| error::Error::Usage(format!("--max-subgraph-nodes: {e}")))?;
            }
            "--max-size-mb" => {
                let mb: u64 = value()?.parse().map_err(|e| error::Error::Usage(format!("--max-size-mb: {e}")))?;
                // saturate rather than wrap: an absurd value should mean "no cap", not a tiny one.
                args.cfg.max_compressed_bytes = mb.saturating_mul(1024 * 1024);
            }
            "-t" | "--threads" => args.cfg.graphaligner_threads = at_least_one("--threads", &value()?)?,
            "-w" | "--workers" => args.cfg.compute_workers = at_least_one("--workers", &value()?)?,
            "--fetch-workers" => args.cfg.fetch_workers = at_least_one("--fetch-workers", &value()?)?,
            "--fetch-attempts" => {
                args.cfg.fetch_attempts =
                    at_least_one("--fetch-attempts", &value()?)?.try_into().map_err(|_| {
                        error::Error::Usage("--fetch-attempts: too large".to_string())
                    })?;
            }
            "--fetch-timeout" => {
                args.cfg.fetch_timeout_secs =
                    at_least_one("--fetch-timeout", &value()?)? as u64;
            }
            "--scratch-dir" => args.scratch_base = Some(PathBuf::from(value()?)),
            "--graphaligner" => args.cfg.graphaligner_bin = value()?,
            // 0 is meaningful (wait forever), so this is a plain parse, not `at_least_one`.
            "--graphaligner-timeout" => {
                args.cfg.graphaligner_timeout_secs = value()?.parse().map_err(|e| {
                    error::Error::Usage(format!("--graphaligner-timeout: {e}"))
                })?;
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
            // TW-04 found this was the only boolean without a positive spelling, while the
            // other six all have one. It is a no-op against the default, but a user scripting
            // the flags should not have to know which direction each one is written in.
            "--indel-alts" | "--no-indel-alts" => {
                args.cfg.allow_indel_alts = flag == "--indel-alts";
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
            "--revalidate" => {
                args.cfg.revalidate = true;
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
            "--max-alts" => {
                let raw = value()?;
                args.cfg.max_alts = raw.parse().map_err(|_| {
                    error::Error::Usage(format!("--max-alts needs a number, got {raw:?}"))
                })?;
                i += 2;
                continue;
            }
            "--rescue-minority" | "--no-rescue-minority" => {
                args.cfg.rescue_minority = flag == "--rescue-minority";
                i += 1;
                continue;
            }
            "--summary-seqs" | "--no-summary-seqs" => {
                args.cfg.summary_seqs = flag == "--summary-seqs";
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

/// Stop glibc's dynamic mmap threshold from turning freed graph arrays into resident slack.
///
/// The parser's five CSR arrays grow by doubling (`stages::parse::build`). glibc services a
/// large allocation with `mmap` only while the request is above its mmap threshold, and that
/// threshold is *dynamic*: each time an mmapped block is freed, glibc raises it towards 32 MiB
/// on the theory that the next allocation of that size should come from the heap instead. One
/// accession per run therefore behaves nothing like twenty. From the second accession on, the
/// arrays grow inside the heap arenas, freed chunks are not returned to the OS, and the steady
/// state measured about 3.17x the compressed size rather than the 1.7x
/// `GRAPH_PEAK_PER_COMPRESSED_BYTE` was calibrated at, which was measured one accession per
/// run.
///
/// Pinning the threshold at 1 MiB makes every array above it an `mmap`/`munmap` pair, so the
/// memory goes back to the OS when the graph is dropped. Measured over 20 cached 34 MB
/// accessions at `-w 1 --fetch-workers 8`: peak RSS 386.5 -> 340.3 MiB, no output change.
///
/// glibc only. musl has no dynamic threshold and no `mallopt` tunable worth setting here, and
/// on other platforms this is a no-op.
///
/// **Deliberately not unit-tested.** glibc exposes no getter for the threshold, and the only
/// observable effect is peak RSS over a multi-accession run, which is not a thing to assert in
/// a unit test. Removing this call changes no output and fails nothing; the evidence is the
/// measurement above, and this comment is here so that a future reader finds the gap stated
/// rather than assumes the guard is covered.
#[cfg(all(target_os = "linux", target_env = "gnu"))]
fn pin_mmap_threshold() {
    // Safe: `mallopt` takes two ints and touches no memory we own. A failure is not actionable
    // (the allocator keeps its default behaviour), so the return value is deliberately ignored.
    unsafe {
        libc::mallopt(libc::M_MMAP_THRESHOLD, 1 << 20);
    }
}

#[cfg(not(all(target_os = "linux", target_env = "gnu")))]
fn pin_mmap_threshold() {}

fn main() -> ExitCode {
    #[cfg(unix)]
    default_sigpipe();
    pin_mmap_threshold();
    let argv: Vec<String> = std::env::args().skip(1).collect();

    if argv.is_empty() || argv.iter().any(|a| a == "-h" || a == "--help") {
        print!("{USAGE}");
        return ExitCode::SUCCESS;
    }
    if argv.iter().any(|a| a == "-V" || a == "--version") {
        crate::outln!("logan-walker {}", env!("CARGO_PKG_VERSION"));
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
/// A line each for the prune outcomes and the provisional results, which `summary.json`
/// records per accession, so a large run's log says how often each happened.
fn decision_lines(results: &[types::AccessionResult], guard: f64) -> Vec<String> {
    use types::PruneResult as P;
    let mut lines = Vec::new();
    let parts: Vec<String> = [
        (P::Applied, "applied"),
        (P::RefusedShorterSpan, "refused (realignment spans less of the reference)"),
        (P::RefusedNoHaplotype, "refused (realignment gave no haplotype)"),
        (P::NothingToRemove, "nothing to remove"),
        (P::WouldEmpty, "skipped (would have removed every unitig)"),
    ]
    .into_iter()
    .filter_map(|(want, label)| {
        let n = results.iter().filter(|r| r.prune.as_ref().map(|p| p.result) == Some(want)).count();
        (n > 0).then(|| format!("{n} {label}"))
    })
    .collect();
    if !parts.is_empty() {
        lines.push(format!("[logan-walker] prune: {}", parts.join(", ")));
    }
    let provisional = results.iter().filter(|r| r.provisional).count();
    if provisional > 0 {
        lines.push(format!(
            "[logan-walker] {provisional} of {} results provisional: neighbourhood above {:.1} \
             unitigs per reference base",
            results.len(),
            guard
        ));
    }
    lines
}

/// WP-02. Say so when the `-q` pool is not one gene.
///
/// Multi-anchor mode ranks every reference against every other and keeps one primary per
/// accession. That is what it is for across strains of a gene; across *different* genes it
/// silently answers a different question than the user asked, because an accession carrying
/// both genes comes back with whichever scored higher, exit 0, and the other one appears only
/// as an unselected candidate in `summary.json`. 0.1.0 printed "One gene per query" for a
/// multi-record `-q`; 0.1.1 dropped the warning and put nothing in its place, while the README
/// goes on listing the 25-gene master target file among the references to start from.
///
/// A warning and not a refusal: a deliberately divergent panel also splits, and the walker
/// cannot know which the user meant. Returns the line to print, or `None` when the pool is one
/// group, so a test can read it without running the binary.
fn mixed_gene_warning(cfg: &config::Config, queries: &[types::Query]) -> Option<String> {
    let seqs: Vec<&[u8]> = queries.iter().map(|q| q.seq.as_slice()).collect();
    let groups = stages::anchor::reference_groups(cfg.k, &seqs, stages::anchor::GROUP_MIN_SHARE);
    if groups.len() < 2 {
        return None;
    }
    // Named by the first record of each group, which is enough to recognise the gene, with
    // the rest of the group as a count: a 25-gene panel must stay one readable line.
    const SHOWN: usize = 6;
    let names: Vec<String> = groups
        .iter()
        .take(SHOWN)
        .map(|g| match g.len() {
            1 => format!("[{}]", queries[g[0]].record_id),
            n => format!("[{} +{}]", queries[g[0]].record_id, n - 1),
        })
        .collect();
    let more = if groups.len() > SHOWN {
        format!(" and {} more", groups.len() - SHOWN)
    } else {
        String::new()
    };
    Some(format!(
        "[warn] multi-anchor: the {} references fall into {} groups sharing under {:.0}% of \
         their {}-mers ({}{}); one primary is kept per accession, so an accession carrying \
         more than one of these genes reports only the winner and the rest stay unselected \
         candidates in summary.json",
        queries.len(),
        groups.len(),
        stages::anchor::GROUP_MIN_SHARE * 100.0,
        cfg.k,
        names.join(" "),
        more
    ))
}

fn check_run(results: &[types::AccessionResult], allow_empty: bool) -> error::Result<()> {
    let total = results.len();
    let success = results.iter().filter(|r| r.outcome == types::Outcome::Success).count();
    let errors = results.iter().filter(|r| r.outcome == types::Outcome::Failed).count();
    if total == 0 {
        return Ok(());
    }
    if errors == total {
        return Err(error::Error::AllErrored { total });
    }
    // Any error is a failure of the run, even beside successes. This used to return `Ok` the
    // moment one accession succeeded, so a 5,000-accession run in which 8 panicked exited 0
    // and the scheduler that launched it saw a green job. The successes are still written;
    // what changes is that the exit code stops claiming they were the whole story.
    if errors > 0 {
        return Err(if success > 0 {
            error::Error::SomeErrored { total, success, errors }
        } else {
            error::Error::NothingReconstructed { total, errors }
        });
    }
    if success > 0 || allow_empty {
        return Ok(());
    }
    Err(error::Error::NothingReconstructed { total, errors })
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
    if args.cfg.revalidate && args.cfg.unitigs_cache.is_none() {
        return Err(error::Error::Usage(
            "--revalidate requires --unitigs-cache <DIR>: there is nothing to revalidate \
             without a cache"
                .into(),
        ));
    }
    // Refused rather than silently ignored. Both flags are about trusting the cache and they
    // want opposite things: --offline forbids the network this needs, so a run with both
    // would quietly do no revalidation at all while its settings claimed otherwise.
    if args.cfg.revalidate && args.cfg.offline {
        return Err(error::Error::Usage(
            "--revalidate and --offline are mutually exclusive: revalidating means asking \
             the bucket, which --offline forbids"
                .into(),
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
        // Try the legacy single-file path first for its n_records semantics. If that file
        // turns out to hold more than one record, the whole pool goes into the anchor union
        // and the align fan-out instead.
        //
        // Both decisions are made from ONE read of the file. Every input the walker takes may
        // be a pipe (`-q <(zcat ref.fa.gz)`, `-q /dev/stdin`), which can only be read once,
        // and this used to call `read_query` and then `read_query_set`, opening the same path
        // twice: a multi-record query through a pipe died with "no FASTA records found" on
        // the second open, naming a file that was perfectly good (WA-01).
        Some(if args.query.len() == 1 {
            let path = args.query[0].clone();
            let bytes = io::input::read_input(&path, "query")?;
            let first = io::input::read_query_from(&path, &bytes)?;
            if first.n_records > 1 {
                io::input::read_query_set_from(&[(path, bytes)])?
            } else {
                vec![first]
            }
        } else {
            io::input::read_query_set(&args.query)?
        })
    };

    // `n_listed` is what the user asked for, before duplicates were dropped; the summary
    // reports the difference rather than letting it disappear (`docs/summaries.md`).
    let (accessions, n_listed) = match (&args.accession, &args.accession_list) {
        // `-a` used to win silently, so a run meant for a whole list did one accession.
        (Some(_), Some(_)) => {
            return Err(error::Error::Usage(
                "give --accession or --accession-list, not both".into(),
            ))
        }
        (Some(a), None) => {
            // The list path has always validated; `-a` did not, though the value reaches the
            // same URL and the same two filesystem paths.
            let acc = io::input::validate_accession(a, "--accession")?;
            // WF-12: `-a srr8209082` used to 404 and be recorded as `no_logan_data`, a
            // biological negative about a run Logan does have.
            if acc != *a {
                eprintln!(
                    "[warn] --accession {a:?} read as {acc:?}: Logan's bucket is                      case-sensitive, and the output files use the uppercase form"
                );
            }
            (vec![acc], 1)
        }
        (None, Some(p)) => io::input::read_accession_list_counted(p)?,
        (None, None) => {
            return Err(error::Error::Usage(
                "provide --accession or --accession-list".into(),
            ))
        }
    };

    // ---- the pipeline -----------------------------------------------------
    // Adjacency spills next to the cache when there is one, because that is a directory the
    // user has already chosen and sized for this run's data. Otherwise the system temp
    // directory. Never the scratch directory, which is a tmpfs by default.
    if args.cfg.spool_dir.is_none() {
        args.cfg.spool_dir = Some(
            args.cfg
                .unitigs_cache
                .clone()
                .unwrap_or_else(std::env::temp_dir),
        );
    }
    if let Some(d) = &args.cfg.spool_dir {
        if d.starts_with("/dev/shm") {
            eprintln!(
                "[warn] --spool-dir {} is a tmpfs, so spilled adjacency stays in RAM and the \
                 memory bound below does not hold",
                d.display()
            );
        }
    }

    // State the memory bound up front. This used to be a comment in `pipeline`, and it was
    // wrong by 5x; a number that moves with -w, --fetch-workers and --max-size-mb is worth
    // printing rather than documenting.
    let gb = |b: u64| b as f64 / (1024.0 * 1024.0 * 1024.0);
    // What the pipeline will actually process: `n_listed` counts the list before duplicates
    // were dropped, and a duplicate never becomes a second in-flight buffer.
    let n = accessions.len();
    let buffers = args.cfg.peak_compressed_buffers_for(n);
    let graphs = args.cfg.peak_graphs_for(n);
    // Where the slots come from. Under --fetch-only there is no queue and no compute pool, so
    // saying "queued + fetching + parsing" there would name two things that cannot happen.
    let shape = if args.cfg.fetch_only {
        format!("{} fetching", args.cfg.fetch_workers)
    } else {
        format!(
            "{} queued + {} fetching + {} parsing",
            args.cfg.channel_capacity, args.cfg.fetch_workers, args.cfg.compute_workers
        )
    };
    // A short list cannot fill the pipeline, and quoting the pipeline's ceiling at a
    // one-accession run overstates it by the worker count. Compared against this mode's own
    // ceiling: under --fetch-only there are no graphs at all, which is the mode, not the list.
    let plural = |k: usize| if k == 1 { "" } else { "s" };
    let capped = if buffers < args.cfg.peak_compressed_buffers()
        || graphs < args.cfg.peak_graphs_for(usize::MAX)
    {
        format!(", capped by {n} accession{}", plural(n))
    } else {
        String::new()
    };
    let graph_clause = if args.cfg.fetch_only {
        // Saying "0 parsed graphs at up to 1.7x = 0.0 GB" invites the reader to wonder what
        // would make it non-zero. Nothing would: --fetch-only never builds one.
        ", and no parsed graphs, because --fetch-only does not parse".to_string()
    } else if args.cfg.graph_bound_is_set_by_the_gene() {
        // The graph term stopped scaling with the file when adjacency began spooling, but
        // this line went on quoting --max-size-mb x 1.7 and leaving the next line to correct
        // it by a factor of fifty. It now states the bound that applies and that the parser
        // enforces; the next line says what makes it true.
        format!(
            ", plus {graphs} parsed graph{} bounded by --max-subgraph-nodes {} at {:.2} GB \
             each = {:.2} GB",
            plural(graphs),
            args.cfg.max_subgraph_nodes,
            gb(args.cfg.max_graph_bytes_per_accession()),
            gb(args.cfg.max_graph_resident_bytes_for(n)),
        )
    } else {
        format!(
            ", plus {graphs} parsed graph{} at up to {:.1}x their compressed size = {:.1} GB",
            plural(graphs),
            config::GRAPH_PEAK_PER_COMPRESSED_BYTE,
            gb(args.cfg.max_graph_resident_bytes_for(n)),
        )
    };
    // The compressed term used to be the whole story and is now, usually, nothing. A
    // download is written to disk as it arrives, so an in-flight accession costs one copy
    // buffer rather than its own size, and `--max-size-mb` stops appearing in the memory
    // bound at all. Stated in the units it actually has: "0.00 GB" invites the reader to
    // wonder what is being hidden, where "2.0 MiB" answers it.
    let mib = |b: u64| b as f64 / (1024.0 * 1024.0);
    let buffer_clause = if args.cfg.downloads_stream_to_disk() {
        format!(
            "{buffers} in-flight accession{} ({shape}{capped}) x a {} KiB stream buffer = \
             {:.1} MiB, because compressed objects go to disk as they arrive",
            plural(buffers),
            config::DOWNLOAD_BUFFER_BYTES / 1024,
            mib(args.cfg.max_compressed_resident_bytes_for(n)),
        )
    } else {
        format!(
            "{buffers} compressed buffer{} ({shape}{capped}) x {} MB cap = {:.1} GB",
            plural(buffers),
            args.cfg.max_compressed_bytes / (1024 * 1024),
            gb(args.cfg.max_compressed_resident_bytes_for(n)),
        )
    };
    // Two decimals, because the number this prints is routinely under 1 GB: at the defaults
    // it is 1.24, and one decimal made the total read 0.1 next to a graph term of 0.16,
    // which looks like an arithmetic error rather than a rounding one.
    eprintln!(
        "[logan-walker] memory bound: {buffer_clause}{graph_clause}; {:.2} GB worst case. \
         Excludes GraphAligner, the scratch GFA and the output buffers.",
        gb(args.cfg.max_resident_bytes_for(n)),
    );
    // This used to restate the graph term with the correct number, because the line above
    // had the wrong one. The line above is now right, so this says what makes it right, and
    // what happens when the thing that makes it right is unavailable.
    if args.cfg.graph_bound_is_set_by_the_gene() && !args.cfg.fetch_only {
        eprintln!(
            "[logan-walker] that bound does not mention --max-size-mb, because nothing in the \
             walk scales with accession size: adjacency spools to {} and the neighbourhood is \
             capped. An accession whose adjacency cannot spill there is refused as too_large \
             rather than parsed in RAM.",
            args.cfg.spool_dir.as_ref().map(|d| d.display().to_string()).unwrap_or_default(),
        );
    }

    // Disk, which neither line above covers and which spooling made load-bearing: the memory
    // win of stage 3 was paid for in bytes that now land on a filesystem nobody was asked to
    // size. Two consumers, named separately rather than summed, because they do not share a
    // shape. The spool is transient, dies with its accession and scales with `-w`; the cache
    // is persistent, outlives the run and scales with the accession count. At a 28 GB cap the
    // spool term is 15 GB at `-w 1` and 120 GB at `-w 8`, while the cache term for 50
    // accessions is 1.4 TB at either. Summing them would suggest that lowering `-w` helps
    // with the larger of the two, and it does not.
    let spool_bytes = args.cfg.max_spool_bytes_for(n);
    let cache_bytes = args.cfg.max_cache_bytes_for(n);
    // The bytes streaming moved off the heap. With a cache they land in it and are counted
    // there; without one they are unlinked files under --spool-dir, one per in-flight
    // accession, which is exactly the term that used to be in the memory line above.
    let download_bytes = args.cfg.max_download_bytes_for(n);
    if spool_bytes > 0 || download_bytes > 0 || cache_bytes.is_some() {
        let mut terms: Vec<String> = Vec::new();
        // Distinct directories, in the order they are named. The spool defaults to the cache,
        // so the common case has one filesystem and should report its free space once.
        let mut dirs: Vec<&std::path::Path> = Vec::new();
        if let (true, Some(d)) = (download_bytes > 0, &args.cfg.spool_dir) {
            terms.push(format!(
                "{buffers} downloaded object{} x {} MB cap = {:.2} GB transient in {}",
                plural(buffers),
                args.cfg.max_compressed_bytes / (1024 * 1024),
                gb(download_bytes),
                d.display(),
            ));
            dirs.push(d.as_path());
        }
        if let (true, Some(d)) = (spool_bytes > 0, &args.cfg.spool_dir) {
            terms.push(format!(
                "{graphs} spool{} x {:.2} GB = {:.2} GB transient in {}",
                plural(graphs),
                gb(args.cfg.max_spool_bytes_per_accession()),
                gb(spool_bytes),
                d.display(),
            ));
            if !dirs.contains(&d.as_path()) {
                dirs.push(d.as_path());
            }
        }
        if let (Some(b), Some(d)) = (cache_bytes, &args.cfg.unitigs_cache) {
            terms.push(format!(
                "{n} accession{} x {} MB = {:.2} GB persistent in {}",
                plural(n),
                args.cfg.max_compressed_bytes / (1024 * 1024),
                gb(b),
                d.display(),
            ));
            if !dirs.contains(&d.as_path()) {
                dirs.push(d.as_path());
            }
        }
        // A filesystem that will not answer is left out rather than reported as zero, the
        // same way the fetch gate disables itself on one.
        let free: Vec<String> = dirs
            .iter()
            .filter_map(|d| {
                types::available_bytes(d).map(|a| format!("{:.1} GB free in {}", gb(a), d.display()))
            })
            .collect();
        let free_clause =
            if free.is_empty() { String::new() } else { format!("; {} now", free.join(", ")) };
        // Only name the axis of a term that is actually printed. `--fetch-only` has no spool
        // and an `--offline` run has no cache, and telling either of them what `-w` does to a
        // number it is not being shown is how a reader learns to stop reading these lines.
        let axes = match (spool_bytes > 0 || download_bytes > 0, cache_bytes.is_some()) {
            (true, true) => {
                " The transient terms move with the worker counts, the cache term with the \
                 list length."
            }
            (true, false) => " The transient terms move with the worker counts.",
            _ => " The cache term moves with the list length, not with the worker counts.",
        };
        eprintln!(
            "[logan-walker] disk bound: {}{free_clause}. Worst case, every accession at the \
             cap.{axes}",
            terms.join("; "),
        );
    }

    // The run's own clock, for `run.started` and `run.wall_seconds` in the summary. Started
    // here rather than at process entry so it measures the work, not the argument parsing,
    // and `n_requested` is taken before the pipeline consumes the list. Above the
    // `--fetch-only` branch because that path writes a summary too.
    let started = std::time::SystemTime::now();
    let run_clock = std::time::Instant::now();
    let n_requested = n_listed;

    // --fetch-only takes the short path: fill the cache and exit. No compute pool, no
    // reconstruction, no query. Success = every accession is in the cache after this.
    //
    // It does write `fetch.summary.json`, because a cache with no record of where its bytes
    // came from is not a frozen input, it is just a directory. Rule 4 of the summary contract
    // applies here as everywhere: the summary is written even when the fetch fails.
    if args.cfg.fetch_only {
        // WP-06, same reason as the walk below: an `-o` that cannot be made is a fact about
        // this command, not about the bucket, so find it before downloading a cohort. No
        // journal here, because the cache entries a fetch writes ARE its durable output: a
        // killed `--fetch-only` leaves everything it fetched, and rerunning skips it.
        std::fs::create_dir_all(&args.out_dir).map_err(error::file_err("create", &args.out_dir))?;
        let summary = pipeline::run_fetch_only(&args.cfg, accessions)?;
        let ok = summary.success;
        let negatives = summary.results.iter().filter(|r| r.outcome.is_negative()).count();
        let errors = summary.results.len() - ok - negatives;
        eprintln!(
            "[logan-walker --fetch-only] {ok} fetched, {negatives} not available \
             (no Logan data or over --max-size-mb), {errors} failed"
        );

        let summary_path = args.out_dir.join("fetch.summary.json");
        // Decided here so the summary can state the code rather than guess it from counts.
        let verdict = check_fetch(ok, negatives, errors, args.allow_empty);
        let run = io::output::RunInfo {
            // GraphAligner never runs here, so `tool.external` stays empty rather than
            // naming a tool that did no work.
            graphaligner: None,
            cfg: &args.cfg,
            queries: &[],
            started,
            wall_seconds: run_clock.elapsed().as_secs_f64(),
            n_requested,
            exit_code: verdict.as_ref().err().map_or(0, |e| e.exit_code()),
        };
        let written = std::fs::create_dir_all(&args.out_dir)
            .map_err(error::file_err("create", &args.out_dir))
            .and_then(|()| {
                io::output::write_summary(&summary, &run, &summary_path)
                    .map_err(|e| e.at("write", &summary_path))
            });
        if let Err(e) = written {
            // A failed write must not mask a failed fetch, and must not be swallowed when
            // the fetch was fine.
            eprintln!("[logan-walker] write failed: {e}");
            verdict?;
            return Err(e);
        }
        eprintln!("[logan-walker --fetch-only] wrote {}", summary_path.display());
        return verdict;
    }

    // Unwrap the query set for the full pipeline. The guard above ensures it is Some here.
    let queries = queries.expect("queries is Some when not --fetch-only");

    // WV-14: the id reaches CHROM and `##contig=<ID=...>` verbatim, so refuse one that
    // cannot legally appear there rather than writing a VCF no parser will read.
    for q in &queries {
        let id = if q.record_id.is_empty() { &q.name } else { &q.record_id };
        if !io::output::valid_contig_id(id) {
            return Err(error::Error::Usage(format!(
                "reference id {id:?} cannot be a VCF contig name. Rename the record in {}: \
                 letters, digits and !#$%&+./:;?@^_|~- are allowed, and * = after the first \
                 character.",
                q.path.display()
            )));
        }
    }

    // Checked before any download. GraphAligner used to be looked up per accession, after
    // its unitigs were fetched and parsed, so a missing binary failed each one in turn.
    let bin = stages::align::resolve_binary(&args.cfg)?;
    let version = stages::align::version(&bin)?;
    eprintln!("[logan-walker] GraphAligner: {version} ({})", bin.display());
    if !version.contains("1.0.17b") {
        eprintln!(
            "[warn] the walker is validated against GraphAligner 1.0.17b; 1.0.19 and later \
             clip the query's last base (docs/notes/graphaligner_end_clipping.md)"
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
        if let Some(line) = mixed_gene_warning(&args.cfg, &queries) {
            eprintln!("{line}");
        }
    }

    // WP-06, the cheap half: touch the output locations BEFORE the cohort runs. An
    // unusable `-o` (a regular file of that name) or an `-f` whose parent does not exist
    // used to be found only after every accession had been fetched and aligned. The probe
    // leaves the filesystem as it found it: `-f` is created only to see whether it can be,
    // and removed again if it was not already there, so the canonical FASTA stays something
    // only a finished run writes.
    std::fs::create_dir_all(&args.out_dir).map_err(error::file_err("create", &args.out_dir))?;
    if let Some(parent) = args.out_fasta.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(error::file_err("create", parent))?;
        }
    }
    let fasta_existed = args.out_fasta.exists();
    std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&args.out_fasta)
        .map_err(error::file_err("open", &args.out_fasta))?;
    if !fasta_existed {
        let _ = std::fs::remove_file(&args.out_fasta);
    }

    // WP-06, the rest: a journal of what has finished, so a run that is killed leaves its
    // reconstructions and a list of what it got through. Removed below once the real outputs
    // are on disk.
    let journal = io::output::Journal::create(&args.out_dir, &queries[0].name)?;
    let (partial_fa, journal_jsonl) = journal.paths();
    // stderr, with the other startup notices. Not stdout: that stream carries the progress
    // lines a pipeline reads, and a run whose reader goes away ends at its first stdout
    // write, so putting a startup line there would move when that happens.
    eprintln!(
        "[logan-walker] progress is journalled to {} and {}; both are removed when the run \
         finishes writing its outputs",
        partial_fa.display(),
        journal_jsonl.display()
    );

    let mut summary = pipeline::run(&args.cfg, accessions, &queries, Some(&journal))?;

    // ---- OUT of the pipeline ----------------------------------------------
    //
    // WV-18: every writer is attempted and the summary is ALWAYS written, with the first
    // error returned at the end. These used to be chained with `?` and the summary came
    // last, so one write failure after the compute (a full disk, a directory in the way)
    // threw away the whole shard's per-accession status and timings although every accession
    // had reconstructed. The run still exits non-zero, so the failure is not hidden.
    let mut first_err: Option<error::Error> = None;
    // The run's verdict, reached before any output is written so the summary can state the
    // exit code instead of inferring it from counts, and returned unchanged at the end.
    let verdict = check_run(&summary.results, args.allow_empty);
    // Scoped so the closure's mutable borrow of `first_err` ends before it is read.
    {
    let mut keep = |r: error::Result<()>| {
        if let Err(e) = r {
            eprintln!("[logan-walker] write failed: {e}");
            first_err.get_or_insert(e);
        }
    };
    keep(io::output::write_fasta(&summary, &args.out_fasta).map_err(|e| e.at("write", &args.out_fasta)));
    let bubbles = args.out_dir.join("bubbles.json");
    keep(io::output::write_bubbles(&summary, &bubbles).map_err(|e| e.at("write", &bubbles)));
    // One alignment pass per accession, shared by both writers.
    if args.cfg.emit_vcf || args.cfg.emit_snv_table {
        let sites = io::output::snv_sites_all(&mut summary, &queries);
        if args.cfg.emit_vcf {
            keep(io::output::write_vcfs(&summary, &queries, &args.out_dir, &sites)
                .map_err(|e| e.at("write VCFs in", &args.out_dir)));
        }
        if args.cfg.emit_snv_table {
            let snvs = args.out_dir.join("snvs.tsv");
            keep(io::output::write_snv_table(&summary, &queries, &snvs, &sites)
                .map_err(|e| e.at("write", &snvs)));
        }
    }
    if args.cfg.emit_alts {
        let alts = args.out_dir.join("alt_haplotypes.fasta");
        keep(io::output::write_alts(&summary, &alts).map_err(|e| e.at("write", &alts)));
    }
    // WV-17: a rerun with a writer switched off used to leave the previous run's file beside
    // fresh results, with nothing to show the mismatch. Remove our own leftovers; a file we
    // did not write is left alone (WV-16).
    io::output::remove_disabled_outputs(&args.cfg, &args.out_dir, &summary);
    let summary_path = args.out_dir.join(format!("{}.summary.json", queries[0].name));
    let bin_display = bin.display().to_string();
    let run = io::output::RunInfo {
        graphaligner: Some((&version, &bin_display)),
        cfg: &args.cfg,
        queries: &queries,
        started,
        wall_seconds: run_clock.elapsed().as_secs_f64(),
        n_requested,
        exit_code: verdict.as_ref().err().map_or(0, |e| e.exit_code()),
    };
    keep(io::output::write_summary(&summary, &run, &summary_path)
        .map_err(|e| e.at("write", &summary_path)));
    // Beside the summary, and on the same terms: written whatever the run's outcome, so a
    // workflow can read it without first deciding whether the run went well enough to have
    // produced one.
    let skipped_path = args.out_dir.join(format!("{}.skipped.tsv", queries[0].name));
    keep(io::output::write_skipped(&summary, &skipped_path)
        .map_err(|e| e.at("write", &skipped_path)));
    }
    if let Some(e) = first_err {
        // The real outputs are incomplete, so the journal is the only complete record of
        // what this run produced. It stays.
        return Err(e);
    }
    // Everything the run produces is on disk, so the salvage copies have nothing left to
    // say. WP-06: they exist only while a run is in flight, or after one died.
    journal.discard();

    crate::outln!(
        "[logan-walker] {}/{} accessions reconstructed",
        summary.success, summary.total
    );
    for line in decision_lines(&summary.results, args.cfg.subgraph_nodes_per_base_guard) {
        crate::outln!("{line}");
    }

    // Per-accession failures are collected rather than propagated, which is what lets one bad
    // accession not sink a 10,000-accession shard. The cost is that a run where *every*
    // accession failed used to write a 0-byte FASTA and exit 0, which any scheduler reads as
    // a good shard; until 0.1.1, --allow-empty still did that for a run where every
    // accession errored. Outputs are already written above, so they stay on disk for
    // inspection.
    verdict
}

#[cfg(test)]
mod flag_tests {
    use super::{error, parse_args, Args, USAGE};

    fn args(extra: &[&str]) -> error::Result<Args> {
        let mut v: Vec<String> = ["-a", "SRR1", "-q", "q.fa"].iter().map(|s| s.to_string()).collect();
        v.extend(extra.iter().map(|s| s.to_string()));
        parse_args(&v)
    }

    /// TW-04. Every long flag the usage text advertises must actually be accepted.
    ///
    /// `parse_args` is a hand-written match over about 40 spellings and nothing checked that
    /// the match and the help agree. A mutant parser that swapped `--min-vaf` with
    /// `--min-cov`, inverted `--prefer-majority`, ignored `--emit-vcf` and skipped unknown
    /// flags passed the entire suite. This test is self-maintaining: a flag added to USAGE
    /// but not to the parser fails here without anyone remembering to extend a list.
    #[test]
    fn every_flag_in_the_usage_text_is_accepted() {
        // Flags that need a value get a plausible one; the rest are bare.
        let value_for = |flag: &str| -> Option<&'static str> {
            match flag {
                "--min-cov" | "--min-vaf" | "--prune-frac" | "--min-abund-frac"
                | "--primary-tie-eps" => Some("0.5"),
                "--min-abund" | "--prune-min-abund" => Some("5"),
                "--max-alts" | "--max-alt-indel" | "--alt-max-detour" | "--hops"
                | "--min-shared-kmers" | "--max-size-mb" | "-w" | "--workers"
                | "--fetch-workers" | "-t" | "--threads" | "-k" | "--graph-k"
                | "--fetch-attempts" | "--fetch-timeout" => Some("2"),
                "--subgraph-nodes-per-base-guard" => Some("3"),
                "--graphaligner-timeout" => Some("60"),
                "--primary-metric" => Some("matched"),
                "-a" | "--accession" => Some("SRR1"),
                "-l" | "--accession-list" | "-q" | "--query" | "-o" | "--out-dir"
                | "-f" | "--out-fasta" | "--scratch-dir" | "--unitigs-cache"
                | "--graphaligner" | "--walker-bin" => Some("x"),
                _ => None,
            }
        };
        let mut checked = 0;
        for token in USAGE.split_whitespace() {
            let flag = token.trim_matches(|c: char| !c.is_ascii_alphanumeric() && c != '-');
            // `--help` and `--version` are answered in `main` before `parse_args` ever runs,
            // so the parser does not know them and should not.
            if !flag.starts_with("--")
                || flag.len() < 4
                || flag.contains('<')
                || matches!(flag, "--help" | "--version")
            {
                continue;
            }
            let mut argv = vec![flag];
            if let Some(v) = value_for(flag) {
                argv.push(v);
            }
            match args(&argv) {
                Err(error::Error::Usage(m)) if m.starts_with("unknown flag") => {
                    panic!("{flag} is in the usage text but the parser rejects it: {m}")
                }
                _ => checked += 1, // accepted, or refused for a reason that is not "unknown"
            }
        }
        assert!(checked > 25, "expected to exercise most of the flags, saw {checked}");
    }

    /// The other half: a flag the parser does not know must be refused, not ignored. The
    /// mutant accepted `--emit-vfc` and exited 0 having silently written no VCF.
    #[test]
    fn an_unknown_flag_is_refused() {
        let err = match args(&["--emit-vfc"]) {
            Err(e) => e,
            Ok(_) => panic!("a misspelled flag must not be accepted"),
        };
        assert!(
            matches!(&err, error::Error::Usage(m) if m.contains("unknown flag --emit-vfc")),
            "a misspelled flag must not be ignored, got {err:?}"
        );
    }

    /// TW-04. Each flag sets the field it names, and the boolean pairs move in opposite
    /// directions. Covers the 0.2.0 additions, which are the ones with no history behind
    /// them.
    #[test]
    fn each_flag_sets_the_field_it_names() {
        let base = args(&[]).unwrap();
        assert_eq!(base.cfg.min_cov, 0.85);
        assert_eq!(base.cfg.min_vaf, 0.10);

        // The pair the mutant swapped.
        let a = args(&["--min-cov", "0.3", "--min-vaf", "0.5"]).unwrap();
        assert_eq!((a.cfg.min_cov, a.cfg.min_vaf), (0.3, 0.5));

        // --max-size-mb is megabytes in, bytes out; the mutant scaled it twice.
        let a = args(&["--max-size-mb", "1"]).unwrap();
        assert_eq!(a.cfg.max_compressed_bytes, 1024 * 1024);

        // Booleans, both directions, including every flag 0.2.0 added.
        for (flag, get) in [
            ("prefer-majority", (|c: &crate::config::Config| c.prefer_majority) as fn(&_) -> bool),
            ("emit-vcf", |c| c.emit_vcf),
            ("emit-snv-table", |c| c.emit_snv_table),
            ("emit-alts", |c| c.emit_alts),
            ("indel-alts", |c| c.allow_indel_alts),
            ("summary-seqs", |c| c.summary_seqs),
            ("rescue-minority", |c| c.rescue_minority),
        ] {
            let on = args(&[&format!("--{flag}")]).unwrap();
            let off = args(&[&format!("--no-{flag}")]).unwrap();
            assert!(get(&on.cfg), "--{flag} must turn it on");
            assert!(!get(&off.cfg), "--no-{flag} must turn it off");
        }

        let a = args(&["--max-alts", "7"]).unwrap();
        assert_eq!(a.cfg.max_alts, 7);
        let a = args(&["--allow-frameshift-alts"]).unwrap();
        assert!(a.cfg.allow_frameshift_alts);
    }
}

#[cfg(test)]
mod range_tests {
    use super::finite_in;

    /// WV-15. These flags were a bare `parse()`, so `inf` and `NaN` were accepted and then
    /// written into `summary.json` as bare tokens that **no JSON parser accepts**: the run
    /// exited 0 and left a summary `json.load` rejects. `--min-cov NaN` additionally disables
    /// the coverage gate, since every comparison against NaN is false.
    #[test]
    fn non_finite_and_out_of_range_values_are_refused() {
        for bad in ["inf", "-inf", "NaN", "nan", "infinity"] {
            assert!(finite_in("--min-cov", bad, 0.0, 1.0).is_err(), "{bad} must be refused");
        }
        // A fraction above 1 can never be met, and a negative one is meaningless.
        assert!(finite_in("--min-vaf", "5", 0.0, 1.0).is_err(), "5 is not a fraction");
        assert!(finite_in("--min-vaf", "-0.1", 0.0, 1.0).is_err());
        assert!(finite_in("--min-abund", "-1", 0.0, f64::MAX).is_err());
    }

    #[test]
    fn ordinary_values_still_parse() {
        assert_eq!(finite_in("--min-cov", "0.85", 0.0, 1.0).unwrap(), 0.85);
        assert_eq!(finite_in("--min-vaf", "0", 0.0, 1.0).unwrap(), 0.0);
        assert_eq!(finite_in("--min-vaf", "1", 0.0, 1.0).unwrap(), 1.0);
        assert_eq!(finite_in("--min-abund", "5", 0.0, f64::MAX).unwrap(), 5.0);
        assert_eq!(finite_in("--min-abund", "25000", 0.0, f64::MAX).unwrap(), 25000.0);
    }
}

#[cfg(test)]
mod tests {
    use super::{at_least_one, check_fetch, check_run, decision_lines, mixed_gene_warning, parse_args};
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

    /// WP-02. Two unrelated genes in the anchor pool must draw a warning naming them, and a
    /// panel of one gene's strains must not. A pool of different genes is not refused: the
    /// walker cannot know whether the user meant a catalogue or a panel, and only one of the
    /// two is a mistake.
    #[test]
    fn a_pool_of_different_genes_is_called_out_and_a_strain_panel_is_not() {
        let seq = |seed: u64, n: usize| -> Vec<u8> {
            let mut x = seed;
            (0..n)
                .map(|_| {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    b"ACGT"[((x >> 33) & 3) as usize]
                })
                .collect()
        };
        let query = |id: &str, s: Vec<u8>| crate::types::Query {
            name: "panel".into(),
            record_id: id.into(),
            seq: s,
            path: std::path::PathBuf::from("panel.fa"),
            n_records: 2,
        };
        let cfg = crate::config::Config::default();

        let gene_a = seq(11, 900);
        let gene_b = seq(22, 900);
        let line = mixed_gene_warning(&cfg, &[
            query("bla_KPC", gene_a.clone()),
            query("bla_NDM", gene_b),
        ])
        .expect("two unrelated genes must warn");
        assert!(line.starts_with("[warn] multi-anchor:"), "{line}");
        assert!(line.contains("2 groups"), "{line}");
        assert!(line.contains("[bla_KPC]") && line.contains("[bla_NDM]"), "{line}");
        assert!(line.contains("21-mers"), "the threshold is in anchor k-mers: {line}");

        // Two strains of one gene: 9 substitutions in 900 bases.
        let mut strain = gene_a.clone();
        for i in (0..strain.len()).step_by(100) {
            strain[i] = if strain[i] == b'A' { b'C' } else { b'A' };
        }
        assert!(
            mixed_gene_warning(&cfg, &[query("wuhan", gene_a), query("jn1", strain)]).is_none(),
            "a strain panel is exactly what multi-anchor mode is for"
        );
    }

    #[test]
    fn decision_lines_count_each_outcome_and_skip_the_absent_ones() {
        use crate::types::{PruneReport, PruneResult};
        let prune = |result| Some(PruneReport {
            result, depth: 100.0, floor: 5.0, unitigs: 10, unitigs_kept: 4,
            unpruned_span: (0, 30), pruned_span: None,
        });
        let mut r = run_of(&[Outcome::Success, Outcome::Success, Outcome::Success, Outcome::NoPath]);
        r[0].prune = prune(PruneResult::Applied);
        r[1].prune = prune(PruneResult::Applied);
        r[2].prune = prune(PruneResult::RefusedShorterSpan);
        r[1].provisional = true;
        assert_eq!(decision_lines(&r, 3.0), vec![
            "[logan-walker] prune: 2 applied, 1 refused (realignment spans less of the reference)",
            "[logan-walker] 1 of 4 results provisional: neighbourhood above 3.0 unitigs per \
reference base",
        ]);
        assert!(decision_lines(&run_of(&[Outcome::NoAnchors]), 0.01).is_empty());
    }

    /// One success no longer excuses an error beside it.
    ///
    /// This reverses the rule the walker shipped until 0.2.0, where `check_run` returned `Ok`
    /// the moment any accession succeeded. Under that rule a 5,000-accession run in which 8
    /// accessions panicked exited 0, and the only way to notice was to read `records[].status`
    /// out of the summary by hand.
    #[test]
    fn an_error_beside_a_success_is_still_a_failed_run() {
        let r = run_of(&[Outcome::Failed, Outcome::Success, Outcome::NoAnchors]);
        for allow_empty in [false, true] {
            let e = check_run(&r, allow_empty).unwrap_err();
            assert!(matches!(e, Error::SomeErrored { total: 3, success: 1, errors: 1 }), "{e}");
            // Neither 1 nor 3: this run reconstructed something, and a wrapper that wants to
            // keep the output it did get has to be able to tell it from a run that produced
            // none and from one where nothing worked at all.
            assert_eq!(e.exit_code(), 4);
        }
    }

    /// Accessions where nothing was found still never fail a run that reconstructed something.
    ///
    /// The strict rule is about errors only. Widening it to the negatives would fail almost
    /// every real cohort, since most accessions simply do not carry the gene.
    #[test]
    fn negatives_beside_a_success_are_a_clean_run() {
        let r = run_of(&[
            Outcome::Success,
            Outcome::NoAnchors,
            Outcome::NoLoganData,
            Outcome::TooLarge,
            Outcome::SubgraphTooLarge,
        ]);
        for allow_empty in [false, true] {
            assert!(check_run(&r, allow_empty).is_ok());
        }
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
