//! Error type for the walker.

use crate::types::Payload;
use std::fmt;
use std::path::{Path, PathBuf};

#[derive(Debug)]
pub enum Error {
    /// Placeholder while the pipeline is a skeleton. Names the stage that was reached,
    /// so a run shows exactly how far the flow gets.
    NotImplemented(&'static str),
    Io(std::io::Error),
    /// A file operation failed, naming the file. A bare `Io` did not: a scratch directory
    /// another user had created failed every accession with an unexplained EACCES.
    File { op: &'static str, path: PathBuf, source: std::io::Error },
    /// Network or S3 failure while streaming unitigs.
    Fetch { accession: String, detail: String },
    /// zstd frame could not be decoded.
    Decompress { accession: String, detail: String },
    /// A unitig header did not match the expected BCALM2 shape.
    Parse { accession: String, detail: String },
    /// GraphAligner exited non-zero or was killed by a signal.
    GraphAligner { accession: String, detail: String },
    /// Bad CLI input.
    Usage(String),
    /// Logan holds no unitigs for this accession. A property of the data, like finding no
    /// anchors, not a failure of the run.
    ///
    /// Two routes reach it, and they are the same fact. The bucket returns **404**, or it
    /// serves an object that decompresses to **nothing**: Logan publishes a 55-byte zstd
    /// frame, valid and empty, for accessions it assembled nothing from, about 0.4% of the
    /// bucket. Treating the second as a parse error reports a cohort of failures it has not
    /// had.
    NoLoganData { accession: String },
    /// Over `--max-size-mb`: skipped by the user's own size cap, not a failure.
    ///
    /// Carries the provenance the refusal was decided on, because this is the one negative
    /// where the walker knows something about the object it did not read: when S3 sent a
    /// `Content-Length` the size is exact, and the response's `ETag` and `Last-Modified`
    /// identify the object that was too big. Without it the size survived only as prose
    /// inside `detail`, which a workflow cannot branch on.
    TooLarge { accession: String, detail: String, payload: Box<Payload> },
    /// No accession was reconstructed. A per-accession failure is not fatal on its own, but
    /// a run that reconstructs nothing has produced no usable output and must not look like
    /// a success to a scheduler. `errors` of the `total` failed outright; the rest found
    /// nothing (no Logan data, over the size cap, no anchors, no path). `--allow-empty`
    /// accepts the second kind only.
    NothingReconstructed { total: usize, errors: usize },
    /// Every accession failed with an error: the run itself did not work. Exit code 3, so a
    /// wrapper can tell this from "nothing found", and `--allow-empty` never hides it.
    AllErrored { total: usize },
    /// Some accessions errored while others reconstructed. **Any error makes the run
    /// non-zero**, so this is not success, but the run did produce usable output and must not
    /// be confused with a run that produced none.
    ///
    /// Its own exit code (4) for that reason: 1 means nothing was reconstructed and 3 means
    /// nothing worked at all, while here `<query>.skipped.tsv` and the FASTA are both real and
    /// a wrapper may want to keep them. Until 0.2.0 this case exited 0, so a 5,000-accession
    /// run in which 8 accessions panicked reported success to the scheduler that launched it.
    SomeErrored { total: usize, success: usize, errors: usize },
    /// `--fetch-only` left accessions out of the cache.
    FetchIncomplete { total: usize, errors: usize, negatives: usize },
}

impl Error {
    /// Say which file an I/O error was about: `Io` becomes `File`, anything else is kept.
    pub fn at(self, op: &'static str, path: &Path) -> Error {
        match self {
            Error::Io(source) => Error::File { op, path: path.to_path_buf(), source },
            e => e,
        }
    }

    /// Process exit code: 3 when every accession errored, 1 for everything else.
    pub fn exit_code(&self) -> u8 {
        match self {
            Error::AllErrored { .. } => 3,
            Error::SomeErrored { .. } => 4,
            _ => 1,
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::NotImplemented(what) => write!(f, "not implemented: {what}"),
            Error::Io(e) => write!(f, "io error: {e}"),
            Error::File { op, path, source } => write!(f, "cannot {op} {}: {source}", path.display()),
            Error::Fetch { accession, detail } => write!(f, "[{accession}] fetch failed: {detail}"),
            Error::Decompress { accession, detail } => {
                write!(f, "[{accession}] decompress failed: {detail}")
            }
            Error::Parse { accession, detail } => write!(f, "[{accession}] parse failed: {detail}"),
            Error::GraphAligner { accession, detail } => {
                write!(f, "[{accession}] GraphAligner failed: {detail}")
            }
            Error::Usage(m) => write!(f, "usage error: {m}"),
            Error::NoLoganData { accession } => {
                write!(f, "[{accession}] no Logan unitigs for this accession")
            }
            Error::TooLarge { accession, detail, .. } => {
                write!(f, "[{accession}] skipped, over --max-size-mb: {detail}")
            }
            Error::NothingReconstructed { total, errors: 0 } => write!(
                f,
                "0/{total} accessions reconstructed, and none errored: nothing was found \
                 (see the log above for why). Pass --allow-empty if an empty result is \
                 acceptable here"
            ),
            Error::NothingReconstructed { total, errors } => write!(
                f,
                "0/{total} accessions reconstructed, and {errors} of them failed with an \
                 error (see the log above). --allow-empty accepts accessions where nothing \
                 was found, not errors"
            ),
            Error::AllErrored { total } => write!(
                f,
                "all {total} accessions failed with an error (see the log above); the run \
                 did not work. --allow-empty does not apply to errors"
            ),
            Error::SomeErrored { total, success, errors } => write!(
                f,
                "{success}/{total} accessions reconstructed, but {errors} of them failed with \
                 an error (see the log above and the skipped table). The output that was \
                 produced is still written; an accession that errored is a failure of the run, \
                 not a thing that was not found, so this is not exit 0"
            ),
            Error::FetchIncomplete { total, errors, negatives } => {
                write!(f, "{} of {total} accessions were not cached", errors + negatives)?;
                if *errors > 0 {
                    write!(f, "; {errors} failed with an error (see the log above)")?;
                }
                if *negatives > 0 {
                    write!(
                        f,
                        "; {negatives} have no Logan data or are over --max-size-mb \
                         (--allow-empty accepts these)"
                    )?;
                }
                Ok(())
            }
        }
    }
}

impl std::error::Error for Error {}

/// `map_err` adapter for a raw `io::Error` about `path`.
pub fn file_err(op: &'static str, path: &Path) -> impl FnOnce(std::io::Error) -> Error {
    let path = path.to_path_buf();
    move |source| Error::File { op, path, source }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}

pub type Result<T> = std::result::Result<T, Error>;
