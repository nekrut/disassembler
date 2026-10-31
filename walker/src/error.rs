//! Error type for the walker.

use std::fmt;

#[derive(Debug)]
pub enum Error {
    /// Placeholder while the pipeline is a skeleton. Names the stage that was reached,
    /// so a run shows exactly how far the flow gets.
    NotImplemented(&'static str),
    Io(std::io::Error),
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
    /// Logan holds no unitigs for this accession (HTTP 404). A property of the data, like
    /// finding no anchors, not a failure of the run.
    NoLoganData { accession: String },
    /// Over `--max-size-mb`: skipped by the user's own size cap, not a failure.
    TooLarge { accession: String, detail: String },
    /// No accession was reconstructed. A per-accession failure is not fatal on its own, but
    /// a run that reconstructs nothing has produced no usable output and must not look like
    /// a success to a scheduler. `errors` of the `total` failed outright; the rest found
    /// nothing (no Logan data, over the size cap, no anchors, no path). `--allow-empty`
    /// accepts the second kind only.
    NothingReconstructed { total: usize, errors: usize },
    /// Every accession failed with an error: the run itself did not work. Exit code 3, so a
    /// wrapper can tell this from "nothing found", and `--allow-empty` never hides it.
    AllErrored { total: usize },
    /// `--fetch-only` left accessions out of the cache.
    FetchIncomplete { total: usize, errors: usize, negatives: usize },
}

impl Error {
    /// Process exit code: 3 when every accession errored, 1 for everything else.
    pub fn exit_code(&self) -> u8 {
        match self {
            Error::AllErrored { .. } => 3,
            _ => 1,
        }
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::NotImplemented(what) => write!(f, "not implemented: {what}"),
            Error::Io(e) => write!(f, "io error: {e}"),
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
                write!(f, "[{accession}] no Logan unitigs for this accession (HTTP 404)")
            }
            Error::TooLarge { accession, detail } => {
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

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}

pub type Result<T> = std::result::Result<T, Error>;
