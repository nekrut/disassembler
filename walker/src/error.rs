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
    /// Every accession failed. A per-accession failure is not fatal on its own, but a run
    /// that reconstructs nothing has produced no usable output and must not look like a
    /// success to a scheduler. `--allow-empty` opts out.
    NothingReconstructed { total: usize },
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
            Error::NothingReconstructed { total } => write!(
                f,
                "0/{total} accessions reconstructed; no output worth keeping. \
                 Check the log above for per-accession errors, or pass --allow-empty \
                 if an empty result is acceptable here"
            ),
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
