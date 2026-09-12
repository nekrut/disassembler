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
