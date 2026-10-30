//! The boundary of the Rust pipeline.
//!
//! ```text
//!  INTO the pipeline (io::input)          OUT of the pipeline (io::output)
//!  ────────────────────────────           ────────────────────────────────
//!  --accession-list  FILE                 --out-fasta  reconstructed.fasta
//!  --accession       ACC                  bubbles.json
//!  --query           reference CDS FASTA  <query>_walker_summary.json
//! ```
//!
//! Everything upstream (kmindex discovery) and downstream (MAFFT, FastTree) stays in the
//! `disassembler` Python package. The file names below are the contract with the commands
//! there that run the walker and read its output, so they must not drift.

pub mod input;
pub mod output;
