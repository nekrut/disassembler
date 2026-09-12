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
//! Everything upstream (kmindex discovery) and downstream (MAFFT, FastTree, HyphAeon)
//! stays in Python. The file names below are the contract with `build_walker_codon_msa.py`
//! and `run_panel_analysis.py`, so they must not drift.

pub mod input;
pub mod output;
