//! One file per pipeline step, in execution order.
//!
//! ```text
//!   stream   S3 + zstd            -> FetchedUnitigs
//!   parse    BCALM2 headers       -> UnitigGraph      (CSR)
//!   anchor   canonical k-mers     -> Anchors
//!   extract  5-hop BFS + GFA      -> Subgraph
//!   align    GraphAligner (C++)   -> GafPath
//!   stitch   path -> sequence     -> Haplotype
//!   phase    forks + VAF          -> Vec<Bubble>
//!   variants on-path forks       -> Vec<AltHaplotype>  (opt-in)
//! ```
//!
//! `phase` reads the `UnitigGraph` and `Subgraph`, not the `GafPath`, so it does not
//! depend on GraphAligner and could run concurrently with `align` if that ever mattered.

pub mod align;
pub mod anchor;
pub mod extract;
pub mod parse;
pub mod phase;
pub mod stitch;
pub mod stream;
pub mod variants;
