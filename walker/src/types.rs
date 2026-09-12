//! Data that flows between stages.
//!
//! Read this file first. Every stage boundary in `pipeline.rs` is one of these types,
//! so the shapes here are the pipeline's contract.

use std::path::PathBuf;

/// An SRA run accession, e.g. `ERR2028820`.
pub type Accession = String;

/// Unitig identifier. BCALM2 ids are dense and contiguous (`ERR2028820`: 1,277,362 unitigs,
/// max id 1,277,361), so 21 bits covers the range and these index arrays directly.
/// Validate density at parse time; fall back to a map if it ever fails.
pub type UnitigId = u32;

/// Strand of one end of a bidirected edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Strand {
    Forward,
    Reverse,
}

/// Traversal direction of a node within a GAF path (`>` forward, `<` reverse-complement).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orient {
    Forward,
    Reverse,
}

/// The reference coding sequence being searched for.
pub struct Query {
    /// File stem, e.g. `query_hiv1_pol`. Used for output filenames, because
    /// `logan_path_walker.py:410` names the summary from `query_fa.stem`. Keeping this
    /// identical is what lets the two implementations be diffed.
    pub name: String,
    /// First token of the FASTA header, e.g. `pfcrt`. Provenance and logging only.
    pub record_id: String,
    /// Uppercased sequence with all whitespace removed.
    pub seq: Vec<u8>,
    /// Path on disk. GraphAligner is handed this directly (`-f`), so it must stay valid
    /// for the life of the run.
    pub path: PathBuf,
}

// ---------------------------------------------------------------------------
// stage 1 output: stream
// ---------------------------------------------------------------------------

/// One accession's unitig file, still compressed.
///
/// Only the compressed bytes cross the channel. `stages::parse` decodes them as a *stream*
/// and never materialises the decompressed text, which for `ERR2028820` would be 145 MB and
/// for the largest accessions in the ERG11 panel about 1.44 GB.
///
/// Everything the rest of the pipeline needs (adjacency, abundance, sequences) is extracted
/// during that single pass, so this buffer can be dropped as soon as `parse::build` returns.
pub struct FetchedUnitigs {
    pub accession: Accession,
    pub compressed: Vec<u8>,
    pub compressed_len: usize,
    /// Seconds spent downloading, recorded by the fetch worker so the compute worker can
    /// report it. Network time belongs to the accession, not to whichever thread fetched it.
    pub fetch_secs: f64,
}

// ---------------------------------------------------------------------------
// stage 2 output: parse
// ---------------------------------------------------------------------------

/// Everything one accession's unitig file contains, in compact form.
///
/// Built in a single streaming pass so the decompressed text is never held. Measured
/// against `ERR2028820` (1,277,362 unitigs, 2,582,206 edge tokens, 68,745,081 bp):
///
/// | field | bytes |
/// |---|---|
/// | `offsets` | 5.1 MB |
/// | `targets` | 10.3 MB |
/// | `abundance` | 5.1 MB |
/// | `seq_2bit` | 17.2 MB |
/// | `seq_offset` | 5.1 MB |
/// | total | **~43 MB** |
///
/// Against ~897 MB for the Python dict-of-sets, and against 145 MB for the decompressed
/// text alone. Holding sequences here rather than re-reading the file is what lets
/// `extract::collect_seqs` be O(nodes kept) instead of O(all unitigs).
pub struct UnitigGraph {
    pub accession: Accession,
    pub n_unitigs: usize,
    /// `offsets[i]..offsets[i+1]` indexes `targets` for unitig `i`. Length `n_unitigs + 1`.
    pub offsets: Vec<u32>,
    /// Packed edge: target id in the low 30 bits, source strand in bit 30, target strand
    /// in bit 31. Use `stages::parse::{pack_edge, unpack_edge}` rather than doing it by hand.
    /// 30 bits caps an accession at 2^30 unitigs, which `parse::build` enforces.
    pub targets: Vec<u32>,
    /// Mean k-mer coverage (`ka:f:`) per unitig.
    pub abundance: Vec<f32>,
    /// All unitig sequences, 2 bits per base, 4 bases per byte.
    ///
    /// 2-bit packing is safe: Logan unitigs are uppercase ACGT only, because BCALM2 never
    /// emits a k-mer containing N. Verified 2026-09-10 by scanning 10 accessions spanning
    /// bacterial, fungal, plant-pathogen, avian-flu and five viral datasets:
    /// **850,835,755 bp across 18,379,016 unitigs, zero non-ACGT characters and no
    /// soft-masked (lowercase) bases.**
    ///
    /// That is strong evidence, not a guarantee. Still validate while parsing and fail
    /// loudly rather than silently corrupting a base; a stray N would otherwise be encoded
    /// as some other nucleotide and propagate into the haplotype.
    pub seq_2bit: Vec<u8>,
    /// Base offset of unitig `i` within `seq_2bit`, in bases (not bytes). Length
    /// `n_unitigs + 1`, so the length of unitig `i` is `seq_offset[i+1] - seq_offset[i]`.
    ///
    /// Packing is contiguous across unitig boundaries, so a sequence may begin mid-byte.
    /// That saves the ~1.9 MB of padding byte-alignment would cost on a 1.28M-unitig
    /// accession, at the price of a shift in `decode_seq`.
    pub seq_offset: Vec<u32>,
}

impl UnitigGraph {
    /// Length of unitig `id` in bases.
    pub fn seq_len(&self, id: UnitigId) -> usize {
        let i = id as usize;
        (self.seq_offset[i + 1] - self.seq_offset[i]) as usize
    }

    /// Decode unitig `id` into `out`, which is cleared first. Reuse one buffer across calls
    /// rather than allocating per unitig.
    pub fn decode_seq(&self, id: UnitigId, out: &mut Vec<u8>) {
        const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let i = id as usize;
        let (start, end) = (self.seq_offset[i] as usize, self.seq_offset[i + 1] as usize);
        out.clear();
        out.reserve(end - start);
        for p in start..end {
            out.push(BASES[((self.seq_2bit[p >> 2] >> (2 * (p & 3))) & 3) as usize]);
        }
    }

    /// 2-bit code of base `j` of unitig `id`, without materialising the sequence.
    /// This is the accessor the k-mer scan should use.
    #[inline]
    pub fn base_at(&self, id: UnitigId, j: usize) -> u64 {
        let p = self.seq_offset[id as usize] as usize + j;
        ((self.seq_2bit[p >> 2] >> (2 * (p & 3))) & 3) as u64
    }

    /// Neighbours of `id` as packed edges.
    #[inline]
    pub fn edges(&self, id: UnitigId) -> &[u32] {
        let i = id as usize;
        &self.targets[self.offsets[i] as usize..self.offsets[i + 1] as usize]
    }
}

// ---------------------------------------------------------------------------
// stage 3 output: anchor
// ---------------------------------------------------------------------------

/// Unitigs judged to overlap the query.
pub struct Anchors {
    pub ids: Vec<UnitigId>,
    /// Shared canonical k-mer count per anchor, parallel to `ids`.
    pub shared_kmers: Vec<u32>,
    /// `ids.len() / n_unitigs`, for diagnostics only.
    ///
    /// This is deliberately NOT the guard signal. Anchor fraction separates a healthy
    /// accession from a pathological one by only about 1.2x at the default threshold, which
    /// is far too little to act on. `extract::subgraph_fraction` separates the same cases by
    /// roughly 315x and is what `pipeline::run_one` actually checks.
    pub fraction_of_graph: f64,
}

// ---------------------------------------------------------------------------
// stage 4 output: extract
// ---------------------------------------------------------------------------

/// The 5-hop neighbourhood around the anchors, serialised as GFA.
pub struct Subgraph {
    pub accession: Accession,
    pub nodes: Vec<UnitigId>,
    /// Sequence per node, parallel to `nodes`.
    pub seqs: Vec<Vec<u8>>,
    /// Where the GFA was written. Should live on `/dev/shm`; it is ~38 KB.
    pub gfa_path: PathBuf,
    pub n_edges: usize,
}

/// Remove the scratch GFA when the subgraph goes out of scope.
///
/// This is the equivalent of the `finally` block at `logan_path_walker.py:327-330`, and it
/// is more reliable: it runs on every exit path, including early returns and errors.
///
/// It matters because `scratch_dir` is tmpfs, so a leaked file occupies RAM that never
/// appears in any process's RSS. At ~1 MB for a large subgraph, a 10,000-accession run
/// would quietly consume gigabytes.
impl Drop for Subgraph {
    fn drop(&mut self) {
        // Set LOGAN_KEEP_SCRATCH=1 to retain the GFA for debugging a specific accession.
        if std::env::var_os("LOGAN_KEEP_SCRATCH").is_none() {
            let _ = std::fs::remove_file(&self.gfa_path);
        }
    }
}

// ---------------------------------------------------------------------------
// stage 5 output: align (GraphAligner)
// ---------------------------------------------------------------------------

/// One step of a GAF path, e.g. `>860013` or `<27101`.
#[derive(Debug, Clone, Copy)]
pub struct PathStep {
    pub orient: Orient,
    pub node: UnitigId,
}

/// The single best alignment GraphAligner reported, as GAF columns.
///
/// Note `path_start` and `path_end` are offsets into the *stitched path*, not into the
/// query. `stitch` slices the assembled sequence with them (`logan_path_walker.py:217-220`),
/// so mistaking them for query coordinates yields a plausible but wrong haplotype.
pub struct GafPath {
    pub steps: Vec<PathStep>,
    /// GAF column 7.
    pub path_len: usize,
    /// GAF column 8.
    pub path_start: usize,
    /// GAF column 9.
    pub path_end: usize,
    /// GAF column 10, residue matches. Also the field the best row is chosen by.
    pub num_matches: usize,
    /// GAF column 11, alignment block length. `pident = num_matches / aln_len`.
    pub aln_len: usize,
}

// ---------------------------------------------------------------------------
// stage 6 output: stitch & detour
// ---------------------------------------------------------------------------

/// Which implementation produces the alternative haplotypes.
///
/// Two independent ones landed in this tree: `stages::detour` came from upstream and
/// `stages::variants` was written here, against the same defect. They differ in three ways
/// that matter, so only one runs at a time and both are kept intact until they can be
/// compared properly on the same data.
///
/// ```text
///                         detour (upstream)        variants (local)
/// gate at a fork          alternative branch only  BOTH branches must be credible
/// length-changing detour  allowed                  rejected
/// records go to           the main FASTA           alt_haplotypes.fasta
/// picks the major allele  ranks candidates by      --prefer-majority rewrites the
///                         mean ka:f                route per site
/// ```
///
/// The gate is the consequential one. On `SRR8209082` the alternative-only rule admits ~57
/// sites where the both-credible rule admits 1, because at 54 of them the branch the route
/// took holds under 10% of local coverage: those are the reference-bias defect, not two
/// alleles co-existing. Under `Detour` they reach the main FASTA, and
/// `build_walker_codon_msa.py` turns every record there into an independent taxon.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum AltEngine {
    /// Neither. One haplotype per accession, as before either feature existed.
    Off,
    /// Upstream's closed-loop detour splicing, `stages::detour`.
    Detour,
    /// The local implementation, `stages::variants`.
    Variants,
}

/// Haplotype emission mode for upstream's `stages::detour`, read only when
/// `AltEngine::Detour` is selected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HaplotypeMode {
    /// Emit only the single highest-coverage (primary/major) haplotype.
    One,
    /// Emit all phased branch detours (quasispecies/alleles) clearing abundance & VAF gates.
    Many,
}

/// A reconstructed full-length coding haplotype.
#[derive(Debug, Clone)]
pub struct Haplotype {
    pub id: String,
    pub seq: Vec<u8>,
    /// Fraction of the reference covered.
    pub cov: f64,
    pub pident: f64,
    pub nodes: usize,
    /// Mean k-mer coverage (`ka:f`) across the unitig path.
    pub mean_ka: f64,
    /// Variant allele frequency at the defining fork (1.0 if primary backbone).
    pub vaf: f64,
}

/// A recovered second allele: the primary haplotype with exactly one substitution.
///
/// The frequencies are what make it interpretable. `alt_vaf` is this allele's share of
/// coverage at its fork, `taken_vaf` the primary's, and both cleared the same gates, so the
/// record asserts "this sample carries both, in roughly this proportion" and nothing more.
pub struct AltHaplotype {
    pub haplotype: Haplotype,
    /// The fork on the walked path where the two alleles diverge.
    pub fork: UnitigId,
    pub fork_strand: Strand,
    /// First unitig of the branch the primary did not take.
    pub alt_node: UnitigId,
    pub alt_strand: Strand,
    pub taken_abundance: f32,
    pub taken_vaf: f64,
    pub alt_abundance: f32,
    pub alt_vaf: f64,
    /// Unitigs traversed before rejoining the walked path.
    pub detour_nodes: usize,
    /// Bases this allele adds (positive) or removes (negative) relative to the primary.
    /// Zero for a substitution, which is the common case.
    pub length_delta: i64,
    /// How much sequence the allele actually replaces, in bases.
    ///
    /// Reported because `max_alt_indel` bounds the *net* length change and not this, and the
    /// two come apart: across the panel one allele nets +21 nt while replacing 194. Small
    /// spans are single-site alleles, and a larger span carrying several linked substitutions
    /// is a phased haplotype block, which is worth having. Measured spans: substitution
    /// alleles median 1 nt (209 of 267 are 1 to 3 nt), indel alleles median 9 nt, max 194.
    /// No cap is applied, because one would discard the phased blocks along with the
    /// reroutings; the number is emitted so a reader can tell them apart.
    pub replaced_span: usize,
}

// ---------------------------------------------------------------------------
// stage 7 output: phase
// ---------------------------------------------------------------------------

/// One branch of a fork, with its coverage and variant allele frequency.
///
/// A branch is an *oriented* unitig. The same unitig entered from its two ends is two
/// different branches, which is exactly the distinction the Python implementation loses.
pub struct Branch {
    pub node: UnitigId,
    pub strand: Strand,
    pub abundance: f32,
    pub vaf: f64,
}

/// A phased bubble: one oriented unitig whose exits diverge and then meet again.
///
/// Gated on `ka:f >= min_abund` and `vaf >= min_vaf`.
///
/// `merge` is meaningful here, unlike in the Python version. There the adjacency drops the
/// strand fields, so a bidirected edge recorded in both endpoints' headers looks like a
/// reconvergence and the fork is reported as its own merge node: `fork == merge` on 93.1%
/// of `ERR2028820`'s bubbles and 58.5% of `SRR8209082`'s. Orientation is what separates
/// "this unitig has a neighbour at each end", which is an ordinary chain, from "this end of
/// this unitig has two alternatives", which is a real variant.
pub struct Bubble {
    pub fork: UnitigId,
    pub fork_strand: Strand,
    pub merge: UnitigId,
    pub merge_strand: Strand,
    pub branches: Vec<Branch>,
}

// ---------------------------------------------------------------------------
// per-accession result
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    Success,
    NoAnchors,
    NoPath,
    Failed,
}

/// Wall time per stage, in seconds. Mirrors the Python `timing` block so the two
/// implementations can be compared directly.
#[derive(Debug, Default, Clone, Copy)]
pub struct Timing {
    pub stream: f64,
    pub parse: f64,
    pub anchor: f64,
    pub extract: f64,
    pub align: f64,
    pub stitch: f64,
    pub phase: f64,
    pub total: f64,
}

/// Everything one accession produced. This is what crosses back to the collector.
pub struct AccessionResult {
    pub accession: Accession,
    pub outcome: Outcome,
    pub n_unitigs: usize,
    pub n_anchors: usize,
    pub subgraph_nodes: usize,
    pub haplotype: Option<Haplotype>,
    /// Everything written to the main FASTA. Under `AltEngine::Detour` this is the major
    /// allele followed by its phased detours; otherwise it holds the primary alone.
    pub haplotypes: Vec<Haplotype>,
    /// Second alleles recovered at biallelic sites, written to their own file. Populated
    /// only under `AltEngine::Variants`. Never combined with each other: one substitution
    /// per record.
    pub alts: Vec<AltHaplotype>,
    pub bubbles: Vec<Bubble>,
    pub timing: Timing,
    pub error: Option<String>,
}

impl AccessionResult {
    pub fn failed(accession: Accession, error: String) -> Self {
        Self {
            accession,
            outcome: Outcome::Failed,
            n_unitigs: 0,
            n_anchors: 0,
            subgraph_nodes: 0,
            haplotype: None,
            haplotypes: Vec::new(),
            alts: Vec::new(),
            bubbles: Vec::new(),
            timing: Timing::default(),
            error: Some(error),
        }
    }
}
