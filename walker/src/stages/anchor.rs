//! Step 3: find the unitigs that carry the query gene.
//!
//! Replaces `logan_path_walker.py:69-99` (`find_anchors_blast`), which writes a BLAST
//! database of every unitig to disk and runs megablast against it. That is the single
//! largest stage in the Python walker, 5.8 s of a 14.6 s accession.
//!
//! Validated 2026-09-10 across 10 accessions, 3 organisms, 3 genes, by running the full
//! downstream on each anchor set and comparing haplotype SHA-256:
//!
//! - byte-identical haplotypes on 5 of the 6 reconstructing accessions for
//!   `min_shared_kmers` in 1..=5; two accessions shift at 10
//! - anchor recall is a poor proxy for correctness. Output was identical at recall 0.375,
//!   because the 5-hop BFS bridges anchors the k-mer pass misses
//! - the low apparent overlap with BLAST is an artifact of BLAST's `-evalue 1e-10`, which
//!   rejects a bare 31 bp exact match. At `1e-5`, 91.5% of k-mer anchors are in BLAST's set
//! - requiring N *consecutive* matching k-mers is a measured no-op. A unitig is contiguous,
//!   so its matching k-mers form one block and the run length equals the total count
//!
//! IN   UnitigGraph, decompressed unitig bytes, Query
//! OUT  Anchors

use crate::config::Config;
use crate::error::{Error, Result};
use crate::types::{Anchors, Query, UnitigGraph, UnitigId};

/// 2-bit code for a nucleotide, or `None` for anything else.
///
/// Logan unitigs are pure ACGT (verified over 850 Mbp), but a reference CDS pulled from
/// NCBI can carry `N` or IUPAC ambiguity codes. Those must be skipped, never coerced:
/// encoding `N` as `A` would invent a k-mer that matches unitigs the query does not contain.
#[inline]
pub fn base_code(b: u8) -> Option<u64> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Canonical k-mer set for the reference CDS. Built once per run, not per accession.
///
/// # Why not `HashSet<u64>`
///
/// `anchor::find` probes this once per base of every unitig, roughly 68 million times for
/// a single mid-sized accession. `HashSet` with the default SipHash costs on the order of
/// 25 ns per probe, which is seconds of pure hashing per accession. The set is tiny
/// (~1,245 entries for PfCRT, ~10 KB) and never changes after construction, so an
/// open-addressed table with multiply-shift hashing sits entirely in L1 and probes in a
/// few nanoseconds.
pub struct QueryIndex {
    /// Open addressing, power-of-two length, linear probing. `EMPTY` marks a free slot.
    table: Vec<u64>,
    mask: u64,
    pub k: usize,
    /// Distinct canonical k-mers stored.
    pub n_kmers: usize,
    /// Positions skipped because the window contained a non-ACGT base.
    pub n_skipped: usize,
}

/// Sentinel for an empty slot.
///
/// Safe because a *canonical* k-mer can never be `u64::MAX`. For k < 32 the top bits are
/// always clear, so the point is trivial. At k = 32 the mask is all ones and an all-T
/// k-mer really does encode as `u64::MAX`, but its reverse complement is all-A, which is
/// 0, and `canonical = min(fwd, rev)` therefore returns 0. Only a value that is `MAX` in
/// both orientations could collide, and no such k-mer exists.
const EMPTY: u64 = u64::MAX;

/// Odd 64-bit constant for Fibonacci hashing.
const MULT: u64 = 0x9E37_79B9_7F4A_7C15;

impl QueryIndex {
    #[inline]
    pub fn contains(&self, kmer: u64) -> bool {
        let mut i = ((kmer.wrapping_mul(MULT)) >> 32) & self.mask;
        loop {
            let slot = self.table[i as usize];
            if slot == kmer {
                return true;
            }
            if slot == EMPTY {
                return false;
            }
            i = (i + 1) & self.mask;
        }
    }

    fn insert(&mut self, kmer: u64) -> bool {
        let mut i = ((kmer.wrapping_mul(MULT)) >> 32) & self.mask;
        loop {
            let slot = self.table[i as usize];
            if slot == kmer {
                return false; // already present
            }
            if slot == EMPTY {
                self.table[i as usize] = kmer;
                return true;
            }
            i = (i + 1) & self.mask;
        }
    }

    pub fn len(&self) -> usize {
        self.n_kmers
    }

    pub fn is_empty(&self) -> bool {
        self.n_kmers == 0
    }
}

/// Build the query k-mer index from a single reference.
///
/// Thin wrapper around [`index_query_set`]. Every current call site holds a single
/// [`Query`]; multi-anchor mode calls the set variant directly with the whole pool.
/// The single-query path stays byte-identical: with one query, the set variant sums
/// windows over `[query]`, which is `query.seq.len() - k + 1` — the same capacity
/// computation the old direct implementation used — and iterates the same bytes in the
/// same order, so the resulting table is the same table.
///
/// Canonical means `min(kmer, revcomp(kmer))`, which makes the later scan strand-agnostic
/// without scanning each unitig twice. Both registers are rolled in lockstep:
///
/// ```text
/// fwd = ((fwd << 2) | code) & mask
/// rev = (rev >> 2) | ((3 - code) << 2(k-1))
/// ```
///
/// A window containing any non-ACGT base is skipped, and the run length resets, so no
/// k-mer spanning an ambiguity code is ever emitted.
pub fn index_query(cfg: &Config, query: &Query) -> Result<QueryIndex> {
    index_query_set(cfg, std::slice::from_ref(query))
}

/// Build one canonical-k-mer index unioning every record in `queries`.
///
/// Multi-anchor mode's anchor pool: pass Wuhan + JN.1 (or a whole strain panel) and every
/// canonical k-mer from any reference is a valid seed. The open-addressed table dedupes
/// insertions, so overlapping k-mers cost table space only once and the `contains` probe
/// stays a single hash.
///
/// Capacity comes from the sum of window counts across queries (rounded up to a power of
/// two, load factor ≤ 0.5). For a single-query call this reduces to the same capacity the
/// legacy single-query path used, and the k-mer stream is iterated in the same order, so
/// [`index_query`] is a byte-identical wrapper.
///
/// Warnings and the "no usable k-mers" error are aggregate: a run whose whole pool yielded
/// nothing is still a hard error, and skipped windows are reported per query so a bad
/// record is identifiable in a wide pool.
pub fn index_query_set(cfg: &Config, queries: &[Query]) -> Result<QueryIndex> {
    let k = cfg.k;
    if k == 0 || k > 32 {
        return Err(Error::Usage(format!(
            "k must be in 1..=32 to fit a u64 (got {k})"
        )));
    }
    if queries.is_empty() {
        return Err(Error::Usage("no queries in anchor pool".into()));
    }
    for q in queries {
        if q.seq.len() < k {
            return Err(Error::Usage(format!(
                "query {:?} is {} bp, shorter than k={k}",
                q.name,
                q.seq.len()
            )));
        }
    }

    let total_windows: usize = queries.iter().map(|q| q.seq.len() - k + 1).sum();
    // Load factor <= 0.5, rounded up to a power of two, so linear probing stays short.
    let mut cap = 16usize;
    while cap < total_windows * 2 {
        cap <<= 1;
    }
    let mut idx = QueryIndex {
        table: vec![EMPTY; cap],
        mask: (cap - 1) as u64,
        k,
        n_kmers: 0,
        n_skipped: 0,
    };

    let shift = 2 * (k as u64 - 1);
    let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };

    for query in queries {
        let (mut fwd, mut rev, mut run) = (0u64, 0u64, 0usize);
        let mut skipped_this_query = 0usize;

        for (pos, &b) in query.seq.iter().enumerate() {
            match base_code(b) {
                Some(code) => {
                    fwd = ((fwd << 2) | code) & mask;
                    rev = (rev >> 2) | ((3 - code) << shift);
                    run += 1;
                }
                None => {
                    // Ambiguity code: drop the window and start counting again after it.
                    run = 0;
                    fwd = 0;
                    rev = 0;
                }
            }
            if pos + 1 >= k {
                if run >= k {
                    let canon = fwd.min(rev);
                    if idx.insert(canon) {
                        idx.n_kmers += 1;
                    }
                } else {
                    skipped_this_query += 1;
                    idx.n_skipped += 1;
                }
            }
        }

        if skipped_this_query > 0 {
            let n_windows = query.seq.len() - k + 1;
            eprintln!(
                "[warn] query {:?}: {} of {} windows skipped for non-ACGT bases",
                query.name, skipped_this_query, n_windows
            );
        }
    }

    if idx.n_kmers == 0 {
        let names: Vec<&str> = queries.iter().map(|q| q.name.as_str()).collect();
        return Err(Error::Usage(format!(
            "anchor pool {names:?} yielded no usable {k}-mers"
        )));
    }
    Ok(idx)
}

/// Scan every unitig, counting shared canonical k-mers.
///
/// A unitig is an anchor when it shares at least `cfg.min_shared_kmers` canonical k-mers
/// with the query. Counting matched *positions*, not distinct k-mers, to match the Python
/// implementation the anchor-recall study validated against.
///
/// Replaces `find_anchors_blast` (`logan_path_walker.py:69-99`), which writes a BLAST
/// database of the whole accession to disk. That is the largest stage in the Python walker
/// at 5.8 s of a 14.6 s accession.
///
/// # The hot loop
///
/// This runs once per base of every unitig: roughly 68 million iterations for a mid-sized
/// accession. It therefore reads the 2-bit store directly instead of calling
/// `graph.base_at`, which would recompute `seq_offset[id]` on every base, and instead of
/// `decode_seq`, which would allocate and expand each unitig back to one byte per base.
///
/// Forward and reverse-complement registers roll in lockstep so each unitig is scanned
/// once rather than twice. No ambiguity handling is needed: `parse::build` has already
/// guaranteed every stored base is ACGT.
pub fn find(cfg: &Config, graph: &UnitigGraph, index: &QueryIndex) -> Result<Anchors> {
    let k = cfg.k;
    if k != index.k {
        return Err(Error::Usage(format!(
            "config k={k} disagrees with the query index k={}",
            index.k
        )));
    }
    let min_shared = cfg.min_shared_kmers.max(1);
    let shift = 2 * (k as u64 - 1);
    let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };

    let mut ids: Vec<UnitigId> = Vec::new();
    let mut shared_kmers: Vec<u32> = Vec::new();

    for id in 0..graph.n_unitigs {
        let start = graph.seq_offset[id] as usize;
        let end = graph.seq_offset[id + 1] as usize;
        if end - start < k {
            continue;
        }

        let (mut fwd, mut rev) = (0u64, 0u64);
        let mut count: u32 = 0;
        for (j, p) in (start..end).enumerate() {
            // Inlined 2-bit read: byte p/4, quarter p%4.
            let code = ((graph.seq_2bit[p >> 2] >> (2 * (p & 3))) & 3) as u64;
            fwd = ((fwd << 2) | code) & mask;
            rev = (rev >> 2) | ((3 - code) << shift);
            if j + 1 >= k && index.contains(fwd.min(rev)) {
                count += 1;
            }
        }

        if count >= min_shared {
            ids.push(id as UnitigId);
            shared_kmers.push(count);
        }
    }

    let fraction_of_graph = if graph.n_unitigs == 0 {
        0.0
    } else {
        ids.len() as f64 / graph.n_unitigs as f64
    };

    Ok(Anchors {
        ids,
        shared_kmers,
        fraction_of_graph,
    })
}

// The anchor-fraction guard that used to live here has been removed. It was calibrated on
// BLAST anchor counts and then applied to k-mer anchors, so at the default
// `min_shared_kmers = 3` it sat 9x below its own threshold and stayed silent on the one
// accession it existed to catch. `Anchors::fraction_of_graph` is still reported for
// diagnostics. The real guard is `extract::exceeds_guard`, which measures the neighbourhood
// after expansion.

#[cfg(test)]
mod find_tests {
    use super::*;
    use crate::stages::parse;
    use std::path::PathBuf;

    fn graph(fa: &str) -> UnitigGraph {
        parse::build(&"T".to_string(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec()))).unwrap()
    }
    fn query(seq: &str) -> Query {
        Query { name: "q".into(), record_id: "q".into(), seq: seq.as_bytes().to_vec(), path: PathBuf::from("q.fa"), n_records: 1 }
    }
    fn cfg(k: usize, m: u32) -> Config {
        Config { k, min_shared_kmers: m, ..Config::default() }
    }

    #[test]
    fn finds_the_unitig_that_contains_the_query() {
        let c = cfg(5, 1);
        let g = graph(">T_0\nAAAAAAAAAA\n>T_1\nACGTTGCAAC\n>T_2\nCCCCCCCCCC\n");
        let idx = index_query(&c, &query("ACGTTGCAAC")).unwrap();
        let a = find(&c, &g, &idx).unwrap();
        assert_eq!(a.ids, vec![1]);
        assert_eq!(a.shared_kmers, vec![6], "10 bases, k=5 -> 6 windows, all matching");
    }

    /// The query may be present on the opposite strand. Canonical k-mers must find it
    /// without scanning each unitig twice.
    #[test]
    fn finds_a_unitig_stored_as_the_reverse_complement() {
        let c = cfg(5, 1);
        // revcomp("ACGTTGCAAC") = "GTTGCAACGT"
        let g = graph(">T_0\nGTTGCAACGT\n");
        let idx = index_query(&c, &query("ACGTTGCAAC")).unwrap();
        let a = find(&c, &g, &idx).unwrap();
        assert_eq!(a.ids, vec![0]);
        assert_eq!(a.shared_kmers, vec![6]);
    }

    #[test]
    fn threshold_filters_weak_matches() {
        // Unitig 1 shares only its first window with the query.
        let g = graph(">T_0\nACGTTGCAAC\n>T_1\nACGTTAAAAA\n");
        let idx = index_query(&cfg(5, 1), &query("ACGTTGCAAC")).unwrap();
        let loose = find(&cfg(5, 1), &g, &idx).unwrap();
        assert_eq!(loose.ids, vec![0, 1]);
        let strict = find(&cfg(5, 3), &g, &idx).unwrap();
        assert_eq!(strict.ids, vec![0], "unitig 1 must fall below the threshold");
    }

    #[test]
    fn unitigs_shorter_than_k_are_skipped_not_counted() {
        let c = cfg(8, 1);
        let g = graph(">T_0\nACGTACG\n>T_1\nACGTACGTTT\n");
        let idx = index_query(&c, &query("ACGTACGTTT")).unwrap();
        let a = find(&c, &g, &idx).unwrap();
        assert_eq!(a.ids, vec![1], "the 7 bp unitig has no 8-mer window");
    }

    #[test]
    fn ids_come_back_ascending() {
        let c = cfg(4, 1);
        let g = graph(">T_0\nACGTA\n>T_1\nTTTTT\n>T_2\nACGTA\n>T_3\nACGTA\n");
        let idx = index_query(&c, &query("ACGTA")).unwrap();
        let a = find(&c, &g, &idx).unwrap();
        assert_eq!(a.ids, vec![0, 2, 3]);
    }

    #[test]
    fn fraction_of_graph_is_reported() {
        let c = cfg(4, 1);
        let g = graph(">T_0\nACGTA\n>T_1\nTTTTT\n>T_2\nGGGGG\n>T_3\nCCCCC\n");
        let idx = index_query(&c, &query("ACGTA")).unwrap();
        let a = find(&c, &g, &idx).unwrap();
        assert_eq!(a.ids.len(), 1);
        assert!((a.fraction_of_graph - 0.25).abs() < 1e-9);
    }

    #[test]
    fn no_match_is_an_empty_anchor_set_not_an_error() {
        let c = cfg(4, 1);
        let g = graph(">T_0\nAAAAAAAA\n");
        let idx = index_query(&c, &query("CGCGCGCG")).unwrap();
        let a = find(&c, &g, &idx).unwrap();
        assert!(a.ids.is_empty());
        assert_eq!(a.fraction_of_graph, 0.0);
    }

    #[test]
    fn mismatched_k_is_rejected() {
        let g = graph(">T_0\nACGTACGTAC\n");
        let idx = index_query(&cfg(5, 1), &query("ACGTACGTAC")).unwrap();
        assert!(find(&cfg(6, 1), &g, &idx).is_err(), "k mismatch must not silently misread");
    }

    #[test]
    fn lower_k_rescues_divergent_haplotype_anchoring() {
        // Reference query of 100 bp
        let ref_seq = "ACGTGGTCAACCGGTTACGTAGGCTAGCTAGCTAGCTGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA";
        // Mutated target with 15 mutations (85% identity), but clustered in loops so 15-mer motifs survive,
        // while 31-mers are all broken.
        // Positions of mutations: 5, 10, 15, 20, 45, 50, 55, 60, 80, 85, 90, 95
        let mut mut_bytes = ref_seq.as_bytes().to_vec();
        for &pos in &[5, 10, 15, 20, 45, 50, 55, 60, 80, 85, 90, 95] {
            mut_bytes[pos] = match mut_bytes[pos] {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            };
        }
        let mut_seq = std::str::from_utf8(&mut_bytes).unwrap();
        let g = graph(&format!(">T_0\n{}\n", mut_seq));

        // k = 31 fails completely (0 matching 31-mers because largest clean stretch is 24 bp [21..45])
        let c31 = cfg(31, 1);
        let idx31 = index_query(&c31, &query(ref_seq)).unwrap();
        let a31 = find(&c31, &g, &idx31).unwrap();
        assert!(a31.ids.is_empty(), "k=31 should find 0 anchors because clean run is <= 24 bp");

        // k = 15 finds anchors in the clean 24 bp stretch
        let c15 = cfg(15, 1);
        let idx15 = index_query(&c15, &query(ref_seq)).unwrap();
        let a15 = find(&c15, &g, &idx15).unwrap();
        assert_eq!(a15.ids, vec![0], "k=15 must detect anchors in the 24 bp clean run");
        assert!(a15.shared_kmers[0] >= 10);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn q(seq: &str) -> Query {
        Query {
            name: "test".into(),
            record_id: "test".into(),
            seq: seq.as_bytes().to_vec(),
            path: PathBuf::from("test.fasta"),
            n_records: 1,
        }
    }
    fn cfg_k(k: usize) -> Config {
        Config { k, ..Config::default() }
    }

    /// ACGT is its own reverse complement, so forward and reverse agree.
    #[test]
    fn palindromic_kmer_is_its_own_canonical() {
        let idx = index_query(&cfg_k(4), &q("ACGT")).unwrap();
        assert_eq!(idx.n_kmers, 1);
    }

    /// A sequence and its reverse complement must produce the identical index. This is the
    /// property that lets each unitig be scanned once instead of twice.
    #[test]
    fn revcomp_of_the_query_gives_the_same_kmers() {
        let fwd = index_query(&cfg_k(5), &q("ACGTTGCAAGGCTA")).unwrap();
        let rev = index_query(&cfg_k(5), &q("TAGCCTTGCAACGT")).unwrap();
        assert_eq!(fwd.n_kmers, rev.n_kmers);
        for slot in &fwd.table {
            if *slot != EMPTY {
                assert!(rev.contains(*slot), "canonical kmer missing from revcomp index");
            }
        }
    }

    #[test]
    fn counts_distinct_windows() {
        // 10 bases, k=4 -> 7 windows, all distinct here.
        let idx = index_query(&cfg_k(4), &q("ACGTACCTGA")).unwrap();
        assert_eq!(idx.n_kmers, 7);
        assert_eq!(idx.n_skipped, 0);
    }

    #[test]
    fn repeated_kmers_are_deduplicated() {
        // "AAAAAA" with k=3 gives 4 windows, all "AAA".
        let idx = index_query(&cfg_k(3), &q("AAAAAA")).unwrap();
        assert_eq!(idx.n_kmers, 1);
    }

    /// An N must not be coerced to a base; every window touching it is dropped.
    ///
    /// "ACGC N TAGG" with k=4 has exactly two valid windows, one either side of the N.
    /// If N were treated as a base, windows like `CGCN` and `GCNT` would invent k-mers the
    /// query does not contain, and those would then match unitigs it has nothing to do with.
    #[test]
    fn windows_spanning_ambiguity_codes_are_skipped() {
        let dirty = index_query(&cfg_k(4), &q("ACGCNTAGG")).unwrap();
        assert_eq!(dirty.n_skipped, 4, "the four windows touching N must be dropped");

        // The result must equal the union of the two clean flanks, and nothing more.
        let left = index_query(&cfg_k(4), &q("ACGC")).unwrap();
        let right = index_query(&cfg_k(4), &q("TAGG")).unwrap();
        assert_eq!(dirty.n_kmers, left.n_kmers + right.n_kmers);
        for side in [&left, &right] {
            for slot in &side.table {
                if *slot != EMPTY {
                    assert!(dirty.contains(*slot), "flank kmer missing");
                }
            }
        }

        let clean = index_query(&cfg_k(4), &q("ACGTACGT")).unwrap();
        assert_eq!(clean.n_skipped, 0);
    }

    #[test]
    fn membership_is_exact() {
        let idx = index_query(&cfg_k(4), &q("ACGTACCTGA")).unwrap();
        // Every stored kmer is found.
        let stored: Vec<u64> = idx.table.iter().copied().filter(|&s| s != EMPTY).collect();
        for s in &stored {
            assert!(idx.contains(*s));
        }
        // A kmer that cannot be in a 4-mer index is not found.
        assert!(!idx.contains(0xDEAD_BEEF_DEAD_BEEF & ((1 << 8) - 1) | 0xF0));
    }

    #[test]
    fn rejects_impossible_k() {
        assert!(index_query(&cfg_k(0), &q("ACGT")).is_err());
        assert!(index_query(&cfg_k(33), &q("ACGT")).is_err());
    }

    #[test]
    fn rejects_query_shorter_than_k() {
        assert!(index_query(&cfg_k(31), &q("ACGTACGT")).is_err());
    }

    #[test]
    fn all_ambiguous_query_is_an_error_not_an_empty_index() {
        assert!(index_query(&cfg_k(4), &q("NNNNNNNN")).is_err());
    }

    #[test]
    fn k32_does_not_overflow_the_mask() {
        let seq = "ACGT".repeat(16); // 64 bases
        let idx = index_query(&cfg_k(32), &q(&seq)).unwrap();
        assert!(idx.n_kmers > 0);
    }

    // ---- index_query_set (multi-anchor pool) ----

    /// The whole point of the wrapper: a single-query `index_query_set` produces exactly the
    /// same table as `index_query`. If this drifts, `bench/check_golden.sh` will fail because
    /// anchor ids depend on the table's capacity and hash slots.
    #[test]
    fn single_query_set_is_bit_identical_to_index_query() {
        let seq = "ACGTGGTCAACCGGTTACGTAGGCTAGCTAGCTAGCTGATCGATCGATCGATC";
        let cfg = cfg_k(11);
        let solo = index_query(&cfg, &q(seq)).unwrap();
        let set = index_query_set(&cfg, std::slice::from_ref(&q(seq))).unwrap();
        assert_eq!(solo.k, set.k);
        assert_eq!(solo.mask, set.mask);
        assert_eq!(solo.n_kmers, set.n_kmers);
        assert_eq!(solo.n_skipped, set.n_skipped);
        assert_eq!(solo.table.len(), set.table.len(), "capacity must match");
        assert_eq!(solo.table, set.table, "slots must be bit-identical");
    }

    /// Two queries union: every k-mer of either query is present, and duplicates counted once.
    #[test]
    fn set_of_two_queries_unions_the_kmer_sets() {
        let cfg = cfg_k(5);
        // Overlap: "ACGTA" appears in both.
        let a = q("ACGTACGTAT");
        let b = q("ACGTAGGGCC");
        let ia = index_query(&cfg, &a).unwrap();
        let ib = index_query(&cfg, &b).unwrap();
        let union = index_query_set(&cfg, &[a, b]).unwrap();

        // Every k-mer of either single-query index must be in the union.
        for src in [&ia, &ib] {
            for &slot in &src.table {
                if slot != EMPTY {
                    assert!(union.contains(slot), "union missing kmer {slot:#x}");
                }
            }
        }
        // The dedupe invariant: union count is at most sum of individual counts, and at least
        // the size of the larger of the two.
        assert!(union.n_kmers <= ia.n_kmers + ib.n_kmers);
        assert!(union.n_kmers >= ia.n_kmers.max(ib.n_kmers));
    }

    /// Aggregate empty-index error names every query so a bad record is identifiable.
    #[test]
    fn set_reports_all_query_names_when_pool_yields_nothing() {
        let cfg = cfg_k(4);
        let err = match index_query_set(&cfg, &[q("NNNN"), q("NNNN")]) {
            Ok(_) => panic!("expected all-ambiguous pool to error"),
            Err(e) => e.to_string(),
        };
        assert!(err.contains("no usable"), "{err}");
    }

    #[test]
    fn set_rejects_empty_pool() {
        let cfg = cfg_k(4);
        assert!(index_query_set(&cfg, &[]).is_err());
    }

    #[test]
    fn set_rejects_a_query_shorter_than_k() {
        let cfg = cfg_k(21);
        assert!(index_query_set(&cfg, &[q("ACGT")]).is_err());
    }
}

/// Cross-validation against the independent Python implementation used for the anchor-recall
/// study (`scratchpad/kmer_anchor.py`, canonical 31-mers via `min(m, revcomp(m))`). The
/// expected counts come from that study's recorded output, so a mismatch means the two
/// implementations disagree about what a canonical k-mer is.
///
/// Run from the crate directory: `cargo test --release -- --ignored --nocapture`
#[cfg(test)]
mod python_concordance {
    use super::*;
    use crate::io::input::read_query;
    use std::path::Path;

    #[test]
    #[ignore]
    fn matches_python_canonical_kmer_counts() {
        // (query file, canonical 31-mer count measured by the Python implementation)
        let cases = [
            ("../queries/pfcrt_only.fasta", 1245usize),
            ("../queries/query_hiv1_pol.fasta", 2709),
            ("../queries/markers/candida_auris_erg11.fasta", 1545),
        ];
        let cfg = Config { k: 31, ..Config::default() };
        assert_eq!(cfg.k, 31, "these counts are for k=31");

        for (path, expected) in cases {
            let p = Path::new(path);
            if !p.exists() {
                eprintln!("SKIP: {path} not found");
                continue;
            }
            let query = read_query(p).expect("read query");
            let idx = index_query(&cfg, &query).expect("index");
            eprintln!(
                "  {:<44} {} bp -> {} canonical 31-mers (python: {})",
                query.name,
                query.seq.len(),
                idx.n_kmers,
                expected
            );
            assert_eq!(idx.n_kmers, expected, "disagreement on {path}");
            // Every window was usable, so the count equals seq_len - k + 1.
            assert_eq!(idx.n_skipped, 0);
            assert_eq!(idx.n_kmers, query.seq.len() - cfg.k + 1);
        }
    }
}

/// Cross-validation against the Python anchor-recall study. Those anchor counts were
/// produced by an independent implementation (`scratchpad/kmer_anchor.py`) and then used to
/// decide that k-mer anchoring is viable, so the two must agree exactly.
///
/// `LOGAN_TEST_ZST=/path/to/ERR2028820.unitigs.fa.zst cargo test --release -- --ignored --nocapture`
#[cfg(test)]
mod python_concordance_find {
    use super::*;
    use crate::io::input::read_query;
    use crate::stages::{parse, stream};
    use crate::types::FetchedUnitigs;
    use std::path::Path;

    #[test]
    #[ignore]
    fn matches_python_anchor_counts_at_every_threshold() {
        let Ok(path) = std::env::var("LOGAN_TEST_ZST") else {
            eprintln!("SKIP: set LOGAN_TEST_ZST to ERR2028820.unitigs.fa.zst");
            return;
        };
        let Ok(compressed) = std::fs::read(&path) else {
            eprintln!("SKIP: cannot read {path}");
            return;
        };
        let qpath = Path::new("../queries/pfcrt_only.fasta");
        if !qpath.exists() {
            eprintln!("SKIP: {} not found", qpath.display());
            return;
        }

        let fetched = FetchedUnitigs {
            index: 0,
            accession: "ERR2028820".to_string(),
            compressed_len: compressed.len(),
            compressed,
            fetch_secs: 0.0,
        };
        let graph = parse::build(&fetched.accession, stream::decoder(&fetched).unwrap()).unwrap();
        let query = read_query(qpath).unwrap();

        // (min_shared_kmers, anchors measured by the Python implementation)
        let expected = [(1u32, 104usize), (2, 83), (3, 77), (5, 59), (10, 36)];
        for (m, want) in expected {
            let cfg = Config { min_shared_kmers: m, ..Config::default() };
            let idx = index_query(&cfg, &query).unwrap();
            let t0 = std::time::Instant::now();
            let a = find(&cfg, &graph, &idx).unwrap();
            eprintln!(
                "  M={:<3} anchors {:>5} (python {:>5})  {:.3}s  fraction {:.6}%",
                m,
                a.ids.len(),
                want,
                t0.elapsed().as_secs_f64(),
                a.fraction_of_graph * 100.0
            );
            assert_eq!(a.ids.len(), want, "anchor count disagreement at M={m}");
            assert_eq!(a.ids.len(), a.shared_kmers.len());
            assert!(a.shared_kmers.iter().all(|&c| c >= m), "a kept anchor is below the threshold");
        }

        // BLAST found 76 anchors on this accession; the k-mer set at M=3 is a comparable size.
        let cfg = Config { k: 31, ..Config::default() };
        let idx = index_query(&cfg, &query).unwrap();
        let a = find(&cfg, &graph, &idx).unwrap();
        assert_eq!(a.ids.len(), 77);
        // fraction_of_graph is diagnostic only; the guard lives in extract, after expansion.
        assert!(a.fraction_of_graph < 0.001);
    }
}
