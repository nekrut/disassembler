//! Step 3: find the unitigs that carry the query gene.
//!
//! Shared canonical k-mers, not a BLAST database of every unitig written to disk and
//! searched with megablast. Under BLAST that was the single largest stage, 5.8 s of a
//! 14.6 s accession.
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
/// with the query. Counting matched *positions*, not distinct k-mers, as the anchor-recall
/// study did.
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
/// Counts shared k-mers as a record streams past, so anchoring stores nothing.
///
/// This is `parse::RecordSink`, driven by `parse::scan`. Fusing anchoring into the decode is
/// what lets the parser avoid storing every unitig's sequence: 68 million bases go past once
/// and only the anchor ids survive. See `count_shared` for the register trick.
pub struct Scanner<'a> {
    k: usize,
    min_shared: u32,
    index: &'a QueryIndex,
    shift: u64,
    mask: u64,
    /// Rolling registers for the record in progress.
    fwd: u64,
    rev: u64,
    /// Bases seen in this record, which is also the k-mer boundary test.
    j: usize,
    count: u32,
    id: UnitigId,
    ids: Vec<UnitigId>,
    shared_kmers: Vec<u32>,
    n_unitigs: usize,
}

impl<'a> Scanner<'a> {
    pub fn new(cfg: &Config, index: &'a QueryIndex) -> Result<Self> {
        let k = cfg.k;
        if k != index.k {
            return Err(Error::Usage(format!(
                "config k={k} disagrees with the query index k={}",
                index.k
            )));
        }
        Ok(Self {
            k,
            min_shared: cfg.min_shared_kmers.max(1),
            index,
            shift: 2 * (k as u64 - 1),
            mask: if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 },
            fwd: 0,
            rev: 0,
            j: 0,
            count: 0,
            id: 0,
            ids: Vec::new(),
            shared_kmers: Vec::new(),
            n_unitigs: 0,
        })
    }

    /// The anchors found, and what share of the graph they are.
    pub fn finish(self) -> Anchors {
        let fraction_of_graph = if self.n_unitigs == 0 {
            0.0
        } else {
            self.ids.len() as f64 / self.n_unitigs as f64
        };
        Anchors {
            ids: self.ids,
            shared_kmers: self.shared_kmers,
            fraction_of_graph,
        }
    }
}

impl crate::stages::parse::RecordSink for Scanner<'_> {
    #[inline]
    fn begin(&mut self, id: UnitigId) {
        self.id = id;
        self.fwd = 0;
        self.rev = 0;
        self.j = 0;
        self.count = 0;
        self.n_unitigs += 1;
    }

    /// The hot loop: once per base of every unitig, roughly 68 million times for a mid-sized
    /// accession. Monomorphised through `parse::scan`, so this inlines into the decode.
    #[inline]
    fn base(&mut self, code: u8) {
        let c = code as u64;
        self.fwd = ((self.fwd << 2) | c) & self.mask;
        self.rev = (self.rev >> 2) | ((3 - c) << self.shift);
        self.j += 1;
        // A unitig shorter than k never reaches this, which is the old loop's explicit skip.
        if self.j >= self.k && self.index.contains(self.fwd.min(self.rev)) {
            self.count += 1;
        }
    }

    #[inline]
    fn end(&mut self) {
        if self.count >= self.min_shared {
            self.ids.push(self.id);
            self.shared_kmers.push(self.count);
        }
    }
}

/// Anchors of a graph that already holds every sequence.
///
/// The production path does not use this: it anchors during `parse::scan` with a [`Scanner`],
/// which is the point of the two-pass parser. Kept because a `build`-made graph is the
/// convenient fixture for tests, and because driving the same `Scanner` here is what makes
/// those tests evidence about the production path rather than about a second implementation.
pub fn find(cfg: &Config, graph: &UnitigGraph, index: &QueryIndex) -> Result<Anchors> {
    use crate::stages::parse::RecordSink;
    let mut scanner = Scanner::new(cfg, index)?;
    let mut seq: Vec<u8> = Vec::new();
    for id in 0..graph.n_unitigs as UnitigId {
        scanner.begin(id);
        graph.decode_seq(id, &mut seq);
        for &b in &seq {
            scanner.base(match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                _ => 3,
            });
        }
        scanner.end();
    }
    Ok(scanner.finish())
}

/// The canonical k-mers of one sequence, sorted and deduplicated.
///
/// Same definition of "canonical" as [`index_query_set`]: `min(kmer, revcomp(kmer))`, with any
/// window holding a non-ACGT base skipped. The rolling registers are written out again here
/// rather than shared, because the indexer carries per-query bookkeeping this has no use for;
/// `the_kmer_set_agrees_with_the_index` pins the two together so they cannot drift.
///
/// Used by [`reference_groups`], which needs each reference's k-mers separately. The index is
/// the union and cannot answer which reference a k-mer came from.
pub fn canonical_kmers(k: usize, seq: &[u8]) -> Vec<u64> {
    let mut out: Vec<u64> = Vec::new();
    if k == 0 || k > 32 || seq.len() < k {
        return out;
    }
    let shift = 2 * (k as u64 - 1);
    let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let (mut fwd, mut rev, mut run) = (0u64, 0u64, 0usize);
    for &b in seq {
        match base_code(b) {
            Some(code) => {
                fwd = ((fwd << 2) | code) & mask;
                rev = (rev >> 2) | ((3 - code) << shift);
                run += 1;
            }
            None => {
                run = 0;
                fwd = 0;
                rev = 0;
            }
        }
        if run >= k {
            out.push(fwd.min(rev));
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// Elements two sorted, deduplicated k-mer lists have in common.
fn shared(a: &[u64], b: &[u64]) -> usize {
    let (mut i, mut j, mut n) = (0usize, 0usize, 0usize);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                n += 1;
                i += 1;
                j += 1;
            }
        }
    }
    n
}

/// The share of the smaller reference's k-mers two references must have in common to count
/// as the same gene.
///
/// 1%. Strains of one gene share far more than this (the SARS-CoV-2 spike panel is one group
/// at any threshold up to 50%), and unrelated genes share far less: two random coding
/// sequences share none, and the 142-record, 25-gene master target file splits into 31
/// groups. The gap is wide enough that the exact value does not matter much, and low enough
/// that a divergent but legitimate panel is not split.
pub const GROUP_MIN_SHARE: f64 = 0.01;

/// Partition the anchor pool into groups of references that look like the same gene.
///
/// Multi-anchor mode ranks every `-q` record against one another and keeps a single
/// primary per accession, which is right for strains of one gene and wrong for a panel of
/// different genes: an accession carrying both `bla_KPC` and `bla_NDM` comes back with one of
/// them, exit 0, and the other visible only as an unselected candidate in the summary. The
/// walker cannot tell which the user meant, so it does not refuse; it says what it sees.
///
/// Two references are linked when they share at least `min_share` of the smaller one's
/// canonical k-mers, and groups are the connected components of that relation. The fraction
/// is of the smaller set so that a short reference and a long one holding it are linked.
///
/// Returned groups are sorted, each ascending by index, so the warning text is deterministic.
///
/// Cost is one k-mer list per reference and a merge per pair. The widest panel in the
/// repository is 142 records, which is about 10,000 merges of about 1,500 k-mers each.
pub fn reference_groups(k: usize, seqs: &[&[u8]], min_share: f64) -> Vec<Vec<usize>> {
    let n = seqs.len();
    let kmers: Vec<Vec<u64>> = seqs.iter().map(|s| canonical_kmers(k, s)).collect();
    // Union-find, by index, with path halving.
    let mut parent: Vec<usize> = (0..n).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    for i in 0..n {
        for j in i + 1..n {
            let smaller = kmers[i].len().min(kmers[j].len());
            if smaller == 0 {
                continue;
            }
            if shared(&kmers[i], &kmers[j]) as f64 >= min_share * smaller as f64 {
                let (a, b) = (root(&mut parent, i), root(&mut parent, j));
                if a != b {
                    parent[a] = b;
                }
            }
        }
    }
    let mut groups: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
    for i in 0..n {
        let r = root(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    let mut out: Vec<Vec<usize>> = groups.into_values().collect();
    out.sort();
    out
}

// There is deliberately no anchor-fraction guard here. One calibrated on BLAST anchor
// counts sits 9x below its own threshold once applied to k-mer anchors, so it stays silent
// on the accession it exists to catch. `Anchors::fraction_of_graph` is reported for
// diagnostics; the guard is `extract::exceeds_guard`, which measures after expansion.

#[cfg(test)]
mod group_tests {
    use super::*;
    use std::path::PathBuf;

    fn q(seq: &str) -> Query {
        Query {
            name: "q".into(),
            record_id: "q".into(),
            seq: seq.as_bytes().to_vec(),
            path: PathBuf::from("q.fa"),
            n_records: 1,
        }
    }

    /// `canonical_kmers` rolls the same registers `index_query_set` does, written out a second
    /// time because the indexer carries bookkeeping this has no use for. If the two ever
    /// disagreed about what a canonical k-mer is, the grouping would be measuring something
    /// the anchoring does not use.
    #[test]
    fn the_kmer_set_agrees_with_the_index() {
        let seq = "ACGTTGCAACGGATCCGTANACGTTTGCAGCATCGGATCAGGCATCAG";
        let c = Config { k: 7, ..Config::default() };
        let idx = index_query(&c, &q(seq)).unwrap();
        let set = canonical_kmers(c.k, seq.as_bytes());
        assert_eq!(set.len(), idx.len(), "same number of distinct canonical k-mers");
        for km in &set {
            assert!(idx.contains(*km), "index is missing a k-mer the set found");
        }
        // The N must have cost both of them the same windows.
        assert!(idx.n_skipped > 0, "the fixture must exercise the ambiguity path");
    }

    /// A sequence and a 1%-mutated copy of it are one gene; two unrelated sequences are two.
    #[test]
    fn related_references_group_and_unrelated_ones_do_not() {
        // Deterministic pseudo-random sequences, so the test cannot flake.
        let gen = |seed: u64, n: usize| -> Vec<u8> {
            let mut x = seed;
            (0..n)
                .map(|_| {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    b"ACGT"[((x >> 33) & 3) as usize]
                })
                .collect()
        };
        let a = gen(1, 900);
        let b = gen(2, 900);
        let mut a_mut = a.clone();
        for i in (0..a_mut.len()).step_by(100) {
            a_mut[i] = if a_mut[i] == b'A' { b'C' } else { b'A' };
        }

        let two = reference_groups(21, &[&a[..], &b[..]], GROUP_MIN_SHARE);
        assert_eq!(two.len(), 2, "unrelated sequences must not be one gene: {two:?}");

        let one = reference_groups(21, &[&a[..], &a_mut[..]], GROUP_MIN_SHARE);
        assert_eq!(one.len(), 1, "a 1%-mutated copy is the same gene: {one:?}");

        // Grouping is transitive through a shared member, which is what makes a strain panel
        // one group even when its extremes are divergent.
        let three = reference_groups(21, &[&a[..], &a_mut[..], &b[..]], GROUP_MIN_SHARE);
        assert_eq!(three, vec![vec![0, 1], vec![2]]);
    }

    /// The repository's own reference files, which are what a user is told to start from: the
    /// strain panel must be one gene and the master target file must not. Skips, loudly,
    /// where the files are not beside the walker (a packaged build).
    #[test]
    fn the_shipped_panels_group_the_way_they_read() {
        let read = |path: &str| -> Option<Vec<Vec<u8>>> {
            let text = std::fs::read_to_string(path).ok()?;
            let mut out: Vec<Vec<u8>> = Vec::new();
            for line in text.lines() {
                if line.starts_with('>') {
                    out.push(Vec::new());
                } else if let Some(last) = out.last_mut() {
                    last.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
                }
            }
            Some(out)
        };

        match read("../panels/sarscov2_spike_11strains.fa") {
            Some(recs) if recs.len() > 1 => {
                let seqs: Vec<&[u8]> = recs.iter().map(|r| r.as_slice()).collect();
                let groups = reference_groups(21, &seqs, GROUP_MIN_SHARE);
                assert_eq!(
                    groups.len(),
                    1,
                    "the spike strain panel is one gene and must not warn: {groups:?}"
                );
            }
            _ => eprintln!("SKIP: panels/sarscov2_spike_11strains.fa not readable from here"),
        }

        match read("../queries/master_planetary_targets.fasta") {
            Some(recs) if recs.len() > 1 => {
                let seqs: Vec<&[u8]> = recs.iter().map(|r| r.as_slice()).collect();
                let groups = reference_groups(21, &seqs, GROUP_MIN_SHARE);
                eprintln!("master_planetary_targets.fasta: {} records, {} groups", recs.len(), groups.len());
                // 25 gene prefixes, but 31 k-mer-sharing groups: some genes split.
                assert_eq!((recs.len(), groups.len()), (142, 31));
                assert!(
                    groups.len() > 1,
                    "the master target file is a catalogue of different genes, not a -q panel"
                );
            }
            _ => eprintln!("SKIP: queries/master_planetary_targets.fasta not readable from here"),
        }
    }
}

#[cfg(test)]
mod find_tests {
    use super::*;
    use crate::stages::parse;
    use std::path::PathBuf;

    fn graph(fa: &str) -> UnitigGraph {
        parse::build(&"T".to_string(), Box::new(std::io::Cursor::new(fa.as_bytes().to_vec())), u64::MAX).unwrap()
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

/// Cross-validation against the independent implementation used for the anchor-recall study
/// (canonical 31-mers via `min(m, revcomp(m))`). The expected counts come from that
/// study's recorded output, so a mismatch means the two implementations disagree about what
/// a canonical k-mer is.
///
/// Run from the crate directory: `cargo test --release -- --ignored --nocapture`
#[cfg(test)]
mod recall_study_kmers {
    use super::*;
    use crate::io::input::read_query;
    use std::path::Path;

    /// Not `#[ignore]`d. It reads only repository files, so there is nothing to skip
    /// for, and while it was ignored it never ran on push.
    #[test]
    fn matches_the_recall_study_kmer_counts() {
        // (query file, canonical 31-mer count the study recorded)
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
                "  {:<44} {} bp -> {} canonical 31-mers (study: {})",
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

/// Cross-validation against the anchor-recall study, whose counts came from an independent
/// implementation and decided that k-mer anchoring is viable, so the two must agree exactly.
///
/// `LOGAN_TEST_ZST=/path/to/ERR2028820.unitigs.fa.zst cargo test --release -- --ignored --nocapture`
#[cfg(test)]
mod recall_study_anchors {
    use super::*;
    use crate::io::input::read_query;
    use crate::stages::{parse, stream};
    use crate::types::FetchedUnitigs;
    use std::path::Path;

    #[test]
    #[ignore]
    fn matches_the_recall_study_anchor_counts_at_every_threshold() {
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
            compressed: crate::types::Compressed::Memory(compressed),
            fetch_secs: 0.0,
            payload: crate::types::Payload::default(),
        };
        let graph = parse::build(&fetched.accession, stream::decoder(&fetched).unwrap(), u64::MAX).unwrap();
        let query = read_query(qpath).unwrap();

        // (min_shared_kmers, anchors the study recorded)
        let expected = [(1u32, 104usize), (2, 83), (3, 77), (5, 59), (10, 36)];
        for (m, want) in expected {
            // The study counted canonical 31-mers; the default -k has since moved to 21.
            let cfg = Config { k: 31, min_shared_kmers: m, ..Config::default() };
            let idx = index_query(&cfg, &query).unwrap();
            let t0 = std::time::Instant::now();
            let a = find(&cfg, &graph, &idx).unwrap();
            eprintln!(
                "  M={:<3} anchors {:>5} (study {:>5})  {:.3}s  fraction {:.6}%",
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
