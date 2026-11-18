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
    /// File stem, e.g. `query_hiv1_pol`. Used for output filenames: the summary is named
    /// from the query file's stem, and that name is part of the output contract.
    pub name: String,
    /// First token of the FASTA header, e.g. `pfcrt`. Provenance and logging only.
    pub record_id: String,
    /// Uppercased sequence with all whitespace removed.
    pub seq: Vec<u8>,
    /// Path on disk. Provenance, warnings and error messages only: nothing re-opens it once
    /// the query has been read, and GraphAligner is never handed it. It picks its parser from
    /// the file name, which is the user's to choose and may be anything (a Galaxy dataset
    /// arrives as `dataset_<uuid>.dat`).
    ///
    /// Until 0.2.0 `align::run` re-read this path on every accession to hand GraphAligner a
    /// byte-for-byte copy. That emptied a piped query after the first accession and passed
    /// the user's header line and any whitespace inside the sequence through to an aligner
    /// that chokes on both; `align::run` now writes `seq` under `record_id` instead.
    pub path: PathBuf,
    /// Records present in the source file. `main` uses it to decide whether a single `-q`
    /// should be read as a whole reference set (multi-anchor mode) rather than as just its
    /// first record.
    pub n_records: usize,
}

// ---------------------------------------------------------------------------
// provenance: what was read, and where it came from
// ---------------------------------------------------------------------------

/// Lowercase hex sha256 of `bytes`.
///
/// One digest helper for the whole walker: `io::output` names a reconstructed sequence with
/// it, and `stages::stream` names the compressed payload an accession was built from.
pub fn sha256_hex(bytes: &[u8]) -> String {
    use sha2::{Digest, Sha256};
    let digest = Sha256::digest(bytes);
    let mut s = String::with_capacity(64);
    for b in digest {
        s.push_str(&format!("{b:02x}"));
    }
    s
}

/// The whole graph's adjacency, held on disk instead of in RAM.
///
/// Adjacency is the last part of the walker that is O(accession) rather than O(neighbourhood):
/// 12 bytes per unitig, 15.4 MB for `ERR2028820`, and it exists from `parse::scan` until
/// `extract::expand` has finished walking it. It is also the term that keeps `--max-size-mb`
/// in the memory bound, because it scales with the file rather than with the gene.
///
/// Two unlinked temp files, read with `pread`. Unlinked at creation, so the data is reclaimed
/// when the handles close even if the process is killed: there is no cleanup path to get
/// wrong and nothing to sweep.
///
/// Spilling is opt-in per accession, past [`SPILL_ADJACENCY_AT`]. A small accession stays in
/// RAM and never touches the disk, which keeps the common case free of a new failure mode, and
/// a spill that cannot be created falls back to memory with a warning rather than failing the
/// accession.
#[derive(Debug)]
pub struct AdjacencySpool {
    offsets: std::fs::File,
    targets: std::fs::File,
    n_offsets: usize,
    n_targets: usize,
}

/// In-memory adjacency above this spills to disk. 4 MiB is about 350,000 unitigs.
///
/// Low enough that any accession big enough for adjacency to matter spills, high enough that
/// the small ones, and every test fixture, stay in RAM.
pub const SPILL_ADJACENCY_AT: usize = 4 << 20;

impl AdjacencySpool {
    /// Two unlinked files under `dir`.
    pub fn create(dir: &std::path::Path) -> std::io::Result<Self> {
        Ok(Self {
            offsets: anon_file(dir, "off")?,
            targets: anon_file(dir, "tgt")?,
            n_offsets: 0,
            n_targets: 0,
        })
    }

    pub fn push_offsets(&mut self, v: &[u32]) -> std::io::Result<()> {
        use std::io::Write;
        self.offsets.write_all(as_bytes(v))?;
        self.n_offsets += v.len();
        Ok(())
    }

    pub fn push_targets(&mut self, v: &[u32]) -> std::io::Result<()> {
        use std::io::Write;
        self.targets.write_all(as_bytes(v))?;
        self.n_targets += v.len();
        Ok(())
    }

    pub fn n_targets(&self) -> usize {
        self.n_targets
    }

    /// Row `id` into `buf`, which is cleared first.
    pub fn edges_into(&self, id: UnitigId, buf: &mut Vec<u32>) -> std::io::Result<()> {
        use std::os::unix::fs::FileExt;
        let i = id as usize;
        if i + 1 >= self.n_offsets {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("unitig {id} is outside the spooled graph ({} rows)", self.n_offsets),
            ));
        }
        let mut bounds = [0u8; 8];
        self.offsets.read_exact_at(&mut bounds, (i * 4) as u64)?;
        let lo = u32::from_ne_bytes([bounds[0], bounds[1], bounds[2], bounds[3]]) as usize;
        let hi = u32::from_ne_bytes([bounds[4], bounds[5], bounds[6], bounds[7]]) as usize;
        if hi < lo || hi > self.n_targets {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("spooled row {id} is [{lo}, {hi}) of {} targets", self.n_targets),
            ));
        }
        buf.clear();
        buf.resize(hi - lo, 0);
        if hi > lo {
            self.targets.read_exact_at(as_bytes_mut(buf), (lo * 4) as u64)?;
        }
        Ok(())
    }
}

/// Bytes an unprivileged process can still write under `dir`, or `None` if it cannot be asked.
///
/// `f_bavail`, not `f_bfree`: they differ by the blocks reserved for root, and the walker is
/// not root. The spool's files are unlinked the moment they are created, so their blocks are
/// still charged to the filesystem and this figure already accounts for every live spool,
/// including other workers' and other processes'. What it cannot see is growth that has not
/// happened yet, which is why callers multiply by the number of concurrent spools.
///
/// Walks up to the nearest existing ancestor, because a spool directory is created lazily by
/// the first spill and may not exist when the bound is printed or the gate is checked. `None`
/// rather than 0 on failure, so a filesystem that will not answer disables the check instead
/// of refusing everything on it.
pub fn available_bytes(dir: &std::path::Path) -> Option<u64> {
    use std::os::unix::ffi::OsStrExt;
    let mut p = dir;
    while !p.exists() {
        p = p.parent()?;
    }
    let c = std::ffi::CString::new(p.as_os_str().as_bytes()).ok()?;
    // SAFETY: `c` is a NUL-terminated path that outlives the call, and `st` is a correctly
    // sized `statvfs` that is only read after the call reports success.
    let mut st: libc::statvfs = unsafe { std::mem::zeroed() };
    if unsafe { libc::statvfs(c.as_ptr(), &mut st) } != 0 {
        return None;
    }
    (st.f_bavail as u64).checked_mul(st.f_frsize as u64)
}

/// A file that exists only as an open handle, created under `dir` and immediately unlinked.
///
/// Two callers: the adjacency spool, and a downloaded object with no `--unitigs-cache` to land
/// in. Both want the same thing, which is bytes that do not live in RAM and do not survive the
/// process. Nothing is left behind by a crash or a signal, and nothing needs sweeping.
pub fn transient_file(dir: &std::path::Path, tag: &str) -> std::io::Result<std::fs::File> {
    anon_file(dir, tag)
}

/// A file that is unlinked the moment it exists, so it needs no cleanup and survives no crash.
fn anon_file(dir: &std::path::Path, tag: &str) -> std::io::Result<std::fs::File> {
    std::fs::create_dir_all(dir)?;
    let ns = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let path = dir.join(format!(".logan-{tag}-{}-{ns}.spool", std::process::id()));
    let f = std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create_new(true)
        .open(&path)?;
    // From here the only reference is the handle. The bytes go back to the filesystem when it
    // closes, including on a signal, and nothing is left behind to sweep.
    std::fs::remove_file(&path)?;
    Ok(f)
}

#[inline]
fn as_bytes(v: &[u32]) -> &[u8] {
    // Safe: u32 has no padding and no invalid bit patterns, and the slice is read-only. The
    // spool is written and read by the same process, so native endianness is consistent.
    unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) }
}

#[inline]
fn as_bytes_mut(v: &mut [u32]) -> &mut [u8] {
    // Safe: as `as_bytes`, and every byte pattern is a valid u32.
    unsafe { std::slice::from_raw_parts_mut(v.as_mut_ptr() as *mut u8, std::mem::size_of_val(v)) }
}

/// Which ingest path produced an accession's bytes.
///
/// This answers "where did these bytes come from", not "what was configured": a run with
/// `--unitigs-cache` and no `--offline` reports `Cache` for the accessions that hit and `S3`
/// for the ones it had to download, which `inputs.source` (`cache+s3`) cannot distinguish.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PayloadSource {
    /// Downloaded from the Logan bucket during this run.
    S3,
    /// Read from `--unitigs-cache`.
    Cache,
    /// No bytes were read: a 404, a transport failure, or an error before the fetch.
    #[default]
    None,
}

impl PayloadSource {
    pub fn as_str(self) -> &'static str {
        match self {
            PayloadSource::S3 => "s3",
            PayloadSource::Cache => "cache",
            PayloadSource::None => "none",
        }
    }
}

/// What one accession's unitig object was, as opposed to what the walker made of it.
///
/// The Logan objects at `logan-pub.s3.amazonaws.com/u/{acc}/` are unversioned and mutable,
/// and they do get replaced: a 5,000-key sample across five ranges found 3.1% of objects
/// rewritten after the 2024-03 bulk load, in waves on 2025-12-10, 2025-12-14 and 2026-02-08.
/// Until 0.2.0 a run recorded nothing about the bytes it consumed, so a rerun after an
/// upstream regeneration produced different results with no detectable difference in any
/// output. This is the record that makes those two runs tell apart.
///
/// `sha256` is the field that works on every path, because it is derived from the bytes
/// rather than from the transport: a cache hit and an `--offline` run have no HTTP response
/// to read headers from. `etag` and `last_modified` are recorded in addition when a live
/// fetch supplied them, because they are comparable with what the bucket serves today.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Payload {
    pub source: PayloadSource,
    /// Compressed size. Known without reading the body when S3 sent a `Content-Length` and
    /// the object was then refused for being over `--max-size-mb`.
    pub bytes: Option<u64>,
    /// sha256 of the compressed bytes exactly as read, which are byte-identical to the S3
    /// object and to the cache file. `None` when no bytes were read.
    pub sha256: Option<String>,
    /// S3 `ETag`, surrounding quotes stripped. For every object below S3's 5 GB single-PUT
    /// limit this is the MD5 of the content, so it is comparable with the live bucket; above
    /// it the value carries a `-<parts>` suffix and is opaque. `--max-size-mb` defaults to
    /// 1024, so in practice every object the walker ingests has the comparable form.
    pub etag: Option<String>,
    /// S3 `Last-Modified`, normalised to `YYYY-MM-DDTHH:MM:SSZ`. Omitted rather than stored
    /// raw when it cannot be parsed, so the field has one format or none.
    pub last_modified: Option<String>,
    /// When the object was read from S3.
    ///
    /// Set on a cache hit too, from the sidecar, since 0.2.0: that is the clock that says how
    /// old the bytes are, which is the question a warm run actually raises. Absent for a cache
    /// entry written before sidecars existed, where nothing records it.
    pub fetched_at: Option<String>,
    /// What `--revalidate` found, when it ran. `None` means it did not.
    ///
    /// - `current`: the bucket confirmed the cached object (HTTP 304).
    /// - `replaced`: the bucket serves different bytes now, and these are the new ones.
    /// - `refetched`: nothing recorded what the cached bytes were a copy of, so they were
    ///   downloaded again unconditionally rather than trusted.
    /// - `unavailable`: the check could not be made, and the cached bytes were used anyway.
    ///
    /// `replaced` is the one worth alerting on: it is a silent substitution made visible, and
    /// it means an earlier run of the same command analysed different data.
    pub revalidated: Option<&'static str>,
}

impl Payload {
    /// Bytes that were read, from `source`, with their digest.
    pub fn read(source: PayloadSource, bytes: &[u8]) -> Self {
        Self {
            source,
            bytes: Some(bytes.len() as u64),
            sha256: Some(sha256_hex(bytes)),
            ..Self::default()
        }
    }

    /// The same for a file whose contents are never held: the digest is taken in chunks and
    /// the handle is rewound, so this costs 128 KiB of stack and one extra pass over bytes the
    /// OS page cache is about to be asked for anyway.
    pub fn read_file(source: PayloadSource, file: &mut std::fs::File) -> std::io::Result<Self> {
        use sha2::{Digest, Sha256};
        use std::io::{Read, Seek, SeekFrom};
        file.seek(SeekFrom::Start(0))?;
        let mut hasher = Sha256::new();
        let mut buf = vec![0u8; 128 * 1024];
        let mut total: u64 = 0;
        loop {
            let n = file.read(&mut buf)?;
            if n == 0 {
                break;
            }
            hasher.update(&buf[..n]);
            total += n as u64;
        }
        file.seek(SeekFrom::Start(0))?;
        let mut hex = String::with_capacity(64);
        for b in hasher.finalize() {
            hex.push_str(&format!("{b:02x}"));
        }
        Ok(Self {
            source,
            bytes: Some(total),
            sha256: Some(hex),
            ..Self::default()
        })
    }
}

/// What a cache entry records about the object it was downloaded from.
///
/// A cached file carries the digest of its own bytes, which says they are intact but not what
/// they are a copy of. Logan's `u/` objects are mutable, so without this a warm run cannot
/// name the bucket object its bytes came from, and `--revalidate` has no `ETag` to send.
///
/// Written beside the object as `<name>.provenance`, one `key value` per line. Not JSON,
/// because the walker deliberately carries no serde (see `Cargo.toml`) and hand-parsing JSON
/// to read back something it hand-wrote is a liability for no gain. Unknown keys are ignored
/// on read, so a later version may add fields without stranding an existing cache.
///
/// Every value this holds is digits, hex, or an ISO 8601 timestamp, none of which contain
/// whitespace. A value that somehow does is dropped on write rather than written as a line
/// that would not read back the way it went in.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Sidecar {
    pub bytes: Option<u64>,
    /// Digest of the object this describes. A sidecar whose digest does not match the file
    /// beside it is describing bytes that have since been replaced, and is ignored.
    pub sha256: Option<String>,
    pub etag: Option<String>,
    pub last_modified: Option<String>,
    pub fetched_at: Option<String>,
}

impl Sidecar {
    pub const SCHEMA: u32 = 1;

    pub fn to_text(&self) -> String {
        let mut out = format!("schema {}\n", Self::SCHEMA);
        let mut put = |k: &str, v: &str| {
            if !v.is_empty() && !v.chars().any(char::is_whitespace) {
                out.push_str(k);
                out.push(' ');
                out.push_str(v);
                out.push('\n');
            }
        };
        if let Some(b) = self.bytes {
            put("bytes", &b.to_string());
        }
        for (k, v) in [
            ("sha256", &self.sha256),
            ("etag", &self.etag),
            ("last_modified", &self.last_modified),
            ("fetched_at", &self.fetched_at),
        ] {
            if let Some(v) = v {
                put(k, v);
            }
        }
        out
    }

    /// Parse, keeping whatever is legible. A sidecar is provenance, not a gate: a line that
    /// does not parse costs that one field, never the cache entry.
    pub fn from_text(text: &str) -> Self {
        let mut out = Self::default();
        for line in text.lines() {
            let Some((k, v)) = line.trim_end().split_once(' ') else {
                continue;
            };
            match k {
                "bytes" => out.bytes = v.parse().ok(),
                "sha256" => out.sha256 = Some(v.to_string()),
                "etag" => out.etag = Some(v.to_string()),
                "last_modified" => out.last_modified = Some(v.to_string()),
                "fetched_at" => out.fetched_at = Some(v.to_string()),
                _ => {}
            }
        }
        out
    }

    /// Whether this describes the bytes that were actually read back.
    pub fn describes(&self, sha256: Option<&str>) -> bool {
        match (&self.sha256, sha256) {
            (Some(a), Some(b)) => a == b,
            _ => false,
        }
    }
}

impl Payload {
    /// The provenance worth persisting beside a cached object.
    pub fn sidecar(&self) -> Sidecar {
        Sidecar {
            bytes: self.bytes,
            sha256: self.sha256.clone(),
            etag: self.etag.clone(),
            last_modified: self.last_modified.clone(),
            fetched_at: self.fetched_at.clone(),
        }
    }

    /// Adopt a sidecar's provenance onto a payload read from the cache.
    ///
    /// Only the fields a cache hit cannot otherwise know: `bytes` and `sha256` stay as they
    /// were measured from the file, because those are facts about the bytes in hand and the
    /// sidecar is only believed at all when it agrees with them.
    pub fn adopt(&mut self, s: &Sidecar) {
        self.etag = s.etag.clone();
        self.last_modified = s.last_modified.clone();
        self.fetched_at = s.fetched_at.clone();
    }
}

/// Where one accession's compressed bytes live while the pipeline works on them.
///
/// Until 0.2.0 this was always a `Vec<u8>`, which made peak memory `(C + F + W) x cap`: the
/// compressed object was resident in the queue, in the fetch worker blocked on `send`, and in
/// the compute worker for the whole parse. 20 of the 33.6 GB the walker printed as its default
/// worst case were these buffers, and none of it was necessary for a cache hit, where the
/// bytes are already a file and `zstd` can stream from the handle with about 128 KB of state.
///
/// Measured over 20 cached 34 MB accessions at `-w 1 --fetch-workers 8`, with the queue depth
/// and allocator changes already in: peak RSS 340.3 -> 49.9 MiB, and 516.4 -> 49.9 against
/// 0.1.1, a 10.4x reduction with the wall clock unchanged and every output byte-identical.
/// At `-w 8 --fetch-workers 8` the floor is the eight resident graphs: 366.9 MiB.
#[derive(Debug)]
pub enum Compressed {
    /// Held in RAM. The test fixtures, and a download by a caller that configured neither
    /// `--unitigs-cache` nor `--spool-dir`, so there was nowhere on disk to put it. No CLI run
    /// reaches that case, because `main` always resolves a spool directory.
    Memory(Vec<u8>),
    /// An open handle to a file holding exactly the compressed object: a cache hit, a download
    /// streamed into the cache, or one streamed into an unlinked file under `--spool-dir`.
    /// Decoded without being read whole.
    File(std::fs::File),
}

impl Compressed {
    /// The bytes, when they are in memory. `None` for a file-backed payload, which is the
    /// point: nothing may quietly read one of those into a `Vec`.
    pub fn as_slice(&self) -> Option<&[u8]> {
        match self {
            Compressed::Memory(v) => Some(v.as_slice()),
            Compressed::File(_) => None,
        }
    }
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
    /// Position of this accession in the input list. Carried the whole way through so the
    /// collector can restore input order; see `AccessionResult::index`.
    pub index: usize,
    pub accession: Accession,
    /// The compressed object, in RAM or as a file handle. See [`Compressed`].
    pub compressed: Compressed,
    pub compressed_len: usize,
    /// Seconds spent downloading, recorded by the fetch worker so the compute worker can
    /// report it. Network time belongs to the accession, not to whichever thread fetched it.
    pub fetch_secs: f64,
    /// Where these bytes came from and what they hash to. Carried through to
    /// `AccessionResult::payload` so it reaches `summary.json` whatever the accession goes
    /// on to do, including failing.
    pub payload: Payload,
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
/// Against ~897 MB measured for a dict-of-sets adjacency, and against 145 MB for the
/// decompressed text alone. Those are the arrays' lengths; at peak, with the slack the arrays
/// keep as they grow, a graph measured 1.24x to 1.62x its compressed file across seven
/// accessions (1.25x here). See `config::GRAPH_PEAK_PER_COMPRESSED_BYTE`. Holding sequences here rather than re-reading the file is what lets
/// `extract::collect_seqs` be O(nodes kept) instead of O(all unitigs).
pub struct UnitigGraph {
    pub accession: Accession,
    pub n_unitigs: usize,
    /// `offsets[i]..offsets[i+1]` indexes `targets`. Length `n_unitigs + 1` after
    /// `parse::scan`, `nodes.len() + 1` after `parse::collect`.
    ///
    /// **Indexed by unitig id until `collect` runs, and by slot afterwards.** `scan` needs the
    /// whole graph's adjacency because `extract::expand` walks it to find the neighbourhood;
    /// once the neighbourhood is known, nothing downstream ever asks about a unitig outside
    /// it, so `collect` rewrites these to the neighbourhood's rows and frees the rest. On a
    /// 1.28M-unitig accession with a 747-node neighbourhood that is 15.4 MB returned to the
    /// allocator while seven other compute workers are still holding theirs.
    ///
    /// Use `edges_by_id` before `collect` and `edges` after. A graph from `build` holds every
    /// unitig, so for that one the two indexings coincide.
    pub offsets: Vec<u32>,
    /// Packed edge: target id in the low 30 bits, source strand in bit 30, target strand
    /// in bit 31. Use `stages::parse::{pack_edge, unpack_edge}` rather than doing it by hand.
    /// 30 bits caps an accession at 2^30 unitigs, which `parse::build` enforces.
    pub targets: Vec<u32>,
    /// Ids whose sequence and abundance this graph holds, ascending.
    ///
    /// Empty between `parse::scan` and `parse::collect`: pass 1 builds adjacency only, which
    /// is all `stages::extract::expand` needs to find the neighbourhood. Pass 2 then decodes a
    /// second time and keeps the per-unitig data for these ids alone.
    ///
    /// Sorted, so `slot` is a binary search. The subgraph is 747 nodes for `ERR2028820`
    /// against 1,277,362 in the graph, and the only production caller that touches a sequence
    /// is `extract::collect_seqs`, once per node, so the log factor is not on any hot path.
    pub nodes: Vec<UnitigId>,
    /// Mean k-mer coverage (`ka:f:`) for `nodes[i]`, parallel to `nodes`.
    pub abundance: Vec<f32>,
    /// Set when the whole graph's adjacency was too large to keep in RAM and went to disk.
    ///
    /// Present only between `parse::scan` and `parse::collect`; `offsets` and `targets` are
    /// empty while it is. `collect` reads the neighbourhood's rows back and drops it, so no
    /// stage after extraction knows or cares whether the accession ever spilled.
    pub spool: Option<AdjacencySpool>,
    /// The sequences of `nodes`, in that order, 2 bits per base, 4 bases per byte.
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
    /// Base offset of `nodes[i]` within `seq_2bit`, in bases (not bytes). Length
    /// `nodes.len() + 1`, so the length of `nodes[i]` is `seq_offset[i+1] - seq_offset[i]`.
    ///
    /// Packing is contiguous across unitig boundaries, so a sequence may begin mid-byte.
    /// That saves the ~1.9 MB of padding byte-alignment would cost on a 1.28M-unitig
    /// accession, at the price of a shift in `decode_seq`.
    pub seq_offset: Vec<u32>,
}

impl UnitigGraph {
    /// Heap this graph holds: full adjacency plus whatever `parse::collect` retained.
    ///
    /// Capacity, not length. The slack doubling leaves behind is resident, and it is most of
    /// what the measured `GRAPH_PEAK_PER_COMPRESSED_BYTE` of 1.7 is made of. This is the
    /// quantity `--max-size-mb x 1.7` bounds, and the single definition both parse passes
    /// measure against.
    pub fn resident_bytes(&self) -> u64 {
        // A spooled graph's adjacency is on disk and costs no resident bytes, which is the
        // whole point of spilling it.
        let w = size_of::<u32>() as u64;
        self.offsets.capacity() as u64 * w
            + self.targets.capacity() as u64 * w
            + self.nodes.capacity() as u64 * w
            + self.seq_offset.capacity() as u64 * w
            + self.abundance.capacity() as u64 * size_of::<f32>() as u64
            + self.seq_2bit.capacity() as u64
    }

    /// Where unitig `id`'s per-unitig data sits, or `None` when this graph does not hold it.
    ///
    /// Absent is a normal answer, not an error: a graph holds the neighbourhood, and a caller
    /// asking about a unitig outside it is asking about something that was deliberately not
    /// kept. Callers that cannot sensibly continue use `expect_slot`, which names the id.
    #[inline]
    pub fn slot(&self, id: UnitigId) -> Option<usize> {
        self.nodes.binary_search(&id).ok()
    }

    /// Whether this graph holds unitig `id`'s sequence and abundance.
    #[inline]
    pub fn has(&self, id: UnitigId) -> bool {
        self.slot(id).is_some()
    }

    #[inline]
    fn expect_slot(&self, id: UnitigId, what: &str) -> usize {
        self.slot(id).unwrap_or_else(|| {
            panic!(
                "{what} for unitig {id}, which is not in this graph's {} retained nodes. \
                 Only the extracted neighbourhood is kept; see UnitigGraph::nodes.",
                self.nodes.len()
            )
        })
    }

    /// Mean k-mer coverage of unitig `id`.
    #[inline]
    pub fn abundance_of(&self, id: UnitigId) -> f32 {
        self.abundance[self.expect_slot(id, "abundance")]
    }

    /// Length of unitig `id` in bases.
    pub fn seq_len(&self, id: UnitigId) -> usize {
        let i = self.expect_slot(id, "sequence length");
        (self.seq_offset[i + 1] - self.seq_offset[i]) as usize
    }

    /// Decode unitig `id` into `out`, which is cleared first. Reuse one buffer across calls
    /// rather than allocating per unitig.
    pub fn decode_seq(&self, id: UnitigId, out: &mut Vec<u8>) {
        const BASES: [u8; 4] = *b"ACGT";
        let i = self.expect_slot(id, "sequence");
        let (start, end) = (self.seq_offset[i] as usize, self.seq_offset[i + 1] as usize);
        out.clear();
        out.reserve(end - start);
        for p in start..end {
            out.push(BASES[((self.seq_2bit[p >> 2] >> (2 * (p & 3))) & 3) as usize]);
        }
    }

    /// 2-bit code of base `j` of unitig `id`, without materialising the sequence.
    #[inline]
    pub fn base_at(&self, id: UnitigId, j: usize) -> u64 {
        let p = self.seq_offset[self.expect_slot(id, "base")] as usize + j;
        ((self.seq_2bit[p >> 2] >> (2 * (p & 3))) & 3) as u64
    }

    /// Neighbours of `id` as packed edges, for a graph whose neighbourhood has been collected.
    ///
    /// Includes edges that leave the neighbourhood: `extract::write_gfa` has to know which
    /// links to drop, so the row is whole even though its targets may not be retained.
    #[inline]
    pub fn edges(&self, id: UnitigId) -> &[u32] {
        let i = self.expect_slot(id, "edges");
        &self.targets[self.offsets[i] as usize..self.offsets[i + 1] as usize]
    }

    /// Neighbours of `id` into `buf`, whether the rows are indexed by id or by slot.
    ///
    /// Between `parse::scan` and `parse::collect` the graph holds a row per unitig and `id`
    /// indexes the arrays directly. After `collect` only the neighbourhood's rows survive and
    /// they are indexed by slot, so the id has to go through `nodes`. Rows are kept whole in
    /// both states, so the row this returns is the same one either way.
    ///
    /// Both states are reachable. `pipeline` expands the neighbourhood a second time when
    /// pruning is rejected, and that runs after `collect`: indexing by id there read a
    /// slot-indexed array, which is shorter than the id space by the whole of the accession.
    /// On a 5,000-accession run all 8 accessions that reached it panicked, and an id that
    /// happened to land inside the shorter array would have silently returned another
    /// unitig's edges instead.
    ///
    /// Fills a caller-owned buffer rather than returning a slice because the rows may be on
    /// disk; pass the same buffer back each time.
    pub fn edges_by_id_into(&self, id: UnitigId, buf: &mut Vec<u32>) -> std::io::Result<()> {
        match &self.spool {
            Some(s) => s.edges_into(id, buf),
            None => {
                // `nodes` is empty exactly while every row is still held, which is also the
                // state `parse::shrink_adjacency` reads the graph in while rebuilding it.
                // Once it is populated it is the id -> slot mapping.
                let i = if self.nodes.is_empty() {
                    id as usize
                } else {
                    self.expect_slot(id, "edges")
                };
                buf.clear();
                buf.extend_from_slice(
                    &self.targets[self.offsets[i] as usize..self.offsets[i + 1] as usize],
                );
                Ok(())
            }
        }
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
    /// True only when this value created `gfa_path`, which in production means
    /// `extract::write_gfa`. `Drop` unlinks the file only then.
    ///
    /// There is no default on purpose. `Subgraph` is built as a struct literal, so a new
    /// production constructor cannot omit this and silently leak tmpfs files; it will not
    /// compile until it says which it is.
    pub owns_gfa: bool,
}

/// Remove the scratch GFA when the subgraph goes out of scope, if this value created it.
///
/// It runs on every exit path, including early returns and errors.
///
/// It matters because `scratch_dir` is tmpfs, so a leaked file occupies RAM that never
/// appears in any process's RSS. At ~1 MB for a large subgraph, a 10,000-accession run
/// would quietly consume gigabytes.
///
/// The `owns_gfa` guard exists because this used to unlink `gfa_path` whoever built the value,
/// and unit tests build `Subgraph` by hand with paths they do not own: `/dev/null` and two
/// fixed `/tmp` names. As an ordinary user those unlinks fail silently, so CI never showed it.
/// As root they succeed, which turned `/dev/null` into a regular file and then failed every
/// `tests/cli.rs` case with `NotFound`, because they open it for child stdio. See TW-01.
impl Drop for Subgraph {
    fn drop(&mut self) {
        // Set LOGAN_KEEP_SCRATCH=1 to retain the GFA for debugging a specific accession.
        if self.owns_gfa && !crate::scratch::keep_scratch() {
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
/// query. `stitch` slices the assembled sequence with them, so mistaking them for query
/// coordinates yields a plausible but wrong haplotype. `query_start` and `query_end` are the
/// query coordinates.
#[derive(Clone)]
pub struct GafPath {
    pub steps: Vec<PathStep>,
    /// GAF column 3: where on the reference the alignment starts.
    pub query_start: usize,
    /// GAF column 4: where on the reference it ends (exclusive).
    pub query_end: usize,
    /// GAF column 7.
    pub path_len: usize,
    /// GAF column 8.
    pub path_start: usize,
    /// GAF column 9.
    pub path_end: usize,
    /// GAF column 10, residue matches, **less the ambiguous reference bases in the aligned
    /// window**. Also the field the best row is chosen by.
    ///
    /// GraphAligner counts `N` as matching any base and each IUPAC code as matching every
    /// base it stands for, so the raw column rewards a reference for the bases it does not
    /// pin down. `align::discount_ambiguous_matches` subtracts them once, immediately after
    /// the GAF is parsed, so every later consumer (`matched`, `pident`, the prune pass, the
    /// candidate report) sees the same corrected number.
    pub num_matches: usize,
    /// GAF column 11, alignment block length. `pident = num_matches / aln_len`.
    pub aln_len: usize,
}

impl GafPath {
    /// Percent identity of this alignment: matches over aligned length.
    ///
    /// WV-24. The formula was spelled out in `stitch` and again in `pipeline`, which is the
    /// kind of duplication that survives until one of the two is changed. `aln_len == 0`
    /// gives 0.0 rather than NaN, which would reach `summary.json` as an invalid token
    /// (WV-15).
    pub fn pident(&self) -> f64 {
        if self.aln_len > 0 {
            (self.num_matches as f64 / self.aln_len as f64) * 100.0
        } else {
            0.0
        }
    }
}

impl GafPath {
    /// Reference bases the alignment spans. Unlike the haplotype's length, an indel in the
    /// sample does not change it; an alignment that stops short does.
    pub fn query_span(&self) -> usize {
        self.query_end.saturating_sub(self.query_start)
    }
}

// ---------------------------------------------------------------------------
// stage 6 output: stitch & variants
// ---------------------------------------------------------------------------

/// A reconstructed full-length coding haplotype.
#[derive(Debug, Clone)]
pub struct Haplotype {
    pub id: String,
    pub seq: Vec<u8>,
    /// Length of this haplotype over the length of the reference.
    ///
    /// **Not what `--min-cov` gates on**, and not a fraction of the reference: an insertion
    /// in the sample pushes it above 1.0 and a deletion below, because it measures what was
    /// walked. See `ref_cov`, which is the fraction the flag's name promises and the number
    /// the gate uses. Both are published because they answer different questions: this one
    /// says how much sequence came out, `ref_cov` says how much of the reference was found.
    pub cov: f64,
    /// Fraction of the reference the alignment spans: `(query_end - query_start) / ref_len`.
    ///
    /// WA-05. This is what `--min-cov` gates on, and what the flag has always been documented
    /// as meaning. Until 0.2.0 the gate used `cov` instead, which measures the walked length,
    /// so an indel in the sample moved the decision: on a 300 nt gene at the default 0.85, a
    /// sample carrying the whole gene with a 60 nt in-frame deletion scored 0.80 and was
    /// published as `no_path`, "the gene is not here", while an alignment that reached only
    /// reference base 240 but walked a 20 nt insertion scored 0.867 and was published as a
    /// primary. The first biases presence calls against deletion alleles, which is the
    /// opposite of what a gene-presence gate should do.
    pub ref_cov: f64,
    /// The reference interval the alignment covers, GAF columns 3 and 4, end exclusive.
    ///
    /// Published so `ref_cov` can be checked and so a short alignment says *where* it was
    /// short: a walk that stops 60 bases from the 3' end is a different finding from one
    /// missing 60 bases in the middle, and the fraction alone cannot tell them apart.
    pub ref_span: (usize, usize),
    /// Identity of the GraphAligner alignment that produced this path: `num_matches /
    /// aln_len` from the GAF.
    ///
    /// When `pident_inherited` is set it describes a DIFFERENT sequence: the route as
    /// aligned, before this haplotype's path was spliced. See that field.
    pub pident: f64,
    /// True when `pident` was carried over from the alignment of another path rather than
    /// measured on this one.
    ///
    /// A spliced path (an alternative allele, or the primary after `prefer_majority` moves
    /// it onto the more abundant allele) has no GAF of its own. `variants::copy_coords`
    /// carries the original's `num_matches`/`aln_len` forward so the coordinates stay usable,
    /// but the resulting identity is the pre-splice one. Measuring it properly would need an
    /// alignment per haplotype, so it is labelled rather than recomputed: publishing it as a
    /// plain `pident=` claimed a number about the wrong sequence.
    pub pident_inherited: bool,
    pub nodes: usize,
    /// Mean k-mer coverage (`ka:f`) across the unitig path.
    ///
    /// Arithmetic mean over path steps, unweighted by unitig length. For the robust statistic
    /// used to decide prune floors, see `stages::extract::path_depth`, which is a median.
    pub mean_ka: f64,
    /// Variant allele frequency at the defining fork (1.0 if primary backbone).
    pub vaf: f64,
    /// Reference record id this haplotype was scored against. `None` on the single-anchor
    /// path (byte-identical output legacy contract); `Some(record_id)` in multi-anchor mode,
    /// where downstream tools need to know which of the pool the `cov` and `pident` numbers
    /// belong to.
    pub anchor_used: Option<String>,
}

/// A recovered second allele: the primary haplotype with exactly one allele changed.
///
/// A substitution in the common case, or an in-frame indel of at most `max_alt_indel` bases
/// when `allow_indel_alts` is on, which it is by default (WV-21).
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
    /// A branch of this fork sits at [`KA_CEILING`], so `taken_vaf` and `alt_vaf` are computed
    /// against a total that is a lower bound and both are overstated (WV-03).
    pub ka_saturated: bool,
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
/// different branches, which is exactly the distinction a strand-less adjacency loses.
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
/// `merge` is meaningful here because orientation is kept. When an adjacency drops the
/// strand fields, a bidirected edge recorded in both endpoints' headers looks like a
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
    /// A branch of this bubble is at [`KA_CEILING`], so every `vaf` on it is computed against
    /// a lower-bound total and is overstated (WV-03).
    pub ka_saturated: bool,
}

/// One reference's candidate at the point `pick_primary` chose the winner.
///
/// Multi-anchor mode fans one alignment per reference, then keeps only the winner. Every
/// non-winning candidate's `cov`/`pident`/nodes were discarded, which is what every
/// after-the-fact "which strain does this look like" analysis had to reconstruct with a
/// post-hoc BLAST against the reference panel. This carries the numbers the walker
/// already computed to the summary writer, so downstream code can re-rank or audit
/// without a second alignment.
///
/// Populated even for candidates that did not clear `min_cov`: `cov` and the alignment
/// counters come from the GAF (which `stitch` inspected before returning `None`), and a
/// user asking "why did nothing win" needs the sub-floor numbers as much as the winning
/// ones. Only single-anchor runs emit an empty vec, so their JSON stays byte-identical.
#[derive(Debug, Clone)]
pub struct CandidateReport {
    /// FASTA record id of the reference, matches [`Haplotype::anchor_used`].
    pub record_id: String,
    /// Order in the anchor pool. `pick_primary` uses this as its final tie-break, so a
    /// consumer wanting to reproduce the winner must have it.
    pub pool_index: usize,
    /// `num_matches / aln_len * 100`, same as [`Haplotype::pident`]. `None` when this
    /// reference produced no GAF row (i.e. GraphAligner found no alignment for it).
    pub pident: Option<f64>,
    /// Length of the stitched span over the reference length. `None` when there was no path
    /// to stitch. As on [`Haplotype::cov`], this is the walked length, not the fraction of
    /// the reference that was found; `ref_cov` is that.
    pub cov: Option<f64>,
    /// Fraction of the reference this candidate's alignment spans, which is what `--min-cov`
    /// gates on (WA-05). `None` when this reference produced no GAF row.
    ///
    /// Present even when `cov` is `None`, because a candidate rejected by the gate has no
    /// stitched span but does have an alignment, and this is the number that rejected it.
    pub ref_cov: Option<f64>,
    /// The reference interval the alignment covers, GAF columns 3 and 4, end exclusive.
    pub ref_span: Option<(usize, usize)>,
    /// Length of the stitched span in bases. `None` when there was no path.
    pub length: Option<usize>,
    pub nodes: Option<usize>,
    /// GAF column 7, path length in bases.
    pub path_len: Option<usize>,
    /// GAF column 10.
    pub num_matches: Option<usize>,
    /// GAF column 11.
    pub aln_len: Option<usize>,
    /// True iff this candidate was the winner `pick_primary` returned.
    pub selected: bool,
}

// ---------------------------------------------------------------------------
// per-accession result
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Outcome {
    Success,
    NoAnchors,
    NoPath,
    /// Logan has no unitigs for this accession (HTTP 404).
    NoLoganData,
    /// Skipped by `--max-size-mb`.
    TooLarge,
    /// The anchor neighbourhood is over `--max-subgraph-nodes`.
    ///
    /// Distinct from `TooLarge`, which is about the file. This is about the *graph*: the
    /// accession downloaded and parsed fine, and then the gene's neighbourhood turned out to
    /// be most of it. That is the one case the two-pass parser does not help with, because an
    /// O(subgraph) design is O(graph) when the subgraph is the graph, and it is the memory
    /// term left after `--max-size-mb` stops governing. Its own status, because the remedy is
    /// different: a bigger cap or a more specific reference, not a bigger machine.
    SubgraphTooLarge,
    Failed,
}

impl Outcome {
    /// Nothing reconstructed, but nothing went wrong either: the data held nothing to find,
    /// or the user's size cap skipped it. `--allow-empty` accepts a run of these; a
    /// `Failed` accession is an error it never hides.
    pub fn is_negative(self) -> bool {
        matches!(
            self,
            Outcome::NoAnchors
                | Outcome::NoPath
                | Outcome::NoLoganData
                | Outcome::TooLarge
                | Outcome::SubgraphTooLarge
        )
    }
}

/// Wall time per stage, in seconds.
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

/// What `stages::variants` considered and why it refused, for one accession.
///
/// Until 0.2.0 every one of these numbers existed only as an `eprintln!`, so a run that
/// recovered nothing because every candidate was refused and a run with no candidates at all
/// produced identical summaries. `docs/summaries.md` rule 3: everything printed is recorded.
#[derive(Debug, Default, Clone, Copy)]
pub struct VariantStats {
    /// Depth and floor the gates were applied at.
    pub depth: f32,
    pub floor: f32,
    /// Forks the walked route passes through, and how many had exactly one viable alternative.
    pub forks: usize,
    pub biallelic: usize,
    /// Alternatives actually emitted.
    pub recovered: usize,
    /// Credible alleles that cleared every gate but did not fit under `--max-alts`
    /// (WV-06). Non-zero means the variant outputs are incomplete for this accession.
    pub beyond_cap: usize,
    /// Of those, how many change length.
    pub indels: usize,
    // Refusals, in the order the stage applies them.
    pub no_rejoin: usize,
    pub not_substitution: usize,
    pub over_indel_cap: usize,
    pub frameshifting: usize,
    pub window_not_shiftable: usize,
    pub identical: usize,
}

/// The value Logan's `ka:f` stops at.
///
/// WV-03. Logan caps abundance, and the published values stop at exactly this number: measured
/// on three cached accessions, `ERR2028820` has 664 unitigs here against 4 at the next value
/// down, `SRR8209082` has 802 against 24, and `SRR12363134` has 1 here and the next is
/// 14,500. A distribution does not pile up on one value and empty out immediately below it;
/// that is a clip.
///
/// It matters because the walker reads `ka:f` as linear depth. Where a branch is clipped the
/// fork total is a lower bound, so the minority allele's share is OVERSTATED: a 3,000x allele
/// beside a true 250,000x major allele is 1.2%, but against a clipped 25,000 it reports 10.7%
/// and clears `--min-vaf`. The same values drive the prune floor, the bubble VAFs and the
/// majority correction. The true depth cannot be recovered from the file, so every output
/// that rests on a saturated fork is flagged rather than dropped.
pub const KA_CEILING: f32 = 25_000.0;

/// What `variants::prefer_majority` did, and what it left behind.
///
/// WV-02. Before 0.2.0 the stage refused a fork outright when the branch the aligner walked
/// was itself below `min_abund` or `min_vaf`, on the reasoning that such a branch is
/// sequencing error and belongs to `--prune-frac`. But pruning removes branches below 1% of
/// *gene depth* while this gate is 10% of the *fork*, so a walked branch between the two was
/// neither pruned, nor swapped, nor reported anywhere. The reference-like minority of a
/// contaminated or mixed sample lands exactly there, and `_primary` published it with
/// `pident=100.0` and nothing to say so.
#[derive(Debug, Default, Clone, Copy)]
pub struct MajorityReport {
    /// Sites where the route was moved onto the majority allele.
    pub swaps: usize,
    /// Of those, sites where the branch the aligner walked was itself a minority of the fork,
    /// so the swap was the WV-02 rescue rather than an ordinary two-allele correction.
    pub minority_rescues: usize,
    /// The correction loop exhausted its `max_alts * 2` pass backstop with swaps still being
    /// applied, so this route was rewritten rather than corrected at a few sites.
    pub hit_pass_cap: bool,
    /// Forks still on the final route where a correction was considered and not made,
    /// counted by reason. A non-zero total is a site where `_primary` knowingly keeps an
    /// allele the stage looked at replacing.
    ///
    /// `exempt_minority` is the WV-02 category, a fork whose walked branch is itself a
    /// minority; the others fire on ordinary two-allele sites as well. The JSON key stays
    /// `minority_left` so the 0.2.0 schema does not churn again.
    pub left: MinorityLeft,
}

/// Why a fork whose untaken branch is an outright majority was left alone. Counts only.
#[derive(Debug, Default, Clone, Copy)]
pub struct MinorityLeft {
    /// The majority branch does not clear `min_abund`: the whole fork is noise.
    pub alt_below_floor: usize,
    /// No single untaken branch holds half the fork, so there is no majority to move onto.
    /// A three-way split with the route on the smallest share lands here.
    pub not_outright_majority: usize,
    /// The branch the aligner walked is a minority of its fork and `--rescue-minority` is off,
    /// which is the default. This is the restored error-branch exemption: the site is reported
    /// but not acted on. See `Config::rescue_minority` for why.
    pub exempt_minority: usize,
    /// A branch of the fork is at [`KA_CEILING`], so the total the majority test divides by is
    /// a lower bound and the test cannot be trusted. Measured on a 5,000-accession slice of
    /// walker20k: without this guard the correction count went from 470 sites to 97,068, with
    /// a median of 64 per corrected accession, which is the loop's own backstop.
    pub saturated: usize,
    /// No path back to the route within `--alt-max-detour`.
    pub no_rejoin: usize,
    /// The length rules refused it: indels off, over `--max-alt-indel`, or a frameshift.
    pub length_refused: usize,
    /// The alignment window could not be shifted onto the spliced route.
    pub window_not_shiftable: usize,
    /// The spliced route does not clear the coverage gate.
    pub coverage_gate: usize,
    /// The detour spells the same sequence, so there is nothing to correct.
    pub identical: usize,
}

impl MinorityLeft {
    pub fn total(&self) -> usize {
        self.alt_below_floor
            + self.not_outright_majority
            + self.exempt_minority
            + self.saturated
            + self.no_rejoin
            + self.length_refused
            + self.window_not_shiftable
            + self.coverage_gate
            + self.identical
    }
}

/// Everything one accession produced. This is what crosses back to the collector.
pub struct AccessionResult {
    /// Position of this accession in the input list, 0-based.
    ///
    /// Results come back in *completion* order, which is decided by S3 and by how long each
    /// graph takes, so without this every output file was ordered differently run to run and
    /// `reconstructed.fasta` was not byte-reproducible. That also silently chose which
    /// haplotype `--sieve` keeps, since the sieve promotes the first record it sees.
    /// `pipeline::run` sorts on this before returning.
    pub index: usize,
    pub accession: Accession,
    pub outcome: Outcome,
    pub n_unitigs: usize,
    pub n_anchors: usize,
    pub subgraph_nodes: usize,
    pub haplotype: Option<Haplotype>,
    /// Everything written to the main FASTA; holds the single primary haplotype.
    pub haplotypes: Vec<Haplotype>,
    /// Second alleles recovered at biallelic sites by `stages::variants`, written to their own
    /// file (`alt_haplotypes.fasta`). Never combined with each other: one allele per record,
    /// a substitution or an in-frame indel.
    pub alts: Vec<AltHaplotype>,
    pub bubbles: Vec<Bubble>,
    pub timing: Timing,
    pub error: Option<String>,
    /// Reference the primary haplotype was scored against. Mirrors [`Haplotype::anchor_used`]
    /// so callers who look at `AccessionResult` directly (e.g. the summary writer) do not
    /// need to reach into `haplotype`.
    pub anchor_used: Option<String>,
    /// One entry per reference in the anchor pool, populated only in multi-anchor mode.
    ///
    /// Single-anchor runs leave this empty so `summary.json` stays byte-identical to the
    /// pre-multi-anchor shape.
    pub candidates: Vec<CandidateReport>,
    /// Share of the graph's unitigs in the extracted neighbourhood
    /// (`extract::subgraph_fraction`). `None` when nothing was extracted.
    pub subgraph_fraction: Option<f64>,
    /// The neighbourhood holds more unitigs per reference base than
    /// `Config::subgraph_nodes_per_base_guard`. Until 0.1.1 this was only a warning on stderr,
    /// which nothing downstream of a large run reads.
    pub provisional: bool,
    /// What the depth-based prune did. `None` when it did not run: `--prune-frac 0`, or no
    /// first-pass haplotype to measure the depth on. Until 0.1.1 this too was stderr only.
    pub prune: Option<PruneReport>,
    /// Sites where `variants::prefer_majority` swapped the route onto the majority allele.
    ///
    /// WV-07. The stage prints this, and prints a WARNING above 11, while `summary.json` said
    /// only `pident_inherited: true`, which cannot tell 1 correction from 27. A primary that
    /// was rerouted through a divergent segment looked identical to one that was not.
    pub majority_corrections: usize,
    /// What the majority correction did and what it left behind. See [`MajorityReport`].
    pub majority: MajorityReport,
    /// What the alt search considered and refused. See [`VariantStats`].
    pub variants: VariantStats,
    /// Intra-host SNV sites called for this accession, or `None` when the primary could not
    /// be mapped to the query within the alignment budget (WV-09/WV-20). `Some(0)` and `None`
    /// are different answers: no variants, against no attempt.
    pub snv_sites: Option<usize>,
    /// Largest `ka:f` anywhere in this accession's extracted neighbourhood.
    pub ka_max: f32,
    /// That maximum has reached [`KA_CEILING`], so this accession's abundances are clipped
    /// and every VAF, the prune floor and the majority calls rest on lower bounds (WV-03).
    pub depth_saturated: bool,
    /// What this accession was built from. Present whatever the outcome: a `too_large` skip
    /// still records the size that got it skipped, and a reconstruction that later failed
    /// still names the bytes it failed on.
    pub payload: Payload,
}

/// How the second, pruned alignment pass ended for one accession.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PruneResult {
    /// Realigned on the pruned subgraph, and that result replaced the first.
    Applied,
    /// Every unitig cleared the floor.
    NothingToRemove,
    /// No unitig cleared the floor, so pruning would have removed the whole subgraph.
    WouldEmpty,
    /// The realignment spanned less of the reference than the first; the first stands.
    RefusedShorterSpan,
    /// The realignment produced no haplotype; the first stands.
    RefusedNoHaplotype,
}

impl PruneResult {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Applied => "applied",
            Self::NothingToRemove => "nothing_to_remove",
            Self::WouldEmpty => "would_empty",
            Self::RefusedShorterSpan => "refused_shorter_span",
            Self::RefusedNoHaplotype => "refused_no_haplotype",
        }
    }
}

/// The prune pass for one accession, as `summary.json` reports it.
#[derive(Debug, Clone)]
pub struct PruneReport {
    pub result: PruneResult,
    /// Median `ka:f` on the first-pass route: the gene's depth.
    pub depth: f32,
    /// Unitigs below this `ka:f` are removed.
    pub floor: f32,
    /// Subgraph unitigs before the prune, and how many cleared the floor.
    pub unitigs: usize,
    pub unitigs_kept: usize,
    /// Reference span `[start, end)` of the first-pass alignment.
    pub unpruned_span: (usize, usize),
    /// Reference span of the realignment, when one produced a haplotype.
    pub pruned_span: Option<(usize, usize)>,
}

impl AccessionResult {
    /// `index` is the accession's position in the input list. It must survive every failure
    /// path, or a run where some accessions failed would come back partly ordered.
    pub fn failed(index: usize, accession: Accession, error: String) -> Self {
        Self::empty(index, accession, Outcome::Failed, Some(error))
    }

    /// A result with nothing reconstructed: a failure, or a negative found before any stage
    /// ran (no Logan data, over the size cap).
    pub fn empty(index: usize, accession: Accession, outcome: Outcome, error: Option<String>) -> Self {
        Self {
            index,
            accession,
            outcome,
            n_unitigs: 0,
            n_anchors: 0,
            subgraph_nodes: 0,
            haplotype: None,
            haplotypes: Vec::new(),
            alts: Vec::new(),
            bubbles: Vec::new(),
            timing: Timing::default(),
            error,
            anchor_used: None,
            candidates: Vec::new(),
            subgraph_fraction: None,
            provisional: false,
            prune: None,
            majority_corrections: 0,
            majority: MajorityReport::default(),
            variants: VariantStats::default(),
            snv_sites: None,
            ka_max: 0.0,
            depth_saturated: false,
            payload: Payload::default(),
        }
    }
}

#[cfg(test)]
mod tests {

    /// The sidecar format must survive a round trip, and a later version's extra keys.
    ///
    /// It is hand-written and hand-parsed, so the round trip is the only thing standing
    /// between a cache's provenance and a quiet misread. Unknown keys are ignored rather
    /// than fatal so a cache written by a later walker still reads here.
    #[test]
    fn a_sidecar_round_trips_and_tolerates_keys_it_does_not_know() {
        use super::Sidecar;
        let full = Sidecar {
            bytes: Some(34_016_719),
            sha256: Some("a".repeat(64)),
            etag: Some("d2dd288b238c3c88b46acdcd78860dce".into()),
            last_modified: Some("2024-03-31T17:38:36Z".into()),
            fetched_at: Some("2026-11-17T22:08:24Z".into()),
        };
        assert_eq!(Sidecar::from_text(&full.to_text()), full);

        // Absent fields stay absent rather than becoming empty strings.
        let bare = Sidecar { sha256: Some("b".repeat(64)), ..Default::default() };
        let back = Sidecar::from_text(&bare.to_text());
        assert_eq!(back, bare);
        assert_eq!(back.etag, None);

        // A key from a later version, and a line that is not a pair at all.
        let text = format!("{}compression zstd-19\ngarbage\n", full.to_text());
        assert_eq!(Sidecar::from_text(&text), full);

        // The schema line is written, and is not mistaken for one of the fields.
        assert!(full.to_text().starts_with("schema 1\n"));
    }

    /// A value with whitespace in it is dropped rather than written unreadably.
    ///
    /// Every real value here is digits, hex or an ISO 8601 timestamp, so this cannot happen
    /// with data from S3. It is enforced anyway because the format's one rule is that a line
    /// is a key, a space, and the rest: a value containing a space would read back short and
    /// silently wrong, which is worse than a field that is simply absent.
    #[test]
    fn a_sidecar_value_that_would_not_read_back_is_not_written() {
        use super::Sidecar;
        let bad = Sidecar {
            etag: Some("two words".into()),
            sha256: Some("c".repeat(64)),
            ..Default::default()
        };
        let text = bad.to_text();
        assert!(!text.contains("two words"), "{text}");
        let back = Sidecar::from_text(&text);
        assert_eq!(back.etag, None, "dropped, not mangled");
        assert_eq!(back.sha256, bad.sha256, "and the good fields survive");
    }

    /// A sidecar only speaks for bytes whose digest it carries.
    #[test]
    fn a_sidecar_describes_only_the_digest_it_names() {
        use super::Sidecar;
        let s = Sidecar { sha256: Some("d".repeat(64)), ..Default::default() };
        assert!(s.describes(Some(&"d".repeat(64))));
        assert!(!s.describes(Some(&"e".repeat(64))));
        assert!(!s.describes(None), "unknown bytes are not a match");
        assert!(!Sidecar::default().describes(Some(&"d".repeat(64))), "and neither is silence");
    }
    use super::*;

    /// Every outcome is deliberately either a negative or a failure.
    ///
    /// The split decides what `--allow-empty` accepts and what a scheduler is told, so a new
    /// outcome that nobody classified would quietly land on the failure side and turn "this
    /// accession had nothing to give" into "this run broke". Written as an exhaustive match so
    /// that adding a variant without deciding which side it is on does not compile.
    #[test]
    fn every_outcome_is_classified_as_found_nothing_or_failed() {
        for o in [
            Outcome::Success,
            Outcome::NoAnchors,
            Outcome::NoPath,
            Outcome::NoLoganData,
            Outcome::TooLarge,
            Outcome::SubgraphTooLarge,
            Outcome::Failed,
        ] {
            let expected = match o {
                // Nothing was found, which is a fact about the accession.
                Outcome::NoAnchors
                | Outcome::NoPath
                | Outcome::NoLoganData
                | Outcome::TooLarge
                | Outcome::SubgraphTooLarge => true,
                // Not negatives: one worked, the other broke.
                Outcome::Success | Outcome::Failed => false,
            };
            assert_eq!(o.is_negative(), expected, "{o:?} is on the wrong side of the split");
        }
    }

    /// `resident_bytes` counts capacity, not length.
    ///
    /// The slack doubling leaves behind is resident memory, and it is most of what the
    /// measured `GRAPH_PEAK_PER_COMPRESSED_BYTE` of 1.7 is made of. Measuring `len` would
    /// under-count by up to half, so a graph could reach twice the budget the startup line
    /// promised before either parse pass refused it. The budget tests cannot tell the two
    /// apart, because their budgets are small enough that either measure trips.
    #[test]
    fn resident_bytes_counts_the_slack_not_only_the_data() {
        let mut g = UnitigGraph {
            accession: "T".into(),
            n_unitigs: 1,
            offsets: Vec::new(),
            targets: Vec::new(),
            nodes: Vec::new(),
            abundance: Vec::new(),
            seq_2bit: Vec::with_capacity(4096),
            seq_offset: Vec::new(),
            spool: None,
        };
        g.seq_2bit.push(0);
        assert!(
            g.resident_bytes() >= 4096,
            "counted {} bytes for one byte held in a 4096-byte allocation",
            g.resident_bytes()
        );
    }

    /// TW-01. `Drop` unlinked `gfa_path` whoever built the value, and unit tests build
    /// `Subgraph` by hand with paths they never created, so as root the suite deleted
    /// `/dev/null` and then failed every `tests/cli.rs` case with `NotFound`.
    #[test]
    fn drop_unlinks_the_gfa_only_when_the_subgraph_owns_it() {
        let dir = std::env::temp_dir().join(format!("logan-types-drop-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let borrowed = dir.join("not-ours.gfa");
        let owned = dir.join("ours.gfa");
        let make = |path: &std::path::Path, owns_gfa: bool| Subgraph {
            accession: "ACC".into(),
            nodes: vec![],
            seqs: vec![],
            gfa_path: path.to_path_buf(),
            n_edges: 0,
            owns_gfa,
        };

        std::fs::write(&borrowed, b"H\tVN:Z:1.0\n").unwrap();
        drop(make(&borrowed, false));
        assert!(borrowed.exists(), "Drop removed a file the Subgraph did not create");

        std::fs::write(&owned, b"H\tVN:Z:1.0\n").unwrap();
        drop(make(&owned, true));
        // TW-08: `scratch::keep_scratch()` ignores the environment under `cfg(test)`, so this
        // asserts one contract rather than branching on the developer's shell. The branch
        // that used to be here was the workaround for the variable leaking into the suite.
        assert!(!owned.exists(), "Drop left a file the Subgraph did create");

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The spool directory is created lazily by the first spill, so both callers of
    /// `available_bytes` routinely ask about a path that is not there yet: the startup line
    /// prints before any accession is fetched, and the fetch gate runs before any is parsed.
    /// Answering `None` for those would disable the check on exactly the first run of a new
    /// `--unitigs-cache`, which is the run most likely to fill a disk.
    #[test]
    fn available_bytes_answers_through_a_directory_that_does_not_exist_yet() {
        let root = std::env::temp_dir();
        let here = available_bytes(&root).expect("the system temp directory must answer");
        assert!(here > 0, "a writable filesystem reported no free space");

        // Three levels of nothing, under a directory that does exist.
        let missing = root.join(format!("logan-absent-{}", std::process::id())).join("a").join("b");
        assert!(!missing.exists());
        let through = available_bytes(&missing).expect("must walk up to an existing ancestor");
        // The same filesystem, so the same answer give or take other processes writing.
        let drift = (through as f64 - here as f64).abs() / here as f64;
        assert!(drift < 0.05, "{through} against {here} is not the same filesystem");

        // A path with no existing ancestor at all cannot be answered, and must not be
        // reported as zero free: the callers read `None` as "do not check" and 0 as "refuse".
        assert_eq!(available_bytes(std::path::Path::new("")), None);
    }
}
