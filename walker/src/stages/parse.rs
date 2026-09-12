//! Step 2: single streaming pass over the unitig file.
//!
//! Replaces the first pass of `logan_path_walker.py:108-126`, which builds three unbounded
//! dicts keyed by every unitig in the accession before the anchors are consulted at all
//! (`:128`). That is the memory defect behind the 21.3 GB out-of-memory event: measured
//! ~700 bytes per unitig, against ~114 bytes on disk.
//!
//! A record looks like:
//!
//! ```text
//! >ERR2028820_2933 ka:f:5.3   L:-:443893:+  L:+:1239391:-
//! CTGTATACACCCTTCCATTTTTAAAACAACGTAAGAATTAATTAGGAAAG…
//! ```
//!
//! `L:` tokens reference the bare integer id, not the accession-prefixed one. Edges are
//! bidirected and each edge appears in *both* endpoints' headers, which is what the strand
//! bits are for. Dropping them is what makes `phase_superbubbles` mistake one edge read
//! twice for a reconvergence.
//!
//! # Why everything happens in this one pass
//!
//! The Python walker traverses the file three times (count at `:62-67`, index at `:108`,
//! sequences at `:138`). Here adjacency, abundance and 2-bit sequences are all extracted
//! while the bytes stream past, so the decompressed text is never held and the file is
//! decoded exactly once. Sequences cost 17.2 MB packed for `ERR2028820`, which is far
//! cheaper than either holding 145 MB of text or decoding a second time later.
//!
//! No regex. The Python parser is already regex-free; the cost there is the data structure
//! and the interpreter, not the tokenizer.
//!
//! IN   accession, streaming zstd decoder
//! OUT  UnitigGraph (~43 MB for 1.28M unitigs)

use crate::error::{Error, Result};
use crate::types::{Accession, UnitigGraph};
use std::io::{BufReader, Read, BufRead};

/// Target id occupies the low 30 bits of a packed edge; the two strand flags sit above it.
pub const ID_BITS: u32 = 30;
pub const ID_MASK: u32 = (1 << ID_BITS) - 1;
const SRC_FWD: u32 = 1 << 30;
const DST_FWD: u32 = 1 << 31;

/// Pack a target id with its source and target strand bits into one `u32`.
#[inline]
pub fn pack_edge(target: u32, src_fwd: bool, dst_fwd: bool) -> u32 {
    debug_assert!(target <= ID_MASK);
    target | if src_fwd { SRC_FWD } else { 0 } | if dst_fwd { DST_FWD } else { 0 }
}

/// Inverse of `pack_edge`, returning `(target, src_fwd, dst_fwd)`.
#[inline]
pub fn unpack_edge(packed: u32) -> (u32, bool, bool) {
    (
        packed & ID_MASK,
        packed & SRC_FWD != 0,
        packed & DST_FWD != 0,
    )
}

/// Consume the decoder and build the complete `UnitigGraph`.
///
/// # Structural assumptions, all enforced
///
/// Verified against real Logan output (`ERR2028820`, 1,277,362 records): ids run strictly
/// 0, 1, 2, ... with no gaps, and every sequence is on a single line. Strict ordering is
/// what makes CSR construction trivial here: a record's `L:` tokens are all the edges
/// incident to that unitig, so appending them in record order already produces the CSR
/// layout, with no edge sort and no temporary edge list.
///
/// These are assumptions about someone else's output format, so each is checked and
/// produces a clear error rather than a silently wrong graph. Multi-line sequences are
/// handled regardless, since that costs nothing.
pub fn build(accession: &Accession, reader: Box<dyn Read + '_>) -> Result<UnitigGraph> {
    let bad = |detail: String| Error::Parse {
        accession: accession.clone(),
        detail,
    };

    let mut rd = BufReader::with_capacity(1 << 18, reader);
    let mut line: Vec<u8> = Vec::with_capacity(512);

    let mut offsets: Vec<u32> = Vec::new();
    let mut targets: Vec<u32> = Vec::new();
    let mut abundance: Vec<f32> = Vec::new();
    let mut seq_2bit: Vec<u8> = Vec::new();
    let mut seq_offset: Vec<u32> = Vec::new();

    let mut expected_id: u64 = 0;
    let mut total_bases: u64 = 0;
    let mut in_record = false;

    loop {
        line.clear();
        let n = rd
            .read_until(b'\n', &mut line)
            .map_err(|e| bad(format!("read failed after {expected_id} records: {e}")))?;
        if n == 0 {
            break;
        }
        while matches!(line.last(), Some(b'\n') | Some(b'\r')) {
            line.pop();
        }
        if line.is_empty() {
            continue;
        }

        if line[0] == b'>' {
            // ---- header: id, abundance, edges ----
            let mut fields = line[1..].split(|b| b.is_ascii_whitespace()).filter(|f| !f.is_empty());
            let name = fields
                .next()
                .ok_or_else(|| bad(format!("record {expected_id}: empty header")))?;
            let id_bytes = match name.iter().rposition(|&b| b == b'_') {
                Some(p) => &name[p + 1..],
                None => name,
            };
            let id = parse_u32(id_bytes).ok_or_else(|| {
                bad(format!(
                    "record {expected_id}: cannot read a unitig id from {:?}",
                    String::from_utf8_lossy(name)
                ))
            })?;
            if id as u64 != expected_id {
                return Err(bad(format!(
                    "unitig ids are not dense and sequential: expected {expected_id}, found {id}. \
                     CSR construction relies on record order; this file needs the fallback path."
                )));
            }
            if id > ID_MASK {
                return Err(bad(format!(
                    "unitig id {id} exceeds the {ID_BITS}-bit packing limit"
                )));
            }

            offsets.push(targets.len() as u32);
            seq_offset.push(
                u32::try_from(total_bases)
                    .map_err(|_| bad("total sequence exceeds u32 addressing".into()))?,
            );

            let mut ab = 0.0f32;
            for f in fields {
                if f.starts_with(b"ka:f:") {
                    ab = std::str::from_utf8(&f[5..])
                        .ok()
                        .and_then(|s| s.parse::<f32>().ok())
                        .ok_or_else(|| {
                            bad(format!(
                                "record {id}: bad abundance {:?}",
                                String::from_utf8_lossy(f)
                            ))
                        })?;
                } else if f.starts_with(b"L:") {
                    // L:<src strand>:<target id>:<dst strand>
                    let mut parts = f[2..].split(|&b| b == b':');
                    let s1 = parts.next().unwrap_or(b"");
                    let tgt = parts.next().unwrap_or(b"");
                    let s2 = parts.next().unwrap_or(b"");
                    let (Some(src_fwd), Some(tid), Some(dst_fwd)) =
                        (strand(s1), parse_u32(tgt), strand(s2))
                    else {
                        return Err(bad(format!(
                            "record {id}: malformed link {:?}",
                            String::from_utf8_lossy(f)
                        )));
                    };
                    if tid > ID_MASK {
                        return Err(bad(format!("record {id}: link target {tid} too large")));
                    }
                    targets.push(pack_edge(tid, src_fwd, dst_fwd));
                }
            }
            abundance.push(ab);
            expected_id += 1;
            in_record = true;
        } else {
            // ---- sequence, possibly wrapped over several lines ----
            if !in_record {
                return Err(bad("sequence data before the first header".into()));
            }
            for &b in &line {
                let code = match b {
                    b'A' => 0u8,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => {
                        return Err(bad(format!(
                            "record {}: non-ACGT base {:?}; 2-bit packing would corrupt it",
                            expected_id - 1,
                            b as char
                        )))
                    }
                };
                if total_bases & 3 == 0 {
                    seq_2bit.push(0);
                }
                let last = seq_2bit.len() - 1;
                seq_2bit[last] |= code << (2 * (total_bases & 3));
                total_bases += 1;
            }
        }
    }

    if expected_id == 0 {
        return Err(bad("no unitig records found".into()));
    }
    let n_unitigs = expected_id as usize;
    offsets.push(targets.len() as u32);
    seq_offset.push(
        u32::try_from(total_bases).map_err(|_| bad("total sequence exceeds u32 addressing".into()))?,
    );

    Ok(UnitigGraph {
        accession: accession.clone(),
        n_unitigs,
        offsets,
        targets,
        abundance,
        seq_2bit,
        seq_offset,
    })
}

#[inline]
fn strand(b: &[u8]) -> Option<bool> {
    match b {
        b"+" => Some(true),
        b"-" => Some(false),
        _ => None,
    }
}

#[inline]
fn parse_u32(b: &[u8]) -> Option<u32> {
    if b.is_empty() {
        return None;
    }
    let mut v: u32 = 0;
    for &c in b {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v.checked_mul(10)?.checked_add((c - b'0') as u32)?;
    }
    Some(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_str(s: &str) -> Result<UnitigGraph> {
        build(&"TEST".to_string(), Box::new(std::io::Cursor::new(s.as_bytes().to_vec())))
    }

    #[test]
    fn edge_packing_round_trips() {
        for &(t, a, b) in &[(0u32, false, false), (1, true, false), (7, false, true), (ID_MASK, true, true)] {
            assert_eq!(unpack_edge(pack_edge(t, a, b)), (t, a, b));
        }
    }

    #[test]
    fn parses_a_minimal_graph() {
        let g = build_str(">T_0 ka:f:5.3   L:-:1:+  L:+:1:-\nACGT\n>T_1 ka:f:9.0   L:+:0:+\nTTTTAA\n").unwrap();
        assert_eq!(g.n_unitigs, 2);
        assert_eq!(g.abundance, vec![5.3, 9.0]);
        // CSR: node 0 has two edges, node 1 has one.
        assert_eq!(g.offsets, vec![0, 2, 3]);
        assert_eq!(g.edges(0).len(), 2);
        assert_eq!(g.edges(1).len(), 1);
        assert_eq!(unpack_edge(g.edges(0)[0]), (1, false, true));
        assert_eq!(unpack_edge(g.edges(0)[1]), (1, true, false));
        assert_eq!(unpack_edge(g.edges(1)[0]), (0, true, true));
    }

    #[test]
    fn sequences_round_trip_through_2bit() {
        let g = build_str(">T_0\nACGT\n>T_1\nTTTTAAAACCCCGGGGA\n").unwrap();
        let mut out = Vec::new();
        g.decode_seq(0, &mut out);
        assert_eq!(out, b"ACGT");
        g.decode_seq(1, &mut out);
        assert_eq!(out, b"TTTTAAAACCCCGGGGA");
        assert_eq!(g.seq_len(0), 4);
        assert_eq!(g.seq_len(1), 17);
    }

    /// Packing is contiguous, so unitig 1 starts mid-byte. This is the case that a
    /// byte-aligned implementation would get right by accident and this one must get right
    /// on purpose.
    #[test]
    fn unaligned_sequence_boundaries_decode_correctly() {
        let g = build_str(">T_0\nACG\n>T_1\nTTGCA\n>T_2\nG\n>T_3\nCCCCCCCC\n").unwrap();
        let mut out = Vec::new();
        for (i, want) in [(0, "ACG"), (1, "TTGCA"), (2, "G"), (3, "CCCCCCCC")] {
            g.decode_seq(i, &mut out);
            assert_eq!(out, want.as_bytes(), "unitig {i}");
        }
        assert_eq!(g.base_at(1, 2), 2, "base_at must agree with decode_seq (G=2)");
    }

    #[test]
    fn handles_multi_line_sequences() {
        let g = build_str(">T_0\nACGT\nACGT\nAC\n").unwrap();
        let mut out = Vec::new();
        g.decode_seq(0, &mut out);
        assert_eq!(out, b"ACGTACGTAC");
    }

    #[test]
    fn handles_crlf() {
        let g = build_str(">T_0 ka:f:2.0   L:+:0:+\r\nACGT\r\n").unwrap();
        assert_eq!(g.n_unitigs, 1);
        assert_eq!(g.edges(0).len(), 1);
    }

    #[test]
    fn record_without_links_or_abundance_is_fine() {
        let g = build_str(">T_0\nACGT\n").unwrap();
        assert_eq!(g.abundance, vec![0.0]);
        assert!(g.edges(0).is_empty());
    }

    /// CSR construction relies on record order. A gap must be a loud error, never a
    /// silently misaligned graph.
    #[test]
    fn non_sequential_ids_are_rejected() {
        let err = build_str(">T_0\nACGT\n>T_2\nACGT\n").err().expect("should fail");
        assert!(err.to_string().contains("not dense and sequential"), "{err}");
    }

    #[test]
    fn non_acgt_base_is_rejected_not_coerced() {
        let err = build_str(">T_0\nACGN\n").err().expect("should fail");
        assert!(err.to_string().contains("non-ACGT"), "{err}");
    }

    #[test]
    fn malformed_inputs_are_errors() {
        assert!(build_str("").is_err(), "empty input");
        assert!(build_str("ACGT\n").is_err(), "sequence before header");
        assert!(build_str(">T_0 L:x:1:+\nACGT\n").is_err(), "bad strand");
        assert!(build_str(">T_0 ka:f:notanumber\nACGT\n").is_err(), "bad abundance");
        assert!(build_str(">T_zero\nACGT\n").is_err(), "non-numeric id");
    }
}

/// Validation against a real Logan file. Ground truth was measured independently in Python:
/// `ERR2028820` has 1,277,362 unitigs, 2,582,206 `L:` tokens and 68,745,081 bp.
///
/// `LOGAN_TEST_ZST=/path/to/ERR2028820.unitigs.fa.zst cargo test --release -- --ignored --nocapture`
#[cfg(test)]
mod real_data {
    use super::*;
    use crate::stages::stream;
    use crate::types::FetchedUnitigs;

    #[test]
    #[ignore]
    fn builds_err2028820_and_matches_measured_ground_truth() {
        let Ok(path) = std::env::var("LOGAN_TEST_ZST") else {
            eprintln!("SKIP: set LOGAN_TEST_ZST");
            return;
        };
        let Ok(compressed) = std::fs::read(&path) else {
            eprintln!("SKIP: cannot read {path}");
            return;
        };
        let fetched = FetchedUnitigs {
            accession: "ERR2028820".to_string(),
            compressed_len: compressed.len(),
            compressed,
            fetch_secs: 0.0,
        };

        let t0 = std::time::Instant::now();
        let g = build(
            &fetched.accession,
            stream::decoder(&fetched).expect("decoder"),
        )
        .expect("build");
        let elapsed = t0.elapsed();

        let mem = g.offsets.len() * 4 + g.targets.len() * 4 + g.abundance.len() * 4
            + g.seq_2bit.len() + g.seq_offset.len() * 4;
        eprintln!(
            "  unitigs {:>10}  edges {:>10}  bases {:>11}  {:.2}s  graph {:.1} MB",
            g.n_unitigs,
            g.targets.len(),
            g.seq_offset[g.n_unitigs],
            elapsed.as_secs_f64(),
            mem as f64 / 1_048_576.0
        );

        assert_eq!(g.n_unitigs, 1_277_362, "unitig count");
        assert_eq!(g.targets.len(), 2_582_206, "edge count");
        assert_eq!(g.seq_offset[g.n_unitigs], 68_745_081, "total bases");
        assert_eq!(g.offsets.len(), g.n_unitigs + 1);
        assert_eq!(g.abundance.len(), g.n_unitigs);

        // Known record: >ERR2028820_2933 ka:f:5.3   L:-:443893:+  L:+:1239391:-
        assert!((g.abundance[2933] - 5.3).abs() < 1e-6, "abundance {}", g.abundance[2933]);
        let e: Vec<_> = g.edges(2933).iter().map(|&p| unpack_edge(p)).collect();
        assert_eq!(e, vec![(443893, false, true), (1239391, true, false)]);

        // Sequence content, decoded back out of the 2-bit store.
        let mut out = Vec::new();
        g.decode_seq(2933, &mut out);
        assert!(
            out.starts_with(b"CTGTATACACCCTTCCATTTTTAAAACAACGTAAGAATTAATTAGGAAAG"),
            "unitig 2933 decoded as {:?}",
            String::from_utf8_lossy(&out[..50.min(out.len())])
        );

        // Every unitig is at least k long and the largest matches the measured maximum.
        let mut min_len = usize::MAX;
        let mut max_len = 0usize;
        for i in 0..g.n_unitigs as u32 {
            let l = g.seq_len(i);
            min_len = min_len.min(l);
            max_len = max_len.max(l);
        }
        eprintln!("  unitig length min {min_len} max {max_len}");
        assert_eq!(min_len, 31, "shortest unitig should be exactly k");
        assert_eq!(max_len, 1588, "longest unitig");
    }
}
