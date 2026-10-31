//! Everything that enters the pipeline.
//!
//! IN   accession list (one per line, `#` comments skipped) or a single accession
//!      reference CDS FASTA, in frame, single record
//! OUT  Vec<Accession>, Query

use crate::error::{Error, Result};
use crate::types::{Accession, Query};
use std::path::{Path, PathBuf};

/// Read accessions, one per line. Blank lines and `#` comments are skipped.
///
/// Duplicates are kept, not collapsed, because collapsing them would change the output
/// record set. They are counted and reported, because each duplicate costs a full download.
///
/// Note the shipped lists carry duds: `DRR000252`, `DRR028850`, `DRR028861`, `DRR006432`
/// and `DRR006471` all fail to reconstruct. Reconstruction yield and throughput are
/// different measurements and should be reported separately.
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
const ZSTD_MAGIC: [u8; 4] = [0x28, 0xb5, 0x2f, 0xfd];

/// The bytes of a user's input file (a query FASTA, an accession list), gunzipped if it is
/// gzip. Decided by the first bytes, never by a `.gz` name, the same rule as the Python
/// side's `disassembler.textio`. zstd is refused with a message rather than read as text.
/// `what` names the input in errors.
pub fn read_input(path: &Path, what: &str) -> Result<Vec<u8>> {
    let bytes = std::fs::read(path)
        .map_err(|e| Error::Usage(format!("cannot read {what} {}: {e}", path.display())))?;
    if bytes.starts_with(&ZSTD_MAGIC) {
        return Err(Error::Usage(format!(
            "{what} {} is zstd-compressed; only plain and gzip-compressed input is read. \
             Decompress it first (zstd -d)",
            path.display()
        )));
    }
    if !bytes.starts_with(&GZIP_MAGIC) {
        return Ok(bytes);
    }
    // MultiGzDecoder, so a file of several gzip members (bgzip, or `cat a.gz b.gz`) is read
    // whole rather than stopping after the first.
    let mut out = Vec::new();
    std::io::Read::read_to_end(&mut flate2::read::MultiGzDecoder::new(&bytes[..]), &mut out)
        .map_err(|e| {
            Error::Usage(format!("{what} {}: corrupt gzip data: {e}", path.display()))
        })?;
    Ok(out)
}

pub fn read_accession_list(path: &Path) -> Result<Vec<Accession>> {
    let bytes = read_input(path, "accession list")?;
    let text = String::from_utf8(bytes).map_err(|e| {
        Error::Usage(format!("cannot read accession list {}: {e}", path.display()))
    })?;

    let mut out = Vec::new();
    for (lineno, raw) in text.lines().enumerate() {
        let line = raw.trim();
        // Trimming first is forgiving and treats an indented comment as a comment.
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        validate_accession(line, &format!("{}:{}", path.display(), lineno + 1))?;
        out.push(line.to_string());
    }

    if out.is_empty() {
        return Err(Error::Usage(format!(
            "no accessions in {} (blank lines and # comments are skipped)",
            path.display()
        )));
    }

    let mut seen = std::collections::HashSet::new();
    let dupes = out.iter().filter(|a| !seen.insert((*a).clone())).count();
    if dupes > 0 {
        eprintln!(
            "[warn] {} duplicate accession(s) in {}; each is fetched again",
            dupes,
            path.display()
        );
    }
    Ok(out)
}

/// Reject anything that cannot be an accession before it becomes a confusing 404 or,
/// worse, a malformed URL.
/// Reject an accession that cannot safely be interpolated into a URL and two filesystem
/// paths, and warn about one that does not look like an SRA run.
///
/// `source` names where the value came from, e.g. `"list.txt:12"` or `"--accession"`. It used
/// to take a path and a line number, which is why only the `--accession-list` path could call
/// it: `-a` therefore accepted anything, and the value reaches `s3://logan-pub/u/<acc>/…`,
/// the scratch GFA name and `<out_dir>/<acc>.vcf`, so `-a ../../foo` escaped both directories.
pub fn validate_accession(s: &str, source: &str) -> Result<()> {
    if s.is_empty() {
        return Err(Error::Usage(format!("{source}: empty accession")));
    }
    if !s.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-') {
        return Err(Error::Usage(format!(
            "{source}: {s:?} is not a usable accession (unexpected characters)"
        )));
    }
    // Logan run accessions are SRR/ERR/DRR followed by digits. Warn rather than reject,
    // since the bucket's naming is not ours to constrain.
    let looks_right = s.len() > 3
        && matches!(&s[..3], "SRR" | "ERR" | "DRR")
        && s[3..].bytes().all(|b| b.is_ascii_digit());
    if !looks_right {
        eprintln!("[warn] {source}: {s:?} does not look like an SRA run accession");
    }
    Ok(())
}

/// Read the reference CDS: the first record, with `n_records` saying how many the file holds.
///
/// It does not warn about the others. `main` sends a file with more than one record to
/// [`read_query_set`] instead, which uses all of them (multi-anchor mode); a warning here
/// that only the first was used was printed on exactly the runs that used every record.
///
/// The sequence is uppercased and stripped of whitespace, so downstream k-mer code needs
/// no case folding.
pub fn read_query(path: &Path) -> Result<Query> {
    let bytes = read_input(path, "query")?;

    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    let mut header: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();

    for line in bytes.split(|&b| b == b'\n') {
        let line = match line.strip_suffix(b"\r") {
            Some(l) => l,
            None => line,
        };
        if line.first() == Some(&b'>') {
            if let Some(h) = header.take() {
                records.push((h, std::mem::take(&mut seq)));
            }
            let text = String::from_utf8_lossy(&line[1..]).to_string();
            header = Some(text.split_whitespace().next().unwrap_or("").to_string());
        } else if header.is_some() {
            seq.extend(line.iter().filter(|b| !b.is_ascii_whitespace()).map(|b| b.to_ascii_uppercase()));
        } else if !line.iter().all(|b| b.is_ascii_whitespace()) {
            return Err(Error::Usage(format!(
                "{}: sequence data before any '>' header; is this FASTA?",
                path.display()
            )));
        }
    }
    if let Some(h) = header.take() {
        records.push((h, seq));
    }

    if records.is_empty() {
        return Err(Error::Usage(format!(
            "{}: no FASTA records found",
            path.display()
        )));
    }
    let n_records = records.len();
    let (record_id, seq) = records.swap_remove(0);
    if seq.is_empty() {
        return Err(Error::Usage(format!(
            "{}: record {record_id:?} has no sequence",
            path.display()
        )));
    }

    // Everything downstream is codon-based, so a frame violation matters. Warn rather than
    // reject: the walker itself is frame-agnostic, and it is the codon MSA (`disassembler
    // msa`) that ultimately depends on this.
    if seq.len() % 3 != 0 {
        eprintln!(
            "[warn] {}: length {} is not a multiple of 3; codon analysis downstream expects a CDS",
            path.display(),
            seq.len()
        );
    }
    let ambiguous = seq.iter().filter(|b| !matches!(b, b'A' | b'C' | b'G' | b'T')).count();
    if ambiguous > 0 {
        eprintln!(
            "[warn] {}: {ambiguous} non-ACGT base(s); k-mers spanning them cannot match \
             Logan unitigs, which are pure ACGT",
            path.display()
        );
    }

    let name = path
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .ok_or_else(|| Error::Usage(format!("{}: cannot derive a name", path.display())))?;

    Ok(Query {
        name,
        record_id,
        seq,
        path: path.to_path_buf(),
        n_records,
    })
}

/// Load every record from every input FASTA as its own single-record [`Query`].
///
/// Used by multi-anchor mode: `read_query` keeps only the first record of a FASTA (matching
/// legacy behaviour and warning about the rest), which is what we want on today's single-`-q`
/// single-record path. `read_query_set` is what runs when the collected reference set is
/// wider than that — repeated `-q file.fa`, or a multi-record FASTA passed intentionally.
/// Each returned [`Query`] holds one record, so the "one record per Query" contract every
/// stage assumes is preserved by splitting up front here. Its `n_records` is still the count
/// of records in its file, as the field promises: it used to be 1, which told
/// `align::run` the file could be handed to GraphAligner as is, so the prune pass's
/// realignment against the winning reference aligned every record in a panel again.
///
/// Naming: the [`Query::name`] is the file stem when a file holds a single record (matching
/// `read_query`), and `<stem>__<record_id>` when a file holds several. That keeps the
/// summary filename derivation (`{name}_walker_summary.json`) unambiguous across files that
/// happen to share a record id.
pub fn read_query_set(paths: &[PathBuf]) -> Result<Vec<Query>> {
    if paths.is_empty() {
        return Err(Error::Usage("no query FASTAs provided".into()));
    }
    let mut out: Vec<Query> = Vec::new();
    for path in paths {
        let bytes = read_input(path, "query")?;

        let mut records: Vec<(String, Vec<u8>)> = Vec::new();
        let mut header: Option<String> = None;
        let mut seq: Vec<u8> = Vec::new();
        for line in bytes.split(|&b| b == b'\n') {
            let line = line.strip_suffix(b"\r").unwrap_or(line);
            if line.first() == Some(&b'>') {
                if let Some(h) = header.take() {
                    records.push((h, std::mem::take(&mut seq)));
                }
                let text = String::from_utf8_lossy(&line[1..]).to_string();
                header = Some(text.split_whitespace().next().unwrap_or("").to_string());
            } else if header.is_some() {
                seq.extend(
                    line.iter()
                        .filter(|b| !b.is_ascii_whitespace())
                        .map(|b| b.to_ascii_uppercase()),
                );
            } else if !line.iter().all(|b| b.is_ascii_whitespace()) {
                return Err(Error::Usage(format!(
                    "{}: sequence data before any '>' header; is this FASTA?",
                    path.display()
                )));
            }
        }
        if let Some(h) = header.take() {
            records.push((h, seq));
        }
        if records.is_empty() {
            return Err(Error::Usage(format!(
                "{}: no FASTA records found",
                path.display()
            )));
        }

        let stem = path
            .file_stem()
            .map(|s| s.to_string_lossy().to_string())
            .ok_or_else(|| Error::Usage(format!("{}: cannot derive a name", path.display())))?;
        let n_in_file = records.len();
        for (record_id, seq) in records {
            if seq.is_empty() {
                return Err(Error::Usage(format!(
                    "{}: record {record_id:?} has no sequence",
                    path.display()
                )));
            }
            if seq.len() % 3 != 0 {
                eprintln!(
                    "[warn] {}: record {record_id:?} length {} is not a multiple of 3; \
                     codon analysis downstream expects a CDS",
                    path.display(),
                    seq.len()
                );
            }
            let ambiguous = seq
                .iter()
                .filter(|b| !matches!(b, b'A' | b'C' | b'G' | b'T'))
                .count();
            if ambiguous > 0 {
                eprintln!(
                    "[warn] {}: record {record_id:?} has {ambiguous} non-ACGT base(s); \
                     k-mers spanning them cannot match Logan unitigs",
                    path.display()
                );
            }
            let name = if n_in_file == 1 {
                stem.clone()
            } else {
                format!("{stem}__{record_id}")
            };
            out.push(Query {
                name,
                record_id,
                seq,
                path: path.to_path_buf(),
                n_records: n_in_file,
            });
        }
    }
    if out.is_empty() {
        return Err(Error::Usage("no FASTA records found across the provided -q files".into()));
    }
    // GraphAligner reports alignments by record id, and the pool is keyed on it from there on.
    // Two references sharing an id therefore collide: the rows of both were read as one
    // reference's, which produced cov=4.000 and two candidates both marked selected.
    let mut seen: std::collections::HashMap<&str, &Path> = std::collections::HashMap::new();
    for q in &out {
        if let Some(first) = seen.insert(q.record_id.as_str(), q.path.as_path()) {
            return Err(Error::Usage(format!(
                "reference record id {:?} appears more than once ({} and {}); every -q record \
                 needs a distinct id",
                q.record_id,
                first.display(),
                q.path.display()
            )));
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn tmp(name: &str, content: &[u8]) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("logan_walker_test_{}_{}", std::process::id(), name));
        let mut f = std::fs::File::create(&p).unwrap();
        f.write_all(content).unwrap();
        p
    }

    // ---- accession list ----

    #[test]
    fn skips_comments_and_blanks_keeps_order() {
        let p = tmp("list.txt", b"# header\n\nSRR8209082\n  \nDRR000252\n   # indented comment\nDRR000253\n");
        let got = read_accession_list(&p).unwrap();
        assert_eq!(got, vec!["SRR8209082", "DRR000252", "DRR000253"]);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn trims_surrounding_whitespace() {
        let p = tmp("ws.txt", b"  SRR8209082  \n\tDRR000252\t\n");
        assert_eq!(read_accession_list(&p).unwrap(), vec!["SRR8209082", "DRR000252"]);
        std::fs::remove_file(p).ok();
    }

    /// Duplicates are kept, because collapsing them would change the
    /// output record set. They are only warned about.
    #[test]
    fn keeps_duplicates() {
        let p = tmp("dup.txt", b"SRR8209082\nSRR8209082\n");
        assert_eq!(read_accession_list(&p).unwrap().len(), 2);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn rejects_entries_that_would_break_a_url() {
        let p = tmp("bad.txt", b"SRR8209082\n../../etc/passwd\n");
        assert!(read_accession_list(&p).is_err());
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn empty_list_is_an_error_not_an_empty_run() {
        let p = tmp("empty.txt", b"# only comments\n\n");
        assert!(read_accession_list(&p).is_err());
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn missing_file_is_a_usage_error() {
        assert!(read_accession_list(Path::new("/nope/missing.txt")).is_err());
    }

    // ---- query ----

    #[test]
    fn reads_single_record_uppercases_and_joins_lines() {
        let p = tmp("q.fasta", b">pfcrt some description\nacgtACGT\nAAAcccGGG\n");
        let q = read_query(&p).unwrap();
        assert_eq!(q.record_id, "pfcrt");
        assert_eq!(q.seq, b"ACGTACGTAAACCCGGG".to_vec());
        std::fs::remove_file(p).ok();
    }

    /// The summary filename is derived from the file stem, so `name` must be the stem, not
    /// the record id.
    #[test]
    fn name_is_the_file_stem_not_the_record_id() {
        let p = tmp("query_hiv1_pol.fasta", b">pol\nACGTAC\n");
        let q = read_query(&p).unwrap();
        assert!(q.name.ends_with("query_hiv1_pol"), "got {}", q.name);
        assert_eq!(q.record_id, "pol");
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn handles_crlf_line_endings() {
        let p = tmp("crlf.fasta", b">pfcrt\r\nACGT\r\nACGT\r\n");
        assert_eq!(read_query(&p).unwrap().seq, b"ACGTACGT".to_vec());
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn multi_record_uses_first_and_warns() {
        let p = tmp("multi.fasta", b">one\nAAA\n>two\nCCC\n");
        let q = read_query(&p).unwrap();
        assert_eq!(q.record_id, "one");
        assert_eq!(q.seq, b"AAA".to_vec());
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn rejects_empty_and_headerless_and_bodyless() {
        let a = tmp("nofasta.fasta", b"ACGTACGT\n");
        assert!(read_query(&a).is_err(), "sequence before any header");
        let b = tmp("norecords.fasta", b"");
        assert!(read_query(&b).is_err(), "empty file");
        let c = tmp("nobody.fasta", b">only_header\n");
        assert!(read_query(&c).is_err(), "header with no sequence");
        for p in [a, b, c] { std::fs::remove_file(p).ok(); }
    }

    #[test]
    fn non_multiple_of_three_is_accepted_with_a_warning() {
        let p = tmp("frame.fasta", b">x\nACGTA\n");
        assert_eq!(read_query(&p).unwrap().seq.len(), 5);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn ambiguous_bases_are_kept_not_silently_dropped() {
        let p = tmp("ambig.fasta", b">x\nACGTNRYACG\n");
        let q = read_query(&p).unwrap();
        assert_eq!(q.seq, b"ACGTNRYACG".to_vec(), "must not alter the query");
        std::fs::remove_file(p).ok();
    }

    // ---- read_query_set (multi-anchor) ----

    #[test]
    fn set_from_two_files_returns_one_query_per_record() {
        let a = tmp("set_a.fasta", b">wuhan\nACGTACGT\n");
        let b = tmp("set_b.fasta", b">jn1\nAAAACCCC\n");
        let qs = read_query_set(&[a.clone(), b.clone()]).unwrap();
        assert_eq!(qs.len(), 2);
        assert_eq!(qs[0].record_id, "wuhan");
        assert_eq!(qs[1].record_id, "jn1");
        // Single-record files use the file stem as the name.
        assert!(qs[0].name.contains("set_a"), "got {}", qs[0].name);
        assert!(qs[1].name.contains("set_b"), "got {}", qs[1].name);
        // Each file holds one record.
        assert!(qs.iter().all(|q| q.n_records == 1));
        for p in [a, b] { std::fs::remove_file(p).ok(); }
    }

    #[test]
    fn set_splits_multi_record_file_into_scoped_queries() {
        let p = tmp("set_multi.fasta", b">wuhan\nAAA\n>jn1\nCCC\n>ba5\nGGG\n");
        let qs = read_query_set(std::slice::from_ref(&p)).unwrap();
        assert_eq!(qs.len(), 3);
        assert_eq!(qs[0].record_id, "wuhan");
        assert_eq!(qs[1].record_id, "jn1");
        assert_eq!(qs[2].record_id, "ba5");
        // Multi-record files scope the Query name to disambiguate identical record ids
        // across files.
        for q in &qs {
            assert!(q.name.contains("set_multi__"), "got {}", q.name);
        }
        assert!(qs[0].name.ends_with("__wuhan"), "got {}", qs[0].name);
        // n_records is the file's record count, so align never hands GraphAligner the whole
        // file for one of its records.
        assert!(qs.iter().all(|q| q.n_records == 3));
        std::fs::remove_file(p).ok();
    }

    /// Alignments come back keyed by record id, so two references sharing one collided:
    /// `cov=4.000` with both candidates marked selected. Same for one file passed twice.
    #[test]
    fn set_rejects_a_record_id_used_twice() {
        let a = tmp("dup_a.fasta", b">spike\nACGTACGTA\n");
        let b = tmp("dup_b.fasta", b">spike\nACGTACGTAACG\n");
        let err = match read_query_set(&[a.clone(), b.clone()]) {
            Ok(_) => panic!("a duplicated record id must be rejected"),
            Err(e) => e.to_string(),
        };
        assert!(err.contains("\"spike\"") && err.contains("more than once"), "got {err}");
        assert!(read_query_set(&[a.clone(), a.clone()]).is_err());
        let multi = tmp("dup_multi.fasta", b">x\nAAA\n>x\nCCC\n");
        assert!(read_query_set(std::slice::from_ref(&multi)).is_err());
        for p in [a, b, multi] { std::fs::remove_file(p).ok(); }
    }

    fn gz(data: &[u8]) -> Vec<u8> {
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        e.write_all(data).unwrap();
        e.finish().unwrap()
    }

    /// Gzip is decided by content: these files are not named `.gz`.
    #[test]
    fn gzipped_inputs_read_like_plain_ones() {
        let list = tmp("gz_list.txt", &gz(b"SRR8209082\nDRR000252\n"));
        assert_eq!(read_accession_list(&list).unwrap(), vec!["SRR8209082", "DRR000252"]);

        let q = tmp("gz_query.fasta", &gz(b">pol\nATGGCTTAA\n"));
        let query = read_query(&q).unwrap();
        assert_eq!((query.record_id.as_str(), &query.seq[..]), ("pol", &b"ATGGCTTAA"[..]));

        let set = tmp("gz_set.fasta", &gz(b">a\nATGGCTTAA\n>b\nATGAAATAA\n"));
        let ids: Vec<String> = read_query_set(std::slice::from_ref(&set))
            .unwrap()
            .into_iter()
            .map(|q| q.record_id)
            .collect();
        assert_eq!(ids, vec!["a", "b"]);
        for p in [list, q, set] { std::fs::remove_file(p).ok(); }
    }

    /// `cat a.gz b.gz`, or bgzip's blocks: every member is read, not just the first.
    #[test]
    fn every_gzip_member_is_read() {
        let mut both = gz(b"SRR8209082\n");
        both.extend(gz(b"DRR000252\n"));
        let p = tmp("gz_members.txt", &both);
        assert_eq!(read_accession_list(&p).unwrap(), vec!["SRR8209082", "DRR000252"]);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn zstd_and_corrupt_gzip_are_named() {
        let z = tmp("zst_query.fasta", &zstd::encode_all(&b">pol\nATG\n"[..], 3).unwrap());
        let err = read_query(&z).err().expect("zstd is refused").to_string();
        assert!(err.contains("zstd-compressed"), "{err}");
        let mut bad = gz(b"SRR8209082\n");
        bad.truncate(12);
        let b = tmp("bad_gz.txt", &bad);
        let err = read_accession_list(&b).expect_err("corrupt gzip").to_string();
        assert!(err.contains("corrupt gzip"), "{err}");
        for p in [z, b] { std::fs::remove_file(p).ok(); }
    }

    #[test]
    fn set_rejects_empty_input_list() {
        assert!(read_query_set(&[]).is_err());
    }

    #[test]
    fn set_preserves_input_ordering_for_tie_break() {
        // The primary picker will tie-break on this order in Phase 4, so it must be stable.
        let a = tmp("first.fasta", b">a\nACGT\n");
        let b = tmp("second.fasta", b">b\nCCCC\n");
        let c = tmp("third.fasta", b">c\nGGGG\n");
        let qs = read_query_set(&[a.clone(), b.clone(), c.clone()]).unwrap();
        assert_eq!(
            qs.iter().map(|q| q.record_id.as_str()).collect::<Vec<_>>(),
            vec!["a", "b", "c"]
        );
        for p in [a, b, c] { std::fs::remove_file(p).ok(); }
    }
}

/// Survey the repository's real reference FASTAs. Ignored by default; run with
/// `cargo test --release -- --ignored --nocapture` from the crate directory.
#[cfg(test)]
mod repo_files {
    use super::*;

    #[test]
    #[ignore]
    fn surveys_real_query_fastas() {
        let dirs = ["../queries/markers", "../queries"];
        let mut n = 0;
        for d in dirs {
            let Ok(entries) = std::fs::read_dir(d) else { continue };
            for e in entries.flatten() {
                let p = e.path();
                if p.extension().and_then(|s| s.to_str()) != Some("fasta") {
                    continue;
                }
                match read_query(&p) {
                    Ok(q) => {
                        n += 1;
                        eprintln!(
                            "  OK   {:<34} id={:<22} {} bp{}",
                            q.name,
                            q.record_id,
                            q.seq.len(),
                            if q.seq.len() % 3 == 0 { "" } else { "  <-- not in frame" }
                        );
                    }
                    Err(err) => eprintln!("  ERR  {}: {err}", p.display()),
                }
            }
        }
        assert!(n > 0, "expected to find some query fastas");
    }
}

#[cfg(test)]
mod validation_tests {
    use super::validate_accession;

    /// The value reaches `s3://logan-pub/u/<acc>/…`, the scratch GFA name and
    /// `<out_dir>/<acc>.vcf`, so anything that can leave a directory must be refused.
    #[test]
    fn path_traversal_and_separators_are_refused() {
        for bad in [
            "../../foo",
            "..",
            ".",
            "SRR8209082/../evil",
            "a/b",
            "a\\b",
            "SRR 123",
            "SRR8209082\n",
            "",
        ] {
            assert!(
                validate_accession(bad, "--accession").is_err(),
                "{bad:?} must be refused"
            );
        }
    }

    /// Accepting is the common case, and the SRA-shape check is only a warning, because the
    /// bucket's naming is not ours to constrain.
    #[test]
    fn real_accessions_are_accepted() {
        for good in ["SRR8209082", "ERR2028820", "DRR012345", "SRR12363134", "odd_but_safe-1"] {
            assert!(
                validate_accession(good, "--accession").is_ok(),
                "{good:?} must be accepted"
            );
        }
    }

    /// The message has to name where the value came from; that is why it takes a `source`
    /// string rather than a path and a line number.
    #[test]
    fn the_error_names_its_source() {
        let e = validate_accession("../x", "--accession").unwrap_err().to_string();
        assert!(e.contains("--accession"), "{e}");
        let e = validate_accession("../x", "list.txt:12").unwrap_err().to_string();
        assert!(e.contains("list.txt:12"), "{e}");
    }
}
