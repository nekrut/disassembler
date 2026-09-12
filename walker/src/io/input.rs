//! Everything that enters the pipeline.
//!
//! IN   accession list (one per line, `#` comments skipped) or a single accession
//!      reference CDS FASTA, in frame, single record
//! OUT  Vec<Accession>, Query

use crate::error::{Error, Result};
use crate::types::{Accession, Query};
use std::path::Path;

/// Read accessions, one per line. Blank lines and `#` comments are skipped.
///
/// Duplicates are kept, not collapsed, to match `logan_path_walker.py:352-354`. They are
/// counted and reported, because each duplicate costs a full download.
///
/// Note the shipped lists carry duds: `DRR000252`, `DRR028850`, `DRR028861`, `DRR006432`
/// and `DRR006471` all fail to reconstruct. Reconstruction yield and throughput are
/// different measurements and should be reported separately.
pub fn read_accession_list(path: &Path) -> Result<Vec<Accession>> {
    let text = std::fs::read_to_string(path).map_err(|e| {
        Error::Usage(format!("cannot read accession list {}: {e}", path.display()))
    })?;

    let mut out = Vec::new();
    for (lineno, raw) in text.lines().enumerate() {
        let line = raw.trim();
        // Python tests `startswith("#")` on the untrimmed line. Trimming first is more
        // forgiving and treats an indented comment as a comment.
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        validate_accession(line, lineno + 1, path)?;
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
fn validate_accession(s: &str, lineno: usize, path: &Path) -> Result<()> {
    if !s.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-') {
        return Err(Error::Usage(format!(
            "{}:{}: {s:?} is not a usable accession (unexpected characters)",
            path.display(),
            lineno
        )));
    }
    // Logan run accessions are SRR/ERR/DRR followed by digits. Warn rather than reject,
    // since the bucket's naming is not ours to constrain.
    let looks_right = s.len() > 3
        && matches!(&s[..3], "SRR" | "ERR" | "DRR")
        && s[3..].bytes().all(|b| b.is_ascii_digit());
    if !looks_right {
        eprintln!(
            "[warn] {}:{}: {s:?} does not look like an SRA run accession",
            path.display(),
            lineno
        );
    }
    Ok(())
}

/// Read the reference CDS.
///
/// Takes the first record, matching `logan_path_walker.py:183`
/// (`next(SeqIO.parse(query_fa, "fasta"))`), but unlike Python it says so out loud when
/// there are more. A multi-record query is a real hazard here: the Galaxy kmindex and
/// LexicMap wrappers silently use only the first sequence too (tools-iuc #8388), so a
/// multi-FASTA query produces results for one gene while looking like it covered several.
///
/// The sequence is uppercased and stripped of whitespace, so downstream k-mer code needs
/// no case folding.
pub fn read_query(path: &Path) -> Result<Query> {
    let bytes = std::fs::read(path)
        .map_err(|e| Error::Usage(format!("cannot read query {}: {e}", path.display())))?;

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
    if records.len() > 1 {
        eprintln!(
            "[warn] {}: {} records present, using only the first ({:?}). \
             One gene per query: see tools-iuc #8388.",
            path.display(),
            records.len(),
            records[0].0
        );
    }

    let (record_id, seq) = records.swap_remove(0);
    if seq.is_empty() {
        return Err(Error::Usage(format!(
            "{}: record {record_id:?} has no sequence",
            path.display()
        )));
    }

    // Everything downstream is codon-based, so a frame violation matters. Warn rather than
    // reject: the walker itself is frame-agnostic, and it is build_walker_codon_msa.py that
    // ultimately depends on this.
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
    })
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

    /// Duplicates are kept, matching Python, because collapsing them would change the
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

    /// The summary filename is derived from the file stem in Python
    /// (`logan_path_walker.py:410`), so `name` must be the stem, not the record id.
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
