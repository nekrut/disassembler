//! Everything that leaves the pipeline.
//!
//! IN   Summary
//! OUT  reconstructed.fasta, bubbles.json, <query>_walker_summary.json
//!
//! `reconstructed.fasta` is consumed by `build_walker_codon_msa.py`, which feeds MAFFT and
//! FastTree and then HyphAeon, so its record format must match what the Python walker emits
//! byte for byte (`logan_path_walker.py:380`).
//!
//! JSON is written by hand rather than through serde. The shapes are small and fixed, and
//! matching Python's `json.dump(..., indent=2)` layout exactly matters more here than
//! convenience, because the two summaries are meant to be diffable.

use crate::error::Result;
use crate::pipeline::Summary;
use crate::types::{AccessionResult, Outcome, Query, Strand};
use std::io::{BufWriter, Write};
use std::path::Path;

fn outcome_str(o: Outcome) -> &'static str {
    match o {
        Outcome::Success => "success",
        Outcome::NoAnchors => "no_anchors",
        Outcome::NoPath => "no_path",
        Outcome::Failed => "error",
    }
}

fn strand_str(s: Strand) -> &'static str {
    match s {
        Strand::Forward => "+",
        Strand::Reverse => "-",
    }
}

/// Minimal JSON string escaping. Accessions and ids are alphanumeric in practice, but an
/// error message can contain anything.
fn esc(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 8);
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

/// One record per reconstructed haplotype, in the order accessions completed.
///
/// Header format is fixed by `logan_path_walker.py:380`:
/// `>{id} cov={cov:.3f} pident={pident:.1f} nodes={n}`, sequence on a single line.
pub fn write_fasta(summary: &Summary, path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    for r in &summary.results {
        for h in &r.haplotypes {
            writeln!(
                w,
                ">{} cov={:.3} pident={:.1} nodes={} ka={:.1} vaf={:.2}",
                h.id, h.cov, h.pident, h.nodes, h.mean_ka, h.vaf
            )?;
            w.write_all(&h.seq)?;
            w.write_all(b"\n")?;
        }
    }
    w.flush()?;
    Ok(())
}

/// Recovered second alleles, one record per biallelic site.
///
/// Deliberately a **separate file** from `reconstructed.fasta`. `build_walker_codon_msa.py`
/// feeds that one straight into MAFFT, FastTree and HyPhy, and every sequence there becomes
/// an independent taxon. Two records from one sample are not two independent observations,
/// and HyPhy has no way to express "this one is 31% of that sample", so folding these in
/// would inflate apparent diversity rather than describe it. Keeping them apart means the
/// variance is recovered and auditable, and using it in an analysis stays a deliberate act.
///
/// The header carries the frequencies, so a reader can tell a 50/50 site from a 90/10 one:
///
/// ```text
/// >ACC_alt1 cov=0.999 pident=98.1 nodes=1281 vaf=0.310 alt_kaf=3600.0 primary_vaf=0.688 fork=49624+
/// ```
pub fn write_alts(summary: &Summary, path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    for r in &summary.results {
        for a in &r.alts {
            let h = &a.haplotype;
            writeln!(
                w,
                ">{} cov={:.3} pident={:.1} nodes={} vaf={:.3} alt_kaf={:.1} \
                 primary_vaf={:.3} fork={}{} alt_node={}{} detour={} dlen={:+} span={}",
                h.id, h.cov, h.pident, h.nodes,
                a.alt_vaf, a.alt_abundance, a.taken_vaf,
                a.fork, strand_str(a.fork_strand),
                a.alt_node, strand_str(a.alt_strand),
                a.detour_nodes,
                a.length_delta,
                a.replaced_span
            )?;
            w.write_all(&h.seq)?;
            w.write_all(b"\n")?;
        }
    }
    w.flush()?;
    Ok(())
}

/// Full per-bubble detail: fork, merge, and every branch with orientation, coverage and VAF.
///
/// This file has no Python counterpart. There, `phase_superbubbles` computes all of this and
/// then `logan_path_walker.py:317` keeps only `len(bubbles)`, so a reported quasispecies
/// figure cannot be checked by anyone. See `work/issue_intra_accession_variance_lost.md`.
pub fn write_bubbles(summary: &Summary, path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    writeln!(w, "[")?;
    for (i, r) in summary.results.iter().enumerate() {
        writeln!(w, "  {{")?;
        writeln!(w, "    \"accession\": \"{}\",", esc(&r.accession))?;
        writeln!(w, "    \"n_bubbles\": {},", r.bubbles.len())?;
        writeln!(w, "    \"bubbles\": [")?;
        for (j, b) in r.bubbles.iter().enumerate() {
            writeln!(w, "      {{")?;
            writeln!(w, "        \"fork\": {}, \"fork_strand\": \"{}\",", b.fork, strand_str(b.fork_strand))?;
            writeln!(w, "        \"merge\": {}, \"merge_strand\": \"{}\",", b.merge, strand_str(b.merge_strand))?;
            writeln!(w, "        \"branches\": [")?;
            for (k, br) in b.branches.iter().enumerate() {
                writeln!(
                    w,
                    "          {{\"node\": {}, \"strand\": \"{}\", \"ka\": {:.1}, \"vaf\": {:.4}}}{}",
                    br.node,
                    strand_str(br.strand),
                    br.abundance,
                    br.vaf,
                    if k + 1 == b.branches.len() { "" } else { "," }
                )?;
            }
            writeln!(w, "        ]")?;
            writeln!(w, "      }}{}", if j + 1 == r.bubbles.len() { "" } else { "," })?;
        }
        writeln!(w, "    ]")?;
        writeln!(w, "  }}{}", if i + 1 == summary.results.len() { "" } else { "," })?;
    }
    writeln!(w, "]")?;
    w.flush()?;
    Ok(())
}

fn write_result(w: &mut impl Write, r: &AccessionResult, last: bool) -> Result<()> {
    writeln!(w, "    {{")?;
    writeln!(w, "      \"acc\": \"{}\",", esc(&r.accession))?;
    writeln!(w, "      \"status\": \"{}\",", outcome_str(r.outcome))?;
    writeln!(w, "      \"n_unitigs\": {},", r.n_unitigs)?;
    writeln!(w, "      \"n_anchors\": {},", r.n_anchors)?;
    writeln!(w, "      \"subgraph_nodes\": {},", r.subgraph_nodes)?;
    if r.haplotypes.is_empty() {
        writeln!(w, "      \"reconstructed\": [],")?;
    } else {
        writeln!(w, "      \"reconstructed\": [")?;
        for (h_idx, h) in r.haplotypes.iter().enumerate() {
            let last_h = h_idx + 1 == r.haplotypes.len();
            writeln!(w, "        {{")?;
            writeln!(w, "          \"haplotype_id\": \"{}\",", esc(&h.id))?;
            writeln!(w, "          \"seq\": \"{}\",", String::from_utf8_lossy(&h.seq))?;
            writeln!(w, "          \"cov\": {},", h.cov)?;
            writeln!(w, "          \"pident\": {},", h.pident)?;
            writeln!(w, "          \"length\": {},", h.seq.len())?;
            writeln!(w, "          \"nodes\": {},", h.nodes)?;
            writeln!(w, "          \"mean_ka\": {:.1},", h.mean_ka)?;
            writeln!(w, "          \"vaf\": {:.2}", h.vaf)?;
            writeln!(w, "        }}{}", if last_h { "" } else { "," })?;
        }
        writeln!(w, "      ],")?;
    }
    writeln!(w, "      \"bubbles\": {},", r.bubbles.len())?;
    if let Some(e) = &r.error {
        writeln!(w, "      \"error\": \"{}\",", esc(e))?;
    }
    writeln!(w, "      \"timing\": {{")?;
    let t = &r.timing;
    writeln!(w, "        \"stream\": {:.3},", t.stream)?;
    writeln!(w, "        \"parse\": {:.3},", t.parse)?;
    writeln!(w, "        \"anchor\": {:.3},", t.anchor)?;
    writeln!(w, "        \"extract\": {:.3},", t.extract)?;
    writeln!(w, "        \"align\": {:.3},", t.align)?;
    writeln!(w, "        \"stitch\": {:.3},", t.stitch)?;
    writeln!(w, "        \"phase\": {:.3},", t.phase)?;
    writeln!(w, "        \"total\": {:.3}", t.total)?;
    writeln!(w, "      }}")?;
    writeln!(w, "    }}{}", if last { "" } else { "," })?;
    Ok(())
}

/// Per-accession counts and stage timings.
///
/// Top-level shape matches `logan_path_walker.py:412`, as do the per-result field names, so
/// the two files can be diffed directly. The `timing` block deliberately does not: Python
/// reports `stream`/`blast`/`subgraph`/`walk`, which do not map one to one onto the seven
/// stages here. `blast` corresponds to `anchor`, and Python's `subgraph` spans what is split
/// here into `parse` and `extract`.
pub fn write_summary(summary: &Summary, path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    writeln!(w, "{{")?;
    writeln!(w, "  \"total\": {},", summary.total)?;
    writeln!(w, "  \"success\": {},", summary.success)?;
    writeln!(w, "  \"results\": [")?;
    for (i, r) in summary.results.iter().enumerate() {
        write_result(&mut w, r, i + 1 == summary.results.len())?;
    }
    writeln!(w, "  ],")?;
    writeln!(w, "  \"bubble_events\": [")?;
    let with_bubbles: Vec<&AccessionResult> =
        summary.results.iter().filter(|r| !r.bubbles.is_empty()).collect();
    for (i, r) in with_bubbles.iter().enumerate() {
        writeln!(
            w,
            "    {{\"acc\": \"{}\", \"bubbles\": {}}}{}",
            esc(&r.accession),
            r.bubbles.len(),
            if i + 1 < with_bubbles.len() { "," } else { "" }
        )?;
    }
    writeln!(w, "  ]")?;
    writeln!(w, "}}")?;
    w.flush()?;
    Ok(())
}

// ---------------------------------------------------------------------------
// VCF: intra-host SNPs with frequency, in query coordinates
// ---------------------------------------------------------------------------

/// Write one VCF per accession listing its intra-host single-nucleotide polymorphisms.
///
/// The sites are the on-path biallelic forks recovered by `stages::variants` (the same set
/// that produces `alt_haplotypes.fasta`), NOT the `phase` stage's `bubbles.json`, most of
/// which sit off the walked path and are not alleles of the gene.
///
/// Coordinates are relative to the query gene, 1-based: `CHROM` is the query record id and
/// `POS` a position within the query CDS, so VCFs for the same gene across accessions share a
/// coordinate frame and `bcftools merge` yields a sites x samples frequency matrix.
///
/// Only substitutions are emitted here; indel alleles are left to `alt_haplotypes.fasta`.
/// Frequencies (`AF`) are ratios of k-mer coverage (`ka:f`) at the fork, a proxy for allele
/// frequency, not read counts, and they are per-site, not phased. Both caveats are stated in
/// the file header.
///
/// One file `<out_dir>/<accession>.vcf` per accession that has at least one recovered SNP.
pub fn write_vcfs(summary: &Summary, query: &Query, out_dir: &Path) -> Result<()> {
    std::fs::create_dir_all(out_dir)?;
    let chrom = chrom_of(query);

    for r in &summary.results {
        let rows = snv_sites(r, query);
        if rows.is_empty() {
            continue;
        }

        let path = out_dir.join(format!("{}.vcf", r.accession));
        let mut w = BufWriter::new(std::fs::File::create(&path)?);
        write_vcf_header(&mut w, chrom, query.seq.len(), &r.accession)?;
        for row in &rows {
            let alt_bases: Vec<String> =
                row.alts.iter().map(|&(b, _, _)| (b as char).to_string()).collect();
            let afs: Vec<String> = row.alts.iter().map(|&(_, af, _)| format!("{af:.4}")).collect();
            let aos: Vec<String> = row.alts.iter().map(|&(_, _, ka)| format!("{ka:.1}")).collect();
            writeln!(
                w,
                "{}\t{}\t.\t{}\t{}\t.\tPASS\tAF={};AO={};RO={:.1};FORK={}\tGT:AF:AO:RO\t.:{}:{}:{:.1}",
                chrom,
                row.pos,
                row.ref_base as char,
                alt_bases.join(","),
                afs.join(","),
                aos.join(","),
                row.ref_ka,
                row.fork,
                afs.join(","),
                aos.join(","),
                row.ref_ka,
            )?;
        }
        w.flush()?;
    }
    Ok(())
}

/// One append-only long ("tidy") table of every intra-host SNP across the whole run, one row
/// per (accession, site, ALT allele).
///
/// This is the cohort output that scales. A merged multi-sample VCF carries one field per
/// sample on every line, so a gene matched by millions of accessions would need millions of
/// columns of an almost-entirely-empty matrix. The long form is sparse-native: total rows are
/// the sum of SNPs across accessions, a few per accession, and it concatenates trivially across
/// cluster shards (each shard writes its own file; `cat` them, or point DuckDB at the glob).
/// Per-site summaries and cross-sample frequency structure are then `GROUP BY` / self-join
/// queries over this table, never a dense matrix.
///
/// Columns: `gene pos ref alt accession af ka_alt ka_ref fork`. Frequencies are `ka:f` ratios,
/// the same values the VCF carries.
pub fn write_snv_table(summary: &Summary, query: &Query, path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let chrom = chrom_of(query);
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    writeln!(w, "gene\tpos\tref\talt\taccession\taf\tka_alt\tka_ref\tfork")?;
    for r in &summary.results {
        for row in snv_sites(r, query) {
            for &(alt, af, ka_alt) in &row.alts {
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.1}\t{:.1}\t{}",
                    chrom,
                    row.pos,
                    row.ref_base as char,
                    alt as char,
                    r.accession,
                    af,
                    ka_alt,
                    row.ref_ka,
                    row.fork,
                )?;
            }
        }
    }
    w.flush()?;
    Ok(())
}

fn chrom_of(query: &Query) -> &str {
    if query.record_id.is_empty() {
        query.name.as_str()
    } else {
        query.record_id.as_str()
    }
}

/// The intra-host SNP sites for one accession, in query coordinates, sorted by position.
///
/// Shared by the per-accession VCF and the long table so the two cannot disagree. Substitutions
/// only; an alt whose length differs from the primary, or that differs at more than one base, is
/// an indel or a multi-substitution block and stays in `alt_haplotypes.fasta`.
fn snv_sites(r: &AccessionResult, query: &Query) -> Vec<VcfRow> {
    let Some(primary) = &r.haplotype else { return Vec::new() };
    if r.alts.is_empty() {
        return Vec::new();
    }
    // Guard: full DP is O(n*m); genes are short. Skip pathological inputs loudly.
    if primary.seq.len() > 20_000 || query.seq.len() > 20_000 {
        eprintln!(
            "[{}] SNP sites skipped: sequence longer than 20 kb (primary {}, query {})",
            r.accession,
            primary.seq.len(),
            query.seq.len()
        );
        return Vec::new();
    }

    // Map each primary base to a query position (None where the primary base aligns to a gap in
    // the query, i.e. an insertion relative to the gene).
    let map = align_primary_to_query(&primary.seq, &query.seq);

    let mut rows: Vec<VcfRow> = Vec::new();
    for a in &r.alts {
        let alt_seq = &a.haplotype.seq;
        if alt_seq.len() != primary.seq.len() {
            continue;
        }
        let mut diff: Option<usize> = None;
        let mut multiple = false;
        for (i, (x, y)) in primary.seq.iter().zip(alt_seq.iter()).enumerate() {
            if x != y {
                if diff.is_some() {
                    multiple = true;
                    break;
                }
                diff = Some(i);
            }
        }
        let (Some(i), false) = (diff, multiple) else { continue };
        let Some(qpos) = map[i] else { continue }; // primary base has no query coordinate

        let ref_base = query.seq[qpos];
        let major = primary.seq[i]; // the primary carries the majority allele
        let minor = alt_seq[i];
        let major_ka = a.taken_abundance;
        let minor_ka = a.alt_abundance;

        // Report every sample allele that differs from the query base as an ALT. At most one of
        // {major, minor} can equal REF (they differ from each other), so this is one ALT in the
        // common case and two only when the gene base is a third allele.
        let mut alts: Vec<(u8, f64, f32)> = Vec::new();
        if minor != ref_base {
            alts.push((minor, a.alt_vaf, minor_ka));
        }
        if major != ref_base {
            alts.push((major, a.taken_vaf, major_ka));
        }
        if alts.is_empty() {
            continue; // neither allele differs from the gene base: nothing to call
        }
        let ref_ka = if ref_base == major {
            major_ka
        } else if ref_base == minor {
            minor_ka
        } else {
            0.0
        };
        rows.push(VcfRow { pos: qpos + 1, ref_base, alts, ref_ka, fork: a.fork });
    }
    rows.sort_by_key(|r| (r.pos, r.fork));
    rows
}

struct VcfRow {
    pos: usize,
    ref_base: u8,
    /// (allele base, frequency, ka:f) for each sample allele that differs from the gene base.
    alts: Vec<(u8, f64, f32)>,
    ref_ka: f32,
    fork: crate::types::UnitigId,
}

fn write_vcf_header<W: Write>(w: &mut W, chrom: &str, len: usize, sample: &str) -> Result<()> {
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(
        w,
        "##source=logan-walker intra-host SNPs. Frequencies are k-mer coverage (ka:f) ratios \
         at a graph fork, a proxy for allele frequency, not read counts. Sites are per-locus \
         and NOT phased: two rows do not assert co-occurrence on one molecule."
    )?;
    writeln!(w, "##contig=<ID={chrom},length={len}>")?;
    writeln!(w, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency, ka:f ratio at the fork\">")?;
    writeln!(w, "##INFO=<ID=AO,Number=A,Type=Float,Description=\"Alternate allele k-mer coverage (ka:f)\">")?;
    writeln!(w, "##INFO=<ID=RO,Number=1,Type=Float,Description=\"Reference allele k-mer coverage (ka:f); 0 if the gene base is not observed\">")?;
    writeln!(w, "##INFO=<ID=FORK,Number=1,Type=Integer,Description=\"Unitig id of the fork on the walked path\">")?;
    writeln!(w, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype (missing: this is a frequency call, not a diploid genotype)\">")?;
    writeln!(w, "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">")?;
    writeln!(w, "##FORMAT=<ID=AO,Number=A,Type=Float,Description=\"Alternate allele k-mer coverage\">")?;
    writeln!(w, "##FORMAT=<ID=RO,Number=1,Type=Float,Description=\"Reference allele k-mer coverage\">")?;
    writeln!(w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}")?;
    Ok(())
}

/// Global (Needleman-Wunsch) alignment of the primary against the query, returning for each
/// primary base the query position it aligns to, or `None` if it aligns to a gap (an
/// insertion relative to the gene).
///
/// Linear gap penalty; match +1, mismatch -1, gap -2. The two sequences are the same gene at
/// high identity, so the alignment is near-diagonal and the exact scoring barely matters; what
/// matters is that a rare indel shifts coordinates correctly rather than being ignored.
fn align_primary_to_query(primary: &[u8], query: &[u8]) -> Vec<Option<usize>> {
    let n = primary.len();
    let m = query.len();
    // Direction matrix: 0 diag, 1 up (gap in query / primary insertion), 2 left (gap in primary).
    let mut dir = vec![0u8; (n + 1) * (m + 1)];
    let mut prev = vec![0i32; m + 1];
    let mut curr = vec![0i32; m + 1];
    const GAP: i32 = -2;
    for j in 0..=m {
        prev[j] = GAP * j as i32;
        dir[j] = 2; // top row: all left
    }
    for i in 1..=n {
        curr[0] = GAP * i as i32;
        dir[i * (m + 1)] = 1; // first column: all up
        for j in 1..=m {
            let s = if primary[i - 1] == query[j - 1] { 1 } else { -1 };
            let diag = prev[j - 1] + s;
            let up = prev[j] + GAP;
            let left = curr[j - 1] + GAP;
            let (mut best, mut d) = (diag, 0u8);
            if up > best {
                best = up;
                d = 1;
            }
            if left > best {
                best = left;
                d = 2;
            }
            curr[j] = best;
            dir[i * (m + 1) + j] = d;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    // Traceback from (n, m).
    let mut map = vec![None; n];
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        let d = dir[i * (m + 1) + j];
        if i > 0 && j > 0 && d == 0 {
            map[i - 1] = Some(j - 1);
            i -= 1;
            j -= 1;
        } else if i > 0 && d == 1 {
            // primary base i-1 aligns to a gap in the query
            map[i - 1] = None;
            i -= 1;
        } else {
            j -= 1;
        }
    }
    map
}

#[cfg(test)]
mod vcf_tests {
    use super::align_primary_to_query;

    #[test]
    fn identical_sequences_map_one_to_one() {
        let s = b"ACGTACGTAC";
        let map = align_primary_to_query(s, s);
        assert_eq!(map, (0..s.len()).map(Some).collect::<Vec<_>>());
    }

    #[test]
    fn a_substitution_keeps_coordinates_aligned() {
        let map = align_primary_to_query(b"ACGTACGT", b"ACGTTCGT");
        assert_eq!(map, (0..8).map(Some).collect::<Vec<_>>());
    }

    #[test]
    fn an_insertion_in_the_primary_has_no_query_coordinate() {
        // primary "ACGTGACGT" has one extra G at index 4 relative to query "ACGTACGT"
        let map = align_primary_to_query(b"ACGTGACGT", b"ACGTACGT");
        assert_eq!(map.len(), 9);
        assert_eq!(map[0], Some(0));
        assert_eq!(map[4], None, "the inserted base maps to no query position");
        assert_eq!(map[8], Some(7), "bases after the insertion still map correctly");
    }

    #[test]
    fn a_deletion_in_the_primary_skips_a_query_position() {
        // primary "ACGACGT" is query "ACGTACGT" missing the T at query index 3
        let map = align_primary_to_query(b"ACGACGT", b"ACGTACGT");
        assert_eq!(map.len(), 7);
        assert_eq!(map[0], Some(0));
        assert_eq!(map[2], Some(2));
        assert_eq!(map[3], Some(4), "query position 3 (the deleted T) is skipped");
        assert_eq!(map[6], Some(7));
    }
}
