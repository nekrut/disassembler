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
///
/// The identity key is `pident=` when it was measured on this haplotype's own alignment and
/// `inherited_pident=` when it came from the pre-splice route (see `Haplotype::pident`), so a
/// consumer reading `pident=` never gets a number about a different sequence.
pub fn write_fasta(summary: &Summary, path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    for r in &summary.results {
        for h in &r.haplotypes {
            // anchor_used is emitted only when set (multi-anchor mode). Single-anchor
            // runs keep the header byte-identical to the pre-multi-anchor shape so all
            // golden fixtures diff clean.
            let anchor_suffix = match &h.anchor_used {
                Some(rid) => format!(" anchor={rid}"),
                None => String::new(),
            };
            writeln!(
                w,
                ">{} cov={:.3} {}={:.1} nodes={} ka={:.1} vaf={:.2}{}",
                h.id,
                h.cov,
                if h.pident_inherited { "inherited_pident" } else { "pident" },
                h.pident,
                h.nodes,
                h.mean_ka,
                h.vaf,
                anchor_suffix,
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
/// >ACC_alt1 cov=0.999 inherited_pident=98.1 nodes=1281 vaf=0.310 alt_kaf=3600.0 primary_vaf=0.688 fork=49624+
///
/// `inherited_pident` is the PRIMARY's alignment identity, not this allele's: an alt is a
/// spliced path with no alignment of its own. It was published as `pident=` per record, which
/// read as a measurement of that allele.
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
            let anchor_suffix = match &h.anchor_used {
                Some(rid) => format!(" anchor={rid}"),
                None => String::new(),
            };
            writeln!(
                w,
                ">{} cov={:.3} inherited_pident={:.1} nodes={} vaf={:.3} alt_kaf={:.1} \
                 primary_vaf={:.3} fork={}{} alt_node={}{} detour={} dlen={:+} span={}{}",
                h.id, h.cov, h.pident, h.nodes,
                a.alt_vaf, a.alt_abundance, a.taken_vaf,
                a.fork, strand_str(a.fork_strand),
                a.alt_node, strand_str(a.alt_strand),
                a.detour_nodes,
                a.length_delta,
                a.replaced_span,
                anchor_suffix,
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
/// figure cannot be checked by anyone. See `docs/issue_intra_accession_variance_lost.md`.
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
            writeln!(w, "          \"pident_inherited\": {},", h.pident_inherited)?;
            writeln!(w, "          \"length\": {},", h.seq.len())?;
            writeln!(w, "          \"nodes\": {},", h.nodes)?;
            writeln!(w, "          \"mean_ka\": {:.1},", h.mean_ka)?;
            // Emit `anchor_used` only when set. Absent field ≡ single-anchor legacy JSON.
            if let Some(rid) = &h.anchor_used {
                writeln!(w, "          \"vaf\": {:.2},", h.vaf)?;
                writeln!(w, "          \"anchor_used\": \"{}\"", esc(rid))?;
            } else {
                writeln!(w, "          \"vaf\": {:.2}", h.vaf)?;
            }
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
/// The intra-host SNP sites for every result, computed once.
///
/// `write_vcfs` and `write_snv_table` both need these, and each used to call `snv_sites`
/// itself, so with both flags on every accession paid for two identical alignments. Measured
/// at 11.1 ms each for pfcrt and 15.7 ms for hiv1_pol, which is 3.1 s of pure duplication on a
/// 200-accession hiv1_pol panel.
pub fn snv_sites_all(summary: &Summary, query: &Query) -> Vec<Vec<VcfRow>> {
    summary.results.iter().map(|r| snv_sites(r, query)).collect()
}

pub fn write_vcfs(
    summary: &Summary,
    query: &Query,
    out_dir: &Path,
    sites: &[Vec<VcfRow>],
) -> Result<()> {
    std::fs::create_dir_all(out_dir)?;
    let chrom = chrom_of(query);

    // The filename is the accession, which `docs/vcf_example.md` documents and the committed
    // examples rely on, so duplicates in the input list land on one path. Their contents are
    // for the same accession and so should agree, but say it rather than silently overwrite.
    {
        let mut seen = std::collections::HashSet::new();
        let dupes: Vec<&str> = summary
            .results
            .iter()
            .filter(|r| !seen.insert(r.accession.as_str()))
            .map(|r| r.accession.as_str())
            .collect();
        if !dupes.is_empty() {
            eprintln!(
                "[warn] {} duplicate accession(s) in the results ({:?}); each writes \
                 <out-dir>/<acc>.vcf, so only the last of each survives",
                dupes.len(),
                dupes
            );
        }
    }

    for (r, rows) in summary.results.iter().zip(sites) {
        if rows.is_empty() {
            continue;
        }

        let path = out_dir.join(format!("{}.vcf", r.accession));
        let mut w = BufWriter::new(std::fs::File::create(&path)?);
        write_vcf_header(&mut w, chrom, query.seq.len(), &r.accession)?;
        for row in rows {
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
pub fn write_snv_table(
    summary: &Summary,
    query: &Query,
    path: &Path,
    sites: &[Vec<VcfRow>],
) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let chrom = chrom_of(query);
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    writeln!(w, "gene\tpos\tref\talt\taccession\taf\tka_alt\tka_ref\tfork")?;
    for (r, rows) in summary.results.iter().zip(sites) {
        for row in rows {
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

    // Coalesce alts at the same query position into one multiallelic record. A triallelic
    // on-path fork produces two `AltHaplotype`s that differ from the primary at the SAME
    // position, so each becomes a row at that pos. Emitting them as separate VCF lines would
    // duplicate POS (and the major allele) with a contradictory RO, which is a malformed
    // record and makes `bcftools merge` warn or build a wrong matrix. Merge them: union the
    // alt alleles by base, and take the reference allele's coverage from whichever row actually
    // observed it (the others report 0 because the gene base is not their allele).
    let mut merged: Vec<VcfRow> = Vec::with_capacity(rows.len());
    for row in rows {
        if let Some(last) = merged.last_mut() {
            if last.pos == row.pos {
                last.ref_ka = last.ref_ka.max(row.ref_ka);
                for (base, af, ka) in row.alts {
                    if !last.alts.iter().any(|&(b, _, _)| b == base) {
                        last.alts.push((base, af, ka));
                    }
                }
                continue;
            }
        }
        merged.push(row);
    }
    merged
}

pub struct VcfRow {
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
/// Map each primary base to its query coordinate, exactly.
///
/// The two sequences are the same gene at high identity, so the optimal path hugs the main
/// diagonal and a band a few dozen cells wide contains it. That matters because the exact
/// implementation below allocates a full `(n+1)*(m+1)` byte traceback matrix: measured at
/// 1.6 MB and 11.1 ms for pfcrt, but **381 MB and 978 ms** at the 20 kb the guard permits,
/// which is a large single allocation on top of the pipeline's compressed-buffer ceiling.
///
/// Banding makes that `(n+1)*w` for a band width `w`, so pfcrt drops from 1.6 MB to ~170 KB.
///
/// # Why this cannot silently return a worse alignment
///
/// A band is only a speedup if the optimum stays inside it. The traceback reports whether it
/// touched a band edge; if it did, the band may have clipped the true optimum, so the band is
/// doubled and the alignment redone. When the band reaches full width this calls `align_full`,
/// so the worst case is exactly the old behaviour and never a wrong answer.
fn align_primary_to_query(primary: &[u8], query: &[u8]) -> Vec<Option<usize>> {
    let (n, m) = (primary.len(), query.len());
    let full = n.max(m);
    // The band has to span the length difference before it can absorb any drift on top.
    let mut band = (n.abs_diff(m) + 16).max(32);
    while band < full {
        if let Some(map) = align_banded(primary, query, band) {
            return map;
        }
        band = band.saturating_mul(2);
    }
    align_full(primary, query)
}

/// One banded pass. `None` means the traceback reached a band edge, so the result may not be
/// optimal and the caller must widen.
///
/// Cells are addressed by their offset from the diagonal, `d = j - i`, which is bounded by
/// `d_min..=d_max`. The band must contain `d = 0` (the start) and `d = m - n` (the end), which
/// is why those two bracket the range before `band` widens it.
fn align_banded(primary: &[u8], query: &[u8], band: usize) -> Option<Vec<Option<usize>>> {
    let (n, m) = (primary.len(), query.len());
    let delta = m as isize - n as isize;
    let d_min = delta.min(0) - band as isize;
    let d_max = delta.max(0) + band as isize;
    let w = (d_max - d_min + 1) as usize;

    const GAP: i32 = -2;
    const NEG: i32 = i32::MIN / 4; // room to add GAP without overflowing

    let mut dir = vec![0u8; (n + 1) * w];
    let mut prev = vec![NEG; w];
    let mut curr = vec![NEG; w];

    let col = |i: usize, j: usize| -> usize { (j as isize - i as isize - d_min) as usize };
    let j_range = |i: usize| -> (usize, usize) {
        let lo = (i as isize + d_min).max(0) as usize;
        let hi = ((i as isize + d_max).min(m as isize)).max(0) as usize;
        (lo, hi)
    };

    // Row 0: only left moves are possible.
    let (lo0, hi0) = j_range(0);
    for j in lo0..=hi0 {
        prev[col(0, j)] = GAP * j as i32;
        dir[col(0, j)] = 2;
    }

    for i in 1..=n {
        let (lo, hi) = j_range(i);
        for c in curr.iter_mut() {
            *c = NEG;
        }
        for j in lo..=hi {
            let c = col(i, j);
            if j == 0 {
                curr[c] = GAP * i as i32;
                dir[i * w + c] = 1;
                continue;
            }
            // diag (i-1, j-1) keeps d; up (i-1, j) is d+1; left (i, j-1) is d-1.
            let diag = if prev[c] == NEG { NEG } else { prev[c] + if primary[i - 1] == query[j - 1] { 1 } else { -1 } };
            let up = if c + 1 < w && prev[c + 1] != NEG { prev[c + 1] + GAP } else { NEG };
            let left = if c > 0 && curr[c - 1] != NEG { curr[c - 1] + GAP } else { NEG };
            let (mut best, mut d) = (diag, 0u8);
            if up > best {
                best = up;
                d = 1;
            }
            if left > best {
                best = left;
                d = 2;
            }
            curr[c] = best;
            dir[i * w + c] = d;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    // Traceback from (n, m), refusing the result if it ever sat on a band edge.
    let mut map = vec![None; n];
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        let d_here = j as isize - i as isize;
        if d_here <= d_min || d_here >= d_max {
            return None;
        }
        let c = col(i, j);
        let d = dir[i * w + c];
        if i > 0 && j > 0 && d == 0 {
            map[i - 1] = Some(j - 1);
            i -= 1;
            j -= 1;
        } else if i > 0 && d == 1 {
            map[i - 1] = None;
            i -= 1;
        } else if j > 0 {
            j -= 1;
        } else {
            return None; // dir says left but there is nothing left: outside the band
        }
    }
    Some(map)
}

/// Linear gap penalty; match +1, mismatch -1, gap -2. The two sequences are the same gene at
/// high identity, so the alignment is near-diagonal and the exact scoring barely matters; what
/// matters is that a rare indel shifts coordinates correctly rather than being ignored.
fn align_full(primary: &[u8], query: &[u8]) -> Vec<Option<usize>> {
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
mod banded_equivalence_tests {
    use super::{align_full, align_primary_to_query};

    fn lcg(seed: &mut u64) -> u64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *seed >> 33
    }

    fn seq(n: usize, seed: &mut u64) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[(lcg(seed) % 4) as usize]).collect()
    }

    /// The band is a speedup, never a different answer. If this ever fails, the widening
    /// loop in `align_primary_to_query` is not detecting a clipped optimum.
    #[test]
    fn banded_matches_the_exact_alignment_on_random_pairs() {
        let mut seed = 0xC0FFEEu64;
        for case in 0..400 {
            let n = 20 + (lcg(&mut seed) % 180) as usize;
            let a = seq(n, &mut seed);
            let mut b = a.clone();

            // substitutions
            for _ in 0..(lcg(&mut seed) % 12) {
                let p = (lcg(&mut seed) as usize) % b.len();
                b[p] = b"ACGT"[(lcg(&mut seed) % 4) as usize];
            }
            // indels, including runs long enough to push the path off the diagonal
            for _ in 0..(lcg(&mut seed) % 4) {
                let p = (lcg(&mut seed) as usize) % b.len();
                let len = 1 + (lcg(&mut seed) % 9) as usize;
                if lcg(&mut seed).is_multiple_of(2) {
                    let ins = seq(len, &mut seed);
                    b.splice(p..p, ins);
                } else {
                    let end = (p + len).min(b.len());
                    b.drain(p..end);
                }
            }
            if b.is_empty() {
                continue;
            }
            assert_eq!(
                align_primary_to_query(&a, &b),
                align_full(&a, &b),
                "case {case}: banded disagreed with the exact DP (n={}, m={})",
                a.len(),
                b.len()
            );
        }
    }

    /// A pair that shares almost nothing forces the band all the way open; the fallback must
    /// still produce the exact answer rather than a clipped one.
    #[test]
    fn a_wildly_divergent_pair_falls_back_to_the_exact_dp() {
        let a = vec![b'A'; 120];
        let mut b = vec![b'T'; 40];
        b.extend(std::iter::repeat_n(b'A', 30));
        assert_eq!(align_primary_to_query(&a, &b), align_full(&a, &b));
    }

    /// Degenerate shapes the band arithmetic has to survive.
    #[test]
    fn edge_shapes_agree() {
        let a = b"ACGTACGTAC".to_vec();
        for b in [
            b"A".to_vec(),
            b"ACGTACGTAC".to_vec(),
            b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec(),
            vec![b'G'; 1],
        ] {
            assert_eq!(align_primary_to_query(&a, &b), align_full(&a, &b), "b len {}", b.len());
        }
    }
}

#[cfg(test)]
mod vcf_tests {
    use super::{align_primary_to_query, snv_sites};
    use crate::types::{
        AccessionResult, AltHaplotype, Haplotype, Outcome, Query, Strand, Timing,
    };
    use std::path::PathBuf;

    fn hap(seq: &[u8]) -> Haplotype {
        Haplotype {
            id: "h".into(),
            seq: seq.to_vec(),
            cov: 1.0,
            pident: 100.0,
            pident_inherited: false,
            nodes: 1,
            mean_ka: 0.0,
            vaf: 1.0,
            anchor_used: None,
        }
    }

    fn alt(seq: &[u8], alt_vaf: f64, alt_ka: f32) -> AltHaplotype {
        AltHaplotype {
            haplotype: hap(seq),
            fork: 7,
            fork_strand: Strand::Forward,
            alt_node: 9,
            alt_strand: Strand::Forward,
            taken_abundance: 100.0,
            taken_vaf: 0.60,
            alt_abundance: alt_ka,
            alt_vaf,
            detour_nodes: 1,
            length_delta: 0,
            replaced_span: 1,
        }
    }

    fn result(primary: &[u8], alts: Vec<AltHaplotype>) -> AccessionResult {
        AccessionResult {
            index: 0,
            accession: "ACC".into(),
            outcome: Outcome::Success,
            n_unitigs: 0,
            n_anchors: 0,
            subgraph_nodes: 0,
            haplotype: Some(hap(primary)),
            haplotypes: vec![hap(primary)],
            alts,
            bubbles: Vec::new(),
            timing: Timing::default(),
            error: None,
            anchor_used: None,
        }
    }

    fn query(seq: &[u8]) -> Query {
        Query { name: "g".into(), record_id: "g".into(), seq: seq.to_vec(), path: PathBuf::from("g.fa"), n_records: 1 }
    }

    /// A triallelic on-path fork yields two alt haplotypes that differ from the primary at the
    /// same position. They must coalesce into ONE multiallelic site, not two records at the same
    /// POS (which would be a malformed VCF and break `bcftools merge`).
    #[test]
    fn triallelic_site_coalesces_into_one_multiallelic_row() {
        let primary = b"TTTATTT"; // major allele A at index 3 (pos 4)
        let q = query(b"TTTATTT"); // gene base == major, so REF = A
        let alts = vec![
            alt(b"TTTCTTT", 0.25, 25.0), // minor C
            alt(b"TTTGTTT", 0.15, 15.0), // minor G
        ];
        let rows = snv_sites(&result(primary, alts), &q);
        assert_eq!(rows.len(), 1, "the two same-position alts must merge to one record");
        assert_eq!(rows[0].pos, 4);
        let bases: Vec<u8> = rows[0].alts.iter().map(|&(b, _, _)| b).collect();
        assert_eq!(bases, vec![b'C', b'G'], "both minor alleles present, deduplicated");
        assert_eq!(rows[0].ref_base, b'A');
        assert_eq!(rows[0].ref_ka, 100.0, "RO is the major (=REF) coverage, once");
    }

    /// Two forks at DIFFERENT positions stay two records.
    #[test]
    fn distinct_positions_stay_separate() {
        let primary = b"TTATTAT"; // A at index 2 and index 5
        let q = query(b"TTATTAT");
        let alts = vec![alt(b"TTCTTAT", 0.20, 20.0), alt(b"TTATTCT", 0.30, 30.0)];
        let rows = snv_sites(&result(primary, alts), &q);
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].pos, 3);
        assert_eq!(rows[1].pos, 6);
    }

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

#[cfg(test)]
mod identity_label_tests {
    use super::{write_alts, write_fasta};
    use crate::pipeline::Summary;
    use crate::types::{AccessionResult, AltHaplotype, Haplotype, Outcome, Strand, Timing};

    fn hap(id: &str, inherited: bool) -> Haplotype {
        Haplotype {
            id: id.into(),
            seq: b"ACGT".to_vec(),
            cov: 1.0,
            pident: 98.1,
            pident_inherited: inherited,
            nodes: 1,
            mean_ka: 12.0,
            vaf: 1.0,
            anchor_used: None,
        }
    }

    fn summary(primary: Haplotype, alts: Vec<AltHaplotype>) -> Summary {
        Summary {
            total: 1,
            success: 1,
            results: vec![AccessionResult {
                index: 0,
                accession: "ACC".into(),
                outcome: Outcome::Success,
                n_unitigs: 0,
                n_anchors: 0,
                subgraph_nodes: 0,
                haplotype: Some(primary.clone()),
                haplotypes: vec![primary],
                alts,
                bubbles: Vec::new(),
                timing: Timing::default(),
                error: None,
                anchor_used: None,
            }],
        }
    }

    fn alt_of(h: Haplotype) -> AltHaplotype {
        AltHaplotype {
            haplotype: h,
            fork: 1,
            fork_strand: Strand::Forward,
            alt_node: 2,
            alt_strand: Strand::Forward,
            taken_abundance: 10.0,
            taken_vaf: 0.6,
            alt_abundance: 5.0,
            alt_vaf: 0.4,
            detour_nodes: 1,
            length_delta: 0,
            replaced_span: 1,
        }
    }

    fn tmp(name: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("logan_ident_{}_{}", std::process::id(), name));
        p
    }

    /// A haplotype whose own alignment produced the number keeps `pident=`; one carrying the
    /// pre-splice route's number must not, or a consumer reads it as a measurement of this
    /// sequence.
    #[test]
    fn primary_identity_key_says_whether_it_was_measured() {
        let p = tmp("measured.fasta");
        write_fasta(&summary(hap("ACC_primary", false), vec![]), &p).unwrap();
        let s = std::fs::read_to_string(&p).unwrap();
        assert!(s.contains("pident=98.1"), "{s}");
        assert!(!s.contains("inherited_pident"), "{s}");
        std::fs::remove_file(&p).ok();

        let p = tmp("inherited.fasta");
        write_fasta(&summary(hap("ACC_primary", true), vec![]), &p).unwrap();
        let s = std::fs::read_to_string(&p).unwrap();
        assert!(s.contains("inherited_pident=98.1"), "{s}");
        assert!(
            !s.split_whitespace().any(|w| w.starts_with("pident=")),
            "a spliced haplotype must not publish a bare pident=: {s}"
        );
        std::fs::remove_file(&p).ok();
    }

    /// Every alt is a spliced path, so its identity is always the primary's.
    #[test]
    fn alt_identity_is_always_labelled_inherited() {
        let p = tmp("alts.fasta");
        let alts = vec![alt_of(hap("ACC_alt1", false))];
        write_alts(&summary(hap("ACC_primary", false), alts), &p).unwrap();
        let s = std::fs::read_to_string(&p).unwrap();
        assert!(s.contains("inherited_pident=98.1"), "{s}");
        assert!(
            !s.split_whitespace().any(|w| w.starts_with("pident=")),
            "alt records must not publish a bare pident=: {s}"
        );
        std::fs::remove_file(&p).ok();
    }
}
