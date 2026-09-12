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
use crate::types::{AccessionResult, Outcome, Strand};
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
