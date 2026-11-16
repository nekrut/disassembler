//! Everything that leaves the pipeline.
//!
//! IN   Summary
//! OUT  reconstructed.fasta, bubbles.json, <query>.summary.json
//!
//! `reconstructed.fasta` is consumed by `disassembler msa`, which feeds MAFFT and FastTree,
//! so its record format is a fixed contract: downstream parsers read it, and the
//! `bench/golden/` fixtures are this file verbatim.
//!
//! JSON is written by hand rather than through serde. The shapes are small and fixed, and
//! keeping one exact layout (two-space indent, fields always in the same order) matters more
//! here than convenience, because summaries are meant to be diffable.

use crate::config::Config;
use crate::error::Result;
use crate::pipeline::Summary;
use crate::types::{AccessionResult, Outcome, Payload, Query, Strand, Timing};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

fn outcome_str(o: Outcome) -> &'static str {
    match o {
        Outcome::Success => "success",
        Outcome::NoAnchors => "no_anchors",
        Outcome::NoPath => "no_path",
        Outcome::NoLoganData => "no_logan_data",
        Outcome::TooLarge => "too_large",
        Outcome::SubgraphTooLarge => "subgraph_too_large",
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

/// One record per reconstructed haplotype, in input order: `pipeline::run` sorts results by
/// `AccessionResult::index` before they reach here, so the file is byte-reproducible.
///
/// Header format is part of the output contract (downstream parsers, `bench/golden/`):
/// `>{id} cov={:.3} pident={:.1} nodes={} ka={:.1} vaf={:.2}`, plus ` anchor={rid}` in
/// multi-anchor mode, sequence on a single line.
///
/// WV-24: this doc stopped at `nodes=` while the writer has emitted `ka=` and `vaf=` for some
/// time, so anyone writing a parser from it got the fields wrong.
///
/// The identity key is `pident=` when it was measured on this haplotype's own alignment and
/// `inherited_pident=` when it came from the pre-splice route (see `Haplotype::pident`), so a
/// consumer reading `pident=` never gets a number about a different sequence.
pub fn write_fasta(summary: &Summary, path: &Path) -> Result<()> {
    let mut w = create(path)?;
    for r in &summary.results {
        for h in &r.haplotypes {
            // anchor_used is emitted only when set (multi-anchor mode). Single-anchor
            // runs keep the header byte-identical to the pre-multi-anchor shape so all
            // golden fixtures diff clean.
            let anchor_suffix = anchor_suffix(h);
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
/// Deliberately a **separate file** from `reconstructed.fasta`. `disassembler msa`
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
/// ```
///
/// WV-24: that fence used to stay open, so the prose below rendered as code.
///
/// `inherited_pident` is the PRIMARY's alignment identity, not this allele's: an alt is a
/// spliced path with no alignment of its own. It was published as `pident=` per record, which
/// read as a measurement of that allele.
///
/// ` ka_sat=1` is appended when the fork is at Logan's abundance ceiling (WV-03), and
/// ` anchor={rid}` in multi-anchor mode.
pub fn write_alts(summary: &Summary, path: &Path) -> Result<()> {
    let mut w = create(path)?;
    for r in &summary.results {
        for a in &r.alts {
            let h = &a.haplotype;
            let anchor_suffix = anchor_suffix(h);
            // WV-03: appended only when the fork is clipped, so an unsaturated record is
            // byte-identical to what 0.1.x wrote.
            let ka_sat_suffix = if a.ka_saturated { " ka_sat=1" } else { "" };
            writeln!(
                w,
                ">{} cov={:.3} inherited_pident={:.1} nodes={} vaf={:.3} alt_kaf={:.1} \
                 primary_vaf={:.3} fork={}{} alt_node={}{} detour={} dlen={:+} span={}{}{}",
                h.id, h.cov, h.pident, h.nodes,
                a.alt_vaf, a.alt_abundance, a.taken_vaf,
                a.fork, strand_str(a.fork_strand),
                a.alt_node, strand_str(a.alt_strand),
                a.detour_nodes,
                a.length_delta,
                a.replaced_span,
                anchor_suffix,
                ka_sat_suffix,
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
/// With only a bubble count, a reported quasispecies figure could not be checked by anyone.
/// See `docs/issue_intra_accession_variance_lost.md`.
pub fn write_bubbles(summary: &Summary, path: &Path) -> Result<()> {
    let mut w = create(path)?;
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
            // WV-03. Written always, including false, so a consumer can tell an unsaturated
            // bubble from one produced before the flag existed.
            writeln!(w, "        \"ka_saturated\": {},", b.ka_saturated)?;
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

/// The first bytes every walker-written VCF starts with.
///
/// WV-16. `write_vcfs` removed or truncated `<out-dir>/<ACC>.vcf` without ever checking who
/// wrote it. Point `-o` at a directory already holding per-sample VCFs from a read-based
/// caller, which use the same `<ACC>.vcf` names, and a run with no sites **deletes them** and
/// exits 0. Observed on a `##source=bcftools` file.
const VCF_SIGNATURE: &[u8] = b"##fileformat=VCFv4.2\n##source=logan-walker";

/// Whether this path is absent or was written by us, so replacing it loses nothing of
/// anyone else's.
fn walker_owned(path: &Path) -> bool {
    match std::fs::read(path) {
        Err(_) => true, // missing, or unreadable: `File::create` will report it properly
        Ok(bytes) => bytes.starts_with(VCF_SIGNATURE),
    }
}

/// A file that appears at its final path only once it is complete.
///
/// WV-18. Every writer opened the destination directly, so a write error after the compute
/// (a full disk, a directory in the way) left a truncated or missing file under the final
/// name, and because `main` used `?` between writers, the **summary was never written at
/// all**: the shard lost its per-accession status and timings even though the accessions had
/// reconstructed. A signal mid-write had the same effect. `stages::stream` already writes its
/// cache this way.
struct AtomicFile {
    tmp: std::path::PathBuf,
    dst: std::path::PathBuf,
    inner: Option<BufWriter<std::fs::File>>,
}

impl AtomicFile {
    fn create(dst: &Path) -> Result<Self> {
        if let Some(parent) = dst.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)?;
            }
        }
        let tmp = dst.with_extension(format!(
            "{}.tmp.{}",
            dst.extension().and_then(|e| e.to_str()).unwrap_or(""),
            std::process::id()
        ));
        let file = std::fs::File::create(&tmp)?;
        Ok(Self { tmp, dst: dst.to_path_buf(), inner: Some(BufWriter::new(file)) })
    }

    fn w(&mut self) -> &mut BufWriter<std::fs::File> {
        self.inner.as_mut().expect("writer used after commit")
    }

    /// Flush, sync and rename. Only now does anything appear at the final path.
    fn commit(mut self) -> Result<()> {
        let mut bw = self.inner.take().expect("commit called twice");
        bw.flush()?;
        bw.into_inner().map_err(|e| e.into_error())?.sync_all()?;
        std::fs::rename(&self.tmp, &self.dst)?;
        Ok(())
    }
}

impl Drop for AtomicFile {
    fn drop(&mut self) {
        // Uncommitted: the write failed or was interrupted, so leave nothing behind.
        if self.inner.is_some() {
            let _ = std::fs::remove_file(&self.tmp);
        }
    }
}

/// Remove this run's own leftovers from a writer that is switched off.
///
/// WV-17. Each optional writer ran only when its flag was set, so a rerun into the same
/// `--out-dir` without `--emit-vcf`, without `--emit-snv-table` or with `--no-emit-alts` left
/// the PREVIOUS run's `<ACC>.vcf`, `snvs.tsv` and `alt_haplotypes.fasta` sitting beside a
/// fresh primary and a fresh summary. Observed on a rerun two seconds later: the three
/// variant files still carried the first run's timestamps, and nothing in the outputs said
/// so. The emit flags are now in the summary too (WV-20), but a stale file is still a stale
/// file.
///
/// Only files this walker wrote are removed. `snvs.tsv` and the VCFs are identified by their
/// header; `alt_haplotypes.fasta` has no signature, so it is removed only when it sits beside
/// a summary we are about to overwrite, which is the same directory contract.
pub fn remove_disabled_outputs(cfg: &Config, out_dir: &Path, summary: &Summary) {
    let drop_it = |path: std::path::PathBuf, what: &str| {
        if path.exists() {
            match std::fs::remove_file(&path) {
                Ok(()) => eprintln!(
                    "[logan-walker] removed stale {} from an earlier run ({what} is off)",
                    path.display()
                ),
                Err(e) => eprintln!("[warn] cannot remove stale {}: {e}", path.display()),
            }
        }
    };
    if !cfg.emit_vcf {
        for r in &summary.results {
            let path = out_dir.join(format!("{}.vcf", r.accession));
            if path.exists() && walker_owned(&path) {
                drop_it(path, "--emit-vcf");
            }
        }
    }
    if !cfg.emit_snv_table {
        let path = out_dir.join("snvs.tsv");
        let ours = std::fs::read(&path)
            .map(|b| b.starts_with(b"gene\tpos\tref\talt\taccession"))
            .unwrap_or(false);
        if ours {
            drop_it(path, "--emit-snv-table");
        }
    }
    if !cfg.emit_alts {
        drop_it(out_dir.join("alt_haplotypes.fasta"), "--emit-alts");
    }
}

/// Open `path` for writing, creating its parent directory.
///
/// WV-24: five writers spelled this block out. The VCF writer uses `AtomicFile` instead,
/// which does the same thing and then some.
fn create(path: &Path) -> Result<BufWriter<std::fs::File>> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    Ok(BufWriter::new(std::fs::File::create(path)?))
}

/// The ` anchor={rid}` suffix a FASTA header carries in multi-anchor mode, or nothing.
///
/// WV-24: `write_fasta` and `write_alts` built this identically.
fn anchor_suffix(h: &crate::types::Haplotype) -> String {
    match &h.anchor_used {
        Some(rid) => format!(" anchor={rid}"),
        None => String::new(),
    }
}

/// A float as JSON, or `null` when it is not finite.
///
/// WV-15. `format!("{:?}")` on an f64 yields the bare tokens `inf` and `NaN`, which **no JSON
/// parser accepts**, so one such value made the whole summary unreadable while the run exited
/// 0. The flags that could carry one are now validated, but a value can still arrive from a
/// division by zero, and a summary that cannot be parsed is worse than one with a null in it.
fn jf(v: f64) -> String {
    if v.is_finite() {
        format!("{v:?}")
    } else {
        "null".to_string()
    }
}

/// Lowercase hex sha256 of a sequence, uppercased first.
///
/// This is how a summary names a sequence (`docs/summaries.md`). It answers the questions the
/// old `seq` field was used for, equality between runs, counting distinct primaries, and
/// joining a record to its FASTA entry, without the summary carrying a second copy of
/// `reconstructed.fasta`.
fn seq_sha256(seq: &[u8]) -> String {
    let upper: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
    crate::types::sha256_hex(&upper)
}

/// `YYYY-MM-DDTHH:MM:SSZ`, second resolution, UTC.
///
/// Written out rather than pulled in: a date crate would be the walker's only dependency that
/// is not doing I/O, and this is 12 lines of arithmetic (Howard Hinnant's `civil_from_days`).
pub(crate) fn iso8601(t: SystemTime) -> String {
    let secs = t.duration_since(UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0);
    let (days, tod) = ((secs / 86_400) as i64, secs % 86_400);
    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z",
        tod / 3600,
        (tod % 3600) / 60,
        tod % 60
    )
}

/// `{"p10": .., "p25": .., "p50": .., "p75": .., "p90": .., "mean": ..}` on one line, or
/// `null` when there is nothing to describe. The same five quantiles and the same key names
/// every stage uses, so two summaries can be compared without a translation table.
fn quantiles(mut v: Vec<f64>) -> String {
    if v.is_empty() {
        return "null".to_string();
    }
    v.sort_by(f64::total_cmp);
    // Linear interpolation between the two neighbouring order statistics (numpy "linear",
    // `statistics.quantiles(method="inclusive")`), which is what the Python side does.
    //
    // PL-15. Both sides used to pick a rounded index, so an even cohort's median was one of
    // the two middle values. Worse for a contract that says `p50` means `p50` on both
    // pathways: Rust's `round` goes half away from zero and Python's goes to even, so the
    // two stages disagreed on exactly the inputs where the choice mattered.
    let at = |q: f64| {
        let pos = (v.len() as f64 - 1.0) * q;
        let lo = pos.floor() as usize;
        let hi = (lo + 1).min(v.len() - 1);
        v[lo] + (v[hi] - v[lo]) * (pos - lo as f64)
    };
    let mean = v.iter().sum::<f64>() / v.len() as f64;
    format!(
        "{{\"p10\": {:.4}, \"p25\": {:.4}, \"p50\": {:.4}, \"p75\": {:.4}, \"p90\": {:.4}, \
         \"mean\": {:.4}}}",
        at(0.10),
        at(0.25),
        at(0.50),
        at(0.75),
        at(0.90),
        mean
    )
}

/// What produced a summary, and what it consumed.
///
/// 0.1.1 changes primaries (`--primary-metric matched`, the header fix, the prune coverage
/// rule), so a summary has to say which walker and which settings wrote it; until 0.1.1 one
/// on disk could not. 0.2.0 adds what it read, how long it took, and which of the requested
/// accessions it never reached.
pub struct RunInfo<'a> {
    /// GraphAligner's own `--version` line, and the binary that produced it. The version alone
    /// does not say which of several installed copies ran.
    ///
    /// `None` for a stage that never runs it, which is `--fetch-only`. Naming a tool that did
    /// no work would be a false provenance claim, so `tool.external` comes out empty instead.
    pub graphaligner: Option<(&'a str, &'a str)>,
    pub cfg: &'a Config,
    pub queries: &'a [Query],
    pub started: SystemTime,
    pub wall_seconds: f64,
    /// Accessions asked for, before any were skipped. `n_requested - n_records` is itself a
    /// discard and is reported as one.
    pub n_requested: usize,
}

fn write_tool(w: &mut impl Write, run: &RunInfo) -> Result<()> {
    writeln!(w, "  \"tool\": {{")?;
    writeln!(w, "    \"name\": \"logan-walker\",")?;
    writeln!(w, "    \"version\": \"{}\",", env!("CARGO_PKG_VERSION"))?;
    // Set by `disassembler walker`; absent when the binary is run directly.
    match std::env::var("DISASSEMBLER_VERSION") {
        Ok(v) if !v.is_empty() => writeln!(w, "    \"disassembler\": \"{}\",", esc(&v))?,
        _ => writeln!(w, "    \"disassembler\": null,")?,
    }
    match run.graphaligner {
        Some((version, path)) => writeln!(
            w,
            "    \"external\": [{{\"name\": \"GraphAligner\", \"version\": \"{}\", \"path\": \"{}\"}}]",
            esc(version),
            esc(path)
        )?,
        None => writeln!(w, "    \"external\": []")?,
    }
    writeln!(w, "  }},")?;
    Ok(())
}

fn write_inputs(w: &mut impl Write, run: &RunInfo) -> Result<()> {
    let c = run.cfg;
    writeln!(w, "  \"inputs\": {{")?;
    // `null`, not `[]`, for a stage with no reference: `docs/summaries.md` reserves the empty
    // array for "this stage takes references and was given none", which cannot happen here.
    // `--fetch-only` is the stage that has none at all.
    if run.queries.is_empty() {
        writeln!(w, "    \"references\": null,")?;
    } else {
        writeln!(w, "    \"references\": [")?;
        for (i, q) in run.queries.iter().enumerate() {
            writeln!(
                w,
                "      {{\"id\": \"{}\", \"length\": {}, \"sha256\": \"{}\"}}{}",
                esc(&q.record_id),
                q.seq.len(),
                seq_sha256(&q.seq),
                if i + 1 == run.queries.len() { "" } else { "," }
            )?;
        }
        writeln!(w, "    ],")?;
    }
    // A cached run and a streamed one can produce the same numbers from different bytes, so
    // the summary says which it was.
    let (kind, id) = match (&c.unitigs_cache, c.offline) {
        (Some(p), true) => ("cache".to_string(), p.display().to_string()),
        (Some(p), false) => ("cache+s3".to_string(), p.display().to_string()),
        (None, _) => ("s3".to_string(), "https://logan-pub.s3.amazonaws.com/u/".to_string()),
    };
    writeln!(
        w,
        "    \"source\": {{\"kind\": \"{}\", \"id\": \"{}\", \"n_requested\": {}}}",
        kind,
        esc(&id),
        run.n_requested
    )?;
    writeln!(w, "  }},")?;
    Ok(())
}

fn write_settings(w: &mut impl Write, run: &RunInfo) -> Result<()> {
    let c = run.cfg;
    // Everything that changes which accessions are processed or what they produce. Workers,
    // timeouts and paths change neither, and are left out.
    let settings: [(&str, String); 24] = [
        ("k", c.k.to_string()),
        ("graph_k", c.graph_k.to_string()),
        ("min_shared_kmers", c.min_shared_kmers.to_string()),
        ("hops", c.hops.to_string()),
        ("subgraph_nodes_per_base_guard", jf(c.subgraph_nodes_per_base_guard)),
        ("max_size_mb", (c.max_compressed_bytes / (1024 * 1024)).to_string()),
        ("min_cov", jf(c.min_cov)),
        ("primary_metric", format!("\"{}\"", c.primary_metric.as_str())),
        ("prune_frac", jf(c.prune_frac)),
        ("prefer_majority", c.prefer_majority.to_string()),
        ("min_abund", jf(c.min_abund as f64)),
        ("min_abund_frac", jf(c.min_abund_frac)),
        ("min_vaf", jf(c.min_vaf)),
        ("alt_max_detour", c.alt_max_detour.to_string()),
        ("max_alts", c.max_alts.to_string()),
        ("indel_alts", c.allow_indel_alts.to_string()),
        ("max_alt_indel", c.max_alt_indel.to_string()),
        ("allow_frameshift_alts", c.allow_frameshift_alts.to_string()),
        ("graphaligner_threads", c.graphaligner_threads.to_string()),
        ("rescue_minority", c.rescue_minority.to_string()),
        ("summary_seqs", c.summary_seqs.to_string()),
        // WV-20: which outputs this run was asked for. Without these, a summary beside a
        // missing VCF cannot say whether the run found no sites or was never asked to look.
        ("emit_alts", c.emit_alts.to_string()),
        ("emit_vcf", c.emit_vcf.to_string()),
        ("emit_snv_table", c.emit_snv_table.to_string()),
    ];
    writeln!(w, "  \"settings\": {{")?;
    for (i, (k, v)) in settings.iter().enumerate() {
        writeln!(w, "    \"{k}\": {v}{}", if i + 1 == settings.len() { "" } else { "," })?;
    }
    writeln!(w, "  }},")?;
    Ok(())
}

fn write_run_block(w: &mut impl Write, s: &Summary, run: &RunInfo) -> Result<()> {
    let n = s.results.len();
    let ok = s.success;
    let status = if n > 0 && ok == n {
        "ok"
    } else if ok == 0 {
        "failed"
    } else {
        "partial"
    };
    let sum = |f: fn(&Timing) -> f64| s.results.iter().map(|r| f(&r.timing)).sum::<f64>();
    writeln!(w, "  \"run\": {{")?;
    writeln!(w, "    \"status\": \"{status}\",")?;
    writeln!(w, "    \"exit_code\": {},", i32::from(status == "failed"))?;
    writeln!(w, "    \"started\": \"{}\",", iso8601(run.started))?;
    writeln!(w, "    \"wall_seconds\": {:.3},", run.wall_seconds)?;
    writeln!(
        w,
        "    \"stage_seconds\": {{\"stream\": {:.3}, \"parse\": {:.3}, \"anchor\": {:.3}, \
         \"extract\": {:.3}, \"align\": {:.3}, \"stitch\": {:.3}, \"phase\": {:.3}}},",
        sum(|t| t.stream),
        sum(|t| t.parse),
        sum(|t| t.anchor),
        sum(|t| t.extract),
        sum(|t| t.align),
        sum(|t| t.stitch),
        sum(|t| t.phase)
    )?;
    writeln!(w, "    \"n_records\": {n},")?;
    writeln!(w, "    \"n_ok\": {ok}")?;
    writeln!(w, "  }},")?;
    Ok(())
}

/// Every rejection the run made, counts only.
///
/// Two groups because there are two things being discarded: whole accessions, and candidate
/// alleles within an accession. Every key is always written, including the zeroes, so a
/// consumer can tell "this did not happen" from "this version did not measure it".
fn write_discards(w: &mut impl Write, s: &Summary, run: &RunInfo) -> Result<()> {
    let count = |o: Outcome| s.results.iter().filter(|r| r.outcome == o).count();
    let v = |f: fn(&crate::types::VariantStats) -> usize| {
        s.results.iter().map(|r| f(&r.variants)).sum::<usize>()
    };
    writeln!(w, "  \"discards\": {{")?;
    writeln!(
        w,
        "    \"accessions\": {{\"never_reached\": {}, \"no_anchors\": {}, \"no_path\": {}, \
         \"no_logan_data\": {}, \"too_large\": {}, \"subgraph_too_large\": {}, \
         \"error\": {}}},",
        run.n_requested.saturating_sub(s.results.len()),
        count(Outcome::NoAnchors),
        count(Outcome::NoPath),
        count(Outcome::NoLoganData),
        count(Outcome::TooLarge),
        count(Outcome::SubgraphTooLarge),
        count(Outcome::Failed)
    )?;
    writeln!(
        w,
        "    \"alleles\": {{\"beyond_max_alts\": {}, \
         \"no_rejoin\": {}, \"not_substitution\": {}, \"over_indel_cap\": {}, \
         \"frameshifting\": {}, \"window_not_shiftable\": {}, \"identical\": {}}},",
        v(|x| x.beyond_cap),
        v(|x| x.no_rejoin),
        v(|x| x.not_substitution),
        v(|x| x.over_indel_cap),
        v(|x| x.frameshifting),
        v(|x| x.window_not_shiftable),
        v(|x| x.identical)
    )?;
    // WV-02: forks the run left with `_primary` on the rarer allele. A third group, because
    // these are not candidate alleles that were refused; they are sites where the published
    // sequence is known to carry the minority base.
    let m = |f: fn(&crate::types::MinorityLeft) -> usize| {
        s.results.iter().map(|r| f(&r.majority.left)).sum::<usize>()
    };
    writeln!(
        w,
        "    \"minority_sites\": {{\"alt_below_floor\": {}, \
         \"not_outright_majority\": {}, \"exempt_minority\": {}, \
         \"saturated\": {}, \
         \"no_rejoin\": {}, \"length_refused\": {}, \"window_not_shiftable\": {}, \
         \"coverage_gate\": {}, \"identical\": {}}}",
        m(|x| x.alt_below_floor),
        m(|x| x.not_outright_majority),
        m(|x| x.exempt_minority),
        m(|x| x.saturated),
        m(|x| x.no_rejoin),
        m(|x| x.length_refused),
        m(|x| x.window_not_shiftable),
        m(|x| x.coverage_gate),
        m(|x| x.identical)
    )?;
    writeln!(w, "  }},")?;
    Ok(())
}

/// What the run produced, as a population rather than a list.
///
/// None of this was in a walker summary before 0.2.0: every population-level number had to be
/// computed by the consumer, while the LexicMap pathway reported its own. `distinct_primaries`
/// is the one that matters most for a cohort, and it is why `seq_sha256` is a record field.
fn write_cohort(w: &mut impl Write, s: &Summary) -> Result<()> {
    let mut hashes: Vec<String> = Vec::new();
    let (mut covs, mut pidents) = (Vec::new(), Vec::new());
    for r in &s.results {
        if let Some(h) = r.haplotypes.first() {
            hashes.push(seq_sha256(&h.seq));
            covs.push(h.cov);
            pidents.push(h.pident);
        }
    }
    hashes.sort();
    let distinct = {
        let mut d = hashes.clone();
        d.dedup();
        d.len()
    };
    let singletons = {
        let mut n = 0usize;
        let mut i = 0usize;
        while i < hashes.len() {
            let j = hashes[i..].partition_point(|x| *x == hashes[i]) + i;
            if j - i == 1 {
                n += 1;
            }
            i = j;
        }
        n
    };
    writeln!(w, "  \"cohort\": {{")?;
    writeln!(w, "    \"reconstructed\": {},", hashes.len())?;
    writeln!(
        w,
        "    \"provisional\": {},",
        s.results.iter().filter(|r| r.provisional).count()
    )?;
    writeln!(
        w,
        "    \"majority_corrected\": {},",
        s.results.iter().filter(|r| r.majority_corrections > 0).count()
    )?;
    writeln!(
        w,
        "    \"majority_corrections\": {},",
        s.results.iter().map(|r| r.majority_corrections).sum::<usize>()
    )?;
    writeln!(
        w,
        "    \"minority_rescues\": {},",
        s.results.iter().map(|r| r.majority.minority_rescues).sum::<usize>()
    )?;
    writeln!(
        w,
        "    \"accessions_with_minority_sites\": {},",
        s.results.iter().filter(|r| r.majority.left.total() > 0).count()
    )?;
    writeln!(
        w,
        "    \"accessions_depth_saturated\": {},",
        s.results.iter().filter(|r| r.depth_saturated).count()
    )?;
    writeln!(w, "    \"distinct_primaries\": {distinct},")?;
    writeln!(w, "    \"singleton_primaries\": {singletons},")?;
    writeln!(
        w,
        "    \"alts_recovered\": {},",
        s.results.iter().map(|r| r.alts.len()).sum::<usize>()
    )?;
    writeln!(
        w,
        "    \"bubbles\": {},",
        s.results.iter().map(|r| r.bubbles.len()).sum::<usize>()
    )?;
    writeln!(w, "    \"cov\": {},", quantiles(covs))?;
    writeln!(w, "    \"pident\": {}", quantiles(pidents))?;
    writeln!(w, "  }},")?;
    Ok(())
}

/// What the accession was built from: `{"source": .., "bytes": .., "sha256": .., "etag": ..,
/// "last_modified": .., "fetched_at": ..}` on one line, absent members omitted.
///
/// Always written, including `{"source": "none"}` for an accession whose bytes were never
/// read, so that "this run did not record provenance" and "there was nothing to record" are
/// different answers. That is the same rule `discards` follows.
///
/// `source` says where the bytes came from, which `inputs.source` cannot: a `cache+s3` run
/// reports `cache` for the accessions that hit and `s3` for the ones it downloaded.
fn write_payload(w: &mut impl Write, p: &Payload) -> Result<()> {
    write!(w, "      \"payload\": {{\"source\": \"{}\"", p.source.as_str())?;
    if let Some(n) = p.bytes {
        write!(w, ", \"bytes\": {n}")?;
    }
    for (key, value) in [
        ("sha256", &p.sha256),
        ("etag", &p.etag),
        ("last_modified", &p.last_modified),
        ("fetched_at", &p.fetched_at),
    ] {
        if let Some(v) = value {
            write!(w, ", \"{key}\": \"{}\"", esc(v))?;
        }
    }
    writeln!(w, "}},")?;
    Ok(())
}

/// One accession. Sub-objects go on one line each, and an empty one is omitted rather than
/// written as null, which keeps a record near the budget in `docs/summaries.md`.
fn write_result(w: &mut impl Write, cfg: &Config, r: &AccessionResult, last: bool) -> Result<()> {
    writeln!(w, "    {{")?;
    writeln!(w, "      \"acc\": \"{}\",", esc(&r.accession))?;
    writeln!(w, "      \"status\": \"{}\",", outcome_str(r.outcome))?;
    write_payload(w, &r.payload)?;
    writeln!(w, "      \"n_unitigs\": {},", r.n_unitigs)?;
    writeln!(w, "      \"n_anchors\": {},", r.n_anchors)?;
    writeln!(w, "      \"subgraph_nodes\": {},", r.subgraph_nodes)?;
    match r.subgraph_fraction {
        Some(f) => writeln!(w, "      \"subgraph_fraction\": {f:.6},")?,
        None => writeln!(w, "      \"subgraph_fraction\": null,")?,
    }
    writeln!(w, "      \"provisional\": {},", r.provisional)?;
    // WV-03. `ka_max` at the ceiling means every VAF, the prune floor and the majority calls
    // on this accession are computed against lower bounds.
    writeln!(w, "      \"ka_max\": {:.1},", r.ka_max)?;
    writeln!(w, "      \"depth_saturated\": {},", r.depth_saturated)?;
    // WV-20: null when the primary could not be mapped to the query at all, which is not
    // the same as finding no sites. A missing VCF reads as the latter.
    match r.snv_sites {
        Some(n) => writeln!(w, "      \"snv_sites\": {n},")?,
        None => writeln!(w, "      \"snv_sites\": null,")?,
    }
    if r.haplotypes.is_empty() {
        writeln!(w, "      \"reconstructed\": [],")?;
    } else {
        writeln!(w, "      \"reconstructed\": [")?;
        for (h_idx, h) in r.haplotypes.iter().enumerate() {
            let last_h = h_idx + 1 == r.haplotypes.len();
            write!(w, "        {{\"haplotype_id\": \"{}\"", esc(&h.id))?;
            write!(w, ", \"seq_sha256\": \"{}\"", seq_sha256(&h.seq))?;
            write!(w, ", \"length\": {}", h.seq.len())?;
            // 0.1.x wrote the sequence here unconditionally. Opt-in from 0.2.0.
            if cfg.summary_seqs {
                write!(w, ", \"seq\": \"{}\"", String::from_utf8_lossy(&h.seq))?;
            }
            // Four decimals, as everywhere else in the summary. 0.1.x wrote `pident` at full
            // f64 precision, 17 digits of it, which is noise in a diff and is finer than the
            // alignment it came from.
            write!(w, ", \"cov\": {:.4}", h.cov)?;
            write!(w, ", \"pident\": {:.4}", h.pident)?;
            write!(w, ", \"pident_inherited\": {}", h.pident_inherited)?;
            write!(w, ", \"nodes\": {}", h.nodes)?;
            write!(w, ", \"mean_ka\": {:.1}", h.mean_ka)?;
            write!(w, ", \"vaf\": {:.2}", h.vaf)?;
            // Emit `anchor_used` only when set. Absent field means single-anchor.
            if let Some(rid) = &h.anchor_used {
                write!(w, ", \"anchor_used\": \"{}\"", esc(rid))?;
            }
            writeln!(w, "}}{}", if last_h { "" } else { "," })?;
        }
        writeln!(w, "      ],")?;
    }
    writeln!(w, "      \"alts\": {},", r.alts.len())?;
    writeln!(w, "      \"bubbles\": {},", r.bubbles.len())?;
    // WV-07. The stage prints this and warns above 11; until 0.2.0 the summary said only
    // `pident_inherited`, which cannot tell 1 correction from 27.
    writeln!(w, "      \"majority_corrections\": {},", r.majority_corrections)?;
    // WV-02. Of those corrections, the ones where the branch the aligner walked was itself a
    // minority of its fork, which 0.1.x refused to touch at all.
    writeln!(w, "      \"minority_rescues\": {},", r.majority.minority_rescues)?;
    // WV-02: the route was still being rewritten when the backstop stopped it.
    writeln!(w, "      \"hit_pass_cap\": {},", r.majority.hit_pass_cap)?;
    // Forks still on the route where `_primary` carries the rarer allele, by why it was left.
    // Omitted when there are none, so an ordinary record does not pay for it.
    let ml = &r.majority.left;
    if ml.total() > 0 {
        writeln!(
            w,
            "      \"minority_left\": {{\"alt_below_floor\": {}, \"not_outright_majority\": {}, \
             \"exempt_minority\": {}, \"saturated\": {}, \
             \"no_rejoin\": {}, \"length_refused\": {}, \
             \"window_not_shiftable\": {}, \
             \"coverage_gate\": {}, \"identical\": {}}},",
            ml.alt_below_floor, ml.not_outright_majority, ml.exempt_minority,
            ml.saturated, ml.no_rejoin,
            ml.length_refused,
            ml.window_not_shiftable, ml.coverage_gate, ml.identical,
        )?;
    }
    let v = &r.variants;
    if v.forks > 0 || v.biallelic > 0 || v.recovered > 0 {
        writeln!(
            w,
            "      \"variants\": {{\"depth\": {:.1}, \"floor\": {:.1}, \"forks\": {}, \
             \"biallelic\": {}, \"recovered\": {}, \"beyond_cap\": {}, \
             \"indels\": {}, \"no_rejoin\": {}, \
             \"not_substitution\": {}, \"over_indel_cap\": {}, \"frameshifting\": {}, \
             \"window_not_shiftable\": {}, \"identical\": {}}},",
            v.depth, v.floor, v.forks, v.biallelic, v.recovered, v.beyond_cap,
            v.indels, v.no_rejoin,
            v.not_substitution, v.over_indel_cap, v.frameshifting, v.window_not_shiftable,
            v.identical
        )?;
    }
    // Per-reference candidate snapshot, multi-anchor only. Omitted in single-anchor runs.
    if !r.candidates.is_empty() {
        writeln!(w, "      \"candidates\": [")?;
        let n = r.candidates.len();
        for (i, c) in r.candidates.iter().enumerate() {
            let num = |x: Option<f64>| x.map_or("null".to_string(), |v| v.to_string());
            let int = |x: Option<usize>| x.map_or("null".to_string(), |v| v.to_string());
            writeln!(
                w,
                "        {{\"record_id\": \"{}\", \"pool_index\": {}, \"cov\": {}, \
                 \"pident\": {}, \"length\": {}, \"nodes\": {}, \"path_len\": {}, \
                 \"num_matches\": {}, \"aln_len\": {}, \"selected\": {}}}{}",
                esc(&c.record_id),
                c.pool_index,
                num(c.cov),
                num(c.pident),
                int(c.length),
                int(c.nodes),
                int(c.path_len),
                int(c.num_matches),
                int(c.aln_len),
                c.selected,
                if i + 1 == n { "" } else { "," }
            )?;
        }
        writeln!(w, "      ],")?;
    }
    if let Some(p) = &r.prune {
        let span = |s: Option<(usize, usize)>| match s {
            Some((a, b)) => format!("[{a}, {b}]"),
            None => "null".to_string(),
        };
        writeln!(
            w,
            "      \"prune\": {{\"outcome\": \"{}\", \"depth\": {:.1}, \"floor\": {:.1}, \
             \"unitigs\": {}, \"unitigs_kept\": {}, \"unpruned_span\": {}, \
             \"pruned_span\": {}}},",
            p.result.as_str(),
            p.depth,
            p.floor,
            p.unitigs,
            p.unitigs_kept,
            span(Some(p.unpruned_span)),
            span(p.pruned_span)
        )?;
    }
    if let Some(e) = &r.error {
        writeln!(w, "      \"error\": \"{}\",", esc(e))?;
    }
    let t = &r.timing;
    writeln!(
        w,
        "      \"timing\": {{\"stream\": {:.3}, \"parse\": {:.3}, \"anchor\": {:.3}, \
         \"extract\": {:.3}, \"align\": {:.3}, \"stitch\": {:.3}, \"phase\": {:.3}, \
         \"total\": {:.3}}}",
        t.stream, t.parse, t.anchor, t.extract, t.align, t.stitch, t.phase, t.total
    )?;
    writeln!(w, "    }}{}", if last { "" } else { "," })?;
    Ok(())
}

/// The summary contract in `docs/summaries.md`, as the walker implements it.
///
/// Seven blocks in a fixed order: `schema`, `tool`, `inputs`, `settings`, `run`, `discards`,
/// `cohort`, `records`. Within a schema version, fields are only ever added.
pub fn write_summary(summary: &Summary, run: &RunInfo, path: &Path) -> Result<()> {
    let mut w = create(path)?;
    writeln!(w, "{{")?;
    writeln!(w, "  \"schema\": 1,")?;
    write_tool(&mut w, run)?;
    write_inputs(&mut w, run)?;
    write_settings(&mut w, run)?;
    write_run_block(&mut w, summary, run)?;
    write_discards(&mut w, summary, run)?;
    write_cohort(&mut w, summary)?;
    writeln!(w, "  \"records\": [")?;
    for (i, r) in summary.results.iter().enumerate() {
        write_result(&mut w, run.cfg, r, i + 1 == summary.results.len())?;
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
pub fn snv_sites_all(summary: &mut Summary, queries: &[Query]) -> Vec<Vec<VcfRow>> {
    let mut all = Vec::with_capacity(summary.results.len());
    for r in summary.results.iter_mut() {
        let (rows, mapped) = snv_sites_counted(r, reference_of(r, queries));
        // WV-20: `Some(0)` is "looked, found nothing"; `None` is "could not look". The
        // summary could not tell them apart, and a missing VCF reads as the former.
        r.snv_sites = mapped.then_some(rows.len());
        all.push(rows);
    }
    all
}

/// The reference an accession's primary was reconstructed against: in multi-anchor mode the
/// winning one, named by the primary's `anchor_used`; otherwise the only one. Its coordinates
/// are the ones the VCF and the long table must use. Writing every accession against the first
/// reference of the pool put a SNV at `gene:601` down as `othergene 438`.
fn reference_of<'a>(r: &AccessionResult, queries: &'a [Query]) -> &'a Query {
    r.anchor_used
        .as_deref()
        .or_else(|| r.haplotype.as_ref().and_then(|h| h.anchor_used.as_deref()))
        .and_then(|rid| queries.iter().find(|q| q.record_id == rid))
        .unwrap_or(&queries[0])
}

pub fn write_vcfs(
    summary: &Summary,
    queries: &[Query],
    out_dir: &Path,
    sites: &[Vec<VcfRow>],
) -> Result<()> {
    std::fs::create_dir_all(out_dir)?;

    // The filename is the accession, which `docs/vcf_example.md` documents and the committed
    // examples rely on. `read_accession_list` drops duplicates, so each path is written once.
    for (r, rows) in summary.results.iter().zip(sites) {
        let path = out_dir.join(format!("{}.vcf", r.accession));
        // WV-16: never touch a file we did not write. The same check guards the removal and
        // the overwrite, because both destroy it.
        if !walker_owned(&path) {
            eprintln!(
                "[warn] refusing to replace {}: it was not written by logan-walker. Use \
                 another --out-dir; this accession's VCF was not written.",
                path.display()
            );
            continue;
        }
        if rows.is_empty() {
            // A VCF an earlier run left in this --out-dir would claim sites this run did not
            // find, beside a primary that no longer has them.
            // WV-19: a stale file we cannot remove is a warning, not a dead run. This used
            // to be fatal, so one EACCES or EISDIR here threw away the summary, the SNV table
            // and the alt FASTA for every accession that had already succeeded.
            if let Err(e) = std::fs::remove_file(&path) {
                if e.kind() != std::io::ErrorKind::NotFound {
                    eprintln!("[warn] cannot remove stale {}: {e}", path.display());
                }
            }
            continue;
        }

        let query = reference_of(r, queries);
        let chrom = chrom_of(query);
        let mut af = AtomicFile::create(&path)?;
        let w = af.w();
        write_vcf_header(w, chrom, query.seq.len(), &r.accession)?;
        for row in rows {
            let alt_bases: Vec<String> =
                row.alts.iter().map(|&(b, _, _)| (b as char).to_string()).collect();
            let afs: Vec<String> = row.alts.iter().map(|&(_, af, _)| format!("{af:.4}")).collect();
            let aos: Vec<String> = row.alts.iter().map(|&(_, _, ka)| format!("{ka:.1}")).collect();
            writeln!(
                w,
                "{}\t{}\t.\t{}\t{}\t.\t{}\tAF={};AO={};RO={:.1};FORK={}{}\tGT:AF:AO:RO\t.:{}:{}:{:.1}",
                chrom,
                row.pos,
                row.ref_base as char,
                alt_bases.join(","),
                // Column 7 is FILTER, after REF and ALT. Getting this one slot wrong wrote
                // the filter into REF, which the CLI test caught.
                if row.ka_saturated { "KaSat" } else { "PASS" },
                afs.join(","),
                aos.join(","),
                row.ref_ka,
                row.fork,
                // WV-08: only multi-site alleles carry one, so single-site rows are
                // byte-identical to what 0.1.x wrote.
                row.block.map_or(String::new(), |b| format!(";BLOCK={b}")),
                afs.join(","),
                aos.join(","),
                row.ref_ka,
            )?;
        }
        af.commit()?;
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
/// Columns: `gene pos ref alt accession af ka_alt ka_ref fork ka_saturated block`. Frequencies are `ka:f` ratios,
/// the same values the VCF carries.
pub fn write_snv_table(
    summary: &Summary,
    queries: &[Query],
    path: &Path,
    sites: &[Vec<VcfRow>],
) -> Result<()> {
    let mut w = create(path)?;
    writeln!(
        w,
        "gene\tpos\tref\talt\taccession\taf\tka_alt\tka_ref\tfork\tka_saturated\tblock"
    )?;
    for (r, rows) in summary.results.iter().zip(sites) {
        let chrom = chrom_of(reference_of(r, queries));
        for row in rows {
            for &(alt, af, ka_alt) in &row.alts {
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.1}\t{:.1}\t{}\t{}\t{}",
                    chrom,
                    row.pos,
                    row.ref_base as char,
                    alt as char,
                    r.accession,
                    af,
                    ka_alt,
                    row.ref_ka,
                    row.fork,
                    // WV-03: `af` here is a lower-bound ratio. A cohort query that averages
                    // these without filtering on this column averages overstated frequencies.
                    u8::from(row.ka_saturated),
                    // WV-08: rows sharing a block are on one molecule. Empty for a
                    // single-site call, so a `GROUP BY block` ignores them naturally.
                    row.block.map_or(String::new(), |b| b.to_string()),
                )?;
            }
        }
    }
    w.flush()?;
    Ok(())
}

/// A reference base as VCF 4.2 allows it to appear in REF: A, C, G, T, or N.
///
/// WV-14. The base was copied verbatim, so an IUPAC ambiguity code in the reference produced
/// `REF=R` and an invalid record.
fn vcf_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        n @ (b'A' | b'C' | b'G' | b'T') => n,
        _ => b'N',
    }
}

/// Whether a reference record id can be a VCF contig name.
///
/// WV-14. CHROM and `##contig=<ID=...>` take the reference's id verbatim, which is the first
/// token of its FASTA header or the file stem. A `,` or `<` in it produces a malformed header
/// line and a record no parser will read. The grammar is VCF 4.3's.
pub fn valid_contig_id(id: &str) -> bool {
    fn ok(c: char, first: bool) -> bool {
        c.is_ascii_alphanumeric()
            || "!#$%&+./:;?@^_|~-".contains(c)
            || (!first && (c == '*' || c == '='))
    }
    let mut chars = id.chars();
    match chars.next() {
        None => false,
        Some(c) if !ok(c, true) => false,
        Some(_) => chars.all(|c| ok(c, false)),
    }
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
    snv_sites_counted(r, query).0
}

/// As `snv_sites`, also saying whether the primary could be mapped to the query at all.
///
/// WV-20: `false` means the mapping was refused (too divergent for the budget), which is a
/// different statement from "no sites" and must not reach `summary.json` as `0`.
fn snv_sites_counted(r: &AccessionResult, query: &Query) -> (Vec<VcfRow>, bool) {
    let Some(primary) = &r.haplotype else { return (Vec::new(), true) };
    if r.alts.is_empty() {
        return (Vec::new(), true);
    }
    // Map each primary base to a query position (None where the primary base aligns to a gap in
    // the query, i.e. an insertion relative to the gene).
    //
    // WV-09: bounded by memory rather than by a 20 kb length test. That test dated from the
    // full O(n*m) DP and survived the move to a banded one, so it dropped every SNV and the
    // whole VCF for any gene over 20 kb (SARS-CoV-2 ORF1ab is 21,290 nt) while the run still
    // exited 0 and the docs say a missing VCF means no sites.
    let Some(map) = align_within(&primary.seq, &query.seq, SNV_ALIGN_BUDGET) else {
        eprintln!(
            "[{}] SNV sites skipped: primary {} nt against query {} nt cannot be mapped \
             within {} MiB; the pair is too divergent for a narrow band. The VCF and the SNV \
             table are empty for this accession, which is NOT the same as having no variants.",
            r.accession,
            primary.seq.len(),
            query.seq.len(),
            SNV_ALIGN_BUDGET >> 20,
        );
        return (Vec::new(), false);
    };

    let mut rows: Vec<VcfRow> = Vec::new();
    let mut block_id = 0usize;
    let mut skipped_blocks = 0usize;
    for a in &r.alts {
        let alt_seq = &a.haplotype.seq;
        if alt_seq.len() != primary.seq.len() {
            continue; // an indel allele; it stays in alt_haplotypes.fasta
        }
        // WV-08. This used to stop at the SECOND differing base and discard the whole allele.
        // In a k=31 graph two substitutions within 30 nt sit on one bubble, so `find_alts`
        // returns them as a single alternative and BOTH sites vanished from the VCF and
        // snvs.tsv with no message. Two changes in one codon are the common case. Measured on
        // the panel, 58 of 267 substitution alleles span more than 3 nt, so per-site counts
        // and cohort matrices were biased against exactly the dense variable regions they are
        // built to find.
        let diffs: Vec<usize> = primary
            .seq
            .iter()
            .zip(alt_seq.iter())
            .enumerate()
            .filter(|(_, (x, y))| x != y)
            .map(|(i, _)| i)
            .collect();
        if diffs.is_empty() || diffs.len() > MAX_LINKED_SNVS {
            // Past the cap this is a rerouted block rather than a few linked substitutions,
            // and calling every base of it as a SNV would be an assertion the graph does not
            // support.
            if diffs.len() > MAX_LINKED_SNVS {
                skipped_blocks += 1;
            }
            continue;
        }
        // Rows from one allele are on one molecule, unlike rows that merely share a position.
        // Only a multi-site allele needs saying so; a single-site one keeps the old output.
        let block = if diffs.len() > 1 {
            block_id += 1;
            Some(block_id)
        } else {
            None
        };
        for &i in &diffs {
        let Some(qpos) = map[i] else { continue }; // primary base has no query coordinate

        // WV-14: REF is copied from the reference, and VCF 4.2 allows only A/C/G/T/N there.
        // An IUPAC code in the query (R, Y, S, ...) produced `REF=R`, which is not a legal
        // record. The ambiguity is in the reference, not in what we observed, so it maps to N.
        let ref_base = vcf_base(query.seq[qpos]);
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
        rows.push(VcfRow {
            pos: qpos + 1,
            ref_base,
            ka_saturated: a.ka_saturated,
            block,
            alts,
            ref_ka,
            fork: a.fork,
        });
        }
    }
    if skipped_blocks > 0 {
        eprintln!(
            "[{}] {} alternative allele(s) carry more than {} linked substitutions and are \
             not called as SNVs; they remain in alt_haplotypes.fasta",
            r.accession, skipped_blocks, MAX_LINKED_SNVS
        );
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
                // OR the flag: if any fork contributing to this record is clipped, every
                // frequency on the record is drawn from a lower-bound total.
                last.ka_saturated |= row.ka_saturated;
                last.block = last.block.or(row.block);
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
    (merged, true)
}

/// Most linked substitutions one alternative may carry and still be called as SNVs.
///
/// WV-08. Below this, an alternative with several substitutions is a short haplotype block
/// worth reporting site by site. Above it, the alternative is a rerouted stretch of graph and
/// calling every base of it as a SNV would assert more than the graph supports; those stay in
/// `alt_haplotypes.fasta` and are counted on stderr.
const MAX_LINKED_SNVS: usize = 8;

pub struct VcfRow {
    pos: usize,
    ref_base: u8,
    /// Set when this row is one of several substitutions carried by a single alternative
    /// haplotype (WV-08). Rows sharing a block ARE on one molecule, unlike two rows that
    /// merely share a position, so the `##source` caveat about co-occurrence does not apply
    /// to them. `None` for an ordinary single-site allele.
    block: Option<usize>,
    /// A branch of this site's fork is at `KA_CEILING`, so AF and RO are lower-bound ratios
    /// and the row is FILTERed `KaSat` rather than `PASS` (WV-03).
    ka_saturated: bool,
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
         and NOT phased: two rows do not assert co-occurrence on one molecule, UNLESS they \
         share a BLOCK, which means one alternative haplotype carries both. ka:f is clipped \
         at 25000, and a row whose fork reaches it is FILTERed KaSat with an overstated AF."
    )?;
    writeln!(w, "##contig=<ID={chrom},length={len}>")?;
    writeln!(w, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency: this allele ka:f over the total of the fork branches that return to the route\">")?;
    writeln!(w, "##INFO=<ID=AO,Number=A,Type=Float,Description=\"Alternate allele k-mer coverage (ka:f)\">")?;
    writeln!(w, "##INFO=<ID=RO,Number=1,Type=Float,Description=\"Reference allele k-mer coverage (ka:f); 0 if the gene base is not observed\">")?;
    writeln!(w, "##INFO=<ID=FORK,Number=1,Type=Integer,Description=\"Unitig id of the fork on the walked path\">")?;
    writeln!(
        w,
        "##INFO=<ID=BLOCK,Number=1,Type=Integer,Description=\"Rows sharing a BLOCK are \
         substitutions carried by ONE alternative haplotype, so they do co-occur on a single \
         molecule. Absent on an ordinary single-site call\">"
    )?;
    writeln!(w, "##FILTER=<ID=PASS,Description=\"All filters passed\">")?;
    writeln!(
        w,
        "##FILTER=<ID=KaSat,Description=\"A branch of this fork is at Logan's ka:f ceiling of \
         25000, so the fork total is a lower bound and AF is overstated. The true depth is not \
         recoverable from the unitig file\">"
    )?;
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
/// A band is only a speedup if the optimum stays inside it, and **the edge check alone does not
/// establish that** (WV-10). When the true optimum leaves the band, the best in-band path can
/// be strictly worse and still never touch an edge, so the traceback accepts it. Reproduced
/// with compensating indels of +33 and -33 nt in a 2.4 kb gene: the banded score was 1367
/// against 2233 for the full alignment, about 940 bases mapped to the wrong reference
/// positions, and an SNV whose true position was 761 was written at 809, with the REF base of
/// whatever sat there.
///
/// So acceptance now needs a score certificate as well. Any path that leaves the band must
/// cross to diagonal `d_max + 1` or `d_min - 1` and come back, which costs at least `g` gaps;
/// a path with `g` gaps scores at most `(n + m - g) / 2 - 2g`, every remaining column being a
/// match. If the in-band optimum is at least that, no out-of-band path can beat it and the
/// result is globally optimal. Otherwise the band doubles and the alignment is redone. At full
/// width this calls `align_full`, so the worst case is the old behaviour.
fn align_primary_to_query(primary: &[u8], query: &[u8]) -> Vec<Option<usize>> {
    align_within(primary, query, usize::MAX).expect("an unbounded budget cannot refuse")
}

/// The traceback matrix this is willing to allocate for one accession's SNV mapping.
///
/// WV-09. What this replaces was a flat 20 kb length test, left over from the full O(n*m) DP
/// and never revisited when the alignment was banded. It refused cheap inputs (a 21 kb gene at
/// high identity costs a few MB banded) while doing nothing about the expensive case, a short
/// but divergent pair that widens the band to full width. A budget measures the thing that
/// actually hurts.
const SNV_ALIGN_BUDGET: usize = 256 << 20;

/// As `align_primary_to_query`, refusing rather than allocating past `budget`.
///
/// `None` means this pair cannot be mapped inside the budget: it is too divergent for a narrow
/// band and too large for a wide one. The caller reports that; it is not the same as "no
/// variants", which is what the old guard's silence looked like.
fn align_within(primary: &[u8], query: &[u8], budget: usize) -> Option<Vec<Option<usize>>> {
    let (n, m) = (primary.len(), query.len());
    let full = n.max(m);
    // The band has to span the length difference before it can absorb any drift on top.
    let mut band = (n.abs_diff(m) + 16).max(32);
    while band < full {
        // The width `align_banded` will allocate for this band, matching its own d_min/d_max.
        let w = n.abs_diff(m) + 2 * band + 1;
        if (n + 1).saturating_mul(w) > budget {
            return None;
        }
        if let Some(map) = align_banded(primary, query, band) {
            return Some(map);
        }
        band = band.saturating_mul(2);
    }
    if (n + 1).saturating_mul(m + 1) > budget {
        return None;
    }
    Some(align_full(primary, query))
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

    // WV-10. `prev` holds row `n` after the final swap, so this is the in-band optimum.
    //
    // The cheapest escape from the band: cross to `d_max + 1` and return to the end diagonal,
    // or cross to `d_min - 1` and return. Each step off the diagonal is one gap.
    let best = prev[col(n, m)] as i64;
    let g_hi = (d_max + 1) + (d_max + 1 - delta);
    let g_lo = (1 - d_min) + (delta - d_min + 1);
    let g = g_hi.min(g_lo) as i64;
    // A path with at least `g` gaps has at most `(n + m - g) / 2` aligned columns, so even
    // with every one of them a match it cannot beat this. If the in-band optimum clears the
    // bound, nothing outside the band can beat it.
    if best < (n as i64 + m as i64 - g) / 2 - 2 * g {
        return None; // the optimum may be outside; widen
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

    /// WV-09. A long gene at high identity is cheap to align banded, and used to be refused
    /// by a flat 20 kb length test left over from the full O(n*m) DP. SARS-CoV-2 ORF1ab is
    /// 21,290 nt, so every SNV and the whole VCF were dropped for it while the run exited 0.
    #[test]
    fn a_long_gene_at_high_identity_is_mapped_not_skipped() {
        let mut seed = 0xA11CEu64;
        let a = seq(21_291, &mut seed);
        let mut b = a.clone();
        b[10_000] = if b[10_000] == b'A' { b'C' } else { b'A' };
        let map = super::align_within(&a, &b, super::SNV_ALIGN_BUDGET)
            .expect("21 kb at high identity is well inside the budget");
        assert_eq!(map.len(), a.len());
        assert_eq!(map[10_000], Some(10_000), "a substitution keeps coordinates aligned");
    }

    /// The other half: the budget refuses rather than allocating. A short but wildly
    /// divergent pair is what widens the band to full width, which the old length test let
    /// through untouched.
    #[test]
    fn a_pass_over_the_budget_is_refused_not_allocated() {
        let a = vec![b'A'; 4_000];
        let b = vec![b'T'; 4_000];
        // 4001 * 4001 bytes is about 16 MB, so a 1 MB budget must refuse.
        assert!(super::align_within(&a, &b, 1 << 20).is_none());
        // And with room, the same pair still produces an answer.
        assert!(super::align_within(&a, &b, super::SNV_ALIGN_BUDGET).is_some());
    }

    /// WV-10. Compensating indels of opposite sign, each longer than the initial band, are
    /// the shape the edge check alone could not catch: the path leaves the band and comes
    /// back, so the best in-band route never touches an edge and was accepted although it was
    /// strictly worse. Before the score certificate, +33/-33 scored 1367 against the exact
    /// 2233, and about 940 bases mapped to the wrong reference positions.
    ///
    /// The initial band is `max(|n - m| + 16, 32)`, and these keep `n == m`, so 33 is already
    /// past it and 40 and 60 are well past it.
    #[test]
    fn compensating_indels_longer_than_the_band_still_match_the_exact_dp() {
        let mut seed = 0x5EEDu64;
        for len in [33usize, 40, 60] {
            let a = seq(2400, &mut seed);
            let mut b = a.clone();
            // An insertion early and a deletion of the same size 900 bases later, so the
            // lengths match and the path has to leave the diagonal and return.
            let ins = seq(len, &mut seed);
            b.splice(300..300, ins);
            b.drain(1200..1200 + len);
            assert_eq!(a.len(), b.len(), "the two indels must compensate");
            assert_eq!(
                align_primary_to_query(&a, &b),
                align_full(&a, &b),
                "compensating +{len}/-{len} disagreed with the exact DP"
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
    use super::{align_primary_to_query, chrom_of, reference_of, snv_sites, snv_sites_all, write_vcfs};
    use crate::pipeline::Summary;
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
            ka_saturated: false,
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
            payload: crate::types::Payload::default(),
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
            candidates: Vec::new(),
            snv_sites: None,
            ka_max: 0.0,
            depth_saturated: false,
            majority_corrections: 0,
            majority: crate::types::MajorityReport::default(),
            variants: crate::types::VariantStats::default(),
            subgraph_fraction: None,
            provisional: false,
            prune: None,
        }
    }

    fn query(seq: &[u8]) -> Query {
        Query { name: "g".into(), record_id: "g".into(), seq: seq.to_vec(), path: PathBuf::from("g.fa"), n_records: 1 }
    }

    /// Multi-anchor: the VCF and the long table are written in the winning reference's
    /// coordinates. They used to be written against the first reference of the pool, so a SNV
    /// at `gene:601` came out as `othergene 438`.
    #[test]
    fn variant_writers_use_the_winning_reference() {
        let first = Query { name: "a".into(), record_id: "first".into(), seq: b"AAAA".to_vec(),
                            path: PathBuf::from("a.fa"), n_records: 1 };
        let winner = Query { name: "b".into(), record_id: "winner".into(), seq: b"CCCCCC".to_vec(),
                             path: PathBuf::from("b.fa"), n_records: 1 };
        let pool = [first, winner];
        let mut r = result(b"CCCCCC", Vec::new());
        assert_eq!(reference_of(&r, &pool).record_id, "first", "single-anchor falls back");
        r.anchor_used = Some("winner".into());
        assert_eq!(reference_of(&r, &pool).record_id, "winner");
        assert_eq!(chrom_of(reference_of(&r, &pool)), "winner");
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

    /// WV-08. Two substitutions on one alternative yield two rows, sharing a BLOCK.
    ///
    /// This used to stop at the second differing base and discard the whole allele, so both
    /// sites vanished from the VCF and snvs.tsv with no message. In a k=31 graph two
    /// substitutions within 30 nt share one bubble and arrive as a single alternative, and
    /// two changes in one codon are the common case, so the loss concentrated in exactly the
    /// dense variable regions these outputs exist to find.
    #[test]
    fn two_linked_substitutions_yield_two_rows_sharing_a_block() {
        let primary = b"TTTATTTATTT";
        let q = query(b"TTTATTTATTT");
        // One alternative differing at index 3 and index 7, four bases apart.
        let rows = snv_sites(&result(primary, vec![alt(b"TTTCTTTGTTT", 0.25, 25.0)]), &q);
        assert_eq!(rows.len(), 2, "both substitutions must be called");
        assert_eq!(rows[0].pos, 4);
        assert_eq!(rows[1].pos, 8);
        // Same block: these two ARE on one molecule, unlike two rows that share a position.
        assert_eq!(rows[0].block, rows[1].block);
        assert!(rows[0].block.is_some());
        assert_eq!(rows[0].alts[0].0, b'C');
        assert_eq!(rows[1].alts[0].0, b'G');
    }

    /// A single-site allele carries no block, so its row is what 0.1.x wrote.
    #[test]
    fn a_single_substitution_carries_no_block() {
        let q = query(b"TTTATTT");
        let rows = snv_sites(&result(b"TTTATTT", vec![alt(b"TTTCTTT", 0.25, 25.0)]), &q);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].block, None);
    }

    /// Past the cap the alternative is a rerouted stretch of graph, not a few linked
    /// substitutions, and calling every base of it would assert more than the graph supports.
    #[test]
    fn an_allele_past_the_linked_cap_is_not_called() {
        let primary = vec![b'T'; 40];
        let mut alt_seq = primary.clone();
        for i in 0..super::MAX_LINKED_SNVS + 1 {
            alt_seq[i * 2] = b'A';
        }
        let q = query(&primary);
        let rows = snv_sites(&result(&primary, vec![alt(&alt_seq, 0.25, 25.0)]), &q);
        assert!(rows.is_empty(), "{} diffs is past the cap", super::MAX_LINKED_SNVS + 1);
    }

    /// WV-13. ALT is every sample allele that differs from the gene base, which is the
    /// MAJORITY allele when the sample's consensus differs from the reference. The docs used
    /// to call ALT "the minor allele" flatly, so a cohort query averaging AF as a
    /// minor-allele frequency silently mixed in values above 0.5.
    #[test]
    fn alt_is_the_majority_allele_when_the_consensus_differs_from_the_reference() {
        // Primary (the sample's majority) carries C at index 3; the gene base there is A.
        let primary = b"TTTCTTT";
        let q = query(b"TTTATTT");
        // The minority allele happens to match the gene base, so only the majority is an ALT.
        let rows = snv_sites(&result(primary, vec![alt(b"TTTATTT", 0.4, 40.0)]), &q);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].ref_base, b'A', "REF is the gene base");
        let alts: Vec<(u8, f64)> = rows[0].alts.iter().map(|&(b, af, _)| (b, af)).collect();
        assert_eq!(alts, vec![(b'C', 0.60)], "ALT is the sample majority, AF above 0.5");
        assert_eq!(rows[0].ref_ka, 40.0, "RO is the gene base's own coverage");
    }

    /// WV-24. One identity formula, used by `stitch` and `pipeline`, including the
    /// `aln_len == 0` case that would otherwise be NaN and make the summary invalid (WV-15).
    #[test]
    fn pident_is_matches_over_aligned_length() {
        use crate::types::GafPath;
        let p = |num_matches, aln_len| GafPath {
            steps: vec![],
            query_start: 0,
            query_end: 0,
            path_len: 0,
            path_start: 0,
            path_end: 0,
            num_matches,
            aln_len,
        };
        assert_eq!(p(95, 100).pident(), 95.0);
        assert_eq!(p(100, 100).pident(), 100.0);
        assert_eq!(p(0, 0).pident(), 0.0, "no alignment is 0.0, never NaN");
        assert!(p(0, 0).pident().is_finite());
    }

    /// WV-14. The reference id becomes CHROM and `##contig=<ID=...>` verbatim, so an id
    /// carrying `,` or `<` writes a malformed header and an unreadable record.
    #[test]
    fn contig_ids_follow_the_vcf_grammar() {
        for good in ["gene", "hiv1_pol", "NC_045512.2", "chr1", "a|b", "x=y", "s*t"] {
            assert!(super::valid_contig_id(good), "{good:?} must be accepted");
        }
        for bad in ["", "gene,with", "a<b>c", "has space", "=leading", "*leading", "tab\there"] {
            assert!(!super::valid_contig_id(bad), "{bad:?} must be refused");
        }
    }

    /// WV-14. VCF 4.2 allows only A/C/G/T/N in REF, and the base was copied verbatim, so an
    /// IUPAC code in the reference produced `REF=R`.
    #[test]
    fn an_ambiguous_reference_base_is_written_as_n() {
        // Gene base R (A or G) at index 3; the sample carries C there.
        let q = query(b"TTTRTTT");
        let rows = snv_sites(&result(b"TTTATTT", vec![alt(b"TTTCTTT", 0.25, 25.0)]), &q);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].ref_base, b'N', "R is not legal in REF");
        // Both sample alleles differ from N, so both are reported, which is the honest answer.
        let bases: Vec<u8> = rows[0].alts.iter().map(|&(b, _, _)| b).collect();
        assert_eq!(bases, vec![b'C', b'A']);
    }

    /// WV-16. A VCF this walker did not write is never removed or overwritten.
    ///
    /// Point `-o` at a directory holding per-sample VCFs from a read-based caller, which use
    /// the same `<ACC>.vcf` names, and a run with no sites used to delete them and exit 0.
    #[test]
    fn a_foreign_vcf_is_left_alone() {
        let dir = std::env::temp_dir().join(format!("logan-vcf-foreign-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let foreign = dir.join("ACC.vcf");
        let body = "##fileformat=VCFv4.2\n##source=bcftools\nnot ours\n";
        std::fs::write(&foreign, body).unwrap();

        let q = [query(b"TTTATTT")];
        // No sites, which is the path that used to remove the file.
        let mut s = Summary { total: 1, success: 1, results: vec![result(b"TTTATTT", vec![])] };
        let sites = snv_sites_all(&mut s, &q);
        write_vcfs(&s, &q, &dir, &sites).unwrap();
        assert_eq!(std::fs::read_to_string(&foreign).unwrap(), body, "must be untouched");

        // And with sites, the overwrite is refused rather than silently replacing it.
        let mut s2 = Summary {
            total: 1,
            success: 1,
            results: vec![result(b"TTTATTT", vec![alt(b"TTTCTTT", 0.25, 25.0)])],
        };
        let sites2 = snv_sites_all(&mut s2, &q);
        write_vcfs(&s2, &q, &dir, &sites2).unwrap();
        assert_eq!(std::fs::read_to_string(&foreign).unwrap(), body, "still untouched");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// WV-18. Nothing appears at the final path until the write is complete, and an
    /// abandoned write leaves no temporary file behind.
    #[test]
    fn an_abandoned_atomic_write_leaves_nothing() {
        let dir = std::env::temp_dir().join(format!("logan-atomic-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let dst = dir.join("out.txt");
        {
            let mut af = super::AtomicFile::create(&dst).unwrap();
            use std::io::Write;
            write!(af.w(), "half a file").unwrap();
            // Dropped without commit, as a `?` on a later write would do.
        }
        assert!(!dst.exists(), "an uncommitted write must not appear at the final path");
        let left: Vec<_> = std::fs::read_dir(&dir).unwrap().filter_map(|e| e.ok()).collect();
        assert!(left.is_empty(), "no temporary file may be left: {:?}", left.len());

        // And a committed one does appear, in full.
        let mut af = super::AtomicFile::create(&dst).unwrap();
        use std::io::Write;
        write!(af.w(), "whole file").unwrap();
        af.commit().unwrap();
        assert_eq!(std::fs::read_to_string(&dst).unwrap(), "whole file");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// WV-03. A site whose fork touches Logan's ceiling is FILTERed, not PASSed, because its
    /// AF is a ratio against a total that is a lower bound.
    ///
    /// The number that makes this matter: a 3,000x minor allele beside a true 250,000x major
    /// allele is 1.2%, but measured against a clipped 25,000 it reports 10.7% and clears the
    /// default `--min-vaf 0.10`. The true depth cannot be recovered from the unitig file, so
    /// the row is published with a flag rather than dropped or silently corrected.
    #[test]
    fn a_saturated_fork_is_flagged_not_passed() {
        let dir = std::env::temp_dir().join(format!("logan-vcf-kasat-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let q = [query(b"TTTATTT")];
        let mut a = alt(b"TTTCTTT", 0.1071, 3000.0);
        a.ka_saturated = true;
        let mut s = Summary { total: 1, success: 1, results: vec![result(b"TTTATTT", vec![a])] };
        let sites = snv_sites_all(&mut s, &q);
        write_vcfs(&s, &q, &dir, &sites).unwrap();
        let vcf = std::fs::read_to_string(dir.join("ACC.vcf")).unwrap();
        let _ = std::fs::remove_dir_all(&dir);
        assert!(vcf.contains("##FILTER=<ID=KaSat,"), "the filter must be declared: {vcf}");
        let row = vcf.lines().find(|l| !l.starts_with('#')).expect("a data row");
        let cols: Vec<&str> = row.split('\t').collect();
        assert_eq!(cols[3], "A", "REF");
        assert_eq!(cols[4], "C", "ALT");
        assert_eq!(cols[6], "KaSat", "FILTER, which is column 7: {row}");
        // The frequency is still published; it is the flag that says not to trust it.
        assert!(cols[7].contains("AF=0.1071"), "{row}");
    }

    /// The unsaturated case stays byte-identical to what 0.1.x wrote, so an ordinary cohort
    /// does not suddenly fail `bcftools view -f PASS`.
    #[test]
    fn an_unsaturated_fork_still_passes() {
        let dir = std::env::temp_dir().join(format!("logan-vcf-pass-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let q = [query(b"TTTATTT")];
        let mut s = Summary {
            total: 1,
            success: 1,
            results: vec![result(b"TTTATTT", vec![alt(b"TTTCTTT", 0.25, 25.0)])],
        };
        let sites = snv_sites_all(&mut s, &q);
        write_vcfs(&s, &q, &dir, &sites).unwrap();
        let vcf = std::fs::read_to_string(dir.join("ACC.vcf")).unwrap();
        let _ = std::fs::remove_dir_all(&dir);
        let row = vcf.lines().find(|l| !l.starts_with('#')).expect("a data row");
        assert_eq!(row.split('\t').nth(6).unwrap(), "PASS", "{row}");
    }

    /// A rerun into the same `--out-dir` that finds no sites for an accession removes the VCF
    /// an earlier run wrote for it, rather than leaving it beside a primary without those sites.
    #[test]
    fn a_rerun_without_sites_removes_the_old_vcf() {
        let dir = std::env::temp_dir().join(format!("logan-vcf-stale-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let q = [query(b"TTTATTT")];
        let run = |alts: Vec<AltHaplotype>| {
            let mut s = Summary { total: 1, success: 1, results: vec![result(b"TTTATTT", alts)] };
            let sites = snv_sites_all(&mut s, &q);
            write_vcfs(&s, &q, &dir, &sites).unwrap();
        };
        run(vec![alt(b"TTTGTTT", 0.25, 25.0)]);
        assert!(dir.join("ACC.vcf").exists(), "the first run found a site");
        run(Vec::new());
        assert!(!dir.join("ACC.vcf").exists(), "the second run found none");
        run(Vec::new()); // and no file to remove is not an error
        std::fs::remove_dir_all(&dir).ok();
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
                payload: crate::types::Payload::default(),
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
                candidates: Vec::new(),
                snv_sites: None,
            ka_max: 0.0,
            depth_saturated: false,
            majority_corrections: 0,
                majority: crate::types::MajorityReport::default(),
                variants: crate::types::VariantStats::default(),
                subgraph_fraction: None,
                provisional: false,
                prune: None,
            }],
        }
    }

    /// `write_summary` with a fixed provenance block.
    fn write_summary(s: &Summary, p: &std::path::Path) -> crate::error::Result<()> {
        let cfg = crate::config::Config::default();
        let queries = [crate::types::Query {
            name: "gene".into(), record_id: "gene".into(), seq: b"ACGT".to_vec(),
            path: std::path::PathBuf::from("gene.fa"), n_records: 1,
        }];
        // A fixed clock so the volatile fields (`run.started`, the timings) are stable and
        // the rest of the file is byte-comparable, the property `docs/summaries.md` promises.
        let run = super::RunInfo {
            graphaligner: Some((
                "GraphAligner bioconda 1.0.17b-0",
                "/opt/conda/bin/GraphAligner",
            )),
            cfg: &cfg,
            queries: &queries,
            started: std::time::UNIX_EPOCH + std::time::Duration::from_secs(1_763_200_000),
            wall_seconds: 1.5,
            n_requested: 1,
        };
        super::write_summary(s, &run, p)
    }

    /// The versions and settings, and the per-accession guard and prune decisions that used
    /// to reach only stderr.
    #[test]
    fn summary_json_records_provenance_and_decisions() {
        use crate::types::{PruneReport, PruneResult};
        let mut s = summary(hap("ACC_primary", false), vec![]);
        s.results[0].subgraph_fraction = Some(0.161189);
        s.results[0].provisional = true;
        s.results[0].payload = crate::types::Payload {
            source: crate::types::PayloadSource::S3,
            bytes: Some(1_861_058),
            sha256: Some("a".repeat(64)),
            etag: Some("d2dd288b238c3c88b46acdcd78860dce".into()),
            last_modified: Some("2024-03-31T17:38:36Z".into()),
            fetched_at: Some("2026-11-16T16:47:50Z".into()),
        };
        s.results[0].prune = Some(PruneReport {
            result: PruneResult::RefusedShorterSpan, depth: 1700.0, floor: 17.0,
            unitigs: 4976, unitigs_kept: 1734, unpruned_span: (0, 3821),
            pruned_span: Some((424, 3821)),
        });
        let p = tmp("summary_provenance.json");
        write_summary(&s, &p).unwrap();
        let out = std::fs::read_to_string(&p).unwrap();
        std::fs::remove_file(&p).ok();
        for want in [
            // The contract's block order and the provenance in them.
            "\"schema\": 1",
            "\"name\": \"logan-walker\"",
            concat!("\"version\": \"", env!("CARGO_PKG_VERSION"), "\""),
            "\"name\": \"GraphAligner\", \"version\": \"GraphAligner bioconda 1.0.17b-0\", \
             \"path\": \"/opt/conda/bin/GraphAligner\"",
            // A reference is identified by hash, not by the path it was read from.
            "{\"id\": \"gene\", \"length\": 4, \"sha256\": \"1dff3e84fe7877e0673b69bbddcf4012\
             4e396e3f9943dd890c91b6a09adb9af0\"}",
            "\"kind\": \"s3\", \"id\": \"https://logan-pub.s3.amazonaws.com/u/\", \
             \"n_requested\": 1",
            "\"primary_metric\": \"matched\"",
            "\"prune_frac\": 0.01",
            "\"min_abund\": 5.0",
            "\"max_size_mb\": 1024",
            "\"summary_seqs\": false",
            // The run block, with the clock the test pinned.
            "\"status\": \"ok\"",
            "\"started\": \"2025-11-15T09:46:40Z\"",
            "\"n_records\": 1",
            "\"n_ok\": 1",
            // Every discard key is present even at zero, so "did not happen" is readable.
            "\"accessions\": {\"never_reached\": 0, \"no_anchors\": 0, \"no_path\": 0, \
             \"no_logan_data\": 0, \"too_large\": 0, \"subgraph_too_large\": 0, \
             \"error\": 0}",
            "\"alleles\": {\"beyond_max_alts\": 0, \
             \"no_rejoin\": 0, \"not_substitution\": 0, \"over_indel_cap\": 0, \
             \"frameshifting\": 0, \"window_not_shiftable\": 0, \"identical\": 0}",
            "\"distinct_primaries\": 1",
            "\"majority_corrections\": 0",
            // Per-accession decisions that used to reach only stderr.
            "\"subgraph_fraction\": 0.161189",
            "\"provisional\": true",
            "\"prune\": {\"outcome\": \"refused_shorter_span\", \"depth\": 1700.0, \"floor\": 17.0, \
             \"unitigs\": 4976, \"unitigs_kept\": 1734, \"unpruned_span\": [0, 3821], \
             \"pruned_span\": [424, 3821]}",
            // What the accession was built from. The Logan objects are unversioned and
            // mutable, so without this a rerun after an upstream regeneration is
            // indistinguishable from a rerun over the same bytes.
            "\"payload\": {\"source\": \"s3\", \"bytes\": 1861058, \"sha256\": \"aaaaaaaaaaaaaaaa\
              aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\", \
              \"etag\": \"d2dd288b238c3c88b46acdcd78860dce\", \
              \"last_modified\": \"2024-03-31T17:38:36Z\", \
              \"fetched_at\": \"2026-11-16T16:47:50Z\"}",
        ] {
            assert!(out.contains(want), "missing {want}: {out}");
        }
        // The sequence is named, not carried: 0.1.x wrote it here unconditionally and it was
        // half the file. `docs/summaries.md`.
        assert!(!out.contains("\"seq\":"), "the literal sequence must be opt-in: {out}");
        assert!(out.contains("\"seq_sha256\":"), "{out}");
    }

    /// A record whose bytes were never read says so, rather than omitting the block. The
    /// same rule `discards` follows: "did not happen" and "this version does not measure it"
    /// have to be different answers.
    #[test]
    fn a_record_with_nothing_to_record_still_carries_a_payload_block() {
        let mut s = summary(hap("ACC_primary", false), vec![]);
        s.results[0].payload = crate::types::Payload::default();
        let p = tmp("payload_none.json");
        write_summary(&s, &p).unwrap();
        let out = std::fs::read_to_string(&p).unwrap();
        std::fs::remove_file(&p).ok();
        assert!(out.contains("\"payload\": {\"source\": \"none\"},"), "{out}");
        // Absent members are omitted, not written as null: a record is on a byte budget.
        for absent in ["\"sha256\": null", "\"etag\": null", "\"bytes\": null"] {
            assert!(!out.contains(absent), "{absent} should be omitted: {out}");
        }
    }

    /// `--fetch-only` writes a summary of its own. It reconstructs nothing and never runs
    /// GraphAligner, so it must not claim a reference or an external tool it did not use,
    /// while still recording what it put in the cache.
    #[test]
    fn a_fetch_only_summary_names_no_reference_and_no_external_tool() {
        let mut s = summary(hap("ACC_primary", false), vec![]);
        s.results[0].haplotypes.clear();
        s.results[0].haplotype = None;
        s.results[0].payload = crate::types::Payload {
            source: crate::types::PayloadSource::Cache,
            bytes: Some(4096),
            sha256: Some("b".repeat(64)),
            ..crate::types::Payload::default()
        };
        let cfg = crate::config::Config::default();
        let run = super::RunInfo {
            graphaligner: None,
            cfg: &cfg,
            queries: &[],
            started: std::time::UNIX_EPOCH + std::time::Duration::from_secs(1_763_200_000),
            wall_seconds: 0.5,
            n_requested: 1,
        };
        let p = tmp("fetch_only.json");
        super::write_summary(&s, &run, &p).unwrap();
        let out = std::fs::read_to_string(&p).unwrap();
        std::fs::remove_file(&p).ok();
        assert!(out.contains("\"references\": null,"), "{out}");
        assert!(out.contains("\"external\": []"), "{out}");
        assert!(!out.contains("GraphAligner"), "no tool that did no work: {out}");
        assert!(
            out.contains("\"payload\": {\"source\": \"cache\", \"bytes\": 4096, \"sha256\": \"bbbb"),
            "{out}"
        );
        // Still a well-formed summary: every contract block, in order.
        let mut at = 0usize;
        for block in ["\"schema\"", "\"tool\"", "\"inputs\"", "\"settings\"", "\"run\"",
                      "\"discards\"", "\"cohort\"", "\"records\""] {
            let i = out[at..].find(block).unwrap_or_else(|| panic!("missing {block}: {out}"));
            at += i + block.len();
        }
    }

    fn alt_of(h: Haplotype) -> AltHaplotype {
        AltHaplotype {
            haplotype: h,
            fork: 1,
            ka_saturated: false,
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

    /// Multi-anchor mode must emit every reference's cov/pident in `summary.json`, not just
    /// the winner's. The whole point of the field is that the panel11 analysis (and every
    /// analysis like it) does not have to re-run BLAST to recover what the walker already
    /// computed.
    #[test]
    fn summary_json_contains_all_candidates_in_multi_anchor_mode() {
        use crate::types::CandidateReport;
        let mut s = summary(hap("ACC_primary", false), vec![]);
        s.results[0].anchor_used = Some("XBB15".into());
        s.results[0].candidates = vec![
            CandidateReport {
                record_id: "Wuhan".into(),
                pool_index: 0,
                pident: Some(97.5),
                cov: Some(0.62),
                length: Some(2350),
                nodes: Some(120),
                path_len: Some(4500),
                num_matches: Some(2288),
                aln_len: Some(2350),
                selected: false,
            },
            CandidateReport {
                record_id: "XBB15".into(),
                pool_index: 1,
                pident: Some(99.6),
                cov: Some(0.99),
                length: Some(3800),
                nodes: Some(160),
                path_len: Some(4600),
                num_matches: Some(3785),
                aln_len: Some(3800),
                selected: true,
            },
        ];
        let p = tmp("summary_candidates.json");
        write_summary(&s, &p).unwrap();
        let out = std::fs::read_to_string(&p).unwrap();
        assert!(out.contains("\"candidates\""), "missing candidates field: {out}");
        assert!(out.contains("\"record_id\": \"Wuhan\""), "loser reference missing: {out}");
        assert!(out.contains("\"record_id\": \"XBB15\""), "winner reference missing: {out}");
        assert!(out.contains("\"selected\": true"), "winner not flagged: {out}");
        assert!(out.contains("\"selected\": false"), "loser not flagged: {out}");
        assert!(out.contains("\"cov\": 0.62"), "loser cov missing: {out}");
        assert!(out.contains("\"pident\": 99.6"), "winner pident missing: {out}");
        std::fs::remove_file(&p).ok();
    }

    /// Single-anchor runs leave `candidates` empty. The JSON writer must OMIT the field
    /// entirely (not emit an empty array), so existing single-anchor summaries diff clean
    /// against the pre-multi-anchor shape.
    #[test]
    fn summary_json_omits_candidates_field_in_single_anchor_mode() {
        let s = summary(hap("ACC_primary", false), vec![]);
        assert!(s.results[0].candidates.is_empty(), "test setup: expected empty candidates");
        let p = tmp("summary_no_candidates.json");
        write_summary(&s, &p).unwrap();
        let out = std::fs::read_to_string(&p).unwrap();
        assert!(
            !out.contains("\"candidates\""),
            "single-anchor JSON must not carry a candidates field: {out}"
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
