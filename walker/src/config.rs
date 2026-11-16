//! Tunables.

use std::path::PathBuf;

/// How multi-anchor mode ranks per-reference candidates in [`stages::stitch::pick_primary`].
///
/// The default is `Matched`: the fraction of the reference's own bases that the alignment
/// matched (`num_matches / reference length`), then identity, then pool order. It cannot
/// exceed 1.0, and it rewards coverage and identity together: a reference loses for every
/// base of its own the walk misses or mismatches.
///
/// `Cov` ranks on the walked length divided by the reference length, then identity. That
/// ratio is not capped, so a reference shorter than the sample's gene (one carrying a
/// deletion the sample lacks) scores above 1.0 and beats an exact match. On the SARS-CoV-2
/// spike 11-strain panel (20,000-accession run) it put 62% of primaries on `Omicron_XBB15`
/// and 21% on `Omicron_JN1`, the two shortest references, picking 6,085 of 9,448 primaries on
/// a ratio above 1.0; `Matched` puts 43% on Delta and 26% on BA.2, the cohort's real mix.
/// The choice changes the sequence, not just the tag: GraphAligner walks the graph
/// differently against each reference.
///
/// `Pident` ranks on identity, then the coverage ratio. It agrees with `Matched` on 96.5% of
/// that run, but it can prefer a truncated walk that matches its reference cleanly over a
/// full-length one.
///
/// `CovXPident` orders on `cov * pident/100`, then cov, then pident.
///
/// Only affects which candidate is tagged as the primary and written to the FASTA. Every
/// candidate's numbers are emitted in `summary.json` under `"candidates"` regardless, so
/// a downstream tool can always re-rank.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PrimaryMetric {
    Matched,
    Cov,
    Pident,
    CovXPident,
}

impl PrimaryMetric {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "matched" => Some(Self::Matched),
            "cov" => Some(Self::Cov),
            "pident" => Some(Self::Pident),
            "cov_x_pident" | "cov-x-pident" => Some(Self::CovXPident),
            _ => None,
        }
    }
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Matched => "matched",
            Self::Cov => "cov",
            Self::Pident => "pident",
            Self::CovXPident => "cov_x_pident",
        }
    }
}

#[derive(Clone)]
pub struct Config {
    // ---- algorithm ----
    /// Logan de Bruijn graph k-mer size. Fixed at 31 by Logan's BCALM2 index construction;
    /// consecutive unitigs overlap by graph_k - 1 = 30 bp.
    pub graph_k: usize,
    /// Anchoring k-mer size used to identify unitigs that carry the reference gene.
    /// Defaults to 21 based on empirical sensitivity simulations (0 FP in 4.64 Mb, rescues 85% divergence cliff).
    pub k: usize,
    /// Hop radius for the local subgraph BFS.
    pub hops: usize,
    /// Minimum `ka:f` for a bubble branch to count as a real variant.
    pub min_abund: f32,
    /// Minimum variant allele frequency for a bubble branch.
    pub min_vaf: f64,
    /// Minimum fraction of the reference a haplotype must cover to be kept.
    pub min_cov: f64,
    /// How multi-anchor mode ranks per-reference candidates. See [`PrimaryMetric`].
    ///
    /// Single-anchor runs ignore this: there is only one candidate to pick from.
    pub primary_metric: PrimaryMetric,
    /// Second-pass prune: drop unitigs below this fraction of the gene's own depth, then
    /// realign. `0.0` disables it and restores the previous behaviour exactly.
    ///
    /// **Defaults to 0.01**, on the strength of a 200-accession panel run (section 16 of
    /// `docs/issue_intra_accession_variance_lost.md`). Run three ways on one build, pruning
    /// never changed whether an accession reconstructed, 187 of 200 in every arm, so it
    /// carries no risk of losing a result. What it does carry:
    ///
    /// ```text
    /// alternatives emitted        540 -> 311, and 94% of the sub-1%-of-depth ones removed
    /// haplotypes changed           84 of 187
    ///   gap-free before/after      52 -> 72       fewer gaps 26, more gaps 5
    ///   identity to reference      higher 44, lower 38
    /// ```
    ///
    /// Twenty more coding sequences stop carrying an indel, which is a hard defect for a
    /// codon-aware alignment. Identity to the reference splits near evenly, which is what
    /// removing a reference bias should look like: the point is to come off the error
    /// branches, not to move toward the reference.
    ///
    /// # The problem it solves
    ///
    /// GraphAligner scores a route by resemblance to the *reference*. Where the sample's
    /// real sequence differs from the reference and a branch carried by a few erroneous
    /// reads happens to match it, the aligner takes the error. This is not a rare event:
    /// on `SRR8209082` the unpruned alignment crosses 118 unitigs below 1% of the gene's
    /// depth, some at 2x in a 14,500x sample, and the emitted sequence carries whatever
    /// those branches spell.
    ///
    /// # Why it must be a second pass
    ///
    /// The floor needs the gene's depth, and that number does not exist before an
    /// alignment. Deriving it from the anchor set was tried and fails badly: `SRR8209082`'s
    /// anchors have a median `ka:f` of 70x against a true 14,500x, because
    /// `min_shared_kmers = 3` admits thousands of shallow incidental matches. Aligning once
    /// to measure, then pruning and realigning, is what works.
    ///
    /// # Measured
    ///
    /// ```text
    /// SRR8209082 (HIV pol, 14,500x)   nodes  aligns    cov   min ka:f  under 1% depth
    ///   unpruned                      34582     3    0.9996      2.0        118
    ///   1% floor (145x)                6930     1    0.9993    147.5          2
    ///   5% floor (725x)                5343     1    0.9993   7000.0          0
    ///
    /// ERR2028820 (PfCRT, 180x)        identical haplotype at every setting, 0 nodes lost
    /// ```
    ///
    /// Three alignments collapsing to one is the intended outcome: they were three routes
    /// through the same noise, differing at 33 positions, not three haplotypes.
    ///
    /// # Effect on variant recovery
    ///
    /// This is what makes `stages::variants` trustworthy. On the unpruned graph that stage
    /// finds 12 apparently biallelic sites on the route, every one of them between 1x and
    /// 90x in a 14,500x sample, which is noise. On the pruned graph it finds exactly one,
    /// at 8,000x against 3,600x, which is a real second population.
    ///
    /// 0.01 is the value to start from. Pruning never regresses an accession: if the second
    /// alignment fails, the first pass's haplotype stands.
    pub prune_frac: f64,
    /// Raise the variant abundance floor to this fraction of the gene's own depth.
    ///
    /// `0.0` disables it, which is the default and leaves `min_abund` alone. Otherwise the
    /// floor used by `stages::variants` is `max(min_abund, min_abund_frac * depth)`, where
    /// depth is the median `ka:f` over the walked path.
    ///
    /// # Why
    ///
    /// `min_abund = 5.0` is absolute, so what it means depends entirely on how deeply the
    /// sample was sequenced. Measured gene depths: `ERR899743` 230x, `ERR2028820` 180x,
    /// `SRR8209082` 15,000x. Five reads is a sane 2% floor on the first two and 0.03% on the
    /// third, where it admits branches carried by a handful of erroneous reads.
    ///
    /// Measured on `SRR8209082`, swept and grouped by walked path (the path itself is not
    /// stable on this accession, see below), every recovered allele falls into one of two
    /// clearly separated classes:
    ///
    /// ```text
    /// ka:f 5.7  VAF 0.182   error, ~6 reads in a 14,500x sample
    /// ka:f 7.7  VAF 0.106   error, ~8 reads
    /// ka:f 3600 VAF 0.310   a real second population
    /// ```
    ///
    /// A 500x gap with nothing in between, so the floor has a wide plateau: anything in
    /// **0.001 to 0.05** drops both error calls and keeps the real one. `min_vaf` alone
    /// cannot do this, because the error calls sit at 18% and 11% of their *local* coverage;
    /// it is only against the gene's depth that they are visible as noise.
    ///
    /// Left at 0.0 regardless, because one accession with genuine variants is not enough to
    /// calibrate a default on. 0.01 is the value to try first.
    ///
    /// # Why the path and not the anchors
    ///
    /// The anchor set is not a usable depth estimate: at `min_shared_kmers = 3` it admits
    /// thousands of shallow incidental matches, and its median `ka:f` on `SRR8209082` is 70x
    /// against a true 15,000x. The walked path *is* the gene, so its median is the thing
    /// being asked about. This is the estimate the abandoned global prune floor needed and
    /// could not get, because it had to run before alignment.
    pub min_abund_frac: f64,
    /// At a site where two alleles are both real, make the primary carry the more abundant
    /// one. **On by default.**
    ///
    /// Without it the emitted primary carries the *minority* allele at 19% of the two-allele
    /// sites it reports (59 of 311 across the panel); with it, 1 of 311. The invariant to
    /// look for is the maximum alternative VAF, which moves from 0.883 to 0.499: above 0.5
    /// would mean the primary had taken the rarer allele.
    ///
    /// It is defaulted on because a consensus carrying the rarer base is simply wrong, and
    /// because the corrections replicate: five independent `p_infestans_avr3a` isolates were
    /// corrected at the same two positions, 119 and 132, all five having matched the
    /// reference exactly beforehand. Independent samples agreeing at identical coordinates
    /// is not a noise process.
    ///
    /// It does change the sequence that feeds the codon alignment and the selection
    /// analysis, so `--no-prefer-majority` is the switch to reach for when reproducing an
    /// older result. 25 of the 26 accessions it touched moved by 1 to 11 sites and cleanly;
    /// the one exception warns, see `stages::variants::prefer_majority`.
    ///
    /// Since 0.2.0 this also moves the route off a branch the aligner walked that is itself a
    /// minority of its fork, when the other side holds an outright majority (WV-02). Before
    /// that, such a fork was skipped as "the error-branch case" and left to `--prune-frac`,
    /// which does not reach it: pruning is 1% of gene depth, the old gate was 10% of the fork,
    /// and a branch between the two was neither pruned nor corrected nor reported.
    ///
    /// Pruning cannot fix these: no coverage floor separates a 35% allele from a 65% one
    /// without also deleting every genuine minority variant, so the aligner keeps choosing
    /// by resemblance to the reference. Measured over the first 68 two-allele sites of a
    /// panel run, with pruning already on, the alternative was the more abundant allele at
    /// 32 of them, and clearly so (0.65 or above) at 18.
    ///
    /// See `stages::variants::prefer_majority`, including why this is not the whole-route
    /// repair that failed earlier.
    pub prefer_majority: bool,
    /// Accept alternatives that change the sequence length, not only substitutions.
    ///
    /// **On by default**, but narrowly: see `max_alt_indel` and `allow_frameshift_alts`,
    /// which together admit only the subset that is safe to publish.
    ///
    /// Length-changing alleles were refused outright until `variants::shifted_window`
    /// existed, because `path_start` and `path_end` are offsets into the stitched sequence
    /// and reusing them after a splice that adds or removes bases slices the wrong window.
    /// The retired upstream `detour` engine had no such check, which was most of why it emitted
    /// alternatives carrying gaps at 42% against this stage's 2%.
    pub allow_indel_alts: bool,
    /// Largest length change, in bases, that still counts as an indel allele.
    ///
    /// Measured over the panel, length-changing calls cluster at +-1 and +-2 with a long
    /// tail at +100, +142 and -282. The tail is not an allele: it is the route being sent
    /// through a different stretch of the graph, rewriting the gene rather than editing it.
    /// 30 keeps a ten-codon indel and excludes all of that.
    /// For `find_alts` this bounds one allele. For `prefer_majority` it bounds **the whole
    /// correction set** on an accession as well, because the corrections compound: on the
    /// WV-02 panel `SRR6007473` took 64 swaps that were each in frame and each inside this
    /// cap, and its ebola GP primary grew from 2,049 nt to 3,129 against a 2,031 nt
    /// reference. Per-swap legality does not bound the total.
    pub max_alt_indel: usize,
    /// Accept indel alleles whose length change is not a multiple of three.
    ///
    /// **Off by default.** Such an allele shifts the reading frame of everything downstream.
    /// Real frameshift alleles exist, but here they cannot be told apart from the
    /// homopolymer-length artifacts a de Bruijn graph produces, and 190 of the 203
    /// length-changing calls on the panel are this case with a median magnitude of one base.
    /// In a codon-aware selection analysis a frameshifted record does more damage than a
    /// missing one, so recovering these is opt-in.
    pub allow_frameshift_alts: bool,
    /// Emit a second record for each genuinely biallelic site on the walked path.
    ///
    /// **On by default**, and safe to be: the records go to `alt_haplotypes.fasta`, not to
    /// the FASTA the rest of the pipeline reads, so nothing downstream changes unless they
    /// are deliberately picked up. `disassembler msa` feeds its input straight into
    /// MAFFT, FastTree and HyPhy, where every sequence becomes an independent taxon, and two
    /// records from one sample are not two independent observations. Keeping them in a
    /// separate file means the variance is recovered and auditable while using it in an
    /// analysis stays a deliberate act. See `stages::variants`.
    pub emit_alts: bool,
    /// Also write a per-accession VCF of intra-host SNPs (the on-path biallelic sites), each
    /// carrying its k-mer-coverage frequency, in query coordinates. Off by default; the
    /// primary FASTA is unaffected. See `io::output::write_vcfs`.
    pub emit_vcf: bool,
    /// Also write one long ("tidy") table of every intra-host SNP across the whole run, one row
    /// per (accession, site, allele), to `<out-dir>/snvs.tsv`. This is the cohort output that
    /// scales: a merged multi-sample VCF needs one column per accession, so a gene matched by
    /// millions of runs would be a million-column near-empty matrix; the long form is
    /// sparse-native and concatenates across cluster shards. See `io::output::write_snv_table`.
    pub emit_snv_table: bool,
    /// Move the route off a branch the aligner walked that is itself a minority of its fork,
    /// when another branch holds an outright majority. **Off by default.**
    ///
    /// This is WV-02's action, and it is off because a 5,000-accession slice of walker20k
    /// showed the correction loop does not converge with it on: 97,068 corrections where the
    /// pre-0.2.0 walker made 470, a median of 64 per corrected accession (the loop's own
    /// `max_alts * 2` backstop) on 1,399 accessions, **7.8 corrections for every net
    /// substitution**, and a net +33 internal stop codons. Sixty-four swaps producing six
    /// bases of change means the swaps are undoing one another.
    ///
    /// WV-02's *finding* stands and its reporting is unconditional: forks where `_primary`
    /// carries the rarer allele are counted in `minority_left` whether or not they are acted
    /// on. What is off is the acting.
    pub rescue_minority: bool,
    /// Write the literal haplotype sequence into the summary, as 0.1.x did.
    ///
    /// Off by default from 0.2.0. The summary carries `seq_sha256` and `length` instead, which
    /// answer the questions the field was used for (equality between runs, counting distinct
    /// primaries, joining a record to its FASTA entry) without duplicating
    /// `reconstructed.fasta`. Measured on three real runs, the sequence was 50.1% to 71.4% of
    /// the file; at walker20k scale that is roughly 67 MB of duplicate FASTA inside a 130 MB
    /// summary. See `docs/summaries.md`.
    pub summary_seqs: bool,
    /// Cap on alternative haplotypes per accession, `--max-alts`.
    ///
    /// WV-06: this used to have no flag and to `break` the 5'-to-3' scan, so an accession
    /// with more credible sites than the cap published the start of its gene and nothing
    /// after, with the stderr counters stopping at the cap so that even they did not reveal
    /// it. The scan now runs to the end, counts what did not fit as `beyond_cap`, and warns.
    ///
    /// Measured site counts are 4, 4 and 0 on the accessions looked at, so this is a
    /// backstop against a pathological neighbourhood rather than a tuning knob.
    pub max_alts: usize,
    /// Unitigs an alternative allele may traverse before it has to rejoin the main route.
    ///
    /// Counted in unitigs, not bases, and the two are not proportional. Pruning removes the
    /// long low-coverage unitigs, so the surviving route is built from many short ones: on
    /// `SRR8209082` the route grows from 1,319 unitigs to 2,479 while the sequence stays
    /// ~2,737 bases. The same allele therefore spans far more unitigs after pruning.
    ///
    /// Measured on that accession's pruned graph, the one real biallelic site is found at
    /// 60, 150 and 400 and missed at 20, so the default sits inside that plateau. A generous
    /// value is safe here because the length check, not this limit, is what rejects
    /// shortcuts: an alternative that rejoins far downstream changes the sequence length and
    /// is discarded whatever its detour.
    pub alt_max_detour: usize,

    // ---- anchoring ----
    /// Minimum shared canonical k-mers for a unitig to be an anchor.
    /// Validated window is 1..=5; haplotypes shift at 10. Total count is used rather than
    /// longest consecutive run: the two agree in 22 of 24 measured cells, because a unitig is
    /// contiguous, and the plain count needs one counter instead of three.
    pub min_shared_kmers: u32,
    /// If the 5-hop neighbourhood holds more than this many unitigs per reference base, the
    /// haplotype is likely under-determined and the result should not be trusted without
    /// review. Checked after `extract::expand`, not on anchor counts.
    ///
    /// Through 0.1.1-dev this was a share of the whole graph (`subgraph_fraction`) at 1%,
    /// calibrated on 10 accessions where every concordant one sat at or below 0.033% and the
    /// single discordant one at 10.383%. That calibration only holds when the target is a
    /// sliver of the sample: on targeted or amplicon data the target *is* most of the graph,
    /// so the share says more about the library than the haplotype. Measured on the
    /// 20,000-accession SARS-CoV-2 run, it flagged 9,440 of the 9,448 successes, and the 8 it
    /// cleared had a *lower* complete-ORF rate (50.0%) than the ones it flagged (67.0%): no
    /// signal at all, in the wrong direction.
    ///
    /// Unitigs per reference base measures how tangled the graph is around the gene, which is
    /// the thing that actually makes a walk ambiguous, and it does not move with the library's
    /// breadth. On the same run, successes below 3 were 89.7% complete-ORF against 51.3% above,
    /// and the cohorts where the walk is easy sit far below the cut (obelisk epsilon's median
    /// is 0.3) rather than all being flagged.
    ///
    /// `subgraph_fraction` is still reported per result; it is a useful number, just not a gate.
    pub subgraph_nodes_per_base_guard: f64,
    /// Refuse an accession whose anchor neighbourhood exceeds this many unitigs.
    ///
    /// The memory term the two-pass parser cannot remove. Everything downstream of
    /// `extract::expand` is proportional to the neighbourhood, at roughly 150 bytes per node
    /// in RAM plus about 157 bytes per node of GFA on scratch, and on amplicon-like data the
    /// neighbourhood is most of the graph rather than a corner of it: a synthetic HIV cohort
    /// reached 1,383,280 nodes, 42% of its graph.
    ///
    /// The default is far above anything real measured here: `ERR2028820`/PfCRT is 747 nodes,
    /// `SRR12363134`/ERG11 is 296, and `SRR8209082`/HIV-1 *pol* is 34,582, which is already
    /// flagged provisional by `subgraph_nodes_per_base_guard`. At 1,000,000 the ceiling is
    /// about 150 MB of graph per compute worker, so it bounds the run without touching any
    /// result this repo has ever produced. 0 disables it.
    ///
    /// This is what `--max-size-mb` can eventually be relaxed against: the compressed size cap
    /// exists to bound memory, and once memory is bounded by the neighbourhood instead, the
    /// size cap is a transfer guardrail rather than a correctness skip.
    pub max_subgraph_nodes: usize,

    // ---- concurrency ----
    /// Threads streaming from S3. Bandwidth saturates around 8 streams at 17.4 MB/s aggregate,
    /// so raising this past that buys nothing.
    pub fetch_workers: usize,
    /// Threads running parse through phase.
    ///
    /// A fixed default of 8, deliberately **not** `available_parallelism()`, though this and
    /// `pipeline`'s diagram both used to say "= cores". Auto-sizing would silently raise the
    /// memory bound on a large machine: each compute worker holds a compressed buffer through
    /// the whole parse, so 32 cores would mean 4 + 8 + 32 = 44 buffers, 44 GB at the default
    /// `--max-size-mb`. See `peak_compressed_buffers`. Set it explicitly with `-w` once you
    /// know the machine can absorb it; the startup line prints what you chose.
    pub compute_workers: usize,
    /// Bounded queue depth between fetch and compute.
    ///
    /// **Zero, a rendezvous channel.** A positive depth buys nothing here and costs one
    /// capped buffer per slot. A fetch worker that has bytes in hand and no compute worker
    /// ready is already a held buffer, counted by the `fetch_workers` term of
    /// `peak_compressed_buffers`; queueing those bytes separately just moves the same
    /// allocation into a second place and adds a slot's worth of cap to the bound. The depth
    /// used to be 4, which at the 1 GB default put 4 GB into the printed worst case for no
    /// throughput: measured over 20 cached 34 MB accessions at `-w 1 --fetch-workers 8`, peak
    /// RSS fell from 516.4 to 386.5 MiB, which is the four slots exactly, with the wall clock
    /// unchanged at 18.0 s.
    ///
    /// Left configurable rather than removed, because a slow network with a fast compute pool
    /// is the one shape where a slot could help, and the startup line prints what is in force.
    pub channel_capacity: usize,
    /// Threads handed to each GraphAligner invocation.
    ///
    /// **One, for determinism.** GraphAligner's multi-threaded search is not reproducible:
    /// on the same graph and the same query it sometimes finds a better alignment and
    /// sometimes does not. Measured on `SRR8209082`'s subgraph, best-alignment score over
    /// five identical invocations:
    ///
    /// ```text
    /// -t 1    2691  2691  2691  2691  2691     deterministic, and the best score
    /// -t 2    2686  2691  2691  2691  2691     varies run to run
    /// -t 8    2686  2686  2686  2686  2686     stable, but consistently worse
    /// ```
    ///
    /// At `-t 2` this reached the pipeline's output: five identical runs produced two
    /// different primary haplotypes, because a run that missed the 2691 alignment emitted
    /// the 2686 one instead. A walker whose result depends on thread scheduling cannot be
    /// pinned by a fixture or trusted downstream, so this is not a tuning knob to raise.
    ///
    /// It costs nothing. On that graph, the largest here at 34,582 nodes, one thread takes
    /// 0.47s against 0.41s for eight. Throughput comes from `compute_workers` running
    /// whole accessions in parallel, which is the axis that actually scales; splitting one
    /// alignment across threads only buys a few hundredths of a second and forfeits
    /// reproducibility to get it.
    ///
    /// Smaller graphs never showed the effect: `ERR2028820` scored 1268 on all five runs at
    /// every thread count. It appears where the alignment is ambiguous enough to search.
    pub graphaligner_threads: usize,

    // ---- fetch ----
    /// Whole-request timeout, including body transfer. A 343 MB accession at the measured
    /// 6.9 MB/s single-stream rate takes ~50 s, so this must be generous.
    pub fetch_timeout_secs: u64,
    /// TCP + TLS connect timeout, kept short so a dead connection fails fast and retries.
    pub fetch_connect_timeout_secs: u64,
    /// Attempts per accession before giving up. Retries apply to transport errors and 5xx
    /// only; a 404 means the accession has no Logan unitigs and is not retried.
    pub fetch_attempts: u32,
    /// Refuse accessions whose compressed size exceeds this.
    ///
    /// This is a safety valve, not a tuning knob. Compressed size is held in memory per
    /// in-flight accession, and decompressed size runs 4-5x higher again. Default 1 GB still
    /// drops multi-gigabyte environmental metagenomes in 50 ms via Content-Length, while
    /// keeping the deep viral runs the earlier 350 MB cap refused (it cost ~7.5% of the panel).
    pub max_compressed_bytes: u64,

    // ---- io ----
    /// Where per-accession GFA, GAF and query files go; payloads are ~38 KB and <1 KB. `main`
    /// points it at the run's private directory (see scratch.rs), so this default is only what
    /// tests get.
    pub scratch_dir: PathBuf,
    /// GraphAligner executable, resolved on PATH if not absolute.
    pub graphaligner_bin: String,

    // ---- offline cache ----
    /// On-disk read-through cache for compressed unitig payloads.
    ///
    /// When set, `Fetcher::fetch` looks for `<dir>/<ACC>.unitigs.fa.zst` before hitting S3.
    /// A miss with `offline == false` falls through to the normal S3 fetch and writes the
    /// downloaded bytes back to the cache atomically (tempfile + rename in the same dir).
    ///
    /// Files are byte-identical to the S3 payload, so `sha256sum` against the original object
    /// verifies. Layout is flat: one file per accession, no manifest.
    pub unitigs_cache: Option<PathBuf>,
    /// Refuse the S3 fallback. Only meaningful with `unitigs_cache`.
    ///
    /// A cache miss becomes a hard error naming the missing accession, rather than a silent
    /// download. Off by default; opt in for reproducibility runs or air-gapped hosts.
    pub offline: bool,
    /// Where adjacency spills when it outgrows `types::SPILL_ADJACENCY_AT`.
    ///
    /// Resolved in `main` as `--spool-dir`, else `--unitigs-cache`, else the system temp
    /// directory. Deliberately **not** `--scratch-dir`, which defaults to `/dev/shm`: spilling
    /// into a tmpfs moves the bytes from one part of RAM to another and measures as a saving
    /// while being none.
    ///
    /// `None` keeps adjacency in memory however large it gets, which is the old behaviour and
    /// what every test uses.
    pub spool_dir: Option<PathBuf>,
    /// Fetch every accession into the cache and exit before the compute stages run.
    ///
    /// Requires `unitigs_cache`. `pipeline::run_fetch_only` handles this shortcut; the compute
    /// pool is not spawned, and the run's success criterion is "every requested accession is
    /// now in the cache". This is what `disassembler fetch` wraps.
    pub fetch_only: bool,
}

impl Config {
    /// Compressed unitig buffers that can be resident at the same time.
    ///
    /// **Not** `channel_capacity`, which is the number this file and `pipeline`'s module doc
    /// both used to claim. Three places hold one concurrently:
    ///
    /// * `channel_capacity` sitting in the bounded queue (zero by default, see the field);
    /// * one per fetch worker, held from the moment `read_to_vec` returns until the blocking
    ///   `send` into a full queue completes (`stages::stream`);
    /// * one per compute worker, because `run_one` keeps `fetched` alive through the whole
    ///   parse and drops it only once the CSR graph exists (`pipeline::run_one`).
    ///
    /// At the defaults that is 0 + 8 + 8 = 16; before the queue depth went to zero it was 20,
    /// five times the figure this file used to document.
    pub fn peak_compressed_buffers(&self) -> usize {
        self.peak_compressed_buffers_for(usize::MAX)
    }

    /// The same bound for a run that has only `n_accessions` of work.
    ///
    /// The pipeline's shape is an upper bound on a long list; it is not what a short one can
    /// reach. A single-accession run cannot fill 16 slots, so quoting 16 GB at it overstates
    /// the requirement by the worker count. That matters because the number is read and acted
    /// on: the Galaxy wrapper derives its memory ask from this formula, and a user sizing a
    /// job reads the startup line.
    ///
    /// `--fetch-only` never starts the compute pool, so it holds neither a queue slot nor a
    /// parsing buffer: its ceiling is one buffer per fetch worker.
    pub fn peak_compressed_buffers_for(&self, n_accessions: usize) -> usize {
        let slots = if self.fetch_only {
            self.fetch_workers
        } else {
            self.channel_capacity + self.fetch_workers + self.compute_workers
        };
        slots.min(n_accessions)
    }

    /// Parsed graphs that can be resident at once: one per compute worker, never more than
    /// the work available, and none at all under `--fetch-only`, which does not parse.
    pub fn peak_graphs_for(&self, n_accessions: usize) -> usize {
        if self.fetch_only {
            0
        } else {
            self.compute_workers.min(n_accessions)
        }
    }

    /// Worst-case bytes resident in compressed buffers alone, if every accession is at the
    /// `--max-size-mb` cap.
    ///
    /// This is a bound, not a prediction: real accessions are far below the cap (`ERR2028820`
    /// is 33 MB, the ERG11 panel reaches 343 MB against a 1 GB default). It excludes the parsed
    /// graphs (`max_graph_resident_bytes`) and any transient growth inside `read_to_vec` while
    /// the receiving `Vec` reallocates.
    pub fn max_compressed_resident_bytes(&self) -> u64 {
        self.max_compressed_resident_bytes_for(usize::MAX)
    }

    /// As above, for a run with `n_accessions` of work.
    pub fn max_compressed_resident_bytes_for(&self, n_accessions: usize) -> u64 {
        (self.peak_compressed_buffers_for(n_accessions) as u64)
            .saturating_mul(self.max_compressed_bytes)
    }

    /// The budget one parsed graph may use, which is what `stages::parse` enforces.
    ///
    /// `max_graph_resident_bytes` is this times the number of compute workers. Splitting it
    /// out makes the per-accession limit the single definition both the printed bound and the
    /// parser work from, so the promise and the enforcement cannot drift apart.
    pub fn max_graph_bytes_per_accession(&self) -> u64 {
        (self.max_compressed_bytes as f64 * GRAPH_PEAK_PER_COMPRESSED_BYTE) as u64
    }

    /// Resident bytes one neighbourhood costs, once adjacency has gone to disk.
    ///
    /// About 150 bytes per node: the retained rows, the 2-bit sequence, the offsets and the
    /// abundance. This is what bounds a run when `spool_dir` is set and
    /// `max_subgraph_nodes` is finite, and it does not mention `--max-size-mb` at all, which
    /// is the point of stage 3: memory stops scaling with the file and starts scaling with the
    /// gene.
    pub fn max_subgraph_resident_bytes(&self) -> Option<u64> {
        if self.spool_dir.is_none() || self.max_subgraph_nodes == 0 {
            return None;
        }
        Some((self.max_subgraph_nodes as u64).saturating_mul(BYTES_PER_SUBGRAPH_NODE))
    }

    /// Worst-case bytes in parsed graphs: one per compute worker, each from an accession at
    /// the cap. See `GRAPH_PEAK_PER_COMPRESSED_BYTE`.
    pub fn max_graph_resident_bytes(&self) -> u64 {
        self.max_graph_resident_bytes_for(usize::MAX)
    }

    /// As above, for a run with `n_accessions` of work.
    pub fn max_graph_resident_bytes_for(&self, n_accessions: usize) -> u64 {
        // A float-to-integer cast saturates, so an unlimited cap stays u64::MAX.
        (self.peak_graphs_for(n_accessions) as f64 * self.max_compressed_bytes as f64
            * GRAPH_PEAK_PER_COMPRESSED_BYTE) as u64
    }

    /// Worst-case resident bytes for the run: compressed buffers plus parsed graphs. What the
    /// startup line prints.
    pub fn max_resident_bytes(&self) -> u64 {
        self.max_resident_bytes_for(usize::MAX)
    }

    /// As above, for a run with `n_accessions` of work. This is what the startup line prints.
    pub fn max_resident_bytes_for(&self, n_accessions: usize) -> u64 {
        self.max_compressed_resident_bytes_for(n_accessions)
            .saturating_add(self.max_graph_resident_bytes_for(n_accessions))
    }
}

/// Peak memory of a parsed `UnitigGraph`, per byte of its compressed unitig file.
///
/// Measured 2026-10-31 as the process's peak RSS, minus a ~7 MB baseline and the compressed
/// buffer, one accession per run: 1.24x to 1.62x across seven accessions of 4.8 to 348 MB
/// (`ERR2028820` 1.25x, `SRR15477033` 1.62x). The parser's arrays grow by doubling and are
/// not shrunk, so the peak includes their slack and the moment a growing array is copied.
/// This replaces a fixed ~43 MB per worker, which was the array lengths of one 33 MB
/// accession; the 348 MB one's graph peaked at 498 MB.
pub const GRAPH_PEAK_PER_COMPRESSED_BYTE: f64 = 1.7;

/// Resident bytes per retained neighbourhood node, once adjacency is spooled.
///
/// The rows, the 2-bit sequence, the offsets and the abundance. Measured as the slope of peak
/// RSS against neighbourhood size rather than derived from the field widths, so it includes
/// the allocator's slack.
pub const BYTES_PER_SUBGRAPH_NODE: u64 = 150;

impl Default for Config {
    fn default() -> Self {
        Self {
            graph_k: 31,
            k: 21,
            hops: 5,
            min_abund: 5.0,
            min_vaf: 0.10,
            min_cov: 0.85,
            primary_metric: PrimaryMetric::Matched,
            prune_frac: 0.01,
            min_abund_frac: 0.0,
            prefer_majority: true,
            allow_indel_alts: true,
            max_alt_indel: 30,
            allow_frameshift_alts: false,
            emit_alts: true,
            emit_vcf: false,
            emit_snv_table: false,
            rescue_minority: false,
            summary_seqs: false,
            max_alts: 32,
            alt_max_detour: 100,
            min_shared_kmers: 3,
            subgraph_nodes_per_base_guard: 3.0,
            max_subgraph_nodes: 1_000_000,
            fetch_workers: 8,
            compute_workers: 8,
            channel_capacity: 0,
            graphaligner_threads: 1,
            fetch_timeout_secs: 900,
            fetch_connect_timeout_secs: 15,
            fetch_attempts: 3,
            max_compressed_bytes: 1024 * 1024 * 1024,
            scratch_dir: std::env::temp_dir().join("logan-walker-test"),
            graphaligner_bin: String::from("GraphAligner"),
            unitigs_cache: None,
            spool_dir: None,
            offline: false,
            fetch_only: false,
        }
    }
}

#[cfg(test)]
mod memory_bound_tests {
    use super::Config;

    /// The documented ceiling was `channel_capacity` alone. Three places hold a compressed
    /// buffer at once, so at the defaults the real figure is 20, not 4.
    #[test]
    fn peak_buffers_counts_all_three_holders() {
        let cfg = Config::default();
        assert_eq!((cfg.channel_capacity, cfg.fetch_workers, cfg.compute_workers), (0, 8, 8));
        assert_eq!(cfg.peak_compressed_buffers(), 16);
        assert_ne!(
            cfg.peak_compressed_buffers(),
            cfg.channel_capacity,
            "the queue depth alone is the claim this check exists to refute"
        );
    }

    /// 20 buffers at the 1 GB default cap is ~20 GB, the figure that made this a finding.
    #[test]
    fn worst_case_bytes_track_the_size_cap() {
        let cfg = Config::default();
        assert_eq!(cfg.max_compressed_bytes, 1024 * 1024 * 1024);
        assert_eq!(cfg.max_compressed_resident_bytes(), 16 * 1024 * 1024 * 1024);

        // The ERG11 panel's largest accession, a real-data case.
        let panel = 343u64 * 1024 * 1024;
        let at_panel_size = cfg.peak_compressed_buffers() as u64 * panel;
        assert!(
            (5.2..5.8).contains(&(at_panel_size as f64 / (1024.0 * 1024.0 * 1024.0))),
            "~5.4 GB expected, got {at_panel_size} bytes"
        );
    }

    /// Every term is user-tunable, which is why the binary prints the number rather than
    /// documenting a constant.
    #[test]
    fn every_term_moves_the_bound() {
        let base = Config::default().peak_compressed_buffers();
        for cfg in [
            Config { channel_capacity: 4, ..Config::default() },
            Config { fetch_workers: 16, ..Config::default() },
            Config { compute_workers: 16, ..Config::default() },
        ] {
            assert!(cfg.peak_compressed_buffers() > base);
        }
    }

    /// An absurd `--max-size-mb` must not wrap the bound to a small number.
    #[test]
    fn an_absurd_cap_saturates_instead_of_wrapping() {
        let cfg = Config { max_compressed_bytes: u64::MAX, ..Config::default() };
        assert_eq!(cfg.max_compressed_resident_bytes(), u64::MAX);
        assert_eq!(cfg.max_graph_resident_bytes(), u64::MAX);
        assert_eq!(cfg.max_resident_bytes(), u64::MAX);
    }

    /// The parsed graphs come on top of the buffers: one per compute worker, up to 1.7x the
    /// cap each, so the defaults' worst case is ~33.6 GB, not the ~20 GB once printed.
    #[test]
    fn the_bound_includes_the_parsed_graphs() {
        let cfg = Config::default();
        let gib = |b: u64| b as f64 / (1024.0 * 1024.0 * 1024.0);
        assert!((gib(cfg.max_graph_resident_bytes()) - 8.0 * 1.7).abs() < 1e-6);
        // 16 GB of buffers plus 13.6 GB of graphs. It was 33.6 until the queue depth went to 0.
        assert!((gib(cfg.max_resident_bytes()) - 29.6).abs() < 1e-6);
        let more = Config { compute_workers: 16, ..Config::default() };
        assert!(more.max_graph_resident_bytes() > cfg.max_graph_resident_bytes());
    }

    /// The pipeline's shape is a ceiling for a long list, not a figure to quote at a short
    /// one. A single-accession run reserving 29.6 GB is what the Galaxy wrapper used to be
    /// told to ask for.
    #[test]
    fn the_bound_is_capped_by_the_work_available() {
        let cfg = Config::default();
        assert_eq!(cfg.peak_compressed_buffers_for(1), 1);
        assert_eq!(cfg.peak_graphs_for(1), 1);
        assert_eq!(cfg.peak_compressed_buffers_for(4), 4);
        // Past the pipeline's own ceiling the accession count stops mattering.
        assert_eq!(cfg.peak_compressed_buffers_for(16), 16);
        assert_eq!(cfg.peak_compressed_buffers_for(10_000), 16);
        assert_eq!(cfg.peak_graphs_for(10_000), 8);
        // 1 x 1 GB + 1 x 1.7 GB, against 29.6 GB for a list long enough to fill the pipeline.
        let gib = |b: u64| b as f64 / (1024.0 * 1024.0 * 1024.0);
        assert!((gib(cfg.max_resident_bytes_for(1)) - 2.7).abs() < 1e-6);
        assert_eq!(cfg.max_resident_bytes_for(usize::MAX), cfg.max_resident_bytes());
    }

    /// `--fetch-only` never starts the compute pool, so it holds neither a queue slot nor a
    /// parsing buffer and builds no graph. Counting those would have the mode reserve memory
    /// for stages it does not run.
    #[test]
    fn fetch_only_counts_neither_the_queue_nor_the_graphs() {
        let cfg = Config { fetch_only: true, channel_capacity: 4, ..Config::default() };
        assert_eq!(cfg.peak_compressed_buffers_for(usize::MAX), cfg.fetch_workers);
        assert_eq!(cfg.peak_graphs_for(usize::MAX), 0);
        assert_eq!(cfg.max_graph_resident_bytes_for(usize::MAX), 0);
        assert_eq!(
            cfg.max_resident_bytes_for(usize::MAX),
            cfg.fetch_workers as u64 * cfg.max_compressed_bytes
        );
        // Still capped by the list.
        assert_eq!(cfg.peak_compressed_buffers_for(2), 2);
    }
}
