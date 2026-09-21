//! Tunables. Defaults mirror `logan_path_walker.py` so the two agree out of the box.

use std::path::PathBuf;

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
    /// analysis, so `--prefer-majority=false` is the switch to reach for when reproducing an
    /// older result. 25 of the 26 accessions it touched moved by 1 to 11 sites and cleanly;
    /// the one exception warns, see `stages::variants::prefer_majority`.
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
    /// are deliberately picked up. `build_walker_codon_msa.py` feeds its input straight into
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
    /// Hard cap on alternative haplotypes per accession.
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
    /// If the 5-hop neighbourhood exceeds this fraction of the whole graph, the haplotype is
    /// likely under-determined and the result should not be trusted without review.
    ///
    /// Checked after `extract::expand`, not on anchor counts. Measured across 10 accessions
    /// at the default `min_shared_kmers`: every accession whose k-mer reconstruction matched
    /// BLAST sat at or below **0.033%**, while the single discordant one (`SRR8209082`) sat
    /// at **10.383%**. A 315x gap with nothing in between, so 1% has an order of magnitude
    /// of margin on both sides.
    pub subgraph_fraction_guard: f64,

    // ---- concurrency ----
    /// Threads streaming from S3. Bandwidth saturates around 8 streams at 17.4 MB/s aggregate,
    /// so raising this past that buys nothing.
    pub fetch_workers: usize,
    /// Threads running parse through phase. Size to cores.
    pub compute_workers: usize,
    /// Bounded queue depth between fetch and compute. This times the largest retained buffer
    /// is the hard ceiling on resident memory, and is the whole point of the bounded channel.
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
    /// Where per-accession GFA and GAF files go. Should be a tmpfs; payloads are ~38 KB and <1 KB.
    pub scratch_dir: PathBuf,
    /// GraphAligner executable, resolved on PATH if not absolute.
    pub graphaligner_bin: String,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            graph_k: 31,
            k: 21,
            hops: 5,
            min_abund: 5.0,
            min_vaf: 0.10,
            min_cov: 0.85,
            prune_frac: 0.01,
            min_abund_frac: 0.0,
            prefer_majority: true,
            allow_indel_alts: true,
            max_alt_indel: 30,
            allow_frameshift_alts: false,
            emit_alts: true,
            emit_vcf: false,
            emit_snv_table: false,
            max_alts: 32,
            alt_max_detour: 100,
            min_shared_kmers: 3,
            subgraph_fraction_guard: 0.01,
            fetch_workers: 8,
            compute_workers: 8,
            channel_capacity: 4,
            graphaligner_threads: 1,
            fetch_timeout_secs: 900,
            fetch_connect_timeout_secs: 15,
            fetch_attempts: 3,
            max_compressed_bytes: 1024 * 1024 * 1024,
            scratch_dir: PathBuf::from("/dev/shm/logan-walker"),
            graphaligner_bin: String::from("GraphAligner"),
        }
    }
}
