//! Tunables.

use std::path::PathBuf;

/// How multi-anchor mode ranks per-reference candidates in `stages::stitch::pick_primary`.
///
/// - `Matched` (default): `num_matches / reference length`, then identity, then pool order.
///   Capped at 1.0, so a reference loses for every base of its own the walk misses.
/// - `Cov`: walked length over reference length, then identity. Uncapped, so a reference
///   shorter than the sample's gene scores above 1.0 and beats an exact match.
/// - `Pident`: identity, then coverage. Can prefer a truncated walk that matches cleanly.
/// - `CovXPident`: `cov * pident/100`, then cov, then pident.
///
/// This decides the sequence, not just the `anchor=` tag, because GraphAligner walks the
/// graph differently against each reference. Every candidate's numbers reach `summary.json`
/// under `"candidates"` either way, so a downstream tool can re-rank.
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
    ///
    /// 21 by simulation: no false positive in a 4.64 Mb background, and it extends recovery
    /// past the divergence cliff k = 31 falls off. See
    /// `bench/simulated_kmer_sensitivity/README.md`.
    pub k: usize,
    /// Hop radius for the local subgraph BFS.
    pub hops: usize,
    /// Minimum `ka:f` for a bubble branch to count as a real variant.
    ///
    /// The **allele** floor, and nothing else: `prune_min_abund` is the prune's own lower
    /// bound, so raising this to be stricter about variants does not also delete unitigs.
    pub min_abund: f32,
    /// Minimum variant allele frequency for a bubble branch.
    pub min_vaf: f64,
    /// Minimum fraction of the reference a haplotype's alignment must span to be kept.
    ///
    /// Measured on the **reference span** (`Haplotype::ref_cov`, GAF columns 3 and 4), not on
    /// the walked length (`Haplotype::cov`). The two differ by exactly the sample's indels, so
    /// gating on the walked length would call a gene with a deletion absent. See
    /// `Haplotype::ref_cov`.
    pub min_cov: f64,
    /// How multi-anchor mode ranks per-reference candidates. See [`PrimaryMetric`].
    ///
    /// Single-anchor runs ignore this: there is only one candidate to pick from.
    pub primary_metric: PrimaryMetric,
    /// Second-pass prune: drop unitigs below this fraction of the gene's own depth, then
    /// realign. `0.0` disables it.
    ///
    /// Defaults to 0.01. GraphAligner scores a route by resemblance to the reference, so
    /// where the sample really differs from it and a branch carried by a few erroneous reads
    /// happens to match, the walk takes the error: on `SRR8209082` the unpruned alignment
    /// crosses 118 unitigs below 1% of the gene's depth, some at 2x in a 14,500x sample.
    ///
    /// It must be a second pass because the floor needs the gene's depth, which does not
    /// exist before an alignment. Deriving it from the anchor set fails: `min_shared_kmers`
    /// admits thousands of shallow incidental matches, putting `SRR8209082`'s anchor median
    /// at 70x against a true 14,500x.
    ///
    /// Pruning never regresses an accession: if the second alignment fails or spans less of
    /// the reference, the first pass's haplotype stands. The panel run behind the default is
    /// section 6 of `docs/notes/intra_accession_variance.md`.
    pub prune_frac: f64,
    /// The prune's lower bound: no `--prune-frac` setting takes the floor below this.
    ///
    /// Defaults to 5.0. A floor is needed because `prune_frac * depth` on a shallow gene is
    /// a fraction of a read: at 60x, 1% is 0.6, which deletes nothing and means nothing. 5.0
    /// is the smallest coverage this project treats as evidence anywhere.
    pub prune_min_abund: f32,
    /// Raise the variant abundance floor to this fraction of the gene's own depth.
    ///
    /// `0.0` disables it, which is the default and leaves `min_abund` alone. Otherwise the
    /// floor used by `stages::variants` is `max(min_abund, min_abund_frac * depth)`, where
    /// depth is the median `ka:f` over the walked path.
    ///
    /// The absolute `min_abund` means different things at different sequencing depths: five
    /// reads is a 2% floor on a 230x gene and 0.03% on a 15,000x one, where it admits
    /// branches carried by a handful of erroneous reads. `min_vaf` cannot substitute, because
    /// an error branch sits at a respectable fraction of its own *local* coverage; only
    /// against the gene's depth is it visible as noise.
    ///
    /// Depth is taken from the walked path rather than the anchors, which are not a usable
    /// estimate: `min_shared_kmers` admits thousands of shallow incidental matches.
    ///
    /// Left at 0.0 because one accession with genuine variants is not enough to calibrate a
    /// default on. 0.01 is the value to try first.
    pub min_abund_frac: f64,
    /// At a site where two alleles are both real, make the primary carry the more abundant
    /// one. **On by default.**
    ///
    /// Without it the emitted primary carries the *minority* allele at 19% of the two-allele
    /// sites it reports (59 of 311 across the panel); with it, 1 of 311. The invariant to
    /// look for is the maximum alternative VAF, which moves from 0.883 to 0.499: above 0.5
    /// would mean the primary had taken the rarer allele.
    ///
    /// On by default because a consensus carrying the rarer base is simply wrong, and the
    /// corrections replicate: five independent `p_infestans_avr3a` isolates were corrected at
    /// the same two positions, all five having matched the reference exactly beforehand.
    ///
    /// Pruning cannot do this instead: no coverage floor separates a 35% allele from a 65%
    /// one without also deleting every genuine minority variant, so the aligner keeps
    /// choosing by resemblance to the reference.
    ///
    /// It changes the sequence that feeds the codon alignment, so `--no-prefer-majority` is
    /// the switch for reproducing an older result. See `stages::variants::prefer_majority`.
    pub prefer_majority: bool,
    /// Accept alternatives that change the sequence length, not only substitutions.
    ///
    /// **On by default**, but narrowly: see `max_alt_indel` and `allow_frameshift_alts`,
    /// which together admit only the subset that is safe to publish.
    ///
    /// `variants::shifted_window` is what makes them safe: `path_start` and `path_end` are
    /// offsets into the stitched sequence, so reusing them after a splice that adds or
    /// removes bases slices the wrong window.
    pub allow_indel_alts: bool,
    /// Largest length change, in bases, that still counts as an indel allele.
    ///
    /// Measured over the panel, length-changing calls cluster at +-1 and +-2 with a long
    /// tail at +100, +142 and -282. The tail is not an allele: it is the route being sent
    /// through a different stretch of the graph, rewriting the gene rather than editing it.
    /// 30 keeps a ten-codon indel and excludes all of that.
    ///
    /// For `find_alts` this bounds one allele. For `prefer_majority` it also bounds **the
    /// whole correction set** on an accession, because corrections compound: one panel
    /// accession took 64 swaps that were each in frame and each inside this cap, and its
    /// ebola GP primary grew from 2,049 nt to 3,129 against a 2,031 nt reference. Per-swap
    /// legality does not bound the total.
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
    /// Off because the correction loop does not converge with it on. Measured over a
    /// 5,000-accession slice: 97,068 corrections against 470 with it off, a median of 64 per
    /// corrected accession, 7.8 corrections for every net substitution, and a net +33
    /// internal stop codons. Swaps that numerous are undoing one another.
    ///
    /// The reporting is unconditional either way: forks where `_primary` carries the rarer
    /// allele are counted in `minority_left` whether or not they are acted on.
    pub rescue_minority: bool,
    /// Write the literal haplotype sequence into the summary, as 0.1.x did.
    ///
    /// Off by default. The summary carries `seq_sha256` and `length` instead, which answer
    /// what the field was used for without duplicating `reconstructed.fasta`: on real runs
    /// the sequence was 50% to 71% of the file. See `docs/summaries.md`.
    pub summary_seqs: bool,
    /// Cap on alternative haplotypes per accession, `--max-alts`.
    ///
    /// The 5'-to-3' scan runs to the end regardless, counting what did not fit as
    /// `beyond_cap` and warning, so reaching the cap does not silently truncate an accession
    /// to the start of its gene. Real site counts are single digits, so this is a backstop
    /// against a pathological neighbourhood rather than a tuning knob.
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
    /// Density, not the neighbourhood's share of the whole graph: on targeted or amplicon
    /// data the target *is* most of the graph, so a share says more about the library than
    /// about the haplotype. Unitigs per reference base measures how tangled the graph is
    /// around the gene, which is what makes a walk ambiguous. On a 20,000-accession
    /// SARS-CoV-2 run, successes below 3 were 89.7% complete-ORF against 51.3% above.
    ///
    /// `subgraph_fraction` is still reported per result; it is a useful number, not a gate.
    pub subgraph_nodes_per_base_guard: f64,
    /// Refuse an accession whose anchor neighbourhood exceeds this many unitigs.
    ///
    /// The memory term the two-pass parser cannot remove: everything downstream of
    /// `extract::expand` is proportional to the neighbourhood, at roughly 150 bytes per node
    /// in RAM plus about 157 bytes of GFA on scratch. On amplicon-like data the neighbourhood
    /// is most of the graph rather than a corner of it.
    ///
    /// The default sits far above anything real here (PfCRT 747 nodes, ERG11 296, HIV-1
    /// *pol* 34,582), so it bounds the run at about 150 MB of graph per compute worker
    /// without touching any result this repo has produced. 0 disables it.
    pub max_subgraph_nodes: usize,

    // ---- concurrency ----
    /// Threads streaming from S3. Bandwidth saturates around 8 streams at 17.4 MB/s aggregate,
    /// so raising this past that buys nothing.
    pub fetch_workers: usize,
    /// Threads running parse through phase.
    ///
    /// A fixed default of 8, deliberately **not** `available_parallelism()`: auto-sizing
    /// would silently raise the memory bound on a large machine, since each worker holds a
    /// graph. Set it with `-w` once you know the machine can absorb it; the startup line
    /// prints what you chose.
    pub compute_workers: usize,
    /// Bounded queue depth between fetch and compute.
    ///
    /// **Zero, a rendezvous channel.** A positive depth buys nothing here and costs one
    /// capped buffer per slot. A fetch worker that has bytes in hand and no compute worker
    /// ready is already a held buffer, counted by the `fetch_workers` term of
    /// `peak_compressed_buffers`; queueing those bytes separately just moves the same
    /// allocation into a second place and adds a slot's worth of cap to the bound. Measured
    /// over 20 cached 34 MB accessions at `-w 1 --fetch-workers 8`, a depth of 4 cost 130 MiB
    /// of peak RSS, exactly the four slots, and bought no wall clock at all.
    ///
    /// Left configurable rather than removed, because a slow network with a fast compute pool
    /// is the one shape where a slot could help, and the startup line prints what is in force.
    pub channel_capacity: usize,
    /// Threads handed to each GraphAligner invocation.
    ///
    /// **One, for determinism.** GraphAligner's multi-threaded search is not reproducible:
    /// on `SRR8209082`'s subgraph, five identical `-t 2` invocations scored 2686 or 2691 and
    /// produced two different primary haplotypes, where `-t 1` scored 2691 every time.
    ///
    /// It costs nothing to pin: on that graph, the largest here at 34,582 nodes, one thread
    /// takes 0.47s against 0.41s for eight. Throughput comes from `compute_workers` running
    /// whole accessions in parallel, which is the axis that scales.
    pub graphaligner_threads: usize,

    /// Deadline for one GraphAligner call, in seconds. `0` waits forever, which is what every
    /// release before 0.2.0 did.
    ///
    /// A call that never returns holds its compute worker until the scheduler's wall-time
    /// kill, which discards every accession the shard had already finished, since the FASTA,
    /// summary and VCF are written after `pipeline::run` returns. With a deadline it costs
    /// one errored accession.
    ///
    /// **Default 3600 s**, three orders of magnitude above anything measured: the largest
    /// subgraph here (34,582 nodes) aligns in 0.47 s. A stall detector, not a budget, set
    /// where it cannot refuse an alignment that was going to succeed.
    ///
    /// Not in the summary's `settings` block, by that block's own rule: a timeout does not
    /// change what an accession produces, only whether it gets to produce it.
    pub graphaligner_timeout_secs: u64,

    // ---- fetch ----
    /// Whole-request timeout, including body transfer. **One attempt's budget, not one
    /// accession's.**
    ///
    /// Without resumption this would be a hard floor on transfer rate, since every attempt
    /// would restart from byte 0 and time out in the same place. An interrupted transfer
    /// **resumes** (`Range` plus `If-Match`), so a short budget costs nothing and buys earlier
    /// detection of a stall. **Default 300 s**, which at the measured 6.9 MB/s single-stream
    /// rate carries 2 GB per attempt.
    ///
    /// Raising it does not make a stalled connection finish; it only delays noticing.
    pub fetch_timeout_secs: u64,
    /// TCP + TLS connect timeout, kept short so a dead connection fails fast and retries.
    pub fetch_connect_timeout_secs: u64,
    /// Failed attempts per accession before giving up. A 404 means the accession has no Logan
    /// unitigs and is not retried; 408, 429, 5xx and transport errors are.
    ///
    /// **Default 6**, with jittered exponential backoff, which covers about half a minute;
    /// `--fetch-attempts 10` covers about four. A budget shorter than a typical network or
    /// DNS interruption discards one accession per fetch worker per retry cycle, and any
    /// error makes the run non-zero.
    ///
    /// Only **failed** attempts count. One that transferred at least
    /// `RESUME_PROGRESS_BYTES` before the connection dropped is a resume, not a failure, so a
    /// large object on a flaky link is not refused for making slow progress.
    pub fetch_attempts: u32,
    /// Refuse accessions whose compressed size exceeds this.
    ///
    /// This is a safety valve, not a tuning knob, and since 0.2.0 it guards disk and time
    /// rather than memory: nothing in the walker scales with accession size any more, so the
    /// cap no longer appears in the memory bound at all.
    ///
    /// **Default 4 GB**, chosen from the bucket rather than from a machine: in a
    /// 57,000-accession sample of Logan object sizes it refuses 0.16%, against 2.26% at 1 GB,
    /// and past it the gains are negligible while the disk cost keeps doubling. The
    /// distribution is tabulated in `docs/walker.md#choosing-the-size-cap`.
    ///
    /// What it costs is disk, worst case with every in-flight accession at the cap: `16 x cap`
    /// of downloaded objects plus `8 x 0.55 x cap` of spool at the defaults. The startup disk
    /// line prints that bound and the free space. `--unitigs-cache` moves the download term
    /// into the persistent cache term.
    ///
    /// Deliberately **not** the structural ceiling, which the 30-bit unitig id packing puts
    /// near 29 GB: the cap is also a relevance filter, and an 11 GB accession takes about 16
    /// minutes to return `no_path`.
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
    /// How close two references must score before route depth decides between them.
    ///
    /// In multi-anchor mode a sample holding two strains that each equal a panel record gives
    /// two candidates at the same score, and before 0.2.0 the one listed first in the panel
    /// won. A 5% strain could therefore be published as `_primary`, carrying its `anchor=`
    /// tag and every VCF coordinate, with the 95% strain absent from the output entirely and
    /// the result depending silently on panel order.
    ///
    /// Candidates within this much of the best score are treated as tied and the deepest of
    /// them wins. Measured on the metric's own value as a fraction, so 0.002 is two parts in
    /// a thousand whichever `--primary-metric` is in use. `0` turns the depth tie-break off
    /// and restores the pre-0.2.0 ordering.
    pub primary_tie_eps: f64,
    /// Confirm each cache hit against the bucket before using it.
    ///
    /// Logan's `u/` objects are mutable and unversioned, so a warm cache can serve bytes the
    /// bucket has since replaced, and nothing in any output says so: two runs of the same
    /// command then analyse different data and agree in every field that records what was
    /// asked for. This turns that into a detected event.
    ///
    /// A conditional `If-None-Match` GET against the `ETag` the sidecar recorded, so a cache
    /// that is still current costs one round trip and no transfer. Off by default: that is a
    /// round trip a run which trusts its cache should not have to pay.
    ///
    /// Requires `unitigs_cache` and is incompatible with `offline`, which forbids the network
    /// this needs.
    pub revalidate: bool,
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
    /// **Not** `channel_capacity`. Three places hold one concurrently:
    ///
    /// * `channel_capacity` sitting in the bounded queue (zero by default, see the field);
    /// * one per fetch worker, from the moment the body is in hand until the blocking `send`
    ///   into a full queue completes (`stages::stream`);
    /// * one per compute worker, because `run_one` keeps `fetched` alive through the whole
    ///   parse and drops it only once the CSR graph exists (`pipeline::run_one`).
    ///
    /// At the defaults that is 0 + 8 + 8 = 16.
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

    /// Whether a downloaded object is written to disk as it arrives instead of being held.
    ///
    /// True when there is anywhere to put it: `--unitigs-cache` (where it is wanted anyway) or
    /// `--spool-dir` (an unlinked file that dies with the process). `main` always resolves a
    /// spool directory, so this is true for every CLI run; it is false for library callers and
    /// tests that configure neither, who keep the old in-memory behaviour.
    ///
    /// This is the single switch that decides whether memory scales with `--max-size-mb`. When
    /// it is true, nothing in the walker does.
    pub fn downloads_stream_to_disk(&self) -> bool {
        self.unitigs_cache.is_some() || self.spool_dir.is_some()
    }

    /// Worst-case bytes resident in compressed buffers alone.
    ///
    /// With `downloads_stream_to_disk`, a download costs one `DOWNLOAD_BUFFER_BYTES` copy
    /// buffer while it is in flight and nothing afterwards, so this is a constant and the cap
    /// does not appear in it. Without it, every in-flight accession can be a full buffer at
    /// the `--max-size-mb` cap.
    ///
    /// Either way a bound, not a prediction: real accessions are far below the cap
    /// (`ERR2028820` is 33 MB, the ERG11 panel reaches 343 MB against a 4 GB default). It
    /// excludes the parsed graphs (`max_graph_resident_bytes`).
    pub fn max_compressed_resident_bytes(&self) -> u64 {
        self.max_compressed_resident_bytes_for(usize::MAX)
    }

    /// As above, for a run with `n_accessions` of work.
    pub fn max_compressed_resident_bytes_for(&self, n_accessions: usize) -> u64 {
        let per = if self.downloads_stream_to_disk() {
            DOWNLOAD_BUFFER_BYTES
        } else {
            self.max_compressed_bytes
        };
        (self.peak_compressed_buffers_for(n_accessions) as u64).saturating_mul(per)
    }

    /// The budget one parsed graph may use, which is what `stages::parse` enforces.
    ///
    /// `max_graph_resident_bytes` is this times the number of compute workers. Splitting it
    /// out makes the per-accession limit the single definition both the printed bound and the
    /// parser work from, so the promise and the enforcement cannot drift apart.
    pub fn max_graph_bytes_per_accession(&self) -> u64 {
        match self.max_subgraph_resident_bytes() {
            // Adjacency spools and the neighbourhood is capped, so the resident graph is a
            // property of the gene and not of the file. This is the number the startup line
            // leads with, and enforcing the cap-derived one instead would mean promising one
            // bound and policing another about fifty times larger. Plus the flush buffers
            // `scan_in` keeps on the way to the spool, so they do not eat the neighbourhood's
            // allowance; they are cleared on every flush and stay near `SPILL_ADJACENCY_AT`.
            //
            // It governs the degraded path too, deliberately: a spool that cannot be created
            // falls back to adjacency in RAM, and holding that to the same budget refuses the
            // accession as `too_large` rather than letting it grow past what the run
            // announced.
            Some(per) => per.saturating_add(SPILL_HEADROOM),
            // No spool directory, or `--max-subgraph-nodes 0`: nothing bounds the graph by the
            // gene, so the only bound left is the one calibrated against compressed size.
            None => (self.max_compressed_bytes as f64 * GRAPH_PEAK_PER_COMPRESSED_BYTE) as u64,
        }
    }

    /// Whether the graph budget is set by the gene rather than by `--max-size-mb`.
    ///
    /// True for every CLI run, since `main` always resolves a spool directory and
    /// `--max-subgraph-nodes` defaults to 1,000,000. The startup lines branch on it, because
    /// the two cases want different sentences and not just different numbers.
    pub fn graph_bound_is_set_by_the_gene(&self) -> bool {
        self.max_subgraph_resident_bytes().is_some()
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
    ///
    /// Built from `max_graph_bytes_per_accession` rather than from its own copy of the
    /// formula, which is what that function's doc comment always claimed and what makes the
    /// printed bound and the parser's refusal the same number by construction. They were two
    /// expressions of the same arithmetic until the per-accession budget stopped being
    /// `cap x 1.7`, at which point the duplicate would have silently disagreed.
    pub fn max_graph_resident_bytes_for(&self, n_accessions: usize) -> u64 {
        (self.peak_graphs_for(n_accessions) as u64)
            .saturating_mul(self.max_graph_bytes_per_accession())
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

    // ---- disk ----
    //
    // Two consumers, and they do not share a shape. The spool is transient, dies with the
    // accession and scales with `compute_workers`; the cache is persistent, outlives the run
    // and scales with the accession count. A single "disk needed" number would hide that one
    // of them is unaffected by `-w` and the other is unaffected by how long the list is.

    /// Spool bytes one accession at the `--max-size-mb` cap would write.
    ///
    /// The disk twin of `max_graph_bytes_per_accession`, and the better-grounded of the two:
    /// see `SPOOL_BYTES_PER_COMPRESSED_BYTE`. Zero without a `spool_dir`, because then
    /// adjacency never leaves RAM and `max_graph_bytes_per_accession` is the only bound that
    /// applies.
    pub fn max_spool_bytes_per_accession(&self) -> u64 {
        if self.spool_dir.is_none() {
            return 0;
        }
        spool_bytes_for_compressed(self.max_compressed_bytes)
    }

    /// Worst-case spool bytes on disk at once, for a run with `n_accessions` of work.
    ///
    /// One spool per compute worker, each from an accession at the cap. `peak_graphs_for` is
    /// the same count the graph bound uses, so the two move together: a spool exists exactly
    /// when a graph is being parsed, and `--fetch-only` has neither.
    pub fn max_spool_bytes_for(&self, n_accessions: usize) -> u64 {
        (self.peak_graphs_for(n_accessions) as u64)
            .saturating_mul(self.max_spool_bytes_per_accession())
    }

    /// Worst-case bytes of downloaded objects sitting on disk at once, outside the cache.
    ///
    /// With `--unitigs-cache` a download lands in the cache and is counted there, permanently;
    /// this is the other case, where it lands in `--spool-dir` as an unlinked file that lives
    /// only as long as the accession does. One per in-flight accession: this term is
    /// precisely the RAM that streaming moved onto disk, and it is the same size.
    pub fn max_download_bytes_for(&self, n_accessions: usize) -> u64 {
        if self.unitigs_cache.is_some() || !self.downloads_stream_to_disk() {
            return 0;
        }
        (self.peak_compressed_buffers_for(n_accessions) as u64)
            .saturating_mul(self.max_compressed_bytes)
    }

    /// Worst-case bytes `--unitigs-cache` can reach over a run of `n_accessions`, or `None`
    /// when nothing will be written to it.
    ///
    /// Every accession at the cap, which no cohort reaches, and it ignores what the directory
    /// already holds: a bound, like the memory line, not a prediction. The shape is what
    /// matters, and it is `n x cap` with no worker term at all, so lowering `-w` does not
    /// move it. `--offline` writes nothing, so there is no bound to state.
    pub fn max_cache_bytes_for(&self, n_accessions: usize) -> Option<u64> {
        if self.unitigs_cache.is_none() || self.offline {
            return None;
        }
        Some((n_accessions as u64).saturating_mul(self.max_compressed_bytes))
    }
}

/// Spool bytes one accession of `compressed_bytes` is predicted to write.
///
/// Free-standing because the fetch gate works from an accession's real `Content-Length`
/// rather than from the cap, which is the whole point of deciding there: `--max-size-mb`
/// bounds the run, this bounds the accession in hand.
pub fn spool_bytes_for_compressed(compressed_bytes: u64) -> u64 {
    (compressed_bytes as f64 * SPOOL_BYTES_PER_COMPRESSED_BYTE) as u64
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

/// Headroom in the graph budget for the adjacency the parser holds before it spills.
///
/// `scan_in` lets adjacency grow to `SPILL_ADJACENCY_AT` before it creates a spool, and the
/// `Vec`s grow by doubling and are never shrunk, so their capacity can reach twice that and
/// then stays at the high-water mark for the rest of the accession. Four times the threshold
/// leaves a factor of two over that worst case.
///
/// It matters most at a small `--max-subgraph-nodes`, where the neighbourhood's own allowance
/// is the smaller term: at 1,000 nodes that allowance is 150 KB, against roughly 8 MiB of
/// retained flush-buffer capacity, so without this every accession would be refused in pass 1
/// for holding its own buffers.
const SPILL_HEADROOM: u64 = 4 * crate::types::SPILL_ADJACENCY_AT as u64;

/// The copy buffer a download uses on its way from the socket to disk.
///
/// The whole of what one in-flight accession costs in RAM once objects stream. 128 KiB to
/// match `Payload::read_file`, large enough that the syscall count is irrelevant next to the
/// network and small enough that 16 of them do not register against a parsed graph.
pub const DOWNLOAD_BUFFER_BYTES: u64 = 128 * 1024;

/// Spooled adjacency bytes, per byte of compressed unitig file.
///
/// The spool is two files of `u32` and nothing else: one offset per unitig plus a terminator,
/// and one entry per edge. No sequences, no abundances, and no allocator slack, which is why
/// this is a tighter constant than `GRAPH_PEAK_PER_COMPRESSED_BYTE`.
///
/// Measured from 32 MB to 11 GB, the ratio holds between 0.477 and 0.529. The term that
/// drifts is edges per unitig, and it drifts *upward* with size because larger assemblies are
/// branchier, so 0.55 is the largest measured ratio rounded up rather than an average: a
/// constant taken from small accessions would under-predict the large ones this exists to
/// bound. The table is in `docs/walker_reference.md#disk`.
pub const SPOOL_BYTES_PER_COMPRESSED_BYTE: f64 = 0.55;

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
            // Its own setting, so that changing the allele floor does not move the prune's.
            prune_min_abund: 5.0,
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
            graphaligner_timeout_secs: 3600,
            fetch_timeout_secs: 300,
            fetch_connect_timeout_secs: 15,
            fetch_attempts: 6,
            max_compressed_bytes: 4096 * 1024 * 1024,
            scratch_dir: std::env::temp_dir().join("logan-walker-test"),
            graphaligner_bin: String::from("GraphAligner"),
            unitigs_cache: None,
            spool_dir: None,
            offline: false,
            primary_tie_eps: 0.002,
            revalidate: false,
            fetch_only: false,
        }
    }
}

#[cfg(test)]
mod memory_bound_tests {
    use super::{
        spool_bytes_for_compressed, Config, DOWNLOAD_BUFFER_BYTES,
        GRAPH_PEAK_PER_COMPRESSED_BYTE,
    };
    use std::path::PathBuf;

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

    /// Written against the cap rather than against a number of gigabytes: the default has
    /// moved twice and the claim is about the multiplier.
    #[test]
    fn worst_case_bytes_track_the_size_cap() {
        let cfg = Config::default();
        assert_eq!(cfg.max_compressed_bytes, 4096 * 1024 * 1024, "the 0.2.0 default");
        // `Config::default()` configures no cache and no spool dir, so downloads are still
        // held and the cap is still the per-buffer size. A CLI run is the streaming case.
        assert!(!cfg.downloads_stream_to_disk());
        assert_eq!(
            cfg.max_compressed_resident_bytes(),
            cfg.peak_compressed_buffers() as u64 * cfg.max_compressed_bytes
        );

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
        // The disk bounds multiply the same cap, so they need the same protection. 0.55 of
        // u64::MAX does not itself overflow a f64-to-u64 cast, but the worker and accession
        // multiplications do.
        let cfg = Config {
            max_compressed_bytes: u64::MAX,
            spool_dir: Some(PathBuf::from("/tmp")),
            unitigs_cache: Some(PathBuf::from("/tmp")),
            ..Config::default()
        };
        assert_eq!(cfg.max_spool_bytes_for(usize::MAX), u64::MAX);
        assert_eq!(cfg.max_cache_bytes_for(1_000_000), Some(u64::MAX));
    }

    /// The graph budget is set by the gene, and the printed bound is the budget.
    ///
    /// Both halves matter: a printed bound the parser does not enforce is a calibration, and
    /// an enforced bound that is not printed cannot be planned against. One expression.
    #[test]
    fn the_graph_budget_is_set_by_the_gene_and_is_what_gets_printed() {
        let gene = Config {
            spool_dir: Some(PathBuf::from("/tmp")),
            ..Config::default()
        };
        assert!(gene.graph_bound_is_set_by_the_gene());
        // Raising the size cap a hundredfold must not move it.
        let big = Config { max_compressed_bytes: gene.max_compressed_bytes * 100, ..gene.clone() };
        assert_eq!(
            big.max_graph_bytes_per_accession(),
            gene.max_graph_bytes_per_accession(),
            "the graph budget must not mention the size cap once adjacency spools"
        );
        // It is the neighbourhood cap that moves it.
        let narrow = Config { max_subgraph_nodes: 1000, ..gene.clone() };
        assert!(narrow.max_graph_bytes_per_accession() < gene.max_graph_bytes_per_accession());
        // Pass 1 holds adjacency in RAM until it reaches `SPILL_ADJACENCY_AT`, and the Vecs
        // keep that capacity for the rest of the accession, so the budget has to leave room
        // for it however small the neighbourhood cap is. Without the headroom,
        // `--max-subgraph-nodes 1000` would budget 150 KB against ~8 MiB of retained buffers
        // and refuse every accession in pass 1.
        let tiny = Config { max_subgraph_nodes: 1, ..gene.clone() };
        assert!(
            tiny.max_graph_bytes_per_accession()
                > 2 * crate::types::SPILL_ADJACENCY_AT as u64,
            "the budget must clear the flush buffers pass 1 holds before it spills, even at \
             the smallest neighbourhood cap"
        );

        // The printed bound is the budget times the graphs, not a second copy of a formula.
        assert_eq!(
            gene.max_graph_resident_bytes_for(usize::MAX),
            gene.peak_graphs_for(usize::MAX) as u64 * gene.max_graph_bytes_per_accession()
        );

        // Without a spool directory, or with the neighbourhood cap off, nothing bounds the
        // graph by the gene and the budget falls back to the compressed size.
        for cfg in [
            Config::default(),
            Config { max_subgraph_nodes: 0, ..gene.clone() },
        ] {
            assert!(!cfg.graph_bound_is_set_by_the_gene());
            assert_eq!(
                cfg.max_graph_bytes_per_accession(),
                (cfg.max_compressed_bytes as f64 * GRAPH_PEAK_PER_COMPRESSED_BYTE) as u64
            );
        }
    }

    /// Streaming downloads is what takes `--max-size-mb` out of the memory bound.
    ///
    /// Before it, every in-flight accession could be a full buffer at the cap, which is why
    /// the cap could not be raised: at 28 GB the defaults would have promised 437 GB of RAM,
    /// and meant it. The test that matters is the last one, that raising the cap a
    /// hundredfold does not move the number at all.
    #[test]
    fn streaming_downloads_take_the_cap_out_of_the_memory_bound() {
        let held = Config::default();
        assert!(!held.downloads_stream_to_disk(), "neither a cache nor a spool dir");
        assert_eq!(
            held.max_compressed_resident_bytes(),
            held.peak_compressed_buffers() as u64 * held.max_compressed_bytes
        );

        for cfg in [
            Config { spool_dir: Some(PathBuf::from("/tmp")), ..Config::default() },
            Config { unitigs_cache: Some(PathBuf::from("/tmp")), ..Config::default() },
        ] {
            assert!(cfg.downloads_stream_to_disk());
            assert_eq!(
                cfg.max_compressed_resident_bytes(),
                cfg.peak_compressed_buffers() as u64 * DOWNLOAD_BUFFER_BYTES
            );
            // 16 x 128 KiB, against 16 GB before.
            assert!(cfg.max_compressed_resident_bytes() < 4 * 1024 * 1024);
            let hundredfold = Config {
                max_compressed_bytes: cfg.max_compressed_bytes * 100,
                ..cfg.clone()
            };
            assert_eq!(
                hundredfold.max_compressed_resident_bytes(),
                cfg.max_compressed_resident_bytes(),
                "the cap must not appear in the memory bound once objects stream"
            );
        }
    }

    /// A downloaded object is on disk in exactly one place, and counted in exactly one term.
    ///
    /// With a cache it lands there and stays; without one it is an unlinked file that dies
    /// with the accession. Counting it in both would double the printed disk bound, and
    /// counting it in neither would hide the bytes that streaming moved off the heap.
    #[test]
    fn a_downloaded_object_is_counted_once_and_in_the_right_term() {
        let dir = || Some(PathBuf::from("/tmp"));
        let cap = Config::default().max_compressed_bytes;

        let spool_only = Config { spool_dir: dir(), ..Config::default() };
        assert_eq!(
            spool_only.max_download_bytes_for(50),
            spool_only.peak_compressed_buffers_for(50) as u64 * cap
        );
        assert_eq!(spool_only.max_cache_bytes_for(50), None);

        let cached = Config { spool_dir: dir(), unitigs_cache: dir(), ..Config::default() };
        assert_eq!(cached.max_download_bytes_for(50), 0, "already counted as cache");
        assert_eq!(cached.max_cache_bytes_for(50), Some(50 * cap));

        let neither = Config::default();
        assert_eq!(neither.max_download_bytes_for(50), 0, "held in RAM, counted as memory");
    }

    /// `SPOOL_BYTES_PER_COMPRESSED_BYTE` is a bound, so it must sit above every accession it
    /// was measured on, and close enough above them to be worth stating.
    ///
    /// The three span a 340x size range. The prediction has to cover the largest ratio,
    /// because an under-prediction is a gate that lets through an accession that then fills
    /// the disk; it must not cover it by so much that the gate refuses feasible runs.
    #[test]
    fn the_spool_constant_bounds_every_accession_it_was_measured_on() {
        // (compressed bytes, spool bytes) as measured.
        let measured: [(u64, u64); 3] = [
            (32_400_000, 15_438_276),          // ERR2028820
            (11_023_800_656, 5_626_407_158),   // DRR240330
            (10_989_106_976, 5_808_943_268),   // DRR212963
        ];
        for (compressed, spool) in measured {
            let predicted = spool_bytes_for_compressed(compressed);
            assert!(
                predicted >= spool,
                "{compressed} compressed: predicted {predicted} under the measured {spool}"
            );
            assert!(
                predicted as f64 <= spool as f64 * 1.2,
                "{compressed} compressed: predicted {predicted} over the measured {spool} by \
                 more than 20%"
            );
        }
    }

    /// The two disk consumers do not share a shape, and the startup line says so because the
    /// numbers do. The spool moves with `-w` and is capped by the pipeline; the cache moves
    /// with the list length and ignores `-w` entirely.
    #[test]
    fn the_spool_scales_with_workers_and_the_cache_with_the_list() {
        let base = Config {
            spool_dir: Some(PathBuf::from("/tmp")),
            unitigs_cache: Some(PathBuf::from("/tmp")),
            ..Config::default()
        };
        let wide = Config { compute_workers: 16, ..base.clone() };
        // Workers move the spool...
        assert!(wide.max_spool_bytes_for(usize::MAX) > base.max_spool_bytes_for(usize::MAX));
        // ...and do not move the cache.
        assert_eq!(wide.max_cache_bytes_for(50), base.max_cache_bytes_for(50));
        // The list moves the cache without limit...
        assert!(base.max_cache_bytes_for(1000) > base.max_cache_bytes_for(50));
        // ...and moves the spool only until the pipeline is full.
        assert_eq!(base.max_spool_bytes_for(1000), base.max_spool_bytes_for(8));
        assert!(base.max_spool_bytes_for(2) < base.max_spool_bytes_for(8));
        // At 8 workers and a 50-accession list: 4.4 x cap of transient spool against 50 x
        // cap of cache. Summing them would hide which one `-w` touches.
        let cap = base.max_compressed_bytes as f64;
        let spool = base.max_spool_bytes_for(50) as f64;
        assert!((spool - 8.0 * 0.55 * cap).abs() / spool < 1e-6);
        assert_eq!(base.max_cache_bytes_for(50), Some(50 * base.max_compressed_bytes));
    }

    /// No spool directory means adjacency never leaves RAM, where the graph budget already
    /// bounds it. Claiming disk for it would have the startup line reserve space nothing
    /// writes, and would arm a fetch gate against a spool that cannot exist.
    #[test]
    fn without_a_spool_directory_there_is_no_disk_bound() {
        let cfg = Config::default();
        assert!(cfg.spool_dir.is_none());
        assert_eq!(cfg.max_spool_bytes_per_accession(), 0);
        assert_eq!(cfg.max_spool_bytes_for(usize::MAX), 0);
        // And no cache means nothing is written to one, `--offline` most of all: it is the
        // mode that reads a cache and never adds to it.
        assert_eq!(cfg.max_cache_bytes_for(50), None);
        let offline = Config {
            unitigs_cache: Some(PathBuf::from("/tmp")),
            offline: true,
            ..Config::default()
        };
        assert_eq!(offline.max_cache_bytes_for(50), None);
    }

    /// `--fetch-only` parses nothing, so it spools nothing, but filling the cache is the
    /// whole mode: the one case where the persistent term exists and the transient one does
    /// not.
    #[test]
    fn fetch_only_spools_nothing_and_caches_everything() {
        let cfg = Config {
            fetch_only: true,
            spool_dir: Some(PathBuf::from("/tmp")),
            unitigs_cache: Some(PathBuf::from("/tmp")),
            ..Config::default()
        };
        assert_eq!(cfg.max_spool_bytes_for(usize::MAX), 0);
        assert_eq!(cfg.max_cache_bytes_for(50), Some(50 * cfg.max_compressed_bytes));
    }

    /// The parsed graphs come on top of the buffers: one per compute worker, up to 1.7x the
    /// cap each.
    #[test]
    fn the_bound_includes_the_parsed_graphs() {
        let cfg = Config::default();
        let cap = cfg.max_compressed_bytes as f64;
        let graphs = cfg.max_graph_resident_bytes() as f64;
        assert!((graphs - 8.0 * 1.7 * cap).abs() / graphs < 1e-6);
        // Buffers plus graphs: 16 x cap + 8 x 1.7 x cap = 29.6 x cap. It was 33.6 until the
        // queue depth went to 0.
        let total = cfg.max_resident_bytes() as f64;
        assert!((total - 29.6 * cap).abs() / total < 1e-6);
        let more = Config { compute_workers: 16, ..Config::default() };
        assert!(more.max_graph_resident_bytes() > cfg.max_graph_resident_bytes());
    }

    /// The pipeline's shape is a ceiling for a long list, not a figure to quote at a short
    /// one: a single-accession run should not be told to reserve the whole pool's worth.
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
        // 1 x cap + 1 x 1.7 x cap, against 29.6 x cap for a list long enough to fill the
        // pipeline.
        let one = cfg.max_resident_bytes_for(1) as f64;
        assert!((one - 2.7 * cfg.max_compressed_bytes as f64).abs() / one < 1e-6);
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
