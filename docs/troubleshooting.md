# Troubleshooting

Failure modes by symptom, with what to do about each. For what an exit code means in general,
see [`walker_reference.md`](walker_reference.md#exit-status) for the walker and
[`codon_msa.md`](codon_msa.md) for the composites.

## The run exits non-zero

| exit | meaning | what to do |
|---|---|---|
| `1` | nothing was reconstructed, or a usage or I/O error | read `<query stem>.skipped.tsv`: it names every accession that produced nothing and why |
| `2` | the command line was not understood | argparse's own code; the message names the argument |
| `3` | every accession errored | the run did not work. `--allow-empty` does not apply, and the cause is in the per-accession `detail` |
| `4` | some accessions errored and others reconstructed | the haplotypes are real and were written; the failures are in `skipped.tsv` |
| `75` | a Galaxy dataset was not ready (`queued`, `new`, `running`, `setting_metadata`) | wait for the dataset and run again. 75 is `EX_TEMPFAIL`, so a caller can retry it; `error` and `discarded` are not retryable and are not reported as 75 |
| `127` | a required external binary could not be found | install it, or name it: `--walker-bin`, `--sieve-bin`, or the walker's own `--graphaligner`. MAFFT and FastTree are taken from `PATH` |
| `141` | standard output was closed early, as in a pipeline ending in `head` | expected; the walker removes its scratch directory on the way out |

## A tool is missing or the wrong version

- **`GraphAligner` not found.** The walker checks it at startup, before downloading anything.
  Put 1.0.17b on `PATH` or name it with `--graphaligner`.
- **GraphAligner is not 1.0.17b.** The run warns and continues. 1.0.19 and later drop the last
  base of every query, so haplotypes come out one base short; see
  [`notes/graphaligner_end_clipping.md`](notes/graphaligner_end_clipping.md).
- **MAFFT 7.525.** `msa` refuses it by name: 7.525 has a bug in `--addfragments`, the path
  `msa` uses. Install `mafft >=7.526`. Channel order decides which you get, so list
  conda-forge before bioconda, or `conda install 'mafft>=7.526'`.
- **The sieve binary dies with SIGILL.** Its Makefile builds with `-march=native`, so a binary
  copied from another machine is illegal on this CPU. Rebuild it here:
  `(cd extern/kmer-diversity-sieve && make serial)`.
- **`logan-walker` is too old.** `disassembler` runs `logan-walker --version` before using it
  and refuses anything below `MIN_WALKER_VERSION`. Update the package, or check what
  `--walker-bin` and `DISASSEMBLER_WALKER_BIN` point at.

## The walker runs but finds nothing

Every accession that produced no haplotype is one row of `<out-dir>/<query stem>.skipped.tsv`,
with its outcome and detail. The outcomes that are not errors:

- `no_anchors`: no unitig shares enough k-mers with the reference. The gene is not in that
  accession, or the reference is too divergent. Lowering `--min-shared-kmers` or `-k` widens
  the net.
- `no_path`: the gene is there but no walk spans `--min-cov` of the reference. Partial
  coverage, a tangled neighbourhood, or a genuinely truncated gene.
- `no_links`: the Logan object carries no `L:` adjacency at all, so nothing could be walked
  across unitigs. A property of the object, not of the sample.
- `no_logan_data`: Logan has no unitigs for the accession.
- `too_large` / `subgraph_too_large`: over `--max-size-mb`, over the parsed graph budget, over
  the free space under `--spool-dir`, or over `--max-subgraph-nodes`. The detail says which.

## Resource and environment problems

- **`/dev/shm` is 64 MB in a default Docker container**, and that is where the walker puts its
  scratch GFA and GAF. Point `--scratch-dir` at a real directory inside the container.
- **The run is using more memory than expected.** The walker prints its bound at startup. It
  is `--max-subgraph-nodes` times about 150 bytes per compute worker, so lower `-w` and
  `--fetch-workers`; they multiply. See
  [`walker_reference.md`](walker_reference.md#memory).
- **The disk filled.** Downloads and the adjacency spool both go under `--spool-dir`, which
  defaults to `--unitigs-cache` when set and the system temp directory otherwise. A write that
  runs out of space is reported as `too_large` for that accession rather than failing the run.
  Point `--spool-dir` somewhere with room; do not point it at `--scratch-dir`, which is tmpfs.
- **A `logan-walker.<pid>.<suffix>` directory was left behind.** A run killed outright cannot
  clean up; the next run with the same `--scratch-dir` sweeps it. Directories left by 0.1.x
  carry no lock and are never swept, so clear those once by hand.

## Cache problems

- **`parse failed: Unknown frame descriptor` on a cached accession.** The cached file is not a
  zstd object, usually a proxy's error page saved under the accession's name. The walker now
  detects this, removes the entry and downloads it again; under `--offline` it refuses and
  names the file for you to inspect.
- **A cache miss is reported as an error under `--offline`.** Either the accession was never
  fetched, or its `.unitigs.absent` note was deleted. Refill the cache with
  `disassembler fetch`.
- **`revalidated: replaced` in the summary.** Logan's `u/` objects are mutable, and the bucket
  now serves different bytes than the cache held. This run used the new object; an earlier run
  of the same command did not. See [`walker.md`](walker.md#--revalidate).

## A run was killed

A run that dies leaves `<query stem>.partial.fasta` and `<query stem>.journal.jsonl` in the
output directory. The FASTA holds every haplotype finished before the run died, byte-identical
to what the final file would have carried; the journal has one line per finished accession. The
alternatives, the VCFs, `bubbles.json` and `summary.json` come only from a run that reaches the
end. Use the journal's `acc` column as the list of accessions a rerun can leave out. Both files
are removed on a clean finish, so finding them means a run is in flight or one died here.

## The results look wrong

- **Every frequency is empty and the prune did nothing.** The unitig file may carry no
  abundance field. The summary reports this as `abundance_missing`, per record and as a cohort
  count. Without abundance, every depth-based stage is inert.
- **`mixed=1` on a primary.** The route crossed two or more two-allele sites, so the published
  sequence is a combination nothing in the data links. The alleles it chose between are in
  `alt_haplotypes.fasta` and the VCF. See [`walker.md`](walker.md#what-the-primary-is).
- **`bias_uncorrected` is true.** The prune ran and its result was thrown away, so the
  first-pass route stands with the aligner's reference pull still in it. Common on amplicon
  data with uneven coverage. See
  [`notes/intra_accession_variance.md`](notes/intra_accession_variance.md).
- **A frequency looks too high.** `ka:f` is clipped at 25,000 in the Logan files, so a fork
  whose major branch is clipped reports every minor allele as more frequent than it is. Those
  rows carry `KaSat` in the VCF FILTER column and `ka_saturated` in `snvs.tsv`.
- **The sieve dropped most of the cohort.** Its threshold is a divergence, not a count, and on
  tiler rows the gap limit acts as a second coverage gate. See [`sieve.md`](sieve.md).
