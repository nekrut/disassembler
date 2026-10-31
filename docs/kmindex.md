# The kmindex route: harvest and stratify

A [kmindex](https://github.com/tlemane/kmindex) query of a reference CDS against the Logan
planetary indexes returns, for every index shard, the SRA accessions that contain enough of the
query's k-mers, with a containment score. On usegalaxy.org that query runs across about 2,869
shards and produces one JSON file per shard. `disassembler harvest` merges those shards into
per-target accession lists, `disassembler stratify` optionally cuts a smaller cohort from them by
score, and the walker reconstructs the cohort (see [`walker.md`](walker.md)).

## `disassembler harvest`

Shards come from one of two places, which produce identical output:

```bash
# A Galaxy dataset collection, over the Galaxy API
disassembler harvest --history-id <HISTORY_ID> --collection-id <COLLECTION_ID> -o candidates/

# A local directory of shard files
disassembler harvest --shards-dir path/to/shards/ -o candidates/
```

In Galaxy API mode, `--download-workers` (default 6) sets the number of parallel downloads. Six
is about the ceiling before usegalaxy.org's rate limiter starts returning HTTP 429; more still
completes, but each throttled shard waits 30 to 60 seconds. Raise it for a self-hosted Galaxy.
`--base-url` selects the server, and the API key comes from `--api-key-file`,
`$GALAXY_API_KEY` or `/tmp/gxy.txt`.

In local mode, every regular file directly in `--shards-dir` is read as a shard, whatever its
name (hidden files are skipped). A file that is not valid JSON counts as a failed shard and is
named in the log. Shards, local or fetched from Galaxy, can be gzip-compressed.

Outputs, written to `-o`:

| file | contents |
|---|---|
| `<target>.accessions.txt` | every accession the query returned for that target, by containment score, highest first |
| `<target>.scores.json` | `{"target": ..., "scores": {accession: containment score}}` |
| `UNION.accessions.txt` | the union across all targets |
| `summary.json` | per-target counts and maximum score |

A target is one record of the kmindex query, named by the first word of its FASTA header, so a
multi-record query gives one set of files per record.

## `disassembler stratify`

`stratify` samples a cohort from harvest score files, by score bin:

```bash
# 30 accessions: 10 from each of three low-score bins
disassembler stratify -o candidates/ \
    --score-bin 0.0:0.05:10 --score-bin 0.05:0.10:10 --score-bin 0.10:0.15:10
```

Each `--score-bin LO:HI:N` selects up to N accessions whose score falls in `[LO, HI)`, spread
evenly across the sorted range.

- **Half-open bins.** Bins can touch without double-counting, and overlapping bins are refused.
  The exception is the top: a bin ending at 1.0 includes 1.0, the highest containment score, so
  `0.9:1.0:N` takes perfect matches too. Before 0.1.1 it did not: in the SARS-CoV-2 spike
  harvest, 274,808 accessions (5.3%) score exactly 1.0, and none could reach the 20,000-accession
  walker cohort drawn with `0.90:1.00`.
- **Nonzero scores only.** Zero-score accessions are excluded, so `0.0:0.05:10` samples from
  `(0, 0.05)`.
- **No backfill.** A bin with fewer than N candidates contributes what it has, and no other bin
  is topped up to compensate.

**Finding the score files.** Each score file names its own target, and `stratify` recognises
score files by their content, never their names (a Galaxy dataset arrives as `dataset_<id>.dat`).
By default it reads every file in the `-o` directory and uses the score files among them;
`--scores PATH` (repeatable) names files or directories instead, plain or gzip-compressed.
`--target NAME` stratifies only that target. A score file written before 0.1.1, or by hand, is a
bare `{accession: score}` with no target inside; name it with `--scores NAME=FILE`, or with
`--target` when it is the only one. Finding no score file is an error, where it used to print a
warning and exit 0.

Outputs, per target, next to the harvest files: `<target>.bins.txt` (the merged list, by score,
highest first) and `<target>.bins.json` (per-bin requested, available and selected counts, and
the accession lists), plus a combined `bins_summary.json`. The list feeds straight into
`disassembler reconstruct -l <target>.bins.txt`.

If download volume is the real constraint, sampling by count is the wrong lever. Unitig file
sizes are very skewed: across 45 randomly sampled accessions the median was 0.6 MB, the mean
17.6 MB and the largest 595 MB, and the top 10% of files carried 90% of the bytes. The walker's
`--max-size-mb` targets the bytes directly, by rejecting an accession on its `Content-Length`
before downloading anything.

## `disassembler from-kmindex`

`from-kmindex` chains harvest, an optional stratify (when any `--score-bin` is given), and
`reconstruct` (walker, optional sieve, codon MSA and tree) for one target:

```bash
disassembler from-kmindex \
    --history-id <HISTORY_ID> --collection-id <COLLECTION_ID> \
    --target zika_ns5 -q queries/markers/zika_ns5.fasta -o out/zika_ns5
```

The shards come from exactly one of three places:

| flags | harvest |
|---|---|
| `--history-id` and `--collection-id` | fetches the shard collection from Galaxy |
| `--shards-dir <DIR>` | reads shard JSON files already on disk, as `harvest --shards-dir` does |
| `--skip-harvest` | none: reuses an earlier harvest in `--cohort-dir` |

Harvest output goes to `--cohort-dir` (default `<out-dir>/cohort`), and `--target` picks which
target's list is reconstructed. `--harvest-workers` is harvest's `--download-workers`. The
walker, sieve and MSA options are the same as `reconstruct`'s (see
[`codon_msa.md`](codon_msa.md#disassembler-reconstruct)). For a run that never touches the
network, harvest from `--shards-dir` (or `--skip-harvest`) and walk with `--offline` and a
pre-filled `--unitigs-cache`:

```bash
disassembler from-kmindex --skip-harvest --cohort-dir out/zika_ns5/cohort \
    --target zika_ns5 -q queries/markers/zika_ns5.fasta -o out/zika_ns5 \
    --unitigs-cache cache/ --offline
```

Before 0.1.1, `from-kmindex` crashed in every harvest with an `AttributeError` (it did not pass
harvest's `--shards-dir` setting), and it required the two Galaxy IDs even with
`--skip-harvest`. Command lines written for that, with placeholder IDs, still run: under
`--skip-harvest` the IDs are ignored with a warning. With `--score-bin`, `from-kmindex` passes its
`--target` to `stratify`, so a hand-made score file in the older shape works without a target
inside.

## Scale

A popular gene can return millions of candidates: a SARS-CoV-2 spike search returned about
5.1 million accessions. Walking all of them is not practical, so harvest, stratify (or
otherwise subset) the list, and walk the subset. On Galaxy the walker is scattered over the list
in shards; see [`galaxy.md`](galaxy.md#scatter).
