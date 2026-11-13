# Walker: fetch, cache, inputs, scratch directory and signals (stream.rs, scratch.rs, io/input.rs, error.rs)

The core of this area works: cache entries get their final name only through a temp-file rename, the size cap is enforced both by Content-Length and by a body limit, every run gets its own 0700 scratch directory, and zstd frame checksums catch corrupted bytes. The weak spots are the scratch sweep and the signal handler, both new in 0.1.1. Together they account for seven of the entries below. The sweep deletes live, kept or user-owned directories. The handler overrides `nohup`, misses several fatal signals, leaves child processes running, and crashes when the walker is PID 1. The cache has a loose trust model: it records nothing about where its bytes came from, does not check that an entry is zstd or belongs to its accession, records no negative results, and treats a failed write as success. Fetching uses fixed policies the user cannot change: a 1.5 s retry budget, a 900 s whole-request timeout that throws away partial downloads, only the CA roots compiled into the binary, and redirects followed to any host. Most fixes share one idea: decide on evidence the walker itself produced (a lock or marker file it wrote, an entry's magic bytes and source, a recorded negative, a signal disposition it inherited), not on guesses from names, PIDs or whether a file exists.

| ID | Severity | Where | Title | Fix effort |
|---|---|---|---|---|
| WF-01 | high | `walker/src/scratch.rs:151` | Scratch sweep judges liveness by PID and deletes live runs' directories across PID namespaces or hosts | medium |
| WF-02 | high | `walker/src/scratch.rs:73` | Walker overrides an inherited SIG_IGN, so a `nohup`'d run dies on hangup (0.1.1 regression) | trivial |
| WF-03 | medium | `walker/src/stages/stream.rs:110` | Fetch then `--offline` turns `no_logan_data` / `too_large` accessions into errors (exit 1 or 3) | small |
| WF-04 | medium | `walker/src/scratch.rs:139` | Scratch sweep recursively deletes any user directory named `logan-walker.<integer>[.*]` | trivial |
| WF-05 | medium | `walker/src/io/output.rs:245` | summary.json records no input identity; cache files are used without any identity check | medium |
| WF-06 | medium | `walker/src/stages/stream.rs:128` | `--fetch-only` reports a failed cache write as "cached" and exits 0 | small |
| WF-07 | medium | `walker/src/stages/stream.rs:161` | 1.5 s fixed retry budget turns brief outages into per-accession errors; 408/429 never retried | small |
| WF-08 | medium | `walker/src/stages/stream.rs:92` | Compressed unitig files held whole in RAM in 4+F+W places; cache hits copied needlessly | small |
| WF-09 | medium | `walker/src/stages/stream.rs:49` | Fixed 900 s whole-request fetch timeout discards partial downloads and restarts from byte 0 | small |
| WF-10 | medium | `walker/src/stages/stream.rs:48` | Fetch trusts only compiled-in Mozilla roots; `SSL_CERT_FILE` ignored behind TLS-inspecting networks | small |
| WF-11 | medium | `walker/src/main.rs:460` | Query FASTA read twice: piped multi-record panel refused, duplicate warnings, two parsers | small |
| WF-12 | low | `walker/src/io/input.rs:107` | Lowercase SRA accession is fetched verbatim, 404s, and is recorded as `no_logan_data` | trivial |
| WF-13 | low | `walker/src/stages/stream.rs:96` | Unitigs cache accepts and trusts non-zstd entries; a bad entry fails forever, even online | small |
| WF-14 | low | `walker/src/scratch.rs:151` | `LOGAN_KEEP_SCRATCH` directory is deleted by the next run's sweep | trivial |
| WF-15 | low | `walker/src/scratch.rs:80` | Signal handler removes scratch and exits without stopping running GraphAligner children | small |
| WF-16 | low | `walker/src/stages/stream.rs:92` | `--fetch-only` reads every cached file whole, validates nothing, counts hits as fetched | small |
| WF-17 | low | `walker/src/scratch.rs:73` | Scratch directory left behind when the run dies by SIGPIPE, SIGQUIT, SIGUSR1/2 or SIGXCPU | small |
| WF-18 | low | `walker/src/error.rs:75` | Runtime failures labelled "usage error"; non-zstd payloads "parse failed"; non-UTF-8 GAF fails every accession | small |
| WF-19 | low | `walker/src/stages/stream.rs:262` | Killed cache writer leaves orphan hidden `.tmp` files; entry not fsynced before rename | small |
| WF-20 | low | `walker/src/scratch.rs:80` | As PID 1 of a container, SIGTERM/SIGINT ends the walker with SIGSEGV (exit 139) | trivial |
| WF-21 | low | `walker/src/stages/stream.rs:48` | Unitig fetcher follows redirects to any origin, including plain http | trivial |
| WF-22 | low | `walker/src/io/input.rs:43` | Inputs read with no size limit (gzip inflated whole); accession length unbounded | small |

### WF-01 · high · Scratch sweep judges liveness by PID and deletes live runs' directories across PID namespaces or hosts

`walker/src/scratch.rs:151` · concurrency · confirmed (reproduced) · fix effort: medium · sources: RA-01, BL2-05, IC-04, R1-01, RB-02, RE-01, RX-04

**Problem.** At startup `sweep_stale` reads the PID out of every `logan-walker.<pid>.*` name in the scratch base. It then runs `if !process_exists(pid) { remove(&entry.path());` (`walker/src/scratch.rs:151`). `process_exists` is `libc::kill(pid, 0)` (`walker/src/scratch.rs:169`), which sees only processes in the caller's PID namespace on the caller's host. A live run in another container, or on another node that shares the base, therefore looks dead, and its directory is deleted. The victim run does not notice: `extract.rs:293` silently recreates the directory with `create_dir_all`, and `align.rs:231` maps a vanished GAF (`ErrorKind::NotFound`) to `Ok(HashMap::new())`, which is recorded as a biological `no_path`.

**Impact.** The trigger is runs of the same user that share a `--scratch-dir` but not a PID namespace. Examples: array jobs on several nodes with scratch on Lustre or NFS, or containers that bind-mount one volume (`docs/walker.md:274` tells Docker users to move `--scratch-dir`). Each new run deletes the in-flight files of all the others. Accessions that are being aligned at that moment are reported as `no_path`, with no warning and exit 0. Not affected: the default `/dev/shm`, the Galaxy wrapper (one scratch directory per job), Apptainer (which uses the host PID namespace), and Docker with the walker as PID 1. A related inverse bug: a SIGKILLed run that was PID 1 leaves `logan-walker.1.*`, which is never swept. A second reviewer rated this medium because the default setups are safe. It stays high because the wrong record is silent (reproduced twice) and a shared scratch path on a cluster filesystem is an ordinary HPC setup.

**Evidence.** Run A used 3 accessions, `-w 1`, and a stand-in aligner that opens its GAF and then works for 4 s. Run B started 2 s later on the same `--scratch-dir` under `unshare --pid --fork --mount-proc`. B printed `removed 1 scratch directory left in ... by runs that were killed`. A's directory reappeared with mode 755 (it was 700). A reported `[1/3] ...: no_path` and `2/3 accessions reconstructed`, and exited 0. In an independent reproduction, 4 of 8 accessions came out `no_path` with no error line and exit 0; the control run gave 8 of 8.

**Proposed fix.** Decide liveness with a kernel lock instead of a PID:
1. In `Scratch::create`, create `<dir>/.lock` with `create_new` and take `flock(LOCK_EX|LOCK_NB)` on it. Then write the kernel's boot id into it (`/proc/sys/kernel/random/boot_id`, or the hostname if that is missing), and keep the `File` open in `Scratch`. Under `LOGAN_KEEP_SCRATCH`, also write `.keep` (WF-14).
2. In `sweep_stale`, consider only names of the form `logan-walker.<digits>.<9 hex>`, which is what `scratch.rs:46` builds (WF-04). Remove a directory only if its `.lock` exists, names this kernel's boot id, there is no `.keep`, and a non-blocking `flock` succeeds. Delete `process_exists`.
3. In the `align::run` NotFound branch, return `Error::GraphAligner` ("removed while GraphAligner ran") when `subgraph.gfa_path` or `query_fa` no longer exists, instead of `no_path`.
4. Stop recreating the directory in production: make the `create_dir_all` at `extract.rs:293` test-only, or keep it (unit tests call `write_gfa` directly) and rely on step 3.

With this patch the reproduction gives 3 of 3 reconstructed and B removes nothing. SIGKILLed directories are still swept, including the PID-1 case. An external `rm -rf` during alignment now exits 3 with an error where 0.1.1 recorded `no_path`. On filesystems where `flock` is missing or works only within one node (Lustre `noflock` or `localflock`, NFS `nolock`), the boot-id check makes the run fail safe: directories leak and are never wrongly deleted. Directories left by 0.1.1 have no `.lock` and are never swept, so the release notes should ask users for a one-time manual cleanup.

**Tests and docs.** Unit tests: a directory whose lock is held survives even under PID 99999999; a lock that names another kernel survives; a free lock is swept even under PID 1. Add a `tests/cli.rs` case that runs the second walker under `unshare --pid` and skips when `unshare` fails. Update `docs/walker.md:266-268` (ended runs on this machine are found by lock, not PID, and other hosts are never touched), `docs/walker.md:274` (a shared base is safe), the module doc at `scratch.rs:9-13`, and the comment at `walker/Cargo.toml:17`.

### WF-02 · high · Walker overrides an inherited SIG_IGN, so a `nohup`'d run dies on hangup (0.1.1 regression)

`walker/src/scratch.rs:73` · reliability · confirmed (reproduced) · fix effort: trivial · sources: RE-02, R1-03, RB-01

**Problem.** `Scratch::remove_on_signal`, called from `main.rs:542`, registers handlers without looking at the dispositions the process inherited: `let mut signals = signal_hook::iterator::Signals::new([SIGINT, SIGTERM, SIGHUP])?;` (`walker/src/scratch.rs:73`). When a signal arrives, the handler thread removes the scratch directory and calls `emulate_default_handler(sig)` (`walker/src/scratch.rs:80`), which kills the process. So a SIGHUP the parent had set to ignored (`nohup`), or a SIGINT that a non-interactive shell ignores for a `&` job, now ends the run. There is a second effect: a caught signal is reset to the default on `exec`, so GraphAligner children no longer inherit the ignore and die from the same hangup. Version 0.1.0 had no signal code: `git grep -n signal v0.1.0 -- walker/src walker/Cargo.toml` finds no handler and no signal-hook dependency.

**Impact.** A run started with `nohup logan-walker ... &` or `nohup disassembler reconstruct ... &` on a login node is killed when the SSH session drops, because bash sends SIGHUP to every job. The walker exits with 129 and `reconstruct` prints `logan-walker (stage 1) failed (exit -1)`. Outputs are written only at the end, so the whole run is lost. Running long jobs under `nohup` is common. tmux, screen and sbatch are not affected.

**Evidence.** An interactive bash on a pseudo-terminal ran `nohup sh -c 'logan-walker ... --offline ...; echo status' &` over 4 cached accessions, with a stand-in GraphAligner that sleeps 2 s, next to a `nohup sleep` control. Closing the terminal's master side simulated the SSH drop. On 0.1.1, `/proc/<pid>/status` showed SIGHUP caught by the walker but ignored by the control. One second after the hangup the walker was dead (exit 129, no output directory) and the control was still alive. With the patch, SIGHUP stayed ignored and all 4 accessions completed with exit 0. Under `nohup disassembler reconstruct` followed by `kill -HUP -<pgid>`, 0.1.1 printed `GraphAligner failed: exited with signal: 1 (SIGHUP)`.

**Proposed fix.** Register only the signals that are not already ignored. The ignore then also passes through to GraphAligner.

```rust
let wanted: Vec<i32> =
    [SIGINT, SIGTERM, SIGHUP].into_iter().filter(|&s| !is_ignored(s)).collect();
if wanted.is_empty() { return Ok(()); }
let mut signals = signal_hook::iterator::Signals::new(&wanted)?;
...
fn is_ignored(sig: libc::c_int) -> bool {
    unsafe {
        let mut old: libc::sigaction = std::mem::zeroed();
        libc::sigaction(sig, std::ptr::null(), &mut old) == 0 && old.sa_sigaction == libc::SIG_IGN
    }
}
```

One behaviour changes: walkers started with `&` from a script keep running when the user presses Ctrl-C, as they did in 0.1.0 and as any POSIX tool does. WF-15, WF-17 and WF-20 change the same handler, so they are best done in one patch.

**Tests and docs.** Add `an_inherited_sighup_ignore_survives_a_hangup` to `walker/tests/cli.rs`. In `pre_exec`, set `signal(SIGHUP, SIG_IGN)` and `process_group(0)`; send `kill(-pgid, SIGHUP)` during alignment; assert exit 0 and no `"error"` in the summary. The current code fails this test; with the fix, all CLI tests pass. In `docs/walker.md:265-266` and the module doc at `scratch.rs:11`, add: a signal the parent set to ignored (SIGHUP under `nohup`) stays ignored.

### WF-03 · medium · Fetch then `--offline` turns `no_logan_data` / `too_large` accessions into errors (exit 1 or 3)

`walker/src/stages/stream.rs:110` · correctness · confirmed (reproduced) · fix effort: small · sources: RB-03, BL1-05, PC-06, R1-06, RC-04

**Problem.** Version 0.1.1 reports HTTP 404 (`no_logan_data`) and the size cap (`too_large`) separately from `error`, and `--fetch-only --allow-empty` accepts both. But nothing is written to the cache for these accessions. The documented follow-up `--offline` run (`docs/walker.md:161-167`) therefore sees an ordinary cache miss: `if self.offline { return Err(Error::Fetch { .. "cache miss under --offline; expected {}" ..}) }` (`walker/src/stages/stream.rs:110-117`). `unfetched` (`pipeline.rs:259-265`) treats only `NoLoganData` and `TooLarge` as negatives, so the accession becomes `Outcome::Failed`, and `check_run` (`main.rs:408-412`) fails the run whatever `--allow-empty` says.

**Impact.** Same list, same flags: the online run says `no_logan_data` and exits 0, while the offline run says `error` and exits 1, or 3 if no accession has data. Summaries mislabel missing data as failures. In a scatter across shards, every shard without a hit fails. 404s are common in SRA lists. No haplotype is wrong.

**Evidence.** Against live S3, `logan-walker --fetch-only --unitigs-cache c -a SRR0000001 --allow-empty` printed `not cached: ... (HTTP 404)`, exited 0 and left the cache empty. The analysis run with `--unitigs-cache c --allow-empty` printed `no_logan_data` and exited 0. The same run with `--offline` printed `ERROR: ... cache miss under --offline` and `all 1 accessions failed with an error`, and exited 3.

**Proposed fix.** Record negatives in the cache as a one-line sidecar file, `<DIR>/<ACC>.unitigs.absent`, holding `no_logan_data` or `too_large <B>`, where B is a lower bound on the compressed size.
1. Change `FetchAttempt::TooLarge` to `TooLarge(u64, String)`, carrying the Content-Length, or `cap+1` when the body limit found the overrun.
2. In `Fetcher::fetch`, when a cache directory is set: on `Missing` or `TooLarge`, write the sidecar (temp file plus rename); on success, remove any sidecar.
3. Before the offline miss error: a `no_logan_data` sidecar gives `Error::NoLoganData`; `too_large B` with `B > max_compressed_bytes` gives `Error::TooLarge`; otherwise (no sidecar, or the cap has since been raised) keep `Error::Fetch`.
4. Online runs ignore sidecars, because Logan may add data later, and refresh them.

In a patched build, `--fetch-only` over `SRR0000001` plus `SRR8209082 --max-size-mb 1` wrote both sidecars. The `--offline` run then reported `no_logan_data` and `too_large` and exited 0; the unpatched binary on the same cache exited 3.

**Tests and docs.** Add the CLI test `an_offline_miss_with_a_negative_sidecar_keeps_its_status`, plus a unit test showing that a garbage sidecar is treated as a plain miss. In the "Offline cache" section of `docs/walker.md`, document `.absent` files and say that deleting one forces a refetch. Update the `--allow-empty` help in `src/disassembler/cmd/fetch.py`.

### WF-04 · medium · Scratch sweep recursively deletes any user directory named `logan-walker.<integer>[.*]`

`walker/src/scratch.rs:139` · data-loss · confirmed (reproduced) · fix effort: trivial · sources: R1-02

**Problem.** The sweep accepts any name that passes `strip_prefix(PREFIX)`, `rest.split('.').next()` and `parse::<i32>()` (`walker/src/scratch.rs:139-143`). It checks only that the entry is a directory owned by the same uid, then calls `remove_dir_all` when `kill(n, 0)` reports no such process (`scratch.rs:148-152`). Nothing shows that the walker created the directory: the name's full shape is not checked and there is no marker file. Any directory of the user's named `logan-walker.<n>` or `logan-walker.<n>.<anything>` is deleted whenever `n` is not a live PID. Dates and Slurm job ids are above `pid_max`, so they never are.

**Impact.** This needs a non-default `--scratch-dir` that also holds user data, for example `--scratch-dir .` with `-o logan-walker.$(date +%Y%m%d)`. Every run then destroys the previous results and reports `removed N scratch directories ... by runs that were killed`. The default `/dev/shm` and Galaxy's per-job scratch never hold such names, and symlinks are not followed. It is an edge case, but it destroys data. A first reviewer rated it high; it is medium because the user has to make two independent choices.

**Evidence.** A run with `--scratch-dir .` deleted `logan-walker.20261110/reconstructed_haplotypes.fasta` and `logan-walker.8812345.results/sub/summary.json`, and exited 0. The controls `logan-walker.run7` and `logan-walker-20261110` survived.

**Proposed fix.** Sweep only exact walker names (the shape `scratch.rs:46` builds) and only directories carrying the walker's marker:

```rust
fn is_run_dir_name(n: &str) -> bool {
    let Some((pid, sfx)) = n.strip_prefix(PREFIX).and_then(|r| r.split_once('.')) else { return false };
    !pid.is_empty() && pid.bytes().all(|b| b.is_ascii_digit())
        && sfx.len() == 9 && sfx.bytes().all(|b| b.is_ascii_hexdigit())
}
```

Then require `.lock` (WF-01) and skip `.keep` (WF-14). If the lock design is deferred, the shape check alone fixes this defect.

**Tests and docs.** Add the unit test `directories_that_only_look_like_ours_are_never_removed`, covering `logan-walker.20261109`, `logan-walker.8812345.results`, a correctly shaped name without `.lock`, and a `.keep` directory. In `docs/walker.md:268`, state that only directories the walker created are ever removed.

### WF-05 · medium · summary.json records no input identity; cache files are used without any identity check

`walker/src/io/output.rs:245` · provenance · confirmed (reproduced) · fix effort: medium · sources: BL1-01, PK2-04, SWR2-02

**Problem.** `write_result` (`walker/src/io/output.rs:245-365`) writes the accession, status, counts, haplotypes and timing, but nothing that identifies the bytes analysed: no URL or cache path, size, hash, ETag or Last-Modified. Settings such as `cache_path` and `offline` are not recorded either, so a summary does not even say whether the input came from S3 or from the cache. On a cache hit, `if path.exists()` (`stream.rs:77`) is the only check, and the doc comment says "A cached file's contents are trusted" (`stream.rs:67-68`). `parse::build` throws away the header prefix (`rposition(|&b| b == b'_')`, `parse.rs:111`), although Logan headers name their run (`>SRR8209082_0`).

**Impact.** A misfiled cache entry, such as a copied or renamed cache or a Galaxy collection staged by element identifier, publishes another run's haplotype under the requested name, with exit 0 and no warning. Logan's `u/` objects are unversioned, so if one is regenerated, a rerun gives different results and no output shows that the input changed. No `u/` rewrite has been observed; sampled objects date from March 2024. Corruption is already caught by zstd frame checksums. A first reviewer rated this high; it is medium because both triggers are unobserved or need manual copying.

**Evidence.** Two cache files swapped between SRR1000001 and SRR1000002 produced `>SRR1000001_primary` carrying SRR1000002's allele, with exit 0. ERR2028820's real file filed as SRR8209082 gave `n_anchors: 220` and no warning. A live-S3 run and an `--offline` run of SRR8209082 produced summaries that are identical apart from timings.

**Proposed fix.**
- Add `InputProvenance {source: "s3"|"cache", location, bytes, sha256, etag, last_modified, fetched_at}` to `FetchedUnitigs`, and an `input` object (or `null`) per accession in summary.json. Attach it to every outcome, including errors.
- In `try_once`, keep the `etag` and `last-modified` headers and hash the body with SHA-256 (`sha2` is already in Cargo.lock). Write `<ACC>.unitigs.fa.zst.meta` next to the cache entry.
- On a cache hit, decode the first 4 KB and refuse the file if its first header has a prefix other than `{ACC}_`. Check the prefix only when there is one: DRR344513's headers are `>0 ka:f:9500.0`. Also refuse a file whose SHA-256 differs from its sidecar, and report the sidecar's ETag and dates.
- Do not add an MD5-against-ETag retry; the zstd checksum already covers that. Strip the quotes from the ETag and write `fetched_at` in RFC 3339.

**Tests and docs.** Add `a_cache_file_holding_another_accession_is_refused`, `a_cache_hit_reports_its_provenance_and_checks_the_sidecar`, and a `write_summary` assertion on `input`. Update `docs/walker.md:116` (summary fields), `:154` ("byte-identical"), `:276` (add a "Logan inputs are mutable" bullet), `config.rs:378`, and the `--unitigs-cache` help in `main.rs`.

### WF-06 · medium · `--fetch-only` reports a failed cache write as "cached" and exits 0

`walker/src/stages/stream.rs:128` · correctness · confirmed (reproduced) · fix effort: small · sources: R1-04, R6-02, RC-03

**Problem.** `Fetcher::fetch` turns a failed cache write into a warning and returns `Ok` (`walker/src/stages/stream.rs:128-137`: `if let Err(e) = write_cache_atomic(...) { eprintln!(... "warning: could not write cache entry: {e}") } ... return Ok(f);`). That is right for analysis runs. But `run_fetch_only` counts every `Ok` as cached (`pipeline.rs:213-219`: `d.0 += 1; println!("... cached ...")`). This breaks the documented contract "success means every accession is now in the cache" (`docs/walker.md:317`). The defect is also present in 0.1.0.

**Impact.** When the cache filesystem is full, over quota, read-only or not writable, every remaining accession is downloaded, thrown away and printed as `cached`. The summary says `0 failed` and the exit status is 0. `disassembler fetch` does not catch this. The later `--offline` run marks the missing accessions `error`, but it still exits 0 if any accession succeeds. A stderr warning does appear for each accession.

**Evidence.** `disassembler fetch -l accs.txt --unitigs-cache DIR --fetch-workers 1` was run with two accessions on a read-only bind mount. It printed `Read-only file system (os error 30)` warnings, then `[2/2] DRR000252: cached (7954205 bytes ...)` and `0 failed`, and exited 0 with the cache empty. On a 3 MB tmpfs it printed `No space left on device`, then `cached`, and exited 0.

**Proposed fix.**
1. Add a flag to `Fetcher`: `require_cache_write = cfg.fetch_only`. In fetch-only mode a write error becomes an error for that accession (a new `Error::CacheWrite { accession, path, source }` or `Error::Fetch`). Analysis runs keep the warning.
2. In `main.rs`, before `run_fetch_only`, call a `probe_cache_dir(dir)` that creates the directory, then writes and removes `.probe.<pid>.tmp`. An unwritable directory then exits 1 before any download.
3. Optional: on `StorageFull`, `QuotaExceeded` or `ReadOnlyFilesystem`, set an `AtomicBool` so workers stop taking accessions, and count the rest as errors. This stops wasted downloads; these `ErrorKind` variants need Rust 1.83.

In a patched build, the read-only case printed `error: cannot write into cache directory ...` and exited 1, and the full-disk case printed `FETCH FAILED: ... could not write the cache entry` and exited non-zero.

**Tests and docs.** Add a unit test: in fetch-only mode, storing into `<file>/cache` returns an error, while analysis mode keeps the bytes. Add a CLI test: `--fetch-only --unitigs-cache <file>/cache` exits 1 before any network access. Update `docs/walker.md:317` and the comments at `stream.rs:64-67`, `stream.rs:129-131` and `pipeline.rs:180-182`.

### WF-07 · medium · 1.5 s fixed retry budget turns brief outages into per-accession errors; 408/429 never retried

`walker/src/stages/stream.rs:161` · reliability · confirmed (reproduced) · fix effort: small · sources: R1-08

**Problem.** A transient failure is retried at most `fetch_attempts = 3` times (`config.rs:476`; no CLI flag sets it) with `std::thread::sleep(Duration::from_millis(500 * attempt as u64))` (`walker/src/stages/stream.rs:161`). That is 0.5 s plus 1 s, with no jitter, no Retry-After handling and no pause for the whole run. `stream.rs:184` (`400..=499 => return Err(FetchAttempt::Permanent(...))`) makes 408 and 429 fail on the first try.

**Impact.** While the network, proxy or DNS is down, each fetch worker discards one accession every 1.5 s. A 60 s interruption with 8 workers marks about 320 accessions `error`. In a full run the exit status is still 0 if any accession reconstructed, so the gap shows only in `summary.json` and the log. A throttling 429 is fatal for that accession immediately. The Python `harvest.py:97-122` already handles Retry-After and uses longer waits.

**Evidence.** `HTTPS_PROXY=http://127.0.0.1:9 logan-walker -l accs.txt --unitigs-cache c3 --fetch-only` with 16 accessions printed 16 × `giving up after 3 attempts: ... Connection refused` and exited 3 after 3.0 s.

**Proposed fix.**
1. Change `FetchAttempt::Transient` to `Transient(String, Option<Duration>)`. In `try_once`, map `408 | 429` to `Transient` carrying the parsed `Retry-After`, capped at 300 s.
2. Sleep `retry_after.unwrap_or_else(|| backoff(attempt))`, where `backoff` is 1, 2, 4, 8, 16 s … with ±50 % jitter, capped at 60 s.
3. Default `fetch_attempts` to 6 and add a `--fetch-attempts N` flag.
4. Optional (medium effort): a shared count of consecutive transport failures in `pipeline.rs`. Above a threshold, workers pause and probe again; if the network never comes back, the run aborts with a clear message instead of working through the rest of the list.

In a patched build, an accession gave up after 32 s instead of 1.5 s. A dead network now costs about 30 s per accession per worker, which the outage gate in step 4 bounds.

**Tests and docs.** Add `throttling_is_transient_and_honours_retry_after` (429 with `Retry-After: 7` gives a 7 s wait; 408 is transient) and `backoff_grows_and_is_capped`. Document the retry policy and `--fetch-attempts` in `docs/walker.md` and `--help`. See also WF-09 and WF-10, which change the same agent.

### WF-08 · medium · Compressed unitig files held whole in RAM in 4+F+W places; cache hits copied needlessly

`walker/src/stages/stream.rs:92` · memory · confirmed (reproduced) · fix effort: small · sources: BL1-02, RF-01

**Problem.** Each in-flight accession's compressed `.zst` file is held as a `Vec<u8>`: in the bounded queue (`channel_capacity: 4`, `config.rs:472`), in each fetch worker blocked in `send`, and in each compute worker for the whole parse (`pipeline.rs:532-536`). `Config::peak_compressed_buffers` = 4 + F + W, where F is the number of fetch workers and W the number of compute workers. At the defaults that is 20 buffers, or 20 GB of the printed 33.6 GB bound. The model is documented honestly (`config.rs:395-409`), but most of it can be avoided. A cache hit runs `let compressed = std::fs::read(&path)` (`walker/src/stages/stream.rs:92`), copying a local file into RAM, although `zstd::stream::read::Decoder` can stream from a `File` with about 128 KB of state. Local reads always outrun compute, so in the documented `--unitigs-cache --offline` workflow the worst case is the steady state.

**Impact.** Peak RSS is (4+F+W) × file size, plus W graphs. Users and the Galaxy wrapper must reserve that memory, so the wrapper lowers `--max-size-mb`, and larger accessions are skipped as `too_large`. No output is wrong.

**Evidence.** 16 hard links to a 46.4 MB synthetic file, run `--offline`: at `-w 1 --fetch-workers 1`, peak RSS was 355.6 MiB shipped against 92–97 MiB patched; at `-w 8 --fetch-workers 8`, it was 1055 MiB against 679 MiB. Outputs were byte-identical.

**Proposed fix.** Stage 1 (about 40 lines):

```rust
// types.rs
pub enum Payload { Memory(Vec<u8>), File(std::fs::File) }
// stream.rs fetch(), cache hit: open first, then check the handle's size against the cap
let file = File::open(&path)?; let len = file.metadata()?.len();
if len > cap { TooLarge } if len == 0 { Fetch error }
return Ok(FetchedUnitigs { payload: Payload::File(file), compressed_len: len as usize, .. });
// stream.rs decoder()
Payload::File(f) => Box::new(zstd::stream::read::Decoder::new(f)?),
Payload::Memory(v) => Box::new(Decoder::with_buffer(Cursor::new(v.as_slice()))?),
```

Stage 2: with `--unitigs-cache`, stream the download into the temp file with `io::copy` and return `Payload::File`. Stage 3: set `channel_capacity` to 0. Then `peak_compressed_buffers` counts only payloads held in memory.

**Tests and docs.** Update the cache tests `a_cache_hit_returns_the_bytes_on_disk` and `a_cached_file_over_the_size_cap_is_too_large` to decode and compare. Add a test that a cache hit returns `Payload::File`. Update `docs/walker.md:231-235`, the module doc at `pipeline.rs:56-80`, `config.rs:395-420`, and the memory formula in the tools-iuc wrapper. See WF-16 for the same read in `--fetch-only`.

### WF-09 · medium · Fixed 900 s whole-request fetch timeout discards partial downloads and restarts from byte 0

`walker/src/stages/stream.rs:49` · reliability · confirmed (reproduced) · fix effort: small · sources: RE-06, R1-07

**Problem.** `Fetcher::new` sets `.timeout_global(Some(Duration::from_secs(cfg.fetch_timeout_secs)))` (`walker/src/stages/stream.rs:49`). In ureq 3.4.1 that timeout covers the whole request, from DNS lookup to the end of the body. The default is 900 s (`config.rs:474`); no flag changes it, and it does not scale with `--max-size-mb`. A body timeout becomes `Transient("body read failed: …")` (`stream.rs:216`), and the retry issues a new GET with no `Range` header, discarding the bytes already received.

**Impact.** A transfer fails when its per-stream rate is below size ÷ 900 s. For a 1 GiB object, the default cap, that is 1.14 MB/s per stream, or 9.1 MB/s in total with 8 workers. The default setup on a fast link is fine. Slow or shared links, nodes far from us-east-1, `--fetch-workers 16`, or a cap above about 2 GB can hit it. A large accession then holds a worker for about 45 minutes and ends as `error`. A stalled connection is detected only after 900 s.

**Evidence.** In a test, a local server sent a 1,000,000-byte object at 250 KB/s with `fetch_timeout_secs: 1`. The fetch failed after 3 connections and 850,000 bytes served (`giving up after 3 attempts: body read failed: timeout: global`). With the resume fix and a server that honours `Range`, it succeeded after 4 connections.

**Proposed fix.**
1. Give `try_once` two extra parameters, `got: &mut Vec<u8>` and `etag: &mut Option<String>`. When `got` is not empty, send `Range: bytes={got.len()}-` and `If-Match: <etag>`. A 200 clears `got`. A 206 must have a `Content-Range` start equal to `got.len()`, or `got` is cleared. A 412 clears `got`. Read with `.take(cap + 1 - got.len()).read_to_end(got)`, which keeps the bytes received before an error.
2. In `fetch()`, keep `got` and `etag` across attempts, and count against `fetch_attempts` only those attempts that added less than 8 MiB.
3. Add `--fetch-timeout <SECS>`, now a per-attempt budget; about 300 s is reasonable once attempts resume.

Do not just raise the timeout: that only delays the detection of stalled connections.

**Tests and docs.** Turn the slow-server probe into a unit test, and add tests for a 200 reply to a Range request and for a 412. Update the doc comment at `config.rs:347`, and the options table in `docs/walker.md` (around lines 197-198).

### WF-10 · medium · Fetch trusts only compiled-in Mozilla roots; `SSL_CERT_FILE` ignored behind TLS-inspecting networks

`walker/src/stages/stream.rs:48` · usability · confirmed (reproduced) · fix effort: small · sources: R1-05

**Problem.** `Fetcher::new` builds the ureq agent with timeouts and `http_status_as_error(false)` only (`walker/src/stages/stream.rs:48-54`). As a result, ureq's default `RootCerts::WebPki`, the Mozilla roots compiled into the binary, is the only trust store. The walker honours `HTTPS_PROXY` and `NO_PROXY`, but no variable or flag adds a site CA. `stream.rs:177` maps every `call()` error to `Transient`, so a certificate error is also retried 3 times.

**Impact.** On corporate, hospital or institutional networks that re-sign TLS, curl, conda and the Python side of disassembler work, but every walker fetch fails with `invalid peer certificate: UnknownIssuer` and the run exits 3. The only workaround is to fill the cache with another tool and run `--offline`. Nothing in the docs mentions this.

**Evidence.** On a network whose direct route to S3 is re-signed by a private CA listed in `SSL_CERT_FILE`, `NO_PROXY=.amazonaws.com logan-walker -a SRR8209082 --unitigs-cache c2 --fetch-only` failed with `UnknownIssuer`, and curl in the same environment succeeded. With the fix, the same command cached a file byte-identical to the proxied download.

**Proposed fix.**

```rust
let roots = std::env::var_os("SSL_CERT_FILE")
    .and_then(|p| std::fs::read(p).ok())
    .map(|pem| parse_pem(&pem).filter_map(|i| match i {
        Ok(PemItem::Certificate(c)) => Some(c), _ => None }).collect::<Vec<_>>())
    .filter(|c| !c.is_empty())
    .map(|c| RootCerts::new_with_certs(&c))
    .unwrap_or(RootCerts::WebPki);
// Fetcher::new: .tls_config(TlsConfig::builder().root_certs(roots).build())
```

Map certificate errors to `Permanent`, with a hint ("behind a TLS-inspecting proxy? set SSL_CERT_FILE"). Log which trust store is in use, and warn rather than fall back silently when the bundle is malformed. The specific roots replace the Mozilla ones, so consider an explicit `--ca-bundle FILE` flag. An alternative is ureq's `platform-verifier` feature, which was not tried.

**Tests and docs.** Add a unit test with a local TLS server and a self-signed CA, setting the bundle through a `Config` field: the fetch fails as `Permanent` without the bundle and succeeds with it. Add a "Network" section to `docs/walker.md` covering `HTTPS_PROXY`, `NO_PROXY` and `SSL_CERT_FILE`.

### WF-11 · medium · Query FASTA read twice: piped multi-record panel refused, duplicate warnings, two parsers

`walker/src/main.rs:460` · correctness · confirmed (reproduced) · fix effort: small · sources: RS-07

**Problem.** With a single `-q`, `main` parses it with `read_query` (`let first = io::input::read_query(&args.query[0])?;`, `walker/src/main.rs:460`). If the file holds more than one record, it opens and parses the same path again with `read_query_set` (`main.rs:462`). The two functions are full copies of the FASTA parser and its per-record checks (`input.rs:124-201` and `input.rs:219-325`). There are three consequences. A multi-record panel given through a pipe or process substitution, which README.md:96-99 says works for every input, is empty on the second read and is refused. The first record's warnings print twice, in two wordings. A parser fix made in one copy only would make acceptance depend on the number of `-q` options. The doc comments at `input.rs:118` and `input.rs:206` contradict each other.

**Impact.** `logan-walker -q <(zcat panel.fa.gz) ...` stops at startup with a false `no FASTA records found`. Only direct CLI users with pipes are affected; Galaxy and the composite commands pass real files.

**Evidence.** With `panel.fasta` holding `>geneA ACGTNNACGT` and `>geneB ATGAAATAG`, `-q panel.fasta` printed four warnings, the two problems each worded twice. `-q <(cat panel.fasta)` printed `error: usage error: /dev/fd/63: no FASTA records found`.

**Proposed fix.** In `main.rs`, replace the `if args.query.len() == 1 { read_query ... }` block with `Some(io::input::read_query_set(&args.query)?)`. This is safe: a lone single-record file still gets `name = stem` and `n_records == 1`, which is the condition for align's byte-identical path (`align.rs:162`). In `input.rs`, factor out `fasta_records(path)`, `check_record(path, id, seq)` and `file_stem(path)`, and reduce `read_query` to a test-only wrapper. This does not fix the single-record pipe case, where `align.rs:169` re-reads the query path (reported in the align section).

**Tests and docs.** Add a CLI test that feeds a multi-record panel through a FIFO, reaches the GraphAligner stage and prints each warning once. Rename `multi_record_uses_first_and_warns` (`input.rs:414`). Fix the doc comments at `input.rs:118-120`, `input.rs:205-207` and `main.rs:448-449`.

### WF-12 · low · Lowercase SRA accession is fetched verbatim, 404s, and is recorded as `no_logan_data`

`walker/src/io/input.rs:107` · usability · confirmed (reproduced) · fix effort: trivial · sources: SWR1-04

**Problem.** `validate_accession` checks the SRA shape case-sensitively (`matches!(&s[..3], "SRR" | "ERR" | "DRR")`, `walker/src/io/input.rs:107-109`) and only prints a generic warning. The value goes into the S3 URL unchanged, the bucket is case-sensitive, and the resulting 404 becomes `NoLoganData` (`pipeline.rs:260`). So `-a srr8209082 --allow-empty` exits 0 with an empty FASTA and `"status": "no_logan_data"`, although Logan has the run: a silent false negative.

**Proposed fix.** Make `validate_accession` return `Result<String>`. If the uppercased value matches `(SRR|ERR|DRR)\d+`, warn once and return the uppercased form. Use the returned value at `main.rs:481` and `input.rs:63`, so deduplication, URL, cache path and output names agree.

**Tests and docs.** Add unit tests: `srr8209082` becomes `SRR8209082`, and `odd_but_safe-1` is unchanged. Add a line to the input section of `docs/walker.md`.

### WF-13 · low · Unitigs cache accepts and trusts non-zstd entries; a bad entry fails forever, even online

`walker/src/stages/stream.rs:96` · reliability · confirmed (reproduced) · fix effort: small · sources: RE-08

**Problem.** A cache hit is accepted if the file is not empty (`walker/src/stages/stream.rs:96-108`), and an online run never refetches a bad entry. `try_once` checks only status, size and non-emptiness, so any 200 body, such as a proxy error page, is cached (`stream.rs:128`). A non-zstd entry fails that accession on every run with `parse failed: ... Unknown frame descriptor`, and `--fetch-only` reports it as `cached` and exits 0. The user has to find and delete the file by hand.

**Proposed fix.** Add `looks_like_zstd(b)`, which accepts frame magic `28 b5 2f fd` or a skippable frame `0x184D2A5?` (little-endian). In `try_once`, treat a non-zstd body as `Transient` and never cache it. On a cache hit with a bad entry: under `--offline`, return an error that names the path ("delete it, or rerun without --offline"); otherwise warn, remove the file and download again.

**Tests and docs.** Add tests that a non-zstd 200 body is refused and that an offline bad entry names the path. Update the test fixture at `stream.rs:508` to start with the zstd magic. Update `docs/walker.md:170` and the comment at `stream.rs:67-69`. Related: WF-16, WF-18 and WF-19.

### WF-14 · low · `LOGAN_KEEP_SCRATCH` directory is deleted by the next run's sweep

`walker/src/scratch.rs:151` · usability · confirmed (reproduced) · fix effort: trivial · sources: RB-09, D1-08, R1-13, RX-07

**Problem.** `LOGAN_KEEP_SCRATCH` (`walker/src/scratch.rs:39`) only stops `Drop` from removing the directory (`scratch.rs:90-91`). Nothing on disk marks the directory as kept, so the next run's sweep (`scratch.rs:151`) removes it as "left ... by runs that were killed". This happens even when that run also sets the variable. Debugging accession A and then accession B in two runs loses A's GFA and GAF, which breaks the workflow in `docs/issue_graphaligner_end_clipping.md:99`. This is a 0.1.1 regression.

**Proposed fix.**

```rust
const KEEP_MARKER: &str = ".keep";
// Scratch::create, after make_private_dir succeeds:
if keep { std::fs::write(path.join(KEEP_MARKER), b"")
              .map_err(file_err("mark the kept scratch directory", &path))?; }
// sweep_stale, after the uid/own-pid check:
if entry.path().join(KEEP_MARKER).exists() { continue; }
// Drop: "LOGAN_KEEP_SCRATCH: kept {} (never swept; remove it yourself)"
```

**Tests and docs.** Add `a_kept_directory_survives_the_next_runs_sweep`. In `docs/walker.md:268` and `scratch.rs:13`, say that a kept directory is never swept and must be removed by hand. WF-01's lock design carries the same marker.

### WF-15 · low · Signal handler removes scratch and exits without stopping running GraphAligner children

`walker/src/scratch.rs:80` · reliability · confirmed (reproduced) · fix effort: small · sources: X1-04

**Problem.** On SIGINT, SIGTERM or SIGHUP, the signal thread only removes the scratch directory and exits (`walker/src/scratch.rs:76-81`). Children started with `cmd.output()` (`align.rs:201`) are not tracked. After `kill <walker pid>` or `timeout N logan-walker`, up to `-w` GraphAligner processes keep using CPU and RAM, reparented to init, on files that no longer exist. Ctrl-C in a terminal, and runners that kill the whole process group or cgroup, are not affected.

**Proposed fix.** Spawn each aligner in its own process group (`CommandExt::process_group(0)`). Record its PID in a static `Mutex<Vec<u32>>` and remove it right after `wait`. In the signal thread, send `kill(-pgid, SIGTERM)` to every recorded group before removing scratch. On Linux, `prctl(PR_SET_PDEATHSIG, SIGKILL)` in `pre_exec` also covers SIGKILL of the walker. Do not use `kill(0, sig)`: that signals the caller's whole process group.

**Tests and docs.** Extend `sigterm_removes_the_scratch_directory` with an aligner that writes its PID and runs `exec sleep 30`, and assert the child is gone within 1 s. Update `docs/walker.md:266`.

### WF-16 · low · `--fetch-only` reads every cached file whole, validates nothing, counts hits as fetched

`walker/src/stages/stream.rs:92` · efficiency · confirmed (reproduced) · fix effort: small · sources: RC-09, R1-12, RE-10, RF-07

**Problem.** `run_fetch_only` calls `fetcher.fetch` (`pipeline.rs:213`) only to print the size. On a hit, `std::fs::read(&path)` (`walker/src/stages/stream.rs:92`) loads the whole file, with up to `--fetch-workers` such files in memory at once. Topping up a 100 GB cache therefore re-reads all of it. Zero-filled files are accepted as "cached", hits are counted as "fetched", and a cached file over the current `--max-size-mb` fails the run with a message about `max_compressed_bytes`. This contradicts `docs/walker.md:170` ("`fetch` skips accessions already in the cache"). On 4 × 150 MB cached files, peak RSS was 537 MB; with the fix it was 12 MB.

**Proposed fix.** Add `Fetcher::ensure_cached(index, acc) -> Result<(u64, bool)>`. It stats the file, rejects an empty file, checks the 4-byte zstd magic (or treats a bad file as a miss when online, as in WF-13), and returns `(len, true)` without reading the rest. Otherwise it calls `fetch`. Print "already cached" for hits, and report "N fetched, M already cached" (`main.rs:518`). Name `--max-size-mb` in the over-cap message.

**Tests and docs.** Add unit tests for a valid file, a zero-filled file and an empty file. Add a CLI test: `--fetch-only --offline` prints "already cached". Related: WF-08.

### WF-17 · low · Scratch directory left behind when the run dies by SIGPIPE, SIGQUIT, SIGUSR1/2 or SIGXCPU

`walker/src/scratch.rs:73` · reliability · confirmed (reproduced) · fix effort: small · sources: R1-10, RB-08

**Problem.** The cleanup thread watches only `[SIGINT, SIGTERM, SIGHUP]` (`walker/src/scratch.rs:73`), and `main.rs:331` resets SIGPIPE to the default. So `logan-walker ... | head -1`, Ctrl-\, or a scheduler's soft limit (SIGXCPU, SIGUSR1/2) kills the process before cleanup runs. The directory stays in `/dev/shm`, using RAM, until the same user starts a run on that host. This contradicts `--help` and `docs/walker.md:266`. Measured: exit statuses 141, 138, 140, 152 and 131 each left the directory; SIGTERM did not.

**Proposed fix.**
1. Also register SIGPIPE, SIGQUIT, SIGUSR1, SIGUSR2, SIGXCPU, SIGXFSZ, SIGALRM, SIGVTALRM and SIGPROF, skipping any signal inherited as ignored (WF-02).
2. In `pipeline::report`, use `let _ = writeln!(stdout(), ...)` instead of `println!`.
3. Remove the production `create_dir_all` at `extract.rs:293`, which an in-flight worker otherwise uses to recreate the directory after the handler removed it.

**Tests and docs.** Extend `a_closed_stdout_ends_the_run_instead_of_hanging` to assert that the base is empty, and add a SIGUSR1 test (status 138, empty base). Update the signal list in the docs.

### WF-18 · low · Runtime failures labelled "usage error"; non-zstd payloads "parse failed"; non-UTF-8 GAF fails every accession

`walker/src/error.rs:75` · usability · confirmed (reproduced) · fix effort: small · sources: SWR1-05

**Problem.** `Error::Usage(m) => write!(f, "usage error: {m}")` (`walker/src/error.rs:75`) is also used for GAF content (`align.rs:269`, `:301`) and for scratch creation (`scratch.rs:51`, `:59`). `stream::decoder` claims to validate the frame header (`stream.rs:289`), but `Decoder::with_buffer` reads nothing, so a non-zstd payload surfaces as `parse failed: read failed after 0 records`, and `Error::Decompress` is never produced. Worse, GraphAligner copies the whole query header into the GAF, so a Latin-1 byte in the header's description (`>gene M\xfcller`) makes `from_utf8` at `align.rs:269` fail for every accession, and the run exits 3.

**Proposed fix.** Parse the GAF with `String::from_utf8_lossy`, and report a bad numeric column as `Error::GraphAligner`. Report scratch failures as `Error::File`. In `decoder`, reject a payload without the zstd magic as `Error::Decompress`, and map read errors in `parse.rs:94` to `Decompress`. Fix the doc comment at `stream.rs:289`.

**Tests and docs.** Add unit tests: `parse_gaf_multi` on a Latin-1 header, `decoder` on `b"<html>"`, and a truncated frame in `parse::build`.

### WF-19 · low · Killed cache writer leaves orphan hidden `.tmp` files; entry not fsynced before rename

`walker/src/stages/stream.rs:262` · reliability · confirmed (reproduced) · fix effort: small · sources: R1-11

**Problem.** `write_cache_atomic` writes `.{ACC}.{pid}.{ns}.tmp` with `std::fs::write` (`walker/src/stages/stream.rs:262`) and renames it (`:264`). A process killed between those two steps (SIGKILL, the OOM killer, or SIGTERM under `--fetch-only`, which installs no handler) leaves a hidden file as large as the entry, up to 1 GB by default. Nothing ever removes these files. There is also no `sync_all` before the rename, so after a power or node crash the entry can be empty or short on some filesystems. Reproduced by delaying the rename with `strace -e inject=rename:delay_enter=8000000` and then sending SIGKILL: the `.tmp` file survived a later successful run.

**Proposed fix.** Use `File::create`, `write_all` and `sync_all()` before the rename; optionally also sync the directory after it. In `Fetcher::new`, remove files whose names exactly match `.<ACC>.<digits>.<digits>.tmp` and whose mtime is more than 24 h old. Use age rather than PID, so a cache shared over NFS stays safe.

**Tests and docs.** Add a sweep test that sets mtimes with `File::set_modified`. Correct the promise at `docs/walker.md:151-152`.

### WF-20 · low · As PID 1 of a container, SIGTERM/SIGINT ends the walker with SIGSEGV (exit 139)

`walker/src/scratch.rs:80` · reliability · confirmed (reproduced) · fix effort: trivial · sources: SWR1-03

**Problem.** After removing scratch, the handler calls `emulate_default_handler(sig)` (`walker/src/scratch.rs:80`). In signal-hook 0.4.5, this raises the signal and calls `abort()` if the raise returns. For PID 1 of a PID namespace, the kernel drops self-sent signals whose disposition is the default, so the raise is ignored, abort's SIGABRT is ignored too, and glibc ends with a SIGSEGV. The `exit(128 + sig)` on the next line is never reached. `docker stop` on a container without `--init` therefore reports a segfault (exit 139), and writes a core file where core dumps are enabled. Reproduced under `unshare -pf --mount-proc`.

**Proposed fix.**

```rust
unsafe { libc::signal(sig, libc::SIG_DFL); libc::raise(sig); }
std::process::exit(128 + sig);
```

**Tests and docs.** Add a CLI test, run as root under `unshare` (skipped otherwise), that expects exit 143. Optionally, document `docker run --init`.

### WF-21 · low · Unitig fetcher follows redirects to any origin, including plain http

`walker/src/stages/stream.rs:48` · security · confirmed (reproduced) · fix effort: trivial · sources: RX-05

**Problem.** The agent keeps ureq 3.4.1's defaults, `max_redirects: 10` and `https_only: false` (`walker/src/stages/stream.rs:48-54`). A 3xx on the fixed `https://logan-pub.s3.amazonaws.com/...` URL therefore moves the download to any host, over cleartext http, and the body is cached under the accession's name. A local 302 test confirmed that the redirect is followed. S3 does not redirect this bucket today, so this is defence in depth.

**Proposed fix.** Add `.max_redirects(0)` and an arm in `try_once` before `_`:

```rust
300..=399 => {
    let loc = resp.headers().get("location").and_then(|v| v.to_str().ok()).unwrap_or("-");
    return Err(FetchAttempt::Permanent(format!("unexpected redirect, HTTP {status} to {loc}; not followed")));
}
```

Without the arm, a 3xx is retried as `Transient`.

**Tests and docs.** Add `redirect_is_not_followed`. In the exit-status table of `docs/walker.md`, note that a redirect is reported as `error`.

### WF-22 · low · Inputs read with no size limit (gzip inflated whole); accession length unbounded

`walker/src/io/input.rs:43` · robustness · confirmed (reproduced) · fix effort: small · sources: RX-03

**Problem.** `read_input` reads plain inputs whole (`input.rs:28`) and inflates gzip inputs without a cap (`walker/src/io/input.rs:43`). `validate_accession` never checks length, and its warning echoes the whole value. A 32 KB gzip list holding one 32 MiB "accession" produced about 300 MB of output across the log, stdout, the summary JSON and `bubbles.json`. Only the submitting user is affected.

**Proposed fix.**

```rust
pub const MAX_INPUT_BYTES: u64 = 256 << 20;
pub const MAX_ACCESSION_LEN: usize = 64;
// read_input: refuse a file over the limit, then
read_to_end(&mut Read::take(MultiGzDecoder::new(&bytes[..]), limit + 1), &mut out)?;
if out.len() as u64 > limit { return Err(too_big()) }
// validate_accession: refuse s.len() > 64, echoing only the first 20 characters
```

**Tests and docs.** Add `an_overlong_accession_is_refused_without_echoing_it` and `a_gzip_bomb_is_refused_at_the_limit`. Document both limits in `docs/walker.md`.
