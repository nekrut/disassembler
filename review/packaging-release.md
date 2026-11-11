# Packaging, release engineering, Bioconda, containers, and claims of the 0.2 backlog and plan section 1

The repository's own pins are sound and documented with evidence: `environment.yml` pins MAFFT 7.526 and GraphAligner 1.0.17b, and `pyproject.toml` lists only the two Python packages the code imports. The weak point is the hand-off from repository to channel. The 0.1.1 Bioconda recipes were published as version-only bumps by Bioconda's autobump bot, so the walker floor, the MAFFT floor and the dependency trim that the repository assumes never shipped. As a result, the only 0.1.1 biocontainer runs logan-walker 0.1.0 under a 0.1.1 label, and nothing in the recipe test or the Python wrapper notices. `docs/packaging.md` and `environment.yml` describe the intended state rather than the published one. Several statements in the 0.2 backlog and the Galaxy plan rest on misread evidence: submodule ancestry, the Logan object format, and single-pass memory. Most fixes share one theme: make the recipe follow the repository through a checked procedure (autobump off, floors and dependencies compared with `pyproject.toml`, inspection after publishing), and have the wrapper enforce its walker floor at run time.

| ID | Severity | Where | Title | Fix effort |
|---|---|---|---|---|
| PK-01 | high | `bioconda:recipes/disassembler/meta.yaml:34` | Biocontainer `disassembler:0.1.1--pyhdfd78af_0` ships logan-walker 0.1.0; nothing checks the walker version | small |
| PK-02 | medium | `docs/packaging.md:40` | Docs hide that the published 0.1.1 can resolve MAFFT 7.525; `environment.yml` gives a bioconda-first install command | trivial |
| PK-03 | medium | `walker/src/stages/stream.rs:48` | logan-walker trusts only compiled-in roots, so every fetch fails behind TLS-inspecting proxies | small |
| PK-04 | medium | `docs/packaging.md:63` | Release procedure lets Bioconda's autobump bot publish version-only recipes | small |
| PK-05 | medium | `docs/packaging.md:27` | packaging.md blames the macOS gap on a fixed CI bug; GraphAligner 1.0.17b makes macOS uninstallable | trivial |
| PK-06 | low | `bioconda:recipes/disassembler/meta.yaml:30` | Published 0.1.1 recipe still pulls zstandard, numpy and pandas | trivial |
| PK-07 | low | `docs/packaging.md:58` | Versions are bumped only at release, so source builds report the previous release | small |
| PK-08 | low | `extern/kmer-diversity-sieve` (gitlink) | Sieve submodule pin was reverted by a docs commit; backlog (e) would move the recipe backwards | trivial |
| PK-09 | low | `docs/packaging.md:39` | Docs attribute the MAFFT 7.525 bug to `--add`; `msa` uses `--addfragments`, which 7.526 fixed | trivial |
| PK-10 | low | `environment.yml:31` | `environment.yml` still installs zstandard, which nothing imports | trivial |
| PK-11 | low | `backlog:61` | Backlog item 0's "stock install has a broken `msa`" holds only with bioconda ahead of conda-forge | trivial |
| PK-12 | low | `walker/src/stages/stream.rs:102` | Backlog 1(a) understates provenance work: cached and offline runs have no headers, manifests are not free | small |
| PK-13 | low | `plan:1(b)` / `walker/src/stages/extract.rs:93` | Plan 1(b): O(subgraph) memory needs two decodes and a spooled adjacency, not one pass | large |
| PK-14 | low | `backlog:1(b)` (line 148) | Backlog calls Logan objects unindexed zstd; they are seekable with a seek table and record index | trivial |
| PK-15 | info | `walker/src/types.rs:97` | Backlog's ~34 B/unitig CSR figure is right for one accession but varies with graph shape | trivial |

### PK-01 · high · Biocontainer `disassembler:0.1.1--pyhdfd78af_0` ships logan-walker 0.1.0; nothing checks the walker version

`bioconda:recipes/disassembler/meta.yaml:34` · packaging · confirmed (reproduced) · fix effort: small · sources: PK2-01, PK1-01

**Problem.** The only published 0.1.1 image, `quay.io/biocontainers/disassembler:0.1.1--pyhdfd78af_0`, contains disassembler 0.1.1 together with `conda-meta/logan-walker-0.1.0-hfa8f182_0.json`. Its `logan-walker` binary reports `logan-walker 0.1.0` and its help shows `--primary-metric ... [default: cov]`. The image was built on 2026-11-02, and logan-walker 0.1.1 reached the channel on 2026-11-05. The recipe floor `- logan-walker >=0.1.0,<0.2` (`meta.yaml:34`) therefore let the solver take 0.1.0. The recipe test `disassembler walker --version | grep -q "logan-walker"` (`meta.yaml:47`) accepts any version. The wrapper only checks that an executable exists: see `resolve_walker` (`src/disassembler/resolve.py:45`), called from `src/disassembler/cmd/walker.py:48` and from the `reconstruct` preflight (`src/disassembler/cmd/pipe_reconstruct.py:320`).

**Impact.** `reconstruct`, `from-kmindex`, `fetch` and `walker` run from this image (docker, apptainer, Nextflow, Snakemake, HPC modules) behave like 0.1.0. The multi-anchor primary is chosen by `cov` rather than `matched`. On a 20,000-accession SARS-CoV-2 run, `cov` tagged 62% of primaries as XBB.1.5, one of the shortest references, where `matched` gives 43% Delta, the cohort's real mix (`docs/walker.md:76-79`). Version 0.1.0 has no prune-span guard and pruning is on by default, so default single-anchor primaries can be truncated without warning (`docs/walker.md:97-99`: 425 nt lost). The 0.1.0 header, `--allow-empty` and SIGPIPE bugs also return. The 0.1.0 summary JSON has no `logan_walker` or `settings` block, so users cannot detect this from their outputs. A fresh `conda install` now gets 0.1.1, and the Galaxy tools in tools-iuc #8800 use the logan-walker image, so neither is affected.

**Evidence.** The image layer was streamed and read without executing anything; it holds the 0.1.0 conda-meta record, the `[default: cov]` help text and `__version__ = "0.1.1"`. The same pairing was rebuilt locally (disassembler 0.1.1 with logan-walker built from tag `logan-walker-v0.1.0`) and run on a two-reference panel: `disassembler walker -q panel.fa -l accs.txt --unitigs-cache cache --offline -o out` exits 0 with `anchor=refShort cov=1.005 pident=99.5` under 0.1.0, and with `anchor=refLong pident=100.0` under 0.1.1. The recipe's test command passes with 0.1.0.

**Proposed fix.**
1. Recipe: build number 1 with `- logan-walker >=0.1.1,<0.2`. Walker 0.1.1 is already on the channel, so the floor is enough. After publishing, read the new image's `conda-meta`, move users' pins to `0.1.1--pyhdfd78af_1`, and mark `_0` as do-not-use. A version regex in the recipe test adds little, because the solver already enforces the floor at build time. The run-time check below is the real guard.
2. Repository: enforce a floor at run time.
```python
# src/disassembler/resolve.py
MIN_WALKER_VERSION = (0, 1, 1)
_WALKER_VER = re.compile(r"^logan-walker (\d+)\.(\d+)\.(\d+)")

def walker_version_error(path: str) -> str | None:
    """None if `path` is a new enough logan-walker, else a message for the user."""
    try:
        out = subprocess.run([path, "--version"], capture_output=True,
                             text=True, timeout=30).stdout.strip()
    except (OSError, subprocess.SubprocessError) as e:
        return f"{path}: cannot run --version ({e})"
    m = _WALKER_VER.match(out)
    if not m:
        return f"{path}: unrecognised --version output {out!r}"
    if tuple(map(int, m.groups())) < MIN_WALKER_VERSION:
        need = ".".join(map(str, MIN_WALKER_VERSION))
        return (f"{path} is logan-walker {'.'.join(m.groups())}; "
                f"this disassembler needs logan-walker >= {need}")
    return None
```
Call it in `_resolved_or_report` (covers the `walker` passthrough, `fetch` and the composites) and in the `pipe_reconstruct.py:320` preflight. In a trial of this change, walker 0.1.0 was refused (`walker` exits 127, `reconstruct` exits 1 with the message above), and walker 0.1.1 ran unchanged.

**Tests and docs.** Add `tests/test_walker_version.py`: an old walker is refused by `reconstruct` and by the passthrough; the floor and newer versions pass; the `walker/Cargo.toml` version meets the floor. The fake walker in `tests/test_composites.py` must answer `--version`, or 10 tests fail. Change `docs/packaging.md:54` to `>=0.1.1,<0.2`. In "Making a release", raise `MIN_WALKER_VERSION` and the recipe floor together and inspect the image's `conda-meta` after publishing (see PK-04). Add a "known bad builds" note for `_0`.

### PK-02 · medium · Docs hide that the published 0.1.1 can resolve MAFFT 7.525; `environment.yml` gives a bioconda-first install command

`docs/packaging.md:40` · docs/packaging · confirmed (by reading) · fix effort: trivial · sources: D2-02, PK1-02, PK2-07

**Problem.** `docs/packaging.md:39-41` says "`environment.yml` and the Galaxy `msa` tool pin 7.526; the `0.1.0` Bioconda recipe does not set a floor yet", which a reader takes to mean that 0.1.1 does set one. It does not: the published disassembler 0.1.1 (`pyhdfd78af_0`, 2026-11-05) depends on a bare `mafft` (`meta.yaml:36`). `environment.yml:3` tells readers to run `conda install -c bioconda disassembler`, whereas `README.md:33` and `docs/packaging.md:20` put conda-forge first. Channels given with `-c` are prepended to the configured ones, so that command ranks bioconda first.

**Impact.** Bioconda's newest mafft is 7.525; conda-forge has 7.526. With bioconda first and flexible channel priority, channel rank beats version, so the solver installs 7.525, which by the repository's own account breaks `disassembler msa` and the end of `reconstruct`. With strict priority the same command fails to solve: bioconda's biopython stops at 1.70, built for Python 3.6 and older. The documented conda-forge-first command gets 7.526 (see PK-11).

**Evidence.** `curl -s https://api.anaconda.org/package/bioconda/disassembler/files` lists the 0.1.1 dependencies as `[..., 'mafft', 'numpy', 'pandas', ..., 'zstandard']`. `.../package/bioconda/mafft` reports 7.525 as latest; `.../package/conda-forge/mafft` has 7.526. No solver was run.

**Proposed fix.**
```diff
--- a/docs/packaging.md
-  `environment.yml` and the Galaxy `msa` tool pin 7.526; the `0.1.0` Bioconda recipe does not
-  set a floor yet.
+  `environment.yml` and the Galaxy `msa` tool pin 7.526. Neither the 0.1.0 nor the published
+  0.1.1 (build 0) Bioconda recipe sets a floor; it arrives with the next recipe build
+  (`mafft >=7.526`). Until then list conda-forge before bioconda, or run
+  `conda install 'mafft>=7.526'` in the environment.
--- a/environment.yml
-# For a packaged install instead, see docs/packaging.md (conda install -c bioconda disassembler).
+# For a packaged install instead, see docs/packaging.md
+# (conda install -c conda-forge -c bioconda disassembler).
```
In the next recipe build (with PK-01 and PK-06), set `mafft >=7.526` and add a recipe test that runs `disassembler msa` on a three-record fixture, so that a bad MAFFT solve fails CI. Fix the `--add` wording on line 39 at the same time (PK-09).

**Tests and docs.** Optional CI grep: every `conda install` line in the repository lists `-c conda-forge` before `-c bioconda`.

### PK-03 · medium · logan-walker trusts only compiled-in roots, so every fetch fails behind TLS-inspecting proxies

`walker/src/stages/stream.rs:48` · packaging · confirmed (reproduced) · fix effort: small · sources: R6-05

**Problem.** `walker/Cargo.toml:23` uses `ureq = "3.4.1"` with its default features (rustls with the bundled webpki root list). `Fetcher::new` (`stream.rs:47-53`) builds its agent without a TLS config, so the only trusted certificate authorities (CAs) are the Mozilla roots compiled into the binary. The walker ignores the OS trust store, `SSL_CERT_FILE` and `SSL_CERT_DIR`, and has no option for adding a CA. Certificate errors are wrapped as `FetchAttempt::Transient` (`stream.rs:177`), so they are retried for nothing.

**Impact.** Some networks re-sign outbound TLS with a site CA (Zscaler and Palo Alto-style inspection, common on institutional and HPC networks). There, curl, Python and conda work, but every accession fails with `invalid peer certificate: UnknownIssuer` after 3 attempts and the run exits 3. This affects `walker`, `fetch`, `reconstruct` and `from-kmindex`. The failure is loud. The only workaround, which no doc mentions, is to download `u/<ACC>/<ACC>.unitigs.fa.zst` with curl into `--unitigs-cache` and run with `--offline`.

**Evidence.** A local intercepting CONNECT proxy was set in `HTTPS_PROXY`, with its CA in `SSL_CERT_FILE`. Through it, curl and Python `urllib` completed TLS and received the proxy's 404. `logan-walker --fetch-only --unitigs-cache cache -a SRR8209082` failed with UnknownIssuer three times and exited 3.

**Proposed fix.** Keep the webpki roots as the default, because minimal containers may have no system store. Let `SSL_CERT_FILE` replace them, as curl and OpenSSL do. No new crates are needed.
```rust
fn tls_config_from_env() -> ureq::tls::TlsConfig {
    let b = ureq::tls::TlsConfig::builder();
    let Some(p) = std::env::var_os("SSL_CERT_FILE") else { return b.build() };
    // read p; collect PemItem::Certificate(c) => c.to_owned() from ureq::tls::parse_pem;
    // on read error or no certificates: warn naming the file, return b.build()
    b.root_certs(ureq::tls::RootCerts::new_with_certs(&certs)).build()
}
// Fetcher::new: Agent::config_builder() ... .tls_config(tls_config_from_env()) ...
```
In a trial build, a fetch through the same proxy completed TLS and reported "not available" (exit 1). With `SSL_CERT_FILE` unset, behaviour was unchanged, and fetches from real S3 still worked. Optionally, add `--ca-bundle FILE` (an unreadable file is a usage error) and classify certificate errors as permanent. Switching to the OS verifier is the weaker choice: it adds a crate and can break busybox-based containers.

**Tests and docs.** Factor out `root_certs_from_pem(&[u8]) -> Option<RootCerts>` and unit-test it with valid, empty and garbage PEM input. In `docs/walker.md`, document `SSL_CERT_FILE` and the curl-into-cache workaround in the offline-cache section and the options table. Add an admin note to `docs/galaxy.md`.

### PK-04 · medium · Release procedure lets Bioconda's autobump bot publish version-only recipes

`docs/packaging.md:63` · packaging · confirmed (reproduced) · fix effort: small · sources: PK1-03, PK2-09

**Problem.** "Making a release" tags first ("3. Tag the release commit ...", `docs/packaging.md:63`) and updates the recipes afterwards. Step 4 (`:64-65`) lists only "version, `sha256` of the new tag tarball, and build number reset to 0". Bioconda's autobump bot reacts to the tag by opening a bump of the current recipe that changes only the version and sha256. If that PR merges first, every recipe-side change the release needs is lost. Step 4 also has these gaps:
- It never mentions the `logan-walker` floor, the run-dependency list, the MAFFT floor or the sieve commit.
- "Reset to 0" is wrong when the version is already published.
- There is no check after publishing.

Step 5 refers to a Galaxy suite that is not yet merged. Backlog item 0.5 calls opening the maintainer PR "before or immediately with" the tag push "the only defence". It misses Bioconda's per-recipe opt-out, `extra: autobump: enable: false`, which bioconda-utils always honours (`bioconda_utils/autobump.py:399-402` `ExcludeDisabled`, installed unconditionally at `cli.py:1661`). None of the three recipes sets it.

**Impact.** This is how 0.1.1 shipped: the published 0.1.0 and 0.1.1 packages have identical dependencies, which caused PK-01, PK-02 and PK-06. The 0.2.0 release repeats this unless the procedure changes.

**Evidence.** `curl -s https://api.anaconda.org/package/bioconda/disassembler/files` shows the same dependency list for 0.1.0 (2026-10-26) and 0.1.1 (2026-11-05): `logan-walker >=0.1.0,<0.2`, a bare `mafft`, `numpy`, `pandas`, `zstandard`. The Bioconda master recipe differs from 0.1.0 only in version and sha256.

**Proposed fix.** Before any 0.2 tag, merge the opt-out into the disassembler and logan-walker recipes (and the sieve recipe while it tracks an untagged commit):
```yaml
extra:
  recipe-maintainers:
    - as042
  autobump:
    enable: false
```
Rewrite "Making a release":
1. Bump the versions and run the checks, including a golden run.
2. Tag.
3. Open the recipe PR at once and check every field:
   - the version;
   - the sha256 of the GitHub tag tarball (`curl -sL https://github.com/nekrut/disassembler/archive/refs/tags/<tag>.tar.gz | sha256sum`, not `git archive`);
   - the build number: 0 for a new version, +1 when that version is already published;
   - `logan-walker >=<new walker>,<next minor>`;
   - run deps = `pyproject.toml` `dependencies` plus the binaries;
   - `mafft >=7.526`;
   - the sieve commit = the submodule pin.
4. After the merge, check that `conda create --dry-run -c conda-forge -c bioconda disassembler=<v>` shows the expected walker and mafft, and that the biocontainer's `conda-meta` contains them.
5. Update the Galaxy version tokens only once tools-iuc #8800 is merged.

Correct backlog 0.5 to name the opt-out.

**Tests and docs.** Optionally add a CI job that compares the recipe's floors and run deps with the repository (see PK-06). With autobump off, a forgotten manual PR means no release reaches Bioconda; step 4 catches that.

### PK-05 · medium · packaging.md blames the macOS gap on a fixed CI bug; GraphAligner 1.0.17b makes macOS uninstallable

`docs/packaging.md:27` · packaging · confirmed (by reading) · fix effort: trivial · sources: PK2-03, D2-01

**Problem.** `docs/packaging.md:27-30` says the macOS builds "were switched off in `0.1.0` because of a Bioconda CI problem (bioconda-recipes#70401, since fixed); re-enabling them is planned for the next release". The pinned GraphAligner 1.0.17b has only two osx-64 builds, `h66d0458_0` and `h8d12793_2`. Both depend on `protobuf 3.14.0`, which conda-forge provides for osx-64 only for Python 3.6 to 3.9. disassembler needs Python >= 3.11. No osx-arm64 GraphAligner exists. Linux solves only because of the protobuf-free build `hdcf5f25_1` (`docs/issue_graphaligner_end_clipping.md:222`). Version 0.1.1 also shipped Linux-only, so the "next release" promise is already stale.

**Impact.** The bioconda-recipes #71068 branch drops `skip: True  # [osx]` from logan-walker. If it merges, an osx-64 logan-walker is published and passes CI, because on its own it pulls Python 3.9. `conda install -c conda-forge -c bioconda disassembler` on osx-64 is then unsatisfiable, and nothing catches it because the noarch disassembler package is tested only on Linux. `conda env create -f environment.yml` fails the same way on a Mac. Users get a solver error, not wrong results, but the docs promise a platform the dependency graph cannot deliver.

**Evidence.** From `https://api.anaconda.org/package/bioconda/graphaligner`, version 1.0.17b has linux-64 builds `hd03093a_0` and `hdcf5f25_1` (no dependencies) and `h21ec9f0_2` (protobuf 3.14.0), and its osx-64 builds both need protobuf 3.14.0. conda-forge's protobuf 3.14.0 for osx-64 exists for py36 to py39 only. No solve was run.

**Proposed fix.**
```diff
--- a/docs/packaging.md
-... The macOS builds were switched off
-in `0.1.0` because of a Bioconda CI problem (...#70401, since fixed); re-enabling them is
-planned for the next release.
+... `logan-walker`, and therefore `disassembler`, are linux-64 only: GraphAligner 1.0.17b has
+no osx-arm64 build, and its osx-64 builds depend on `protobuf 3.14.0`, which exists only for
+Python <= 3.9, while `disassembler` needs Python >= 3.11. macOS returns once a protobuf-free
+osx-64 GraphAligner 1.0.17b build exists (as `hdcf5f25_1` is for Linux).
```
In #71068, keep the osx skip on logan-walker and give its comment the real reason. Re-enable osx only for kmer-diversity-sieve. An alternative is to bump the graphaligner 1.0.17b recipe's build number first, so that it gets a protobuf-free osx-64 build. Add "linux-64 only (see docs/packaging.md#platforms)" to the `environment.yml` header.

**Tests and docs.** Add a release step: `conda create --dry-run --platform osx-64 -c conda-forge -c bioconda disassembler` must solve before a platform is enabled. Update backlog item 0.2, which blames the gap on the skip alone.

### PK-06 · low · Published 0.1.1 recipe still pulls zstandard, numpy and pandas

`bioconda:recipes/disassembler/meta.yaml:30` · packaging · confirmed (reproduced) · fix effort: trivial · sources: PK1-05

**Problem.** `meta.yaml:30-32` lists `zstandard`, `numpy` and `pandas` as run requirements. `pyproject.toml:12-18` dropped them in 0.1.1 because nothing in `src/` imports them. Every conda install and the 0.1.1 image download about 27 MB more (numpy 9.2 MB, pandas 16.9 MB, zstandard 0.5 MB), and the solver gets more constraints next to GraphAligner's `libzlib <1.3` cap. With the three modules blocked through `sys.modules[name] = None`, all 17 `disassembler` modules still import.

**Proposed fix.** In the build-1 recipe (together with PK-01 and PK-02):
```diff
-    - zstandard
-    - numpy
-    - pandas
-    - logan-walker >=0.1.0,<0.2
+    - logan-walker >=0.1.1,<0.2
     - kmer-diversity-sieve >=0.1.0,<0.2
-    - mafft
+    - mafft >=7.526
```

**Tests and docs.** Add a release check (about 20 lines with `tomllib`) that the recipe's Python run deps equal `[project].dependencies`. Change "plus Python deps" in `docs/packaging.md:8` to "plus biopython and requests".

### PK-07 · low · Versions are bumped only at release, so source builds report the previous release

`docs/packaging.md:58` · provenance · confirmed (by reading) · fix effort: small · sources: SWK-02

**Problem.** The checklist bumps `__version__` and `walker/Cargo.toml` only in the release commit. The walker's version is `env!("CARGO_PKG_VERSION")` (`walker/src/main.rs:345`, `walker/src/io/output.rs:201`), with no commit id. Builds of `main` between releases therefore claim to be the previous release. For example, at 7e5653d `Cargo.toml` says 0.1.0 but `--primary-metric` defaults to `matched`, while v0.1.0 defaulted to `cov`: two binaries both reporting "logan-walker 0.1.0" pick different primaries. This affects users who build from source (`README.md:46-49`), not Bioconda users.

**Proposed fix.** Add a post-tag step that bumps to `0.1.2-dev` in Cargo (regenerating `Cargo.lock`) and to `__version__ = "0.1.2.dev0"`. Optionally, add a `walker/build.rs` that emits `cargo:rustc-env=LOGAN_WALKER_GIT=<git describe --tags --always --dirty>`, or an empty string when there is no `.git`. Print it in `--version` and write it as `logan_walker_commit` in `summary.json`.

**Tests and docs.** Add a unit test that `--version` starts with `CARGO_PKG_VERSION`. Document the new field in `docs/walker.md`.

### PK-08 · low · Sieve submodule pin was reverted by a docs commit; backlog (e) would move the recipe backwards

`extern/kmer-diversity-sieve` (gitlink) · packaging · confirmed (reproduced) · fix effort: trivial · sources: K2, PK2-08

**Problem.** History of the submodule pin:
- f0d9285 pinned `39d90c0`.
- 00e3b97 advanced it to `70022da`, the commit that the Bioconda sieve recipe builds and that `docs/packaging.md:42` names.
- b1b6acf, "Docs changes", set it back to `39d90c0`.

Upstream `main` is `70022da`, and `39d90c0` is its parent; the two differ only by 56 README lines. Backlog (e), lines 203-204, misreads this as upstream rewriting history and proposes moving the recipe to `39d90c0`. That would move the package backwards. There is no functional effect today.

**Evidence.** `git ls-tree HEAD extern/kmer-diversity-sieve` prints `39d90c0`; `git log --format=%h -- extern/kmer-diversity-sieve` shows the three pins; `git ls-remote https://github.com/stevenweaver/kmer-diversity-sieve.git` shows `main` at `70022da`.

**Proposed fix.** In a separate commit, run `git -C extern/kmer-diversity-sieve checkout 70022da && git add extern/kmer-diversity-sieve`. Leave the recipe alone. In backlog (e), delete the "rewrote history" sentences and the options "move the recipe" and "move the submodule back".

**Tests and docs.** Add a CI step that compares `git ls-tree HEAD extern/kmer-diversity-sieve` with the commit named in `docs/packaging.md`.

### PK-09 · low · Docs attribute the MAFFT 7.525 bug to `--add`; `msa` uses `--addfragments`, which 7.526 fixed

`docs/packaging.md:39` · docs · confirmed (by reading) · fix effort: trivial · sources: PK2-11

**Problem.** `docs/packaging.md:39` and `docs/galaxy.md:22` say MAFFT 7.525 segfaults "on the `--add` path" that `msa` uses. `src/disassembler/cmd/msa.py:284` passes `--addfragments`. The only functional change from 7.525 to 7.526 is upstream commit b4cc76f0, "Fixed a bug in --addfragments" (`core/pairlocalalign.c`). Someone re-testing the pin with plain `--add` would not hit the bug and might drop the floor. Evidence: `curl -sS 'https://gitlab.com/api/v4/projects/sysimm%2Fmafft/repository/compare?from=v7.525&to=v7.526'`.

**Proposed fix.**
```diff
-- **MAFFT 7.526.** MAFFT 7.525 segfaults on the `--add` path that `disassembler msa` uses.
+- **MAFFT 7.526.** MAFFT 7.525 has a bug in `--addfragments` (the pairlocalalign stage), which
+  `disassembler msa` uses; it is fixed in 7.526 (upstream commit b4cc76f0).
```
Use the same wording in `docs/galaxy.md:22`, in the recipe comment for the floor, in plan 1(f) and in backlog line 59.

**Tests and docs.** Optionally, make `tests/test_msa_real_tools.py` skip with a clear message on MAFFT 7.525.

### PK-10 · low · `environment.yml` still installs zstandard, which nothing imports

`environment.yml:31` · packaging · confirmed (reproduced) · fix effort: trivial · sources: PK1-08

**Problem.** Commit 25c27f2 removed zstandard from `pyproject.toml` ("nothing imports zstandard at all", `pyproject.toml:14`). `environment.yml:31` still lists it under "Python deps for the disassembler package and the benchmark scripts". The development environment therefore has a module that the Bioconda install lacks, so code that imports it would pass locally and fail in the package.

**Proposed fix.**
```diff
-  # Python deps for the disassembler package and the benchmark scripts.
+  # Python deps: the package's runtime imports (pyproject `dependencies`)...
   - biopython
   - requests
-  - zstandard
+  # ...and the benchmark scripts (pyproject `bench` extra).
   - numpy
```
Keep `zstd=1.5.6`, which is the command-line tool.

**Tests and docs.** Optionally, add a unittest asserting that the unpinned Python packages in the yml equal `dependencies` plus the `bench` extra.

### PK-11 · low · Backlog item 0's "stock install has a broken `msa`" holds only with bioconda ahead of conda-forge

`backlog:61` · backlog-claim · confirmed (by reading) · fix effort: trivial · sources: PK2-02

**Problem.** Backlog line 61 ("A stock install of the published package has a broken `msa`") and the backlog description on line 3 rest on a dry-run with `-c bioconda -c conda-forge`. The outcome depends on channel order and priority:

| Command | Channel priority | Result |
|---|---|---|
| `-c conda-forge -c bioconda` (README, `docs/packaging.md:20`) | any | mafft 7.526 |
| `-c bioconda -c conda-forge` | flexible | mafft 7.525 |
| `-c bioconda -c conda-forge` | strict | does not solve (bioconda's biopython 1.70 shadows conda-forge's) |

Inspection of the published images found conda-forge mafft 7.526 (re-checked for the 0.1.1 disassembler image, see PK-01), which answers todo 4. The argument in todo 3 that "a bad solve passes CI" does not hold, because Bioconda's test solve puts conda-forge first.

**Proposed fix.** Reword lines 3 and 61: "With bioconda ahead of conda-forge and flexible priority, the unpinned recipe resolves mafft 7.525 and `msa` can fail. The documented order and the published images get 7.526." Keep `mafft >=7.526` as hardening (PK-02). Call the `msa` recipe test a smoke test; it needs `test: source_files:` for its fixtures. Strike todo 4.

**Tests and docs.** Backlog text only.

### PK-12 · low · Backlog 1(a) understates provenance work: cached and offline runs have no headers, manifests are not free

`walker/src/stages/stream.rs:102` · backlog-claim · confirmed (reproduced) · fix effort: small · sources: BL1-10, PK2-10

**Problem.** Backlog 1(a) scopes per-accession provenance as reading ETag and Last-Modified in `try_once`, "~10 lines". The Galaxy plan 1(a) says a manifest gives release-level provenance "for free". Neither holds:
- The cache-hit branch returns at `stream.rs:102` without any HTTP response, and `write_cache_atomic` (`stream.rs:249`) stores only the bytes. So `--unitigs-cache --offline`, the reproducible route, would record nothing.
- The data must be carried from `FetchedUnitigs` through to `AccessionResult` and `write_result`.
- The `stats/logan_accessions_v*` manifests are 5.1 GB and 11.8 GB SRA-metadata tables, and `u/` has no versioned prefix, so naming a manifest pins none of the bytes the walker reads.

HEAD requests showed that each `u/` object's ETag equals the MD5 of its bytes, and that Last-Modified is a bulk-upload date (31 Mar 2024).

**Proposed fix.** Rescope 1(a). Then:
- Add a `Provenance { source, md5, sha256, etag, last_modified, fetched_at }` field to `FetchedUnitigs` and carry it to `AccessionResult`.
- Hash the bytes on both the fetch path and the cache path.
- Write a `<ACC>.unitigs.fa.zst.meta.json` sidecar atomically next to the cache file.
- Emit a `"source"` object in `summary.json`.
- Record one `ListObjectsV2 prefix=stats/logan_accessions_` listing in the run settings.

This is about 150 to 250 lines with tests. Multipart ETags (`-N` suffix) are not MD5s.

**Tests and docs.** Unit-test a cache hit with and without the sidecar (MD5 of the on-disk bytes, null or present ETag), and a local test server returning an ETag that must reach `summary.json`. Update the backlog 1(a) estimate, plan 1(a) and the `summary.json` schema in `docs/walker.md`.

### PK-13 · low · Plan 1(b): O(subgraph) memory needs two decodes and a spooled adjacency, not one pass

`walker/src/stages/extract.rs:93` · backlog-claim · confirmed (reproduced) · fix effort: large · sources: BL1-04

**Problem.** Galaxy plan 1(b) (lines 54-56) says that decoding incrementally and reducing to the anchor-proximal subgraph "during the pass" cuts memory from O(all unitigs) to O(subgraph). A single pass cannot do this:
- `extract::expand` follows `graph.edges(node)` (`extract.rs:93`) to arbitrary unitig ids, and Logan's unitig order is unrelated to graph distance. A unitig that streams past before a kept neighbour asks for it is lost.
- Bounded memory therefore needs two decodes: pass 1 keeps the adjacency (spooled to disk) and the anchors, a BFS runs between the passes, and pass 2 collects the sequences.
- Between passes the compressed object must be kept or fetched again; today it is dropped at `walker/src/pipeline.rs:536`.

If the adjacency stays in RAM, memory is still O(all unitigs), only smaller.

**Evidence.** A simulated best single pass, compared with the walker's own GFA, kept 18,424 of 34,582 nodes on SRR8209082/hiv1_pol and 289 of 747 on ERR2028820/pfcrt. A two-decode prototype produced byte-identical GFAs with 12.6 MiB peak RSS when the adjacency was spooled.

**Proposed fix.** Rewrite plan 1(b) and backlog (b) as a two-decode design, and add scratch disk (up to 1 GB per accession) to the resource model. Code changes:
- Split `parse::build` into a scan pass and a collect pass.
- Give downstream stages a subgraph struct instead of `&UnitigGraph`.
- Keep the unpruned node list for the refused-prune branch (`pipeline.rs:843`).
- Gate the change on GFA byte identity with the old path.

**Tests and docs.** A data-backed (`#[ignore]`) test comparing GFA bytes between the old and new extract paths on the CLI fixtures, ERR2028820 and SRR8209082, plus `bench/check_golden.sh`. Update the memory-model notes in `extract.rs`, `pipeline.rs`, `config.rs` and `docs/walker.md`.

### PK-14 · low · Backlog calls Logan objects unindexed zstd; they are seekable with a seek table and record index

`backlog:1(b)` (line 148) · backlog-claim · confirmed (reproduced) · fix effort: trivial · sources: R2-07

**Problem.** Backlog 1(b) rules out range fetches because "the object is a flat unindexed `.zst`". The objects are multi-frame zstd: about 128 MiB of decompressed data per frame, with frames ending on record boundaries. Two trailing skippable frames follow:
- frame `0x184D2A5E` holds a standard zstd seek table;
- frame `0x184D2A5F`, ending in `XDIF`, lists each frame's first record name.

The conclusion still holds, since there is no content index. But the stated reason rules out options the format allows: frame-level random access, parallel decoding, knowing the exact decompressed size before decoding, and detecting truncation. Evidence: a ranged fetch of the last 400 bytes (`curl -r -400` on `https://logan-pub.s3.amazonaws.com/u/ERR2028820/ERR2028820.unitigs.fa.zst`) shows 2 frames, with first records `ERR2028820_0` and `ERR2028820_1164640`.

**Proposed fix.** Correct the sentence and list those options in 1(b). The layout is not a published Logan contract, so any code that uses it must fall back to a full stream when the tail is missing.

**Tests and docs.** Backlog text only; if 1(b) uses the seek table, add a fixture with a synthetic multi-frame object and seek table.

### PK-15 · info · Backlog's ~34 B/unitig CSR figure is right for one accession but varies with graph shape

`walker/src/types.rs:97` · backlog-claim · confirmed (reproduced) · fix effort: trivial · sources: R2-10

**Problem.** The memory held by the CSR arrays is 12 + 4 × (links per unitig) + (bases per unitig)/4 bytes per unitig (`types.rs:101-125`). That is 33.5 B on ERR2028820 but 38.2 B on SRR8209082, a link-dense viral graph. On SRR8209082 the arrays alone take 1.67 × the compressed size, close to `GRAPH_PEAK_PER_COMPRESSED_BYTE = 1.7` (`walker/src/config.rs:445`), a constant that is meant to include Vec growth slack as well.

**Proposed fix.** In backlog 1(b), write "~34 B/unitig on ERR2028820 (12 + 4 × links/unitig + bases/4; 38 B on a link-dense viral graph)". Measure peak RSS on a large link-dense accession before relying on 1.7 as a bound.

**Tests and docs.** None; backlog wording only.
