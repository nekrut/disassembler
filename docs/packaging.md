# Packaging and releases

The repository builds three independent artifacts, a Python package, a Rust binary and a C
binary, and each is its own Bioconda package:

| package | built from | provides | depends on |
|---|---|---|---|
| `disassembler` | `src/disassembler/` (`noarch: python`) | the `disassembler` CLI | `logan-walker`, `kmer-diversity-sieve`, `mafft`, `fasttree`, plus biopython and requests |
| `logan-walker` | `walker/` | `logan-walker` | `graphaligner ==1.0.17b` |
| `kmer-diversity-sieve` | [stevenweaver/kmer-diversity-sieve](https://github.com/stevenweaver/kmer-diversity-sieve) | `l1_pipeline_serial` | |

```
                 disassembler
           /        |          \          \
   logan-walker  kmer-diversity-sieve  mafft  fasttree
        |
   graphaligner
```

`conda install -c conda-forge -c bioconda disassembler` therefore installs everything the CLI
calls. The walker and the sieve are packaged separately because each is useful on its own, and
so that a Galaxy tool can depend on only the piece it runs.

## Platforms

All three build for `linux-64`, and `kmer-diversity-sieve` also for `linux-aarch64` and
`osx-arm64`. GraphAligner has no ARM builds, so `logan-walker` cannot have them either.

`logan-walker`, and therefore `disassembler`, are **linux-64 only**, and not because of the
Bioconda CI problem the `0.1.0` recipe cited
([bioconda-recipes#70401](https://github.com/bioconda/bioconda-recipes/issues/70401), since
closed). GraphAligner 1.0.17b has no `osx-arm64` build at all, and both of its `osx-64` builds
(`h66d0458_0`, `h8d12793_2`) depend on `protobuf 3.14.0`, which conda-forge ships for `osx-64`
only up to Python 3.9, while `disassembler` needs Python >= 3.11. Linux solves only because of
the protobuf-free build `hdcf5f25_1`. Publishing an `osx-64` `logan-walker` would pass its own
CI, because on its own it can take Python 3.9, and still leave `disassembler` unsatisfiable on
macOS. macOS returns once an equivalent protobuf-free `osx-64` GraphAligner 1.0.17b build
exists.

Before enabling a platform, check that it actually solves:
`conda create --dry-run --platform osx-64 -c conda-forge -c bioconda disassembler`.

## Known bad builds

- **The biocontainer `disassembler:0.1.1--pyhdfd78af_0` bundles `logan-walker 0.1.0`.** Its
  build-time solve ran seconds after `logan-walker 0.1.1` was uploaded, saw a stale channel
  index, and the recipe's floor at the time (`>=0.1.0,<0.2`) allowed it. The `0.1.1` build 0
  conda package is not itself wrong, it only carries that permissive floor, so a fresh
  `conda install` today gets `logan-walker 0.1.1`; it is the image that has 0.1.0 baked in. Anything running `walker`, `fetch`, `reconstruct` or
  `from-kmindex` from that image gets 0.1.0 behaviour, the `cov` primary metric and no
  prune-span guard, and 0.1.0's `summary.json` has no `logan_walker` block to reveal it. Use
  build `_1` or later. A published build cannot be withdrawn, only superseded.

## Pinned tools

- **GraphAligner 1.0.17b.** 1.0.19 and 1.0.20 drop the last query base of every alignment, and
  1.0.18 cannot be installed alongside Python 3.11 or later. The evidence and the conditions for
  revisiting the pin are in
  [`issue_graphaligner_end_clipping.md`](issue_graphaligner_end_clipping.md). The golden
  fixtures in `bench/golden/` were captured under 1.0.17b.
- **MAFFT 7.526.** MAFFT 7.525 has a bug in `--addfragments` (the pairlocalalign stage), which
  is the path `disassembler msa` uses; it is fixed in 7.526 (upstream commit `b4cc76f0`). Note
  the flag precisely: re-testing with a plain `--add` does not reproduce it. `environment.yml`
  and the Galaxy `msa` tool pin 7.526. Neither the `0.1.0` recipe nor the published `0.1.1`
  build 0 sets a floor; `mafft >=7.526` arrives with `0.1.1` build 1. Until that is published,
  list conda-forge before bioconda, or run `conda install 'mafft>=7.526'` in the environment.
  Channel order decides this: `-c` prepends, so `-c conda-forge -c bioconda` gets 7.526 while
  `-c bioconda -c conda-forge` with flexible priority gets 7.525.
- **kmer-diversity-sieve 0.1.0** is upstream commit `70022da`, because upstream has no tagged
  releases. The recipe overrides the Makefile's `-march=native`, so the binary runs on any
  x86-64 machine; the portable build gave byte-identical output on the test fixture.

## Versions and tags

- `disassembler`'s version is `__version__` in `src/disassembler/__init__.py`; `disassembler
  --version` prints it.
- `logan-walker`'s version is in `walker/Cargo.toml`; `logan-walker --version` prints it.
- The two have separate tag series: `v<version>` for `disassembler`, and
  `logan-walker-v<version>` for the walker. The Bioconda recipes download those tags'
  tarballs.
- `disassembler` depends on a minor version range of the walker and sieve: the walker at
  `>=0.1.1,<0.2` from `0.1.1` build 1 onward, the sieve at `>=0.1.0,<0.2`.

## Making a release

1. Bump the version in `src/disassembler/__init__.py`, and in `walker/Cargo.toml` (and
   `walker/Cargo.lock`) if the walker changed.
2. Run the checks: `cargo clippy --all-targets --release -- -D warnings` and
   `cargo test --release` in `walker/`, the Python tests in `tests/`, and the golden gate
   (`bench/check_golden.sh`, which needs GraphAligner 1.0.17b and network access).
3. Tag the release commit `v<version>`, and `logan-walker-v<version>` if the walker changed.
4. **Open the recipe PR immediately**, in the same sitting as the tag. Bioconda's autobump bot
   opens its own version-only PR within minutes of a tag appearing, and that PR cannot be
   extended: its branch lives in `bioconda/bioconda-recipes`, where a recipe maintainer has
   `pull` only. If it merges first, the release ships with the previous recipe body and the
   fixes then need a build-number bump on an already-published version. All three recipes
   carry `extra: autobump: enable: False` to prevent this; check it is still there.

   Check every line of the recipe, not just the version:

   - `version` and the `sha256` of the **GitHub tag tarball**, not a `git archive`:
     `curl -sL https://github.com/nekrut/disassembler/archive/refs/tags/<tag>.tar.gz | sha256sum`
   - build number: `0` for a version never published, **one higher** than the published build
     if that version already exists in the channel. It is not "reset to 0".
   - `logan-walker >=<new walker version>,<next minor>`
   - run dependencies equal to `pyproject.toml`'s `[project] dependencies`, plus the binaries
     the subcommands shell out to. Nothing else.
   - `mafft >=7.526`
   - the sieve commit equals the `extern/kmer-diversity-sieve` submodule pin
     (`git ls-tree HEAD extern/kmer-diversity-sieve`)
5. Bump `main` to the next development version (`__version__ = "<next>.dev0"`, and the walker's
   `Cargo.toml` if it moved), so builds from `main` stop claiming the released version.
6. **After the recipe PR merges and the packages build**, verify what was actually published:

   - `conda create --dry-run -c conda-forge -c bioconda disassembler=<version>` resolves the
     expected `logan-walker` and `mafft`.
   - the biocontainer carries them too. Read its `conda-meta` records rather than trusting the
     solve; the image is built by its own solve, which can race the channel index.
7. Update the version tokens in the Galaxy suite's `macros.xml` (`@TOOL_VERSION@`,
   `@WALKER_VERSION@`), once tools-iuc #8800 is merged. Until then a version bump there is not
   just a token change.
