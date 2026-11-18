# Packaging

The repository builds three independent artifacts, a Python package, a Rust binary and a C
binary, and each is its own Bioconda package:

| package | built from | provides | depends on |
|---|---|---|---|
| `disassembler` | `src/disassembler/` (`noarch: python`) | the `disassembler` CLI | `logan-walker`, `kmer-diversity-sieve`, `mafft`, `fasttree`, plus biopython and requests |
| `logan-walker` | `walker/` | `logan-walker` | `graphaligner ==1.0.17b` |
| `kmer-diversity-sieve` | [stevenweaver/kmer-diversity-sieve](https://github.com/stevenweaver/kmer-diversity-sieve) | `l1_pipeline_serial` | |

```text
                 disassembler
           /        |          \          \
   logan-walker  kmer-diversity-sieve  mafft  fasttree
        |
   graphaligner
```

So `conda install -c conda-forge -c bioconda disassembler` installs everything the CLI calls.
The walker and the sieve are packaged separately because each is useful on its own, and so that
a Galaxy tool can depend on only the piece it runs.

## Platforms

`kmer-diversity-sieve` builds for `linux-64` and `linux-aarch64`. `logan-walker`, and therefore
`disassembler`, are **linux-64 only**, because of GraphAligner:

- 1.0.17b has no `osx-arm64` build at all;
- both of its `osx-64` builds (`h66d0458_0`, `h8d12793_2`) depend on `protobuf 3.14.0`, which
  conda-forge ships for `osx-64` only up to Python 3.9, while `disassembler` needs Python 3.11
  or later. Linux solves only because of the protobuf-free build `hdcf5f25_1`.

An `osx-64` `logan-walker` would pass its own CI, since on its own it can take Python 3.9, and
still leave `disassembler` unsatisfiable on macOS. macOS returns once an equivalent
protobuf-free `osx-64` GraphAligner 1.0.17b build exists. The `kmer-diversity-sieve` recipe
skips macOS too, by choice: it is packaged for `disassembler`, which cannot be installed there.
Before enabling a platform, check that it solves:

```bash
conda create --dry-run --platform osx-64 -c conda-forge -c bioconda disassembler
```

## Pinned tools

- **GraphAligner 1.0.17b.** 1.0.19 and 1.0.20 drop the last query base of every alignment, and
  1.0.18 cannot be installed alongside Python 3.11 or later. The evidence and the conditions for
  revisiting the pin are in
  [`notes/graphaligner_end_clipping.md`](notes/graphaligner_end_clipping.md). The golden
  fixtures in `bench/golden/` were captured under 1.0.17b.
- **MAFFT 7.526.** 7.525 has a bug in `--addfragments` (the pairlocalalign stage), which is the
  path `disassembler msa` uses; it is fixed in 7.526 (upstream commit `b4cc76f0`). The flag
  matters: re-testing with a plain `--add` does not reproduce it. `environment.yml`, the recipe
  and the Galaxy `msa` tool all pin `mafft >=7.526`, and `msa` runs `mafft --version` before any
  work and refuses 7.525 by name.
- **kmer-diversity-sieve 0.1.0** is upstream tag `v0.1.0` (commit `70022da`). The recipe
  overrides the Makefile's `-march=native`, so the binary runs on any x86-64 machine; the
  portable build gives byte-identical output on the test fixture.

The walker floor is enforced at run time as well as by the recipe: the package runs
`logan-walker --version` before using it and refuses anything below `MIN_WALKER_VERSION`
(`src/disassembler/resolve.py`), in the `walker` passthrough and in the composites' preflight.
That catches a stale conda environment, a `--walker-bin` pointed at an old build and a
hand-built binary on `PATH`, none of which a recipe floor constrains. Raise the constant and the
recipe floor together.

## Versions and tags

- `disassembler`'s version is `__version__` in `src/disassembler/__init__.py`; `disassembler
  --version` prints it.
- `logan-walker`'s version is in `walker/Cargo.toml`; `logan-walker --version` prints it.
- The two have separate tag series, `v<version>` for `disassembler` and
  `logan-walker-v<version>` for the walker. The Bioconda recipes download those tags' tarballs.
- `kmer-diversity-sieve` is tagged `v<version>` in its own repository, and its recipe downloads
  that tag's tarball.
- `disassembler` depends on a minor version range of each binary: the walker at
  `>=<current minor>,<next minor>`, the sieve at `>=0.1.0,<0.2`.

## Making a release

1. Bump the version in `src/disassembler/__init__.py`, and in `walker/Cargo.toml` (and
   `walker/Cargo.lock`) if the walker changed.
2. Run the checks: `cargo clippy --all-targets --release -- -D warnings` and
   `cargo test --release` in `walker/`, the Python tests in `tests/`, and the golden gate
   (`bench/check_golden.sh`, which needs GraphAligner 1.0.17b and network access). The gate also
   runs in CI on pushes to `main` and on release tags (`.github/workflows/golden.yml`); record
   the URL of a green run, because a local run leaves no record.
3. Tag the release commit `v<version>`, and `logan-walker-v<version>` if the walker changed.
4. Update the Bioconda recipes. Bioconda's autobump bot opens a pull request for a new tag
   that changes only the version and checksum, and it merges unreviewed after three days if CI
   passes. The three recipes are set up differently:

   - `disassembler` has the bot switched off (`extra: autobump: enable: False`), because its
     run dependencies and walker floor change between releases. Nothing is published until you
     open the recipe PR yourself.
   - `logan-walker` and `kmer-diversity-sieve` leave the bot on. A `logan-walker-v<version>`
     tag here, or a `v<version>` tag in the sieve's repository, is published by the bot without
     further action. That is right when only the source changed. If the recipe body has to
     change as well (a different GraphAligner pin, build flags, a renamed binary), open your own
     PR in the same sitting as the tag, before the bot's merges.

   Any edit to a recipe whose current build is already published needs a build number bump,
   even a comment-only edit; Bioconda's linter rejects it otherwise. For a recipe PR you write
   by hand, check every line, not just the version:

   - `version`, and the `sha256` of the **GitHub tag tarball** rather than a `git archive`:
     `curl -sL https://github.com/nekrut/disassembler/archive/refs/tags/<tag>.tar.gz | sha256sum`
   - build number: `0` for a version never published, **one higher** than the published build
     if that version already exists in the channel. It is not reset to 0.
   - `logan-walker >=<new walker version>,<next minor>`
   - run dependencies equal to `pyproject.toml`'s `[project] dependencies`, plus the binaries
     the subcommands shell out to, and nothing else
   - `mafft >=7.526`
   - the sieve version's tag points at the `extern/kmer-diversity-sieve` submodule pin
     (`git ls-tree HEAD extern/kmer-diversity-sieve`)
5. Bump `main` to the next development version (`__version__ = "<next>.dev0"`, and the walker's
   `Cargo.toml` if it moved), so builds from `main` stop claiming the released version.
6. Once the recipe PR has merged and the packages have built, verify what was published:
   `conda create --dry-run -c conda-forge -c bioconda disassembler=<version>` resolves the
   expected `logan-walker` and `mafft`, and the biocontainer carries them too. Read the image's
   `conda-meta` records rather than trusting a fresh solve, since the image is built by its own
   solve, which can race the channel index.
7. Update the version tokens in the Galaxy suite's `macros.xml` (`@TOOL_VERSION@`,
   `@WALKER_VERSION@`).
