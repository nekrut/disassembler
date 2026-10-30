# Packaging and releases

The repository builds three independent artifacts, a Python package, a Rust binary and a C
binary, and each is its own Bioconda package:

| package | built from | provides | depends on |
|---|---|---|---|
| `disassembler` | `src/disassembler/` (`noarch: python`) | the `disassembler` CLI | `logan-walker`, `kmer-diversity-sieve`, `mafft`, `fasttree`, plus Python deps |
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

All three build for `linux-64`, and `kmer-diversity-sieve` also for `linux-aarch64`. GraphAligner
has no ARM builds, so `logan-walker` cannot have them either. The macOS builds were switched off
in `0.1.0` because of a Bioconda CI problem
([bioconda-recipes#70401](https://github.com/bioconda/bioconda-recipes/issues/70401), since
fixed); re-enabling them is planned for the next release.

## Pinned tools

- **GraphAligner 1.0.17b.** 1.0.19 and 1.0.20 drop the last query base of every alignment, and
  1.0.18 cannot be installed alongside Python 3.11 or later. The evidence and the conditions for
  revisiting the pin are in
  [`issue_graphaligner_end_clipping.md`](issue_graphaligner_end_clipping.md). The golden
  fixtures in `bench/golden/` were captured under 1.0.17b.
- **MAFFT 7.526.** MAFFT 7.525 segfaults on the `--add` path that `disassembler msa` uses.
  `environment.yml` and the Galaxy `msa` tool pin 7.526; the `0.1.0` Bioconda recipe does not
  set a floor yet.
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
- `disassembler` depends on a minor version range of the walker and sieve (`>=0.1.0,<0.2`).

## Making a release

1. Bump the version in `src/disassembler/__init__.py`, and in `walker/Cargo.toml` (and
   `walker/Cargo.lock`) if the walker changed.
2. Run the checks: `cargo clippy --all-targets --release -- -D warnings` and
   `cargo test --release` in `walker/`, the Python tests in `tests/`, and the golden gate
   (`bench/check_golden.sh`, which needs GraphAligner 1.0.17b and network access).
3. Tag the release commit `v<version>`, and `logan-walker-v<version>` if the walker changed.
4. Update the Bioconda recipes (`recipes/disassembler`, `recipes/logan-walker`): version,
   `sha256` of the new tag tarball, and build number reset to 0.
5. Update the version tokens in the Galaxy suite's `macros.xml` (`@TOOL_VERSION@`,
   `@WALKER_VERSION@`).
