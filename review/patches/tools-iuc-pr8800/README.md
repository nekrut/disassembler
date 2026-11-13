# Fix proposals for the Galaxy wrappers of tools-iuc PR #8800

These four patches are for the wrapper files of galaxyproject/tools-iuc pull request #8800, not for this repository. They are rooted at the tools-iuc repository root (`tools/disassembler/...`) and each applies on its own to the PR head `b9a16ad` with `git apply -p1` (checked with `git apply --check` against the PR's files). They are proposals for review, not merged changes. The findings they belong to are in [`../../galaxy-wrappers.md`](../../galaxy-wrappers.md).

| Patch | Finding | What it changes |
|---|---|---|
| `GX-01-walker-stage-reference-in-own-directory.patch` | GX-01 | Stages the reference dataset in its own directory in `disassembler_walker.xml`, so a reference named like an output file can no longer be overwritten. |
| `GX-02-sieve-max-gap-frac.patch` | GX-02 | Exposes `--max-gap-frac` in `disassembler_sieve.xml`, so LexicMap cohorts keep partly covered rows the composite keeps. |
| `GX-04-walker-primary-metric-matched-default.patch` | GX-04 | Makes `matched` the default `--primary-metric` in `disassembler_walker.xml` and `macros.xml`, keeps `cov` selectable, and adds a test. |
| `GX-10-bump-to-0.1.1.patch` | GX-10 | The changes a bump to 0.1.1 needs beyond the version tokens: stratify and walker wrappers, `macros.xml`, and the two score-file test fixtures. |

GX-04 and GX-10 both edit `disassembler_walker.xml` and `macros.xml`; apply GX-10 first and rebase GX-04 onto it.
