# Archived: Python-walker tooling

These scripts targeted the original Python reference walker (`python/logan_path_walker.py`),
which has been removed now that the Rust `logan-walker` is the sole implementation and the
golden gate (`../check_golden.sh`) verifies it against committed fixtures without needing the
Python walker.

They are kept for the record, not as working tools — each depends on the removed walker:

| script | what it did | replaced by |
|---|---|---|
| `verify_golden.py` | ran `capture_golden.py` and diffed against `../golden/` | `../check_golden.sh` (Rust-native, byte-exact on ERR2028820 + SRR12363134) |
| `capture_golden.py` | captured a golden fixture from the Python walker | fixtures are now frozen; regenerate from the Rust walker if ever needed |
| `measure_fork_loss.py` | measured the reference-bias fork loss (a copy of the Python stages) | findings written up in `../../docs/issue_intra_accession_variance_lost.md` |
| `run_benchmark.sh` | timed the Python walker over 50 accessions | benchmark the Rust binary directly (`logan-walker -l ... -q ...`) |

The golden fixtures themselves (`../golden/*.primary.fasta`, `*.golden.json`) remain live — they
are the ground truth the Rust gate checks against.
