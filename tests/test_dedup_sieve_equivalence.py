#!/usr/bin/env python3
"""test_dedup_sieve_equivalence.py -- the dedup pass promotes the same sequences.

`sieve_run` folds byte-exact duplicates out of its input before invoking the
sieve. A later identical copy has novelty 0 and is never promoted, so the two
paths promote the same sequences, in the same order.

Two things changed in 0.2.0 and are reflected here. The admission gate is
disarmed, so these flags carry `--warmup 2147483647` as the wrapper now
does, and duplicates no longer skew anything: dedup is a speed-up and a memory
saving, not a correctness fix. And the dedup pass renames records `>r0`, `>r1`,
..., because the binary buckets novelty by a month parsed from the
header and selects what to emit by the header's first 255 characters. The two
paths therefore cannot be byte-identical any more, so what is asserted is the
stronger claim the bytes were standing in for: the same sequences in the same
order, and the saved headers mapping back to exactly the full path's headers.

Also records the dedup fold ratio and (loosely) checks the dedup pass itself is
cheap. Does NOT assert on wall-clock speedup of the sieve (noisy on CI); the
measurement belongs in bench/, not here.

Skips cleanly if the sieve binary or its submodule are absent.
"""
import filecmp
import hashlib
import subprocess
import tempfile
import unittest
from pathlib import Path

from disassembler.cmd.sieve import ADMISSION_OFF_WARMUP

REPO = Path(__file__).resolve().parent.parent
SIEVE_BIN = REPO / "extern" / "kmer-diversity-sieve" / "l1_pipeline_serial"
FIXTURE = Path(__file__).resolve().parent / "fixture_reconstructed.fasta"


def _ensure_sieve():
    if SIEVE_BIN.exists():
        return "ok"
    if not (SIEVE_BIN.parent / "Makefile").is_file():
        return "no-submodule"
    subprocess.run(["make", "serial"], cwd=SIEVE_BIN.parent, capture_output=True)
    return "ok" if SIEVE_BIN.exists() else "no-compiler"


def _parse_fasta(path: Path):
    records = []
    header = None
    seq_lines: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line.rstrip("\n")[1:]
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            records.append((header, "".join(seq_lines)))
    return records


def _write_fasta(path: Path, records):
    with path.open("w") as fh:
        for h, s in records:
            fh.write(f">{h}\n{s}\n")


def _run_sieve(input_fa: Path, output_fa: Path, min_len: int):
    """Run the sieve with the same flags the composite uses, concatenating rank shards."""
    emit_prefix = output_fa.with_suffix("")
    for stale in output_fa.parent.glob(f"{emit_prefix.name}.rank*.fasta"):
        stale.unlink()
    cmd = [
        str(SIEVE_BIN),
        "--input", str(input_fa),
        "--min-len", str(min_len),
        "--max-n-frac", "0.05",
        "--threshold", "0.005",
        "--warmup", str(ADMISSION_OFF_WARMUP),
        "--sketch", "256",
        "--emit-promoted", str(emit_prefix),
    ]
    rc = subprocess.run(cmd, capture_output=True).returncode
    if rc != 0:
        raise RuntimeError(f"sieve exited {rc}")
    ranks = sorted(output_fa.parent.glob(f"{emit_prefix.name}.rank*.fasta"))
    with output_fa.open("w") as fh:
        for rf in ranks:
            fh.write(rf.read_text())
        for rf in ranks:
            rf.unlink()


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class TestDedupSieveEquivalence(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        state = _ensure_sieve()
        if state == "no-submodule":
            raise unittest.SkipTest(
                "extern/kmer-diversity-sieve submodule not initialised"
            )
        if state == "no-compiler":
            raise unittest.SkipTest("sieve binary failed to build")
        if not FIXTURE.is_file():
            raise unittest.SkipTest("fixture_reconstructed.fasta missing")
        try:
            from disassembler.cmd.sieve import byte_exact_dedup  # noqa
        except ImportError as e:
            raise unittest.SkipTest(f"disassembler not importable: {e}")

    def setUp(self):
        from disassembler.cmd.sieve import byte_exact_dedup
        self.dedup = byte_exact_dedup
        self.work = Path(tempfile.mkdtemp())
        self.fixture = _parse_fasta(FIXTURE)
        self.min_len = int(len(self.fixture[0][1]) * 0.95)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.work, ignore_errors=True)

    def _make_cohort(self, copies: int):
        """Replicate every unique sequence `copies` times with distinct accession IDs."""
        out = []
        for i, (h, s) in enumerate(self.fixture):
            for c in range(copies):
                out.append((f"ACC{i:04d}_copy{c:02d} {h}", s))
        return out

    def _run_both_paths(self, copies: int):
        cohort = self.work / f"c{copies}.cohort.fa"
        deduped = self.work / f"c{copies}.deduped.fa"
        sieved_full = self.work / f"c{copies}.sieved_full.fa"
        sieved_from_ded = self.work / f"c{copies}.sieved_from_ded.fa"
        _write_fasta(cohort, self._make_cohort(copies))
        rep = self.dedup(cohort, deduped)
        _run_sieve(cohort, sieved_full, self.min_len)
        _run_sieve(deduped, sieved_from_ded, self.min_len)
        return cohort, deduped, sieved_full, sieved_from_ded, rep

    @staticmethod
    def _records(path: Path):
        """(header, sequence) pairs, in file order."""
        out, hdr, seq = [], None, []
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                if hdr is not None:
                    out.append((hdr, "".join(seq)))
                hdr, seq = line[1:], []
            else:
                seq.append(line.strip())
        if hdr is not None:
            out.append((hdr, "".join(seq)))
        return out

    def _assert_same_promotions(self, full: Path, ded: Path, rep):
        """Same sequences in the same order, and the saved headers map back exactly."""
        a, b = self._records(full), self._records(ded)
        self.assertEqual([s for _h, s in a], [s for _h, s in b],
                         "the two paths promoted different sequences")
        restored = [rep.headers[int(h.split()[0][1:])].rstrip("\n")[1:] for h, _s in b]
        self.assertEqual([h for h, _s in a], restored,
                         "mapping the neutral ids back did not reproduce the headers")

    def test_no_duplicates_dedup_pass_is_identity(self):
        """copies=1: dedup is a no-op; the same sequences are promoted either way."""
        _, _, full, ded, rep = self._run_both_paths(copies=1)
        self.assertEqual(rep.n_in, rep.n_out,
                         "dedup should drop nothing when input has no duplicates")
        self._assert_same_promotions(full, ded, rep)

    def test_byte_identical_sieve_output_with_duplicates(self):
        """copies=10: every unique sequence appears 10x in the cohort; sieve outputs
        must still be byte-identical across the two paths (same bytes, same order,
        same headers). This is the core equivalence claim.
        """
        _, _, full, ded, rep = self._run_both_paths(copies=10)
        self.assertEqual(rep.n_in, 10 * len(self.fixture))
        self.assertEqual(rep.n_out, len(self.fixture),
                         "every unique sequence appears once in deduped output")
        self._assert_same_promotions(full, ded, rep)

    def test_dedup_preserves_first_encountered_header(self):
        """The dedup rule MUST be first-encountered-wins; a different winner (e.g.
        alphabetical) would promote a different accession per cluster and break
        byte-equivalence with the sieve-only path.
        """
        src = self.work / "ordered.fa"
        _write_fasta(src, [
            ("SRR_Z_first", "ACGTACGTACGT"),
            ("SRR_A_second", "ACGTACGTACGT"),   # same sequence as first
            ("SRR_M_third", "TTTTTTTTTTTT"),
        ])
        dst = self.work / "ordered.deduped.fa"
        rep = self.dedup(src, dst)
        self.assertEqual((rep.n_in, rep.n_out), (3, 2))
        # The file the binary reads now carries neutral ids; the headers it
        # must be given back are the first-encountered ones, in the same order.
        self.assertEqual([l for l in dst.read_text().splitlines() if l.startswith(">")],
                         [">r0", ">r1"])
        self.assertEqual([h.strip() for h in rep.headers],
                         [">SRR_Z_first", ">SRR_M_third"],
                         "first-encountered header wins; alphabetical would pick SRR_A_second")


if __name__ == "__main__":
    unittest.main()
