#!/usr/bin/env python3
"""test_sieve_integration.py -- exercise `disassembler reconstruct --sieve` end to end.

The bash test (test_sieve.sh) covers the sieve BINARY. This covers the Python SEAM
that test can not: it stubs stage 1 (the walker, which needs GraphAligner + S3 and
so can not run on every host) with a fake `--walker-bin`, and it stubs stage 2
(codon MSA + FastTree) by monkey-patching `disassembler.cmd.msa.run`. Then it runs
`pipe_reconstruct.run(...)` in-process and asserts:

  - the sieve stage ran on the walker's reconstructed.fasta,
  - stage 2 received the SIEVED file (not the raw one),
  - the `--sieve-sketch` flag reaches the sieve binary CLI (a regression that
    would let the flag look-like-it-works while the binary ran at its default of 64,
    which is precisely the failure mode the gene-scale calibration doc characterises).

Skips cleanly if the sieve binary or its submodule are absent.
"""
import contextlib
import io
import os
import shutil
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SIEVE_BIN = REPO / "extern" / "kmer-diversity-sieve" / "l1_pipeline_serial"
FIXTURE = Path(__file__).resolve().parent / "fixture_reconstructed.fasta"
QUERY = Path(__file__).resolve().parent / "fixture_query.fasta"


def ensure_sieve():
    if SIEVE_BIN.exists():
        return "ok"
    if not (SIEVE_BIN.parent / "Makefile").is_file():
        return "no-submodule"
    print("[integration] building sieve (make serial)...")
    subprocess.run(["make", "serial"], cwd=SIEVE_BIN.parent, capture_output=True, text=True)
    return "ok" if SIEVE_BIN.exists() else "no-compiler"


class TestReconstructSieveSeam(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        state = ensure_sieve()
        if state == "no-submodule":
            raise unittest.SkipTest(
                f"extern/kmer-diversity-sieve is empty (run `git submodule update "
                f"--init {SIEVE_BIN.parent.relative_to(REPO)}`)"
            )
        if state == "no-compiler":
            raise unittest.SkipTest("sieve binary failed to build (no C compiler?)")
        if not FIXTURE.is_file() or not QUERY.is_file():
            raise unittest.SkipTest("fixtures missing")
        try:
            from disassembler.cmd import pipe_reconstruct  # noqa: F401
        except ImportError as e:  # pragma: no cover
            raise unittest.SkipTest(f"disassembler package not importable: {e}")

    def setUp(self):
        # Import inside setUp so the SkipTest above triggers before we touch the module.
        from disassembler.cmd import msa as msa_cmd
        from disassembler.cmd import pipe_reconstruct
        self.msa_cmd = msa_cmd
        self.pipe_reconstruct = pipe_reconstruct
        self.work = Path(tempfile.mkdtemp())

        # A fake walker: writes the fixture as reconstructed.fasta into -o <dir>, honoring -f.
        self.fake_walker = self.work / "fake_walker.sh"
        self.fake_walker.write_text(
            "#!/usr/bin/env bash\n"
            "out=''; fasta=''\n"
            'while [ $# -gt 0 ]; do case "$1" in\n'
            '  -o) out="$2"; shift 2;; -f) fasta="$2"; shift 2;; *) shift;; esac; done\n'
            f'cp "{FIXTURE}" "$fasta"\n'
        )
        self.fake_walker.chmod(0o755)

        # Preflight in pipe_reconstruct.run checks that mafft and fasttree are on PATH.
        # Satisfy that with stubs; the actual MSA/tree call is monkey-patched below and
        # never runs the real binaries.
        self.binstub = self.work / "bin"
        self.binstub.mkdir()
        for tool in ("mafft", "fasttree"):
            (self.binstub / tool).write_text("#!/usr/bin/env bash\nexit 0\n")
            (self.binstub / tool).chmod(0o755)
        self._old_path = os.environ.get("PATH", "")
        os.environ["PATH"] = f"{self.binstub}{os.pathsep}{self._old_path}"

        # Monkey-patch msa.run to record its --input and skip real work.
        self._real_msa_run = self.msa_cmd.run
        self.stage2_input = None

        def stub_msa_run(ns):
            self.stage2_input = ns.input
            Path(ns.out_aln).write_text(">stub\nACGT\n")
            Path(ns.out_tree).write_text("(stub);\n")
            return 0

        self.msa_cmd.run = stub_msa_run

    def tearDown(self):
        self.msa_cmd.run = self._real_msa_run
        os.environ["PATH"] = self._old_path
        shutil.rmtree(self.work, ignore_errors=True)

    def _make_ns(self, out_dir, **overrides):
        # Default namespace that mirrors what the argparse subparser would produce.
        ns = types.SimpleNamespace(
            accession="STUBACC",
            accession_list=None,
            query=str(QUERY),
            out_dir=str(out_dir),
            workers=None,
            fetch_workers=None,
            mafft_threads=1,
            msa_min_cov=0.70,
            emit_vcf=False,
            emit_snv_table=False,
            no_emit_alts=False,
            prune_frac=None,
            max_size_mb=None,
            no_prefer_majority=False,
            min_abund=None,
            min_vaf=None,
            walker_bin=str(self.fake_walker),
            skip_tree=False,
            sieve=True,
            sieve_threshold=0.005,
            sieve_sketch=256,
            sieve_min_len=None,
            sieve_bin=str(SIEVE_BIN),
            unitigs_cache=None,
            offline=False,
            walker_args=[],
        )
        for k, v in overrides.items():
            setattr(ns, k, v)
        return ns

    def test_on_sieve_stage2_receives_the_sieved_file(self):
        out = self.work / "out"
        ns = self._make_ns(out)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = self.pipe_reconstruct.run(ns)
        stdout = buf.getvalue()

        self.assertEqual(rc, 0, f"reconstruct exit code 0\n---\n{stdout}")
        recon = out / "reconstructed.fasta"
        sieved = out / "reconstructed.sieved.fasta"
        self.assertTrue(recon.is_file(), "reconstructed.fasta written by (fake) stage 1")
        self.assertTrue(sieved.is_file(), "reconstructed.sieved.fasta produced by sieve")
        self.assertEqual(self.stage2_input, str(sieved),
                         "stage 2 --input must be the sieved file, not the raw one")
        n = sum(1 for l in sieved.open() if l.startswith(">"))
        self.assertTrue(1 <= n < 33, f"sieve thinned 33 -> {n} (1 <= n < 33)")

    def test_sieve_sketch_reaches_the_binary_cli(self):
        out = self.work / "out_sketch"
        ns = self._make_ns(out, sieve_sketch=128)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = self.pipe_reconstruct.run(ns)
        stdout = buf.getvalue()
        self.assertEqual(rc, 0, f"reconstruct exit code 0\n---\n{stdout}")
        self.assertIn("--sketch 128", stdout,
                      "--sieve-sketch forwards --sketch 128 to the sieve binary CLI")
        self.assertIn("sketch=128", stdout,
                      "--sieve-sketch echoed in the sieve summary log line")

    def test_off_sieve_stage2_receives_the_raw_reconstruction(self):
        out = self.work / "out_off"
        ns = self._make_ns(out, sieve=False)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = self.pipe_reconstruct.run(ns)
        stdout = buf.getvalue()
        self.assertEqual(rc, 0, f"reconstruct exit code 0\n---\n{stdout}")
        self.assertEqual(self.stage2_input, str(out / "reconstructed.fasta"),
                         "off-path: stage 2 --input is the raw reconstructed.fasta")


if __name__ == "__main__":
    unittest.main()
