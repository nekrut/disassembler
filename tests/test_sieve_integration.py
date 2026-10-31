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
        # Parsed by the real `reconstruct` parser, so every option `run` reads exists with
        # its real default. A hand-built namespace went stale each time an option changed
        # (-q became repeatable, --msa-reference was added) and failed before the seam.
        from disassembler.cli import build_parser
        ns = build_parser().parse_args([
            "reconstruct",
            "-a", "STUBACC",
            "-q", str(QUERY),
            "-o", str(out_dir),
            "--walker-bin", str(self.fake_walker),
            "--sieve",
            "--sieve-bin", str(SIEVE_BIN),
            "--mafft-threads", "1",
        ])
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


def _write(path: Path, records):
    with path.open("w") as fh:
        for h, s in records:
            fh.write(f">{h}\n{s}\n")


def _count(path: Path) -> int:
    return sum(1 for l in path.open() if l.startswith(">")) if path.is_file() else 0


def _gapped(seq: str, frac: float) -> str:
    """`seq` with its last `frac` of positions turned to '-' (a reference-grid partial)."""
    n = int(len(seq) * frac)
    return seq[:len(seq) - n] + "-" * n


class TestSieveRun(unittest.TestCase):
    """`sieve_run` against the real binary: deduplication and the completeness gate."""

    @classmethod
    def setUpClass(cls):
        state = ensure_sieve()
        if state != "ok":
            raise unittest.SkipTest(f"sieve binary unavailable ({state})")
        from disassembler.cmd import sieve as sieve_cmd
        cls.sieve = sieve_cmd
        cls.fixture = []
        header = None
        for line in FIXTURE.read_text().splitlines():
            if line.startswith(">"):
                header = line[1:].split()[0]
            else:
                cls.fixture.append((header, line.strip()))

    def setUp(self):
        self.work = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.work, ignore_errors=True)

    def _run(self, records, **kw):
        src = self.work / "in.fasta"
        out = self.work / "out.fasta"
        _write(src, records)
        kw.setdefault("min_len", 100)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
            rc = self.sieve.sieve_run(src, out, sieve_bin=str(SIEVE_BIN), **kw)
        return rc, _count(out), buf.getvalue()

    def test_duplicates_past_the_warmup_do_not_reject_new_lineages(self):
        # 2100 copies of four sequences, then 20 that each differ from them at ~30 sites.
        # Every copy adds a zero to the binary's singleton-count statistics, so once its
        # 2000-record warmup is over each new lineage looks like an outlier and is
        # rejected as an error. Deduplicated, the cohort is 24 records and all survive.
        import random
        rng = random.Random(7)
        base = "".join(rng.choice("ACGT") for _ in range(900))

        def mutate(n):
            s = list(base)
            for p in rng.sample(range(len(s)), n):
                s[p] = rng.choice([b for b in "ACGT" if b != s[p]])
            return "".join(s)

        common = [mutate(3) for _ in range(4)]
        records = [(f"dup{i}", common[i % 4]) for i in range(2100)]
        records += [(f"new{j}", mutate(30)) for j in range(20)]

        raw = self.work / "raw.fasta"
        _write(raw, records)
        subprocess.run([str(SIEVE_BIN), "--input", str(raw), "--min-len", "850",
                        "--max-n-frac", "0.05", "--threshold", "0.005", "--warmup", "2000",
                        "--sketch", "256", "--emit-promoted", str(self.work / "raw")],
                       capture_output=True, check=True)
        self.assertLess(_count(self.work / "raw.rank0.fasta"), 24,
                        "premise: on the raw cohort the binary rejects new lineages")

        rc, n, log = self._run(records, min_len=850)
        self.assertEqual(rc, 0, log)
        self.assertIn("dedup: 2120 -> 24 records", log)
        self.assertEqual(n, 24, "all four common and all 20 new sequences are promoted")

    def test_a_gzipped_input_gives_the_same_output(self):
        import gzip
        rc, n_plain, log = self._run(self.fixture)
        self.assertEqual(rc, 0, log)
        plain = (self.work / "out.fasta").read_text()
        src = self.work / "in.dat"
        src.write_bytes(gzip.compress("".join(f">{h}\n{s}\n" for h, s in self.fixture).encode()))
        out = self.work / "out_gz.fasta"
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
            rc = self.sieve.sieve_run(src, out, min_len=100, sieve_bin=str(SIEVE_BIN))
        self.assertEqual(rc, 0, buf.getvalue())
        self.assertEqual(out.read_text(), plain)

    def test_duplicates_keep_the_first_copy(self):
        seq = self.fixture[0][1]
        rc, n, log = self._run([("first", seq), ("second", seq)])
        self.assertEqual((rc, n), (0, 1), log)
        self.assertIn(">first", (self.work / "out.fasta").read_text())

    def test_gapped_rows_fail_the_default_gate_and_pass_a_raised_one(self):
        rows = [(h, _gapped(s, 0.20)) for h, s in self.fixture]
        rc, n, log = self._run(rows)
        self.assertEqual((rc, n), (1, 0), log)
        self.assertIn("completeness gate drops 33 of 33: 0 shorter than 100 nt (--min-len), "
                      "33 with more than 5% gaps", log)

        rc, n, log = self._run(rows, max_gap_frac=0.25)
        self.assertEqual(rc, 0, log)
        self.assertGreaterEqual(n, 1)
        self.assertIn("--max-n-frac 0.25", log)
        self.assertNotIn("completeness gate drops", log)

    def test_completeness_failures_agrees_with_the_binary(self):
        s = self.fixture[0][1]
        records = [
            ("full", s),
            ("short", s[:1000]),
            ("gap10", _gapped(self.fixture[1][1], 0.10)),
            ("gap04", _gapped(self.fixture[2][1], 0.04)),
            ("n10", "N" * 274 + self.fixture[3][1][274:]),
            ("lower", self.fixture[4][1].lower()),
        ]
        src = self.work / "mixed.fasta"
        _write(src, records)
        short, gappy = self.sieve.completeness_failures(src, 2000, 0.05)
        self.assertEqual((short, gappy), (1, 2))
        out = subprocess.run([str(SIEVE_BIN), "--input", str(src), "--min-len", "2000",
                              "--max-n-frac", "0.05", "--emit-promoted",
                              str(self.work / "mixed")],
                             capture_output=True, text=True, check=True).stdout
        self.assertRegex(out, r"LOW-QUAL \(partial\) : 3 ")


class TestAutoMinLen(unittest.TestCase):
    """The auto --min-len is 95% of the shortest reference record, not of record 1."""

    def test_sized_on_the_shortest_record(self):
        from disassembler.cmd import sieve as sieve_cmd
        with tempfile.TemporaryDirectory() as d:
            ref = Path(d) / "panel.fasta"
            ref.write_text(">long\n" + "A" * 1000 + "\n>short\n" + "C" * 400 + "\n" + "C" * 400 + "\n")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                self.assertEqual(sieve_cmd.auto_min_len(ref), 760)
            self.assertIn("2 records; sizing --min-len on the shortest (short, 800 nt)",
                          buf.getvalue())
            ref.write_text(">only\n" + "A" * 1000 + "\n")
            self.assertEqual(sieve_cmd.auto_min_len(ref), 950)
            import gzip
            ref.write_bytes(gzip.compress(b">only\n" + b"A" * 1000 + b"\n"))
            self.assertEqual(sieve_cmd.auto_min_len(ref), 950)
            ref.write_text(">empty\n")
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertIsNone(sieve_cmd.auto_min_len(ref))


class TestFromLexicmapGates(unittest.TestCase):
    """from-lexicmap derives the sieve's gap limit and the MSA gate from its cohort gate.

    Every stage is stubbed: lexicmap writes a two-row cohort, and the sieve and MSA record
    the values they were given.
    """

    def setUp(self):
        from disassembler.cmd import lexicmap as lexicmap_cmd
        from disassembler.cmd import msa as msa_cmd
        from disassembler.cmd import pipe_from_lexicmap
        from disassembler.cmd import sieve as sieve_cmd
        self.mods = (lexicmap_cmd, msa_cmd, sieve_cmd)
        self.real = (lexicmap_cmd.run, msa_cmd.run, sieve_cmd.sieve_run)
        self.pipe = pipe_from_lexicmap
        self.work = Path(tempfile.mkdtemp())
        self.seen = {}

        def fake_lexicmap(ns):
            for name in ("clean", "clean_expanded"):
                _write(Path(f"{ns.out_prefix}.{name}.accessions.fasta"),
                       [("a", "ACGTAC"), ("b", "ACG---")])
            return 0

        def fake_sieve(src, dst, **kw):
            self.seen["sieve_input"] = Path(src).name
            self.seen["max_gap_frac"] = kw["max_gap_frac"]
            Path(dst).write_text(Path(src).read_text())
            return 0

        def fake_msa(ns):
            self.seen["msa_min_cov"] = ns.msa_min_cov
            return 0

        lexicmap_cmd.run, msa_cmd.run, sieve_cmd.sieve_run = fake_lexicmap, fake_msa, fake_sieve

    def tearDown(self):
        lexicmap_cmd, msa_cmd, sieve_cmd = self.mods
        lexicmap_cmd.run, msa_cmd.run, sieve_cmd.sieve_run = self.real
        shutil.rmtree(self.work, ignore_errors=True)

    def _run(self, *extra):
        from disassembler.cli import build_parser
        ns = build_parser().parse_args([
            "from-lexicmap", "-q", str(QUERY), "-f", str(QUERY),
            "-o", str(self.work / "out"), "--sieve", "--sieve-bin", sys.executable, *extra,
        ])
        with contextlib.redirect_stdout(io.StringIO()):
            rc = self.pipe.run(ns)
        self.assertEqual(rc, 0)
        return self.seen

    def test_default_cohort(self):
        seen = self._run()
        self.assertEqual(seen["sieve_input"], "lexicmap.clean.accessions.fasta")
        self.assertEqual((seen["max_gap_frac"], seen["msa_min_cov"]), (0.25, 0.70))

    def test_full_length_cohort(self):
        seen = self._run("--require-full-length")
        self.assertEqual((seen["max_gap_frac"], seen["msa_min_cov"]), (0.05, 0.70))

    def test_partials_cohort(self):
        seen = self._run("--include-partials")
        self.assertEqual(seen["sieve_input"], "lexicmap.clean_expanded.accessions.fasta")
        self.assertEqual((seen["max_gap_frac"], seen["msa_min_cov"]), (0.55, 0.45))

    def test_lower_min_cov_moves_both_gates(self):
        seen = self._run("--min-cov", "0.6", "--min-cov-partial", "0.5")
        self.assertEqual((seen["max_gap_frac"], seen["msa_min_cov"]), (0.45, 0.55))

    def test_explicit_values_win(self):
        seen = self._run("--sieve-max-gap-frac", "0.1", "--msa-min-cov", "0.9")
        self.assertEqual((seen["max_gap_frac"], seen["msa_min_cov"]), (0.1, 0.9))


if __name__ == "__main__":
    unittest.main()
