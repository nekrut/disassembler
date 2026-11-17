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
import random
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
            # PK-01: the composite checks the walker's version before running it.
            'if [ "$1" = --version ]; then echo "logan-walker 0.1.1"; exit 0; fi\n'
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

    def test_a_partly_errored_walker_still_builds_the_tree_and_still_fails(self):
        """Walker exit 4: some accessions errored, others reconstructed.

        The haplotypes it produced are real, so the remaining stages run on them, and the
        composite then carries the walker's code to its own exit instead of reporting
        success. Dying at stage 1 instead would discard a whole cohort's alignment and tree
        because one accession in thousands hit a network error, which contradicts the
        walker's own rule that one bad accession does not stop the run.
        """
        partial = self.work / "partial_walker.sh"
        partial.write_text(self.fake_walker.read_text() + "exit 4\n")
        partial.chmod(0o755)

        out = self.work / "out-partial"
        ns = self._make_ns(out, walker_bin=str(partial))
        buf, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(err):
            rc = self.pipe_reconstruct.run(ns)
        stdout, stderr = buf.getvalue(), err.getvalue()

        self.assertEqual(rc, 4, f"the walker's code must reach the composite's exit\n{stdout}")
        self.assertIsNotNone(self.stage2_input, f"the MSA stage must still run\n{stdout}")
        self.assertTrue((out / "reconstructed.fasta").is_file(), stdout)
        self.assertIn("skipped.tsv", stderr)

    def test_a_walker_that_failed_outright_stops_the_run_and_keeps_its_code(self):
        """Only exit 4 is carried past stage 1; the rest stop, with their own code.

        The walker's 3 ("no accession worked at all") and 1 ("nothing was reconstructed")
        say different things, and 127 says the binary was not found. `_die` used to report
        every one of them as 1, so a caller learned only that something went wrong.
        """
        for code in (1, 3, 127):
            with self.subTest(code=code):
                dead = self.work / f"dead_walker_{code}.sh"
                dead.write_text(self.fake_walker.read_text() + f"exit {code}\n")
                dead.chmod(0o755)
                self.stage2_input = None
                out = self.work / f"out-dead-{code}"
                ns = self._make_ns(out, walker_bin=str(dead))
                buf, err = io.StringIO(), io.StringIO()
                with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(err):
                    rc = self.pipe_reconstruct.run(ns)
                self.assertEqual(rc, code, buf.getvalue())
                self.assertIsNone(self.stage2_input, "stage 2 must not run")
                self.assertIn("failed (exit", err.getvalue())

    def test_a_stage_that_died_on_a_signal_is_reported_as_one_not_as_a_masked_byte(self):
        """A non-positive code means the stage did not normalise a signal death.

        Passing it on would let the shell mask it to an unrelated byte (PM-19's `-4` read as
        252), so it becomes 1, which is at least true.
        """
        self.assertEqual(self.pipe_reconstruct._stage_code(-4), 1)
        self.assertEqual(self.pipe_reconstruct._stage_code(0), 1)
        self.assertEqual(self.pipe_reconstruct._stage_code(3), 3)

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

    def test_only_the_output_is_left_behind(self):
        # The binary's out.rank*.fasta shards used to stay beside out.fasta, a second copy
        # of every promoted record.
        rc, n, log = self._run(self.fixture)
        self.assertEqual(rc, 0, log)
        self.assertGreater(n, 0)
        # out.summary.json joined these in 0.2.0; the shards and the deduped temporary
        # copy must still be gone.
        self.assertEqual(sorted(p.name for p in self.work.iterdir()),
                         ["in.fasta", "out.fasta", "out.summary.json"])

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


class TestSieveAdmissionAndHeaders(unittest.TestCase):
    """PM-01, PM-04, PM-05, PM-07, PM-13, PM-14, PM-15, PM-19, PM-20 against the real binary."""

    @classmethod
    def setUpClass(cls):
        state = ensure_sieve()
        if state != "ok":
            raise unittest.SkipTest(f"sieve binary unavailable ({state})")
        from disassembler.cmd import sieve as sieve_cmd
        cls.sieve = sieve_cmd

    def setUp(self):
        self.work = Path(tempfile.mkdtemp())
        self.rng = random.Random(11)
        self.base = "".join(self.rng.choice("ACGT") for _ in range(1500))
        # One SNV in 1,500 nt changes about 1% of the 15-mers, which is over the 0.5%
        # default threshold, so a long base is needed where near-duplicates are the point.
        self.long_base = "".join(self.rng.choice("ACGT") for _ in range(30000))

    def tearDown(self):
        shutil.rmtree(self.work, ignore_errors=True)

    def mutate(self, seq, n):
        s = list(seq)
        for pos in self.rng.sample(range(len(s)), n):
            s[pos] = self.rng.choice([b for b in "ACGT" if b != s[pos]])
        return "".join(s)

    def run_sieve(self, records, out_name="out.fasta", src_bytes=None, **kw):
        src = self.work / "in.fasta"
        out = self.work / out_name
        if src_bytes is not None:
            src.write_bytes(src_bytes)
        else:
            _write(src, records)
        kw.setdefault("min_len", 1000)
        buf, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(err):
            rc = self.sieve.sieve_run(src, out, sieve_bin=str(SIEVE_BIN), **kw)
        return rc, out, buf.getvalue(), err.getvalue()

    @staticmethod
    def headers_of(path):
        return [l[1:] for l in path.read_text().splitlines() if l.startswith(">")] \
            if path.is_file() else []

    def test_a_late_divergent_lineage_survives_a_large_unique_cohort(self):
        # PM-01. 2,100 unique haplotypes, then one 2% away. Past the binary's 2,000-record
        # warmup that record used to be rejected as an admission outlier, which means never
        # sketched, never promoted and never added to the k-mer table, so the rest of its
        # clade went with it: the lineage vanished from the MSA and tree at exit 0, and
        # whether it survived depended on where in the file it sat.
        #
        # 45 SNVs in 1,500 nt. Measured against this binary on this fixture, 30 SNVs (2%)
        # is admitted and 45 (3%) is rejected; on a cohort built with a different random
        # stream the boundary sat between 20 and 30. The admission band is
        # `mean + 5*sd + 0.5` over the running singleton counts, and record 1 contributes
        # about 1,486 novel k-mers against every later record's ~30, so the band's width is
        # a property of the FILE, not of the biology: there is no divergence a user can
        # stay under. That is a sharper statement of the defect than a single cut-off, and
        # the direction is the alarming one, since it is the genuinely novel lineages that
        # were dropped, which are the ones the sieve exists to keep.
        cohort = [(f"u{i}", self.mutate(self.base, 2)) for i in range(2100)]
        divergent = ("late_lineage", self.mutate(self.base, 45))   # 45/1500 = 3%
        rc, out, log, _err = self.run_sieve(cohort + [divergent])
        self.assertEqual(rc, 0, log)
        self.assertIn("late_lineage", self.headers_of(out),
                      "a divergent lineage after record 2,000 must still be promoted")
        m = self.sieve._REJECT_RE.search(log)
        self.assertIsNotNone(m, log)
        self.assertEqual(m.group(1), "0", "nothing may be rejected on admission")

    def test_the_same_lineage_first_was_always_kept(self):
        # The control that made this order-dependent: placed first, it survived even before.
        cohort = [(f"u{i}", self.mutate(self.base, 2)) for i in range(2100)]
        divergent = ("early_lineage", self.mutate(self.base, 45))
        rc, out, log, _err = self.run_sieve([divergent] + cohort)
        self.assertEqual(rc, 0, log)
        self.assertIn("early_lineage", self.headers_of(out))

    def test_dated_headers_do_not_partition_novelty(self):
        # PM-04: the binary folds a month parsed from `|YYYY-MM` into every LSH band key, so
        # a record was compared only with promoted records of the same month. Six copies of
        # one sequence, one SNV apart, in six different months, all survived.
        seqs = [self.mutate(self.long_base, 1) for _ in range(6)]
        dated = [(f"s{i}|2021-{i+1:02d}-15 desc", s) for i, s in enumerate(seqs)]
        plain = [(f"s{i}", s) for i, s in enumerate(seqs)]
        rc_d, out_d, log_d, _e = self.run_sieve(dated, out_name="dated.fasta")
        rc_p, out_p, _log, _e = self.run_sieve(plain, out_name="plain.fasta")
        self.assertEqual((rc_d, rc_p), (0, 0), log_d)
        self.assertEqual(len(self.headers_of(out_d)), len(self.headers_of(out_p)),
                         "the header date must not decide what is novel")
        self.assertEqual(self.headers_of(out_d), ["s0|2021-01-15 desc"],
                         "and the original header comes back intact")

    def test_a_long_header_survives_intact(self):
        # PM-13: the binary keys its emit set on the first 255 characters of the header and
        # writes that truncated copy out, so a from-lexicmap header with a long `indels=`
        # list lost its tail.
        long_hdr = "acc1 " + "indels=" + ",".join(f"{i}_{i+1}insA" for i in range(40))
        self.assertGreater(len(long_hdr), 255)
        rc, out, log, _e = self.run_sieve([(long_hdr, self.base)])
        self.assertEqual(rc, 0, log)
        self.assertEqual(self.headers_of(out), [long_hdr])

    def test_two_records_sharing_a_header_are_not_over_emitted(self):
        # PM-13: only one was promoted, but the binary emitted every record whose truncated
        # header matched, so `PROMOTED 3` came out as `4 -> 4`.
        a, b = self.mutate(self.long_base, 1), self.mutate(self.long_base, 8000)
        rc, out, log, _e = self.run_sieve([("same", a), ("same", b)])
        self.assertEqual(rc, 0, log)
        self.assertEqual(len(self.headers_of(out)), 2, "two distinct sequences, two records")
        a2 = self.mutate(a, 1)   # one SNV in 30 kb: a near-duplicate at tau=0.005
        rc, out, log, _e = self.run_sieve([("same", a), ("same", a2)],
                                          out_name="near.fasta")
        self.assertEqual(rc, 0, log)
        self.assertEqual(len(self.headers_of(out)), 1,
                         "one promoted means one emitted, whatever the headers say")

    def test_a_byte_order_mark_does_not_drop_the_first_record(self):
        # PM-07: the BOM survived decoding, so the first header failed `startswith(">")` and
        # its lines were dropped while `header is None`: the record vanished and was not
        # even counted in n_in.
        records = [(f"r{i}", self.mutate(self.base, 2)) for i in range(5)]
        body = "".join(f">{h}\n{s}\n" for h, s in records)
        rc, out, log, _e = self.run_sieve(None, src_bytes=b"\xef\xbb\xbf" + body.encode())
        self.assertEqual(rc, 0, log)
        self.assertIn("dedup: 5 -> 5", log)
        self.assertIn("r0", self.headers_of(out))

    def test_trailing_whitespace_does_not_count_as_bad_bases(self):
        # PM-07: the binary strips only \n and \r and counts every other byte as an N,
        # while the wrapper's own predictor stripped all whitespace, so records at 4% N were
        # dropped by the completeness gate with no "completeness gate drops" line and the
        # loss read as thinning.
        records = [(f"r{i}", self.mutate(self.base, 2)) for i in range(5)]
        wrapped = "".join(
            ">" + h + "\n" + "".join(f"{s[i:i+60]}   \n" for i in range(0, len(s), 60))
            for h, s in records)
        rc, out, log, _e = self.run_sieve(None, src_bytes=wrapped.encode(),
                                          max_gap_frac=0.05)
        clean_rc, clean_out, _l, _e = self.run_sieve(records, out_name="clean.fasta",
                                                     max_gap_frac=0.05)
        self.assertEqual((rc, clean_rc), (0, 0), log)
        self.assertEqual(len(self.headers_of(out)), len(self.headers_of(clean_out)))

    def test_text_before_the_first_header_is_refused(self):
        rc, _out, _log, err = self.run_sieve(
            None, src_bytes=b"; a comment\n>r0\n" + self.base.encode() + b"\n")
        self.assertEqual(rc, 1)
        self.assertIn("text before the first", err)

    def test_sibling_files_and_glob_characters_in_the_name_are_safe(self):
        # PM-05: the shard glob was `{stem}.rank*.fasta`, unescaped, in the user's own
        # directory. It deleted look-alikes before the input was read, merged any match into
        # the output, and missed the real shard when the stem held a bracket.
        for name in ("out.ranked_by_abundance.fasta", "out.rank1_manual.fasta",
                     "out.notes.fasta"):
            (self.work / name).write_text(">keep\nACGT\n")
        records = [(f"r{i}", self.mutate(self.base, 2)) for i in range(5)]
        rc, out, log, _e = self.run_sieve(records)
        self.assertEqual(rc, 0, log)
        for name in ("out.ranked_by_abundance.fasta", "out.rank1_manual.fasta",
                     "out.notes.fasta"):
            self.assertTrue((self.work / name).is_file(), f"{name} was deleted")
            self.assertEqual((self.work / name).read_text(), ">keep\nACGT\n")
        self.assertNotIn("keep", out.read_text(), "a look-alike was merged into the output")

        rc, out2, log, _e = self.run_sieve(records, out_name="cohort[v2].fasta")
        self.assertEqual(rc, 0, log)
        self.assertEqual(len(self.headers_of(out2)), 5,
                         "a glob character in the name must not hide the shards")
        # Only the two files the test planted; the run leaves no shard of its own beside
        # the output, and no `.sieve_*` staging directory either.
        self.assertEqual(sorted(f.name for f in self.work.glob("*.rank*.fasta")),
                         ["out.rank1_manual.fasta", "out.ranked_by_abundance.fasta"])
        self.assertEqual([f.name for f in self.work.glob(".sieve_*")], [])

    def test_a_run_that_promotes_nothing_leaves_the_previous_output_alone(self):
        # PM-14: the output was opened for writing before anything was counted, so a run
        # that promoted nothing truncated the earlier result, and the hint led with
        # --threshold, which cannot be the cause: an empty promoted set accepts everything.
        out = self.work / "out.fasta"
        out.write_text(">sentinel\nACGT\n")
        records = [(f"r{i}", self.mutate(self.base, 2)) for i in range(5)]
        rc, out, log, err = self.run_sieve(records, min_len=5000)
        self.assertEqual(rc, 1)
        self.assertEqual(out.read_text(), ">sentinel\nACGT\n", "the earlier result was lost")
        self.assertNotIn("--threshold", err)
        self.assertIn("--min-len", err)
        self.assertEqual([f.name for f in self.work.glob(".sieve_*")], [])

    def test_empty_records_are_not_counted_as_duplicates(self):
        # PM-15
        records = [("a", self.mutate(self.base, 2)), ("empty1", ""),
                   ("b", self.mutate(self.base, 2)), ("empty2", "")]
        rc, _out, log, _e = self.run_sieve(records)
        self.assertEqual(rc, 0, log)
        self.assertIn("dedup: 4 -> 2 records (0 exact duplicates removed, first copy kept, "
                      "2 empty records dropped)", log)

    def test_a_binary_killed_by_a_signal_is_named(self):
        # PM-19: a signal death is a negative returncode; it was printed as "exited -4" and
        # returned to sys.exit, so the shell saw 252. The sieve Makefile builds with
        # -march=native, so a binary moved to another CPU dies exactly this way.
        killer = self.work / "killer.sh"
        killer.write_text("#!/bin/sh\nkill -ILL $$\n")
        killer.chmod(0o755)
        records = [(f"r{i}", self.mutate(self.base, 2)) for i in range(3)]
        src = self.work / "in.fasta"
        _write(src, records)
        err = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
            rc = self.sieve.sieve_run(src, self.work / "o.fasta", min_len=1000,
                                      sieve_bin=str(killer))
        self.assertEqual(rc, 128 + 4)
        self.assertIn("SIGILL", err.getvalue())
        self.assertIn("-march=native", err.getvalue())

    def test_an_unwritable_output_directory_is_not_blamed_on_the_input(self):
        # PM-20: `byte_exact_dedup` both reads the input and writes the scratch copy.
        ro = self.work / "ro"
        ro.mkdir()
        src = self.work / "in.fasta"
        _write(src, [("r0", self.base)])
        ro.chmod(0o500)
        try:
            err = io.StringIO()
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
                rc = self.sieve.sieve_run(src, ro / "out.fasta", min_len=1000,
                                          sieve_bin=str(SIEVE_BIN))
            if os.geteuid() == 0:
                self.skipTest("root writes anywhere")
            self.assertEqual(rc, 1)
            self.assertIn("cannot write beside the output", err.getvalue())
        finally:
            ro.chmod(0o700)

    def test_an_input_that_would_hang_the_binary_is_refused(self):
        # PM-02: the promoted set's LSH index is a fixed 2^20-slot open-addressing table
        # that is never grown, and lookup stops only at an empty slot, so past about 65,536
        # admitted records the binary can loop forever at 100% CPU. That is upstream's code
        # and is not patched here; the wrapper refuses rather than hanging a Galaxy job
        # until its walltime. The guard is checked without building such a cohort.
        self.assertEqual(self.sieve.LSH_TABLE_LIMIT, 65536)
        rep = self.sieve.DedupReport(n_in=70000, n_out=70000, n_empty=0, headers=[],
                                     short=0, gappy=0)
        self.assertGreater(rep.n_out - rep.short - rep.gappy, self.sieve.LSH_TABLE_LIMIT)

    def test_the_admission_gate_is_disarmed_in_the_command_line(self):
        # The flag the binary is actually given, since everything above depends on it.
        self.assertEqual(self.sieve.ADMISSION_OFF_WARMUP, 2 ** 31 - 1)
        records = [(f"r{i}", self.mutate(self.base, 2)) for i in range(3)]
        _rc, _out, log, _e = self.run_sieve(records)
        self.assertIn(f"--warmup {2 ** 31 - 1}", log)


if __name__ == "__main__":
    unittest.main()
