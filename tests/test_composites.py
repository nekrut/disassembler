#!/usr/bin/env python3
"""test_composites.py -- the three composites, end to end through the real CLI parser.

Each composite builds its stages' arguments with `stage_args`, from each atomic's own
parser. These tests run `disassembler.cli.main` on real arguments, so an option that a
composite fails to hand a stage (the bug that crashed every `from-kmindex` harvest) fails
here. External tools are stubs: a fake `logan-walker` that writes one reference-length
haplotype per accession, and MAFFT / FastTree stubs on PATH. Nothing touches the network.
"""
import contextlib
import io
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

REF = "ATGGCTGAATAA"


def _run(argv):
    from disassembler.cli import main
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            rc = main(argv)
        except SystemExit as e:  # argparse usage errors
            rc = e.code
    return rc, out.getvalue() + err.getvalue()


class _Stubs(unittest.TestCase):
    """A scratch dir, a reference, a fake walker, and MAFFT / FastTree stubs on PATH."""

    def setUp(self):
        self.work = Path(tempfile.mkdtemp())
        self.ref = self.work / "gene.fasta"
        self.ref.write_text(f">gene\n{REF}\n")

        # Writes the first record of the first -q as every accession's primary, and records
        # its argv.
        self.walker = self.work / "fake_walker.sh"
        self.walker.write_text(
            "#!/usr/bin/env bash\n"
            # The composites check the walker's version before running it, because
            # the published 0.1.1 container pairs this package with walker 0.1.0. A stub
            # that cannot answer --version is, as far as that check is concerned, exactly
            # the broken install it exists to catch.
            'if [ "$1" = --version ]; then echo "logan-walker 0.1.1"; exit 0; fi\n'
            'echo "$@" > "$(dirname "$0")/walker_argv.txt"\n'
            "fasta=''; q=''; acc=''; list=''\n"
            'while [ $# -gt 0 ]; do case "$1" in\n'
            '  -f) fasta="$2"; shift 2;; -q) [ -z "$q" ] && q="$2"; shift 2;;\n'
            '  -a) acc="$2"; shift 2;; -l) list="$2"; shift 2;; *) shift;; esac; done\n'
            "seq=$(gzip -dcf \"$q\" | awk '/^>/{n++; next} n==1' | tr -d '\\n')\n"
            '{ [ -n "$acc" ] && echo "$acc"; [ -n "$list" ] && gzip -dcf "$list"; } |\n'
            '  while read -r a; do [ -n "$a" ] && printf ">%s_primary\\n%s\\n" "$a" "$seq"; done'
            ' > "$fasta"\n'
        )
        bindir = self.work / "bin"
        bindir.mkdir()
        # mafft --auto --keeplength --thread N --addfragments FRAGS REF: echo both back.
        (bindir / "mafft").write_text(
            '#!/usr/bin/env bash\n'
            'ref="${@: -1}"\n'
            'while [ $# -gt 0 ]; do [ "$1" = --addfragments ] && frags="$2"; shift; done\n'
            'cat "$ref" "$frags"\n'
        )
        (bindir / "fasttree").write_text("#!/usr/bin/env bash\ncat >/dev/null\necho '(Reference);'\n")
        for exe in (self.walker, bindir / "mafft", bindir / "fasttree"):
            exe.chmod(0o755)
        self._old_path = os.environ.get("PATH", "")
        os.environ["PATH"] = f"{bindir}{os.pathsep}{self._old_path}"

    def tearDown(self):
        os.environ["PATH"] = self._old_path
        shutil.rmtree(self.work, ignore_errors=True)

    def walker_argv(self) -> str:
        return (self.work / "walker_argv.txt").read_text()


class TestStageArgs(unittest.TestCase):

    def test_unset_options_take_the_stage_default(self):
        from disassembler.cmd import harvest
        from disassembler.stage import stage_args
        ns = stage_args(harvest, history_id="H", collection_id="C")
        self.assertIsNone(ns.shards_dir)
        self.assertEqual(ns.download_workers, 6)
        self.assertIs(ns.func, harvest.run)

    def test_an_unknown_option_is_refused(self):
        from disassembler.cmd import harvest
        from disassembler.stage import stage_args
        with self.assertRaisesRegex(TypeError, "no option.*workers"):
            stage_args(harvest, workers=3)

    def test_a_missing_required_option_is_refused(self):
        from disassembler.cmd import msa
        from disassembler.stage import stage_args
        # `--out-tree` is no longer argparse-required: an alignment can be produced without
        # FastTree, under `--no-tree`. Exactly one of the two is still required, and
        # msa itself enforces that, so it is checked below rather than here.
        with self.assertRaisesRegex(TypeError, "missing required input, out_aln"):
            stage_args(msa, ref="r.fa")

    def test_msa_requires_a_tree_path_or_no_tree(self):
        # In a temp directory: `msa` writes a summary on every exit path, so a relative
        # --out-aln would drop one into whatever the working directory happens to be.
        import contextlib, io, tempfile
        from pathlib import Path as _P
        from disassembler.cmd import msa
        from disassembler.stage import stage_args
        with tempfile.TemporaryDirectory() as d:
            aln = str(_P(d) / "a.fa")
            for kw, want in (({}, "--no-tree"),
                             ({"out_tree": str(_P(d) / "t.nwk"), "no_tree": True},
                              "mutually exclusive")):
                err = io.StringIO()
                with contextlib.redirect_stdout(io.StringIO()), \
                     contextlib.redirect_stderr(err):
                    rc = msa.run(stage_args(msa, ref="r.fa", input="i.fa", out_aln=aln, **kw))
                self.assertEqual(rc, 2, err.getvalue())
                self.assertIn(want, err.getvalue())

    def test_a_required_either_or_group_is_enforced(self):
        from disassembler.cmd import pipe_reconstruct
        from disassembler.stage import stage_args
        with self.assertRaisesRegex(TypeError, "accession/accession_list"):
            stage_args(pipe_reconstruct, query=["q.fa"], out_dir="o")


class TestFromKmindex(_Stubs):

    def setUp(self):
        super().setUp()
        self.shards = self.work / "shards"
        self.shards.mkdir()
        (self.shards / "dataset_1.dat").write_text(
            json.dumps({"idx": {"gene": {"SRR1": 0.9, "SRR2": 0.4}}}))
        self.out = self.work / "out"

    def _argv(self, *extra):
        return ["from-kmindex", "--target", "gene", "-q", str(self.ref),
                "-o", str(self.out), "--walker-bin", str(self.walker), *extra]

    def test_shards_dir_runs_end_to_end(self):
        rc, log = _run(self._argv("--shards-dir", str(self.shards)))
        self.assertEqual(rc, 0, log)
        self.assertIn(f"-l {self.out}/cohort/gene.accessions.txt", self.walker_argv())
        aln = (self.out / "gene_codon_msa.fasta").read_text()
        self.assertEqual(aln.count(">"), 3, "reference plus two primaries")
        self.assertTrue((self.out / "gene_fasttree.nwk").is_file())

    def test_a_target_with_no_hits_in_this_harvest_is_refused(self):
        """The cohort is picked by file name, and harvest writes a file per target
        of this run. Re-running into the same -o with a collection in which --target has
        no hits found the previous run's file and reconstructed the previous collection's
        accessions, with exit 0 and a summary listing only the new targets.
        """
        rc, log = _run(self._argv("--shards-dir", str(self.shards), "--skip-tree"))
        self.assertEqual(rc, 0, log)
        old_cohort = (self.out / "cohort" / "gene.accessions.txt").read_text()
        self.assertIn("SRR1", old_cohort)

        # A second collection that has a different target entirely.
        other = self.work / "shards2"
        other.mkdir()
        (other / "dataset_1.dat").write_text(
            json.dumps({"idx": {"HA": {"SRR9": 0.9}}}))
        rc, log = _run(self._argv("--shards-dir", str(other), "--skip-tree"))
        self.assertEqual(rc, 1, log)
        self.assertIn("has no hits in the harvest", log)
        self.assertIn("HA", log, "the message names what the harvest did find")
        # And the stale cohort is gone, so a third run cannot pick it up either.
        self.assertFalse((self.out / "cohort" / "gene.accessions.txt").exists())
        self.assertFalse((self.out / "cohort" / "gene.scores.json").exists())

    def test_a_target_with_no_hits_is_refused_under_score_bins_too(self):
        """The same hazard through `<target>.bins.txt`: stratify scanned the whole cohort
        directory, so a stale `<target>.scores.json` was binned as well."""
        rc, log = _run(self._argv("--shards-dir", str(self.shards), "--skip-tree",
                                  "--score-bin", "0.0:1.01:10"))
        self.assertEqual(rc, 0, log)
        other = self.work / "shards2"
        other.mkdir()
        (other / "dataset_1.dat").write_text(
            json.dumps({"idx": {"HA": {"SRR9": 0.9}}}))
        rc, log = _run(self._argv("--shards-dir", str(other), "--skip-tree",
                                  "--score-bin", "0.0:1.01:10"))
        self.assertEqual(rc, 1, log)
        self.assertIn("has no hits in the harvest", log)

    def test_score_bins_feed_the_bins_list(self):
        rc, log = _run(self._argv("--shards-dir", str(self.shards),
                                  "--score-bin", "0.5:1.01:5", "--skip-tree"))
        self.assertEqual(rc, 0, log)
        self.assertIn(f"-l {self.out}/cohort/gene.bins.txt", self.walker_argv())
        self.assertEqual((self.out / "reconstructed.fasta").read_text().count(">"), 1)

    def test_skip_harvest_needs_no_galaxy_ids(self):
        rc, log = _run(self._argv("--shards-dir", str(self.shards), "--skip-tree"))
        self.assertEqual(rc, 0, log)
        rc, log = _run(self._argv("--skip-harvest", "--skip-tree"))
        self.assertEqual(rc, 0, log)
        self.assertIn("harvest skipped: reusing the earlier harvest in", log)

    def test_galaxy_ids_reach_the_harvest(self):
        # The path that crashed before reaching Galaxy. The fetch itself is stubbed.
        from disassembler.cmd import harvest
        seen = {}

        def fake_collection(coll_id, history_id, base_url, api_key, out_dir, max_workers):
            seen.update(coll=coll_id, hist=history_id, url=base_url, key=api_key,
                        workers=max_workers)
            return harvest.harvest_local(self.shards, out_dir)

        key = self.work / "gxy.txt"
        key.write_text("secret\n")
        real = harvest.harvest_collection
        harvest.harvest_collection = fake_collection
        try:
            rc, log = _run(self._argv("--history-id", "H1", "--collection-id", "C1",
                                      "--api-key-file", str(key), "--harvest-workers", "3",
                                      "--skip-tree"))
        finally:
            harvest.harvest_collection = real
        self.assertEqual(rc, 0, log)
        self.assertEqual(seen, dict(coll="C1", hist="H1", url="https://usegalaxy.org",
                                    key="secret", workers=3))

    def test_no_harvest_source_is_a_usage_error(self):
        rc, log = _run(self._argv())
        self.assertEqual(rc, 2, log)
        self.assertIn("either --shards-dir, or both --history-id and --collection-id", log)
        self.assertFalse((self.work / "walker_argv.txt").exists())

    def test_history_without_collection_is_refused(self):
        rc, log = _run(self._argv("--history-id", "H1"))
        self.assertEqual(rc, 2, log)
        self.assertIn("must be given together", log)

    def test_shards_dir_with_galaxy_ids_is_refused(self):
        rc, log = _run(self._argv("--shards-dir", str(self.shards),
                                  "--history-id", "H1", "--collection-id", "C1"))
        self.assertEqual(rc, 2, log)
        self.assertIn("mutually exclusive", log)

    def test_skip_harvest_ignores_a_source_with_a_warning(self):
        # Before 0.1.1 --history-id/--collection-id were required even with --skip-harvest,
        # so existing scripts pass placeholders; refusing them broke those scripts.
        rc, log = _run(self._argv("--shards-dir", str(self.shards), "--skip-tree"))
        self.assertEqual(rc, 0, log)
        rc, log = _run(self._argv("--skip-harvest", "--history-id", "placeholder",
                                  "--collection-id", "placeholder", "--skip-tree"))
        self.assertEqual(rc, 0, log)
        self.assertIn("--history-id, --collection-id ignored under --skip-harvest", log)

    def test_skip_harvest_with_a_hand_made_pre_0_1_1_score_file(self):
        # The SARS-CoV-2 walker run's shape: a flat {accession: score} file written by a
        # script, placeholder Galaxy IDs, and a top bin ending at 1.0.
        cohort = self.out / "cohort"
        cohort.mkdir(parents=True)
        (cohort / "gene.scores.json").write_text(
            json.dumps({"SRR1": 1.0, "SRR2": 0.95, "SRR3": 0.5}))
        rc, log = _run(self._argv("--skip-harvest", "--history-id", "placeholder",
                                  "--collection-id", "placeholder",
                                  "--score-bin", "0.90:1.00:10", "--skip-tree"))
        self.assertEqual(rc, 0, log)
        self.assertEqual((cohort / "gene.bins.txt").read_text().split(), ["SRR1", "SRR2"],
                         "the perfect score is sampled")


class TestReconstruct(_Stubs):

    def test_accession_list_to_tree(self):
        accs = self.work / "accs.txt"
        accs.write_text("SRR1\nSRR2\n")
        out = self.work / "out"
        rc, log = _run(["reconstruct", "-l", str(accs), "-q", str(self.ref),
                        "-o", str(out), "--walker-bin", str(self.walker)])
        self.assertEqual(rc, 0, log)
        self.assertEqual((out / "gene_codon_msa.fasta").read_text().count(">"), 3)

    def test_a_walker_that_writes_nothing_does_not_align_an_earlier_runs_cohort(self):
        """The reconstruction is found by a fixed name, and nothing checked that
        this run wrote it. A walker flag that makes the walker exit 0 without writing `-f`
        (`--fetch-only`, `-h`, `--version`) therefore left the composite aligning whatever
        an earlier run had put at that path, reporting it as this run's, with exit 0.
        """
        out = self.work / "out"
        out.mkdir()
        stale = out / "reconstructed.fasta"
        stale.write_text(">OLD1_primary\n" + REF + "\n>OLD2_primary\n" + REF + "\n")
        # A walker that answers --version and otherwise does nothing at all.
        quiet = self.work / "quiet_walker.sh"
        quiet.write_text(
            "#!/usr/bin/env bash\n"
            'if [ "$1" = --version ]; then echo "logan-walker 0.1.1"; exit 0; fi\n'
            "exit 0\n"
        )
        quiet.chmod(0o755)
        rc, log = _run(["reconstruct", "-a", "SRR1", "-q", str(self.ref),
                        "-o", str(out), "--walker-bin", str(quiet)])
        self.assertEqual(rc, 1, log)
        self.assertIn("wrote no", log)
        self.assertFalse((out / "gene_codon_msa.fasta").exists(),
                         "an MSA of the old cohort was built")
        self.assertFalse(stale.exists(), "the stale handoff survived the run")

    def test_walker_flags_that_skip_the_reconstruction_are_refused(self):
        """Each of these has a better form of its own, and
        `docs/codon_msa.md` invites any walker flag after `--`."""
        out = self.work / "out"
        for flag in ("--fetch-only", "-h", "--help", "-V", "--version"):
            rc, log = _run(["reconstruct", "-a", "SRR1", "-q", str(self.ref),
                            "-o", str(out), "--walker-bin", str(self.walker), "--", flag])
            self.assertEqual(rc, 1, f"{flag}: {log}")
            self.assertIn("cannot be passed through", log)
            self.assertIn(flag, log)

    def test_skip_tree_removes_an_earlier_runs_derived_outputs(self):
        # A full run, then --skip-tree into the same directory: the earlier MSA, tree and
        # sieved copy no longer match the reconstruction the walker has just rewritten.
        accs = self.work / "accs.txt"
        accs.write_text("SRR1\nSRR2\n")
        out = self.work / "out"
        argv = ["reconstruct", "-l", str(accs), "-q", str(self.ref), "-o", str(out),
                "--walker-bin", str(self.walker)]
        rc, log = _run(argv)
        self.assertEqual(rc, 0, log)
        stale = [out / "gene_codon_msa.fasta", out / "gene_fasttree.nwk",
                 out / "reconstructed.sieved.fasta"]
        stale[2].write_text(f">old\n{REF}\n")
        self.assertTrue(all(p.is_file() for p in stale))
        rc, log = _run(argv + ["--skip-tree"])
        self.assertEqual(rc, 0, log)
        self.assertTrue((out / "reconstructed.fasta").is_file())
        self.assertEqual([p.name for p in stale if p.exists()], [])

    def test_gzipped_query_and_accession_list(self):
        import gzip
        q = self.work / "gene.fa"  # gzipped, but not named .gz
        q.write_bytes(gzip.compress(f">gene\n{REF}\n".encode()))
        accs = self.work / "accs.txt"
        accs.write_bytes(gzip.compress(b"SRR1\nSRR2\n"))
        out = self.work / "out"
        rc, log = _run(["reconstruct", "-l", str(accs), "-q", str(q),
                        "-o", str(out), "--walker-bin", str(self.walker)])
        self.assertEqual(rc, 0, log)
        self.assertEqual((out / "gene_codon_msa.fasta").read_text().count(">"), 3)

    def test_a_multi_record_query_is_multi_anchor(self):
        # One -q file with two records: the notice used to need two -q flags, and the sieve's
        # floor was 95% of record 1 whatever the other references' lengths.
        import sys
        from disassembler.cmd import sieve
        panel = self.work / "panel.fasta"
        panel.write_text(f">gene\n{REF}\n>short\n{REF[:9]}\n")
        seen = {}

        def fake_sieve(src, dst, **kw):
            seen.update(kw)
            Path(dst).write_text(Path(src).read_text())
            return 0

        real = sieve.sieve_run
        sieve.sieve_run = fake_sieve
        try:
            rc, log = _run(["reconstruct", "-a", "SRR1", "-q", str(panel),
                            "-o", str(self.work / "out"), "--walker-bin", str(self.walker),
                            "--sieve", "--sieve-bin", sys.executable])
        finally:
            sieve.sieve_run = real
        self.assertEqual(rc, 0, log)
        self.assertIn("multi-anchor: 2 references. The MSA aligns against the first record of "
                      "the first -q (gene in", log)
        # The walker's own gate, not a flat 95%. The walker publishes a primary
        # from 0.85 of the reference, so thinning at 0.95 removed haplotypes both
        # neighbouring stages accept.
        self.assertEqual(seen["min_len"], int(9 * 0.85),
                         "the walker's --min-cov over the shortest reference")

    def test_the_sieve_gate_follows_a_passed_through_min_cov(self):
        """`-- --min-cov 0.6` widens what the walker publishes, and the sieve's
        completeness gate has to widen with it or it removes exactly the rows that flag
        was passed to keep."""
        import sys
        from disassembler.cmd import sieve
        panel = self.work / "panel.fasta"
        panel.write_text(f">gene\n{REF}\n")
        seen = {}

        def fake_sieve(src, dst, **kw):
            seen.update(kw)
            Path(dst).write_text(Path(src).read_text())
            return 0

        real = sieve.sieve_run
        sieve.sieve_run = fake_sieve
        try:
            rc, log = _run(["reconstruct", "-a", "SRR1", "-q", str(panel),
                            "-o", str(self.work / "out"), "--walker-bin", str(self.walker),
                            "--sieve", "--sieve-bin", sys.executable,
                            "--", "--min-cov", "0.6"])
        finally:
            sieve.sieve_run = real
        self.assertEqual(rc, 0, log)
        self.assertEqual(seen["min_len"], int(len(REF) * 0.6))
        # An explicit --sieve-min-len still wins: this is a default, not a clamp.
        seen.clear()
        sieve.sieve_run = fake_sieve
        try:
            rc, log = _run(["reconstruct", "-a", "SRR1", "-q", str(panel),
                            "-o", str(self.work / "out2"), "--walker-bin", str(self.walker),
                            "--sieve", "--sieve-bin", sys.executable,
                            "--sieve-min-len", "3"])
        finally:
            sieve.sieve_run = real
        self.assertEqual(rc, 0, log)
        self.assertEqual(seen["min_len"], 3)


class TestFromLexicmap(_Stubs):

    def test_the_cohort_reaches_the_msa_unmoved(self):
        """The tiler writes each base at the reference position its HSP gives it, so
        the cohort is already an alignment. It used to go through MAFFT anyway, which
        discards the gaps and places the bases again; where several placements score the
        same, a SNV slides across an uncovered stretch and the codon MSA and tree disagree
        with the tiler's own haplotype table, at exit 0.

        The MAFFT stub on PATH here concatenates whatever it is given, so a run that still
        called it would pass a naive row-count check. This asserts the rows themselves.
        """
        table = self.work / "hits.tsv"
        # One full-length hit and one with a hole: HSPs 1-3 and 7-12, a deletion between.
        table.write_text(
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            f"SRR10\t100.0\t100.0\t1\t12\t{REF}\t{REF}\n"
            f"SRR11\t100.0\t25.0\t1\t3\t{REF[:3]}\t{REF[:3]}\n"
            f"SRR11\t100.0\t50.0\t7\t12\t{REF[6:]}\t{REF[6:]}\n"
        )
        out = self.work / "out"
        rc, log = _run(["from-lexicmap", "-q", str(self.ref), "-f", str(table),
                        "-o", str(out), "--msa-min-cov", "0.5"])
        self.assertEqual(rc, 0, log)
        self.assertIn("MAFFT is not run", log)

        def rows(path):
            out_rows, rid = {}, None
            for line in Path(path).read_text().splitlines():
                if line.startswith(">"):
                    rid = line[1:].split()[0]
                else:
                    out_rows[rid] = out_rows.get(rid, "") + line.strip()
            return out_rows

        cohort = rows(out / "tiler.clean.accessions.fasta")
        msa = rows(out / "gene_codon_msa.fasta")
        for rid, seq in cohort.items():
            self.assertEqual(msa[rid], seq,
                             f"{rid} moved between the cohort and the MSA")

    def test_a_tiler_directory_in_the_way_is_refused_not_worked_around(self):
        """`tile` turns a prefix that is an existing directory into
        `<dir>/<ref>.*`, so the cohort lands inside `out/tiler/` while this composite
        goes on reading `out/tiler.clean.accessions.fasta`. With an earlier run's file
        at that path, the MSA and tree were built from the old cohort, with exit 0.

        The shape that gets a user here is `disassembler tile --out-dir results/tiler`
        followed by `from-lexicmap -o results`.
        """
        table = self.work / "hits.tsv"
        table.write_text(
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            f"SRR10\t100.0\t100.0\t1\t12\t{REF}\t{REF}\n"
        )
        out = self.work / "out"
        out.mkdir()
        (out / "tiler").mkdir()
        stale = out / "tiler.clean.accessions.fasta"
        stale.write_text(f">OLD1\n{REF}\n")
        rc, log = _run(["from-lexicmap", "-q", str(self.ref), "-f", str(table),
                        "-o", str(out)])
        self.assertEqual(rc, 1, log)
        self.assertIn("is a directory", log)
        self.assertFalse((out / "gene_codon_msa.fasta").exists())

    def test_a_cohort_this_run_did_not_write_is_not_aligned(self):
        """The cohort files this run reads are removed before
        the tiler runs, so one left by an earlier run cannot stand in for this run's.

        the tiler is stubbed to the exact shape that caused it: exit 0 having written no
        cohort, which is what a prefix-is-a-directory run does.
        """
        from disassembler.cmd import tiler as tiler_cmd
        table = self.work / "hits.tsv"
        table.write_text(
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            f"SRR10\t100.0\t100.0\t1\t12\t{REF}\t{REF}\n"
        )
        out = self.work / "out"
        out.mkdir()
        stale = out / "tiler.clean.accessions.fasta"
        stale.write_text(f">OLD1\n{REF}\n>OLD2\n{REF}\n")
        real = tiler_cmd.run
        tiler_cmd.run = lambda ns: 0
        try:
            rc, log = _run(["from-lexicmap", "-q", str(self.ref), "-f", str(table),
                            "-o", str(out)])
        finally:
            tiler_cmd.run = real
        self.assertEqual(rc, 1, log)
        self.assertIn("not found", log)
        msa = out / "gene_codon_msa.fasta"
        self.assertFalse(msa.exists() and "OLD1" in msa.read_text(),
                         "the old cohort was aligned as if it were this run's")

    def test_table_to_tree(self):
        table = self.work / "hits.tsv"
        table.write_text(
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            f"SRR10\t100.0\t100.0\t1\t12\t{REF}\t{REF}\n"
            "SRR11\t91.7\t100.0\t1\t12\tATGGCTGAATAA\tATGGTTGAATAA\n"
        )
        out = self.work / "out"
        rc, log = _run(["from-lexicmap", "-q", str(self.ref), "-f", str(table),
                        "-o", str(out)])
        self.assertEqual(rc, 0, log)
        self.assertEqual((out / "gene_codon_msa.fasta").read_text().count(">"), 3)
        self.assertTrue((out / "gene_fasttree.nwk").is_file())

        # The manifest ties the stage summaries together. Before 0.2.0 a composite wrote
        # nothing of its own, so a run left per-stage files with no record of the sequence,
        # the arguments, or which stage failed.
        manifest = json.loads((out / "from-lexicmap.summary.json").read_text())
        self.assertEqual(manifest["run"]["status"], "ok")
        self.assertEqual([r["stage"] for r in manifest["records"]],
                         ["tiler", "sieve", "msa"])
        by_stage = {r["stage"]: r for r in manifest["records"]}
        self.assertEqual(by_stage["sieve"]["status"], "skipped")
        # Each stage that wrote a summary names it, and the file is really there.
        for record in manifest["records"]:
            if "summary" in record:
                self.assertTrue((out / record["summary"]).is_file(), record)
        self.assertIn("summary", by_stage["tiler"])
        self.assertIn("summary", by_stage["msa"])

    def test_a_run_without_sieve_removes_an_earlier_runs_sieved_cohorts(self):
        # the tiler rewrites both cohorts every run, so sieved copies an earlier --sieve run
        # left are stale; they used to stay beside the fresh MSA.
        out = self.work / "out"
        out.mkdir()
        stale = [out / "tiler.clean.accessions.sieved.fasta",
                 out / "tiler.clean_expanded.accessions.sieved.fasta"]
        for p in stale:
            p.write_text(f">old\n{REF}\n")
        table = self.work / "hits.tsv"
        table.write_text(
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            f"SRR10\t100.0\t100.0\t1\t12\t{REF}\t{REF}\n"
        )
        rc, log = _run(["from-lexicmap", "-q", str(self.ref), "-f", str(table),
                        "-o", str(out)])
        self.assertEqual(rc, 0, log)
        self.assertTrue((out / "gene_codon_msa.fasta").is_file())
        self.assertEqual([p.name for p in stale if p.exists()], [])

    def test_allow_empty_accepts_an_empty_cohort(self):
        table = self.work / "hits.tsv"
        out = self.work / "out"
        argv = ["from-lexicmap", "-q", str(self.ref), "-f", str(table), "-o", str(out)]
        # An earlier run into the same directory made an MSA and a tree.
        table.write_text(
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            f"SRR10\t100.0\t100.0\t1\t12\t{REF}\t{REF}\n"
        )
        rc, log = _run(argv)
        self.assertEqual(rc, 0, log)
        self.assertTrue((out / "gene_codon_msa.fasta").is_file())
        # Now one partial sequence: the tiler's clean cohort (full + partial) is not empty, but
        # the full-length cohort from-lexicmap aligns is. --allow-empty used to exit 1 anyway.
        table.write_text(
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "SRR20\t100.0\t50.0\t1\t6\tATGGCT\tATGGCT\n"
        )
        rc, log = _run(argv)
        self.assertEqual(rc, 1, log)
        self.assertIn("cohort is empty", log)
        rc, log = _run(argv + ["--allow-empty"])
        self.assertEqual(rc, 0, log)
        # The earlier run's MSA and tree no longer describe this cohort, and must not sit
        # beside an exit 0 as if this run had made them.
        self.assertFalse((out / "gene_codon_msa.fasta").exists())
        self.assertFalse((out / "gene_fasttree.nwk").exists())

    def test_a_queued_dataset_stops_before_a_stale_cohort(self):
        # the tiler used to return 0 for a queued dataset having written nothing, and
        # from-lexicmap then aligned whatever cohort an earlier run had left.
        from disassembler.cmd import tiler
        out = self.work / "out"
        out.mkdir()
        (out / "tiler.clean.accessions.fasta").write_text(f">stale\n{REF}\n")
        real = tiler.resolve_dataset_state
        tiler.resolve_dataset_state = lambda *a, **k: "queued"
        try:
            rc, log = _run(["from-lexicmap", "-q", str(self.ref), "-d", "f2db41e1fa331b3e",
                            "-o", str(out)])
        finally:
            tiler.resolve_dataset_state = real
        self.assertEqual(rc, tiler.EX_TEMPFAIL, log)
        self.assertIn("is 'queued', not finished yet", log)
        self.assertFalse((out / "gene_codon_msa.fasta").exists())


if __name__ == "__main__":
    unittest.main()
