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
        with self.assertRaisesRegex(TypeError, "missing required input, out_aln, out_tree"):
            stage_args(msa, ref="r.fa")

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
        self.assertEqual(seen["min_len"], int(9 * 0.95), "95% of the shortest reference")


class TestFromLexicmap(_Stubs):

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

    def test_allow_empty_accepts_an_empty_cohort(self):
        # One partial sequence: lexicmap's clean cohort (full + partial) is not empty, but the
        # full-length cohort from-lexicmap aligns is. --allow-empty used to exit 1 here anyway.
        table = self.work / "hits.tsv"
        table.write_text(
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "SRR20\t100.0\t50.0\t1\t6\tATGGCT\tATGGCT\n"
        )
        out = self.work / "out"
        argv = ["from-lexicmap", "-q", str(self.ref), "-f", str(table), "-o", str(out)]
        rc, log = _run(argv)
        self.assertEqual(rc, 1, log)
        self.assertIn("cohort is empty", log)
        rc, log = _run(argv + ["--allow-empty"])
        self.assertEqual(rc, 0, log)
        self.assertFalse((out / "gene_codon_msa.fasta").exists())

    def test_a_queued_dataset_stops_before_a_stale_cohort(self):
        # lexicmap used to return 0 for a queued dataset having written nothing, and
        # from-lexicmap then aligned whatever cohort an earlier run had left.
        from disassembler.cmd import lexicmap
        out = self.work / "out"
        out.mkdir()
        (out / "lexicmap.clean.accessions.fasta").write_text(f">stale\n{REF}\n")
        real = lexicmap.resolve_dataset_state
        lexicmap.resolve_dataset_state = lambda *a, **k: "queued"
        try:
            rc, log = _run(["from-lexicmap", "-q", str(self.ref), "-d", "f2db41e1fa331b3e",
                            "-o", str(out)])
        finally:
            lexicmap.resolve_dataset_state = real
        self.assertEqual(rc, lexicmap.EX_TEMPFAIL, log)
        self.assertIn("is 'queued', not finished yet", log)
        self.assertFalse((out / "gene_codon_msa.fasta").exists())


if __name__ == "__main__":
    unittest.main()
