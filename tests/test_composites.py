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

        # Writes the reference as every accession's primary and records its argv.
        self.walker = self.work / "fake_walker.sh"
        self.walker.write_text(
            "#!/usr/bin/env bash\n"
            'echo "$@" > "$(dirname "$0")/walker_argv.txt"\n'
            "fasta=''; q=''; acc=''; list=''\n"
            'while [ $# -gt 0 ]; do case "$1" in\n'
            '  -f) fasta="$2"; shift 2;; -q) q="$2"; shift 2;;\n'
            '  -a) acc="$2"; shift 2;; -l) list="$2"; shift 2;; *) shift;; esac; done\n'
            'seq=$(grep -v ">" "$q" | tr -d "\\n")\n'
            '{ [ -n "$acc" ] && echo "$acc"; [ -n "$list" ] && cat "$list"; } |\n'
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
        self.assertIn("harvest skipped: reusing 1 *.scores.json", log)

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

    def test_skip_harvest_with_a_source_is_refused(self):
        rc, log = _run(self._argv("--skip-harvest", "--shards-dir", str(self.shards)))
        self.assertEqual(rc, 2, log)
        self.assertIn("--shards-dir cannot be combined with --skip-harvest", log)


class TestReconstruct(_Stubs):

    def test_accession_list_to_tree(self):
        accs = self.work / "accs.txt"
        accs.write_text("SRR1\nSRR2\n")
        out = self.work / "out"
        rc, log = _run(["reconstruct", "-l", str(accs), "-q", str(self.ref),
                        "-o", str(out), "--walker-bin", str(self.walker)])
        self.assertEqual(rc, 0, log)
        self.assertEqual((out / "gene_codon_msa.fasta").read_text().count(">"), 3)


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


if __name__ == "__main__":
    unittest.main()
