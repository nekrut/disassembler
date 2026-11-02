#!/usr/bin/env python3
"""test_msa_real_tools.py -- `disassembler msa` end to end against real MAFFT and FastTree.

`test_msa.py` replaces both tools with stubs, which is what makes its coverage-gate and
strand-orientation checks deterministic, but it means nothing in the suite ever runs the real
MAFFT command line. A stub that echoes its input cannot catch a flag this version of MAFFT
stopped honouring, and `--keeplength` is exactly the kind of flag whose silent failure would
only show up as mis-sized rows far downstream.

Skipped when MAFFT or FastTree is missing, so a checkout without them still passes. The CI job
that installs them treats a skip here as a failure, the same way the data-backed walker tests
are handled.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
QUERY = HERE / "fixture_query.fasta"
INPUT = HERE / "fixture_reconstructed.fasta"


def have_tools():
    mafft = shutil.which("mafft")
    tree = shutil.which("fasttree") or shutil.which("FastTree")
    return mafft and tree


def read_fasta(path):
    recs, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs[name] = []
        elif name:
            recs[name].append(line.strip())
    return {k: "".join(v) for k, v in recs.items()}


@unittest.skipUnless(have_tools(), "SKIP: real mafft and fasttree are not both on PATH")
class TestMsaWithRealTools(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.work = Path(tempfile.mkdtemp())
        cls.aln = cls.work / "aln.fasta"
        cls.tree = cls.work / "tree.nwk"
        # src/ goes on PYTHONPATH as an absolute path, so this works both against an installed
        # package (CI) and against a bare checkout (the local `PYTHONPATH=../src` habit).
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join(
            [str(REPO_ROOT / "src")] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
        cls.proc = subprocess.run(
            [sys.executable, "-m", "disassembler", "msa",
             "--ref", str(QUERY), "--input", str(INPUT),
             "--out-aln", str(cls.aln), "--out-tree", str(cls.tree),
             "--mafft-threads", "2"],
            capture_output=True, text=True, cwd=REPO_ROOT, env=env,
        )

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def test_it_exits_zero(self):
        self.assertEqual(self.proc.returncode, 0,
                         f"stdout:\n{self.proc.stdout}\nstderr:\n{self.proc.stderr}")

    def test_every_row_is_reference_length(self):
        # This is what --keeplength buys, and the only check that would notice MAFFT quietly
        # ignoring it: a row longer than the reference means insertions were kept, which breaks
        # the codon coordinates every downstream count relies on.
        ref_len = len(next(iter(read_fasta(QUERY).values())))
        rows = read_fasta(self.aln)
        self.assertEqual(ref_len % 3, 0, "fixture reference is not a whole number of codons")
        widths = sorted({len(s) for s in rows.values()})
        self.assertEqual(widths, [ref_len], f"row widths {widths}, reference is {ref_len}")

    def test_the_reference_and_every_input_row_are_present(self):
        rows = read_fasta(self.aln)
        self.assertEqual(len(rows), len(read_fasta(INPUT)) + 1,
                         "expected one row per input sequence plus the reference")

    def test_the_tree_has_a_tip_for_every_aligned_row(self):
        newick = self.tree.read_text().strip()
        self.assertTrue(newick.endswith(";"), f"not a newick tree: {newick[:60]!r}")
        # A tip label is one that follows `(` or `,`; a label after `)` is FastTree's support
        # value for that node. Comparing the names themselves, not a count, also catches a
        # mangled label.
        tips = set(re.findall(r"[(,]([^(),:;]+):", newick))
        self.assertEqual(tips, set(read_fasta(self.aln)),
                         "tree tips do not match the alignment's rows")


if __name__ == "__main__":
    unittest.main()
