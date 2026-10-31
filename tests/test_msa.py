#!/usr/bin/env python3
"""test_msa.py -- the `msa` coverage gate counts bases, not characters.

MAFFT and FastTree are replaced by stubs on PATH: the MAFFT stub echoes the reference and
the fragments it was given (every row here is already reference-length), so the alignment
lists exactly the rows that passed the gate.
"""
import contextlib
import io
import os
import shutil
import tempfile
import types
import unittest
from pathlib import Path

REF = "ATG" + "GCA" * 32 + "TAA"  # 102 nt


def _write(path: Path, records):
    with path.open("w") as fh:
        for h, s in records:
            fh.write(f">{h}\n{s}\n")


class TestMsaCoverageGate(unittest.TestCase):

    def setUp(self):
        from disassembler.cmd import msa as msa_cmd
        self.msa = msa_cmd
        self.work = Path(tempfile.mkdtemp())
        bindir = self.work / "bin"
        bindir.mkdir()
        # mafft --auto --keeplength --thread N --addfragments FRAGS REF
        (bindir / "mafft").write_text(
            '#!/usr/bin/env bash\n'
            'ref="${@: -1}"\n'
            'while [ $# -gt 0 ]; do [ "$1" = --addfragments ] && frags="$2"; shift; done\n'
            'cat "$ref" "$frags"\n'
        )
        (bindir / "fasttree").write_text("#!/usr/bin/env bash\ncat >/dev/null\necho '(Reference);'\n")
        for tool in ("mafft", "fasttree"):
            (bindir / tool).chmod(0o755)
        self._old_path = os.environ.get("PATH", "")
        os.environ["PATH"] = f"{bindir}{os.pathsep}{self._old_path}"

    def tearDown(self):
        os.environ["PATH"] = self._old_path
        shutil.rmtree(self.work, ignore_errors=True)

    def _aligned_ids(self, records, min_cov=0.70):
        ref = self.work / "ref.fasta"
        inp = self.work / "in.fasta"
        _write(ref, [("ref", REF)])
        _write(inp, records)
        ns = types.SimpleNamespace(
            ref=str(ref), input=str(inp), out_aln=str(self.work / "aln.fasta"),
            out_tree=str(self.work / "tree.nwk"), mafft_threads=1, msa_min_cov=min_cov,
        )
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            rc = self.msa.run(ns)
        self.assertEqual(rc, 0)
        return [l[1:].strip() for l in (self.work / "aln.fasta").read_text().splitlines()
                if l.startswith(">")]

    def test_gaps_and_n_do_not_count_toward_coverage(self):
        ids = self._aligned_ids([
            ("full", REF),
            ("gap20", REF[:81] + "-" * 21),  # 79% bases
            ("gap40", REF[:60] + "-" * 42),  # 59% bases: reference-length, but under 0.70
            ("n40", "N" * 42 + REF[42:]),
            ("lower", REF.lower()),
        ])
        self.assertEqual(ids, ["Reference", "full", "gap20", "lower"])

    def test_gate_of_zero_keeps_every_row(self):
        ids = self._aligned_ids([("full", REF), ("gap40", REF[:60] + "-" * 42)], min_cov=0.0)
        self.assertEqual(ids, ["Reference", "full", "gap40"])

    def test_real_bases(self):
        self.assertEqual(self.msa.real_bases("ACGTacgt-NnRY."), 8)


if __name__ == "__main__":
    unittest.main()
