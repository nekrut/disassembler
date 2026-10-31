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
import unittest
from pathlib import Path

REF = "ATG" + "GCA" * 32 + "TAA"  # 102 nt
# A reference with no internal repeats, for the orientation and identity checks.
GENE = ("ATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTAT"
        "GATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAAC"
        "ATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTAA")


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

    def _msa(self, records, ref_records=None, min_cov=0.70):
        """Run msa; return (exit code, log, {id: aligned row} or None if not written)."""
        ref = self.work / "ref.fasta"
        inp = self.work / "in.fasta"
        aln = self.work / "aln.fasta"
        aln.unlink(missing_ok=True)
        _write(ref, ref_records or [("ref", REF)])
        _write(inp, records)
        from disassembler.stage import stage_args
        ns = stage_args(
            self.msa, ref=str(ref), input=str(inp), out_aln=str(aln),
            out_tree=str(self.work / "tree.nwk"), mafft_threads=1, msa_min_cov=min_cov,
        )
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = self.msa.run(ns)
        rows = None
        if aln.is_file():
            rows, rid = {}, None
            for line in aln.read_text().splitlines():
                if line.startswith(">"):
                    rid = line[1:].strip()
                else:
                    rows[rid] = rows.get(rid, "") + line.strip()
        return rc, out.getvalue() + err.getvalue(), rows

    def _aligned_ids(self, records, min_cov=0.70):
        rc, log, rows = self._msa(records, min_cov=min_cov)
        self.assertEqual(rc, 0, log)
        return list(rows)

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

    def test_rows_on_the_other_strand_are_reverse_complemented(self):
        rc, log, rows = self._msa([("fwd", GENE), ("rev", self.msa.revcomp(GENE))],
                                  ref_records=[("gene", GENE)])
        self.assertEqual(rc, 0, log)
        self.assertEqual(rows["rev"], GENE)
        self.assertIn("1 of 2 sequences were on the reference's other strand", log)

    def test_a_reversed_reference_flips_every_row_and_says_so(self):
        rc, log, rows = self._msa([("a", GENE), ("b", GENE)],
                                  ref_records=[("gene_rc", self.msa.revcomp(GENE))])
        self.assertEqual(rc, 0, log)
        self.assertEqual(rows["a"], self.msa.revcomp(GENE))
        self.assertIn("every sequence was on the other strand", log)

    def test_an_alignment_that_is_mostly_noise_is_not_written(self):
        # The stub MAFFT passes rows through as given, so an unrelated row stays unrelated.
        import random
        rng = random.Random(1)
        noise = "".join(rng.choice("ACGT") for _ in GENE)
        rc, log, rows = self._msa([("x", noise), ("y", noise[::-1])], ref_records=[("gene", GENE)])
        self.assertEqual(rc, 1, log)
        self.assertIsNone(rows, "a failed identity check must not write the alignment")
        self.assertIn("under --min-identity 0.60", log)

    def test_a_multi_record_reference_uses_the_first_and_says_so(self):
        rc, log, rows = self._msa([("a", GENE)],
                                  ref_records=[("gene", GENE), ("other", REF)])
        self.assertEqual(rc, 0, log)
        self.assertIn("2 records; aligning against the first (gene", log)

    def test_gzipped_reference_and_input(self):
        import gzip
        ref = self.work / "ref.dat"  # not named .gz: gzip is decided by content
        inp = self.work / "in.dat"
        ref.write_bytes(gzip.compress(f">gene\n{GENE}\n".encode()))
        inp.write_bytes(gzip.compress(f">a\n{GENE}\n>b\n{GENE}\n".encode()))
        from disassembler.stage import stage_args
        ns = stage_args(self.msa, ref=str(ref), input=str(inp),
                        out_aln=str(self.work / "aln.fasta"),
                        out_tree=str(self.work / "tree.nwk"), mafft_threads=1)
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(self.msa.run(ns), 0)
        self.assertEqual((self.work / "aln.fasta").read_text().count(">"), 3)

    def test_orientation(self):
        fwd = self.msa.kmer_set(GENE)
        rev = self.msa.kmer_set(self.msa.revcomp(GENE))
        self.assertEqual(self.msa.orientation(GENE, fwd, rev), 1)
        self.assertEqual(self.msa.orientation(self.msa.revcomp(GENE), fwd, rev), -1)
        self.assertEqual(self.msa.orientation("ACGT" * 3, fwd, rev), 0)
        self.assertEqual(self.msa.revcomp("ACGTRYN-a"), "t-NRYACGT")

    def test_identity_ignores_gaps_and_ambiguity(self):
        self.assertEqual(self.msa.identity_to_reference("ACGTACGT", "ACGA--NT"), 4 / 5)

    def test_tools_are_found_on_path_only(self):
        # A `codon` env under $MAMBA_ROOT_PREFIX used to be a silent fallback; it belonged to
        # one old setup and could pick up a MAFFT environment.yml does not pin.
        prefix = self.work / "mamba"
        codon_bin = prefix / "envs" / "codon" / "bin"
        codon_bin.mkdir(parents=True)
        (codon_bin / "clustalo").write_text("#!/bin/sh\n")
        (codon_bin / "clustalo").chmod(0o755)
        old = os.environ.get("MAMBA_ROOT_PREFIX")
        os.environ["MAMBA_ROOT_PREFIX"] = str(prefix)
        try:
            with self.assertRaisesRegex(FileNotFoundError, "'clustalo' not found on PATH"):
                self.msa.find_bin("clustalo")
        finally:
            if old is None:
                del os.environ["MAMBA_ROOT_PREFIX"]
            else:
                os.environ["MAMBA_ROOT_PREFIX"] = old
        self.assertTrue(self.msa.find_bin("mafft").endswith("/mafft"), "the PATH stub is found")

    def test_real_bases(self):
        self.assertEqual(self.msa.real_bases("ACGTacgt-NnRY."), 8)


if __name__ == "__main__":
    unittest.main()
