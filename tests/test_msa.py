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

    def _msa_aligned(self, records, ref_records=None, min_cov=0.70):
        """`msa --aligned`, with MAFFT removed from PATH so a call to it would fail."""
        ref = self.work / "ref.fasta"
        inp = self.work / "in.fasta"
        aln = self.work / "aln_grid.fasta"
        aln.unlink(missing_ok=True)
        _write(ref, ref_records or [("ref", REF)])
        _write(inp, records)
        from disassembler.stage import stage_args
        ns = stage_args(
            self.msa, ref=str(ref), input=str(inp), out_aln=str(aln),
            out_tree=str(self.work / "tree.nwk"), mafft_threads=1, msa_min_cov=min_cov,
            aligned=True,
        )
        out, err = io.StringIO(), io.StringIO()
        mafft = self.work / "bin" / "mafft"
        hidden = mafft.with_suffix(".hidden")
        mafft.rename(hidden)
        try:
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                rc = self.msa.run(ns)
        finally:
            hidden.rename(mafft)
        rows = None
        if aln.is_file():
            rows, rid = {}, None
            for line in aln.read_text().splitlines():
                if line.startswith(">"):
                    rid = line[1:].strip()
                else:
                    rows[rid] = rows.get(rid, "") + line.strip()
        return rc, out.getvalue() + err.getvalue(), rows

    def test_aligned_rows_are_published_exactly_as_given(self):
        """PC-01. A lexicmap row carries every base at the reference position its HSP gave
        it, and `-` where no HSP reached, so the cohort is already an alignment. MAFFT
        discards those gaps and places the bases again, and where several placements score
        the same it slides bases across the uncovered stretch: on real HIV-1 pol rows a
        C1416G call moved to T1424C, taking a codon with it. `--aligned` takes the rows.
        """
        # Two HSPs with a hole between them, and a substitution just before the hole.
        hole = REF[:60] + "T" + "-" * 6 + REF[67:]
        rc, log, rows = self._msa_aligned([("full", REF), ("holed", hole)], min_cov=0.5)
        self.assertEqual(rc, 0, log)
        self.assertEqual(rows["holed"], hole, "the row was re-placed, not taken as given")
        self.assertEqual(rows["full"], REF)
        self.assertEqual(list(rows), ["Reference", "full", "holed"])
        self.assertIn("MAFFT is not run", log)

    def test_aligned_refuses_a_row_that_is_not_on_the_grid(self):
        """The one thing `--aligned` must not do is take a row that is not on the grid and
        publish it as an alignment: the length check exists because nothing else would
        catch it."""
        rc, log, rows = self._msa_aligned([("full", REF), ("short", REF[:-9])], min_cov=0.5)
        self.assertEqual(rc, 1, log)
        self.assertIn("are not 102 nt", log)
        self.assertIn("short (93 nt)", log)
        self.assertIsNone(rows, "nothing may be written when a row is off the grid")

    def test_aligned_keeps_the_coverage_gate_and_the_identity_check(self):
        """Everything except the placement is unchanged: a row under the coverage gate is
        still dropped before the alignment."""
        rc, log, rows = self._msa_aligned([
            ("full", REF),
            ("thin", REF[:30] + "-" * 72),  # 29% of the reference in bases
        ])
        self.assertEqual(rc, 0, log)
        self.assertEqual(list(rows), ["Reference", "full"])

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

    def test_temp_files_are_removed_when_mafft_is_interrupted(self):
        # The MAFFT inputs were removed only on the normal path, so a Ctrl-C during MAFFT left
        # .tmp_* files beside the output.
        from unittest import mock
        with mock.patch.object(self.msa.subprocess, "run", side_effect=KeyboardInterrupt):
            with self.assertRaises(KeyboardInterrupt):
                self._msa([("full", REF)])
        self.assertEqual(sorted(p.name for p in self.work.glob(".tmp_*")), [])

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
        #
        # PATH is cut down to the stub directory for this one check, so the result does not
        # depend on what the host has installed: a machine with its own /usr/bin/clustalo
        # would otherwise find that and never reach the fallback this test is about.
        from unittest import mock
        prefix = self.work / "mamba"
        codon_bin = prefix / "envs" / "codon" / "bin"
        codon_bin.mkdir(parents=True)
        (codon_bin / "clustalo").write_text("#!/bin/sh\n")
        (codon_bin / "clustalo").chmod(0o755)
        env = {"MAMBA_ROOT_PREFIX": str(prefix), "PATH": str(self.work / "bin")}
        with mock.patch.dict(os.environ, env, clear=False):
            with self.assertRaisesRegex(FileNotFoundError, "'clustalo' not found on PATH"):
                self.msa.find_bin("clustalo")
            self.assertTrue(self.msa.find_bin("mafft").endswith("/mafft"),
                            "the PATH stub is found")

    def test_real_bases(self):
        self.assertEqual(self.msa.real_bases("ACGTacgt-NnRY."), 8)


class TestMsaUnits(unittest.TestCase):
    """PM-03, PM-12, PM-22, PM-26: the pieces, without running either tool."""

    @classmethod
    def setUpClass(cls):
        from disassembler.cmd import msa as msa_cmd
        cls.msa = msa_cmd

    def test_ids_are_not_truncated(self):
        # PM-26: ids were cut to 40 characters, which is exactly where a GISAID header
        # carries the accession that tells two records apart.
        raw = ["hCoV-19/USA/CA-CDC-QDX12345678/2021|EPI_ISL_1000001|x",
               "hCoV-19/USA/CA-CDC-QDX12345678/2021|EPI_ISL_1000002|x"]
        ids, changed = self.msa.assign_ids(raw)
        self.assertTrue(changed)
        self.assertEqual(len(set(ids)), 2)
        self.assertIn("EPI_ISL_1000001", ids[0])
        self.assertIn("EPI_ISL_1000002", ids[1])

    def test_a_clean_unique_id_keeps_its_name_whatever_the_order(self):
        # PM-26: `SRR1.1` used to be renamed to `SRR1_1`, taking the name of the record that
        # was already called `SRR1_1`, which was then renamed `SRR1_1_1`.
        for raw in (["SRR1.1", "SRR1_1"], ["SRR1_1", "SRR1.1"]):
            ids, _changed = self.msa.assign_ids(raw)
            with self.subTest(raw=raw):
                self.assertEqual(ids[raw.index("SRR1_1")], "SRR1_1",
                                 "an id that is already clean and unique keeps itself")
                self.assertEqual(len(set(ids)), 2)

    def test_reference_is_reserved(self):
        ids, _c = self.msa.assign_ids(["Reference", "Reference"])
        self.assertNotIn("Reference", ids)
        self.assertEqual(len(set(ids)), 2)

    def test_many_ids_sharing_a_prefix_are_not_quadratic(self):
        # PM-22: the suffix search restarted at 1 for every record, so n ids with the same
        # prefix cost about n squared probes; 10,000 took 31.8 s in renaming alone.
        import time
        raw = ["Severe_acute_respiratory_syndrome_coronavirus_2_isolate"] * 20000
        t0 = time.perf_counter()
        ids, _c = self.msa.assign_ids(raw)
        self.assertEqual(len(set(ids)), 20000)
        self.assertLess(time.perf_counter() - t0, 5.0)

    def test_codon_qc_names_what_keeplength_did(self):
        # PM-03: the only check after MAFFT was a row-length test that --keeplength
        # guarantees, printed under the banner "Validating codon lengths...".
        ref = "ATGGCTGAACCCTTTAAA"
        clean = self.msa.codon_qc(ref, self.msa.real_bases(ref))
        self.assertFalse(clean["frameshift"])
        self.assertEqual(clean["split_codons"], [])

        # a 1 nt deletion: one gap, splitting a codon, and a frameshift by base count
        row = ref[:7] + "-" + ref[8:]
        q = self.msa.codon_qc(row, self.msa.real_bases(ref))
        self.assertEqual(q["trimmed"], 1, "one base went missing against the input")
        self.assertTrue(q["frameshift"])
        self.assertEqual(q["split_codons"], [3])
        self.assertEqual(q["internal_gap_runs"], [1])

        # an in-frame 3 nt deletion on a codon boundary: no frameshift, no split codon
        row3 = ref[:6] + "---" + ref[9:]
        q3 = self.msa.codon_qc(row3, self.msa.real_bases(ref))
        self.assertFalse(q3["frameshift"])
        self.assertEqual(q3["split_codons"], [])

        # an internal stop
        stop = "ATGTAAGAACCCTTTAAA"
        self.assertEqual(self.msa.codon_qc(stop, self.msa.real_bases(stop))["internal_stops"], [2])

    def test_an_insertion_that_keeplength_deleted_is_reported(self):
        # The row MAFFT returns is identical to the reference, which is why nothing noticed.
        ref = "ATGGCTGAACCCTTTAAA"
        q = self.msa.codon_qc(ref, self.msa.real_bases(ref) + 1)
        self.assertEqual(q["trimmed"], 1)
        self.assertTrue(q["frameshift"])

    def test_a_divergent_reverse_row_is_oriented_by_the_fallback(self):
        # PM-12: the exact 15-mers find nothing in common once a row is ~30% divergent, and
        # a tie was written out as given, so a minus-strand row aligned at ~0.43 identity
        # and stayed in the tree as noise.
        import random
        rng = random.Random(5)
        ref = "".join(rng.choice("ACGT") for _ in range(900))
        fwd, rev = self.msa.kmer_set(ref), self.msa.kmer_set(self.msa.revcomp(ref))
        fallback = (self.msa.kmer_set(ref, 1, self.msa.ORIENT_FALLBACK_K),
                    self.msa.kmer_set(self.msa.revcomp(ref), 1, self.msa.ORIENT_FALLBACK_K))
        s = list(self.msa.revcomp(ref))
        for i in range(0, len(s), 12):          # a change every 12 nt: ~8% apart
            s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
        row = "".join(s)
        self.assertEqual(self.msa.orientation(row, fwd, rev), 0, "the 15-mers tie, as before")
        self.assertEqual(self.msa.orientation(row, fwd, rev, fallback), -1,
                         "the 11-mers settle it: this row is the other strand")

    def test_a_genuinely_undecidable_row_is_still_undecided(self):
        import random
        rng = random.Random(6)
        ref = "".join(rng.choice("ACGT") for _ in range(900))
        other = "".join(rng.choice("ACGT") for _ in range(900))
        fwd, rev = self.msa.kmer_set(ref), self.msa.kmer_set(self.msa.revcomp(ref))
        fallback = (self.msa.kmer_set(ref, 1, self.msa.ORIENT_FALLBACK_K),
                    self.msa.kmer_set(self.msa.revcomp(ref), 1, self.msa.ORIENT_FALLBACK_K))
        self.assertEqual(self.msa.orientation(other, fwd, rev, fallback), 0)

    def test_a_bad_mafft_version_is_refused(self):
        # PM-17: docs/packaging.md has recorded since 0.1.0 that 7.525 segfaults on the
        # --addfragments path, and nothing ever ran `mafft --version`.
        with tempfile.TemporaryDirectory() as d:
            for text, bad in (("#!/bin/sh\necho 'v7.525 (2024)' >&2\n", True),
                              ("#!/bin/sh\necho 'v7.526 (2024)' >&2\n", False),
                              ("#!/bin/sh\necho nonsense\n", False)):
                stub = Path(d) / "mafft"
                stub.write_text(text)
                stub.chmod(0o755)
                self.msa.tool_version.cache_clear()
                if bad:
                    with self.assertRaisesRegex(RuntimeError, "7.525"):
                        self.msa.check_mafft(str(stub))
                else:
                    self.msa.check_mafft(str(stub))

    def test_publish_writes_into_a_fifo_instead_of_replacing_it(self):
        # PM-09: a rename swaps a directory entry, so a FIFO output became a regular file
        # and its reader blocked forever; /dev/null became a regular file as root.
        import threading
        with tempfile.TemporaryDirectory() as d:
            fifo = Path(d) / "out.fifo"
            os.mkfifo(fifo)
            src = Path(d) / "src"
            src.write_text("hello\n")
            got = []
            # daemon: when `publish` renames over the FIFO instead of writing into it,
            # which is the defect, nothing ever opens the write end and this thread blocks
            # forever. A non-daemon thread would then hang the whole test run rather than
            # failing this one test.
            reader = threading.Thread(target=lambda: got.append(fifo.read_text()),
                                      daemon=True)
            reader.start()
            self.msa.publish(src, fifo)
            reader.join(timeout=10)
            self.assertEqual(got, ["hello\n"])
            self.assertTrue(fifo.is_fifo(), "the FIFO must survive as a FIFO")

    def test_publish_follows_a_symlink_instead_of_replacing_it(self):
        with tempfile.TemporaryDirectory() as d:
            target = Path(d) / "target"
            target.write_text("old\n")
            link = Path(d) / "link"
            link.symlink_to(target)
            src = Path(d) / "src"
            src.write_text("new\n")
            self.msa.publish(src, link)
            self.assertTrue(link.is_symlink(), "the link must still be a link")
            self.assertEqual(target.read_text(), "new\n")

    def test_publish_renames_onto_a_plain_file(self):
        with tempfile.TemporaryDirectory() as d:
            dest = Path(d) / "dest"
            dest.write_text("old\n")
            src = Path(d) / "src"
            src.write_text("new\n")
            self.msa.publish(src, dest)
            self.assertEqual(dest.read_text(), "new\n")
            self.assertFalse(src.exists(), "a plain file is still a rename, not a copy")


if __name__ == "__main__":
    unittest.main()
