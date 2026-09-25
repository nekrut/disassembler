#!/usr/bin/env python3
"""test_lexicmap.py: Unit and integration tests for the lexicmap subcommand.

Tests:
  - Multi-HSP coordinate tiling & overlap resolution
  - In-frame translation and premature stop-codon detection
  - 4-tier cohort partitioning
  - Online haplotype collapsing and mutation classification
  - Streaming from file and stdin
"""

import sys
import unittest
import tempfile
from pathlib import Path

try:
    from disassembler.cmd import lexicmap as lms
except ImportError as e:  # pragma: no cover - environment-dependent
    print(f"SKIP: disassembler.cmd.lexicmap not importable ({e})")
    sys.exit(0)


class TestLexicMapStreamer(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.dir_path = Path(self.test_dir.name)

        # 12 nt CDS: ATG GCT GAA TAA (Met Ala Glu Stop) -> 4 codons
        self.ref_fasta = self.dir_path / "ref.fasta"
        self.ref_fasta.write_text(">test_cds len=12\nATGGCTGAATAA\n")

    def tearDown(self):
        self.test_dir.cleanup()

    def test_translation_and_stops(self):
        self.assertEqual(lms.translate_nt("ATGGCTGAATAA"), "MAE*")
        self.assertEqual(lms.find_internal_stops("ATGGCTGAATAA", 12), [])

        seq_stop = "ATGTAGGAATAA"
        self.assertEqual(lms.translate_nt(seq_stop), "M*E*")
        self.assertEqual(lms.find_internal_stops(seq_stop, 12), [2])

        self.assertEqual(lms.find_internal_stops("ATGGCTGAATGA", 12), [])

    def test_multi_hsp_overlap_resolution(self):
        hsp1 = (100.0, 50.0, 1, 6, "ATGGCT", "ATGGCT")
        hsp2 = (90.0, 66.7, 5, 12, "CTGAATAA", "TTGAATAA")

        merged, cov_bases, mean_pid, single_ge_cov, single_full = lms.merge_accession_hsps(
            [hsp1, hsp2], ref_len=12, min_pident=60.0, min_coverage=0.80
        )

        self.assertEqual(merged, "ATGGCTGAATAA")
        self.assertEqual(cov_bases, 12)
        self.assertAlmostEqual(cov_bases / 12, 1.0)
        self.assertFalse(single_ge_cov)
        self.assertFalse(single_full)

    def test_end_to_end_streaming(self):
        tsv_content = (
            "query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\tcls\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tqseq\tsseq\n"
            "test_cds\t12\t2\tSRR001\tu1\t100.0\t1\t1\t50.0\t6\t100.0\t0\t1\t6\t1\t6\tATGGCT\tATGGCT\n"
            "test_cds\t12\t2\tSRR001\tu2\t100.0\t1\t2\t50.0\t6\t100.0\t0\t7\t12\t1\t6\tGAATAA\tGAATAA\n"
            "test_cds\t12\t1\tSRR002\tu3\t100.0\t1\t1\t100.0\t12\t91.7\t0\t1\t12\t1\t12\tATGGCTGAATAA\tATGGTTGAATAA\n"
            "test_cds\t12\t1\tSRR003\tu4\t100.0\t1\t1\t100.0\t12\t91.7\t0\t1\t12\t1\t12\tATGGCTGAATAA\tATGTAGGAATAA\n"
            "test_cds\t12\t1\tSRR004\tu5\t50.0\t1\t1\t50.0\t6\t100.0\t0\t1\t6\t1\t6\tATGGCT\tATGGCT\n"
            "test_cds\t12\t1\tSRR005\tu6\t33.3\t1\t1\t33.3\t4\t100.0\t0\t1\t4\t1\t4\tATGG\tATGG\n"
        )
        tsv_path = self.dir_path / "lexicmap_search.tsv"
        tsv_path.write_text(tsv_content)

        out_prefix = str(self.dir_path / "output")
        summary = lms.stream_lexicmap_to_msa(
            ref_path=str(self.ref_fasta),
            file_path=str(tsv_path),
            out_prefix=out_prefix
        )

        breakdown = summary["cohort_breakdown"]
        self.assertEqual(breakdown["clean_full_length"], 2)
        self.assertEqual(breakdown["clean_partial"], 1)
        self.assertEqual(breakdown["combined_clean_cohort"], 3)
        self.assertEqual(breakdown["flagged_premature_stops"], 1)
        self.assertEqual(breakdown["flagged_low_coverage"], 1)

        self.assertEqual(summary["recovery_metrics"]["multi_hsp_rescued_accessions"], 1)

        self.assertTrue(Path(f"{out_prefix}.clean.msa.fasta").exists())
        self.assertTrue(Path(f"{out_prefix}.clean.haplotypes.tsv").exists())
        self.assertTrue(Path(f"{out_prefix}.cohort_ledger.tsv").exists())
        self.assertTrue(Path(f"{out_prefix}.flagged.tsv").exists())

        tsv_lines = Path(f"{out_prefix}.clean.haplotypes.tsv").read_text().strip().split("\n")
        header = tsv_lines[0].split("\t")
        class_idx = header.index("mutation_class")

        classes = [line.split("\t")[class_idx] for line in tsv_lines[1:]]
        self.assertIn("WT", classes)
        self.assertIn("Missense", classes)

    def test_mean_pident_excludes_gap_positions(self):
        merged, cov_bases, mean_pid, _, _ = lms.merge_accession_hsps(
            [(90.0, 100.0, 1, 6, "ATGGCT", "AT-GCT")], ref_len=6, min_pident=60.0, min_coverage=0.80
        )
        self.assertEqual(cov_bases, 5)
        self.assertLessEqual(mean_pid, 100.0)
        self.assertAlmostEqual(mean_pid, 90.0)

    def test_ungrouped_input_is_rejected(self):
        tsv = (
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "SRR001\t100.0\t50.0\t1\t6\tATGGCT\tATGGCT\n"
            "SRR002\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
            "SRR001\t100.0\t50.0\t7\t12\tGAATAA\tGAATAA\n"
        )
        p = self.dir_path / "ungrouped.tsv"
        p.write_text(tsv)
        with self.assertRaises(ValueError):
            lms.stream_lexicmap_to_msa(
                ref_path=str(self.ref_fasta), file_path=str(p),
                out_prefix=str(self.dir_path / "ung"),
            )

    def test_empty_stream_does_not_crash(self):
        p = self.dir_path / "empty.tsv"
        p.write_text("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n")
        summary = lms.stream_lexicmap_to_msa(
            ref_path=str(self.ref_fasta), file_path=str(p),
            out_prefix=str(self.dir_path / "emp"),
        )
        self.assertEqual(summary["total_accessions_evaluated"], 0)
        self.assertEqual(summary["cohort_breakdown"]["clean_full_length"], 0)

    def test_multi_record_reference_uses_only_the_first(self):
        multi = self.dir_path / "multi.fasta"
        multi.write_text(">first_cds\nATGGCTGAATAA\n>second_cds\nATGAAACCCGGGTTTTAA\n")
        name, header, seq = lms.load_reference(str(multi))
        self.assertEqual(name, "first_cds")
        self.assertEqual(seq, "ATGGCTGAATAA", "record 2 must not be concatenated on")
        self.assertEqual(len(seq), 12)

        name, header, seq = lms.load_reference(str(self.ref_fasta))
        self.assertEqual((name, seq), ("test_cds", "ATGGCTGAATAA"))

    def test_full_length_only_requires_complete_coverage(self):
        tsv = (
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "SRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
            "SRR011\t100.0\t50.0\t1\t6\tATGGCT\tATGGCT\n"
            "SRR012\t100.0\t91.7\t1\t11\tATGGCTGAATA\tATGGCTGAATA\n"
        )
        p = self.dir_path / "fl.tsv"
        p.write_text(tsv)

        on = lms.stream_lexicmap_to_msa(
            ref_path=str(self.ref_fasta), file_path=str(p),
            out_prefix=str(self.dir_path / "fl_on"), full_length_only=True,
        )
        off = lms.stream_lexicmap_to_msa(
            ref_path=str(self.ref_fasta), file_path=str(p),
            out_prefix=str(self.dir_path / "fl_off"), full_length_only=False,
        )

        self.assertEqual(on["cohort_breakdown"]["clean_full_length"], 1)
        self.assertEqual(off["cohort_breakdown"]["clean_full_length"], 2)
        self.assertNotEqual(
            on["cohort_breakdown"]["clean_full_length"],
            off["cohort_breakdown"]["clean_full_length"],
            "fixture must discriminate, or deleting the flag's effect would leave this green",
        )


class TestMalformedInput(unittest.TestCase):
    """One bad field used to abort the stream after the cohort files were partly written."""

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.dir_path = Path(self.test_dir.name)
        self.ref = self.dir_path / "ref.fasta"
        self.ref.write_text(">test_cds\nATGGCTGAATAA\n")
        self.header = "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
        self.good = "SRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"

    def tearDown(self):
        self.test_dir.cleanup()

    def _run(self, tsv, prefix):
        p = self.dir_path / f"{prefix}.tsv"
        p.write_text(tsv)
        return lms.stream_lexicmap_to_msa(
            ref_path=str(self.ref), file_path=str(p),
            out_prefix=str(self.dir_path / prefix),
        )

    def test_a_bad_numeric_field_skips_the_row_and_finishes(self):
        tsv = self.header + self.good + "SRR011\tNA\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
        summary = self._run(tsv, "na")
        self.assertEqual(summary["cohort_breakdown"]["clean_full_length"], 1)
        self.assertTrue((self.dir_path / "na.summary.json").is_file())

    def test_an_embedded_header_row_is_skipped(self):
        tsv = self.header + self.good + self.header.strip() + "\n" + self.good.replace("SRR010", "SRR012")
        summary = self._run(tsv, "hdr")
        self.assertEqual(summary["total_accessions_evaluated"], 2)

    def test_every_row_malformed_is_an_error_not_an_empty_cohort(self):
        tsv = self.header + "SRR010\tNA\tNA\tNA\tNA\tATG\tATG\n"
        with self.assertRaises(ValueError) as cm:
            self._run(tsv, "allbad")
        self.assertIn("unparseable numeric fields", str(cm.exception))


class TestReferenceGuards(unittest.TestCase):
    """load_reference is the choke point for every coverage figure."""

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.dir_path = Path(self.test_dir.name)

    def tearDown(self):
        self.test_dir.cleanup()

    def test_an_empty_reference_is_named_not_a_zero_division(self):
        p = self.dir_path / "empty.fasta"
        p.write_text("")
        with self.assertRaises(ValueError) as cm:
            lms.load_reference(str(p))
        self.assertIn("no sequence", str(cm.exception))

    def test_a_header_only_reference_is_named(self):
        p = self.dir_path / "hdr.fasta"
        p.write_text(">just_a_header\n")
        with self.assertRaises(ValueError):
            lms.load_reference(str(p))

    def test_a_directory_is_a_file_not_found_not_an_isadirectory(self):
        with self.assertRaises(FileNotFoundError):
            lms.load_reference(str(self.dir_path))

    def test_a_real_reference_still_loads(self):
        p = self.dir_path / "ok.fasta"
        p.write_text(">cds desc here\nATGGCTGAATAA\n")
        name, header, seq = lms.load_reference(str(p))
        self.assertEqual((name, seq), ("cds", "ATGGCTGAATAA"))


class TestQueryColumn(unittest.TestCase):
    """A multi-target search emits rows for every target in one file."""

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.dir_path = Path(self.test_dir.name)
        self.ref = self.dir_path / "ref.fasta"
        self.ref.write_text(">geneA\nATGGCTGAATAA\n")

    def tearDown(self):
        self.test_dir.cleanup()

    def _run(self, tsv, prefix, **kw):
        p = self.dir_path / f"{prefix}.tsv"
        p.write_text(tsv)
        return lms.stream_lexicmap_to_msa(
            ref_path=str(self.ref), file_path=str(p),
            out_prefix=str(self.dir_path / prefix), **kw)

    def test_rows_for_another_gene_are_not_tiled_onto_this_reference(self):
        tsv = (
            "query\tsgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "geneA\tSRR010\t100.0\t50.0\t1\t6\tATGGCT\tATGGCT\n"
            "geneB\tSRR010\t100.0\t50.0\t7\t12\tTTTTTT\tTTTTTT\n"
        )
        summary = self._run(tsv, "multi", min_coverage=0.5)
        self.assertEqual(summary["total_accessions_evaluated"], 1)
        hap = (self.dir_path / "multi.clean.haplotypes.tsv").read_text()
        self.assertIn("ATGGCT------", hap, "only geneA's half should be tiled")
        self.assertNotIn("TTTTTT", hap, "geneB's bases must not reach geneA's grid")

    def test_a_file_sorted_by_query_no_longer_trips_the_grouping_check(self):
        tsv = (
            "query\tsgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "geneA\tSRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
            "geneA\tSRR011\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
            "geneB\tSRR010\t100.0\t100.0\t1\t12\tTTTTTTTTTTTT\tTTTTTTTTTTTT\n"
        )
        summary = self._run(tsv, "sorted")
        self.assertEqual(summary["total_accessions_evaluated"], 2, "both geneA accessions")

    def test_no_row_for_this_query_is_an_error_not_an_empty_cohort(self):
        tsv = (
            "query\tsgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "geneB\tSRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
        )
        with self.assertRaises(ValueError) as cm:
            self._run(tsv, "none")
        self.assertIn("no row matched query 'geneA'", str(cm.exception))

    def test_a_file_without_the_column_still_works(self):
        tsv = (
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "SRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
        )
        summary = self._run(tsv, "nocol")
        self.assertEqual(summary["cohort_breakdown"]["clean_full_length"], 1)

    def test_a_qlen_disagreement_is_refused(self):
        tsv = (
            "query\tqlen\tsgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "geneA\t999\tSRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
        )
        with self.assertRaises(ValueError) as cm:
            self._run(tsv, "qlen")
        self.assertIn("reference length disagreement", str(cm.exception))

    def test_a_matching_qlen_passes(self):
        tsv = (
            "query\tqlen\tsgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "geneA\t12\tSRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
        )
        summary = self._run(tsv, "qlenok")
        self.assertEqual(summary["cohort_breakdown"]["clean_full_length"], 1)


class TestMutationClassification(unittest.TestCase):
    """`Complex` was unreachable, and a partially covered match was labelled WT."""

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.dir_path = Path(self.test_dir.name)
        self.ref = self.dir_path / "ref.fasta"
        self.ref.write_text(">cds\nATGGCTGAATAA\n")

    def tearDown(self):
        self.test_dir.cleanup()

    def _classes(self, rows, prefix, **kw):
        header = "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
        p = self.dir_path / f"{prefix}.tsv"
        p.write_text(header + rows)
        lms.stream_lexicmap_to_msa(
            ref_path=str(self.ref), file_path=str(p),
            out_prefix=str(self.dir_path / prefix), **kw)
        out = (self.dir_path / f"{prefix}.clean.haplotypes.tsv").read_text().splitlines()
        cols = out[0].split("\t")
        i = cols.index("mutation_class")
        return [r.split("\t")[i] for r in out[1:]]

    def test_complex_is_gone_and_missense_is_reported(self):
        rows = "SRR010\t100.0\t100.0\t1\t12\tATGGGTGAATAA\tATGGGTGAATAA\n"
        self.assertEqual(self._classes(rows, "mis"), ["Missense"])

    def test_synonymous_and_nonsense_still_classify(self):
        syn = "SRR010\t100.0\t100.0\t1\t12\tATGGCTGAGTAA\tATGGCTGAGTAA\n"
        self.assertEqual(self._classes(syn, "syn"), ["Synonymous"])
        non = "SRR010\t100.0\t100.0\t1\t12\tATGTGAGAATAA\tATGTGAGAATAA\n"
        self.assertEqual(self._classes(non, "non", max_internal_stops=1), ["Nonsense"])

    def test_an_exact_match_is_WT(self):
        rows = "SRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
        self.assertEqual(self._classes(rows, "wt"), ["WT"])

    def test_a_partially_covered_match_is_not_called_WT(self):
        rows = "SRR010\t100.0\t75.0\t1\t9\tATGGCTGAA\tATGGCTGAA\n"
        classes = self._classes(rows, "part", min_coverage=0.5)
        self.assertEqual(classes, ["WT_partial"])
        self.assertNotIn("WT", [c for c in classes if c == "WT"])


if __name__ == "__main__":
    unittest.main()
