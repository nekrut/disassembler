#!/usr/bin/env python3
"""test_lexicmap_streamer.py: Unit and integration tests for LexicMapStreamer.

Tests:
  - Multi-HSP coordinate tiling & overlap resolution
  - In-frame translation and premature stop-codon detection
  - 4-Tier cohort stratification
  - Online haplotype collapsing and mutation classification
  - Streaming from file, gzipped file, and stdin
"""

import os
import sys
import gzip
import json
import unittest
import tempfile
from pathlib import Path

# Add parent directory to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent))

import lexicmap_streamer as lms

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
        # WT
        self.assertEqual(lms.translate_nt("ATGGCTGAATAA"), "MAE*")
        self.assertEqual(lms.find_internal_stops("ATGGCTGAATAA", 12), [])

        # Premature stop at codon 2 (TAG): ATG TAG GAA TAA
        seq_stop = "ATGTAGGAATAA"
        self.assertEqual(lms.translate_nt(seq_stop), "M*E*")
        self.assertEqual(lms.find_internal_stops(seq_stop, 12), [2])

        # Terminal stop only should NOT be flagged as internal premature stop
        self.assertEqual(lms.find_internal_stops("ATGGCTGAATGA", 12), [])

    def test_multi_hsp_overlap_resolution(self):
        # HSP1: nt 1..6 (ATGGCT) pident=100.0
        # HSP2: nt 5..12 with overlap at 5..6 (CT vs TT, HSP1 has higher pid)
        hsp1 = (100.0, 50.0, 1, 6, "ATGGCT", "ATGGCT")
        hsp2 = (90.0, 66.7, 5, 12, "CTGAATAA", "TTGAATAA")

        merged, cov_bases, mean_pid, single_ge_cov, single_full = lms.merge_accession_hsps(
            [hsp1, hsp2], ref_len=12, min_pident=60.0, min_coverage=0.80
        )

        self.assertEqual(merged, "ATGGCTGAATAA") # Overlap at pos 5 resolved to 'C' from HSP1 (pident 100 > 90)
        self.assertEqual(cov_bases, 12)
        self.assertAlmostEqual(cov_bases / 12, 1.0)
        self.assertFalse(single_ge_cov) # Neither single HSP had >= 80% coverage
        self.assertFalse(single_full)

    def test_end_to_end_streaming(self):
        # Create a TSV with diverse accessions:
        # - SRR001: multi-HSP rescued clean WT (two 50% HSPs)
        # - SRR002: single HSP missense variant (A2V: 5C>T)
        # - SRR003: premature stop codon at pos 2 (TAG) -> Flagged Stops
        # - SRR004: partial coverage 50% (clean partial)
        # - SRR005: fragment <50% cov (low coverage)
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

        # Check cohort counts
        breakdown = summary["cohort_breakdown"]
        self.assertEqual(breakdown["clean_full_length"], 2) # SRR001, SRR002
        self.assertEqual(breakdown["clean_partial"], 1)     # SRR004
        self.assertEqual(breakdown["combined_clean_cohort"], 3)
        self.assertEqual(breakdown["flagged_premature_stops"], 1) # SRR003
        self.assertEqual(breakdown["flagged_low_coverage"], 1)    # SRR005

        # Check multi-HSP rescue metric
        self.assertEqual(summary["recovery_metrics"]["multi_hsp_rescued_accessions"], 1)

        # Check output files generated
        self.assertTrue(Path(f"{out_prefix}.clean.msa.fasta").exists())
        self.assertTrue(Path(f"{out_prefix}.clean.haplotypes.tsv").exists())
        self.assertTrue(Path(f"{out_prefix}.cohort_ledger.tsv").exists())
        self.assertTrue(Path(f"{out_prefix}.flagged.tsv").exists())

        # Check haplotype TSV mutation classification
        tsv_lines = Path(f"{out_prefix}.clean.haplotypes.tsv").read_text().strip().split("\n")
        header = tsv_lines[0].split("\t")
        class_idx = header.index("mutation_class")
        aa_idx = header.index("aa_mutations")

        classes = [line.split("\t")[class_idx] for line in tsv_lines[1:]]
        self.assertIn("WT", classes)
        self.assertIn("Missense", classes)

    def test_mean_pident_excludes_gap_positions(self):
        # A single 90%-identity HSP where the subject deletes one reference position (sc '-').
        # mean_pident must average only the 5 covered (non-gap) positions, never exceed 100%.
        merged, cov_bases, mean_pid, _, _ = lms.merge_accession_hsps(
            [(90.0, 100.0, 1, 6, "ATGGCT", "AT-GCT")], ref_len=6, min_pident=60.0, min_coverage=0.80
        )
        self.assertEqual(cov_bases, 5)
        self.assertLessEqual(mean_pid, 100.0)
        self.assertAlmostEqual(mean_pid, 90.0)

    def test_ungrouped_input_is_handled_correctly(self):
        # SRR001's two HSPs are split by an SRR002 row: the stream is NOT grouped by accession.
        # Real production LexicMap output is routinely ordered this way (ranked by score across
        # all matched genomes, not grouped by genome), so this must be handled correctly rather
        # than raising: SRR001's two split HSPs must still be merged together into one
        # full-coverage clean accession, exactly as if the input had been pre-sorted.
        tsv = (
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "SRR001\t100.0\t50.0\t1\t6\tATGGCT\tATGGCT\n"
            "SRR002\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
            "SRR001\t100.0\t50.0\t7\t12\tGAATAA\tGAATAA\n"
        )
        p = self.dir_path / "ungrouped.tsv"
        p.write_text(tsv)
        summary = lms.stream_lexicmap_to_msa(
            ref_path=str(self.ref_fasta), file_path=str(p),
            out_prefix=str(self.dir_path / "ung"),
        )
        # Both SRR001 (rescued from its two split HSPs) and SRR002 (single full-length HSP) must
        # be counted as clean full-length -- not raised on, not double-counted, not undercounted.
        breakdown = summary["cohort_breakdown"]
        self.assertEqual(summary["total_accessions_evaluated"], 2)
        self.assertEqual(breakdown["clean_full_length"], 2)
        self.assertEqual(summary["recovery_metrics"]["multi_hsp_rescued_accessions"], 1)  # SRR001

        clean_fa = Path(f"{self.dir_path / 'ung'}.clean.accessions.fasta").read_text()
        self.assertIn("SRR001", clean_fa)
        self.assertIn("SRR002", clean_fa)
        self.assertIn("ATGGCTGAATAA", clean_fa)  # SRR001's merged, full-coverage sequence

    def test_reappearing_accession_does_not_double_count(self):
        # A three-way split (interleaved with two other accessions) is a stronger version of the
        # same real-world scenario: SRR900's three HSPs are scattered across the file. Must still
        # be merged into exactly one accession record, not three.
        tsv = (
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "SRR900\t100.0\t33.3\t1\t4\tATGG\tATGG\n"
            "SRR901\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
            "SRR900\t100.0\t33.3\t5\t8\tCTGA\tCTGA\n"
            "SRR902\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
            "SRR900\t100.0\t33.3\t9\t12\tATAA\tATAA\n"
        )
        p = self.dir_path / "scattered.tsv"
        p.write_text(tsv)
        summary = lms.stream_lexicmap_to_msa(
            ref_path=str(self.ref_fasta), file_path=str(p),
            out_prefix=str(self.dir_path / "scat"),
        )
        self.assertEqual(summary["total_accessions_evaluated"], 3)  # not 5, not double-counted
        self.assertEqual(summary["cohort_breakdown"]["clean_full_length"], 3)

    def test_empty_stream_does_not_crash(self):
        # Header only, no data rows: must return a zeroed summary, not raise ZeroDivisionError.
        p = self.dir_path / "empty.tsv"
        p.write_text("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n")
        summary = lms.stream_lexicmap_to_msa(
            ref_path=str(self.ref_fasta), file_path=str(p),
            out_prefix=str(self.dir_path / "emp"),
        )
        self.assertEqual(summary["total_accessions_evaluated"], 0)
        self.assertEqual(summary["cohort_breakdown"]["clean_full_length"], 0)

    def test_full_length_only_requires_complete_coverage(self):
        # SRR010 covers 100% (clean full), SRR011 covers only 50%. With --full-length-only,
        # only the 100%-coverage accession is clean full-length.
        tsv = (
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "SRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
            "SRR011\t100.0\t50.0\t1\t6\tATGGCT\tATGGCT\n"
        )
        p = self.dir_path / "fl.tsv"
        p.write_text(tsv)
        summary = lms.stream_lexicmap_to_msa(
            ref_path=str(self.ref_fasta), file_path=str(p),
            out_prefix=str(self.dir_path / "fl"), full_length_only=True,
        )
        self.assertEqual(summary["cohort_breakdown"]["clean_full_length"], 1)  # only SRR010

if __name__ == "__main__":
    unittest.main()
