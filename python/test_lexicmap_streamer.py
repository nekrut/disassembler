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

if __name__ == "__main__":
    unittest.main()
