#!/usr/bin/env python3
"""test_lexicmap.py: Unit and integration tests for the lexicmap subcommand.

Tests:
  - Multi-HSP coordinate tiling & overlap resolution
  - In-frame translation and premature stop-codon detection
  - 4-tier cohort partitioning
  - Online haplotype collapsing and mutation classification
  - Streaming from file and stdin
  - Input opening (gzip by content, pipes), Galaxy dataset state, API key order
"""

import gzip
import os
import json
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

        # Straight at `tile_accession_hsps`: the `merge_accession_hsps` wrapper these two
        # tests used existed only for them, and a test-only entry point is one more thing
        # that can agree with the tests while the production path does not (PL-28).
        tiling = lms.tile_accession_hsps(
            [hsp1, hsp2], ref_len=12, min_pident=60.0, min_coverage=0.80
        )
        cov_bases = sum(1 for c in tiling.seq if c != '-')

        self.assertEqual(tiling.seq, "ATGGCTGAATAA")
        self.assertEqual(cov_bases, 12)
        self.assertAlmostEqual(cov_bases / 12, 1.0)
        self.assertFalse(tiling.single_ge_min_cov)
        self.assertFalse(tiling.single_full)

    def test_the_cohort_ledger_is_opt_in_and_matches_the_records(self):
        """`--cohort-ledger` renders records[] as a TSV, and the two must agree.

        The ledger used to be written unconditionally, which meant every accession was on
        disk twice once records[] existed. It is kept for Galaxy and for spreadsheets, so
        when it is asked for it has to carry exactly the same rows.
        """
        tsv_path = self.dir_path / "ledger_in.tsv"
        tsv_path.write_text(
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "SRR001\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
            "SRR002\t91.7\t100.0\t1\t12\tATGGCTGAATAA\tATGGTTGAATAA\n"
        )
        prefix = str(self.dir_path / "ledger_out")
        lms.stream_lexicmap_to_msa(ref_path=str(self.ref_fasta), file_path=str(tsv_path),
                                   out_prefix=prefix, write_ledger=True)

        ledger_path = Path(f"{prefix}.cohort_ledger.tsv")
        self.assertTrue(ledger_path.is_file())
        lines = ledger_path.read_text().splitlines()
        header = lines[0].split("\t")
        rows = [dict(zip(header, l.split("\t"))) for l in lines[1:]]
        records = {r["acc"]: r for r in
                   json.loads(Path(f"{prefix}.summary.json").read_text())["records"]}

        self.assertEqual(len(rows), len(records))
        for row in rows:
            record = records[row["accession"]]
            self.assertEqual(row["status"], record["status"])
            self.assertEqual(int(row["tier"]), record["tier"])
            self.assertEqual(float(row["coverage"]), record["coverage"])
            self.assertEqual(float(row["mean_pident"]), record["mean_pident"])
            self.assertEqual(int(row["hsps_merged"]), record["hsps_merged"])
            self.assertEqual(int(row["n_stops"]), record["n_stops"])
            self.assertEqual(row["indels"], record["indels"])

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

        breakdown = summary["cohort"]
        self.assertEqual(breakdown["clean_full_length"], 2)
        self.assertEqual(breakdown["clean_partial"], 1)
        self.assertEqual(breakdown["combined_clean_cohort"], 3)
        # The flagged tiers are discards, and live with the rest of them.
        flagged = summary["discards"]["accessions"]
        self.assertEqual(flagged["flagged_premature_stops"], 1)
        self.assertEqual(flagged["flagged_low_coverage"], 1)
        self.assertEqual(breakdown["flagged"], 2)
        # Counted but never recorded before 0.2.0.
        # The whole dict, not a subset: a row that is dropped has to be counted somewhere,
        # and a new reason for dropping one has to appear here too (PL-14, PL-20).
        self.assertEqual(summary["discards"]["rows"], {
            "malformed": 0, "other_query": 0, "qlen_mismatch": 0,
            "too_few_columns": 0, "inconsistent_qseq": 0})

        self.assertEqual(summary["cohort"]["multi_hsp_rescued_accessions"], 1)

        self.assertTrue(Path(f"{out_prefix}.clean.msa.fasta").exists())
        self.assertTrue(Path(f"{out_prefix}.clean.haplotypes.tsv").exists())
        # Off by default: records[] already carries every one of its rows.
        self.assertFalse(Path(f"{out_prefix}.cohort_ledger.tsv").exists())
        # flagged.tsv was a strict subset of the ledger; records[] carries every accession
        # with its status, so it is gone rather than duplicated a third time.
        self.assertFalse(Path(f"{out_prefix}.flagged.tsv").exists())
        # Read records from the file on disk: the returned dict deliberately omits them,
        # since they are unbounded in the cohort size and stream straight to the summary.
        on_disk = json.loads(Path(f"{out_prefix}.summary.json").read_text())
        self.assertEqual(len(on_disk["records"]), 5)

        tsv_lines = Path(f"{out_prefix}.clean.haplotypes.tsv").read_text().strip().split("\n")
        header = tsv_lines[0].split("\t")
        class_idx = header.index("mutation_class")

        classes = [line.split("\t")[class_idx] for line in tsv_lines[1:]]
        self.assertIn("WT", classes)
        self.assertIn("Missense", classes)

    def test_out_dir_that_does_not_exist_yet_is_a_directory(self):
        # `--out-dir` was a third spelling of --out-prefix, so `--out-dir out/new` wrote
        # out/new.clean.msa.fasta beside a directory that was never made.
        from disassembler import cli
        tsv = self.dir_path / "hits.tsv"
        tsv.write_text(
            "query\tqlen\thits\tsgenome\tsseqid\tqcovGnm\tcls\thsp\tqcovHSP\talenHSP\tpident\tgaps\tqstart\tqend\tsstart\tsend\tqseq\tsseq\n"
            "test_cds\t12\t1\tSRR002\tu3\t100.0\t1\t1\t100.0\t12\t100.0\t0\t1\t12\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n")
        out = self.dir_path / "out" / "new"
        rc = cli.main(["lexicmap", "-q", str(self.ref_fasta), "-f", str(tsv),
                       "--out-dir", str(out)])
        self.assertEqual(rc, 0)
        self.assertTrue(out.is_dir())
        self.assertTrue(any(p.name.endswith(".clean.msa.fasta") for p in out.iterdir()))
        self.assertFalse((self.dir_path / "out" / "new.clean.msa.fasta").exists())
        # -o keeps its prefix meaning.
        rc = cli.main(["lexicmap", "-q", str(self.ref_fasta), "-f", str(tsv),
                       "-o", str(self.dir_path / "pfx")])
        self.assertEqual(rc, 0)
        self.assertTrue((self.dir_path / "pfx.clean.msa.fasta").exists())

    def test_gzip_is_detected_by_content_not_by_name(self):
        # Galaxy hands every dataset over as dataset_<uuid>.dat, compressed or not.
        tsv = (
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            "SRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
        )
        gzip_named_dat = self.dir_path / "dataset_1.dat"
        with gzip.open(gzip_named_dat, "wt") as fh:
            fh.write(tsv)
        plain_named_gz = self.dir_path / "dataset_2.tsv.gz"
        plain_named_gz.write_text(tsv)

        for path in (gzip_named_dat, plain_named_gz):
            with self.subTest(path=path.name):
                summary = lms.stream_lexicmap_to_msa(
                    ref_path=str(self.ref_fasta),
                    file_path=str(path),
                    out_prefix=str(self.dir_path / path.name),
                )
                self.assertEqual(summary["cohort"]["clean_full_length"], 1)

    def test_mean_pident_excludes_gap_positions(self):
        tiling = lms.tile_accession_hsps(
            [(90.0, 100.0, 1, 6, "ATGGCT", "AT-GCT")], ref_len=6, min_pident=60.0, min_coverage=0.80
        )
        self.assertEqual(sum(1 for c in tiling.seq if c != '-'), 5)
        self.assertLessEqual(tiling.mean_pident, 100.0)
        self.assertAlmostEqual(tiling.mean_pident, 90.0)

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
        self.assertEqual(summary["run"]["n_records"], 0)
        self.assertEqual(summary["cohort"]["clean_full_length"], 0)

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

        self.assertEqual(on["cohort"]["clean_full_length"], 1)
        self.assertEqual(off["cohort"]["clean_full_length"], 2)
        self.assertNotEqual(
            on["cohort"]["clean_full_length"],
            off["cohort"]["clean_full_length"],
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
        self.assertEqual(summary["cohort"]["clean_full_length"], 1)
        self.assertTrue((self.dir_path / "na.summary.json").is_file())

    def test_an_embedded_header_row_is_skipped(self):
        tsv = self.header + self.good + self.header.strip() + "\n" + self.good.replace("SRR010", "SRR012")
        summary = self._run(tsv, "hdr")
        self.assertEqual(summary["run"]["n_records"], 2)

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
        self.assertEqual(summary["run"]["n_records"], 1)
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
        self.assertEqual(summary["run"]["n_records"], 2, "both geneA accessions")

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
        self.assertEqual(summary["cohort"]["clean_full_length"], 1)

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
        self.assertEqual(summary["cohort"]["clean_full_length"], 1)


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
        # qseq is the REFERENCE over the span, sseq the sample: these fixtures used to
        # put the mutant in both, which is what PL-20 found nothing was checking.
        rows = "SRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGGTGAATAA\n"
        self.assertEqual(self._classes(rows, "mis"), ["Missense"])

    def test_synonymous_and_nonsense_still_classify(self):
        syn = "SRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAGTAA\n"
        self.assertEqual(self._classes(syn, "syn"), ["Synonymous"])
        non = "SRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGTGAGAATAA\n"
        self.assertEqual(self._classes(non, "non", max_internal_stops=1), ["Nonsense"])

    def test_an_exact_match_is_WT(self):
        rows = "SRR010\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
        self.assertEqual(self._classes(rows, "wt"), ["WT"])

    def test_a_partially_covered_match_is_not_called_WT(self):
        rows = "SRR010\t100.0\t75.0\t1\t9\tATGGCTGAA\tATGGCTGAA\n"
        classes = self._classes(rows, "part", min_coverage=0.5)
        self.assertEqual(classes, ["WT_partial"])
        self.assertNotIn("WT", [c for c in classes if c == "WT"])


# 60 nt CDS: ATG, 18 sense codons, TAA. Long enough to place an indel more than
# INDEL_EDGE_WINDOW aligned bases from either end of its HSP.
CDS60 = "ATG" + "GCTGAAAAACCCGGTTTCCATGACATCAGCCTGTGGCGTACCGTACAGAACTAC" + "TAA"


def _ins(after, bases):
    """(qseq, sseq) for CDS60 with `bases` inserted after 0-based position `after`."""
    return (CDS60[:after + 1] + "-" * len(bases) + CDS60[after + 1:],
            CDS60[:after + 1] + bases + CDS60[after + 1:])


def _del(start, length):
    """(qseq, sseq) for CDS60 with `length` bases deleted from 0-based `start`."""
    return CDS60, CDS60[:start] + "-" * length + CDS60[start + length:]


class TestIndels(unittest.TestCase):
    """Insertions and deletions read from the aligned qseq/sseq. Before this, the tiling dropped
    inserted bases silently and the stop scan skipped codons containing a deletion, so a 1-nt
    insertion came out CLEAN_FULL_LENGTH / WT."""

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.dir_path = Path(self.test_dir.name)
        self.ref = self.dir_path / "cds60.fasta"
        self.ref.write_text(f">cds60\n{CDS60}\n")
        self.assertEqual(len(CDS60), 60)
        self.assertEqual(lms.find_internal_stops(CDS60, 60), [])

    def tearDown(self):
        self.test_dir.cleanup()

    def _run(self, rows, name, **kw):
        """rows: list of (accession, pident, qseq, sseq), each one full-length HSP."""
        lines = ["sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq"]
        for acc, pid, qseq, sseq in rows:
            lines.append(f"{acc}\t{pid}\t100.0\t1\t60\t{qseq}\t{sseq}")
        tsv = self.dir_path / f"{name}.tsv"
        tsv.write_text("\n".join(lines) + "\n")
        prefix = str(self.dir_path / name)
        summary = lms.stream_lexicmap_to_msa(ref_path=str(self.ref), file_path=str(tsv),
                                             out_prefix=prefix, **kw)
        # records[] in the summary is the per-accession layer now; cohort_ledger.tsv is an
        # opt-in rendering of these same rows (--cohort-ledger).
        on_disk = json.loads(Path(f"{prefix}.summary.json").read_text())
        ledger = {r["acc"]: r for r in on_disk["records"]}
        with open(f"{prefix}.clean.haplotypes.tsv") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            haps = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]
        return summary, ledger, haps, prefix

    def test_notation(self):
        self.assertEqual(lms.Indel("ins", 685, 1, "T", False).notation(), "685_686insT")
        self.assertEqual(lms.Indel("del", 2315, 6, "ACGTAC", False).notation(), "2316_2321del")
        self.assertEqual(lms.Indel("del", 708, 1, "C", True).notation(), "709del(edge)")

    def test_one_nt_insertion_is_a_frameshift_not_wild_type(self):
        qseq, sseq = _ins(29, "G")
        _, ledger, _, _ = self._run([("SRR1", 98.4, qseq, sseq)], "ins1")
        self.assertEqual(ledger["SRR1"]["status"], "FLAGGED_FRAMESHIFT")
        self.assertEqual(ledger["SRR1"]["tier"], 5)
        self.assertEqual(ledger["SRR1"]["indels"], "30_31insG")

    def test_one_nt_deletion_is_a_frameshift(self):
        qseq, sseq = _del(30, 1)
        _, ledger, _, _ = self._run([("SRR1", 98.3, qseq, sseq)], "del1")
        self.assertEqual(ledger["SRR1"]["status"], "FLAGGED_FRAMESHIFT")
        self.assertEqual(ledger["SRR1"]["indels"], "31del")

    def test_in_frame_deletion_counts_as_covered_under_require_full_length(self):
        # The TAC1B case: a real two-codon deletion used to be read as missing coverage, so
        # --require-full-length demoted every carrier to the partial tier.
        qseq, sseq = _del(30, 6)
        _, ledger, haps, _ = self._run([("SRR1", 100.0, qseq, sseq)], "del6",
                                       full_length_only=True)
        self.assertEqual(ledger["SRR1"]["status"], "CLEAN_FULL_LENGTH")
        self.assertEqual(ledger["SRR1"]["coverage"], 1.0)
        self.assertEqual(ledger["SRR1"]["indels"], "31_36del")
        self.assertEqual(haps[0]["mutation_class"], "InFrameIndel")
        self.assertIn("31_36del", haps[0]["nt_mutations"].split(","))
        self.assertEqual(haps[0]["coverage_pct"], "100.00")

    def test_in_frame_insertion_is_its_own_haplotype(self):
        qseq, sseq = _ins(29, "GCA")
        summary, ledger, haps, prefix = self._run(
            [("WT1", 100.0, CDS60, CDS60), ("INS", 100.0, qseq, sseq), ("WT2", 100.0, CDS60, CDS60)],
            "ins3")
        self.assertEqual(ledger["INS"]["status"], "CLEAN_FULL_LENGTH")
        self.assertEqual(summary["cohort"]["clean_full_length"], 3)
        self.assertEqual(summary["cohort"]["clean_unique_haplotypes"], 2)
        self.assertEqual(summary["cohort"]["wildtype_count"], 2)
        by_class = {h["mutation_class"]: h for h in haps}
        self.assertEqual(by_class["WT"]["count"], "2")
        self.assertEqual(by_class["InFrameIndel"]["count"], "1")
        self.assertEqual(by_class["InFrameIndel"]["nt_mutations"], "30_31insGCA")
        # The rows are identical (an insertion has no column); the header tells them apart.
        msa = Path(f"{prefix}.clean.msa.fasta").read_text().split(">")[1:]
        self.assertEqual(len(msa), 3)
        rows = [rec.split("\n", 1)[1].strip() for rec in msa[1:]]
        self.assertEqual(rows[0], rows[1])
        self.assertTrue(any("indels=30_31insGCA" in rec.split("\n", 1)[0] for rec in msa))

    def test_stop_codon_formed_by_an_insertion_is_found(self):
        qseq, sseq = _ins(29, "TAA")
        _, ledger, _, _ = self._run([("SRR1", 95.2, qseq, sseq)], "insstop")
        self.assertEqual(ledger["SRR1"]["status"], "FLAGGED_PREMATURE_STOPS")
        self.assertEqual(ledger["SRR1"]["stop_codons"], "11")

    def test_edge_indel_is_reported_but_not_trusted(self):
        qseq, sseq = _ins(4, "T")   # 5 aligned bases before it, under INDEL_EDGE_WINDOW
        summary, ledger, _, _ = self._run(
            [("EDGE", 99.0, qseq, sseq), ("WT", 100.0, CDS60, CDS60)], "edge")
        self.assertEqual(ledger["EDGE"]["status"], "CLEAN_FULL_LENGTH")
        self.assertEqual(ledger["EDGE"]["indels"], "5_6insT(edge)")
        self.assertEqual(summary["cohort"]["clean_unique_haplotypes"], 1)
        self.assertEqual(summary["cohort"]["wildtype_count"], 2)

    def test_allow_frameshift_alts_waives_only_the_frameshift_tier(self):
        qseq, sseq = _ins(29, "G")
        stop = CDS60[:30] + "TAG" + CDS60[33:]
        _, ledger, haps, _ = self._run(
            [("FS", 98.4, qseq, sseq), ("STOP", 98.3, CDS60, stop)], "allowfs",
            allow_frameshifts=True)
        self.assertEqual(ledger["FS"]["status"], "CLEAN_FULL_LENGTH")
        self.assertEqual(haps[0]["mutation_class"], "Frameshift")
        self.assertEqual(ledger["STOP"]["status"], "FLAGGED_PREMATURE_STOPS")

    @staticmethod
    def _ins_then_del(after, base, del_at):
        """(qseq, sseq): `base` inserted after `after`, and the base at `del_at` deleted."""
        q = CDS60[:after + 1] + "-" + CDS60[after + 1:]
        s = CDS60[:after + 1] + base + CDS60[after + 1:del_at] + "-" + CDS60[del_at + 1:]
        return q, s

    def test_split_by_frame(self):
        ins, dele = lms.Indel("ins", 30, 1, "G", False), lms.Indel("del", 35, 1, "G", False)
        self.assertEqual(lms.split_by_frame([ins, dele]), ([ins, dele], []))
        far = lms.Indel("del", 30 + lms.FRAMESHIFT_MIN_SPAN, 1, "G", False)
        self.assertEqual(lms.split_by_frame([ins, far]), ([], [ins, far]))
        self.assertEqual(lms.split_by_frame([ins]), ([], [ins]))
        three = lms.Indel("del", 40, 3, "GTG", False)
        self.assertEqual(lms.split_by_frame([three]), ([three], []))
        self.assertEqual(lms.split_by_frame([lms.Indel("ins", 30, 1, "G", True)]), ([], []))

    def test_indels_that_restore_the_frame_quickly_are_not_a_frameshift(self):
        # The H. pylori babA case: +1 then -1 a few bases later. The aligner's way of writing a
        # divergent stretch, not a broken gene.
        qseq, sseq = self._ins_then_del(29, "G", 35)
        _, ledger, haps, _ = self._run([("SRR1", 98.3, qseq, sseq)], "comp")
        self.assertEqual(ledger["SRR1"]["status"], "CLEAN_FULL_LENGTH")
        self.assertEqual(ledger["SRR1"]["indels"], "30_31insG,36del")
        self.assertEqual(haps[0]["mutation_class"], "InFrameIndel")

    def test_stop_inside_a_briefly_shifted_stretch_is_found(self):
        # Inserting T before AGC reads TAG in the sample's own frame.
        qseq, sseq = self._ins_then_del(29, "T", 35)
        _, ledger, _, _ = self._run([("SRR1", 98.3, qseq, sseq)], "compstop")
        self.assertEqual(ledger["SRR1"]["status"], "FLAGGED_PREMATURE_STOPS")
        self.assertEqual(ledger["SRR1"]["stop_codons"], "11")

    def test_a_shift_restored_only_after_the_span_is_a_frameshift(self):
        qseq, sseq = self._ins_then_del(11, "G", 12 + lms.FRAMESHIFT_MIN_SPAN)
        _, ledger, _, _ = self._run([("SRR1", 98.3, qseq, sseq)], "longshift")
        self.assertEqual(ledger["SRR1"]["status"], "FLAGGED_FRAMESHIFT")

    def test_an_indel_counts_only_where_its_hsp_supplied_the_row(self):
        qseq, sseq = _del(30, 1)
        wt = (99.0, 100.0, 1, 60, CDS60, CDS60)
        gapped = (90.0, 100.0, 1, 60, qseq, sseq)
        t = lms.tile_accession_hsps([wt, gapped], 60, 60.0, 0.8)
        self.assertEqual(t.indels, [])
        self.assertEqual(t.seq, CDS60)
        t = lms.tile_accession_hsps([(99.0, 100.0, 1, 60, qseq, sseq), (90.0, 100.0, 1, 60, CDS60, CDS60)],
                                    60, 60.0, 0.8)
        self.assertEqual([e.notation() for e in t.indels], ["31del"])
        self.assertEqual(t.covered, 60)



TABLE = ("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
         "SRR10\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
         "SRR11\t91.7\t100.0\t1\t12\tATGGCTGAATAA\tATGGTTGAATAA\n")


class TestOpenTable(unittest.TestCase):
    """One opener (`textio.text_stream`) for files, stdin, pipes and Galaxy downloads; gzip
    decided by content."""

    def read(self, raw):
        from disassembler.textio import text_stream
        return text_stream(raw).read()

    def test_plain_and_gzip(self):
        import io
        self.assertEqual(self.read(io.BytesIO(TABLE.encode())), TABLE)
        self.assertEqual(self.read(io.BytesIO(gzip.compress(TABLE.encode()))), TABLE)

    def test_a_pipe_loses_nothing_to_the_gzip_check(self):
        # The old check reopened the path to look at the first bytes, which emptied a pipe
        # (`-f <(zcat x)`) before it was read, compressed or not.
        import os
        for payload in (TABLE.encode(), gzip.compress(TABLE.encode())):
            r, w = os.pipe()
            os.write(w, payload)
            os.close(w)
            with os.fdopen(r, "rb") as fh:
                self.assertEqual(self.read(fh), TABLE)

    def test_a_stream_without_peek_is_buffered(self):
        import io

        class Raw(io.RawIOBase):  # like an HTTP response body
            def __init__(self, data):
                self.data = io.BytesIO(data)

            def readable(self):
                return True

            def readinto(self, b):
                return self.data.readinto(b)

        self.assertEqual(self.read(Raw(gzip.compress(TABLE.encode()))), TABLE)

    def test_open_text_names_the_file_it_refuses(self):
        from disassembler.textio import open_text
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.fasta"
            p.write_bytes(b"\x28\xb5\x2f\xfd" + b"\0" * 16)
            with self.assertRaisesRegex(ValueError, "x.fasta: zstd-compressed"):
                open_text(p)

    def test_zstd_is_named_not_garbled(self):
        import io
        from disassembler.textio import text_stream
        with self.assertRaisesRegex(ValueError, "zstd-compressed"):
            text_stream(io.BytesIO(b"\x28\xb5\x2f\xfd" + b"\0" * 16))


class TestGalaxyDataset(unittest.TestCase):
    """`-d`: only an `ok` dataset is read. The state lookup and the download are stubbed."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.ref = self.dir / "gene.fasta"
        self.ref.write_text(">gene\nATGGCTGAATAA\n")
        self.real = (lms.resolve_dataset_state, lms.get_galaxy_stream)
        self.downloaded = False

        def fake_stream(*a, **k):
            import io
            self.downloaded = True
            return io.BytesIO(gzip.compress(TABLE.encode()))

        lms.get_galaxy_stream = fake_stream

    def tearDown(self):
        lms.resolve_dataset_state, lms.get_galaxy_stream = self.real
        self.tmp.cleanup()

    def run_with_state(self, state):
        import contextlib
        import io
        from disassembler.stage import stage_args
        if isinstance(state, Exception):
            def fake_state(*a, **k):
                raise state
        else:
            def fake_state(*a, **k):
                return state
        lms.resolve_dataset_state = fake_state
        ns = stage_args(lms, ref=str(self.ref), dataset="f2db41e1fa331b3e",
                        out_prefix=str(self.dir / "out") + "/")
        err = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
            rc = lms.run(ns)
        return rc, err.getvalue()

    def test_an_ok_dataset_is_read_gzipped_or_not(self):
        rc, err = self.run_with_state("ok")
        self.assertEqual(rc, 0, err)
        clean = (self.dir / "out" / "gene.clean.accessions.fasta").read_text()
        self.assertEqual(clean.count(">"), 2)

    def test_a_gzipped_reference(self):
        self.ref.write_bytes(gzip.compress(b">gene\nATGGCTGAATAA\n"))
        rc, err = self.run_with_state("ok")
        self.assertEqual(rc, 0, err)
        self.assertTrue((self.dir / "out" / "gene.clean.accessions.fasta").is_file())

    def test_a_queued_dataset_is_a_temporary_failure_not_a_success(self):
        # Galaxy's own `Dataset.non_ready_states`. `paused` and `deferred` were in this list
        # and are not states that become `ok` on their own (PL-23).
        for state in ("queued", "running", "new", "upload", "setting_metadata"):
            rc, err = self.run_with_state(state)
            self.assertEqual(rc, lms.EX_TEMPFAIL, state)
            self.assertIn(f"is '{state}', not finished yet", err)
        self.assertFalse(self.downloaded, "nothing may be streamed from an unfinished dataset")
        self.assertFalse((self.dir / "out").exists())

    def test_a_dataset_that_needs_a_person_is_an_error_not_a_retry(self):
        # A script retrying on 75, which docs/lexicmap.md invites, looped until its own
        # limit: a deferred dataset is terminal and a paused job waits for a person.
        for state, word in (("deferred", "materialize"), ("paused", "resume")):
            rc, err = self.run_with_state(state)
            self.assertEqual(rc, 1, f"{state}: {err}")
            self.assertIn(word, err)
        self.assertFalse(self.downloaded, "nothing may be streamed from either")

    def test_a_failed_dataset_is_an_error(self):
        rc, err = self.run_with_state("error")
        self.assertEqual(rc, 1)
        self.assertIn("in state 'error'; only an 'ok' dataset can be read", err)
        self.assertFalse(self.downloaded)

    def test_an_unreadable_state_is_an_error_not_ignored(self):
        rc, err = self.run_with_state(RuntimeError("could not read the state: HTTP Error 403"))
        self.assertEqual(rc, 1)
        self.assertIn("HTTP Error 403", err)
        self.assertFalse(self.downloaded)


class TestGalaxyApiKeyOrder(unittest.TestCase):
    """--api-key, then --api-key-file, then $GALAXY_API_KEY: the help's order, and harvest's."""

    def test_order(self):
        import os
        from disassembler.cmd import harvest
        with tempfile.TemporaryDirectory() as d:
            key_file = Path(d) / "key.txt"
            key_file.write_text("from-file\n")
            old = os.environ.get("GALAXY_API_KEY")
            os.environ["GALAXY_API_KEY"] = "from-env"
            try:
                self.assertEqual(lms.get_galaxy_api_key("literal", str(key_file)), "literal")
                self.assertEqual(lms.get_galaxy_api_key(None, str(key_file)), "from-file")
                self.assertEqual(harvest.get_galaxy_api_key(str(key_file)), "from-file")
                self.assertEqual(lms.get_galaxy_api_key(None, None), "from-env")
                self.assertEqual(harvest.get_galaxy_api_key(None), "from-env")
            finally:
                if old is None:
                    del os.environ["GALAXY_API_KEY"]
                else:
                    os.environ["GALAXY_API_KEY"] = old


class TestCoverageTakesPrecedence(unittest.TestCase):
    """A fragment is reported as a fragment, whatever its frame and stop codons say.

    Through 0.1.1-dev the frame and stop rules were checked before coverage, so an accession
    with almost none of the gene tiled still came out FLAGGED_FRAMESHIFT or
    FLAGGED_PREMATURE_STOPS. On the SARS-CoV-2 strict pass that was 61,897 of 64,034
    "frameshifts", at a median coverage of 0.18, which is far too little of a gene to tell a
    real shift from a gap in the tiling.
    """

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.dir_path = Path(self.test_dir.name)
        self.ref = self.dir_path / "cds60.fasta"
        self.ref.write_text(f">cds60\n{CDS60}\n")

    def tearDown(self):
        self.test_dir.cleanup()

    def _status(self, name, qstart, qend, qseq, sseq):
        """One HSP covering reference bases qstart..qend (1-based, inclusive)."""
        tsv = self.dir_path / f"{name}.tsv"
        tsv.write_text(
            "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
            f"SRR1\t98.0\t{100.0 * (qend - qstart + 1) / len(CDS60):.1f}\t{qstart}\t{qend}\t{qseq}\t{sseq}\n")
        prefix = str(self.dir_path / name)
        lms.stream_lexicmap_to_msa(ref_path=str(self.ref), file_path=str(tsv), out_prefix=prefix)
        on_disk = json.loads(Path(f"{prefix}.summary.json").read_text())
        return on_disk["records"][0]["status"]

    def test_a_fragment_with_a_frameshift_is_a_fragment(self):
        # 20 of 60 bases tiled (0.33, under the 0.50 partial floor) with a 1-nt insertion.
        frag = CDS60[:20]
        qseq = frag[:10] + "-" + frag[10:]
        sseq = frag[:10] + "G" + frag[10:]
        self.assertEqual(self._status("frag_shift", 1, 20, qseq, sseq), "LOW_COVERAGE_FRAGMENT")

    def test_a_fragment_with_a_premature_stop_is_a_fragment(self):
        # Same 20 bases, with codon 3 (0-based 6..8) turned into TAA.
        frag = CDS60[:20]
        self.assertEqual(self._status("frag_stop", 1, 20, frag, frag[:6] + "TAA" + frag[9:]),
                         "LOW_COVERAGE_FRAGMENT")

    def test_a_well_covered_frameshift_still_flags(self):
        # The control: the same insertion over the whole gene is a frameshift, not a fragment.
        qseq, sseq = _ins(10, "G")
        self.assertEqual(self._status("full_shift", 1, 60, qseq, sseq), "FLAGGED_FRAMESHIFT")

    def test_a_well_covered_premature_stop_still_flags(self):
        self.assertEqual(
            self._status("full_stop", 1, 60, CDS60, CDS60[:6] + "TAA" + CDS60[9:]),
            "FLAGGED_PREMATURE_STOPS")


class TestCoverageIndependentHaplotypes(unittest.TestCase):
    """PL-02: diversity must not move when an accession's contigs stop short.

    Tier 1 rows are collapsed on the whole reference-length row, which carries `-` at every
    untiled position, and with the default `--min-cov 0.80` a tier 1 row can be a fifth `-`.
    Two accessions identical wherever they are covered therefore became two haplotypes, so
    `clean_unique_haplotypes`, `shannon_entropy_bits` and `wildtype_fraction` measured
    assembly fragmentation as much as sequence diversity: the numbers tracked sequencing
    depth, and two cohorts or two targets could not be compared.

    The original fields are unchanged on purpose, so runs already done stay comparable; what
    is asserted here is that the fields beside them do not move with coverage.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def cohort(self, name, ref, rows):
        """`rows` are (accession, qstart, qend, sseq) with qseq taken from the reference."""
        ref_fa = self.dir / f"{name}.fasta"
        ref_fa.write_text(f">{name}\n{ref}\n")
        tsv = self.dir / f"{name}.tsv"
        with tsv.open("w") as fh:
            fh.write("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n")
            for acc, qs, qe, sseq in rows:
                qseq = ref[qs - 1:qe]
                pid = 100.0 * sum(a == b for a, b in zip(qseq, sseq)) / len(qseq)
                cov = 100.0 * (qe - qs + 1) / len(ref)
                fh.write(f"{acc}\t{pid:.2f}\t{cov:.2f}\t{qs}\t{qe}\t{qseq}\t{sseq}\n")
        prefix = str(self.dir / name)
        lms.stream_lexicmap_to_msa(ref_path=str(ref_fa), file_path=str(tsv), out_prefix=prefix)
        return json.loads(Path(f"{prefix}.summary.json").read_text())["cohort"]

    def test_coverage_alone_splits_the_old_fields_and_not_the_new_ones(self):
        # The review's fixture: a 600 nt CDS, five full-length wild types, three wild types
        # whose contigs stop short, and one full-length accession carrying 300T>A.
        ref = "ATG" + "GCT" * 198 + "TAA"
        snv = ref[:299] + "A" + ref[300:]
        rows = [(f"WTFULL{i}", 1, 600, ref) for i in range(5)]
        rows += [("WTSHORT1", 1, 570, ref[:570]),
                 ("WTSHORT2", 31, 600, ref[30:]),
                 ("WTSHORT3", 1, 540, ref[:540]),
                 ("SNVFULL", 1, 600, snv)]
        c = self.cohort("cds600", ref, rows)

        # What the cohort really holds: two sequences, eight accessions carrying one of them.
        # The old fields see five, because three of the eight stop short in three places.
        self.assertEqual(c["clean_full_length"], 9)
        self.assertEqual(c["clean_unique_haplotypes"], 5)
        self.assertAlmostEqual(c["shannon_entropy_bits"], 1.8800, places=3)
        self.assertAlmostEqual(c["wildtype_fraction"], 0.5556, places=3)

        # Full-coverage subset: nothing inferred, only the rows that cover the reference.
        self.assertEqual(c["full_coverage_rows"], 6)
        self.assertEqual(c["clean_unique_haplotypes_full_coverage"], 2)
        self.assertAlmostEqual(c["shannon_entropy_bits_full_coverage"], 0.6500, places=3)

        # Folded: each short row carries 300T reference, which the SNV haplotype contradicts,
        # so all three fold to wild type and nothing is left ambiguous.
        self.assertEqual((c["partial_rows_folded"], c["partial_rows_unresolved"]), (3, 0))
        self.assertEqual(c["clean_unique_haplotypes_folded"], 2)
        self.assertAlmostEqual(c["shannon_entropy_bits_folded"], 0.5033, places=3)
        self.assertEqual(c["wildtype_compatible_count"], 8)
        self.assertAlmostEqual(c["wildtype_fraction_folded"], 0.8889, places=3)

    def test_four_wild_types_with_four_contig_ends_are_one_haplotype(self):
        rows = [("A", 1, 60, CDS60), ("B", 1, 57, CDS60[:57]),
                ("C", 4, 60, CDS60[3:]), ("D", 2, 60, CDS60[1:])]
        c = self.cohort("ends", CDS60, rows)
        self.assertEqual(c["clean_unique_haplotypes"], 4, "the old field splits on the ends")
        self.assertEqual(c["clean_unique_haplotypes_full_coverage"], 1)
        self.assertEqual(c["clean_unique_haplotypes_folded"], 1)
        self.assertEqual(c["shannon_entropy_bits_folded"], 0.0)
        self.assertEqual(c["partial_rows_unresolved"], 0)
        self.assertEqual(c["wildtype_fraction_folded"], 1.0)

    def test_a_partial_that_could_be_either_haplotype_is_never_guessed(self):
        # 55T>C sits at reference position 55; the short row stops at 50, so nothing in it
        # says which of the two haplotypes it is. Folding on a variant signature alone would
        # call it wild type, which is exactly the wrong answer half the time.
        snv = CDS60[:54] + "C" + CDS60[55:]
        rows = [("WT", 1, 60, CDS60), ("SNV", 1, 60, snv), ("SHORT", 1, 50, CDS60[:50])]
        c = self.cohort("ambiguous", CDS60, rows)
        self.assertEqual(c["clean_unique_haplotypes_full_coverage"], 2)
        self.assertEqual((c["partial_rows_folded"], c["partial_rows_unresolved"]), (0, 1))
        # The unresolved row is left out rather than counted as a third haplotype.
        self.assertEqual(c["clean_unique_haplotypes_folded"], 2)
        self.assertEqual(c["wildtype_compatible_count"], 1)

    def test_a_trusted_deletion_counts_as_covered_and_folds_separately(self):
        # A deletion is '-' in the row like an untiled position, so coverage and folding both
        # have to read it as covered, and a row carrying one is not the same haplotype as a
        # row missing those bases through fragmentation.
        qseq = CDS60[:30] + "---" + CDS60[33:]
        sseq = CDS60[:30] + "---" + CDS60[33:]
        rows = [("WT", 1, 60, CDS60),
                ("DEL", 1, 60, CDS60[:30] + "---" + CDS60[33:]),
                ("WTSHORT", 1, 54, CDS60[:54])]
        c = self.cohort("withdel", CDS60, rows)
        self.assertEqual(c["full_coverage_rows"], 2, "the deleted bases still count as covered")
        self.assertEqual(c["clean_unique_haplotypes_full_coverage"], 2)
        # The short wild type carries no deletion, so it can only fold into the wild type.
        self.assertEqual((c["partial_rows_folded"], c["partial_rows_unresolved"]), (1, 0))
        self.assertEqual(c["wildtype_compatible_count"], 2)


CDS300 = "ATG" + "GCT" * 98 + "TAA"


class TestMixedHsps(unittest.TestCase):
    """PL-01: an accession whose own HSPs disagree is reported, not quietly resolved.

    Logan leaves an unpopped SNP bubble as two separate contigs. Both 61 nt branches clear
    LexicMap's default `-l 50 -i 70`, and the reference branch is 100% identical while the
    alternative one is not. The tiler resolved overlaps by whole-HSP `pident`, so the
    reference branch won every position it touched regardless of which allele the sample
    carries: the variant vanished, and a deletion on the losing HSP was dropped as "not
    owned", taking the premature stop or frameshift with it. The accession came out
    CLEAN_FULL_LENGTH wild type with nothing in any column saying so.

    Each fixture below is the same shape: one contig spanning the gene and carrying the real
    change, and one short reference-matching bubble branch sitting on top of it.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.ref = self.dir / "cds300.fasta"
        self.ref.write_text(f">cds300\n{CDS300}\n")

    def tearDown(self):
        self.tmp.cleanup()

    def run_hsps(self, name, rows, **kw):
        """`rows` are (qstart, qend, qseq, sseq); pident is computed from the pair."""
        tsv = self.dir / f"{name}.tsv"
        with tsv.open("w") as fh:
            fh.write("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n")
            for qs, qe, qseq, sseq in rows:
                pid = 100.0 * sum(a == b for a, b in zip(qseq, sseq)) / len(qseq)
                cov = 100.0 * (qe - qs + 1) / len(CDS300)
                fh.write(f"SRR1\t{pid:.2f}\t{cov:.2f}\t{qs}\t{qe}\t{qseq}\t{sseq}\n")
        prefix = str(self.dir / name)
        lms.stream_lexicmap_to_msa(ref_path=str(self.ref), file_path=str(tsv),
                                   out_prefix=prefix, **kw)
        s = json.loads(Path(f"{prefix}.summary.json").read_text())
        return s["records"][0], s

    # The bubble branch: 61 reference bases over 120-180, exactly what Logan emits for an
    # unpopped SNP bubble and comfortably over LexicMap's -l 50.
    BUBBLE = (120, 180, CDS300[119:180], CDS300[119:180])

    def test_a_reference_bubble_no_longer_overrides_the_variant(self):
        snv = CDS300[:149] + "A" + CDS300[150:]
        rec, s = self.run_hsps("snv", [(1, 300, CDS300, snv), self.BUBBLE])

        # The variant survives: position 150 lies 149 bases inside the spanning contig and
        # 29 inside the bubble, so depth decides it, not identity to the reference.
        self.assertEqual(rec["status"], "FLAGGED_MIXED_HSPS")
        self.assertEqual(rec["tier"], 6)
        self.assertEqual(rec["n_conflicts"], 1)
        self.assertEqual(rec["conflict_sites"], "150")
        self.assertEqual(s["cohort"]["clean_full_length"], 0)
        self.assertEqual(s["discards"]["accessions"]["flagged_mixed_hsps"], 1)

    def test_allow_mixed_puts_it_back_in_the_clean_tier_with_the_variant(self):
        snv = CDS300[:149] + "A" + CDS300[150:]
        rec, s = self.run_hsps("snv_allowed", [(1, 300, CDS300, snv), self.BUBBLE],
                               allow_mixed=True)
        self.assertEqual(rec["status"], "CLEAN_FULL_LENGTH")
        self.assertEqual(rec["n_conflicts"], 1, "waived, not unrecorded")
        clean = (self.dir / "snv_allowed.clean.accessions.fasta").read_text().splitlines()[1]
        self.assertEqual(clean[149], "A", "the allele the spanning contig carries")
        self.assertEqual(s["settings"]["allow_mixed"], True)

    def test_a_masked_premature_stop_reaches_its_own_tier(self):
        stopped = CDS300[:150] + "TAA" + CDS300[153:]
        rec, _s = self.run_hsps("stop", [(1, 300, CDS300, stopped), self.BUBBLE])
        # Stops are judged before mixed HSPs, so this is tier 3 and not tier 6: the stop is
        # the more specific finding, and it is the one that used to disappear entirely.
        self.assertEqual(rec["status"], "FLAGGED_PREMATURE_STOPS")
        self.assertEqual(rec["tier"], 3)
        self.assertEqual(rec["n_stops"], 1)

    def test_a_masked_frameshift_reaches_its_own_tier(self):
        rec, _s = self.run_hsps(
            "del", [(1, 300, CDS300, CDS300[:150] + "-" + CDS300[151:]), self.BUBBLE])
        # The deletion is kept because the spanning contig now owns the position it removes;
        # under the old rule the bubble owned it, so the indel was dropped as "not owned".
        self.assertEqual(rec["status"], "FLAGGED_FRAMESHIFT")
        self.assertEqual(rec["tier"], 5)
        self.assertEqual(rec["indels"], "151del")

    def test_a_ragged_end_overlap_is_not_a_conflict(self):
        # Two contigs that merely meet, disagreeing 3 bases from one end of the first and 2
        # from the start of the second. An aligner running out of sequence is not two
        # contigs disagreeing, and flagging it would flag most multi-HSP accessions.
        a_sseq = CDS300[:60]
        b_sseq = CDS300[54:120]
        b_sseq = b_sseq[:2] + ("A" if b_sseq[2] != "A" else "C") + b_sseq[3:]
        rec, _s = self.run_hsps("ragged", [(1, 60, CDS300[:60], a_sseq),
                                           (55, 120, CDS300[54:120], b_sseq)])
        self.assertEqual(rec["n_conflicts"], 0)
        self.assertEqual(rec["conflict_sites"], "None")

    def test_an_interior_insertion_one_contig_lacks_is_a_conflict(self):
        # An insertion occupies no reference position, so it cannot show up as two symbols at
        # a position; it needs its own rule or a contig carrying one and a contig spanning
        # the same junction without one would read as agreement.
        ins_q = CDS300[:150] + "---" + CDS300[150:]
        ins_s = CDS300[:150] + "GGG" + CDS300[150:]
        rec, _s = self.run_hsps("ins", [(1, 300, ins_q, ins_s), self.BUBBLE])
        self.assertGreaterEqual(rec["n_conflicts"], 1)
        self.assertEqual(rec["status"], "FLAGGED_MIXED_HSPS")

    def test_one_hsp_alone_is_never_mixed(self):
        snv = CDS300[:149] + "A" + CDS300[150:]
        rec, _s = self.run_hsps("single", [(1, 300, CDS300, snv)])
        self.assertEqual(rec["status"], "CLEAN_FULL_LENGTH")
        self.assertEqual(rec["n_conflicts"], 0)

    def test_the_ledger_carries_the_conflict_sites(self):
        snv = CDS300[:149] + "A" + CDS300[150:]
        self.run_hsps("ledger", [(1, 300, CDS300, snv), self.BUBBLE], write_ledger=True)
        rows = (self.dir / "ledger.cohort_ledger.tsv").read_text().splitlines()
        self.assertEqual(rows[0].split("\t")[-1], "conflict_sites")
        self.assertEqual(rows[1].split("\t")[-1], "150")


CDS120 = "ATG" + "GCT" * 38 + "TAA"


class TestIndelsAcrossOverlappingHsps(unittest.TestCase):
    """PL-06 and PL-07: an indel seen by two contigs is judged by both of them.

    Adjacent Logan contigs overlap by k-1 = 30 nt, which is exactly the window in which an
    indel can sit 10 bases from one contig's end and deep inside the next.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def tile(self, ref, rows):
        """`rows` are (pident, qstart, qend, qseq, sseq) exactly as the table carries them."""
        hsps = [(pid, 100.0 * (qe - qs + 1) / len(ref), qs, qe, q, s)
                for pid, qs, qe, q, s in rows]
        return lms.tile_accession_hsps(hsps, len(ref), 60.0, 0.80)

    def run_rows(self, name, ref, rows, **kw):
        ref_fa = self.dir / f"{name}.fasta"
        ref_fa.write_text(f">{name}\n{ref}\n")
        tsv = self.dir / f"{name}.tsv"
        with tsv.open("w") as fh:
            fh.write("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n")
            for pid, qs, qe, q, s in rows:
                fh.write(f"SRR1\t{pid:.2f}\t{100.0*(qe-qs+1)/len(ref):.2f}\t{qs}\t{qe}\t{q}\t{s}\n")
        prefix = str(self.dir / name)
        lms.stream_lexicmap_to_msa(ref_path=str(ref_fa), file_path=str(tsv),
                                   out_prefix=prefix, **kw)
        return json.loads(Path(f"{prefix}.summary.json").read_text())["records"][0]

    def deletion_in_both_contigs(self):
        """One 1-nt deletion at reference 71, carried by two overlapping contigs."""
        a_q, a_s = CDS120[:75], CDS120[:70] + "-" + CDS120[71:75]
        b_q, b_s = CDS120[45:], CDS120[45:70] + "-" + CDS120[71:]
        return [(98.67, 1, 75, a_q, a_s), (98.67, 46, 120, b_q, b_s)]

    def test_a_real_deletion_is_not_an_edge_artefact_in_any_hsp_order(self):
        # 71del is 4 bases from the first contig's end and 25 inside the second. Judged from
        # the first contig alone it is "(edge)", so the frame rule ignores it and a
        # frameshifted gene enters the clean cohort; which contig is consulted used to depend
        # on input order and on whose pident was higher.
        a, b = self.deletion_in_both_contigs()
        for n, rows in enumerate(([a, b], [b, a],
                                  [(99.5,) + a[1:], b], [a, (99.5,) + b[1:]])):
            with self.subTest(order=n):
                rec = self.run_rows(f"order{n}", CDS120, rows)
                self.assertEqual(rec["indels"], "71del", "trusted, not (edge)")
                self.assertEqual(rec["status"], "FLAGGED_FRAMESHIFT")

    def test_one_contig_alone_still_sees_it_as_an_edge_artefact(self):
        # The control: with nothing to corroborate it, a deletion 4 bases from the only
        # contig's end stays untrusted, which is the behaviour INDEL_EDGE_WINDOW exists for.
        a, _b = self.deletion_in_both_contigs()
        rec = self.run_rows("alone", CDS120, [a])
        self.assertEqual(rec["indels"], "71del(edge)")
        self.assertEqual(rec["status"], "CLEAN_PARTIAL", "74 of 120 positions tiled")

    def test_an_edge_deletion_does_not_hide_a_base_another_contig_supplied(self):
        # PL-07's fixture: the higher-identity contig closes with a forced 1-nt gap inside
        # the overlap, which is the C. auris case. The gap entered the per-position contest
        # like a base and stayed even after being ruled an artefact, so the row showed '-',
        # coverage was 59/60, and the accession became a haplotype of its own.
        a = (99.0, 1, 60, CDS60, CDS60[:54] + "-" + CDS60[55:])
        b = (98.0, 31, 60, CDS60[30:], CDS60[30:])
        for n, rows in enumerate(([a, b], [b, a])):
            with self.subTest(order=n):
                tiling = self.tile(CDS60, rows)
                self.assertEqual(tiling.seq, CDS60, "the base B supplied is in the row")
                self.assertEqual(tiling.covered, 60)
                self.assertEqual(tiling.indels, [], "the artefact is gone, not recorded")

    def test_that_contig_alone_keeps_the_edge_deletion(self):
        a = (99.0, 1, 60, CDS60, CDS60[:54] + "-" + CDS60[55:])
        tiling = self.tile(CDS60, [a])
        self.assertEqual([e.notation() for e in tiling.indels], ["55del(edge)"])
        self.assertEqual(tiling.covered, 59)

    def test_a_trusted_deletion_still_wins_over_a_lower_depth_base(self):
        # The guard on the refill: only an *edge* deletion is put back. A deletion the tiler
        # trusts is a statement about the sample and must stay in the row.
        a = (99.0, 1, 60, CDS60, CDS60[:30] + "---" + CDS60[33:])
        b = (98.0, 25, 40, CDS60[24:40], CDS60[24:40])
        tiling = self.tile(CDS60, [a, b])
        self.assertEqual([e.notation() for e in tiling.indels], ["31_33del"])
        self.assertEqual(tiling.seq[30:33], "---")


class TestGeneticCodeAndReferenceChecks(unittest.TestCase):
    """PL-04, PL-05, PL-09: what the reference is, and which code reads it."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def run_hits(self, name, ref, samples, **kw):
        ref_fa = self.dir / f"{name}.fasta"
        ref_fa.write_text(f">{name}\n{ref}\n")
        tsv = self.dir / f"{name}.tsv"
        with tsv.open("w") as fh:
            fh.write("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n")
            for i, (qs, qe, qseq, sseq) in enumerate(samples):
                pid = 100.0 * sum(a == b for a, b in zip(qseq, sseq)) / len(qseq)
                fh.write(f"SRR{i}\t{pid:.2f}\t{100.0*(qe-qs+1)/len(ref):.2f}"
                         f"\t{qs}\t{qe}\t{qseq}\t{sseq}\n")
        prefix = str(self.dir / name)
        s = lms.stream_lexicmap_to_msa(ref_path=str(ref_fa), file_path=str(tsv),
                                       out_prefix=prefix, **kw)
        table = (self.dir / f"{name}.clean.haplotypes.tsv")
        rows = [r.split("\t") for r in table.read_text().splitlines()[1:]] if table.is_file() else []
        return json.loads(Path(f"{prefix}.summary.json").read_text()), rows

    def test_a_reference_with_its_own_stops_is_refused(self):
        # Before: the reference's own stops were counted against every sample, so the whole
        # cohort landed in FLAGGED_PREMATURE_STOPS and the run exited 1 saying the clean
        # cohort was empty. Right exit code, wrong cause, and nothing pointed at the input.
        ref = "ATGTAGGAATAA"   # TAG at codon 2
        with self.assertRaisesRegex(ValueError, "in-frame stop codon"):
            self.run_hits("refstop", ref, [(1, 12, ref, ref)])

    def test_non_coding_switches_the_coding_rules_off(self):
        ref = "ATGTAGGAATAA"
        s, _rows = self.run_hits("nc", ref, [(1, 12, ref, ref)], non_coding=True)
        self.assertEqual(s["cohort"]["clean_full_length"], 1)
        self.assertFalse(s["settings"]["reference_is_cds"])

    def test_a_reference_in_another_genetic_code(self):
        # TGA is Trp in code 4: under code 11 this reference cannot be used at all.
        ref = "ATGTGAGCTTAA"
        with self.assertRaisesRegex(ValueError, "genetic code 11"):
            self.run_hits("code11", ref, [(1, 12, ref, ref)])
        s, _rows = self.run_hits("code4", ref, [(1, 12, ref, ref)], genetic_code=4)
        self.assertEqual(s["cohort"]["clean_full_length"], 1)
        self.assertEqual(s["settings"]["genetic_code"], 4)

    def test_the_ctg_clade_code_changes_the_residue_called(self):
        # The shipped candida_auris_erg11 marker is a CTG-clade yeast. CTG is Leu in code 11
        # and Ser in code 12, so the same change is Synonymous under one and Missense under
        # the other, and the residue label itself was wrong.
        ref = "ATGTTGGAATAA"            # codon 2 TTG = Leu in both codes
        alt = "ATGCTGGAATAA"            # CTG: Leu in 11, Ser in 12
        _s, rows = self.run_hits("c11", ref, [(1, 12, ref, alt)])
        self.assertEqual(rows[0][10], "Synonymous")
        _s, rows = self.run_hits("c12", ref, [(1, 12, ref, alt)], genetic_code=12)
        self.assertEqual(rows[0][10], "Missense")
        self.assertEqual(rows[0][9], "L2S")

    def test_an_unknown_genetic_code_names_the_ids_that_exist(self):
        with self.assertRaisesRegex(ValueError, "unknown genetic code 99"):
            lms.codon_table_for(99)

    def test_an_ambiguity_code_is_not_a_mutation(self):
        # `N` used to be reported as a substitution and, through codon_table.get(codon, 'X'),
        # as an X residue classed Missense.
        ref = "ATGGCTGAATAA"
        _s, rows = self.run_hits("amb", ref, [(1, 12, ref, "ATGGNTGAATAA")])
        self.assertEqual(rows[0][6], "0", f"nt mutations: {rows[0][7]}")
        self.assertEqual(rows[0][8], "0", f"aa mutations: {rows[0][9]}")

    def test_the_last_codon_is_scanned_when_the_reference_has_no_terminal_stop(self):
        # Six shipped references are polyprotein segments ending in a sense codon
        # (zika_ns5 ends CTG). A nonsense change in the last codon passed the tier rule while
        # the haplotype table called the same haplotype Nonsense.
        self.assertEqual(lms.find_internal_stops("ATGGCTGAATAA", 12, skip_last=False), [4])
        self.assertEqual(lms.find_internal_stops("ATGGCTGAATAA", 12), [])

        ref = "ATGGCTGAACTG"            # ends CTG, no terminal stop
        s, _rows = self.run_hits("lastcodon", ref, [(1, 12, ref, "ATGGCTGAATAA")])
        rec = s["records"][0]
        self.assertEqual(rec["status"], "FLAGGED_PREMATURE_STOPS")
        self.assertEqual(rec["stop_codons"], "4")
        self.assertFalse(s["settings"]["reference_terminal_stop"])

    def test_a_terminal_stop_is_still_not_a_premature_stop(self):
        ref = "ATGGCTGAATAA"
        s, _rows = self.run_hits("termstop", ref, [(1, 12, ref, ref)])
        self.assertEqual(s["records"][0]["status"], "CLEAN_FULL_LENGTH")
        self.assertTrue(s["settings"]["reference_terminal_stop"])

    def test_a_non_coding_reference_is_not_frameshifted(self):
        # PL-09: for an rRNA, ITS or UTR marker the stop scan was already off and the
        # NonCoding_* classes already existed, but any accession with a 1 nt indel lasting 30
        # bases was still filed FLAGGED_FRAMESHIFT, so that allele could never be clean.
        ref = "".join("ACGT"[(i * 7) % 4] for i in range(100))   # 100 nt, not a multiple of 3
        qseq, sseq = ref, ref[:40] + "-" + ref[41:]
        s, rows = self.run_hits("nc100", ref, [(1, 100, qseq, sseq)])
        rec = s["records"][0]
        self.assertEqual(rec["status"], "CLEAN_FULL_LENGTH")
        self.assertEqual(rec["indels"], "41del")
        self.assertEqual(s["discards"]["accessions"]["flagged_frameshift"], 0)
        self.assertEqual(rows[0][10], "NonCoding_Variant")


class TestReferenceWhitespaceAndNaming(unittest.TestCase):
    """PL-08 and PL-10: what a reference file may contain, and what its id may name."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_whitespace_inside_a_sequence_line_is_not_a_reference_position(self):
        # `validators.fasta_record_lengths` and Biopython both drop it, so one file had two
        # lengths in a single from-lexicmap run. Here it made a 12 nt CDS 13 nt, which turned
        # the coding rules off and let a premature stop into the clean cohort at exit 0.
        for text in (">g\nATGGCT GAATAA\n", ">g\nATGGCT\tGAATAA\n", ">g\nATGGCT\nGAATAA\n"):
            fa = self.dir / "sp.fasta"
            fa.write_text(text)
            name, _header, seq = lms.load_reference(str(fa))
            self.assertEqual(len(seq), 12, repr(text))
            self.assertEqual(seq, "ATGGCTGAATAA")

    def test_text_before_the_first_header_is_refused(self):
        fa = self.dir / "lead.fasta"
        fa.write_text("# a comment\n>g\nATGGCTGAATAA\n")
        with self.assertRaisesRegex(ValueError, "text before the first"):
            lms.load_reference(str(fa))

    def test_a_spaced_reference_still_flags_a_premature_stop(self):
        fa = self.dir / "sp2.fasta"
        fa.write_text(">g\nATGGCT GAATAA\n")
        tsv = self.dir / "sp2.tsv"
        stop = "ATGTAGGAATAA"
        tsv.write_text("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
                       f"SRR1\t91.7\t100.0\t1\t12\tATGGCTGAATAA\t{stop}\n")
        prefix = str(self.dir / "sp2")
        lms.stream_lexicmap_to_msa(ref_path=str(fa), file_path=str(tsv), out_prefix=prefix)
        s = json.loads(Path(f"{prefix}.summary.json").read_text())
        self.assertEqual(s["discards"]["accessions"]["flagged_premature_stops"], 1)

    def names(self):
        return [
            ("../../../victim/lexicmap", "_.._.._victim_lexicmap"),
            ("/abs/dir/x", "_abs_dir_x"),
            ("hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124", "hCoV-19_Wuhan_WIV04_2019|EPI_ISL_402124"),
            ("A/California/07/2009(H1N1)", "A_California_07_2009(H1N1)"),
            ("lcl|NC_045512.2_cds_YP_009724390.1_3", "lcl|NC_045512.2_cds_YP_009724390.1_3"),
        ]

    def test_an_id_that_is_a_path_becomes_one_safe_component(self):
        for raw, want in self.names():
            got = lms.file_safe_name(raw)
            self.assertEqual(got, want, raw)
            # The property the rest of it rests on: one path component, and not a name that
            # resolves to the directory itself or its parent.
            self.assertNotIn("/", got)
            self.assertNotIn("\\", got)
            self.assertNotIn(got, (".", "..", ""))
        self.assertEqual(len(lms.file_safe_name("x" * 300)), 200)
        # Nothing usable left: the caller falls back to the file stem, then to "lexicmap".
        self.assertEqual(lms.file_safe_name("..."), "")
        self.assertEqual(lms.file_safe_name(".."), "")
        # Separators alone still leave a usable, harmless component.
        self.assertEqual(lms.file_safe_name("/"), "_")

    def test_outputs_stay_inside_the_output_directory(self):
        # The crafted id used to write all twelve outputs into a directory elsewhere, and
        # because the cohort files are opened before the table is checked, even a failing run
        # truncated whatever was already there.
        outside = self.dir / "victim"
        outside.mkdir()
        out = self.dir / "sandbox" / "out"
        fa = self.dir / "evil.fasta"
        fa.write_text(">../../victim/lexicmap\nATGGCTGAATAA\n")
        tsv = self.dir / "evil.tsv"
        tsv.write_text("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
                       "SRR1\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n")
        lms.stream_lexicmap_to_msa(ref_path=str(fa), file_path=str(tsv),
                                   out_prefix=str(out) + "/")
        self.assertEqual(list(outside.iterdir()), [], "nothing may be written outside")
        self.assertTrue(any(out.iterdir()))
        self.assertTrue((out / "_.._victim_lexicmap.summary.json").is_file())

    def test_a_gisaid_style_id_no_longer_stops_the_run(self):
        fa = self.dir / "gisaid.fasta"
        fa.write_text(">hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124\nATGGCTGAATAA\n")
        tsv = self.dir / "gisaid.tsv"
        tsv.write_text("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
                       "SRR1\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n")
        out = self.dir / "g"
        s = lms.stream_lexicmap_to_msa(ref_path=str(fa), file_path=str(tsv),
                                       out_prefix=str(out) + "/")
        self.assertEqual(s["cohort"]["clean_full_length"], 1)


class TestQuantilesAndRegime(unittest.TestCase):
    """PL-15 and PL-16: what the median is, and what "no data" is called."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_an_even_cohorts_median_is_between_its_middle_values(self):
        from disassembler import summary as summary_io
        self.assertIsNone(summary_io.quantiles([]))
        self.assertEqual(summary_io.quantiles([0.5])["p50"], 0.5)
        self.assertEqual(summary_io.quantiles([0.55, 0.65])["p50"], 0.6)
        self.assertEqual(summary_io.quantiles([0.5, 0.55, 0.65, 0.7])["p50"], 0.6)
        self.assertEqual(summary_io.quantiles([1, 2, 3, 4, 5, 6])["p50"], 3.5)
        # The old rounded index returned one of the two middle values, which one depending on
        # banker's rounding, so these two cohorts used to straddle the divergent threshold.
        self.assertEqual(summary_io.quantiles([0.55, 0.65])["p25"], 0.575)

    def run_covs(self, name, covs):
        ref = CDS60
        fa = self.dir / f"{name}.fasta"
        fa.write_text(f">{name}\n{ref}\n")
        tsv = self.dir / f"{name}.tsv"
        with tsv.open("w") as fh:
            fh.write("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n")
            for i, cov in enumerate(covs):
                n = int(round(cov * len(ref)))
                fh.write(f"SRR{i}\t100.00\t{cov*100:.2f}\t1\t{n}\t{ref[:n]}\t{ref[:n]}\n")
        prefix = str(self.dir / name)
        lms.stream_lexicmap_to_msa(ref_path=str(fa), file_path=str(tsv), out_prefix=prefix,
                                   min_coverage_partial=0.0)
        return json.loads(Path(f"{prefix}.summary.json").read_text())["cohort"]

    def test_the_regime_call_follows_the_interpolated_median(self):
        c = self.run_covs("two", [0.55, 0.65])
        self.assertEqual(c["cov"]["p50"], 0.6)
        self.assertEqual(c["regime"], "conserved", "0.60 is not below 0.60")

    def test_a_search_that_found_nothing_is_not_divergent(self):
        # A header-only table used to be labelled `divergent` with "median identity is 0.0%",
        # on quantiles invented for an empty list, and under --allow-empty a fan-out over
        # many targets exited 0 carrying that record.
        fa = self.dir / "none.fasta"
        fa.write_text(f">none\n{CDS60}\n")
        tsv = self.dir / "none.tsv"
        tsv.write_text("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n")
        prefix = str(self.dir / "none")
        lms.stream_lexicmap_to_msa(ref_path=str(fa), file_path=str(tsv), out_prefix=prefix)
        c = json.loads(Path(f"{prefix}.summary.json").read_text())["cohort"]
        self.assertEqual(c["regime"], "no_data")
        self.assertIsNone(c["cov"])
        self.assertIsNone(c["pident"])
        self.assertIsNone(c["regime_warning"])

    def test_a_cohort_under_min_pident_says_so(self):
        fa = self.dir / "low.fasta"
        fa.write_text(f">low\n{CDS60}\n")
        alt = "".join("A" if c != "A" else "C" for c in CDS60)
        tsv = self.dir / "low.tsv"
        tsv.write_text("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
                       f"SRR1\t50.00\t100.00\t1\t60\t{CDS60}\t{alt}\n")
        prefix = str(self.dir / "low")
        lms.stream_lexicmap_to_msa(ref_path=str(fa), file_path=str(tsv), out_prefix=prefix)
        c = json.loads(Path(f"{prefix}.summary.json").read_text())["cohort"]
        self.assertEqual(c["regime"], "divergent")
        self.assertIn("--min-pident", c["regime_warning"])
        self.assertIsNone(c["pident"], "no HSP was tiled, so there is no identity to report")


class TestFrameSpanAndCoverageGates(unittest.TestCase):
    """PL-17, PL-18, PL-19: what counts as out of frame, as covered, and as the gate."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    @staticmethod
    def ind(kind, pos, length):
        return lms.Indel(kind, pos, length, "A" * length, False)

    def test_a_deletions_own_bases_are_not_read_out_of_frame(self):
        # The sample does not have those bases, so the frame shifts only after them. Opening
        # the cluster at the deletion's start made `31_64del,70_71del` (5 bases genuinely out
        # of frame) a frameshift while `31_34del,40_41del`, the same 5 bases, stayed clean.
        long_pair = [self.ind("del", 30, 34), self.ind("del", 69, 2)]
        short_pair = [self.ind("del", 30, 4), self.ind("del", 39, 2)]
        for name, pair in (("long", long_pair), ("short", short_pair)):
            with self.subTest(name):
                in_frame, shifted = lms.split_by_frame(pair)
                self.assertEqual(shifted, [], f"{name}: 5 bases is under FRAMESHIFT_MIN_SPAN")
                self.assertEqual(len(in_frame), 2)

    def test_a_stretch_genuinely_out_of_frame_is_still_a_frameshift(self):
        # The control, placed one base past the limit so the boundary itself is pinned.
        pair = [self.ind("del", 100, 34),
                self.ind("del", 100 + 34 + lms.FRAMESHIFT_MIN_SPAN, 2)]
        _in_frame, shifted = lms.split_by_frame(pair)
        self.assertEqual(len(shifted), 2)

    def test_an_insertions_bases_still_count(self):
        # An insertion's bases ARE in the sample and are read in the shifted frame, so only
        # deletions are discounted.
        pair = [self.ind("ins", 100, 1), self.ind("del", 100 + lms.FRAMESHIFT_MIN_SPAN, 2)]
        _in_frame, shifted = lms.split_by_frame(pair)
        self.assertEqual(len(shifted), 2)

    def run_hits(self, name, ref, rows, **kw):
        fa = self.dir / f"{name}.fasta"
        fa.write_text(f">{name}\n{ref}\n")
        tsv = self.dir / f"{name}.tsv"
        with tsv.open("w") as fh:
            fh.write("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n")
            for acc, pid, qs, qe, qseq, sseq in rows:
                fh.write(f"{acc}\t{pid:.2f}\t{100.0*(qe-qs+1)/len(ref):.2f}"
                         f"\t{qs}\t{qe}\t{qseq}\t{sseq}\n")
        prefix = str(self.dir / name)
        lms.stream_lexicmap_to_msa(ref_path=str(fa), file_path=str(tsv),
                                   out_prefix=prefix, **kw)
        return json.loads(Path(f"{prefix}.summary.json").read_text())

    def test_an_accession_with_nothing_tiled_is_a_fragment_at_any_threshold(self):
        # PL-18: with both coverage gates at 0, which the Galaxy wrapper allows, an accession
        # whose HSPs all failed --min-pident became CLEAN_PARTIAL, and its all-gap row became
        # a WT_partial haplotype counted in the frequencies and the entropy.
        alt = "".join("A" if c != "A" else "C" for c in CDS60)
        s = self.run_hits("nothing", CDS60, [("SRR1", 20.0, 1, 60, CDS60, alt)],
                          min_coverage=0.0, min_coverage_partial=0.0)
        rec = s["records"][0]
        self.assertEqual(rec["status"], "LOW_COVERAGE_FRAGMENT")
        self.assertEqual(rec["coverage"], 0.0)
        self.assertEqual(s["cohort"]["clean_full_length"], 0)
        self.assertEqual(s["cohort"]["clean_partial"], 0)
        clean = (self.dir / "nothing.clean_expanded.accessions.fasta").read_text()
        self.assertNotIn("-" * 60, clean, "an all-gap row is not a clean sequence")

    def test_require_full_length_is_the_gate_that_min_cov_partial_is_checked_against(self):
        # PL-19: `--require-full-length --min-cov-partial 0.9` is well formed, because the
        # tier 1 gate is 1.0, but the check compared against --min-cov and refused it.
        from disassembler.stage import stage_args
        import contextlib, io as _io
        ns = stage_args(lms, ref="x", file="y", require_full_length=True,
                        min_cov=0.8, min_cov_partial=0.9)
        err = _io.StringIO()
        with contextlib.redirect_stderr(err):
            rc = lms.run(ns)
        self.assertNotEqual(rc, 2, f"refused a valid combination: {err.getvalue()}")

    def test_the_rescue_count_is_judged_against_the_active_gate(self):
        # `multi_hsp_rescued_accessions` counts accessions no single HSP could have carried.
        # Under --require-full-length the gate is 1.0, but tile_accession_hsps was handed
        # --min-cov, so a single 90% HSP looked like it had cleared the gate on its own.
        rows = [("SRR1", 100.0, 1, 54, CDS60[:54], CDS60[:54]),
                ("SRR1", 100.0, 55, 60, CDS60[54:], CDS60[54:])]
        s = self.run_hits("rescue", CDS60, rows, full_length_only=True)
        self.assertEqual(s["records"][0]["status"], "CLEAN_FULL_LENGTH")
        self.assertEqual(s["cohort"]["multi_hsp_rescued_accessions"], 1,
                         "90% is not full length, so this accession was rescued by tiling")


class TestRowLoopAndOutputs(unittest.TestCase):
    """PL-12, PL-13, PL-14, PL-20, PL-21, PL-22: the stream's edges, and what lands on disk."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.ref = self.dir / "ref.fasta"
        self.ref.write_text(f">cds60\n{CDS60}\n")

    def tearDown(self):
        self.tmp.cleanup()

    HEADER = "sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"

    def row(self, acc, qs, qe, sseq=None):
        q = CDS60[qs - 1:qe]
        s = sseq if sseq is not None else q
        pid = 100.0 * sum(a == b for a, b in zip(q, s)) / len(q)
        return f"{acc}\t{pid:.2f}\t{100.0*(qe-qs+1)/60:.2f}\t{qs}\t{qe}\t{q}\t{s}\n"

    def run_text(self, name, text, **kw):
        tsv = self.dir / f"{name}.tsv"
        tsv.write_text(text)
        prefix = str(self.dir / name)
        lms.stream_lexicmap_to_msa(ref_path=str(self.ref), file_path=str(tsv),
                                   out_prefix=prefix, **kw)
        return json.loads(Path(f"{prefix}.summary.json").read_text())

    def test_max_records_finishes_the_accession_it_is_in_the_middle_of(self):
        # PL-13: the limit was tested after a row joined the open block and the partial block
        # was then tiered, so the last accession of a preview was judged from a fraction of
        # its HSPs: a clean three-HSP accession came out LOW_COVERAGE_FRAGMENT. Galaxy
        # exposes this as a "preview" knob, where that is the wrong thing to show.
        text = self.HEADER + "".join(
            [self.row("SRR2", 1, 20), self.row("SRR2", 21, 40), self.row("SRR2", 41, 60),
             self.row("SRR3", 1, 60)])
        s = self.run_text("preview", text, max_records=2, min_coverage_partial=0.0)
        self.assertEqual(len(s["records"]), 1, "the next accession is not opened")
        rec = s["records"][0]
        self.assertEqual(rec["acc"], "SRR2")
        self.assertEqual(rec["hsps_merged"], 3, "its block was finished")
        self.assertEqual(rec["status"], "CLEAN_FULL_LENGTH")

    def test_a_short_row_is_counted_and_named(self):
        # PL-14: a silent `continue`, with no counter and no message, unlike every other
        # skip. A table cut inside its last row demoted its final accession at exit 0.
        text = self.HEADER + self.row("SRR1", 1, 60) + "SRR2\t100.0\t50.0\n"
        s = self.run_text("short", text)
        self.assertEqual(s["discards"]["rows"]["too_few_columns"], 1)
        self.assertEqual(len(s["records"]), 1)

    def test_a_cut_last_row_is_refused(self):
        text = self.HEADER + self.row("SRR1", 1, 60) + "SRR2\t100.0\t50.0\t1"
        with self.assertRaisesRegex(ValueError, "looks truncated"):
            self.run_text("cut", text)

    def test_a_row_whose_qseq_is_not_the_reference_is_skipped(self):
        # PL-20: nothing compared qseq with the reference slice it is projected onto, so a
        # reference carrying the search's query name and length but a different sequence
        # tiled every accession against the wrong sequence at exit 0.
        wrong = "".join("A" if c != "A" else "C" for c in CDS60[:30])
        text = self.HEADER + self.row("SRR1", 1, 60) + \
            f"SRR2\t100.00\t50.00\t1\t30\t{wrong}\t{wrong}\n"
        s = self.run_text("wrongq", text)
        self.assertEqual(s["discards"]["rows"]["inconsistent_qseq"], 1)
        self.assertEqual(len(s["records"]), 1)

    def test_a_whole_table_against_the_wrong_reference_is_refused(self):
        wrong = "".join("A" if c != "A" else "C" for c in CDS60)
        text = self.HEADER + f"SRR1\t100.00\t100.00\t1\t60\t{wrong}\t{wrong}\n"
        with self.assertRaisesRegex(ValueError, "wrong reference for this search"):
            self.run_text("wrongref", text)

    def test_a_table_without_the_alignment_columns_names_the_flag_to_use(self):
        # PL-21: a KeyError traceback for the likeliest mistake there is, running
        # `lexicmap search` without -a/--all.
        text = "sgenome\tpident\tqcovHSP\tqstart\tqend\nSRR1\t100.0\t100.0\t1\t60\n"
        with self.assertRaises(ValueError) as cm:
            self.run_text("noall", text)
        self.assertIn("-a/--all", str(cm.exception))
        self.assertIn("qseq", str(cm.exception))

    def test_the_sort_advice_keeps_the_header_where_it_is(self):
        # PL-12: the advice used to be a bare `sort`, which sorts the header in with the
        # data; under LC_ALL=C the header goes last and the re-run dies in the schema check.
        text = self.HEADER + self.row("SRR1", 1, 30) + self.row("SRR2", 1, 60) + \
            self.row("SRR1", 31, 60)
        with self.assertRaises(ValueError) as cm:
            self.run_text("ungrouped", text)
        msg = str(cm.exception)
        self.assertIn("head -n1", msg)
        self.assertIn("tail -n +2", msg)
        self.assertIn("-k1,1", msg, "the real column, not a placeholder")

    def test_a_failed_run_leaves_nothing_half_written(self):
        # PL-22: the cohort files were opened under their final names before a row was read.
        text = self.HEADER + self.row("SRR1", 1, 30) + self.row("SRR2", 1, 60) + \
            self.row("SRR1", 31, 60)
        with self.assertRaises(ValueError):
            self.run_text("part", text)
        leftovers = sorted(f.name for f in self.dir.iterdir() if f.name.startswith("part"))
        self.assertEqual([f for f in leftovers if f.endswith(".part")], [])
        self.assertEqual([f for f in leftovers if f.endswith(".fasta")], [],
                         f"a failed run published cohort files: {leftovers}")

    def test_a_rerun_without_a2m_does_not_leave_the_old_one(self):
        text = self.HEADER + self.row("SRR1", 1, 60)
        self.run_text("rerun", text)
        a2m = self.dir / "rerun.clean.a2m"
        self.assertTrue(a2m.is_file())
        first = a2m.read_text()
        self.run_text("rerun", text + self.row("SRR2", 1, 60, CDS60[:30] + CDS60[30:]),
                      write_a2m=False)
        self.assertFalse(a2m.is_file(),
                         f"the previous run's A2M survived a --no-a2m rerun:\n{first[:80]}")


class TestStreamHandling(unittest.TestCase):
    """PL-21 and PL-26: a damaged gzip, and whose stdin it is."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.ref = self.dir / "ref.fasta"
        self.ref.write_text(f">cds60\n{CDS60}\n")
        self.table = ("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
                      f"SRR1\t100.00\t100.00\t1\t60\t{CDS60}\t{CDS60}\n")

    def tearDown(self):
        self.tmp.cleanup()

    def test_a_truncated_gzip_names_the_file_instead_of_raising_eoferror(self):
        # gzip.GzipFile raises EOFError for a truncated member and zlib.error for a corrupt
        # one, and neither is caught by the `except (OSError, ValueError)` that lexicmap,
        # msa, sieve, stratify, reconstruct and from-lexicmap all use (PL-21).
        from disassembler.textio import open_text
        blob = gzip.compress(self.table.encode())
        cut = self.dir / "cut.tsv.gz"
        cut.write_bytes(blob[: len(blob) // 2])
        with self.assertRaisesRegex(ValueError, "truncated or corrupt gzip"):
            open_text(cut).read()

    def test_a_corrupt_gzip_names_the_file_too(self):
        from disassembler.textio import open_text
        blob = bytearray(gzip.compress(self.table.encode()))
        blob[len(blob) // 2] ^= 0xFF
        bad = self.dir / "bad.tsv.gz"
        bad.write_bytes(bytes(blob))
        with self.assertRaisesRegex(ValueError, "truncated or corrupt gzip"):
            open_text(bad).read()

    def test_a_whole_gzip_still_reads(self):
        from disassembler.textio import open_text
        ok = self.dir / "ok.tsv.gz"
        ok.write_bytes(gzip.compress(self.table.encode()))
        self.assertEqual(open_text(ok).read(), self.table)

    def test_reading_from_stdin_leaves_the_process_stdin_open(self):
        # PL-26: the TextIOWrapper around sys.stdin.buffer was neither closed nor detached,
        # so IOBase.__del__ closed the process's stdin when it was collected and a library
        # caller reading stdin afterwards got "I/O operation on closed file". The CLI never
        # noticed, because fd 0 stays open either way.
        import contextlib, io as _io
        real = sys.stdin
        buf = _io.BufferedReader(_io.BytesIO(self.table.encode()))
        sys.stdin = _io.TextIOWrapper(buf, encoding="utf-8")
        try:
            with contextlib.redirect_stdout(_io.StringIO()):
                lms.stream_lexicmap_to_msa(ref_path=str(self.ref), file_path="-",
                                           out_prefix=str(self.dir / "in"))
            import gc
            gc.collect()
            self.assertFalse(sys.stdin.buffer.closed, "the caller's stdin was closed")
        finally:
            sys.stdin = real


# 126 nt, MKPREFNETKPREFNET...*, with a repeating motif so an indel's effect is easy to read.
CDS126 = "ATG" + "AAACCTCGTGAATTTAATGAAACT" * 5 + "TAA"


class TestIndelProteinChanges(unittest.TestCase):
    """PL-29: a haplotype whose only change is an in-frame indel used to report none.

    Amino-acid changes were called on the reference-GRID row, where an insertion is absent
    and deleted bases are `-`, and `translate_nt` returns `-` for any codon holding a gap. So
    every codon an indel touched was skipped, inserted residues never appeared, and the table
    printed `num_aa_muts 0` and `aa_mutations NA` for a haplotype that `nt_mutations` and the
    `InFrameIndel` class both said had changed.

    Each expectation below is checked against the protein Biopython gets from the sample's
    own spliced sequence, so the test does not simply restate the implementation.
    """

    REF_AA = None

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.ref_aa = lms.translate_nt(CDS126)

    def tearDown(self):
        self.tmp.cleanup()

    @staticmethod
    def protein_of(row, indels):
        from Bio.Seq import Seq
        sample, _s, _p = lms.apply_indels(row, len(CDS126), indels)
        return str(Seq(sample.replace("N", "")).translate(table=11))

    def call(self, indels, row):
        changes, covered = lms.indel_protein_changes(CDS126, row, self.ref_aa, indels)
        return changes, covered, self.protein_of(row, indels)

    def test_a_codon_aligned_deletion_is_one_residue_deleted(self):
        d = lms.Indel("del", 30, 3, CDS126[30:33], False)
        changes, _cov, protein = self.call([d], CDS126[:30] + "---" + CDS126[33:])
        self.assertEqual(changes, ["P11del"])
        self.assertEqual(protein, self.ref_aa[:10] + self.ref_aa[11:],
                         "the reference protein minus residue 11")

    def test_a_codon_aligned_insertion_is_one_residue_inserted(self):
        i = lms.Indel("ins", 30, 3, "GCT", False)
        changes, _cov, protein = self.call([i], CDS126)
        self.assertEqual(changes, ["K10_P11insA"])
        self.assertEqual(protein, self.ref_aa[:10] + "A" + self.ref_aa[10:])

    def test_two_distant_indels_are_two_changes(self):
        # The window is computed per cluster; one global prefix/suffix comparison would merge
        # these into a single delins spanning 20 residues.
        far = [lms.Indel("del", 30, 3, CDS126[30:33], False),
               lms.Indel("del", 90, 3, CDS126[90:93], False)]
        row = CDS126[:30] + "---" + CDS126[33:90] + "---" + CDS126[93:]
        changes, _cov, _p = self.call(far, row)
        self.assertEqual(changes, ["P11del", "N31del"])

    def test_a_pair_that_restores_the_frame_reports_the_stretch_between_them(self):
        # The code deliberately keeps a +1/-1 pair inside FRAMESHIFT_MIN_SPAN in the clean
        # tiers. The residues between them differ in the real protein, and the grid shows
        # them as reference, so up to ten substitutions were hidden.
        pair = [lms.Indel("ins", 40, 1, "G", False),
                lms.Indel("del", 52, 1, CDS126[52], False)]
        row = CDS126[:52] + "-" + CDS126[53:]
        changes, covered, protein = self.call(pair, row)
        self.assertEqual(len(changes), 1)
        self.assertIn("delins", changes[0])
        self.assertTrue(covered, "the codons it describes are marked, so they are not "
                                 "reported a second time from the grid")
        self.assertNotEqual(protein[13:18], self.ref_aa[13:18])

    def test_an_edge_indel_is_not_described(self):
        # Untrusted indels are left as the row shows them, here as everywhere else.
        d = lms.Indel("del", 30, 3, CDS126[30:33], True)
        changes, covered, _p = self.call([d], CDS126[:30] + "---" + CDS126[33:])
        self.assertEqual((changes, covered), ([], set()))

    def test_a_frameshift_is_not_described_as_an_in_frame_change(self):
        d = lms.Indel("del", 30, 1, CDS126[30], False)
        changes, _cov, _p = self.call([d], CDS126[:30] + "-" + CDS126[31:])
        self.assertEqual(changes, [], "a lasting frameshift is reported by its own tier")

    def test_the_haplotype_table_carries_the_indels_protein_change(self):
        fa = self.dir / "cds126.fasta"
        fa.write_text(f">cds126\n{CDS126}\n")
        row = CDS126[:30] + "---" + CDS126[33:]
        tsv = self.dir / "cds126.tsv"
        tsv.write_text("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n"
                       f"SRR1\t97.62\t100.00\t1\t126\t{CDS126}\t{row}\n")
        prefix = str(self.dir / "cds126")
        lms.stream_lexicmap_to_msa(ref_path=str(fa), file_path=str(tsv), out_prefix=prefix)
        rows = (self.dir / "cds126.clean.haplotypes.tsv").read_text().splitlines()
        cols = rows[0].split("\t")
        rec = dict(zip(cols, rows[1].split("\t")))
        self.assertEqual(rec["mutation_class"], "InFrameIndel")
        self.assertEqual(rec["aa_mutations"], "P11del", "was NA, with num_aa_muts 0")
        self.assertEqual(rec["num_aa_muts"], "1")


class TestSharedHelpers(unittest.TestCase):
    """PL-24 and PL-27: helpers that existed twice, or in a form that could not scale."""

    def test_counted_quantiles_agree_with_the_list_based_ones(self):
        # PL-24 reports coverage quantiles from a Counter rather than one float per
        # accession. The two have to be the same number: `cohort.cov` is compared between
        # runs and across cohorts.
        import random
        from collections import Counter
        from disassembler import summary as summary_io
        random.seed(3)
        for n in (1, 2, 3, 4, 7, 50, 999):
            vals = [random.randrange(0, 601) for _ in range(n)]
            want = summary_io.quantiles([v / 600 for v in vals])
            got = summary_io.quantiles_from_counts(Counter(vals), 1.0 / 600)
            self.assertEqual(got, want, f"n={n}")
        self.assertIsNone(summary_io.quantiles_from_counts(Counter()))

    def test_one_galaxy_api_key_implementation(self):
        # PL-27: `get_galaxy_api_key` existed in both harvest and lexicmap and had drifted.
        # harvest's copy took no api_key argument, so --api-key was honoured by one command
        # and quietly ignored by the other.
        from disassembler.cmd import harvest as hv
        from disassembler import galaxy
        self.assertIs(lms.get_galaxy_api_key.__wrapped__
                      if hasattr(lms.get_galaxy_api_key, "__wrapped__")
                      else galaxy.api_key, galaxy.api_key)
        with tempfile.TemporaryDirectory() as d:
            kf = Path(d) / "key.txt"
            kf.write_text("from-the-file\n")
            env = os.environ.get("GALAXY_API_KEY")
            os.environ["GALAXY_API_KEY"] = "from-the-env"
            try:
                # The named file beats the environment, in both commands.
                self.assertEqual(lms.get_galaxy_api_key(None, str(kf)), "from-the-file")
                self.assertEqual(hv.get_galaxy_api_key(str(kf)), "from-the-file")
                # And an explicit key beats everything, in both.
                self.assertEqual(lms.get_galaxy_api_key("literal", str(kf)), "literal")
                self.assertEqual(hv.get_galaxy_api_key(str(kf), "literal"), "literal")
            finally:
                if env is None:
                    del os.environ["GALAXY_API_KEY"]
                else:
                    os.environ["GALAXY_API_KEY"] = env

    def test_one_message_for_a_missing_binary(self):
        # PL-27: four hand-written strings. `sieve` said "put it on PATH" without naming the
        # file it looks for, so a user who installed something called
        # `kmer-diversity-sieve` still failed, and the composites exited 1 where the stages
        # exited 127.
        from disassembler.resolve import not_found_message, EX_TOOL_NOT_FOUND
        msg = not_found_message("sieve", "from-lexicmap")
        self.assertTrue(msg.startswith("from-lexicmap: "))
        self.assertIn("l1_pipeline_serial", msg)
        self.assertIn("DISASSEMBLER_SIEVE_BIN", msg)
        self.assertIn("--sieve-bin", msg)
        self.assertIn("make serial", msg)
        self.assertEqual(EX_TOOL_NOT_FOUND, 127)
        walker = not_found_message("walker", "reconstruct")
        self.assertIn("logan-walker", walker)
        self.assertIn("DISASSEMBLER_WALKER_BIN", walker)

    def test_a_missing_sieve_exits_127_from_the_stage_and_the_composite(self):
        # The composite exited 1 where the stage it wraps exits 127, so a caller could not
        # tell "the tool is missing" from "the data was bad" by status alone (PL-27).
        # `resolve_sieve` is patched rather than pointed at a bad path, because it falls
        # through an unusable --sieve-bin to PATH and then to the checkout, and this machine
        # has the binary built.
        import contextlib, io as _io
        from disassembler.stage import stage_args
        from disassembler.cmd import sieve as sv
        from disassembler.cmd import pipe_from_lexicmap as pfl
        with tempfile.TemporaryDirectory() as d:
            fa = Path(d) / "q.fasta"
            fa.write_text(f">q\n{CDS60}\n")
            tsv = Path(d) / "hits.tsv"
            tsv.write_text("sgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n")
            real = (sv.resolve_sieve, pfl.resolve_sieve)
            sv.resolve_sieve = pfl.resolve_sieve = lambda *a, **k: None
            try:
                for mod, ns, name in (
                    (sv, stage_args(sv, input=str(fa), output=str(Path(d) / "o.fa"),
                                    min_len=10), "sieve"),
                    (pfl, stage_args(pfl, ref=str(fa), file=str(tsv),
                                     out_dir=str(Path(d) / "out"), sieve=True),
                     "from-lexicmap"),
                ):
                    err = _io.StringIO()
                    with contextlib.redirect_stdout(_io.StringIO()), \
                         contextlib.redirect_stderr(err):
                        rc = mod.run(ns)
                    self.assertEqual(rc, 127, f"{name}: {err.getvalue()}")
                    self.assertIn("l1_pipeline_serial", err.getvalue(), name)
                    self.assertIn("DISASSEMBLER_SIEVE_BIN", err.getvalue(), name)
            finally:
                sv.resolve_sieve, pfl.resolve_sieve = real


if __name__ == "__main__":
    unittest.main()
