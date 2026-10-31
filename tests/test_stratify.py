#!/usr/bin/env python3
"""test_stratify.py -- unit tests for the cohort stratifier.

--score-bin: arbitrary half-open [LO, HI) ranges (a bin ending at 1.0 includes it), no
backfill, nonzero scores only. Score files are found by content and name their own target.
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

try:
    from disassembler.cmd import stratify as hk
except ImportError as e:  # pragma: no cover - environment-dependent
    print(f"SKIP: disassembler.cmd.stratify not importable ({e})")
    sys.exit(0)


class TestStratifiedSample(unittest.TestCase):
    """A zero quota is legitimate: int(total_target * 0.20) is 0 for --diversity under 5."""

    def test_a_zero_quota_returns_nothing_instead_of_dividing_by_zero(self):
        items = [("a", 0.99), ("b", 0.98)]
        self.assertEqual(hk.stratified_sample(items, 0), [])
        self.assertEqual(hk.stratified_sample(items, -1), [])

    def test_normal_quotas_are_unchanged(self):
        items = [(c, 1.0) for c in "abcdefgh"]
        self.assertEqual(len(hk.stratified_sample(items, 4)), 4)
        self.assertEqual(hk.stratified_sample(items, 99), items)


class TestBinStratification(unittest.TestCase):
    """--score-bin gives the caller arbitrary [LO, HI) ranges with per-range counts, no
    backfill, and only nonzero-score accessions.
    """

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.out = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def _write_scores(self, scores):
        (self.out / "t.scores.json").write_text(json.dumps({"target": "t", "scores": scores}))

    def _stratify(self, bins):
        hk.stratify_bins(self.out, bins, hk.find_score_files([self.out]))

    def _txt(self):
        return [x for x in (self.out / "t.bins.txt").read_text().splitlines() if x]

    def _summary(self):
        return json.loads((self.out / "t.bins.json").read_text())

    def test_parse_bin_spec_roundtrips_and_rejects_bad_input(self):
        self.assertEqual(hk.parse_bin_spec("0.0:0.05:10"), (0.0, 0.05, 10))
        for bad in ("0.0:0.05", "a:b:c", "0.05:0.0:10", "0.0:0.05:-1"):
            with self.assertRaises(ValueError, msg=f"'{bad}' should be rejected"):
                hk.parse_bin_spec(bad)

    def test_assert_bins_disjoint_flags_overlap(self):
        hk.assert_bins_disjoint([(0.0, 0.05, 1), (0.05, 0.10, 1)])  # boundary touch is OK
        with self.assertRaises(ValueError):
            hk.assert_bins_disjoint([(0.0, 0.10, 1), (0.05, 0.15, 1)])

    def test_each_bin_gets_up_to_its_requested_count(self):
        # 100 accessions in [0.0, 0.05), 100 in [0.05, 0.10), 100 in [0.10, 0.15).
        scores = {}
        for base, tag in ((0.0, "a"), (0.05, "b"), (0.10, "c")):
            for i in range(100):
                scores[f"{tag}{i:03d}"] = base + 1e-9 + i * (0.04999 / 100)
        self._write_scores(scores)
        bins = [(0.0, 0.05, 10), (0.05, 0.10, 10), (0.10, 0.15, 10)]
        self._stratify(bins)
        got = self._txt()
        self.assertEqual(len(got), 30, f"expected 30 picks, got {len(got)}")
        buckets = {"a": 0, "b": 0, "c": 0}
        for acc in got:
            buckets[acc[0]] += 1
        self.assertEqual(buckets, {"a": 10, "b": 10, "c": 10})

    def test_zero_score_accessions_are_excluded(self):
        # Ten zeros and ten non-zeros in the same bin range: only the non-zeros are candidates.
        scores = {f"z{i:02d}": 0.0 for i in range(10)}
        scores.update({f"n{i:02d}": 0.02 + i * 1e-4 for i in range(10)})
        self._write_scores(scores)
        self._stratify([(0.0, 0.05, 15)])
        got = self._txt()
        self.assertEqual(len(got), 10, "should get all 10 nonzeros, no zeros")
        self.assertTrue(all(a.startswith("n") for a in got))
        s = self._summary()
        self.assertEqual(s["n_available_nonzero"], 10)
        self.assertEqual(s["n_total_input"], 20)

    def test_no_backfill_short_bin_contributes_what_it_has(self):
        # First bin: 3 candidates for a request of 10. Others plentiful.
        scores = {f"a{i}": 0.01 + i * 1e-3 for i in range(3)}
        scores.update({f"b{i:02d}": 0.06 + i * 1e-4 for i in range(50)})
        self._write_scores(scores)
        self._stratify([(0.0, 0.05, 10), (0.05, 0.10, 10)])
        s = self._summary()
        self.assertEqual(s["bins"][0]["selected"], 3, "short bin gives what it has")
        self.assertEqual(s["bins"][0]["available"], 3)
        self.assertEqual(s["bins"][1]["selected"], 10, "the other bin is NOT topped up above 10")
        self.assertEqual(s["total_selected"], 13)

    def test_output_txt_is_sorted_by_descending_score(self):
        scores = {f"x{i:03d}": 0.01 + i * 1e-3 for i in range(100)}
        self._write_scores(scores)
        self._stratify([(0.0, 0.15, 5)])
        txt = self._txt()
        s = self._summary()
        score_map = {}
        for b in s["bins"]:
            for row in b["accessions"]:
                score_map[row["acc"]] = row["score"]
        ordered = [score_map[a] for a in txt]
        self.assertEqual(ordered, sorted(ordered, reverse=True))

    def test_half_open_boundary_hi_is_exclusive(self):
        scores = {"lo": 0.05, "mid": 0.075, "at_hi": 0.10, "above": 0.12}
        self._write_scores(scores)
        # Bin [0.05, 0.10) must include lo and mid, exclude at_hi and above.
        self._stratify([(0.05, 0.10, 5)])
        self.assertEqual(sorted(self._txt()), ["lo", "mid"])

    def test_a_bin_ending_at_one_includes_a_perfect_score(self):
        # `0.9:1.0` used to miss every score of exactly 1.0, the strongest hits.
        self._write_scores({"perfect": 1.0, "high": 0.95, "edge": 0.9, "low": 0.5})
        self._stratify([(0.9, 1.0, 10)])
        self.assertEqual(sorted(self._txt()), ["edge", "high", "perfect"])
        self.assertTrue(hk.in_bin(1.0, 0.9, 1.01), "an upper bound above 1 still works")
        self.assertFalse(hk.in_bin(0.5, 0.0, 0.5), "below 1.0 bins stay half-open")

    def test_a_bin_ending_at_one_overlaps_one_starting_there(self):
        with self.assertRaises(ValueError):
            hk.assert_bins_disjoint([(0.9, 1.0, 1), (1.0, 1.01, 1)])
        hk.assert_bins_disjoint([(0.5, 0.9, 1), (0.9, 1.0, 1)])

    def test_multi_target_writes_per_target_files(self):
        (self.out / "geneA.scores.json").write_text(json.dumps(
            {"target": "geneA", "scores": {f"a{i}": 0.02 + i * 1e-4 for i in range(20)}}))
        (self.out / "geneB.scores.json").write_text(json.dumps(
            {"target": "geneB", "scores": {f"b{i}": 0.07 + i * 1e-4 for i in range(20)}}))
        self._stratify([(0.0, 0.05, 5), (0.05, 0.10, 5)])
        for target in ("geneA", "geneB"):
            self.assertTrue((self.out / f"{target}.bins.txt").is_file())
            self.assertTrue((self.out / f"{target}.bins.json").is_file())
        combined = json.loads((self.out / "bins_summary.json").read_text())
        self.assertEqual(set(combined.keys()), {"geneA", "geneB"})



class TestScoreFiles(unittest.TestCase):
    """Score files are recognised by content, and the target comes from inside them."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def write(self, name, obj, gz=False):
        import gzip
        data = json.dumps(obj).encode()
        (self.dir / name).write_bytes(gzip.compress(data) if gz else data)
        return self.dir / name

    def run_cli(self, *argv):
        import contextlib
        import io
        from disassembler.cli import build_parser
        ns = build_parser().parse_args(["stratify", *argv])
        err = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
            rc = ns.func(ns)
        return rc, err.getvalue()

    def test_any_name_and_gzip_and_other_files_passed_over(self):
        self.write("dataset_17.dat", {"target": "ns5", "scores": {"SRR1": 0.9}}, gz=True)
        self.write("summary.json", {"ns5": {"n": 1, "max_score": 0.9}})
        self.write("ns5.bins.json", {"target": "ns5", "bins": []})
        (self.dir / "ns5.accessions.txt").write_text("SRR1\n")
        found = hk.find_score_files([self.dir])
        self.assertEqual([(f.target, f.scores) for f in found], [("ns5", {"SRR1": 0.9})])

    def test_finding_nothing_is_an_error(self):
        (self.dir / "notes.txt").write_text("nothing here\n")
        rc, err = self.run_cli("-o", str(self.dir), "--score-bin", "0:1:5")
        self.assertEqual(rc, 1)
        self.assertIn("no harvest score files", err)
        self.assertIn("1 other file(s) there are not score files", err)

    def test_a_pre_0_1_1_file_needs_its_target_named(self):
        old = self.write("ns5.scores.json", {"SRR1": 0.9, "SRR2": 0.4})
        rc, err = self.run_cli("-o", str(self.dir), "--score-bin", "0:1:5")
        self.assertEqual(rc, 1)
        self.assertIn("no target inside", err)
        out = self.dir / "out"
        rc, err = self.run_cli("-o", str(out), "--scores", f"ns5={old}", "--score-bin", "0:1:5")
        self.assertEqual(rc, 0, err)
        self.assertEqual((out / "ns5.bins.txt").read_text().split(), ["SRR1", "SRR2"])
        rc, err = self.run_cli("-o", str(self.dir), "--target", "gyrA", "--score-bin", "0:1:5")
        self.assertEqual(rc, 0, err)
        self.assertTrue((self.dir / "gyrA.bins.txt").is_file(), "the one unnamed file is gyrA")

    def test_target_filters_and_must_exist(self):
        self.write("a.json", {"target": "geneA", "scores": {"SRR1": 0.9}})
        self.write("b.json", {"target": "geneB", "scores": {"SRR2": 0.9}})
        rc, err = self.run_cli("-o", str(self.dir), "--target", "geneB", "--score-bin", "0:1:5")
        self.assertEqual(rc, 0, err)
        self.assertEqual(json.loads((self.dir / "bins_summary.json").read_text()).keys(),
                         {"geneB"})
        rc, err = self.run_cli("-o", str(self.dir), "--target", "geneC", "--score-bin", "0:1:5")
        self.assertEqual(rc, 1)
        self.assertIn("no score file for target 'geneC'", err)

    def test_bad_targets_are_refused(self):
        self.write("a.json", {"target": "../evil", "scores": {"SRR1": 0.9}})
        with self.assertRaisesRegex(ValueError, "cannot be used as a file name"):
            hk.find_score_files([self.dir])
        (self.dir / "a.json").unlink()
        self.write("a.json", {"target": "g", "scores": {"SRR1": 0.9}})
        self.write("b.json", {"target": "g", "scores": {"SRR2": 0.9}})
        with self.assertRaisesRegex(ValueError, "in two score files"):
            hk.find_score_files([self.dir])

    def test_a_named_file_that_is_not_scores_is_an_error_not_skipped(self):
        p = self.write("summary.json", {"ns5": {"n": 1}})
        with self.assertRaisesRegex(ValueError, "not a harvest score file"):
            hk.find_score_files([p])


if __name__ == "__main__":
    unittest.main()
