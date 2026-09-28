#!/usr/bin/env python3
"""test_stratify.py -- unit tests for the cohort stratifier.

--score-bin: arbitrary half-open [LO, HI) ranges, no backfill, nonzero scores only.
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
        (self.out / "t.scores.json").write_text(json.dumps(scores))

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
        hk.stratify_bins(self.out, bins)
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
        hk.stratify_bins(self.out, [(0.0, 0.05, 15)])
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
        hk.stratify_bins(self.out, [(0.0, 0.05, 10), (0.05, 0.10, 10)])
        s = self._summary()
        self.assertEqual(s["bins"][0]["selected"], 3, "short bin gives what it has")
        self.assertEqual(s["bins"][0]["available"], 3)
        self.assertEqual(s["bins"][1]["selected"], 10, "the other bin is NOT topped up above 10")
        self.assertEqual(s["total_selected"], 13)

    def test_output_txt_is_sorted_by_descending_score(self):
        scores = {f"x{i:03d}": 0.01 + i * 1e-3 for i in range(100)}
        self._write_scores(scores)
        hk.stratify_bins(self.out, [(0.0, 0.15, 5)])
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
        hk.stratify_bins(self.out, [(0.05, 0.10, 5)])
        self.assertEqual(sorted(self._txt()), ["lo", "mid"])

    def test_multi_target_writes_per_target_files(self):
        (self.out / "geneA.scores.json").write_text(json.dumps({f"a{i}": 0.02 + i * 1e-4 for i in range(20)}))
        (self.out / "geneB.scores.json").write_text(json.dumps({f"b{i}": 0.07 + i * 1e-4 for i in range(20)}))
        hk.stratify_bins(self.out, [(0.0, 0.05, 5), (0.05, 0.10, 5)])
        for target in ("geneA", "geneB"):
            self.assertTrue((self.out / f"{target}.bins.txt").is_file())
            self.assertTrue((self.out / f"{target}.bins.json").is_file())
        combined = json.loads((self.out / "bins_summary.json").read_text())
        self.assertEqual(set(combined.keys()), {"geneA", "geneB"})


if __name__ == "__main__":
    unittest.main()
