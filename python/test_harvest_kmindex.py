#!/usr/bin/env python3
"""test_harvest_kmindex.py -- unit tests for the Galaxy kmindex shard harvester.

Covers the failure accounting: a shard fetch that exhausts its retries must be
distinguishable from a shard that legitimately matched nothing, and a run in which every
shard failed must not report success.

Skips cleanly when `requests` is unavailable, since harvest_kmindex imports it at module
scope and CI does not install third-party packages.
"""
import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

try:
    import harvest_kmindex as hk
except ImportError as e:  # pragma: no cover - environment-dependent
    print(f"SKIP: harvest_kmindex not importable ({e})")
    sys.exit(0)


class FakeResponse:
    def __init__(self, status_code=200, text=""):
        self.status_code = status_code
        self.text = text


class FakeSession:
    """Records calls and replays a scripted sequence of responses or exceptions."""

    def __init__(self, script):
        self.script = list(script)
        self.calls = 0

    def get(self, url, params=None, timeout=None):
        self.calls += 1
        item = self.script[min(self.calls - 1, len(self.script) - 1)]
        if isinstance(item, Exception):
            raise item
        return item


class TestFetchShard(unittest.TestCase):
    def setUp(self):
        self.slept = []
        self._real_sleep = hk.time.sleep
        hk.time.sleep = self.slept.append

    def tearDown(self):
        hk.time.sleep = self._real_sleep

    def test_a_failed_fetch_is_none_not_an_empty_dict(self):
        # {} is what a shard with no hits legitimately returns, so exhausted retries must
        # not produce it: that conflation is what let an expired API key look like a
        # successful harvest of nothing.
        sess = FakeSession([FakeResponse(403, "denied")])
        got = hk.fetch_shard(("shard1", "abc"), "http://x", sess)
        self.assertIsNone(got)
        self.assertEqual(sess.calls, hk.ATTEMPTS)

    def test_an_empty_shard_is_distinguishable_from_a_failure(self):
        sess = FakeSession([FakeResponse(200, "{}")])
        got = hk.fetch_shard(("shard1", "abc"), "http://x", sess)
        self.assertEqual(got, {})
        self.assertIsNotNone(got)
        self.assertEqual(sess.calls, 1, "a good response must not be retried")

    def test_no_backoff_after_the_final_attempt(self):
        # The last sleep bought nothing and cost ~2.5 s per shard.
        sess = FakeSession([FakeResponse(500, "boom")])
        hk.fetch_shard(("shard1", "abc"), "http://x", sess)
        self.assertEqual(len(self.slept), hk.ATTEMPTS - 1,
                         f"one sleep per retry gap, got {self.slept}")
        self.assertAlmostEqual(sum(self.slept), 5.0, places=6)

    def test_exceptions_are_retried_then_reported(self):
        sess = FakeSession([ConnectionError("no route")])
        self.assertIsNone(hk.fetch_shard(("s", "d"), "http://x", sess))
        self.assertEqual(sess.calls, hk.ATTEMPTS)

    def test_a_late_success_still_returns_data(self):
        sess = FakeSession([
            FakeResponse(500, "boom"),
            FakeResponse(200, '{"shard": {"target": {"SRR1": 0.9}}}'),
        ])
        got = hk.fetch_shard(("s", "d"), "http://x", sess)
        self.assertEqual(got, {"shard": {"target": {"SRR1": 0.9}}})
        self.assertEqual(sess.calls, 2)


class FakeCollectionResponse:
    """The /contents/dataset_collections payload harvest_collection walks."""

    status_code = 200

    def __init__(self, n):
        self._n = n

    def json(self):
        return {
            "elements": [
                {"element_identifier": f"shard{i}", "object": {"id": f"id{i}",
                                                               "history_content_type": "dataset"}}
                for i in range(self._n)
            ]
        }


class _NullSession:
    headers = {}
    closed = False

    def mount(self, *a, **kw):
        pass

    def close(self):
        type(self).closed = True

    def get(self, *a, **kw):
        raise AssertionError("fetch_shard is stubbed in these tests")


class TestHarvestAccounting(unittest.TestCase):
    """A harvest that fetched nothing must not report success."""

    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.out = Path(self._tmp.name)
        self._real_get = hk.requests.get
        self._real_session = hk.requests.Session
        self._real_fetch = hk.fetch_shard
        self._real_adapter = hk.HTTPAdapter
        hk.requests.get = lambda *a, **kw: FakeCollectionResponse(3)
        hk.requests.Session = _NullSession
        hk.HTTPAdapter = lambda *a, **kw: None

    def tearDown(self):
        hk.requests.get = self._real_get
        hk.requests.Session = self._real_session
        hk.fetch_shard = self._real_fetch
        hk.HTTPAdapter = self._real_adapter
        self._tmp.cleanup()

    def _harvest(self):
        return hk.harvest_collection("c1", "h1", "http://x", "key", self.out, max_workers=2)

    def test_every_shard_failing_raises_instead_of_reporting_success(self):
        # The reported behaviour: "Harvest complete: 0 unique SRA runs across 0 targets",
        # a UNION file holding a single newline, and exit 0.
        hk.fetch_shard = lambda *a, **kw: None
        with self.assertRaises(RuntimeError) as cm:
            self._harvest()
        self.assertIn("all 3 shard fetches failed", str(cm.exception))

    def test_all_shards_legitimately_empty_is_not_an_error(self):
        # No hits anywhere is a real answer, and must stay distinct from "nothing fetched".
        hk.fetch_shard = lambda *a, **kw: {}
        per_target = self._harvest()
        self.assertEqual(dict(per_target), {})
        self.assertTrue((self.out / "summary.json").is_file())

    def test_the_session_is_closed_when_the_harvest_finishes(self):
        # 32 pooled connections against one host outlived the harvest.
        _NullSession.closed = False
        hk.fetch_shard = lambda *a, **kw: {}
        self._harvest()
        self.assertTrue(_NullSession.closed)

    def test_the_session_is_closed_even_when_the_harvest_raises(self):
        _NullSession.closed = False
        hk.fetch_shard = lambda *a, **kw: (_ for _ in ()).throw(RuntimeError("boom"))
        with self.assertRaises(Exception):
            self._harvest()
        self.assertTrue(_NullSession.closed)

    def test_a_partial_failure_still_returns_what_was_harvested(self):
        calls = {"n": 0}

        def flaky(*a, **kw):
            calls["n"] += 1
            if calls["n"] == 1:
                return {"shardA": {"gene": {"SRR1": 0.91}}}
            return None

        hk.fetch_shard = flaky
        per_target = self._harvest()
        self.assertEqual(dict(per_target), {"gene": {"SRR1": 0.91}})



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


class TestStratifyDiversity(unittest.TestCase):
    """--diversity below 5 used to raise ZeroDivisionError before writing anything."""

    def setUp(self):
        import json
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.out = Path(self._tmp.name)
        (self.out / "t.scores.json").write_text(
            json.dumps({"a": 0.995, "b": 0.991, "c": 0.85, "d": 0.6})
        )

    def tearDown(self):
        self._tmp.cleanup()

    def test_small_cohorts_no_longer_crash(self):
        for k in (1, 2, 3, 4, 10, 2500):
            with self.subTest(diversity=k):
                hk.stratify_diversity(self.out, total_target=k)
                written = self.out / f"t.diversity_{k}.txt"
                self.assertTrue(written.is_file(), f"--diversity {k} wrote no cohort")
                n = len([x for x in written.read_text().splitlines() if x])
                self.assertEqual(n, min(k, 4), f"--diversity {k} selected {n}")


class TestDivergentBracket(unittest.TestCase):
    """The fourth bracket is an `else`, so it inherits the search's own threshold.

    The kmindex search is launched at containment >= 0.30 to reach divergent homologs and
    filters server-side, so a second floor here can only discard what the search was told to
    find. The chain used to end at `elif sc >= 0.50` with nothing catching the remainder.
    """

    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.out = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def _run(self, scores, target):
        import json
        (self.out / "t.scores.json").write_text(json.dumps(scores))
        hk.stratify_diversity(self.out, total_target=target)
        return [x for x in (self.out / f"t.diversity_{target}.txt").read_text().split() if x]

    def test_nothing_falls_through_when_there_is_room(self):
        scores = {"canon": 0.995, "high": 0.85, "mod": 0.62,
                  "div_a": 0.47, "div_b": 0.38, "div_c": 0.31}
        got = self._run(scores, 2500)
        self.assertEqual(sorted(got), sorted(scores),
                         "every harvested accession must reach the cohort when it fits")

    def test_the_divergent_band_is_represented_under_a_binding_cap(self):
        # 40 near-identical accessions and 10 divergent ones, cohort capped at 10.
        # Taking the top 10 by score would return zero divergent.
        scores = {f"hi{i}": 0.99 for i in range(40)}
        scores.update({f"div{i}": 0.31 + i * 0.01 for i in range(10)})
        got = self._run(scores, 10)
        self.assertEqual(len(got), 10)
        n_div = sum(1 for a in got if a.startswith("div"))
        self.assertGreaterEqual(n_div, 1, f"no divergent accession selected: {got}")

    def test_the_summary_counts_the_bracket_it_names(self):
        import json
        self._run({"a": 0.995, "b": 0.85, "c": 0.62, "d": 0.4, "e": 0.33}, 2500)
        smry = json.load(open(self.out / "diversity_2500_summary.json"))["t"]
        self.assertEqual(smry["n_canonical_ge_99"], 1)
        self.assertEqual(smry["n_high_homology_80_99"], 1)
        self.assertEqual(smry["n_moderate_50_80"], 1)
        self.assertEqual(smry["n_divergent_below_50"], 2)
        self.assertEqual(smry["n_available"], 5)
        # The old key named the 0.50-0.80 bracket "divergent" while the genuinely
        # divergent accessions were being discarded.
        self.assertNotIn("n_divergent_50_80", smry)

    def test_the_cohort_floor_is_the_searchs_floor_not_ours(self):
        scores = {"a": 0.99, "b": 0.72, "c": 0.30}
        got = self._run(scores, 2500)
        self.assertIn("c", got, "an accession at the search threshold must survive")


class TestSafeTargetName(unittest.TestCase):
    """Target names come from shard JSON: untrusted, used as a dict key AND a filename."""

    def test_a_blank_name_is_rejected_rather_than_raising(self):
        # raw.split()[0] used to raise IndexError here, which escaped as_completed and tore
        # down the executor mid-harvest with nothing written.
        for blank in ("", "   ", "\t\n"):
            self.assertIsNone(hk.safe_target_name(blank))

    def test_phix174_collapsing_to_empty_is_rejected(self):
        # Reachable without a malformed shard: the phiX174 strip empties this name, which
        # then wrote dotfiles.
        self.assertIsNone(hk.safe_target_name("phiX174"))
        self.assertEqual(hk.safe_target_name("phiX174_gp1"), "gp1")

    def test_path_separators_and_traversal_cannot_escape_out_dir(self):
        self.assertIsNone(hk.safe_target_name(".."))
        self.assertIsNone(hk.safe_target_name("."))
        self.assertEqual(hk.safe_target_name("a/b"), "a_b")
        self.assertEqual(hk.safe_target_name("/etc/passwd"), "_etc_passwd")
        for name in ("..", ".", "a/b", "/etc/passwd", "../../x"):
            got = hk.safe_target_name(name)
            if got is not None:
                self.assertNotIn("/", got)
                self.assertNotIn("\\", got)
                self.assertNotEqual(got, "..")

    def test_ordinary_names_survive_intact(self):
        self.assertEqual(hk.safe_target_name("pfcrt"), "pfcrt")
        self.assertEqual(hk.safe_target_name("hiv1_pol  extra words"), "hiv1_pol")
        self.assertEqual(hk.safe_target_name("candida_auris_erg11"), "candida_auris_erg11")


class TestBinStratification(unittest.TestCase):
    """--bin gives the caller arbitrary [LO, HI) ranges with per-range counts, no backfill,
    and only nonzero-score accessions. It exists so users can drive walker experiments off
    their own coverage histograms without being tied to the four --diversity brackets.
    """

    def setUp(self):
        import tempfile
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
