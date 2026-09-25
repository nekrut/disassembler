#!/usr/bin/env python3
"""test_harvest.py -- unit tests for the Galaxy kmindex shard harvester.

Covers the failure accounting: a shard fetch that exhausts its retries must be
distinguishable from a shard that legitimately matched nothing, and a run in which
every shard failed must not report success.
"""
import sys
import tempfile
import unittest
from pathlib import Path

try:
    from disassembler.cmd import harvest as hk
except ImportError as e:  # pragma: no cover - environment-dependent
    print(f"SKIP: disassembler.cmd.harvest not importable ({e})")
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


if __name__ == "__main__":
    unittest.main()
