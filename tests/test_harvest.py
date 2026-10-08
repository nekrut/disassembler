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


class TestHarvestLocal(unittest.TestCase):
    """Local-files input mode: read shard JSONs from disk, same merge semantics as API."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.shards = Path(self._tmp.name) / "shards"
        self.shards.mkdir()
        self.out = Path(self._tmp.name) / "out"

    def tearDown(self):
        self._tmp.cleanup()

    def _write_shard(self, name, payload):
        import json as _json
        (self.shards / name).write_text(_json.dumps(payload))

    def test_merges_scores_per_target_across_shard_files(self):
        self._write_shard("s1.json", {"idx": {"pfcrt": {"SRR1": 0.4, "SRR2": 0.9}}})
        self._write_shard("s2.json", {"idx": {"pfcrt": {"SRR1": 0.7}}})
        per_target = hk.harvest_local(self.shards, self.out)
        self.assertEqual(dict(per_target["pfcrt"]), {"SRR1": 0.7, "SRR2": 0.9})
        self.assertTrue((self.out / "pfcrt.accessions.txt").is_file())
        self.assertTrue((self.out / "pfcrt.scores.json").is_file())
        self.assertTrue((self.out / "UNION.accessions.txt").is_file())
        self.assertTrue((self.out / "summary.json").is_file())

    def test_accession_order_matches_api_mode(self):
        # Sorted by score descending, same as harvest_collection.
        self._write_shard("s.json", {"idx": {"g": {"A": 0.1, "B": 0.9, "C": 0.5}}})
        hk.harvest_local(self.shards, self.out)
        got = (self.out / "g.accessions.txt").read_text().splitlines()
        self.assertEqual(got, ["B", "C", "A"])

    def test_missing_directory_raises(self):
        with self.assertRaises(RuntimeError) as cm:
            hk.harvest_local(self.shards / "nope", self.out)
        self.assertIn("not found", str(cm.exception))

    def test_empty_directory_raises(self):
        with self.assertRaises(RuntimeError) as cm:
            hk.harvest_local(self.shards, self.out)
        self.assertIn("no .json shard files", str(cm.exception))

    def test_malformed_shard_is_counted_as_failure_not_lost_silently(self):
        self._write_shard("good.json", {"idx": {"g": {"A": 0.5}}})
        (self.shards / "bad.json").write_text("{not json")
        per_target = hk.harvest_local(self.shards, self.out)
        self.assertEqual(dict(per_target["g"]), {"A": 0.5})

    def test_every_shard_malformed_raises(self):
        (self.shards / "a.json").write_text("{")
        (self.shards / "b.json").write_text("<")
        with self.assertRaises(RuntimeError) as cm:
            hk.harvest_local(self.shards, self.out)
        self.assertIn("all 2 shard files failed", str(cm.exception))

    def test_empty_shards_are_a_valid_answer_not_a_failure(self):
        self._write_shard("a.json", {})
        self._write_shard("b.json", {})
        per_target = hk.harvest_local(self.shards, self.out)
        self.assertEqual(dict(per_target), {})
        self.assertTrue((self.out / "summary.json").is_file())

    def test_non_json_files_in_the_directory_are_ignored(self):
        self._write_shard("s.json", {"idx": {"g": {"A": 0.5}}})
        (self.shards / "README").write_text("not a shard")
        (self.shards / "s.json.bak").write_text("stale")
        per_target = hk.harvest_local(self.shards, self.out)
        self.assertEqual(dict(per_target["g"]), {"A": 0.5})


class TestModeValidation(unittest.TestCase):
    """`run()` must route by mode and reject mixed or incomplete flag sets."""

    class Args:
        history_id = None
        collection_id = None
        shards_dir = None

    def test_shards_dir_alone_selects_local_mode(self):
        a = self.Args()
        a.shards_dir = "/tmp/x"
        self.assertEqual(hk._validate_mode(a), ("local", None))

    def test_history_and_collection_together_select_galaxy_mode(self):
        a = self.Args()
        a.history_id = "h"
        a.collection_id = "c"
        self.assertEqual(hk._validate_mode(a), ("galaxy", None))

    def test_shards_dir_with_galaxy_flags_is_rejected(self):
        a = self.Args()
        a.shards_dir = "/tmp/x"
        a.history_id = "h"
        a.collection_id = "c"
        mode, err = hk._validate_mode(a)
        self.assertIsNone(mode)
        self.assertIn("mutually exclusive", err)

    def test_partial_galaxy_flags_rejected(self):
        a = self.Args()
        a.history_id = "h"
        mode, err = hk._validate_mode(a)
        self.assertIsNone(mode)
        self.assertIn("together", err)

    def test_no_flags_rejected(self):
        a = self.Args()
        mode, err = hk._validate_mode(a)
        self.assertIsNone(mode)
        self.assertIn("required", err)


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
