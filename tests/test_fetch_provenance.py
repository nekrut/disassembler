#!/usr/bin/env python3
"""`disassembler fetch` puts the walker's manifest with the cache it describes.

The walker writes `fetch.summary.json` into `--out-dir`, which records the size, sha256 and
(for anything downloaded) the S3 ETag and Last-Modified of every object it cached. The Logan
objects are unversioned and mutable, so that file is what turns a filled cache into a frozen
input that can be cited rather than a directory of bytes.

The wrapper never passed `-o`, so the walker's default would apply and the manifest would land
in `./walker_output/` beside whatever the user happened to be standing in, describing a cache
somewhere else entirely.
"""
import argparse
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from disassembler.cmd import fetch as fetch_cmd


def parse(argv: list[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser()
    fetch_cmd.add_parser(p.add_subparsers(dest="cmd"))
    return p.parse_args(["fetch", *argv])


class TestFetchArgv(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.cache = str(Path(self.tmp.name) / "unitigs")

    def tearDown(self):
        self.tmp.cleanup()

    def argv_for(self, extra: list[str]) -> list[str]:
        """The argv `run` would hand the walker binary."""
        seen: list[list[str]] = []
        real = fetch_cmd.walker_cmd.run_walker
        fetch_cmd.walker_cmd.run_walker = lambda argv, walker_bin=None: (seen.append(argv), 0)[1]
        try:
            rc = fetch_cmd.run(parse(["-a", "SRR000001", "--unitigs-cache", self.cache, *extra]))
        finally:
            fetch_cmd.walker_cmd.run_walker = real
        self.assertEqual(rc, 0)
        self.assertEqual(len(seen), 1)
        return seen[0]

    def test_the_manifest_defaults_to_the_cache_directory(self):
        argv = self.argv_for([])
        self.assertIn("-o", argv)
        self.assertEqual(argv[argv.index("-o") + 1], self.cache)
        # And it is still a fetch-only run against that same cache.
        self.assertIn("--fetch-only", argv)
        self.assertEqual(argv[argv.index("--unitigs-cache") + 1], self.cache)

    def test_out_dir_overrides_it(self):
        argv = self.argv_for(["--out-dir", "/tmp/elsewhere"])
        self.assertEqual(argv[argv.index("-o") + 1], "/tmp/elsewhere")
        self.assertEqual(argv[argv.index("--unitigs-cache") + 1], self.cache)

    def test_no_walker_output_directory_is_implied(self):
        """The walker's own default (./walker_output) must never be what fetch relies on."""
        self.assertNotIn("walker_output", " ".join(self.argv_for([])))


if __name__ == "__main__":
    unittest.main()
