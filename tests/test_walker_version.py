#!/usr/bin/env python3
"""PK-01: this package refuses a logan-walker older than it needs.

The only published 0.1.1 biocontainer, `disassembler:0.1.1--pyhdfd78af_0`, carries
disassembler 0.1.1 beside **logan-walker 0.1.0**: the recipe floor was `>=0.1.0,<0.2` and the
image was built three days before walker 0.1.1 reached the channel. Everything run from that
image behaves like 0.1.0, which picks the multi-anchor primary by `cov` rather than `matched`
(62% of a 20,000-accession SARS-CoV-2 run tagged XBB.1.5, one of the shortest references,
against 43% Delta, the cohort's real mix) and has no prune-span guard. The 0.1.0 summary JSON
has no `settings` block, so nothing in a user's outputs says which walker produced them.

A recipe floor only constrains the build. This is the stronger guard: it also catches a stale
conda environment, a `--walker-bin` pointed at an old build, and a hand-built binary on PATH.
"""
import contextlib
import io
import os
import tempfile
import unittest
from pathlib import Path

from disassembler.resolve import MIN_WALKER_VERSION, walker_version_error

REPO_ROOT = Path(__file__).resolve().parents[1]


def stub(dir_path: Path, text: str) -> str:
    p = dir_path / "logan-walker"
    p.write_text(text)
    p.chmod(0o755)
    return str(p)


class TestWalkerVersionFloor(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def say(self, version: str) -> str | None:
        return walker_version_error(
            stub(self.dir, f'#!/bin/sh\necho "logan-walker {version}"\n'))

    def test_the_version_in_the_bad_container_is_refused(self):
        err = self.say("0.1.0")
        self.assertIsNotNone(err, "walker 0.1.0 must be refused")
        self.assertIn("0.1.0", err)
        self.assertIn("pyhdfd78af_0", err, "name the container this is really about")

    def test_the_floor_itself_and_anything_newer_passes(self):
        for v in (".".join(map(str, MIN_WALKER_VERSION)), "0.1.2", "0.2.0", "1.0.0"):
            with self.subTest(version=v):
                self.assertIsNone(self.say(v))

    def test_an_older_minor_and_major_are_refused(self):
        for v in ("0.0.9", "0.1.0"):
            with self.subTest(version=v):
                self.assertIsNotNone(self.say(v))

    def test_output_that_is_not_a_walker_is_refused_rather_than_guessed_at(self):
        err = walker_version_error(stub(self.dir, '#!/bin/sh\necho "GraphAligner 1.0.17b"\n'))
        self.assertIsNotNone(err)
        self.assertIn("unrecognised", err)

    def test_a_binary_that_cannot_run_is_refused(self):
        err = walker_version_error(str(self.dir / "does-not-exist"))
        self.assertIsNotNone(err)
        self.assertIn("cannot run --version", err)

    def test_the_walker_in_this_checkout_meets_the_floor(self):
        # The source of truth for what this package ships against.
        cargo = (REPO_ROOT / "walker" / "Cargo.toml").read_text()
        line = next(l for l in cargo.splitlines() if l.startswith("version"))
        ver = tuple(int(x) for x in line.split('"')[1].split("-")[0].split("."))
        self.assertGreaterEqual(ver, MIN_WALKER_VERSION,
                                f"walker/Cargo.toml is {ver}, below the floor this package "
                                f"enforces at run time")

    def test_the_walker_passthrough_refuses_an_old_binary(self):
        # `run_walker` is the entry the composites use; the top-level `walker` subcommand
        # execv's through the same resolver.
        from disassembler.cmd import walker as walker_cmd
        old = stub(self.dir, '#!/bin/sh\necho "logan-walker 0.1.0"\n')
        err = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
            rc = walker_cmd.run_walker(["--version"], walker_bin=old)
        self.assertNotEqual(rc, 0)
        self.assertIn("0.1.0", err.getvalue())

    def test_reconstruct_refuses_before_it_fetches_anything(self):
        from disassembler.cmd import pipe_reconstruct as pr
        from disassembler.stage import stage_args
        old = stub(self.dir, '#!/bin/sh\necho "logan-walker 0.1.0"\n')
        q = self.dir / "q.fasta"
        q.write_text(">g\nATGGCTGAATAA\n")
        out = self.dir / "out"
        ns = stage_args(pr, query=[str(q)], accession="SRR1",
                        out_dir=str(out), walker_bin=old, walker_args=[])
        err = io.StringIO()
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
            rc = pr.run(ns)
        self.assertNotEqual(rc, 0)
        self.assertIn("logan-walker 0.1.0", err.getvalue())
        self.assertFalse(any(out.rglob("*.fasta")) if out.exists() else False,
                         "nothing may be produced from a walker this package refuses")


if __name__ == "__main__":
    unittest.main()
