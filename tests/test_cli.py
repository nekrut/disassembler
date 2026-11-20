"""test_cli.py -- the top-level `disassembler` dispatcher.

Covers `--version`: it must work with no subcommand (it is what a Galaxy wrapper's
version_command runs), and the version it prints must be the one the package is
installed under.

Covers the 0.1.x subcommand names (`walker`, `lexicmap`): they still dispatch to `walk`
and `tile`, say so on stderr, and stay out of `--help`.
"""

import contextlib
import io
import sys
import tomllib
import unittest
from pathlib import Path

try:
    import disassembler
    from disassembler import cli
except ImportError as e:  # pragma: no cover - environment-dependent
    print(f"SKIP: disassembler not importable ({e})")
    sys.exit(0)

REPO_ROOT = Path(__file__).resolve().parents[1]


class TestVersion(unittest.TestCase):

    def test_version_prints_the_package_version_and_exits_0(self):
        out = io.StringIO()
        with contextlib.redirect_stdout(out), self.assertRaises(SystemExit) as cm:
            cli.main(["--version"])
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(out.getvalue().strip(), f"disassembler {disassembler.__version__}")

    def test_installed_metadata_reads_the_same_version_string(self):
        # A second copy of the number in pyproject.toml would let `--version` and the
        # installed package disagree after a bump that touched only one of them.
        project = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text())
        self.assertNotIn("version", project["project"])
        self.assertIn("version", project["project"]["dynamic"])
        self.assertEqual(
            project["tool"]["setuptools"]["dynamic"]["version"],
            {"attr": "disassembler.__version__"},
        )


class TestRenamedSubcommands(unittest.TestCase):

    def _main(self, argv):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            try:
                rc = cli.main(argv)
            except SystemExit as e:
                rc = e.code
        return rc, out.getvalue(), err.getvalue()

    def test_lexicmap_still_runs_as_tile_and_says_so(self):
        # `-h` is enough to show which parser the old name reached.
        rc_old, out_old, err_old = self._main(["lexicmap", "-h"])
        rc_new, out_new, err_new = self._main(["tile", "-h"])
        self.assertEqual((rc_old, rc_new), (0, 0))
        self.assertEqual(out_old, out_new)
        self.assertIn("usage: disassembler tile", out_new)
        self.assertIn("`lexicmap` is now `tile`", err_old)
        self.assertEqual(err_new, "")

    def test_walker_still_runs_as_walk_and_says_so(self):
        from disassembler.cmd import walker
        seen = []
        real = walker.run
        walker.run = lambda ns: seen.append(ns.passthrough) or 0
        try:
            rc_old, _, err_old = self._main(["walker", "--version"])
            rc_new, _, err_new = self._main(["walk", "--version"])
        finally:
            walker.run = real
        self.assertEqual((rc_old, rc_new), (0, 0))
        self.assertEqual(seen, [["--version"], ["--version"]])
        self.assertIn("`walker` is now `walk`", err_old)
        self.assertEqual(err_new, "")

    def test_help_lists_only_the_current_names(self):
        rc, out, _ = self._main(["-h"])
        self.assertEqual(rc, 0)
        listed = out.split("positional arguments:")[-1]
        for name in ("walk ", "tile "):
            self.assertIn(f"    {name}", listed)
        for name in ("walker ", "lexicmap "):
            self.assertNotIn(f"    {name}", listed)


if __name__ == "__main__":
    unittest.main()
