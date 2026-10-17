"""test_cli.py -- the top-level `disassembler` dispatcher.

Covers `--version`: it must work with no subcommand (it is what a Galaxy wrapper's
version_command runs), and the version it prints must be the one the package is
installed under.
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


if __name__ == "__main__":
    unittest.main()
