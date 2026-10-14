"""Locate the external binaries this suite shells out to.

Two binaries live outside the Python package:

  * `logan-walker` (the Rust binary in `walker/target/release/`)
  * `l1_pipeline_serial` (the k-mer diversity sieve in `extern/kmer-diversity-sieve/`)

Resolution order, in both cases:

  1. An explicit path on the caller (e.g. `--walker-bin`).
  2. An environment override (`DISASSEMBLER_WALKER_BIN`, `DISASSEMBLER_SIEVE_BIN`).
  3. `PATH`.
  4. In-tree build under the repo root, discovered by walking up from this file.

Returns a runnable path or None. `Path.is_file()` accepted files without the execute
bit and let the run die deep in subprocess with a bare PermissionError; `os.access(...,
X_OK)` is the same check `shutil.which` uses on PATH entries.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path


def _repo_root() -> Path:
    # src/disassembler/resolve.py -> repo root is three parents up.
    return Path(__file__).resolve().parents[2]


def _executable(p: str | os.PathLike | None) -> str | None:
    if not p:
        return None
    q = Path(p)
    if q.is_file() and os.access(q, os.X_OK):
        return str(q)
    return None


def _from_path(name: str) -> str | None:
    return shutil.which(name)


def resolve_walker(explicit: str | None = None) -> str | None:
    """Locate `logan-walker`. See module docstring for resolution order."""
    return (
        _executable(explicit)
        or _executable(os.environ.get("DISASSEMBLER_WALKER_BIN"))
        or _from_path("logan-walker")
        or _executable(_repo_root() / "walker" / "target" / "release" / "logan-walker")
    )


def resolve_sieve(explicit: str | None = None) -> str | None:
    """Locate `l1_pipeline_serial`. See module docstring for resolution order."""
    return (
        _executable(explicit)
        or _executable(os.environ.get("DISASSEMBLER_SIEVE_BIN"))
        or _from_path("l1_pipeline_serial")
        or _executable(
            _repo_root() / "extern" / "kmer-diversity-sieve" / "l1_pipeline_serial"
        )
    )
