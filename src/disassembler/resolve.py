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
import re
import shutil
import subprocess
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


# What each tool is called, where it comes from, and how to build it. The four sites that
# report a missing binary used to write their own message: `sieve` said "put it on PATH"
# without naming the file it looks for (`l1_pipeline_serial`, not `kmer-diversity-sieve`,
# which is what a user would reasonably install) or the environment override, and the
# composites exited 1 where the stages exited 127.
TOOLS = {
    "walker": {
        "file": "logan-walker",
        "env": "DISASSEMBLER_WALKER_BIN",
        "flag": "--walker-bin",
        "build": "cd walker && cargo build --release",
        "conda": "conda install -c conda-forge -c bioconda logan-walker",
    },
    "sieve": {
        "file": "l1_pipeline_serial",
        "env": "DISASSEMBLER_SIEVE_BIN",
        "flag": "--sieve-bin",
        "build": "cd extern/kmer-diversity-sieve && make serial",
        "conda": "conda install -c conda-forge -c bioconda kmer-diversity-sieve",
    },
}

EX_TOOL_NOT_FOUND = 127


class ToolNotFound(RuntimeError):
    """A required external binary could not be resolved. Exit 127, as a shell would."""


def not_found_message(tool: str, prog: str = "disassembler") -> str:
    """One message naming the file, the resolution order, and how to get it."""
    t = TOOLS[tool]
    return (
        f"{prog}: {t['file']} not found. It is looked for as {t['flag']}, then "
        f"${t['env']}, then on PATH, then in this checkout. Build it "
        f"(`{t['build']}`) or install it (`{t['conda']}`)."
    )


# The published biocontainer `disassembler:0.1.1--pyhdfd78af_0` carries disassembler
# 0.1.1 beside **logan-walker 0.1.0**: the recipe floor was `>=0.1.0,<0.2` and the image was
# built three days before walker 0.1.1 reached the channel. Everything run from that image
# behaves like 0.1.0, which picks the multi-anchor primary by `cov` rather than `matched` (62%
# of a 20,000-accession SARS-CoV-2 run tagged XBB.1.5, one of the shortest references, against
# 43% Delta, the cohort's real mix) and has no prune-span guard. The 0.1.0 summary JSON has no
# `settings` block, so nothing in a user's outputs reveals which walker produced them.
#
# A recipe floor only constrains the build. This is the stronger guard: it also catches a stale
# conda environment, an `--walker-bin` pointed at an old build, and a hand-built binary on PATH.
MIN_WALKER_VERSION = (0, 1, 1)
_WALKER_VER = re.compile(r"^logan-walker (\d+)\.(\d+)\.(\d+)")


def walker_version_error(path: str) -> str | None:
    """None if `path` is a new enough `logan-walker`, else a message naming what is wrong."""
    try:
        out = subprocess.run([path, "--version"], capture_output=True, text=True,
                             stdin=subprocess.DEVNULL, timeout=30).stdout.strip()
    except (OSError, subprocess.SubprocessError) as e:
        return f"{path}: cannot run --version ({e})"
    m = _WALKER_VER.match(out)
    if not m:
        return f"{path}: unrecognised --version output {out!r}"
    if tuple(map(int, m.groups())) < MIN_WALKER_VERSION:
        need = ".".join(map(str, MIN_WALKER_VERSION))
        return (f"{path} is logan-walker {'.'.join(m.groups())}; this disassembler needs "
                f"logan-walker >= {need}. The biocontainer "
                f"disassembler:0.1.1--pyhdfd78af_0 ships 0.1.0: use a later build, or "
                f"`conda install 'logan-walker>={need}'`.")
    return None
