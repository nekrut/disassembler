"""Galaxy access shared by the stages that read a dataset.

`get_galaxy_api_key` existed twice, in `harvest` and in `tile`, and had already
drifted: harvest's copy took no `api_key` argument at all, so `--api-key` was honoured by one
command and not the other, and both were pinned by their own tests, which is how the drift
survived. One implementation, in the order the help documents.
"""

from __future__ import annotations

import os
from pathlib import Path

DEFAULT_KEY_FILE = Path("/tmp/gxy.txt")


def api_key(api_key: str = None, key_file: str = None) -> str:
    """`--api-key`, then `--api-key-file`, then `$GALAXY_API_KEY`, then `/tmp/gxy.txt`.

    The environment variable used to win over an explicitly named file, which made
    `--api-key-file` silently inert wherever the variable happened to be set.
    """
    if api_key:
        return api_key.strip()
    if key_file and Path(key_file).exists():
        return Path(key_file).read_text().strip()
    if os.environ.get("GALAXY_API_KEY"):
        return os.environ["GALAXY_API_KEY"].strip()
    if DEFAULT_KEY_FILE.exists():
        return DEFAULT_KEY_FILE.read_text().strip()
    raise RuntimeError(
        "Galaxy API key not found. Provide --api-key or --api-key-file, set GALAXY_API_KEY, "
        "or write /tmp/gxy.txt."
    )
