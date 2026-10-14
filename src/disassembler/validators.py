"""Argparse `type=` converters that reject out-of-range values at parse time.

These map 1:1 onto Galaxy's `<param type="float|integer" min=".." max="..">` and
`<validator type="in_range">`; wiring them into every `add_argument` here means
the Galaxy XML wrapper can rely on the tool to validate ranges and just pass
values through. On a bad value argparse prints `error: argument --foo: <msg>`
and exits 2, which is a clean, predictable Galaxy failure.
"""

from __future__ import annotations

import argparse


def fraction(s: str) -> float:
    """A float in [0.0, 1.0]. For coverage fractions, VAF, prune fractions, LSH tau."""
    try:
        v = float(s)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{s!r} is not a number")
    if not (0.0 <= v <= 1.0):
        raise argparse.ArgumentTypeError(f"{v} out of range [0.0, 1.0]")
    return v


def percentage(s: str) -> float:
    """A float in [0.0, 100.0]. For percent-identity gates."""
    try:
        v = float(s)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{s!r} is not a number")
    if not (0.0 <= v <= 100.0):
        raise argparse.ArgumentTypeError(f"{v} out of range [0.0, 100.0]")
    return v


def positive_int(s: str) -> int:
    """An integer >= 1. For worker counts, sketch sizes, sample caps, row limits."""
    try:
        v = int(s)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{s!r} is not an integer")
    if v < 1:
        raise argparse.ArgumentTypeError(f"{v} must be >= 1")
    return v


def non_negative_int(s: str) -> int:
    """An integer >= 0. For counters that legitimately allow zero (e.g. max stops)."""
    try:
        v = int(s)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{s!r} is not an integer")
    if v < 0:
        raise argparse.ArgumentTypeError(f"{v} must be >= 0")
    return v


def positive_float(s: str) -> float:
    """A float > 0.0. For thresholds that must be strictly positive (e.g. min-abund)."""
    try:
        v = float(s)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{s!r} is not a number")
    if v <= 0.0:
        raise argparse.ArgumentTypeError(f"{v} must be > 0.0")
    return v


def fasta_has_record(path: str) -> int:
    """Count `>` lines in a FASTA file. Zero means the file exists but is header-empty."""
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def validate_unitigs_cache(path: str | None, *, offline: bool) -> str | None:
    """Return an error string if `--unitigs-cache DIR` is misconfigured, else None.

    Semantics:
      - None (unset): fine unless offline (walker fetches straight from S3).
      - Exists and is a directory: fine.
      - Exists but is not a directory: hard fail (the walker treats this as a dir).
      - Does not exist and offline: hard fail (nothing to read from).
      - Does not exist and not offline: fine; walker will create it on write.
    """
    if path is None:
        if offline:
            return "--offline requires --unitigs-cache <DIR> (nothing to read from)."
        return None
    from pathlib import Path
    p = Path(path)
    if p.exists() and not p.is_dir():
        return (
            f"--unitigs-cache {path!r} exists but is not a directory. The walker "
            f"writes per-accession files inside this path; point it at a directory."
        )
    if not p.exists() and offline:
        return (
            f"--unitigs-cache {path!r} does not exist, and --offline forbids fetching "
            f"from S3. Warm the cache with `disassembler fetch --unitigs-cache {path}` "
            f"first, or drop --offline."
        )
    return None


def validate_api_key_file(path: str | None) -> str | None:
    """Return an error string if `--api-key-file` names an unusable file, else None.

    Only runs when the caller actually intends to use the file: an unset path is
    considered fine here (the caller's fallback chain, e.g. $GALAXY_API_KEY, may
    still find a key). If the file exists it must have non-whitespace content; a
    zero-byte or whitespace-only key file is the common footgun this catches.
    """
    if path is None:
        return None
    from pathlib import Path
    p = Path(path)
    if not p.exists():
        return (
            f"--api-key-file {path!r} does not exist. Point it at a file whose "
            f"contents are the Galaxy API key, or unset it and let $GALAXY_API_KEY "
            f"or /tmp/gxy.txt take over."
        )
    if not p.is_file():
        return f"--api-key-file {path!r} is not a regular file."
    try:
        content = p.read_text().strip()
    except OSError as e:
        return f"--api-key-file {path!r} could not be read: {e}"
    if not content:
        return (
            f"--api-key-file {path!r} is empty (or whitespace-only). The Galaxy API "
            f"key must be a non-empty token."
        )
    return None
