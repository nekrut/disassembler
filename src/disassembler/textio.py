"""Read input files as text, gzip-compressed or not.

Every file disassembler reads from a user (FASTA references and haplotypes, LexicMap tables,
kmindex shards) goes through here, so any of them can be gzipped. Compression is decided by
the first bytes, never by a `.gz` name: a Galaxy dataset arrives as `dataset_<uuid>.dat`
either way. The bytes are peeked, not read, so a pipe (`<(zcat x)`) or standard input, which
cannot be opened twice, loses nothing to the check. zstd is recognised and refused with a
message rather than decoded as garbage.
"""

from __future__ import annotations

import gzip
import io
import zlib
from pathlib import Path

GZIP_MAGIC = b"\x1f\x8b"
ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"
ZSTD_REFUSED = ("zstd-compressed; only plain and gzip-compressed input is read. "
                "Decompress it first (zstd -d)")


class _Gunzip(io.RawIOBase):
    """A gzip stream whose read errors are `ValueError`, naming the source.

    PL-21. `gzip.GzipFile` raises `EOFError` for a truncated member and `zlib.error` for a
    corrupt one. Neither is caught by `except (OSError, ValueError)`, which is what
    `lexicmap`, `msa`, `sieve`, `stratify`, `reconstruct` and `from-lexicmap` all use, so a
    half-downloaded `.gz` ended in a traceback from every one of them. Fixed here rather than
    in six `except` clauses, because here is where the cause is known.

    `BadGzipFile` (a flipped byte, so the CRC fails) is an `OSError` and was already caught,
    but the message it carries is "CRC check failed 0x... != 0x..." and names no file. It is
    wrapped for that reason alone.
    """

    def __init__(self, gz, name):
        self._gz, self._name = gz, name

    def readable(self):
        return True

    def readinto(self, b):
        try:
            return self._gz.readinto(b)
        except (EOFError, zlib.error, gzip.BadGzipFile) as e:
            raise ValueError(f"{self._name}: truncated or corrupt gzip data ({e})") from e

    def close(self):
        try:
            self._gz.close()
        finally:
            super().close()


def text_stream(raw, name: str = None) -> io.TextIOWrapper:
    """UTF-8 text from a binary stream, gunzipped when it starts with the gzip magic bytes.
    A stream without `peek` (an HTTP body) is buffered first."""
    buf = raw if hasattr(raw, "peek") else io.BufferedReader(raw)
    head = buf.peek(4)[:4]
    if head.startswith(ZSTD_MAGIC):
        raise ValueError(f"the input is {ZSTD_REFUSED}")
    if head.startswith(GZIP_MAGIC):
        label = name or getattr(raw, "name", None) or "the input"
        buf = io.BufferedReader(_Gunzip(gzip.GzipFile(fileobj=buf), label))
    return io.TextIOWrapper(buf, encoding="utf-8")


def open_text(path) -> io.TextIOWrapper:
    """`open(path)` for reading text, except that a gzipped file is decompressed."""
    raw = open(path, "rb")
    try:
        return text_stream(raw)
    except ValueError:
        raw.close()
        raise ValueError(f"{path}: {ZSTD_REFUSED}") from None


def read_text(path) -> str:
    """`Path(path).read_text()`, except that a gzipped file is decompressed."""
    with open_text(Path(path)) as fh:
        return fh.read()


def decode_bytes(data: bytes) -> str:
    """Bytes already in memory (a downloaded file) as text, gunzipped if they are gzip."""
    return text_stream(io.BytesIO(data)).read()
