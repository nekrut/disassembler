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
from pathlib import Path

GZIP_MAGIC = b"\x1f\x8b"
ZSTD_MAGIC = b"\x28\xb5\x2f\xfd"
ZSTD_REFUSED = ("zstd-compressed; only plain and gzip-compressed input is read. "
                "Decompress it first (zstd -d)")


def text_stream(raw) -> io.TextIOWrapper:
    """UTF-8 text from a binary stream, gunzipped when it starts with the gzip magic bytes.
    A stream without `peek` (an HTTP body) is buffered first."""
    buf = raw if hasattr(raw, "peek") else io.BufferedReader(raw)
    head = buf.peek(4)[:4]
    if head.startswith(ZSTD_MAGIC):
        raise ValueError(f"the input is {ZSTD_REFUSED}")
    if head.startswith(GZIP_MAGIC):
        buf = gzip.GzipFile(fileobj=buf)
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
