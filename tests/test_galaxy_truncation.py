"""A Galaxy download that ends early must not read as a complete table.

Run against a real local HTTP server rather than a stubbed stream, because the defect lives
in the socket layer: `http.client.HTTPResponse.read1` returns `b""` when the connection
closes before Content-Length bytes have arrived, which CPython does deliberately for
compatibility, and `io.TextIOWrapper` passes that on as an ordinary end of file. A stub built
on `BytesIO` cannot reproduce it, which is why `TestGalaxyDataset` never caught this.

Four framings, because they fail in three different ways:

  * Content-Length, complete        -> exit 0
  * Content-Length, body cut short  -> silent EOF, caught by the bytes-owed check
  * chunked, terminator missing     -> http.client.IncompleteRead
  * close-delimited gzip, cut       -> EOFError out of the gzip reader

The last two were tracebacks before 0.2.0; `run` now catches both.
"""

import contextlib
import gzip
import io
import json
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from disassembler.cmd import tiler as lms
from disassembler.stage import stage_args

REF = ">gene\nATGGCTGAATAA\n"

# Enough accessions that a cut lands well inside the table rather than inside the header.
TABLE = "query\tsgenome\tpident\tqcovHSP\tqstart\tqend\tqseq\tsseq\n" + "".join(
    f"gene\tSRR{i:05d}\t100.0\t100.0\t1\t12\tATGGCTGAATAA\tATGGCTGAATAA\n"
    for i in range(200)
)


class Handler(BaseHTTPRequestHandler):
    """Serves the two endpoints `-d` uses, framing the body as `server.mode` says."""

    protocol_version = "HTTP/1.1"

    def log_message(self, *a):  # keep the test output clean
        pass

    def do_GET(self):
        if not self.path.endswith("/display"):
            body = json.dumps({"state": "ok"}).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return

        mode = self.server.mode
        full = TABLE.encode()

        if mode in ("complete", "short"):
            sent = full if mode == "complete" else full[: len(full) // 3]
            self.send_response(200)
            # Declared in full either way; `short` simply stops writing and hangs up.
            self.send_header("Content-Length", str(len(full)))
            self.end_headers()
            self.wfile.write(sent)
        elif mode == "chunked_cut":
            self.send_response(200)
            self.send_header("Transfer-Encoding", "chunked")
            self.end_headers()
            part = full[: len(full) // 3]
            self.wfile.write(b"%x\r\n" % len(part) + part + b"\r\n")
            # and no terminating 0-length chunk
        elif mode == "gzip_cut":
            # No Content-Length and no chunking: the body ends when the connection does, so
            # nothing is owed and only the gzip reader can notice the cut.
            blob = gzip.compress(full)
            self.send_response(200)
            self.send_header("Connection", "close")
            self.end_headers()
            self.wfile.write(blob[: len(blob) // 3])
        else:  # pragma: no cover - a typo in a test's mode
            raise AssertionError(f"unknown mode {mode!r}")

        self.wfile.flush()
        self.close_connection = True


class TestTruncatedGalaxyDownload(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        cls.server.mode = "complete"
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.base = f"http://127.0.0.1:{cls.server.server_address[1]}"

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join(timeout=5)

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.ref = self.dir / "gene.fasta"
        self.ref.write_text(REF)

    def tearDown(self):
        self.tmp.cleanup()

    def run_tiler(self, mode, **kw):
        type(self).server.mode = mode
        ns = stage_args(
            lms,
            ref=str(self.ref),
            dataset="f2db41e1fa331b3e",
            base_url=self.base,
            api_key="not-a-real-key",
            out_prefix=str(self.dir / "out") + "/",
            **kw,
        )
        err, out = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = lms.run(ns)
        return rc, err.getvalue()

    def summary(self):
        p = self.dir / "out" / "gene.summary.json"
        return json.loads(p.read_text()) if p.is_file() else None

    def test_a_complete_download_is_read_in_full(self):
        rc, err = self.run_tiler("complete")
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.summary()["run"]["status"], "ok")
        self.assertEqual(self.summary()["cohort"]["clean_full_length"], 200)

    def test_a_body_cut_short_of_content_length_is_refused(self):
        # The defect: this used to exit 0 with a biased prefix of the cohort and
        # "status": "ok", with nothing in any output saying the table was half a table.
        rc, err = self.run_tiler("short")
        self.assertEqual(rc, 1, err)
        self.assertIn("short of its Content-Length", err)
        self.assertEqual(self.summary()["run"]["status"], "failed")
        # No cohort block at all: a half-delivered table has no cohort to report, and the
        # summary names the cause instead of describing a prefix of one as if it were whole.
        self.assertNotIn("cohort", self.summary())
        self.assertIn("short of its Content-Length", self.summary()["run"]["error"])

    def test_a_chunked_body_without_its_terminator_is_refused(self):
        rc, err = self.run_tiler("chunked_cut")
        self.assertEqual(rc, 1, err)
        self.assertIn("IncompleteRead", err)
        self.assertNotIn("Traceback", err)

    def test_a_cut_gzip_body_is_refused(self):
        rc, err = self.run_tiler("gzip_cut")
        self.assertEqual(rc, 1, err)
        self.assertIn("Compressed file ended", err)
        self.assertNotIn("Traceback", err)

    def test_max_records_stops_the_loop_without_tripping_the_check(self):
        # A preview deliberately leaves bytes unread, so owing bytes is expected here and
        # must not be reported as a truncated download.
        rc, err = self.run_tiler("complete", max_records=10)
        self.assertEqual(rc, 0, err)
        self.assertNotIn("Content-Length", err)


if __name__ == "__main__":
    unittest.main()
