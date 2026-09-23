#!/usr/bin/env python3
"""test_sieve_integration.py -- exercise disassemble.py's --sieve wiring end to end.

The bash test (test_sieve.sh) covers the sieve BINARY. This covers the Python SEAM
that test could not: it stubs stage 1 (the walker, which needs GraphAligner + S3 and
so cannot run on every host) with a trivial script that just drops a known
`reconstructed.fasta` into the out-dir, then runs the real disassemble.py with
--sieve --skip-tree... except --skip-tree now skips the sieve by design. So instead
we stub BOTH ends: a fake walker for stage 1 and a fake codon-MSA/tree for stage 2,
and assert disassemble.py sieved the reconstructions and handed the SIEVED file to
stage 2.

Dependency-free: standard library only. Skips cleanly if the sieve binary can't build.
"""
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DISASSEMBLE = REPO / "python" / "disassemble.py"
SIEVE_BIN = REPO / "extern" / "kmer-diversity-sieve" / "l1_pipeline_serial"
FIXTURE = Path(__file__).resolve().parent / "fixture_reconstructed.fasta"
QUERY = Path(__file__).resolve().parent / "fixture_query.fasta"

FAILURES = []


def check(cond, msg):
    if cond:
        print(f"  PASS: {msg}")
    else:
        print(f"  FAIL: {msg}")
        FAILURES.append(msg)


def ensure_sieve():
    if SIEVE_BIN.exists():
        return "ok"
    if not (SIEVE_BIN.parent / "Makefile").is_file():
        return "no-submodule"
    print("[integration] building sieve (make serial)...")
    subprocess.run(["make", "serial"], cwd=SIEVE_BIN.parent, capture_output=True, text=True)
    return "ok" if SIEVE_BIN.exists() else "no-compiler"


def main():
    state = ensure_sieve()
    if state == "no-submodule":
        print(f"SKIP: extern/kmer-diversity-sieve is empty (run `git submodule update --init "
              f"{SIEVE_BIN.parent.relative_to(REPO)}`) -- not a failure")
        return 0
    if state == "no-compiler":
        print("SKIP: sieve binary failed to build (no C compiler?) -- not a failure")
        return 0
    if not FIXTURE.is_file() or not QUERY.is_file():
        print("SKIP: fixtures missing")
        return 0

    work = Path(tempfile.mkdtemp())
    try:
        # a fake walker: writes the fixture as reconstructed.fasta into -o <dir>, honoring -f
        fake_walker = work / "fake_walker.sh"
        fake_walker.write_text(
            "#!/usr/bin/env bash\n"
            "out=''; fasta=''\n"
            'while [ $# -gt 0 ]; do case "$1" in\n'
            '  -o) out="$2"; shift 2;; -f) fasta="$2"; shift 2;; *) shift;; esac; done\n'
            f'cp "{FIXTURE}" "$fasta"\n'
        )
        fake_walker.chmod(0o755)

        # disassemble.py preflights the MSA-stage tools on PATH (a real check); satisfy them
        # with stubs on a PATH we prepend for the subprocess. The actual MSA/tree step is our
        # fake_msa, so these are never executed, they only have to resolve.
        # Both tools, not just mafft: the preflight covers fasttree as well.
        binstub = work / "bin"
        binstub.mkdir()
        for tool in ("mafft", "fasttree"):
            (binstub / tool).write_text("#!/usr/bin/env bash\nexit 0\n")
            (binstub / tool).chmod(0o755)
        env = dict(os.environ, PATH=f"{binstub}{os.pathsep}{os.environ.get('PATH', '')}")

        # a fake codon-MSA/tree: records which --input it received, writes stub outputs
        fake_msa = work / "fake_msa.py"
        fake_msa.write_text(
            "import sys, pathlib\n"
            "args = sys.argv[1:]\n"
            "inp = args[args.index('--input') + 1]\n"
            "aln = args[args.index('--out-aln') + 1]\n"
            "tree = args[args.index('--out-tree') + 1]\n"
            # Recorded next to the alignment, so each run writes into its own -o directory.
            # A single shared path let run 2 read run 1's file, which made the
            # "stage 2 never ran" branch unreachable and misreported it as a wrong --input.
            "pathlib.Path(aln).with_name('stage2_input.txt').write_text(inp)\n"
            "pathlib.Path(aln).write_text('>stub\\nACGT\\n')\n"
            "pathlib.Path(tree).write_text('(stub);\\n')\n"
        )

        out = work / "out"
        # Point disassemble.py's CODON_MSA at our stub by monkeypatching via env is not
        # supported; instead we run with a patched copy that swaps CODON_MSA. Simplest:
        # copy disassemble.py and rewrite the CODON_MSA constant to the stub.
        patched = work / "disassemble_patched.py"
        src = DISASSEMBLE.read_text().replace(
            'CODON_MSA = ROOT / "python" / "build_walker_codon_msa.py"',
            f'CODON_MSA = pathlib.Path("{fake_msa}")',
        )
        # ensure pathlib is imported in the patched module (it already imports Path; add module)
        src = src.replace("from pathlib import Path", "import pathlib\nfrom pathlib import Path")
        patched.write_text(src)

        cmd = [
            sys.executable, str(patched),
            "-a", "STUBACC",
            "-q", str(QUERY),
            "-o", str(out),
            "--walker-bin", str(fake_walker),
            "--sieve-bin", str(SIEVE_BIN),   # patched copy lives in /tmp, so ROOT-relative default won't resolve
            "--sieve",
        ]
        r = subprocess.run(cmd, capture_output=True, text=True, env=env)
        print(r.stdout)
        if r.returncode != 0:
            print(r.stderr)

        recon = out / "reconstructed.fasta"
        sieved = out / "reconstructed.sieved.fasta"
        stage2_input = out / "stage2_input.txt"

        check(r.returncode == 0, "disassemble.py --sieve exits 0")
        check(recon.is_file(), "reconstructed.fasta written by (fake) stage 1")
        check(sieved.is_file(), "reconstructed.sieved.fasta produced by the sieve stage")
        # the seam: stage 2 must receive the SIEVED file, not the raw reconstructions
        if stage2_input.is_file():
            got = stage2_input.read_text().strip()
            check(got == str(sieved), f"stage 2 --input is the sieved file (got {Path(got).name})")
        else:
            check(False, "stage 2 was invoked (stage2_input.txt recorded)")
        # the sieve actually thinned the 33 near-duplicate fixture haplotypes. The bound is
        # generous on purpose: the wrapper's defaults were bumped from the sieve binary's own
        # archive-scale defaults (tau=0.15 sketch=64 -> 33 -> 1) to gene-scale ones
        # (tau=0.005 sketch=256 -> 33 -> ~28), and asserting any specific count would couple
        # this seam test to a calibration decision it does not own. What matters here is that
        # something got thinned and the sieved file reached stage 2, which the other assertions
        # cover.
        if sieved.is_file():
            n = sum(1 for l in sieved.open() if l.startswith(">"))
            check(1 <= n < 33, f"sieve thinned 33 -> {n} (1 <= n < 33)")

        # `--sieve-sketch` reaches the sieve binary and shows up in the summary log.
        # A silent regression here would let the sketch flag look-like-it-works while
        # the binary ran at its default of 64, which is precisely the failure mode the
        # gene-scale calibration doc characterises.
        out_sk = work / "out_sketch"
        cmd_sk = [c if c != str(out) else str(out_sk) for c in cmd] + ["--sieve-sketch", "128"]
        r_sk = subprocess.run(cmd_sk, capture_output=True, text=True, env=env)
        if r_sk.returncode != 0:
            print(r_sk.stderr)
        check(r_sk.returncode == 0, "--sieve-sketch: disassemble.py exits 0")
        check("--sketch 128" in r_sk.stdout,
              "--sieve-sketch forwards --sketch 128 to the sieve binary CLI")
        check("sketch=128" in r_sk.stdout,
              "--sieve-sketch echoed in the disassemble.py summary log line")

        # off-path: WITHOUT --sieve, stage 2 must receive the RAW reconstructed.fasta
        out2 = work / "out2"
        cmd_off = [
            sys.executable, str(patched),
            "-a", "STUBACC", "-q", str(QUERY), "-o", str(out2),
            "--walker-bin", str(fake_walker),
        ]
        r_off = subprocess.run(cmd_off, capture_output=True, text=True, env=env)
        if r_off.returncode != 0:
            print(r_off.stderr)
        check(r_off.returncode == 0, "off-path: disassemble.py (no --sieve) exits 0")
        stage2_off = out2 / "stage2_input.txt"
        if stage2_off.is_file():
            got_off = stage2_off.read_text().strip()
            check(got_off == str(out2 / "reconstructed.fasta"),
                  "off-path: stage 2 --input is the raw reconstructed.fasta (no sieve)")
        else:
            check(False, "off-path: stage 2 was invoked")

    finally:
        shutil.rmtree(work, ignore_errors=True)

    print()
    if FAILURES:
        print(f"integration: {len(FAILURES)} FAILED")
        return 1
    print("integration: ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
