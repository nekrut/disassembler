#!/usr/bin/env python3
"""Verify a walker build against the golden-truth fixtures (issue #2 criteria 2 and 3).

Usage: verify_golden.py [--fixture-dir DIR] [--workdir DIR]
Exit 0 = concordant, 1 = regression.

Asserts only fields proven stable across PYTHONHASHSEED. The bubble 'merge' node is
deliberately NOT asserted (see _determinism in each fixture).
"""
import sys, json, argparse, subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fixture-dir", default=str(HERE / "golden"))
    ap.add_argument("--workdir", default="/tmp/verify_golden")
    args = ap.parse_args()

    fixtures = sorted(Path(args.fixture_dir).glob("*.golden.json"))
    if not fixtures:
        print(f"no fixtures in {args.fixture_dir}"); return 1

    fails = 0
    for fp in fixtures:
        exp = json.loads(fp.read_text())
        acc = exp["accession"]
        qmap = {"query_hiv1_pol.fasta": HERE.parent / "queries/query_hiv1_pol.fasta",
                "pfcrt_only.fasta": HERE.parent / "queries/pfcrt_only.fasta"}
        q = qmap[exp["query"]]
        out = subprocess.run(
            [sys.executable, str(HERE / "capture_golden.py"), "--accession", acc,
             "--query", str(q), "--workdir", f"{args.workdir}/{acc}"],
            capture_output=True, text=True)
        if out.returncode != 0:
            print(f"FAIL {acc}: capture failed\n{out.stderr[-500:]}"); fails += 1; continue
        got = json.loads(out.stdout)

        checks = []
        checks.append(("graph", exp["graph"], got["graph"]))
        checks.append(("primary_sha256",
                       [p["seq_sha256"] for p in exp["primary"]],
                       [p["seq_sha256"] for p in got["primary"]]))
        checks.append(("n_bubbles", exp["quasispecies"]["n_bubbles"],
                       got["quasispecies"]["n_bubbles"]))
        norm = lambda bs: sorted([{"fork": b["fork"], "branches": b["branches"]} for b in bs],
                                 key=lambda x: (x["fork"], json.dumps(x["branches"])))
        checks.append(("quasispecies_branches", exp["quasispecies"]["bubbles"],
                       norm(got["quasispecies"]["bubbles"])))

        bad = [n for n, e, g in checks if e != g]
        if bad:
            print(f"FAIL {acc}: mismatch in {', '.join(bad)}"); fails += 1
        else:
            print(f"PASS {acc}  ({exp['quasispecies']['n_bubbles']} bubbles, "
                  f"primary {exp['primary'][0]['length']} nt)")
    return 1 if fails else 0

sys.exit(main())
