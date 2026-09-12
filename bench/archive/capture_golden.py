#!/usr/bin/env python3
"""Capture golden-truth fixtures. Imports logan_path_walker read-only; modifies nothing."""
import sys, json, hashlib, subprocess, argparse
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))
import logan_path_walker as W

def canon_bubbles(bubbles):
    """Canonicalise: sort branches within a bubble, then sort bubbles. Removes
    set-iteration order dependence so the fixture is hash-seed stable."""
    out = []
    for b in bubbles:
        branches = sorted([[str(n), round(float(ab), 4), round(float(vaf), 6)]
                           for n, ab, vaf in b["branches"]])
        out.append({"fork": str(b["fork"]), "merge": str(b["merge"]), "branches": branches})
    out.sort(key=lambda x: (x["fork"], x["merge"], json.dumps(x["branches"])))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--accession", required=True)
    ap.add_argument("--query", required=True)
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()

    acc, qfa = args.accession, Path(args.query)
    wd = Path(args.workdir); wd.mkdir(parents=True, exist_ok=True)
    uni, gfa, gaf = wd/f"{acc}.unitigs.fa", wd/f"{acc}.sub.gfa", wd/f"{acc}.sub.gaf"
    try:
        n_unitigs = W.stream_unitigs(acc, uni)
        anchors = W.find_anchors_blast(uni, qfa, wd, acc)
        sub = W.extract_local_subgraph(uni, anchors, gfa, hops=5)
        recon = W.align_and_reconstruct(gfa, qfa, gaf, sub["seqs"], acc, min_cov=0.85, threads=2)
        bubbles = W.phase_superbubbles(sub, acc, min_abund=5.0)
        fx = {
            "accession": acc,
            "query": qfa.name,
            "query_md5": hashlib.md5(qfa.read_bytes()).hexdigest(),
            "params": {"hops": 5, "min_abund": 5.0, "min_vaf": 0.10, "min_cov": 0.85, "k": 31},
            "graph": {"n_unitigs": n_unitigs, "n_anchors": len(anchors),
                      "subgraph_nodes": sub["num_nodes"]},
            "primary": ([{"haplotype_id": r["haplotype_id"], "length": len(r["seq"]),
                          "cov": round(r["cov"], 6), "pident": round(r["pident"], 4),
                          "seq_sha256": hashlib.sha256(r["seq"].encode()).hexdigest(),
                          "seq": r["seq"]} for r in recon]),
            "quasispecies": {"n_bubbles": len(bubbles), "bubbles": canon_bubbles(bubbles)},
        }
        print(json.dumps(fx))
    finally:
        for p in (uni, gfa, gaf):
            if p.exists(): p.unlink()

main()
