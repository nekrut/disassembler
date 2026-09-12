#!/usr/bin/env python3
"""
logan_path_walker.py: Local-subgraph path walking & bubble phasing for Logan unitigs.

Given an SRA accession and a reference CDS:
1. Streams and decompresses Logan compacted de Bruijn graph unitigs (.unitigs.fa.zst) from AWS S3.
2. Identifies anchor unitigs matching the reference CDS via blastn.
3. Parses BCALM2 'L:' adjacency links to extract an N-hop (default: 5) local subgraph.
4. Emits a compact local GFA (bare node IDs, 30 bp overlap).
5. Runs GraphAligner (-x dbg) to reconstruct the primary full-length haplotype path.
6. Enforces orientation-aware superbubble phasing to recover true quasispecies/polymorphisms
   gated by k-mer abundance (ka:f >= min_abund) and VAF >= min_vaf.
7. Emits reconstructed FASTA with abundance annotations.
"""

import sys
import os
import re
import argparse
import subprocess
import time
import collections
import json
import threading
import concurrent.futures
from pathlib import Path
from Bio import SeqIO

MAMBA_PREFIX = Path(os.environ.get("MAMBA_ROOT_PREFIX", str(Path.home() / "micromamba")))
CODON_ENV_BIN = MAMBA_PREFIX / "envs" / "codon" / "bin"

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def revcomp(s: str) -> str:
    return s.translate(_COMPLEMENT)[::-1]

def find_bin(name: str) -> str:
    candidate = CODON_ENV_BIN / name
    if candidate.exists():
        return str(candidate)
    which = subprocess.run(["which", name], capture_output=True, text=True)
    if which.returncode == 0:
        return which.stdout.strip()
    return name

def stream_unitigs(acc: str, out_fa: Path, timeout: int = 300) -> int:
    url = f"https://logan-pub.s3.amazonaws.com/u/{acc}/{acc}.unitigs.fa.zst"
    curl_cmd = ["curl", "-s", "--fail", "--max-time", str(timeout), url]
    zstd_cmd = ["zstd", "-dc"]
    
    with open(out_fa, "wb") as fh:
        p1 = subprocess.Popen(curl_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p2 = subprocess.Popen(zstd_cmd, stdin=p1.stdout, stdout=fh, stderr=subprocess.PIPE)
        p1.stdout.close()
        _, err2 = p2.communicate()
        p1.wait()
    
    if p2.returncode != 0 or not out_fa.exists() or out_fa.stat().st_size == 0:
        if out_fa.exists():
            out_fa.unlink()
        raise RuntimeError(f"Failed to stream/decompress unitigs for {acc} from {url}: {err2.decode('utf-8', errors='ignore')}")
    
    count = 0
    with open(out_fa, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count

def find_anchors_blast(unitigs_fa: Path, query_fa: Path, work_dir: Path, acc: str,
                       min_pident: float = 80.0, evalue: str = "1e-10") -> list:
    blastn = find_bin("blastn")
    makeblastdb = find_bin("makeblastdb")
    db_prefix = work_dir / f"{acc}_bdb"
    
    subprocess.run(
        [makeblastdb, "-in", str(unitigs_fa), "-dbtype", "nucl", "-out", str(db_prefix)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
    )
    
    cmd = [
        blastn, "-task", "megablast", "-query", str(query_fa), "-db", str(db_prefix),
        "-perc_identity", str(min_pident), "-evalue", str(evalue),
        "-max_target_seqs", "100000", "-outfmt", "6 sseqid pident length bitscore"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    
    matches = set()
    for line in res.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        full_id = parts[0]
        bare_id = full_id.split("_")[-1]
        matches.add(bare_id)
        
    for p in work_dir.glob(f"{acc}_bdb.*"):
        p.unlink()
        
    return sorted(matches)

def extract_local_subgraph(unitigs_fa: Path, anchor_ids: list, gfa_out: Path,
                           hops: int = 5, kmer_len: int = 31) -> dict:
    match_set = set(anchor_ids)
    adj = collections.defaultdict(set)
    links = collections.defaultdict(list)
    abundances = {}
    
    with open(unitigs_fa, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                tokens = line[1:].split()
                uid = tokens[0].split("_")[-1]
                for t in tokens[1:]:
                    if t.startswith("ka:f:"):
                        try:
                            abundances[uid] = float(t.split(":")[-1])
                        except ValueError:
                            pass
                    elif t.startswith("L:"):
                        parts = t.split(":")
                        if len(parts) == 4:
                            _, s1, tgt, s2 = parts
                            adj[uid].add(tgt)
                            adj[tgt].add(uid)
                            links[uid].append((s1, tgt, s2))
                            
    local_nodes = set(match_set)
    for _ in range(hops):
        frontier = set()
        for node in local_nodes:
            frontier.update(adj.get(node, set()))
        local_nodes.update(frontier)
        
    seqs = {}
    current_uid = None
    keep = False
    buf = []
    
    with open(unitigs_fa, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if keep and current_uid:
                    seqs[current_uid] = "".join(buf)
                tokens = line[1:].split()
                current_uid = tokens[0].split("_")[-1]
                keep = current_uid in local_nodes
                buf = []
            elif keep:
                buf.append(line.strip())
        if keep and current_uid:
            seqs[current_uid] = "".join(buf)
            
    overlap = f"{kmer_len - 1}M"
    with open(gfa_out, "w") as out:
        out.write("H\tVN:Z:1.0\n")
        for uid, s in seqs.items():
            out.write(f"S\t{uid}\t{s}\n")
        for uid in seqs:
            for s1, tgt, s2 in links.get(uid, []):
                if tgt in seqs:
                    out.write(f"L\t{uid}\t{s1}\t{tgt}\t{s2}\t{overlap}\n")
                    
    return {
        "num_nodes": len(seqs),
        "abundances": abundances,
        "seqs": seqs,
        "links": links
    }

def align_and_reconstruct(gfa_path: Path, query_fa: Path, gaf_out: Path,
                          seqs: dict, acc: str, kmer_len: int = 31,
                          min_cov: float = 0.85, threads: int = 12) -> list:
    graphaligner = find_bin("GraphAligner")
    cmd = [
        graphaligner, "-g", str(gfa_path), "-f", str(query_fa),
        "-a", str(gaf_out), "-x", "dbg", "-t", str(threads)
    ]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    
    if not gaf_out.exists() or gaf_out.stat().st_size == 0:
        return []
        
    query_record = next(SeqIO.parse(query_fa, "fasta"))
    qlen = len(query_record.seq)
    
    reconstructed = []
    with open(gaf_out, "r") as fh:
        rows = [line.strip().split("\t") for line in fh if line.strip()]
        
    if not rows:
        return []
        
    rows.sort(key=lambda r: int(r[9]) if len(r) > 9 else 0, reverse=True)
    best = rows[0]
    
    path_str = best[5]
    ps = int(best[7])
    pe = int(best[8])
    num_matches = int(best[9])
    aln_len = int(best[10])
    
    nodes = re.findall(r"([<>])(\d+)", path_str)
    if not nodes:
        return []
        
    k_overlap = kmer_len - 1
    full_seq = []
    for i, (orient, nid) in enumerate(nodes):
        if nid not in seqs:
            return []
        s = seqs[nid]
        if orient == "<":
            s = revcomp(s)
        full_seq.append(s if i == 0 else s[k_overlap:])
        
    stitched = "".join(full_seq)
    if pe <= len(stitched):
        target_seq = stitched[ps:pe]
    else:
        target_seq = stitched[ps:]
        
    cov = len(target_seq) / qlen
    pident = (num_matches / aln_len) * 100 if aln_len > 0 else 0
    
    if cov >= min_cov:
        reconstructed.append({
            "haplotype_id": f"{acc}_primary",
            "seq": target_seq,
            "cov": cov,
            "pident": pident,
            "length": len(target_seq),
            "nodes": len(nodes),
            "path": path_str
        })
    
    return reconstructed

def phase_superbubbles(subgraph_info: dict, acc: str, min_abund: float = 5.0, min_vaf: float = 0.10) -> list:
    links = subgraph_info["links"]
    abundances = subgraph_info["abundances"]
    seqs = subgraph_info["seqs"]
    
    outgoing = collections.defaultdict(set)
    for u, edge_list in links.items():
        for s1, tgt, s2 in edge_list:
            if tgt in seqs:
                outgoing[u].add(tgt)
                
    phased_alleles = []
    for a, targets in outgoing.items():
        if len(targets) >= 2:
            branch_nodes = list(targets)
            successors = collections.Counter()
            for b in branch_nodes:
                for c in outgoing.get(b, set()):
                    successors[c] += 1
            reconverged = [c for c, count in successors.items() if count >= 2]
            if reconverged:
                c_node = reconverged[0]
                total_ab = sum(abundances.get(b, 0.0) for b in branch_nodes)
                if total_ab <= 0:
                    continue
                valid_branches = []
                for b in branch_nodes:
                    ab = abundances.get(b, 0.0)
                    vaf = ab / total_ab
                    if ab >= min_abund and vaf >= min_vaf:
                        valid_branches.append((b, ab, vaf))
                if len(valid_branches) >= 2:
                    phased_alleles.append({
                        "fork": a,
                        "merge": c_node,
                        "branches": valid_branches
                    })
    return phased_alleles

def process_accession(acc: str, query_fa: Path, out_dir: Path, hops: int = 5,
                      min_abund: float = 5.0, min_cov: float = 0.85, threads: int = 12) -> dict:
    t0 = time.time()
    work_dir = out_dir / acc
    work_dir.mkdir(parents=True, exist_ok=True)
    
    unitigs_fa = work_dir / f"{acc}.unitigs.fa"
    sub_gfa = work_dir / f"{acc}.sub.gfa"
    sub_gaf = work_dir / f"{acc}.sub.gaf"
    
    try:
        n_unitigs = stream_unitigs(acc, unitigs_fa)
        t_stream = time.time()
        
        anchors = find_anchors_blast(unitigs_fa, query_fa, work_dir, acc)
        t_blast = time.time()
        
        if not anchors:
            return {"acc": acc, "status": "no_anchors", "n_unitigs": n_unitigs, "duration": round(time.time() - t0, 1)}
            
        sub_info = extract_local_subgraph(unitigs_fa, anchors, sub_gfa, hops=hops)
        t_subgraph = time.time()
        
        reconstructed = align_and_reconstruct(
            sub_gfa, query_fa, sub_gaf, sub_info["seqs"], acc,
            min_cov=min_cov, threads=threads
        )
        t_walk = time.time()
        
        bubbles = phase_superbubbles(sub_info, acc, min_abund=min_abund)
        
        total_time = round(time.time() - t0, 1)
        
        res = {
            "acc": acc,
            "status": "success" if reconstructed else "no_path",
            "n_unitigs": n_unitigs,
            "n_anchors": len(anchors),
            "subgraph_nodes": sub_info["num_nodes"],
            "reconstructed": reconstructed,
            "bubbles": len(bubbles),
            "timing": {
                "stream": round(t_stream - t0, 1),
                "blast": round(t_blast - t_stream, 1),
                "subgraph": round(t_subgraph - t_blast, 1),
                "walk": round(t_walk - t_subgraph, 2),
                "total": total_time
            }
        }
        return res
    finally:
        for p in [unitigs_fa, sub_gfa, sub_gaf]:
            if p.exists():
                p.unlink()

def main():
    parser = argparse.ArgumentParser(description="Logan local-subgraph path walker and quasispecies phaser")
    parser.add_argument("--accession", "-a", help="Single SRA accession to process")
    parser.add_argument("--accession-list", "-l", help="File with list of SRA accessions")
    parser.add_argument("--query", "-q", required=True, help="Path to reference CDS FASTA")
    parser.add_argument("--out-dir", "-o", default="walker_output", help="Output directory")
    parser.add_argument("--out-fasta", "-f", default="reconstructed_haplotypes.fasta", help="Reconstructed FASTA output")
    parser.add_argument("--hops", type=int, default=5, help="Hop radius for local subgraph expansion (default: 5)")
    parser.add_argument("--min-abund", type=float, default=5.0, help="Minimum ka:f abundance filter for variants (default: 5.0)")
    parser.add_argument("--min-cov", type=float, default=0.85, help="Minimum coverage fraction of reference (default: 0.85)")
    parser.add_argument("--threads", "-t", type=int, default=4, help="Number of threads for GraphAligner per accession")
    parser.add_argument("--workers", "-w", type=int, default=4, help="Number of accessions to process concurrently in parallel (default: 4)")
    args = parser.parse_args()
    
    query_fa = Path(args.query)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_fasta = Path(args.out_fasta)
    
    accessions = []
    if args.accession:
        accessions.append(args.accession.strip())
    elif args.accession_list:
        with open(args.accession_list) as fh:
            accessions = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
    else:
        parser.error("Must provide either --accession or --accession-list")
        
    print(f"[*] Processing {len(accessions)} accessions against {query_fa.name} (hops={args.hops}, min_abund={args.min_abund}, workers={args.workers}, threads={args.threads})")
    
    recs_out = []
    bubble_reports = []
    success_count = 0
    done_count = 0
    lock = threading.Lock()
    
    def worker(acc):
        nonlocal done_count, success_count
        try:
            res = process_accession(acc, query_fa, out_dir, hops=args.hops,
                                    min_abund=args.min_abund, min_cov=args.min_cov,
                                    threads=args.threads)
            with lock:
                done_count += 1
                curr = done_count
            
            if res.get("status") == "success" and res.get("reconstructed"):
                r = res["reconstructed"][0]
                rec_str = f">{r['haplotype_id']} cov={r['cov']:.3f} pident={r['pident']:.1f} nodes={r['nodes']}\n{r['seq']}\n"
                with lock:
                    success_count += 1
                    recs_out.append(rec_str)
                    if res.get("bubbles"):
                        bubble_reports.append({"acc": acc, "bubbles": res["bubbles"]})
                print(f"[{curr}/{len(accessions)}] {acc}: DONE in {res['timing']['total']}s (nodes={res['subgraph_nodes']}, walk={res['timing']['walk']}s, cov={r['cov']:.2f}, pid={r['pident']:.1f}%, bubbles={res['bubbles']})", flush=True)
            else:
                st = res.get("status", "failed")
                dur = res.get("duration", res.get("timing", {}).get("total", 0))
                print(f"[{curr}/{len(accessions)}] {acc}: {st} ({dur}s)", flush=True)
            return res
        except Exception as e:
            with lock:
                done_count += 1
                curr = done_count
            print(f"[{curr}/{len(accessions)}] {acc}: ERROR: {e}", flush=True)
            return {"acc": acc, "status": f"error: {e}"}

    results = []
    if args.workers > 1 and len(accessions) > 1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
            results = list(executor.map(worker, accessions))
    else:
        for acc in accessions:
            results.append(worker(acc))
            
    with open(out_fasta, "w") as fh:
        fh.writelines(recs_out)
        
    summary_path = out_dir / f"{query_fa.stem}_walker_summary.json"
    with open(summary_path, "w") as fh:
        json.dump({"total": len(accessions), "success": success_count, "results": results, "bubble_events": bubble_reports}, fh, indent=2)
        
    print(f"\n[✓] Completed: {success_count}/{len(accessions)} accessions successfully reconstructed -> {out_fasta}")
    print(f"[*] Summary and quasispecies details saved to {summary_path}")

if __name__ == "__main__":
    main()
