#!/usr/bin/env python3
"""
Simulated k-mer Sensitivity and Specificity Benchmark for logan-walker.

Evaluates how varying k-mer size (k) and min_shared_kmers threshold affects:
1. Target Anchor Recall / Sensitivity
2. Gene Recovery via 5-hop BFS Subgraph Expansion
3. Background False Positive Anchors
4. Background Subgraph Expansion (5-hop noise snowball)
5. Subgraph Fraction Guard Violations (1% threshold)
6. Theoretical vs Empirical limits across sequence divergence (100% down to 60%)
"""

import sys
import os
import math
import json
import time
import random
import argparse
from pathlib import Path
from collections import deque, defaultdict
import numpy as np
import pandas as pd

# Paths
HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent.parent
CDS_DIR = REPO_ROOT / "work" / "cds" / "individual"
ECOLI_PATH = REPO_ROOT / "work" / "galaxytools" / "tools" / "gotohscan" / "test-data" / "NC_000913.fna"

# Reverse complement translation table
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]

def canonical_kmer(kmer: str) -> str:
    rc = revcomp(kmer)
    return kmer if kmer < rc else rc

def load_fasta_seq(path: Path) -> str:
    lines = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                lines.append(line.upper())
    return "".join(lines)

# -----------------------------------------------------------------------------
# Mutation Models
# -----------------------------------------------------------------------------

def mutate_uniform(seq: str, divergence: float, rng: random.Random) -> str:
    """Uniform random point substitutions."""
    if divergence <= 0.0:
        return seq
    bases = ['A', 'C', 'G', 'T']
    seq_list = list(seq)
    n_mutations = int(round(len(seq) * divergence))
    indices = rng.sample(range(len(seq)), n_mutations)
    for idx in indices:
        orig = seq_list[idx]
        alt = rng.choice([b for b in bases if b != orig])
        seq_list[idx] = alt
    return "".join(seq_list)

def mutate_codon_biased(seq: str, divergence: float, rng: random.Random) -> str:
    """
    Biology-informed codon substitutions:
    70% at codon pos 3 (wobble / synonymous),
    20% at codon pos 1,
    10% at codon pos 2.
    Transitions (A<->G, C<->T) twice as likely as transversions (Ti/Tv = 2.0).
    """
    if divergence <= 0.0:
        return seq
    seq_list = list(seq)
    n_mutations = int(round(len(seq) * divergence))
    
    # Weights for each base position in the CDS
    weights = []
    for i in range(len(seq)):
        pos_in_codon = i % 3
        if pos_in_codon == 2:   # 3rd position (0-indexed: 2)
            weights.append(0.70)
        elif pos_in_codon == 0: # 1st position
            weights.append(0.20)
        else:                   # 2nd position
            weights.append(0.10)
    
    weights = np.array(weights)
    weights /= weights.sum()
    
    # Choose mutation positions without replacement
    chosen_indices = rng.choices(range(len(seq)), weights=weights, k=n_mutations * 2)
    seen = set()
    indices = []
    for idx in chosen_indices:
        if idx not in seen:
            seen.add(idx)
            indices.append(idx)
            if len(indices) == n_mutations:
                break
    if len(indices) < n_mutations:
        remaining = list(set(range(len(seq))) - seen)
        rng.shuffle(remaining)
        indices.extend(remaining[:n_mutations - len(indices)])
        
    transitions = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
    transversions = {
        'A': ['C', 'T'], 'G': ['C', 'T'],
        'C': ['A', 'G'], 'T': ['A', 'G']
    }
    
    for idx in indices:
        orig = seq_list[idx]
        if orig not in transitions:
            continue
        # Ti/Tv = 2: 2/3 Ti, 1/3 Tv
        if rng.random() < 0.667:
            seq_list[idx] = transitions[orig]
        else:
            seq_list[idx] = rng.choice(transversions[orig])
            
    return "".join(seq_list)

def mutate_clustered_motifs(seq: str, divergence: float, rng: random.Random) -> str:
    """
    Realistic structural evolution:
    Active-site catalytic motifs (~20% of sequence) are protected / conserved (<5% divergence),
    while loop/surface regions absorb the majority of mutations.
    """
    if divergence <= 0.0:
        return seq
    
    seq_list = list(seq)
    L = len(seq)
    
    # Define protected regions (active site motifs in beta-lactamases)
    protected = [
        (int(0.20 * L), int(0.20 * L) + 30),
        (int(0.45 * L), int(0.45 * L) + 30),
        (int(0.75 * L), int(0.75 * L) + 30),
    ]
    
    is_protected = [False] * L
    for start, end in protected:
        for p in range(start, min(end, L)):
            is_protected[p] = True
            
    weights = np.array([0.05 if is_protected[i] else 1.0 for i in range(L)])
    weights /= weights.sum()
    
    n_mutations = int(round(L * divergence))
    chosen = rng.choices(range(L), weights=weights, k=n_mutations * 3)
    seen = set()
    indices = []
    for idx in chosen:
        if idx not in seen:
            seen.add(idx)
            indices.append(idx)
            if len(indices) == n_mutations:
                break
    if len(indices) < n_mutations:
        remaining = list(set(range(L)) - seen)
        rng.shuffle(remaining)
        indices.extend(remaining[:n_mutations - len(indices)])
        
    bases = ['A', 'C', 'G', 'T']
    for idx in indices:
        orig = seq_list[idx]
        seq_list[idx] = rng.choice([b for b in bases if b != orig])
        
    return "".join(seq_list)

# -----------------------------------------------------------------------------
# Unitig Fragmentation & Graph Representation
# -----------------------------------------------------------------------------

class UnitigGraphSim:
    """Simulated Unitig Graph supporting BFS expansion."""
    def __init__(self):
        self.unitigs = []          # list of sequence strings
        self.is_target = []        # bool: True if target unitig, False if background
        self.target_idx = []       # relative index among target unitigs
        self.adj = defaultdict(set) # node_id -> set of neighbor node_ids
        
    def add_unitig(self, seq: str, is_target: bool, target_idx: int = -1) -> int:
        uid = len(self.unitigs)
        self.unitigs.append(seq)
        self.is_target.append(is_target)
        self.target_idx.append(target_idx)
        return uid
        
    def add_edge(self, u1: int, u2: int):
        self.adj[u1].add(u2)
        self.adj[u2].add(u1)
        
    def bfs_expansion(self, seed_nodes: list, hops: int = 5) -> set:
        """Runs hops-radius BFS from seeds, returning all reached node IDs."""
        visited = set(seed_nodes)
        frontier = list(seed_nodes)
        for _ in range(hops):
            next_frontier = []
            for u in frontier:
                for v in self.adj[u]:
                    if v not in visited:
                        visited.add(v)
                        next_frontier.append(v)
            frontier = next_frontier
            if not frontier:
                break
        return visited

def fragment_target(seq: str, mode: str, k_dbg: int = 31, rng: random.Random = None) -> list:
    """
    Fragments target sequence into unitigs overlapping by (k_dbg - 1).
    Modes:
    - 'single': 1 unitig (~800 bp)
    - 'moderate': 4-5 unitigs (lengths 120-250 bp)
    - 'fragmented': 8-10 unitigs (lengths 60-140 bp)
    """
    L = len(seq)
    overlap = k_dbg - 1
    
    if mode == 'single':
        return [seq]
    elif mode == 'moderate':
        n_segs = 4
    elif mode == 'fragmented':
        n_segs = 8
    else:
        raise ValueError(f"Unknown mode: {mode}")
        
    min_len = overlap + 20
    step = (L - overlap) / n_segs
    cuts = [0]
    for i in range(1, n_segs):
        cut = int(i * step + rng.randint(-15, 15))
        cut = max(cuts[-1] + min_len, min(cut, L - (n_segs - i) * min_len))
        cuts.append(cut)
    cuts.append(L)
    
    unitigs = []
    for i in range(n_segs):
        start = cuts[i]
        seg_start = max(0, start - (overlap if i > 0 else 0))
        seg_end = cuts[i+1]
        unitigs.append(seq[seg_start:seg_end])
        
    return unitigs

# -----------------------------------------------------------------------------
# Fast Rolling Canonical k-mer Engine (Concordant with stages/anchor.rs)
# -----------------------------------------------------------------------------

def build_query_index(query_seq: str, k: int) -> set:
    """Builds set of canonical k-mers for query sequence."""
    n_windows = len(query_seq) - k + 1
    if n_windows <= 0:
        return set()
    kmers = set()
    for i in range(n_windows):
        sub = query_seq[i:i+k]
        kmers.add(canonical_kmer(sub))
    return kmers

def count_shared_kmers(unitig: str, query_kmers: set, k: int) -> int:
    """
    Counts matching canonical k-mer windows in unitig.
    Exactly mirrors `anchor::find`: counts matching positions.
    """
    if len(unitig) < k:
        return 0
    count = 0
    for i in range(len(unitig) - k + 1):
        window = unitig[i:i+k]
        if canonical_kmer(window) in query_kmers:
            count += 1
    return count

# -----------------------------------------------------------------------------
# Background Graph Construction
# -----------------------------------------------------------------------------

def build_background_graph(ecoli_seq: str, n_unitigs: int = 30000, rng: random.Random = None) -> UnitigGraphSim:
    """
    Constructs a realistic bacterial background graph of `n_unitigs`
    derived from E. coli K-12 with realistic de Bruijn graph connectivity.
    """
    sim = UnitigGraphSim()
    L_genome = len(ecoli_seq)
    
    lengths = np.random.lognormal(mean=5.0, sigma=0.4, size=n_unitigs).astype(int)
    lengths = np.clip(lengths, 45, 600)
    
    current_pos = 0
    overlap = 30
    for i in range(n_unitigs):
        length = lengths[i]
        if current_pos + length > L_genome:
            current_pos = rng.randint(0, 10000)
        seq = ecoli_seq[current_pos : current_pos + length]
        sim.add_unitig(seq, is_target=False)
        
        if i > 0 and rng.random() < 0.92:
            sim.add_edge(i - 1, i)
        if i > 10 and rng.random() < 0.08:
            rand_prev = rng.randint(max(0, i - 50), i - 1)
            sim.add_edge(rand_prev, i)
            
        current_pos += max(20, length - overlap)
        
    return sim

# -----------------------------------------------------------------------------
# Experiment Execution
# -----------------------------------------------------------------------------

def run_experiment(
    query_name: str,
    query_seq: str,
    bg_graph: UnitigGraphSim,
    k_grid: list,
    min_shared_grid: list,
    divergence_grid: list,
    mutation_models: list,
    fragmentation_modes: list,
    n_replicates: int = 30,
    seed: int = 42
):
    rng = random.Random(seed)
    np.random.seed(seed)
    
    print(f"=== Starting Experiment for {query_name} ({len(query_seq)} bp) ===")
    print(f"Background unitigs: {len(bg_graph.unitigs)}")
    print(f"Testing k in {k_grid}")
    print(f"Testing min_shared_kmers in {min_shared_grid}")
    print(f"Testing identities: {[round(100*(1-d), 1) for d in divergence_grid]}%")
    print(f"Replicates per cell: {n_replicates}")
    
    print("\nPre-scanning background graph for false positives...")
    bg_kmer_counts = {}
    for k in k_grid:
        t0 = time.time()
        q_kmers = build_query_index(query_seq, k)
        counts = [count_shared_kmers(u, q_kmers, k) for u in bg_graph.unitigs]
        bg_kmer_counts[k] = counts
        t1 = time.time()
        print(f"  k={k:2d}: query={len(q_kmers):4d} kmers | max bg matches={max(counts)} | time={t1-t0:.2f}s")
        
    results = []
    
    total_runs = len(k_grid) * len(min_shared_grid) * len(divergence_grid) * len(mutation_models) * len(fragmentation_modes)
    print(f"\nSimulating {total_runs} conditions x {n_replicates} replicates = {total_runs * n_replicates} trials...")
    
    trial_count = 0
    start_time = time.time()
    
    for k in k_grid:
        q_kmers = build_query_index(query_seq, k)
        bg_counts = bg_kmer_counts[k]
        
        for m in min_shared_grid:
            bg_fp_ids = [uid for uid, c in enumerate(bg_counts) if c >= m]
            n_fp_anchors = len(bg_fp_ids)
            
            if n_fp_anchors > 0:
                bg_subgraph_nodes = bg_graph.bfs_expansion(bg_fp_ids, hops=5)
                n_fp_subgraph = len(bg_subgraph_nodes)
            else:
                n_fp_subgraph = 0
                
            for div in divergence_grid:
                ident = round(1.0 - div, 4)
                
                for model in mutation_models:
                    for frag_mode in fragmentation_modes:
                        
                        target_anchors_detected = []
                        target_recalls = []
                        gene_recovered_flags = []
                        subgraph_target_recalls = []
                        total_subgraph_sizes = []
                        guard_violations = []
                        snrs = []
                        
                        for rep in range(n_replicates):
                            if model == 'uniform':
                                mutated = mutate_uniform(query_seq, div, rng)
                            elif model == 'codon':
                                mutated = mutate_codon_biased(query_seq, div, rng)
                            elif model == 'cluster':
                                mutated = mutate_clustered_motifs(query_seq, div, rng)
                            else:
                                raise ValueError(f"Unknown model: {model}")
                                
                            target_u_seqs = fragment_target(mutated, frag_mode, k_dbg=31, rng=rng)
                            n_target_u = len(target_u_seqs)
                            
                            target_counts = [count_shared_kmers(tu, q_kmers, k) for tu in target_u_seqs]
                            target_anchors = [i for i, c in enumerate(target_counts) if c >= m]
                            n_ta = len(target_anchors)
                            
                            recall = n_ta / n_target_u
                            target_recalls.append(recall)
                            target_anchors_detected.append(n_ta)
                            
                            gene_found = 1 if n_ta > 0 else 0
                            gene_recovered_flags.append(gene_found)
                            
                            if n_ta == 0:
                                recovered_target_nodes = set()
                            else:
                                recovered_target_nodes = set()
                                for a in target_anchors:
                                    for delta in range(-5, 6):
                                        idx = a + delta
                                        if 0 <= idx < n_target_u:
                                            recovered_target_nodes.add(idx)
                            subgraph_target_recall = len(recovered_target_nodes) / n_target_u
                            subgraph_target_recalls.append(subgraph_target_recall)
                            
                            total_subgraph = len(recovered_target_nodes) + n_fp_subgraph
                            total_subgraph_sizes.append(total_subgraph)
                            
                            total_graph_nodes = len(bg_graph.unitigs) + n_target_u
                            guard_violated = 1 if (total_subgraph / total_graph_nodes) > 0.01 else 0
                            guard_violations.append(guard_violated)
                            
                            snr = n_ta / (n_fp_anchors + 1e-9)
                            snrs.append(snr)
                            
                            trial_count += 1
                            
                        row = {
                            "gene": query_name,
                            "k": k,
                            "min_shared_kmers": m,
                            "divergence": div,
                            "identity": ident,
                            "identity_pct": round(ident * 100, 1),
                            "mutation_model": model,
                            "frag_mode": frag_mode,
                            "mean_target_recall": float(np.mean(target_recalls)),
                            "std_target_recall": float(np.std(target_recalls)),
                            "mean_target_anchors": float(np.mean(target_anchors_detected)),
                            "gene_recovery_rate": float(np.mean(gene_recovered_flags)),
                            "mean_subgraph_target_recall": float(np.mean(subgraph_target_recalls)),
                            "fp_anchors": n_fp_anchors,
                            "fp_subgraph_nodes": n_fp_subgraph,
                            "mean_total_subgraph": float(np.mean(total_subgraph_sizes)),
                            "guard_violation_rate": float(np.mean(guard_violations)),
                            "mean_snr": float(np.mean(snrs)),
                            "theory_kmer_match_prob": float((1.0 - div) ** k) if div < 1.0 else 0.0
                        }
                        results.append(row)
                        
    elapsed = time.time() - start_time
    print(f"\nCompleted {trial_count} trials in {elapsed:.2f}s ({trial_count/elapsed:.1f} trials/s)")
    return results

# -----------------------------------------------------------------------------
# Main Entry Point
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Run simulated k-mer sensitivity experiment.")
    parser.add_argument("--replicates", type=int, default=30, help="Number of replicates per condition")
    parser.add_argument("--bg-unitigs", type=int, default=30000, help="Number of background unitigs")
    parser.add_argument("--out-dir", type=str, default=str(HERE), help="Output directory")
    args = parser.parse_args()
    
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    vim_path = CDS_DIR / "bla_vim.fasta"
    kpc_path = CDS_DIR / "bla_kpc.fasta"
    
    vim_seq = load_fasta_seq(vim_path)
    kpc_seq = load_fasta_seq(kpc_path)
    
    print(f"Loaded bla_vim: {len(vim_seq)} bp")
    print(f"Loaded bla_kpc: {len(kpc_seq)} bp")
    
    if not ECOLI_PATH.exists():
        print(f"Warning: E. coli genome not found at {ECOLI_PATH}. Generating synthetic background.")
        rng = random.Random(1337)
        ecoli_seq = "".join(rng.choices("ACGT", weights=[0.24, 0.26, 0.26, 0.24], k=4600000))
    else:
        print(f"Loading real E. coli K-12 genome from {ECOLI_PATH}...")
        ecoli_seq = load_fasta_seq(ECOLI_PATH)
        print(f"E. coli genome: {len(ecoli_seq)} bp")
        
    print(f"Building background graph with {args.bg_unitigs} unitigs...")
    rng = random.Random(42)
    bg_graph = build_background_graph(ecoli_seq, n_unitigs=args.bg_unitigs, rng=rng)
    print(f"Background graph built: {len(bg_graph.unitigs)} unitigs, {sum(len(v) for v in bg_graph.adj.values())//2} edges.")
    
    k_grid = [15, 19, 21, 25, 27, 31]
    min_shared_grid = [1, 2, 3, 5, 8]
    divergence_grid = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40] # 100% to 60%
    mutation_models = ['uniform', 'codon', 'cluster']
    fragmentation_modes = ['single', 'moderate', 'fragmented']
    
    all_results = []
    
    res_vim = run_experiment(
        query_name="bla_vim",
        query_seq=vim_seq,
        bg_graph=bg_graph,
        k_grid=k_grid,
        min_shared_grid=min_shared_grid,
        divergence_grid=divergence_grid,
        mutation_models=mutation_models,
        fragmentation_modes=fragmentation_modes,
        n_replicates=args.replicates,
        seed=101
    )
    all_results.extend(res_vim)
    
    res_kpc = run_experiment(
        query_name="bla_kpc",
        query_seq=kpc_seq,
        bg_graph=bg_graph,
        k_grid=k_grid,
        min_shared_grid=min_shared_grid,
        divergence_grid=divergence_grid,
        mutation_models=mutation_models,
        fragmentation_modes=fragmentation_modes,
        n_replicates=args.replicates,
        seed=202
    )
    all_results.extend(res_kpc)
    
    json_path = out_dir / "simulation_results.json"
    csv_path = out_dir / "simulation_results.csv"
    
    with open(json_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved JSON results to {json_path}")
    
    df = pd.DataFrame(all_results)
    df.to_csv(csv_path, index=False)
    print(f"Saved CSV results to {csv_path} ({len(df)} rows)")

if __name__ == "__main__":
    main()
