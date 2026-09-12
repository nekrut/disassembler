#!/usr/bin/env python3
"""
Generate publication-quality figures and summary statistics
from simulated k-mer sensitivity experiment results.
"""

import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

HERE = Path(__file__).resolve().parent
CSV_PATH = HERE / "simulation_results.csv"
FIGS_DIR = HERE / "figures"
FIGS_DIR.mkdir(parents=True, exist_ok=True)

def set_style():
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    plt.rcParams.update({
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 14,
        'lines.linewidth': 2.0,
        'lines.markersize': 6,
        'grid.alpha': 0.4,
    })

def plot_sensitivity_vs_identity(df: pd.DataFrame):
    """Plot Anchor Recall vs Sequence Identity across k and min_shared_kmers."""
    # Filter: bla_kpc, moderate fragmentation, uniform vs codon
    sub = df[(df['gene'] == 'bla_kpc') & (df['frag_mode'] == 'moderate')]
    
    fig, axes = plt.subplots(2, 3, figsize=(16, 10), sharex=True, sharey=True)
    axes = axes.flatten()
    
    m_values = [1, 2, 3, 5, 8]
    k_values = [15, 19, 21, 25, 27, 31]
    colors = {15: '#2ca02c', 19: '#1f77b4', 21: '#ff7f0e', 25: '#9467bd', 27: '#8c564b', 31: '#d62728'}
    
    for ax_idx, m in enumerate(m_values):
        ax = axes[ax_idx]
        for k in k_values:
            cell = sub[(sub['k'] == k) & (sub['min_shared_kmers'] == m) & (sub['mutation_model'] == 'uniform')]
            cell = cell.sort_values('identity_pct', ascending=False)
            ax.plot(cell['identity_pct'], cell['mean_target_recall'] * 100, 
                    marker='o', color=colors[k], label=f'k={k}')
            
        ax.set_title(f"min_shared_kmers = {m}")
        ax.set_xlabel("Sequence Identity (%)")
        ax.set_ylabel("Anchor Recall / Sensitivity (%)")
        ax.set_xlim(100, 60)
        ax.set_ylim(-2, 105)
        ax.grid(True, linestyle='--', alpha=0.5)
        if ax_idx == 0:
            ax.legend(loc='lower left', frameon=True)
            
    # Hide the 6th subplot
    axes[5].axis('off')
    
    fig.suptitle("Anchor Sensitivity vs Sequence Identity (bla_kpc, Moderate Fragmentation)", fontsize=15, y=0.98)
    plt.tight_layout()
    out_file = FIGS_DIR / "anchor_sensitivity_vs_identity.png"
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"Saved {out_file}")

def plot_gene_recovery_and_bfs_bridge(df: pd.DataFrame):
    """Plot Gene Recovery (BFS 5-hop) vs Anchor Recall, demonstrating bridging power."""
    sub = df[(df['gene'] == 'bla_kpc') & (df['frag_mode'] == 'moderate') & (df['min_shared_kmers'] == 3)]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    
    k_values = [15, 19, 21, 25, 27, 31]
    colors = {15: '#2ca02c', 19: '#1f77b4', 21: '#ff7f0e', 25: '#9467bd', 27: '#8c564b', 31: '#d62728'}
    
    # Panel 1: Raw Anchor Recall
    for k in k_values:
        cell = sub[(sub['k'] == k) & (sub['mutation_model'] == 'uniform')].sort_values('identity_pct', ascending=False)
        ax1.plot(cell['identity_pct'], cell['mean_target_recall'] * 100,
                 marker='s', color=colors[k], label=f'k={k}')
    ax1.set_title("Step 3: Direct Anchor Recall (min_shared=3)")
    ax1.set_xlabel("Sequence Identity (%)")
    ax1.set_ylabel("Fraction (%)")
    ax1.set_xlim(100, 60)
    ax1.set_ylim(-2, 105)
    ax1.legend(loc='lower left')
    ax1.grid(True, linestyle='--', alpha=0.5)
    
    # Panel 2: Full Gene Recovery via 5-hop BFS
    for k in k_values:
        cell = sub[(sub['k'] == k) & (sub['mutation_model'] == 'uniform')].sort_values('identity_pct', ascending=False)
        ax2.plot(cell['identity_pct'], cell['mean_subgraph_target_recall'] * 100,
                 marker='o', color=colors[k], label=f'k={k}')
    ax2.set_title("Step 4: Subgraph Target Coverage (After 5-Hop BFS)")
    ax2.set_xlabel("Sequence Identity (%)")
    ax2.set_xlim(100, 60)
    ax2.grid(True, linestyle='--', alpha=0.5)
    
    fig.suptitle("The 5-Hop BFS Bridging Effect: Subgraph Expansion Recovers Partially Anchored Genes", fontsize=14)
    plt.tight_layout()
    out_file = FIGS_DIR / "gene_recovery_bfs_bridging.png"
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"Saved {out_file}")

def plot_false_positives_and_subgraph_blowup(df: pd.DataFrame):
    """Plot False Positive Anchors and Subgraph Node Expansion vs k and min_shared_kmers."""
    # Dedup background FP rows (gene independent since bg is constant)
    sub = df[(df['gene'] == 'bla_kpc') & (df['identity_pct'] == 100.0) & 
             (df['mutation_model'] == 'uniform') & (df['frag_mode'] == 'moderate')].copy()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    m_values = [1, 2, 3, 5, 8]
    k_values = [15, 19, 21, 25, 27, 31]
    markers = {1: 'o', 2: 's', 3: '^', 5: 'd', 8: 'v'}
    
    for m in m_values:
        cell = sub[sub['min_shared_kmers'] == m].sort_values('k')
        ax1.plot(cell['k'], cell['fp_anchors'], marker=markers[m], label=f'min_shared={m}')
        ax2.plot(cell['k'], cell['fp_subgraph_nodes'], marker=markers[m], label=f'min_shared={m}')
        
    ax1.set_title("False Positive Anchors in 30,000 Unitig Background")
    ax1.set_xlabel("k-mer size (k)")
    ax1.set_ylabel("False Positive Unitigs (Count)")
    ax1.set_xticks(k_values)
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.5)
    
    ax2.set_title("5-Hop Expanded False Subgraph Nodes (Noise Snowball)")
    ax2.set_xlabel("k-mer size (k)")
    ax2.set_ylabel("Background Nodes in Subgraph")
    ax2.set_xticks(k_values)
    # Draw 1% guard line (30,000 * 0.01 = 300 nodes)
    ax2.axhline(300, color='red', linestyle='--', linewidth=1.5, label='1% Subgraph Guard (300 nodes)')
    ax2.legend()
    ax2.grid(True, linestyle='--', alpha=0.5)
    
    fig.suptitle("Specificity & False Positive Subgraph Explosion in E. coli Background (4.6 Mb)", fontsize=14)
    plt.tight_layout()
    out_file = FIGS_DIR / "false_positives_and_subgraph_expansion.png"
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"Saved {out_file}")

def plot_theory_vs_empirical(df: pd.DataFrame):
    """Plot theoretical P(k-mer match) = (1-d)^k vs empirical anchor retention across mutation models."""
    sub = df[(df['gene'] == 'bla_kpc') & (df['frag_mode'] == 'single') & (df['min_shared_kmers'] == 1)]
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
    
    models = ['uniform', 'codon', 'cluster']
    titles = ['Uniform Point Mutations', 'Codon-Biased (70% 3rd, Ti/Tv=2)', 'Domain Motifs Conserved (Active Site)']
    
    div_continuous = np.linspace(0.0, 0.40, 100)
    k_to_plot = [15, 21, 31]
    colors = {15: '#2ca02c', 21: '#ff7f0e', 31: '#d62728'}
    
    for ax_idx, (model, title) in enumerate(zip(models, titles)):
        ax = axes[ax_idx]
        
        # Plot theoretical lines
        for k in k_to_plot:
            theory_p = (1.0 - div_continuous) ** k
            ax.plot((1.0 - div_continuous) * 100, theory_p * 100, '--', color=colors[k], alpha=0.6,
                    label=f'Theory k={k}: (1-d)^{k}' if ax_idx == 0 else "")
            
        # Plot empirical
        for k in k_to_plot:
            cell = sub[(sub['k'] == k) & (sub['mutation_model'] == model)].sort_values('identity_pct', ascending=False)
            ax.plot(cell['identity_pct'], cell['mean_target_recall'] * 100, 'o-', color=colors[k],
                    label=f'Empirical k={k}' if ax_idx == 0 else "")
            
        ax.set_title(title)
        ax.set_xlabel("Sequence Identity (%)")
        if ax_idx == 0:
            ax.set_ylabel("k-mer / Anchor Retention (%)")
            ax.legend(loc='lower left')
        ax.set_xlim(100, 60)
        ax.grid(True, linestyle='--', alpha=0.5)
        
    fig.suptitle("Theoretical Limit (1 - divergence)^k vs Empirical Anchor Retention", fontsize=14)
    plt.tight_layout()
    out_file = FIGS_DIR / "theoretical_vs_empirical_decay.png"
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"Saved {out_file}")

def main():
    if not CSV_PATH.exists():
        print(f"Error: {CSV_PATH} not found. Run harness.py first.")
        sys.exit(1)
        
    set_style()
    df = pd.read_csv(CSV_PATH)
    print(f"Loaded {len(df)} simulation records from {CSV_PATH}")
    
    plot_sensitivity_vs_identity(df)
    plot_gene_recovery_and_bfs_bridge(df)
    plot_false_positives_and_subgraph_blowup(df)
    plot_theory_vs_empirical(df)
    print("All figures successfully generated.")

if __name__ == "__main__":
    main()
