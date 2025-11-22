#!/usr/bin/env python3
import argparse
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from rich_argparse import RichHelpFormatter

# -----------------------------
# Helper functions
# -----------------------------
def clean_seq(seq: str) -> str:
    """Return uppercase sequence containing only A/C/G/T."""
    return ''.join(b for b in seq.upper() if b in "ACGT")

def base_index(b: str) -> int:
    return {'A':0, 'C':1, 'G':2, 'T':3}[b]

def per_cycle_base_freqs(seq: str, S: int, max_cycles: int) -> np.ndarray:
    """Expected base fractions (A,C,G,T) for cycles 1..max_cycles, pooling stagger s=0..S."""
    freqs = np.zeros((max_cycles, 4), dtype=float)
    for c in range(max_cycles):
        counts = np.zeros(4, dtype=float)
        for s in range(S + 1):
            if c <= s:
                counts += 0.25
            else:
                pos = c - s
                if 0 <= pos < len(seq):
                    counts[base_index(seq[pos])] += 1.0
        freqs[c, :] = counts / (S + 1)
    return freqs

def shannon_entropy_bits(p: np.ndarray) -> np.ndarray:
    """H = -Σ p_i log2 p_i (bits)."""
    eps = 1e-12
    return -np.sum(p * np.log2(p + eps), axis=-1)

# -----------------------------
# Main
# -----------------------------
def main():
    # Custom formatter for rich help messages
    class MyFormatter(RichHelpFormatter):
        styles = {
            "argparse.prog": "green",           # program name
            "argparse.args": "cyan",            # positional arguments
            "argparse.option": "",              # options like --flag
            "argparse.metavar": "dark_magenta", # meta variable (actual function argument name)
            "argparse.help": "blue",            # help text
            "argparse.text": "green",           # normal text in help message
            "argparse.groups": "red",           # group titles
            "argparse.description": "",         # description at the top
            "argparse.epilog": "",              # ... -h; epilog at the bottom
            "argparse.syntax": "white",         # []
        }
    parser = argparse.ArgumentParser(
        description="Plot expected per-cycle nucleotide diversity (Shannon entropy) for staggered forward primers.",
        formatter_class=MyFormatter
    )
    parser.add_argument("--seq", required=True,
                        help="Sequence string or path to FASTA/plain-text file.")
    parser.add_argument("--S_values", nargs="+", type=int, default=[0, 2, 4, 6, 8],
                        help="List of max stagger values (e.g. 0 2 4 6 8). Each pool includes staggers 0..S inclusive.")
    parser.add_argument("--max_cycles", type=int, default=None,
                        help="Number of sequencing cycles. Defaults to full sequence length.")
    parser.add_argument("--out_csv", type=str, default=None,
                        help="Optional path to save CSV with per-cycle stats.")
    parser.add_argument("--out_plot", type=str, default=None,
                        help="Optional path to save the plot image file.")
    parser.add_argument("--out_facet", type=str, default=None,
                        help="Optional path to save the base-fractions facet plot (A/C/G/T).")
    args = parser.parse_args()

    seq = clean_seq(args.seq)
    if len(seq) == 0:
        raise ValueError("Sequence appears empty or invalid.")
    max_cycles = args.max_cycles or len(seq)

    records = []
    for S in args.S_values:
        freqs = per_cycle_base_freqs(seq, S, max_cycles)
        H = shannon_entropy_bits(freqs)
        for c in range(max_cycles):
            records.append({
                "cycle": c + 1,
                "S_max": S,
                "entropy_bits": float(H[c]),
                "A": float(freqs[c, 0]),
                "C": float(freqs[c, 1]),
                "G": float(freqs[c, 2]),
                "T": float(freqs[c, 3]),
                "min_base_fraction": float(freqs[c].min()),
            })

    df = pd.DataFrame.from_records(records)

    # Save CSV if requested
    if args.out_csv:
        df.to_csv(args.out_csv, index=False)
        print(f"Saved results to {args.out_csv}")

    # Plot
    plt.figure(figsize=(9, 5))
    for S in args.S_values:
        d = df[df["S_max"] == S]
        plt.plot(d["cycle"], d["entropy_bits"], label=f"S={S}")
    plt.xlabel("Sequencing cycle")
    plt.ylabel("Shannon entropy (bits)")
    plt.title("Per-cycle nucleotide diversity with forward-primer stagger (0..S Ns)")
    plt.legend()
    plt.tight_layout()
    if args.out_plot:
        plt.savefig(args.out_plot, dpi=300)
        print(f"Saved plot to {args.out_plot}")

    # Facet grid plot: facet by S (stagger); within each facet plot A, C, G, T fractions together
    unique_S = list(dict.fromkeys(args.S_values))  # preserve CLI order
    unique_S.remove(0)  # skip S=0 for base-fraction facets
    n = len(unique_S)
    ncols = min(2, n)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 3.5*nrows), sharex=True, sharey=True)
    # Normalize axes to 2D array
    if isinstance(axes, np.ndarray):
        axes = np.atleast_2d(axes)
    else:
        axes = np.array([[axes]])

    # Plot per S facet
    legend_handles, legend_labels = None, None
    for idx, S in enumerate(unique_S):
        r, c = divmod(idx, ncols)
        ax = axes[r, c]
        d = df[df["S_max"] == S]
        lineA, = ax.plot(d["cycle"], d["A"], label="A")
        lineC, = ax.plot(d["cycle"], d["C"], label="C")
        lineG, = ax.plot(d["cycle"], d["G"], label="G")
        lineT, = ax.plot(d["cycle"], d["T"], label="T")
        ax.set_title(f"S = {S}")
        ax.grid(False)
        if legend_handles is None:
            legend_handles = [lineA, lineC, lineG, lineT]
            legend_labels = ["A", "C", "G", "T"]

    # Hide any unused subplots
    total_axes = nrows * ncols
    for j in range(n, total_axes):
        r, c = divmod(j, ncols)
        axes[r, c].axis('off')

    # Set x-labels on bottom row
    for c in range(ncols):
        axes[nrows-1, c].set_xlabel("Sequencing cycle")
    
    # Set y-labels on left column
    for r in range(nrows):
        axes[r, 0].set_ylabel("Base fraction")

    # Single legend above all facets
    if legend_handles is not None:
        fig.legend(legend_handles, legend_labels, loc="lower center", ncol=4)

    fig.suptitle("Per-cycle base fractions by forward-primer stagger")
    fig.tight_layout(rect=(0, 0, 1, 1))
    if args.out_facet:
        fig.savefig(args.out_facet, dpi=300)
        print(f"Saved facet plot to {args.out_facet}")
    plt.show()

if __name__ == "__main__":
    main()