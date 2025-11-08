#!/usr/bin/env python3
"""
Bias-factor curves for clustered editing.

Model:
- Let N be edits per cell. Start from Poisson with lambda λ = -ln(1 - p_any),
  where p_any is the overall per-cell probability of ≥1 edit.
- For k >= 1, up-weight the mass by b^(k-1) (b >= 1). This captures that a
  cell with one edit is b× more likely to gain an additional edit than the
  Poisson baseline would predict, and similarly for higher k.
- Keep P(N=0) = 1 - p_any, and renormalize the k>=1 mass to sum to p_any.

We then compute P(N=2) and report P(two edits on different chromatids) = 0.5*P(N=2)
under symmetric chromatids.
"""

from __future__ import annotations
import math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from rich_argparse import RichHelpFormatter
from rich import print as rprint


# ---------- Core helpers ----------
def poisson_pmf(lam: float, kmax: int) -> np.ndarray:
    """PMF for Poisson(lam) for k=0..kmax."""
    ks = np.arange(kmax + 1)
    pmf = np.exp(-lam) * np.power(lam, ks) / np.array([math.factorial(k) for k in ks])
    return pmf


def biased_pmf_from_p_any(p_any: float, b: float, kmax: int = 20) -> np.ndarray:
    """
    Build the overdispersed PMF with bias factor b at a fixed overall per-cell
    edit probability p_any = P(N>=1).
    - Baseline λ = -ln(1 - p_any)
    - For k>=1, weight Poisson mass by b^(k-1).
    - Set P0 = 1 - p_any, renormalize the rest to sum to p_any.
    """
    if not (0.0 < p_any < 1.0):
        raise ValueError("p_any must be in (0, 1).")
    if b < 1.0:
        raise ValueError("b should be ≥ 1 for 'more likely' clustering (b=1 is Poisson).")

    lam = -math.log(1.0 - p_any)
    base = poisson_pmf(lam, kmax)

    weights = np.zeros_like(base)
    # weight k>=1 by b^(k-1)
    for k in range(1, kmax + 1):
        weights[k] = b ** (k - 1)

    weighted = base * weights
    S = weighted[1:].sum()

    pmf = np.zeros_like(base)
    pmf[0] = 1.0 - p_any
    pmf[1:] = (p_any * weighted[1:] / S) if S > 0 else 0.0
    return pmf


def p_two_diff(p_any: float, b: float, kmax: int = 20) -> float:
    """Return P(exactly two edits on different chromatids) under the biased PMF."""
    pmf = biased_pmf_from_p_any(p_any, b, kmax=kmax)
    return 0.5 * float(pmf[2])  # split equally across chromatids under symmetry


# ---------- Plotting ----------
def plot_p_two_diff_vs_p_any(
    p_any_grid: np.ndarray,
    b_list: list[float],
    kmax: int = 20,
    out_png: str | Path = "p_two_diff_vs_pany_multi_bias.png",
    out_csv: str | Path = "p_two_diff_vs_pany_multi_bias.csv",
) -> tuple[Path, Path]:
    rows = []
    plt.figure(figsize=(8, 5))
    for b in b_list:
        ys = [p_two_diff(p_any, b, kmax=kmax) for p_any in p_any_grid]
        plt.plot(p_any_grid, ys, marker="o", label=f"b={b:g}")
        rows.extend({"per_cell_edit_prob": x, "bias_b": b, "p_two_diff": y} for x, y in zip(p_any_grid, ys))

    plt.xlabel("Per-cell edit probability (≥1 mutation)")
    plt.ylabel("P(two edits on different chromatids)")
    plt.title("Effect of clustering across overall editing rates")
    plt.legend(title="Bias factor b")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()

    out_png = Path(out_png)
    out_csv = Path(out_csv)
    plt.savefig(out_png, dpi=150)
    plt.close()

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    return out_png, out_csv

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
        description="Generate bias-factor curves for clustered editing and save plot + CSV.",
        formatter_class=MyFormatter,
    )
    parser.add_argument(
        "--pmin", type=float, default=0.01,
        help="Minimum per-cell edit probability (default: 0.01)",
    )
    parser.add_argument(
        "--pmax", type=float, default=0.50,
        help="Maximum per-cell edit probability (default: 0.50)",
    )
    parser.add_argument(
        "--points", type=int, default=50,
        help="Number of points in the per-cell edit probability grid (default: 50)",
    )
    parser.add_argument(
        "--bias-list", type=float, nargs="+", default=[1.0, 1.5, 2.0, 3.0, 4.0],
        help="List of bias factors b to plot (default: 1.0 1.5 2.0 3.0 4.0)",
    )
    parser.add_argument(
        "--kmax", type=int, default=20, help="Maximum k to include in PMF (default: 20)"
    )
    parser.add_argument(
        "--out-png", type=str, default="p_two_diff_vs_pany_multi_bias.png",
        help="Output PNG file path (default: p_two_diff_vs_pany_multi_bias.png)",
    )
    parser.add_argument(
        "--out-csv", type=str, default="p_two_diff_vs_pany_multi_bias.csv",
        help="Output CSV file path (default: p_two_diff_vs_pany_multi_bias.csv)",
    )

    args = parser.parse_args()

    # Sweep overall per-cell editing efficiency
    p_any_grid = np.linspace(args.pmin, args.pmax, args.points)

    # Run the plot generation
    png, csv = plot_p_two_diff_vs_p_any(
        p_any_grid=p_any_grid,
        b_list=args.bias_list,
        kmax=args.kmax,
        out_png=args.out_png,
        out_csv=args.out_csv,
    )

    print(f"Saved plot to {png}")
    print(f"Saved CSV to {csv}")

if __name__ == "__main__":
    main()