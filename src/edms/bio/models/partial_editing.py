#!/usr/bin/env python3
"""
plot_editing_compression.py

Plot |log2FC_obs| = |log2(1 + e*(FC_edit - 1))|
as a function of editing fraction e and edited-only |log2FC|.

Y-axis: edited-only |log2FC|
Color:  observed  |log2FC|
Optional: vertical dashed lines at specified editing fractions.
"""

from __future__ import annotations

import argparse
from typing import Sequence, Optional
from rich_argparse import RichHelpFormatter
import numpy as np
import matplotlib.pyplot as plt


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
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

    p = argparse.ArgumentParser(
        description="Visualize deterministic compression of log2FC by partial editing.",
        formatter_class=MyFormatter,
    )

    p.add_argument("-m",
        "--max_log2FC_edit_mag",
        type=float,
        default=2.0,
        help="Maximum edited-only |log2FC| shown on y-axis (default: 2.0).",
    )
    p.add_argument(
        "-c",
        "--cmap",
        type=str,
        default="Reds",
        help="Matplotlib colormap name (default: Reds).",
    )
    p.add_argument(
        "-ef",
        "--editing_fractions",
        type=float,
        nargs="*",
        default=[],
        help="One or more editing fractions e in [0,1] to draw as vertical dashed lines.",
    )
    p.add_argument(
        "-el",
        "--edited_only_log2FCs",
        type=float,
        nargs="*",
        default=[],
        help="One or more edited-only |log2FC| values to draw as horizontal dashed lines.",
    )
    p.add_argument(
        "-ol",
        "--observed_log2FCs",
        type=float,
        nargs="*",
        default=[],
        help="One or more |observed log2FC| values to draw as curved contour lines.",
    )

    # Nice-to-have options
    p.add_argument(
        "-gs",
        "--grid_size",
        type=int,
        default=300,
        help="Grid resolution for the heatmap (default: 300).",
    )
    p.add_argument(
        "-fs",
        "--figsize",
        type=float,
        nargs=2,
        default=(7.5, 5.5),
        metavar=("W", "H"),
        help="Figure size in inches (default: 7.5 5.5).",
    )
    p.add_argument(
        "-t",
        "--title",
        type=str,
        default="Deterministic Compression of log2FC by Partial Editing",
        help="Plot title.",
    )
    p.add_argument(
        "-o",
        "--out",
        type=str,
        default=None,
        help="Output file path (e.g. plot.png). If omitted, shows interactively.",
    )
    p.add_argument(
        "-d",
        "--dpi",
        type=int,
        default=600,
        help="DPI for saved figure (default: 600).",
    )

    args = p.parse_args(argv)

    # Validate
    if args.max_log2FC_edit_mag <= 0:
        p.error("--max_log2FC_edit_mag must be > 0")

    bad = [x for x in args.editing_fractions if (x < 0 or x > 1)]
    if bad:
        p.error(f"--editing_fractions must be in [0,1]. Invalid: {bad}")

    bad = [x for x in args.edited_only_log2FCs if x < 0]
    if bad:
        p.error(f"--edited_only_log2FCs must be >= 0. Invalid: {bad}")

    bad = [x for x in args.observed_log2FCs if x < 0]
    if bad:
        p.error(f"--observed_log2FCs must be >= 0. Invalid: {bad}")

    if args.grid_size < 50:
        p.error("--grid_size should be >= 50 for a reasonable plot")

    return args


def make_plot(
    max_log2fc_edit_mag: float = 2.0,
    cmap: str = "Reds",
    editing_fractions: Sequence[float] = (),
    edited_only_log2fcs: Sequence[float] = (),
    observed_log2fcs: Sequence[float] = (),
    grid_size: int = 300,
    figsize: tuple[float, float] = (7.5, 5.5),
    title: str = "Deterministic Compression of log2FC by Partial Editing",
) -> plt.Figure:
    # Axes ranges
    e_vals = np.linspace(0, 1, grid_size)
    log2fc_edit_mag = np.linspace(0, max_log2fc_edit_mag, grid_size)

    # Convert edited-only |log2FC| magnitude to FC (>1); sign is symmetric
    fc_edit = 2 ** log2fc_edit_mag

    # Grid
    E, FC = np.meshgrid(e_vals, fc_edit)

    # Mixture model: observed log2FC
    # log2FC_obs = log2(1 + e*(FC_edit - 1))
    log2fc_obs = np.log2(1 + E * (FC - 1))

    # Plot
    fig = plt.figure(figsize=figsize)
    ax = plt.gca()

    im = ax.imshow(
        np.abs(log2fc_obs),
        origin="lower",
        aspect="auto",
        extent=[e_vals.min(), e_vals.max(), log2fc_edit_mag.min(), log2fc_edit_mag.max()],
        cmap=cmap,
    )

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("|Observed log2FC| = |log2(1 + e*(FC_edit - 1))|")

    ax.set_xlabel("Editing fraction (e)")
    ax.set_ylabel("|Edited-only log2FC|")
    ax.set_title(title)

    # Vertical dashed lines at specified editing fractions
    # (no explicit color set, per your request; matplotlib will choose default cycle color)
    for e in editing_fractions:
        ax.axvline(e, linestyle="--", linewidth=1.5, color="black")

    # Horizontal dashed lines at specified edited-only |log2FC| values
    for y in edited_only_log2fcs:
        if y < 0:
            continue
        ax.axhline(y, linestyle="--", linewidth=1.5, color="black")

    # Curved contour lines where |observed log2FC| equals specified value(s)
    obs_levels = [v for v in observed_log2fcs if v >= 0]
    if obs_levels:
        # Build contour grid in axis units
        X, Y = np.meshgrid(e_vals, log2fc_edit_mag)
        Z = np.abs(log2fc_obs)
        cs = ax.contour(X, Y, Z, levels=sorted(obs_levels), colors="black", linewidths=1.5)
        ax.clabel(cs, inline=True, fontsize=8, fmt=lambda v: f"{v:g}")

    fig.tight_layout()
    return fig


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)

    fig = make_plot(
        max_log2fc_edit_mag=args.max_log2FC_edit_mag,
        cmap=args.cmap,
        editing_fractions=args.editing_fractions,
        edited_only_log2fcs=args.edited_only_log2FCs,
        observed_log2fcs=args.observed_log2FCs,
        grid_size=args.grid_size,
        figsize=tuple(args.figsize),
        title=args.title,
    )

    if args.out:
        fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    else:
        plt.show()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())