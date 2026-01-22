#!/usr/bin/env python3
from __future__ import annotations

import argparse
from rich_argparse import RichHelpFormatter
from rich import print as rprint
from pathlib import Path
from typing import Iterable

import numpy as np
import matplotlib.pyplot as plt


def moi_curves(moi: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns (infected_total, infected_once, infected_multiple) as FRACTIONS in [0,1]
    using the standard Poisson MOI model.
    """
    infected_total = 1.0 - np.exp(-moi)                 # P(k>=1)
    infected_once = moi * np.exp(-moi)                  # P(k=1)
    infected_multiple = infected_total - infected_once  # P(k>=2)
    return infected_total, infected_once, infected_multiple


def parse_float_list(s: str) -> list[float]:
    """
    Accepts comma/space-separated floats: "0.1,0.3, 1.0"
    """
    if s is None:
        return []
    parts = [p.strip() for p in s.replace(";", ",").replace(" ", ",").split(",")]
    parts = [p for p in parts if p]
    out: list[float] = []
    for p in parts:
        try:
            out.append(float(p))
        except ValueError as e:
            raise argparse.ArgumentTypeError(f"Could not parse float from '{p}'") from e
    return out


def ensure_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def save_figure(fig: plt.Figure, out: Path, dpi: int = 300) -> None:
    ensure_dir(out)
    # dpi only matters for raster formats; matplotlib ignores it for vector formats like pdf/svg.
    fig.savefig(out, bbox_inches="tight", dpi=dpi)


def plot_moi(
    *,
    max_moi: float = 5.0,
    n_points: int = 500,
    vlines: Iterable[float] = (),
    percent: bool = True,
    title: str = "MOI vs Cell Infection Outcomes",
    show: bool = False,
    out_files: list[Path] | None = None,
    dpi: int = 300,
) -> None:
    if max_moi <= 0:
        raise ValueError("--max-moi must be > 0")
    if n_points < 2:
        raise ValueError("--n-points must be >= 2")

    moi = np.linspace(0, max_moi, n_points)
    infected_total, infected_once, infected_multiple = moi_curves(moi)

    # Convert to %
    if percent:
        infected_total *= 100.0
        infected_once *= 100.0
        infected_multiple *= 100.0
        ylab = "% of cells"
        ylim = (0, 100)
    else:
        ylab = "Fraction of cells"
        ylim = (0, 1)

    fig, ax = plt.subplots(figsize=(7.0, 4.6))

    ax.plot(moi, infected_total, label="Infected (≥1)", linewidth=2)
    ax.plot(moi, infected_once, label="Infected once", linewidth=2)
    ax.plot(moi, infected_multiple, label="Infected >1", linewidth=2)

    # Vertical reference lines
    for x in vlines:
        if x < 0:
            continue
        ax.axvline(x=x, linestyle="--", linewidth=1)
        # small label near top
        ax.text(x, ylim[1] * 0.98, f"{x:g}", rotation=90, va="top", ha="right")

    ax.set_xlabel("MOI")
    ax.set_ylabel(ylab)
    ax.set_title(title)
    ax.set_xlim(0, max_moi)
    ax.set_ylim(*ylim)

    # Minor ticks for decimal resolution
    from matplotlib.ticker import AutoMinorLocator

    # X-axis: decimal MOI minor ticks
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    # Y-axis: finer resolution (percent or fraction)
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    ax.tick_params(axis="both", which="minor", length=3, width=0.6)

    ax.grid(True, linewidth=0.5, alpha=0.4)
    ax.legend(frameon=False)

    fig.tight_layout()

    # Save outputs
    if out_files:
        for out in out_files:
            save_figure(fig, out, dpi=dpi)

    if show:
        plt.show()

    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
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
        description="Plot MOI Poisson infection curves: total infected, infected once, infected >1.",
        formatter_class=MyFormatter
    )
    p.add_argument("-m", "--max-moi", type=float, default=5.0, help="Maximum MOI to plot (default: 5).")
    p.add_argument("-n", "--n-points", type=int, default=500, help="Number of x points (default: 500).")

    p.add_argument(
        "-v", "--vlines",
        type=parse_float_list,
        default=[],
        help='Comma/space-separated MOIs for vertical dashed lines, e.g. "0.1,0.3,1".',
    )

    p.add_argument("-p", "--percent", action="store_true", help="Use percent y-axis (0–100). (default)")
    p.add_argument("-f", "--fraction", dest="percent", action="store_false", help="Use fraction y-axis (0–1).")
    p.set_defaults(percent=True)

    p.add_argument("-t", "--title", type=str, default="MOI vs Cell Infection Outcomes", help="Plot title.")
    p.add_argument("-d", "--dpi", type=int, default=300, help="DPI for PNG output (default: 300).")

    # Outputs: allow multiple
    p.add_argument(
        "-o", "--out",
        action="append",
        default=[],
        help="Output file path. Can be provided multiple times. Supports .pdf .svg .png",
    )

    p.add_argument("-s", "--show", action="store_true", help="Display interactively.")
    return p


def main() -> None:
    p = build_parser()
    args = p.parse_args()

    out_files = [Path(x) for x in (args.out or [])]
    plot_moi(
        max_moi=args.max_moi,
        n_points=args.n_points,
        vlines=args.vlines,
        percent=args.percent,
        title=args.title,
        show=args.show,
        out_files=out_files if out_files else None,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()