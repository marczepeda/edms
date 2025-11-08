import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from rich_argparse import RichHelpFormatter
from rich import print as rprint


def build_parser():
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
        description="Plot expected errors (EE) vs average Phred Q for multiple read lengths.",
        formatter_class=MyFormatter
    )
    p.add_argument("--qmin", type=float, default=10, help="Minimum Q value (default: 10)")
    p.add_argument("--qmax", type=float, default=45, help="Maximum Q value (default: 45)")
    p.add_argument("--qstep", type=float, default=0.5, help="Step size for Q grid (default: 0.5)")
    p.add_argument(
        "--cycles", type=int, nargs="+", default=[100, 200, 300, 400, 500, 600],
        help="Read lengths (cycles) to plot; space-separated list (default: 100 200 300 400 500 600)"
    )
    p.add_argument(
        "--ee-lines", type=float, nargs="+", default=[1, 2, 3, 4, 5],
        help="Vertical EE reference lines (default: 1 2 3 4 5)"
    )
    p.add_argument(
        "--q-lines", type=float, nargs="+", default=[40, 35, 30, 25, 20],
        help="Horizontal Q reference lines (default: 40 35 30 25 20)"
    )
    p.add_argument("--width", type=float, default=9, help="Figure width in inches (default: 9)")
    p.add_argument("--height", type=float, default=2.5, help="Figure height in inches (default: 2.5)")
    p.add_argument("--title", type=str, default="Expected errors vs. Average Phred Q for different read lengths",
                   help="Plot title")
    p.add_argument("--no-legend", action="store_true", help="Disable legend")
    p.add_argument("--no-show", action="store_true", help="Do not display plot window")
    p.add_argument("--csv-out", type=str, default="", help="Optional path to write CSV (EE vs Q table)")
    p.add_argument("--png-out", type=str, default="", help="Optional path to save the plot as PNG")
    return p


def p_error_from_Q(Q):
    """Per-base error probability from Phred Q."""
    return 10 ** (-Q / 10.0)


def expected_errors_from_avgQ(Q, n_cycles):
    """EE using a single average Q for the whole read."""
    return n_cycles * p_error_from_Q(Q)


def expected_errors_from_percycle_Q(Q_list):
    """EE from a list/array of per-cycle Q values."""
    Q_arr = np.asarray(Q_list, dtype=float)
    return np.sum(10 ** (-Q_arr / 10.0))


def main():
    args = build_parser().parse_args()

    # Build dense Q grid and EE columns
    Q_vals = np.arange(args.qmin, args.qmax + args.qstep, args.qstep)
    df = pd.DataFrame({"Q": Q_vals})
    df["p_error_per_base"] = p_error_from_Q(df["Q"])

    for n in args.cycles:
        df[f"EE_n{n}"] = expected_errors_from_avgQ(df["Q"], n)

    # Optional CSV output (drop p column to keep it similar to your previous export)
    if args.csv_out:
        df_out = df.drop(columns=["p_error_per_base"]) if "p_error_per_base" in df.columns else df.copy()
        df_out.to_csv(args.csv_out, index=False)

    # Plot EE vs Q for each read length
    plt.figure(figsize=(args.width, args.height))
    for n in args.cycles:
        plt.plot(df[f"EE_n{n}"], df["Q"], label=f"{n} cycles")

    # EE reference lines
    for ee in args.ee_lines:
        plt.axvline(ee, linestyle="--", linewidth=1)

    # Q reference lines
    for q in args.q_lines:
        plt.axhline(q, linestyle="--", linewidth=1)

    plt.xscale("log")
    plt.ylabel("Average Phred Q")
    plt.xlabel("Expected errors per read (EE)")
    plt.title(args.title)
    if not args.no_legend:
        plt.legend(title="Read length", bbox_to_anchor=(1.01, 1), loc='upper left')
    plt.tight_layout()

    if args.png_out:
        plt.savefig(args.png_out, dpi=300, bbox_inches="tight")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()