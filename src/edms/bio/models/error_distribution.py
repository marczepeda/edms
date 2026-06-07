#!/usr/bin/env python

import argparse
from rich_argparse import RichHelpFormatter
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def error_distribution_fastq(fastq_path, max_errors=None):
    counts = Counter()
    n_reads = 0

    with open(fastq_path) as f:
        while True:
            header = f.readline()
            if not header:
                break

            seq = f.readline().rstrip()
            plus = f.readline()
            qual = f.readline().rstrip()

            n_reads += 1

            q = np.array([ord(c) - 33 for c in qual])
            p_err = 10 ** (-q / 10)

            probs = np.array([1.0])

            for p in p_err:
                probs = np.convolve(probs, [1 - p, p])

            if max_errors is not None:
                probs = probs[: max_errors + 1]

            for k, prob in enumerate(probs):
                counts[k] += prob

    dist = {k: v / n_reads for k, v in sorted(counts.items())}
    return dist


def dist_to_dataframe(dist):
    return pd.DataFrame({
        "n_errors": list(dist.keys()),
        "fraction_reads": list(dist.values()),
        "percent_reads": [100 * v for v in dist.values()],
    })


def save_error_distribution(dist, output_csv):
    df = dist_to_dataframe(dist)
    df.to_csv(output_csv, index=False)
    return df


def plot_error_distribution(
    dist,
    save=None,
    cumulative=False,
    logy=True,
):
    k = np.array(list(dist.keys()))
    p = np.array(list(dist.values()))

    if cumulative:
        p = 1 - np.cumsum(p) + p
        ylabel = "Reads with ≥ k Errors (%)"
        title = "Cumulative Read Error Distribution"
    else:
        ylabel = "Reads (%)"
        title = "Read Error Distribution"

    p *= 100

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(k, p)

    ax.set_xlabel("Number of Sequencing Errors")
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if logy:
        ax.set_yscale("log")

    plt.tight_layout()

    if save:
        plt.savefig(save, bbox_inches="tight", dpi=300)

    return fig, ax


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
        description=(
            "Predict the percentage of reads with 0, 1, 2, ... errors "
            "from Phred+33 FASTQ quality scores."
        ),
        formatter_class=MyFormatter
    )

    parser.add_argument(
        "fastq",
        help="Input FASTQ file with Phred+33 quality scores.",
    )

    parser.add_argument(
        "-o",
        "--output_csv",
        default="error_distribution.csv",
        help="Output CSV file. Default: error_distribution.csv",
    )

    parser.add_argument(
        "-p",
        "--plot",
        default="error_distribution.png",
        help="Output plot file. Default: error_distribution.png",
    )

    parser.add_argument(
        "--max_errors",
        type=int,
        default=None,
        help="Maximum number of errors to report.",
    )

    parser.add_argument(
        "--cumulative",
        action="store_true",
        help="Plot percentage of reads with >= k errors.",
    )

    parser.add_argument(
        "--no_logy",
        action="store_true",
        help="Disable log-scaled y-axis.",
    )

    parser.add_argument(
        "--show",
        action="store_true",
        help="Show plot interactively.",
    )

    args = parser.parse_args()

    dist = error_distribution_fastq(
        args.fastq,
        max_errors=args.max_errors,
    )

    df = save_error_distribution(
        dist,
        args.output_csv,
    )

    fig, ax = plot_error_distribution(
        dist,
        save=args.plot,
        cumulative=args.cumulative,
        logy=not args.no_logy,
    )

    if args.show:
        plt.show()
    else:
        plt.close(fig)

    print(f"Saved distribution: {args.output_csv}")
    print(f"Saved plot: {args.plot}")

    mean_errors = sum(k * v for k, v in dist.items())
    print(f"Mean expected errors/read: {mean_errors:.4f}")


if __name__ == "__main__":
    main()