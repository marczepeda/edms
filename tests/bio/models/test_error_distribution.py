"""
Tests for edms.bio.models.error_distribution — predicts the percentage of
reads with 0, 1, 2, ... sequencing errors from Phred+33 FASTQ quality scores.

error_distribution_fastq() is cross-checked against an independent
hand-computed convolution of per-base error probabilities; main() gets a
light smoke test over a tiny synthetic FASTQ file.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from edms.bio.models import error_distribution as ed


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def _manual_error_count_dist(quals):
    """Independent reimplementation: probability of exactly k errors per read via convolution."""
    n_reads = len(quals)
    from collections import Counter
    counts = Counter()
    for qual in quals:
        q = np.array([ord(c) - 33 for c in qual])
        p_err = 10 ** (-q / 10)
        probs = np.array([1.0])
        for p in p_err:
            probs = np.convolve(probs, [1 - p, p])
        for k, prob in enumerate(probs):
            counts[k] += prob
    return {k: v / n_reads for k, v in sorted(counts.items())}


@pytest.fixture()
def mini_fastq(tmp_path):
    # read1: quals '!'=Q0 (p_err=1.0), 'I'=Q40 (p_err=1e-4)
    # read2: quals 'I','I' both Q40
    content = (
        "@read1\nAC\n+\n!I\n"
        "@read2\nGT\n+\nII\n"
    )
    pt = tmp_path / "test.fastq"
    pt.write_text(content)
    return str(pt), ["!I", "II"]


def test_error_distribution_fastq_matches_manual_convolution(mini_fastq):
    fastq_path, quals = mini_fastq
    dist = ed.error_distribution_fastq(fastq_path)
    expected = _manual_error_count_dist(quals)
    assert set(dist.keys()) == set(expected.keys())
    for k in expected:
        assert dist[k] == pytest.approx(expected[k])


def test_error_distribution_fastq_max_errors_truncates(mini_fastq):
    fastq_path, _ = mini_fastq
    dist = ed.error_distribution_fastq(fastq_path, max_errors=0)
    assert set(dist.keys()) == {0}


def test_dist_to_dataframe_shape_and_values():
    dist = {0: 0.9, 1: 0.1}
    df = ed.dist_to_dataframe(dist)
    assert list(df["n_errors"]) == [0, 1]
    assert list(df["fraction_reads"]) == [0.9, 0.1]
    assert list(df["percent_reads"]) == pytest.approx([90.0, 10.0])


def test_save_error_distribution_writes_csv(tmp_path):
    dist = {0: 0.9, 1: 0.1}
    out_csv = tmp_path / "dist.csv"
    df = ed.save_error_distribution(dist, str(out_csv))
    assert out_csv.exists()
    reloaded = pd.read_csv(out_csv)
    assert list(reloaded["n_errors"]) == [0, 1]


def test_plot_error_distribution_returns_fig_ax():
    dist = {0: 0.9, 1: 0.09, 2: 0.01}
    fig, ax = ed.plot_error_distribution(dist, cumulative=False, logy=True)
    assert fig is not None
    assert ax.get_xlabel() == "Number of Sequencing Errors"


def test_plot_error_distribution_cumulative_mode():
    dist = {0: 0.9, 1: 0.09, 2: 0.01}
    fig, ax = ed.plot_error_distribution(dist, cumulative=True, logy=False)
    assert "Cumulative" in ax.get_title()


def test_main_smoke(mini_fastq, tmp_path, monkeypatch):
    fastq_path, _ = mini_fastq
    out_csv = tmp_path / "dist.csv"
    out_png = tmp_path / "dist.png"
    monkeypatch.setattr(sys, "argv", [
        "prog", fastq_path, "-o", str(out_csv), "-p", str(out_png),
    ])
    ed.main()
    assert out_csv.exists()
    assert out_png.exists()
