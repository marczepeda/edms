"""
Tests for edms.bio.models.mutation_rate — mutation frequency distributions
(truncated-normal models for total/single-mutation edit rates).

truncnorm_pdf() is a thin wrapper around scipy.stats.truncnorm; it is
cross-checked directly against an independently constructed scipy call.
main() gets a light smoke test with plotting disabled.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest
from scipy.stats import truncnorm

from edms.bio.models import mutation_rate as mr


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def test_truncnorm_pdf_matches_scipy_directly():
    mu, sigma = 0.2, 0.05
    x = np.array([0.05, 0.1, 0.2, 0.3, 0.4])
    a = (0.0 - mu) / sigma
    b = (1.0 - mu) / sigma
    expected = truncnorm(a, b, loc=mu, scale=sigma).pdf(x)
    got = mr.truncnorm_pdf(x, mu, sigma)
    assert np.allclose(got, expected)


def test_truncnorm_pdf_integrates_to_one():
    mu, sigma = 0.5, 0.1
    x = np.linspace(0, 1, 20001)
    pdf = mr.truncnorm_pdf(x, mu, sigma)
    trapezoid = getattr(np, "trapezoid", None) or np.trapz
    integral = trapezoid(pdf, x)
    assert integral == pytest.approx(1.0, abs=1e-3)


def test_truncnorm_pdf_respects_custom_bounds():
    mu, sigma, lo, hi = 5.0, 2.0, 0.0, 10.0
    x = np.array([1.0, 5.0, 9.0])
    a = (lo - mu) / sigma
    b = (hi - mu) / sigma
    expected = truncnorm(a, b, loc=mu, scale=sigma).pdf(x)
    got = mr.truncnorm_pdf(x, mu, sigma, lo=lo, hi=hi)
    assert np.allclose(got, expected)


def test_main_smoke_no_plot(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "prog", "--no-plot", "--n-mut", "100", "--mean", "0.1",
    ])
    mr.main()  # should run to completion without raising and without plotting


def test_main_smoke_with_plot_and_amplicon_len(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", [
        "prog", "--n-mut", "100", "--mean", "0.1", "--amplicon-len", "150",
    ])
    mr.main()
