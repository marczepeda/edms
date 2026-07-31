"""
Tests for edms.bio.models.possion_bias — bias-factor curves for clustered
editing (overdispersed Poisson model for multiple edits per cell).

Pure-math helpers are cross-checked against scipy.stats.poisson and hand
derivations; plot_p_two_diff_vs_p_any()/main() get light smoke tests.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest
from scipy.stats import poisson

from edms.bio.models import possion_bias as pb


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def test_poisson_pmf_matches_scipy():
    lam = 2.5
    kmax = 10
    got = pb.poisson_pmf(lam, kmax)
    expected = poisson.pmf(np.arange(kmax + 1), lam)
    assert np.allclose(got, expected)


def test_poisson_pmf_sums_close_to_one_for_large_kmax():
    got = pb.poisson_pmf(3.0, kmax=30)
    assert got.sum() == pytest.approx(1.0, abs=1e-6)


def test_biased_pmf_b1_reduces_to_plain_poisson():
    # At b=1 (no bias), lam is chosen so that base[0] == 1 - p_any exactly
    # (lam = -ln(1-p_any) => exp(-lam) = 1-p_any), and the renormalization of
    # k>=1 mass is a no-op since sum(base[1:]) == p_any exactly by construction.
    p_any = 0.3
    kmax = 20
    lam = -np.log(1 - p_any)
    base = poisson.pmf(np.arange(kmax + 1), lam)

    biased = pb.biased_pmf_from_p_any(p_any, b=1.0, kmax=kmax)
    assert biased[0] == pytest.approx(1 - p_any)
    assert biased[0] == pytest.approx(base[0])
    assert np.allclose(biased[1:], base[1:])


def test_biased_pmf_sums_to_one():
    biased = pb.biased_pmf_from_p_any(0.4, b=2.5, kmax=20)
    assert biased.sum() == pytest.approx(1.0)


def test_biased_pmf_higher_bias_shifts_mass_to_higher_k():
    # With more bias (b>1), P(k>=2 | at least one edit) should be relatively
    # larger than under the unbiased (b=1) baseline, since higher k is
    # up-weighted by b^(k-1).
    p_any = 0.3
    kmax = 20
    unbiased = pb.biased_pmf_from_p_any(p_any, b=1.0, kmax=kmax)
    biased = pb.biased_pmf_from_p_any(p_any, b=3.0, kmax=kmax)
    frac_ge2_unbiased = unbiased[2:].sum() / unbiased[1:].sum()
    frac_ge2_biased = biased[2:].sum() / biased[1:].sum()
    assert frac_ge2_biased > frac_ge2_unbiased


def test_biased_pmf_invalid_p_any_raises():
    with pytest.raises(ValueError):
        pb.biased_pmf_from_p_any(0.0, b=1.0)
    with pytest.raises(ValueError):
        pb.biased_pmf_from_p_any(1.0, b=1.0)


def test_biased_pmf_invalid_b_raises():
    with pytest.raises(ValueError):
        pb.biased_pmf_from_p_any(0.3, b=0.5)


def test_p_two_diff_is_half_of_pmf_at_k2():
    p_any, b, kmax = 0.3, 2.0, 20
    pmf = pb.biased_pmf_from_p_any(p_any, b, kmax=kmax)
    assert pb.p_two_diff(p_any, b, kmax=kmax) == pytest.approx(0.5 * pmf[2])


def test_plot_p_two_diff_vs_p_any_writes_outputs(tmp_path):
    out_png = tmp_path / "p.png"
    out_csv = tmp_path / "p.csv"
    png, csv = pb.plot_p_two_diff_vs_p_any(
        p_any_grid=np.linspace(0.05, 0.3, 5), b_list=[1.0, 2.0],
        out_png=out_png, out_csv=out_csv,
    )
    assert png.exists()
    assert csv.exists()


def test_main_smoke(tmp_path, monkeypatch):
    out_png = tmp_path / "p.png"
    out_csv = tmp_path / "p.csv"
    monkeypatch.setattr(sys, "argv", [
        "prog", "--points", "5", "--out-png", str(out_png), "--out-csv", str(out_csv),
    ])
    pb.main()
    assert out_png.exists()
    assert out_csv.exists()
