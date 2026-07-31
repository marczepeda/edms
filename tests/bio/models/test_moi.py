"""
Tests for edms.bio.models.moi — MOI (multiplicity of infection) Poisson
infection-outcome curves.

moi_curves() is cross-checked against scipy.stats.poisson independently;
parse_float_list is tested with hand-picked inputs; plot_moi()/main() get
smoke tests.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest
from scipy.stats import poisson

from edms.bio.models import moi


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def test_moi_curves_matches_scipy_poisson():
    moi_vals = np.array([0.0, 0.5, 1.0, 2.0, 5.0])
    total, once, multiple = moi.moi_curves(moi_vals)

    assert np.allclose(total, 1 - poisson.pmf(0, moi_vals))
    assert np.allclose(once, poisson.pmf(1, moi_vals))
    assert np.allclose(multiple, total - once)
    # sanity: at moi=0 nobody is infected
    assert total[0] == pytest.approx(0.0)
    assert once[0] == pytest.approx(0.0)


def test_moi_curves_multiple_equals_poisson_ge2():
    moi_vals = np.array([0.3, 1.5, 3.0])
    total, once, multiple = moi.moi_curves(moi_vals)
    expected_multiple = 1 - poisson.pmf(0, moi_vals) - poisson.pmf(1, moi_vals)
    assert np.allclose(multiple, expected_multiple)


def test_parse_float_list_comma_and_space_separated():
    assert moi.parse_float_list("0.1,0.3, 1.0") == [0.1, 0.3, 1.0]
    assert moi.parse_float_list("0.1 0.3 1.0") == [0.1, 0.3, 1.0]
    assert moi.parse_float_list("0.1;0.3;1.0") == [0.1, 0.3, 1.0]


def test_parse_float_list_none_returns_empty():
    assert moi.parse_float_list(None) == []


def test_parse_float_list_invalid_raises():
    import argparse
    with pytest.raises(argparse.ArgumentTypeError):
        moi.parse_float_list("abc")


def test_ensure_dir_creates_parent(tmp_path):
    target = tmp_path / "nested" / "dir" / "file.png"
    moi.ensure_dir(target)
    assert target.parent.is_dir()


def test_plot_moi_raises_on_invalid_max_moi():
    with pytest.raises(ValueError):
        moi.plot_moi(max_moi=0)


def test_plot_moi_raises_on_invalid_n_points():
    with pytest.raises(ValueError):
        moi.plot_moi(n_points=1)


def test_plot_moi_saves_figure(tmp_path):
    out = tmp_path / "moi.png"
    moi.plot_moi(max_moi=3.0, n_points=50, out_files=[out], show=False)
    assert out.exists()


def test_build_parser_defaults():
    parser = moi.build_parser()
    args = parser.parse_args([])
    assert args.max_moi == 5.0
    assert args.percent is True


def test_main_smoke(tmp_path, monkeypatch):
    out_png = tmp_path / "moi.png"
    monkeypatch.setattr(sys, "argv", ["prog", "-o", str(out_png), "--max-moi", "3"])
    moi.main()
    assert out_png.exists()
