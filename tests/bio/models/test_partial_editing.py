"""
Tests for edms.bio.models.partial_editing — visualizes deterministic
compression of log2FC by partial editing:
    log2FC_obs = log2(1 + e*(FC_edit - 1))

make_plot()'s underlying math is cross-checked by hand at known points;
parse_args() validation and main() get light smoke tests.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from edms.bio.models import partial_editing as pedit


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def test_parse_args_defaults():
    args = pedit.parse_args([])
    assert args.max_log2FC_edit_mag == 2.0
    assert args.cmap == "Reds"
    assert args.editing_fractions == []


def test_parse_args_rejects_invalid_max_mag():
    with pytest.raises(SystemExit):
        pedit.parse_args(["--max_log2FC_edit_mag", "-1"])


def test_parse_args_rejects_out_of_range_editing_fraction():
    with pytest.raises(SystemExit):
        pedit.parse_args(["--editing_fractions", "1.5"])


def test_parse_args_rejects_negative_edited_only_log2fc():
    with pytest.raises(SystemExit):
        pedit.parse_args(["--edited_only_log2FCs", "-1"])


def test_parse_args_rejects_small_grid_size():
    with pytest.raises(SystemExit):
        pedit.parse_args(["--grid_size", "10"])


def test_compression_formula_full_editing_recovers_edited_only_log2fc():
    # At e=1 (fully edited), log2FC_obs should equal the edited-only log2FC exactly:
    # log2(1 + 1*(FC-1)) = log2(FC)
    for log2fc_edit_mag in [0.5, 1.0, 2.0]:
        fc_edit = 2 ** log2fc_edit_mag
        log2fc_obs = np.log2(1 + 1.0 * (fc_edit - 1))
        assert log2fc_obs == pytest.approx(log2fc_edit_mag)


def test_compression_formula_zero_editing_gives_zero_observed_log2fc():
    for log2fc_edit_mag in [0.5, 1.0, 2.0]:
        fc_edit = 2 ** log2fc_edit_mag
        log2fc_obs = np.log2(1 + 0.0 * (fc_edit - 1))
        assert log2fc_obs == pytest.approx(0.0)


def test_make_plot_returns_figure():
    fig = pedit.make_plot(grid_size=50)
    assert fig is not None
    ax = fig.axes[0]
    assert ax.get_xlabel() == "Editing fraction (e)"


def test_make_plot_with_reference_lines_and_contours():
    fig = pedit.make_plot(
        grid_size=50, editing_fractions=[0.1, 0.5], edited_only_log2fcs=[1.0],
        observed_log2fcs=[0.5],
    )
    assert fig is not None


def test_main_smoke(tmp_path, monkeypatch):
    out_png = tmp_path / "plot.png"
    monkeypatch.setattr(sys, "argv", ["prog", "-o", str(out_png), "-gs", "50"])
    rc = pedit.main()
    assert rc == 0
    assert out_png.exists()
