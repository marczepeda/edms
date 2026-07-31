"""
Tests for edms.gen.plot

plot.py is huge (~4k lines); this suite focuses on the small pure-logic
helpers (hand-computed expectations) plus happy-path smoke tests for the
public graph functions (scat, cat, dist, heat, stack, vol) - asserting they
run without error on valid input, return the documented (fig, axes) types,
and write files to disk when a `file`/`dir` is given.
"""
import math

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
import seaborn as sns

from edms.gen import plot as p


# --------------------------------------------------------------------------- #
# Small pure-logic helpers
# --------------------------------------------------------------------------- #

def test_re_un_cap_basic():
    assert p.re_un_cap("hello_world") == "Hello World"


def test_re_un_cap_no_underscore():
    assert p.re_un_cap("alreadycapitalized") == "Alreadycapitalized"


def test_re_un_cap_multiple_underscores():
    assert p.re_un_cap("a_b_c") == "A B C"


def test_round_up_pow_10():
    # round_up_pow_10() rounds up to the *next power of ten* (10**ceil(log10(n))),
    # not to a "nice" round number near n - e.g. 45 and 99 both round up all
    # the way to 100, not to 50/100 respectively.
    assert p.round_up_pow_10(45) == 100
    assert p.round_up_pow_10(150) == 1000
    assert p.round_up_pow_10(0) == 0
    assert p.round_up_pow_10(100) == 100
    assert p.round_up_pow_10(999) == 1000


def test_round_down_pow_10():
    assert p.round_down_pow_10(45) == 40
    assert p.round_down_pow_10(150) == 100
    assert p.round_down_pow_10(0) == 0
    assert p.round_down_pow_10(999) == 900


def test_log10_elementwise_clipped_at_one():
    # log10() clips values below 1 up to 1 before taking log10, so results
    # are always >= 0 (it is an elementwise transform, not "log10 of the
    # series maximum" despite the one-line docstring).
    out = p.log10([0.001, 1, 10, 100])
    np.testing.assert_allclose(out, [0.0, 0.0, 1.0, 2.0])


def test_extract_pivots():
    df = pd.DataFrame({
        "x": ["c1", "c1", "c2", "c2"],
        "y": ["r1", "r2", "r1", "r2"],
        "variable": ["v1", "v1", "v1", "v1"],
        "value": [1.0, 2.0, 3.0, 4.0],
    })
    pivots = p.extract_pivots(df, x="x", y="y")
    assert set(pivots.keys()) == {"v1"}
    piv = pivots["v1"]
    assert piv.loc["r1", "c1"] == 1.0
    assert piv.loc["r2", "c2"] == 4.0


def test_autoscale_limits_basic():
    lo, hi = p.autoscale_limits([0.0, 10.0], buffer=0.1)
    assert lo == pytest.approx(-1.0)
    assert hi == pytest.approx(11.0)


def test_autoscale_limits_all_identical_uses_min_range():
    lo, hi = p.autoscale_limits([5.0, 5.0, 5.0], min_range=2.0)
    assert lo == pytest.approx(4.0)
    assert hi == pytest.approx(6.0)


def test_autoscale_limits_empty_returns_min_range_centered_on_zero():
    lo, hi = p.autoscale_limits([], min_range=4.0)
    assert lo == pytest.approx(-2.0)
    assert hi == pytest.approx(2.0)


def test_repeat_palette_cmap_matplotlib_cmap():
    cmap = p.repeat_palette_cmap("viridis", repeats=2)
    assert isinstance(cmap, matplotlib.colors.ListedColormap)


def test_repeat_palette_cmap_invalid_repeats_raises():
    with pytest.raises(ValueError):
        p.repeat_palette_cmap("viridis", repeats=0)
    with pytest.raises(ValueError):
        p.repeat_palette_cmap("viridis", repeats=-1)


def test_repeat_palette_cmap_seaborn_palette_name():
    # 'deep' is a documented, valid seaborn palette name (see
    # seaborn_palettes() below); it should return a repeated ListedColormap
    # built from sns.palettes.SEABORN_PALETTES['deep'], not crash.
    cmap = p.repeat_palette_cmap("deep", repeats=2)
    assert isinstance(cmap, matplotlib.colors.ListedColormap)
    base_len = len(sns.palettes.SEABORN_PALETTES["deep"])
    assert cmap.N == base_len * 2


def test_repeat_palette_cmap_invalid_name_returns_input_unchanged(capsys):
    # Not a seaborn palette name and not a matplotlib colormap -> prints a
    # message and returns the input unchanged (no crash referencing an
    # unassigned variable).
    out = p.repeat_palette_cmap("not_a_real_cmap", repeats=2)
    assert out == "not_a_real_cmap"
    captured = capsys.readouterr()
    assert "not a valid matplotlib color map" in captured.out


def test_matplotlib_cmaps_smoke():
    # Just confirm it runs to completion without error (Agg backend, so
    # plt.show() is a no-op).
    p.matplotlib_cmaps()
    plt.close("all")


def test_seaborn_palettes_smoke():
    p.seaborn_palettes()
    plt.close("all")


# --------------------------------------------------------------------------- #
# save_fig()
# --------------------------------------------------------------------------- #

def test_save_fig_png(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    p.save_fig(file="plot.png", dir=str(tmp_path), fig=fig, dpi=72)
    assert (tmp_path / "plot.png").exists()
    plt.close(fig)


def test_save_fig_all_formats(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    p.save_fig(file="plot.all", dir=str(tmp_path), fig=fig, dpi=72)
    for ext in ("png", "pdf", "svg"):
        assert (tmp_path / f"plot.{ext}").exists()
    plt.close(fig)


def test_save_fig_html(tmp_path):
    fig, ax = plt.subplots()
    ax.scatter([0, 1, 2], [0, 1, 2])
    p.save_fig(file="plot.html", dir=str(tmp_path), fig=fig)
    assert (tmp_path / "plot.html").exists()
    plt.close(fig)


def test_save_fig_none_file_is_noop(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    # file=None -> check_outpath returns (None, None); should return quietly.
    p.save_fig(file=None, dir=str(tmp_path), fig=fig)
    assert list(tmp_path.iterdir()) == []
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Graph methods: happy-path smoke tests
# --------------------------------------------------------------------------- #

@pytest.fixture
def scat_df():
    rng = np.random.default_rng(0)
    return pd.DataFrame({
        "x": rng.normal(size=20),
        "y": rng.normal(size=20),
        "grp": ["a", "b"] * 10,
    })


def test_scat_basic_returns_fig_axes(scat_df):
    fig, axes = p.scat("scat", scat_df, x="x", y="y", show=False)
    assert isinstance(fig, plt.Figure)
    assert axes.shape == (1, 1)
    plt.close(fig)


def test_scat_with_color_column_and_save(tmp_path, scat_df):
    fig, axes = p.scat("scat", scat_df, x="x", y="y", cols="grp",
                        file="scat.png", dir=str(tmp_path), show=False)
    assert (tmp_path / "scat.png").exists()
    plt.close(fig)


def test_scat_line_graph_type(scat_df):
    fig, axes = p.scat("line", scat_df.sort_values("x"), x="x", y="y", show=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_scat_faceted(scat_df):
    df = scat_df.copy()
    df["facet"] = ["f1", "f2"] * 10
    fig, axes = p.scat("scat", df, x="x", y="y", facetx="facet", show=False)
    assert axes.shape == (1, 2)
    plt.close(fig)


@pytest.fixture
def cat_df():
    return pd.DataFrame({
        "grp": ["a", "a", "a", "b", "b", "b"],
        "val": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    })


def test_cat_bar(cat_df):
    fig, axes = p.cat("bar", cat_df, x="grp", y="val", show=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_cat_box_and_save(tmp_path, cat_df):
    fig, axes = p.cat("box", cat_df, x="grp", y="val", file="cat.png", dir=str(tmp_path), show=False)
    assert (tmp_path / "cat.png").exists()
    plt.close(fig)


def test_cat_count_x_only(cat_df):
    # x-only (y left at its default '') should now autoscale correctly
    # instead of trying to look up df[''].
    fig, axes = p.cat("count", cat_df, x="grp", show=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_cat_count_y_only(cat_df):
    fig, axes = p.cat("count", cat_df, y="grp", show=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_cat_count_with_both_x_and_y_raises(cat_df):
    with pytest.raises(ValueError):
        p.cat("count", cat_df, x="grp", y="val", show=False)


@pytest.fixture
def dist_df():
    rng = np.random.default_rng(1)
    return pd.DataFrame({
        "val": rng.normal(size=50),
        "grp": ["a", "b"] * 25,
    })


def test_dist_hist(dist_df):
    fig, axes = p.dist("hist", dist_df, x="val", show=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_dist_kde_with_cols_and_save(tmp_path, dist_df):
    fig, axes = p.dist("kde", dist_df, x="val", cols="grp", file="dist.png", dir=str(tmp_path), show=False)
    assert (tmp_path / "dist.png").exists()
    plt.close(fig)


def test_heat_matrix_input():
    # Single-pivot input path (x=y=vals=None): a plain pre-pivoted matrix.
    matrix = pd.DataFrame(
        [[1.0, 2.0], [3.0, 4.0]],
        index=["r1", "r2"],
        columns=["c1", "c2"],
    )
    fig, axes = p.heat(matrix, show=False)
    assert isinstance(fig, plt.Figure)
    assert axes.shape == (1, 1)
    plt.close(fig)


def test_heat_tidy_input_and_save(tmp_path):
    df = pd.DataFrame({
        "x": ["c1", "c1", "c2", "c2"],
        "y": ["r1", "r2", "r1", "r2"],
        "vals": [1.0, 2.0, 3.0, 4.0],
    })
    fig, axes = p.heat(df, x="x", y="y", vals="vals", file="heat.png", dir=str(tmp_path), show=False)
    assert (tmp_path / "heat.png").exists()
    plt.close(fig)


@pytest.fixture
def stack_df():
    return pd.DataFrame({
        "x": ["g1", "g1", "g2", "g2"],
        "y": [10.0, 20.0, 15.0, 25.0],
        "cols": ["cat1", "cat2", "cat1", "cat2"],
    })


def test_stack_basic(stack_df):
    fig, axes = p.stack(stack_df, x="x", y="y", cols="cols", show=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_stack_save(tmp_path, stack_df):
    fig, axes = p.stack(stack_df, x="x", y="y", cols="cols", file="stack.png", dir=str(tmp_path), show=False)
    assert (tmp_path / "stack.png").exists()
    plt.close(fig)


@pytest.fixture
def vol_df():
    rng = np.random.default_rng(2)
    return pd.DataFrame({
        "log2fc": rng.normal(size=30),
        "neglog10p": np.abs(rng.normal(size=30)),
    })


def test_vol_return_df_true_returns_dataframe(vol_df):
    out = p.vol(vol_df, x="log2fc", y="neglog10p", show=False, return_df=True)
    assert isinstance(out, pd.DataFrame)
    assert "log2fc" in out.columns


def test_vol_return_df_false_returns_fig_axes(vol_df):
    fig, axes = p.vol(vol_df, x="log2fc", y="neglog10p", show=False, return_df=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_vol_save(tmp_path, vol_df):
    out = p.vol(vol_df, x="log2fc", y="neglog10p", file="vol.png", dir=str(tmp_path), show=False, return_df=False)
    assert (tmp_path / "vol.png").exists()
