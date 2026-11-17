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
from matplotlib.collections import PolyCollection
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


def _visible_ticks(ax, axis="x", minor=False):
    """Tick locations inside the axis limits, which is what the reader actually sees."""
    lo, hi = ax.get_xlim() if axis == "x" else ax.get_ylim()
    locs = ax.get_xticks(minor=minor) if axis == "x" else ax.get_yticks(minor=minor)
    return [loc for loc in locs if lo <= loc <= hi]


def _visible_tick_labels(ax, axis="x"):
    """Major tick label text for the ticks inside the axis limits.

    A locator may hand back a tick just past each end; those are never drawn, but
    get_xticklabels() still lists them, so they are filtered out by position here.
    """
    lo, hi = ax.get_xlim() if axis == "x" else ax.get_ylim()
    labels = ax.get_xticklabels() if axis == "x" else ax.get_yticklabels()
    pos = 0 if axis == "x" else 1
    return [t.get_text() for t in labels if lo <= t.get_position()[pos] <= hi]


def test_log_ticks_majors_on_powers_of_ten():
    fig, ax = plt.subplots()
    ax.set_xscale("log")
    ax.set_xlim(1e2, 1e7)
    p.log_ticks(ax, "x")
    assert _visible_ticks(ax, "x") == pytest.approx([1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
    plt.close(fig)


def test_log_ticks_minors_span_the_decade():
    fig, ax = plt.subplots()
    ax.set_xscale("log")
    ax.set_xlim(1e3, 1e4)
    p.log_ticks(ax, "x")
    # 0.1x10^n coincides with a major tick, so matplotlib leaves 0.2 through 0.9 showing
    assert _visible_ticks(ax, "x", minor=True) == pytest.approx(
        [sub * 1e4 for sub in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)]
    )
    plt.close(fig)


def test_log_ticks_minors_are_unlabeled():
    fig, ax = plt.subplots()
    ax.set_yscale("log")
    ax.set_ylim(1e0, 1e3)
    p.log_ticks(ax, "y")
    fig.canvas.draw()
    assert all(t.get_text() == "" for t in ax.get_yticklabels(minor=True))
    labels = _visible_tick_labels(ax, "y")
    assert all("10^{" in label for label in labels)
    assert any("10^{1}" in label for label in labels)
    plt.close(fig)


def test_log_ticks_both_axes():
    fig, ax = plt.subplots()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e1, 1e3)
    ax.set_ylim(1e1, 1e3)
    p.log_ticks(ax)
    assert _visible_ticks(ax, "x") == pytest.approx([1e1, 1e2, 1e3])
    assert _visible_ticks(ax, "y") == pytest.approx([1e1, 1e2, 1e3])
    plt.close(fig)


def test_log_ticks_leaves_a_linear_axis_alone():
    fig, ax = plt.subplots()
    ax.set_xlim(0, 10)
    before = list(ax.get_xticks())
    p.log_ticks(ax, "x")
    assert list(ax.get_xticks()) == before
    assert _visible_ticks(ax, "x", minor=True) == []
    plt.close(fig)


def test_log_ticks_custom_subs():
    fig, ax = plt.subplots()
    ax.set_xscale("log")
    ax.set_xlim(1e2, 1e3)
    p.log_ticks(ax, "x", subs=(0.5,))
    assert _visible_ticks(ax, "x", minor=True) == pytest.approx([5e2])
    plt.close(fig)


def test_log_ticks_invalid_axis_raises():
    fig, ax = plt.subplots()
    with pytest.raises(ValueError):
        p.log_ticks(ax, "z")
    plt.close(fig)


def test_scat_log_scale_gets_decade_ticks():
    df = pd.DataFrame({"x": [1e2, 1e4, 1e6], "y": [1e1, 1e2, 1e3]})
    fig, axes = p.scat("scat", df, x="x", y="y",
                       x_axis_scale="log", y_axis_scale="log",
                       x_axis_dims=(1e2, 1e6), y_axis_dims=(1e1, 1e3), show=False)
    ax = axes[0, 0]
    assert _visible_ticks(ax, "x") == pytest.approx([1e2, 1e3, 1e4, 1e5, 1e6])
    assert len(_visible_ticks(ax, "x", minor=True)) == 8 * 4  # 0.2 to 0.9 per decade
    assert _visible_ticks(ax, "y") == pytest.approx([1e1, 1e2, 1e3])
    plt.close(fig)


def test_log10_ticks_majors_on_whole_decades():
    fig, ax = plt.subplots()
    ax.set_xlim(3.0, 6.0)  # linear axis holding log10 values
    p.log10_ticks(ax, "x")
    fig.canvas.draw()
    assert _visible_ticks(ax, "x") == pytest.approx([3.0, 4.0, 5.0, 6.0])
    assert _visible_tick_labels(ax, "x") == [
        f"$\\mathdefault{{10^{{{n}}}}}$" for n in (3, 4, 5, 6)
    ]
    plt.close(fig)


def test_log10_ticks_minors_sit_at_log10_of_the_subs():
    fig, ax = plt.subplots()
    ax.set_xlim(4.0, 5.0)
    p.log10_ticks(ax, "x")
    # 0.1x10^5 is 10^4, which coincides with a major, leaving 0.2 through 0.9 of that decade
    assert _visible_ticks(ax, "x", minor=True) == pytest.approx(
        [5 + math.log10(sub) for sub in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)]
    )
    plt.close(fig)


def test_log10_ticks_handles_negative_decades():
    fig, ax = plt.subplots()
    ax.set_xlim(-2.0, 0.0)
    p.log10_ticks(ax, "x")
    fig.canvas.draw()
    assert _visible_ticks(ax, "x") == pytest.approx([-2.0, -1.0, 0.0])
    assert _visible_tick_labels(ax, "x")[0] == "$\\mathdefault{10^{-2}}$"
    plt.close(fig)


def test_log10_ticks_leaves_a_log_axis_alone():
    fig, ax = plt.subplots()
    ax.set_xscale("log")
    ax.set_xlim(1e2, 1e4)
    p.log10_ticks(ax, "x")
    assert _visible_ticks(ax, "x") == pytest.approx([1e2, 1e3, 1e4])  # not 2, 3, 4
    plt.close(fig)


def test_log10_ticks_invalid_axis_raises():
    fig, ax = plt.subplots()
    with pytest.raises(ValueError):
        p.log10_ticks(ax, "z")
    plt.close(fig)


@pytest.mark.parametrize("graph", ["kde", "hist_kde", "offset", "half_offset"])
def test_dist_log_scale_log10_graphs_get_decade_ticks(graph):
    rng = np.random.default_rng(3)
    df = pd.DataFrame({"v": 10 ** rng.normal(4, 0.5, 300), "g": ["a", "b", "c"] * 100})
    fig, axes = p.dist(graph, df, x="v", cols="g", x_axis_scale="log", show=False)
    ax = axes[0, 0]
    assert ax.get_xscale() == "linear"  # the values were log10 transformed for the kde
    assert all(float(t).is_integer() for t in _visible_ticks(ax, "x"))
    assert all("10^{" in label for label in _visible_tick_labels(ax, "x"))
    assert len(_visible_ticks(ax, "x", minor=True)) > 0
    assert "log10" not in ax.get_xlabel().lower()  # the ticks already say the decade
    plt.close(fig)


def test_dist_hist_log_scale_uses_a_real_log_axis():
    rng = np.random.default_rng(4)
    df = pd.DataFrame({"v": 10 ** rng.normal(4, 0.5, 300)})
    fig, axes = p.dist("hist", df, x="v", x_axis_scale="log", show=False)
    ax = axes[0, 0]
    assert ax.get_xscale() == "log"
    assert all("10^{" in label for label in _visible_tick_labels(ax, "x"))
    plt.close(fig)


def test_dist_linear_scale_is_untouched():
    rng = np.random.default_rng(5)
    df = pd.DataFrame({"v": rng.normal(0, 1, 200)})
    fig, axes = p.dist("half_offset", df, x="v", show=False)
    ax = axes[0, 0]
    assert _visible_ticks(ax, "x", minor=True) == []
    assert ax.get_xlabel() == "V"
    plt.close(fig)


def test_scat_linear_scale_keeps_default_ticks():
    df = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0]})
    fig, axes = p.scat("scat", df, x="x", y="y", show=False)
    assert _visible_ticks(axes[0, 0], "x", minor=True) == []
    plt.close(fig)


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


def test_repeat_palette_cmap_color_list():
    # A plain list of colors is repeated as-is (no LinearSegmentedColormap needed).
    cmap = p.repeat_palette_cmap(["#4EC569", "#FA5757"], repeats=3)
    assert isinstance(cmap, matplotlib.colors.ListedColormap)
    assert cmap.N == 6


def test_is_color_list():
    assert p._is_color_list(["#4EC569", "red", (0, 0, 1)])
    assert not p._is_color_list("colorblind")
    assert not p._is_color_list([])
    assert not p._is_color_list(["not_a_color"])
    assert not p._is_color_list({"a": "#4EC569"})


def test_palette_colors_list_cycles_and_truncates():
    colors = ["#4EC569", "#FA5757"]
    assert len(p.palette_colors(colors)) == 2
    # cycles when more colors are requested than given
    assert p.palette_colors(colors, n_colors=5) == [matplotlib.colors.to_rgba(c)
                                                    for c in colors * 3][:5]
    # truncates when fewer are requested
    assert p.palette_colors(colors, n_colors=1) == [matplotlib.colors.to_rgba(colors[0])]


def test_palette_colors_name_and_colormap():
    assert len(p.palette_colors("colorblind", n_colors=4)) == 4
    assert len(p.palette_colors("viridis", n_colors=7)) == 7
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        "custom", ["#4EC569", "#FA5757"], N=10
    )
    assert len(p.palette_colors(cmap)) == 10
    assert len(p.palette_colors(cmap, n_colors=3)) == 3


def test_palette_colors_invalid_name_falls_back(capsys):
    colors = p.palette_colors("not_a_real_palette", n_colors=3)
    assert len(colors) == 3
    assert "used colorblind instead" in capsys.readouterr().out


def test_category_color_map_accepts_color_list():
    lut = p._category_color_map(["a", "b"], palette_or_cmap=["#4EC569", "#FA5757"])
    assert lut["a"] == matplotlib.colors.to_rgba("#4EC569")
    assert lut["b"] == matplotlib.colors.to_rgba("#FA5757")


def test_category_color_map_passes_dict_through():
    lut = {"a": "#4EC569", "b": "#FA5757"}
    assert p._category_color_map(["a", "b"], palette_or_cmap=lut) is lut


def test_as_cmap_inputs():
    # list of colors and seaborn palette names both become colormaps; matplotlib
    # colormap names are left for matplotlib/seaborn to resolve downstream.
    assert isinstance(p._as_cmap(["#4EC569", "#FA5757"]), matplotlib.colors.ListedColormap)
    assert isinstance(p._as_cmap("colorblind"), matplotlib.colors.ListedColormap)
    assert p._as_cmap("Reds") == "Reds"


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
    p.save_fig(file=str(tmp_path / "plot.png"), fig=fig, dpi=72)
    assert (tmp_path / "plot.png").exists()
    plt.close(fig)


def test_save_fig_all_formats(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    p.save_fig(file=str(tmp_path / "plot.all"), fig=fig, dpi=72)
    for ext in ("png", "pdf", "svg"):
        assert (tmp_path / f"plot.{ext}").exists()
    plt.close(fig)


def test_save_fig_html(tmp_path):
    fig, ax = plt.subplots()
    ax.scatter([0, 1, 2], [0, 1, 2])
    p.save_fig(file=str(tmp_path / "plot.html"), fig=fig)
    assert (tmp_path / "plot.html").exists()
    plt.close(fig)


def test_save_fig_none_file_is_noop(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    # file=None -> check_outpath returns None; should return quietly.
    p.save_fig(file=None, fig=fig)
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
                        file=str(tmp_path / "scat.png"), show=False)
    assert (tmp_path / "scat.png").exists()
    plt.close(fig)


def test_scat_line_graph_type(scat_df):
    fig, axes = p.scat("line", scat_df.sort_values("x"), x="x", y="y", show=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_scat_scat_line_connects_raw_points_per_color(scat_df):
    df = pd.concat([scat_df, scat_df.assign(y=scat_df["y"] + 1)])  # replicate x values would be averaged by line_scat
    fig, axes = p.scat("scat_line", df, x="x", y="y", cols="grp", show=False)
    ax = axes[0, 0]
    lines = [ln for ln in ax.get_lines() if len(ln.get_xdata()) > 0]
    assert len(lines) == 2  # one line per color
    assert sorted(len(ln.get_xdata()) for ln in lines) == [20, 20]  # every raw point, no aggregation
    assert not any(isinstance(c, PolyCollection) for c in ax.collections)  # no error band
    assert len(ax.collections) >= 1  # scatter points
    plt.close(fig)


def test_scat_scat_line_without_color_column(scat_df):
    fig, axes = p.scat("scat_line", scat_df, x="x", y="y", show=False)
    lines = [ln for ln in axes[0, 0].get_lines() if len(ln.get_xdata()) > 0]
    assert len(lines) == 1
    assert len(lines[0].get_xdata()) == len(scat_df)
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
    fig, axes = p.cat("box", cat_df, x="grp", y="val", file=str(tmp_path / "cat.png"), show=False)
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
    fig, axes = p.dist("kde", dist_df, x="val", cols="grp", file=str(tmp_path / "dist.png"), show=False)
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
    fig, axes = p.heat(df, x="x", y="y", vals="vals", file=str(tmp_path / "heat.png"), show=False)
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
    fig, axes = p.stack(stack_df, x="x", y="y", cols="cols", file=str(tmp_path / "stack.png"), show=False)
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
    out = p.vol(vol_df, x="log2fc", y="neglog10p", file=str(tmp_path / "vol.png"), show=False, return_df=False)
    assert (tmp_path / "vol.png").exists()


# --------------------------------------------------------------------------- #
# palette_or_cmap / cmap input types across the graph methods
# --------------------------------------------------------------------------- #

CUSTOM_COLORS = ["#4EC569", "#2FCFE4", "#A969D6", "#FA5757"]
CUSTOM_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list(
    "custom", CUSTOM_COLORS, N=4
)


@pytest.mark.parametrize("palette", [
    CUSTOM_COLORS,                      # list of colors
    CUSTOM_CMAP,                        # colormap object
    "colorblind",                       # seaborn palette name
    "viridis",                          # matplotlib colormap name
    {"a": "#4EC569", "b": "#FA5757"},   # {category: color}
])
def test_cat_palette_input_types_without_cols(cat_df, palette):
    fig, axes = p.cat("bar", cat_df, x="grp", y="val", palette_or_cmap=palette, show=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


@pytest.mark.parametrize("palette", [CUSTOM_COLORS, CUSTOM_CMAP, "colorblind"])
def test_cat_palette_input_types_with_cols(cat_df, palette):
    df = cat_df.copy()
    df["hue"] = ["x", "y"] * 3
    fig, axes = p.cat("bar", df, x="grp", y="val", cols="hue",
                      palette_or_cmap=palette, show=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


def test_cat_color_list_maps_onto_cats_order(cat_df):
    # Colors are assigned to categories in cats_order (seaborn desaturates bar
    # faces by 0.75, which is its own long-standing default).
    colors = ["#4EC569", "#FA5757"]
    fig, axes = p.cat("bar", cat_df, x="grp", y="val", cats_order=["a", "b"],
                      palette_or_cmap=colors, show=False)
    faces = [matplotlib.colors.to_hex(patch.get_facecolor())
             for patch in axes[0, 0].patches][:2]
    assert faces == [matplotlib.colors.to_hex(sns.desaturate(c, 0.75)) for c in colors]
    plt.close(fig)


@pytest.mark.parametrize("palette", [CUSTOM_COLORS, CUSTOM_CMAP])
def test_scat_and_dist_palette_input_types(scat_df, dist_df, palette):
    fig, axes = p.scat("scat", scat_df, x="x", y="y", cols="grp",
                       palette_or_cmap=palette, show=False)
    plt.close(fig)
    fig, axes = p.dist("hist", dist_df, x="val", cols="grp",
                       palette_or_cmap=palette, show=False)
    plt.close(fig)


@pytest.mark.parametrize("palette", [CUSTOM_COLORS, CUSTOM_CMAP])
def test_stack_palette_input_types(stack_df, palette):
    fig, axes = p.stack(stack_df, x="x", y="y", cols="cols",
                        palette_or_cmap=palette, show=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


@pytest.mark.parametrize("cmap", [CUSTOM_COLORS, CUSTOM_CMAP, "colorblind", "Reds"])
def test_heat_cmap_input_types(cmap):
    df = pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0]}, index=["r1", "r2"])
    fig, axes = p.heat(df, cmap=cmap, show=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Offset distributions & grid() composition
# --------------------------------------------------------------------------- #

def test_dist_offset_labels_baselines_with_cols(dist_df):
    # One baseline per cols value, spaced a full trace height apart, first value on top.
    fig, axes = p.dist("offset", dist_df, x="val", cols="grp", show=False)
    ax = axes[0, 0]
    assert list(ax.get_yticks()) == [1.0, 0.0]
    assert [t.get_text() for t in ax.get_yticklabels()] == ["A", "B"]
    plt.close(fig)


def test_dist_half_offset_halves_the_spacing(dist_df):
    fig, axes = p.dist("half_offset", dist_df, x="val", cols="grp", show=False)
    assert list(axes[0, 0].get_yticks()) == [0.5, 0.0]
    plt.close(fig)


def test_dist_offset_ticks_false_makes_a_legend(dist_df):
    fig, axes = p.dist("offset", dist_df, x="val", cols="grp", offset_ticks=False, show=False)
    assert list(axes[0, 0].get_yticks()) == []
    assert len(fig.legends) == 1
    plt.close(fig)


def _vlines(ax):
    """x positions of the dashed reference lines, which are the only dashed lines dist() draws."""
    return sorted(line.get_xdata()[0] for line in ax.lines if line.get_linestyle() == "--")


def test_dist_line_draws_a_vertical_line(dist_df):
    fig, axes = p.dist("hist", dist_df, x="val", line=0.5, show=False)
    assert _vlines(axes[0, 0]) == pytest.approx([0.5])
    plt.close(fig)


def test_dist_line_accepts_a_list(dist_df):
    fig, axes = p.dist("hist", dist_df, x="val", line=[-1.0, 0.0, 1.0], show=False)
    assert _vlines(axes[0, 0]) == pytest.approx([-1.0, 0.0, 1.0])
    plt.close(fig)


def test_dist_line_none_draws_nothing(dist_df):
    fig, axes = p.dist("hist", dist_df, x="val", show=False)
    assert _vlines(axes[0, 0]) == []
    plt.close(fig)


def test_dist_line_is_drawn_on_every_facet_panel(dist_df):
    fig, axes = p.dist("hist", dist_df, x="val", facetx="grp", line=0.5, show=False)
    assert all(_vlines(ax) == pytest.approx([0.5]) for ax in axes.ravel())
    plt.close(fig)


def test_dist_line_moves_into_log10_space_on_a_log_scale():
    # kde & offset graphs on a log scale are drawn as log10 values on a linear axis.
    df = pd.DataFrame({"val": [10.0, 100.0, 1000.0, 5000.0] * 5, "grp": ["a", "b"] * 10})
    fig, axes = p.dist("half_offset", df, x="val", cols="grp", x_axis_scale="log",
                       line=1000, show=False)
    assert _vlines(axes[0, 0]) == pytest.approx([3.0]) # log10(1000), where the data sits
    plt.close(fig)


def test_dist_line_stays_in_x_units_on_a_real_log_axis():
    # hist keeps a genuine log axis, so the line needs no transform.
    df = pd.DataFrame({"val": [10.0, 100.0, 1000.0, 5000.0] * 5})
    fig, axes = p.dist("hist", df, x="val", x_axis_scale="log", line=1000, show=False)
    assert axes[0, 0].get_xscale() == "log"
    assert _vlines(axes[0, 0]) == pytest.approx([1000.0])
    plt.close(fig)


def test_dist_line_skips_non_positive_values_on_a_log_scale():
    df = pd.DataFrame({"val": [10.0, 100.0, 1000.0, 5000.0] * 5, "grp": ["a", "b"] * 10})
    fig, axes = p.dist("half_offset", df, x="val", cols="grp", x_axis_scale="log",
                       line=[0, -5, 100], show=False)
    assert _vlines(axes[0, 0]) == pytest.approx([2.0]) # Only log10(100) survives
    plt.close(fig)


def test_dist_line_also_marks_the_stat_panel(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", x_axis_scale="log",
                       offset_stat="gmean", line=500, show=False)
    dist_ax, stat_ax = fig.axes
    assert _vlines(dist_ax) == pytest.approx([math.log10(500)])
    assert _vlines(stat_ax) == pytest.approx([500.0]) # The bar panel keeps a real log axis
    plt.close(fig)


def test_dist_line_does_not_stretch_the_stat_panel(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", x_axis_scale="log",
                       offset_stat="gmean", offset_stat_dims=(100, 10000), line=10 ** 9, show=False)
    assert fig.axes[1].get_xlim() == (100, 10000)
    plt.close(fig)


def test_dist_offset_single_value_category_does_not_raise():
    # gaussian_kde needs spread; a lone value falls back to a narrow spike.
    df = pd.DataFrame({"val": [1.0, 2.0, 3.0, 5.0], "grp": ["a", "a", "a", "b"]})
    fig, axes = p.dist("offset", df, x="val", cols="grp", show=False)
    assert isinstance(fig, plt.Figure)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Offset distribution summary statistic bar panels
# --------------------------------------------------------------------------- #

@pytest.fixture
def offset_stat_df():
    rng = np.random.default_rng(7)
    rows = []
    for line in ["L1", "L2"]:
        for block in ["top", "bottom"]:
            for grp, center in [("a", 100.0), ("b", 1000.0)]:
                for rep in ["r1", "r2", "r3"]:
                    for value in rng.lognormal(mean=math.log(center), sigma=0.3, size=40):
                        rows.append({"val": value, "grp": grp, "line": line, "block": block,
                                     "rep": f"{line}-{block}-{grp}-{rep}"})
    return pd.DataFrame(rows)


def _bars(ax):
    """(center, right edge) of each horizontal bar, in drawing order."""
    return [(rect.get_y() + rect.get_height() / 2, rect.get_x() + rect.get_width())
            for rect in ax.patches]


def test_offset_stats_mean_matches_a_groupby(offset_stat_df):
    stats = p._offset_stats(offset_stat_df, "val", cols="grp", cats=["a", "b"], stat="mean")
    expected = offset_stat_df.groupby("grp")["val"].mean()
    assert stats["a"]["center"] == pytest.approx(expected["a"])
    assert stats["b"]["center"] == pytest.approx(expected["b"])
    # Without replicates there is no spread to show.
    assert stats["a"]["lower"] is None and stats["a"]["upper"] is None


def test_offset_stats_gmean_is_the_exponentiated_log_mean(offset_stat_df):
    stats = p._offset_stats(offset_stat_df, "val", cols="grp", cats=["a"], stat="gmean")
    values = offset_stat_df.loc[offset_stat_df["grp"] == "a", "val"]
    assert stats["a"]["center"] == pytest.approx(float(np.exp(np.log(values).mean())))
    # The geometric mean sits below the arithmetic mean for a spread distribution.
    assert stats["a"]["center"] < values.mean()


def test_offset_stats_gmean_drops_non_positive_values():
    df = pd.DataFrame({"val": [0.0, -5.0, 1.0, 100.0], "grp": ["a"] * 4})
    stats = p._offset_stats(df, "val", cols="grp", cats=["a"], stat="gmean")
    assert stats["a"]["center"] == pytest.approx(10.0) # sqrt(1 * 100), with 0 & -5 left out


def test_offset_stats_without_cols_summarizes_every_row(offset_stat_df):
    stats = p._offset_stats(offset_stat_df, "val", stat="mean")
    assert list(stats) == [None]
    assert stats[None]["center"] == pytest.approx(offset_stat_df["val"].mean())


def test_offset_stats_empty_group_is_none():
    df = pd.DataFrame({"val": [1.0, 2.0], "grp": ["a", "a"]})
    stats = p._offset_stats(df, "val", cols="grp", cats=["a", "b"], stat="mean")
    assert stats["b"] is None


def test_offset_stats_invalid_stat_raises(offset_stat_df):
    with pytest.raises(ValueError, match="Invalid stat"):
        p._offset_stats(offset_stat_df, "val", cols="grp", stat="median")


def test_offset_stats_invalid_errorbar_raises(offset_stat_df):
    with pytest.raises(ValueError, match="Invalid errorbar"):
        p._offset_stats(offset_stat_df, "val", cols="grp", errorbar="ci")


def test_offset_stats_reps_center_and_sd_come_from_the_replicates(offset_stat_df):
    stats = p._offset_stats(offset_stat_df, "val", cols="grp", cats=["a"], stat="mean",
                            reps="rep", errorbar="sd")
    rep_means = offset_stat_df[offset_stat_df["grp"] == "a"].groupby("rep")["val"].mean()
    assert stats["a"]["n"] == len(rep_means) # n counts replicates, not values
    assert stats["a"]["center"] == pytest.approx(rep_means.mean())
    spread = rep_means.std(ddof=1)
    assert stats["a"]["lower"] == pytest.approx(rep_means.mean() - spread)
    assert stats["a"]["upper"] == pytest.approx(rep_means.mean() + spread)


def test_offset_stats_reps_sem_divides_the_sd_by_root_n(offset_stat_df):
    sd = p._offset_stats(offset_stat_df, "val", cols="grp", cats=["a"], stat="mean",
                         reps="rep", errorbar="sd")["a"]
    sem = p._offset_stats(offset_stat_df, "val", cols="grp", cats=["a"], stat="mean",
                          reps="rep", errorbar="sem")["a"]
    assert sem["center"] == pytest.approx(sd["center"])
    assert (sem["upper"] - sem["center"]) == pytest.approx((sd["upper"] - sd["center"]) / math.sqrt(sd["n"]))


def test_offset_stats_reps_gmean_spread_is_multiplicative(offset_stat_df):
    stats = p._offset_stats(offset_stat_df, "val", cols="grp", cats=["a"], stat="gmean",
                            reps="rep", errorbar="sd")["a"]
    rep_gmeans = (offset_stat_df[offset_stat_df["grp"] == "a"]
                  .groupby("rep")["val"].apply(lambda s: float(np.exp(np.log(s).mean()))))
    logs = np.log(rep_gmeans)
    assert stats["center"] == pytest.approx(float(np.exp(logs.mean())))
    assert stats["lower"] == pytest.approx(float(np.exp(logs.mean() - logs.std(ddof=1))))
    # Asymmetric in linear space, equal factors either way.
    assert (stats["center"] / stats["lower"]) == pytest.approx(stats["upper"] / stats["center"])


def test_offset_stats_reps_errorbar_none_leaves_out_the_bounds(offset_stat_df):
    stats = p._offset_stats(offset_stat_df, "val", cols="grp", cats=["a"], stat="mean",
                            reps="rep", errorbar=None)["a"]
    assert stats["lower"] is None and stats["upper"] is None
    rep_means = offset_stat_df[offset_stat_df["grp"] == "a"].groupby("rep")["val"].mean()
    assert stats["center"] == pytest.approx(rep_means.mean())


def test_offset_stats_single_replicate_has_no_spread():
    df = pd.DataFrame({"val": [1.0, 3.0, 5.0], "grp": ["a"] * 3, "rep": ["r1"] * 3})
    stats = p._offset_stats(df, "val", cols="grp", cats=["a"], stat="mean", reps="rep")["a"]
    assert stats["n"] == 1
    assert stats["lower"] == stats["upper"] == pytest.approx(3.0)


def test_dist_offset_stat_adds_one_bar_panel_per_facet_column(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp",
                       facetx="line", facety="block", offset_stat="mean", show=False)
    # The returned axes stay the facet grid; the bar panels are the odd figure axes.
    assert axes.shape == (2, 2)
    assert len(fig.axes) == 8
    assert all(ax not in list(axes.ravel()) for ax in fig.axes[1::2])
    plt.close(fig)


def test_dist_without_offset_stat_has_no_extra_panels(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp",
                       facetx="line", facety="block", show=False)
    assert axes.shape == (2, 2)
    assert len(fig.axes) == 4
    plt.close(fig)


def test_dist_offset_stat_bars_rest_on_the_baselines(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp",
                       cols_order=["a", "b"], offset_stat="mean", offset_stat_graph="bar", show=False)
    dist_ax, stat_ax = fig.axes
    # First cols value on top, half a trace height apart, matching the offset graph's ticks.
    assert list(dist_ax.get_yticks()) == [0.5, 0.0]
    # The default shift is half the bar height, so a bar starts on its baseline the way its trace does.
    assert [rect.get_y() for rect in stat_ax.patches] == pytest.approx([0.5, 0.0])
    values = [value for _, value in _bars(stat_ax)]
    expected = offset_stat_df.groupby("grp")["val"].mean()
    assert values[0] == pytest.approx(expected["a"])
    assert values[1] == pytest.approx(expected["b"])
    assert len(stat_ax.collections) == 0 # No replicates, so no error bars
    plt.close(fig)


def test_dist_offset_stat_shift_zero_centers_bars_on_the_baselines(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_graph="bar", offset_stat_shift=0, show=False)
    centers = [center for center, _ in _bars(fig.axes[1])]
    assert centers == pytest.approx(list(fig.axes[0].get_yticks()))
    plt.close(fig)


def test_dist_offset_stat_shift_is_in_offset_spacing_units(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_graph="bar", offset_stat_shift=1, show=False)
    # 1 whole spacing of lift puts each bar on the next baseline up.
    assert [center for center, _ in _bars(fig.axes[1])] == pytest.approx([1.0, 0.5])
    plt.close(fig)


def test_dist_offset_stat_bars_share_the_trace_colors(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp",
                       cols_order=["a", "b"], offset_stat="mean", offset_stat_graph="bar",
                       palette_or_cmap=["#4EC569", "#FF0000"], show=False)
    bars = fig.axes[1].patches
    assert [matplotlib.colors.to_hex(rect.get_facecolor()) for rect in bars] == ["#4ec569", "#ff0000"]
    plt.close(fig)


def test_dist_offset_stat_auto_picks_gmean_on_a_log_scale(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       x_axis_scale="log", offset_stat="auto", offset_stat_graph="bar", show=False)
    values = [value for _, value in _bars(fig.axes[1])]
    group = offset_stat_df.loc[offset_stat_df["grp"] == "a", "val"]
    assert values[0] == pytest.approx(float(np.exp(np.log(group).mean())))
    assert fig.axes[1].get_xlabel() == "Geometric mean"
    plt.close(fig)


def test_dist_offset_stat_auto_picks_mean_on_a_linear_scale(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat=True, offset_stat_graph="bar", show=False)
    values = [value for _, value in _bars(fig.axes[1])]
    assert values[0] == pytest.approx(offset_stat_df.loc[offset_stat_df["grp"] == "a", "val"].mean())
    assert fig.axes[1].get_xlabel() == "Mean"
    plt.close(fig)


def test_dist_offset_stat_log_panel_gets_whole_decade_ticks(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", x_axis_scale="log",
                       offset_stat="gmean", offset_stat_graph="bar", offset_stat_dims=(0, 0),
                       show=False)
    stat_ax = fig.axes[1]
    assert stat_ax.get_xscale() == "log"
    # Limits land on powers of 10 so a bar starting at the axis still has length to read.
    lower, upper = stat_ax.get_xlim()
    assert math.log10(lower) == pytest.approx(round(math.log10(lower)))
    assert math.log10(upper) == pytest.approx(round(math.log10(upper)))
    assert all(rect.get_x() == pytest.approx(lower) for rect in stat_ax.patches)
    assert all(rect.get_width() > 0 for rect in stat_ax.patches)
    assert _visible_ticks(stat_ax, "x") == [10 ** n for n in range(round(math.log10(lower)),
                                                                  round(math.log10(upper)) + 1)]
    plt.close(fig)


def test_dist_offset_stat_linear_panel_starts_its_bars_at_zero(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp",
                       offset_stat="mean", offset_stat_graph="bar", offset_stat_dims=(0, 0),
                       show=False)
    stat_ax = fig.axes[1]
    assert stat_ax.get_xscale() == "linear"
    assert stat_ax.get_xlim()[0] == 0
    assert all(rect.get_x() == 0 for rect in stat_ax.patches)
    plt.close(fig)


def test_dist_offset_stat_panels_share_one_set_of_limits(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", facetx="line",
                       facety="block", offset_stat="mean", offset_stat_dims=(0, 0), show=False)
    limits = {ax.get_xlim() for ax in fig.axes[1::2]}
    assert len(limits) == 1
    plt.close(fig)


def test_dist_offset_stat_dims_match_the_offset_graph_on_a_linear_scale(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp",
                       offset_stat="mean", offset_stat_graph="bar",
                       offset_stat_dims="match", show=False)
    dist_ax, stat_ax = fig.axes
    assert stat_ax.get_xlim() == pytest.approx(dist_ax.get_xlim())


def test_dist_offset_stat_dims_match_in_x_units_on_a_log_scale(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp",
                       x_axis_scale="log", offset_stat="gmean", offset_stat_graph="bar",
                       offset_stat_dims="match", show=False)
    dist_ax, stat_ax = fig.axes
    # The offset graph holds log10 values on a linear axis; the bar panel is a real log axis.
    assert stat_ax.get_xscale() == "log"
    assert list(stat_ax.get_xlim()) == pytest.approx([10 ** limit for limit in dist_ax.get_xlim()])
    plt.close(fig)


def test_dist_offset_stat_dims_match_follows_per_panel_limits(offset_stat_df):
    # 1 set of x limits per facetx column turns panel sharing off, so each pair stands alone.
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", facetx="line",
                       x_axis_dims=[(0, 500), (0, 5000)], offset_stat="mean",
                       offset_stat_graph="bar", offset_stat_dims="match", show=False)
    assert fig.axes[1].get_xlim() == pytest.approx((0, 500))
    assert fig.axes[3].get_xlim() == pytest.approx((0, 5000))
    plt.close(fig)


def test_dist_offset_stat_width_defaults_to_two_thirds_one_third(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp",
                       offset_stat="mean", show=False)
    dist_width = axes[0, 0].get_position().width
    stat_width = fig.axes[1].get_position().width
    assert stat_width == pytest.approx(dist_width / 2, rel=0.05)
    assert dist_width / (dist_width + stat_width) == pytest.approx(2 / 3, rel=0.05)
    plt.close(fig)


def test_dist_offset_stat_invalid_dims_raise(offset_stat_df):
    with pytest.raises(ValueError, match="Invalid offset_stat_dims"):
        p.dist("half_offset", offset_stat_df, x="val", cols="grp", offset_stat="mean",
               offset_stat_dims="auto", show=False)


def test_dist_offset_stat_dims_are_respected(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp",
                       offset_stat="mean", offset_stat_dims=(0, 5000), show=False)
    assert fig.axes[1].get_xlim() == (0, 5000)
    plt.close(fig)


def test_dist_offset_stat_width_sets_the_panel_ratio(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp",
                       offset_stat="mean", offset_stat_width=0.25, show=False)
    dist_width = axes[0, 0].get_position().width
    assert fig.axes[1].get_position().width == pytest.approx(dist_width * 0.25, rel=0.05)
    plt.close(fig)


def test_dist_offset_stat_height_sets_the_bar_thickness(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp",
                       offset_stat="mean", offset_stat_height=0.5, show=False)
    assert all(rect.get_height() == pytest.approx(0.25) for rect in fig.axes[1].patches) # 0.5 spacing
    plt.close(fig)


def test_dist_offset_stat_axis_overrides_the_label(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp",
                       offset_stat="gmean", offset_stat_axis="gMFI", show=False)
    assert fig.axes[1].get_xlabel() == "gMFI"
    plt.close(fig)


def test_dist_offset_stat_hides_the_bar_panel_of_an_empty_facet(offset_stat_df):
    df = offset_stat_df[~((offset_stat_df["line"] == "L2") & (offset_stat_df["block"] == "bottom"))]
    fig, axes = p.dist("half_offset", df, x="val", cols="grp",
                       facetx="line", facety="block", offset_stat="mean", show=False)
    blank = [ax for ax in fig.axes if not ax.get_visible()]
    assert len(blank) == 2 # The distribution panel & its bar panel both go
    plt.close(fig)


def test_dist_offset_stat_without_cols_draws_one_bar(offset_stat_df):
    fig, axes = p.dist("offset", offset_stat_df, x="val", offset_stat="mean",
                       offset_stat_graph="bar", show=False)
    bars = _bars(fig.axes[1])
    assert len(bars) == 1
    assert bars[0][1] == pytest.approx(offset_stat_df["val"].mean())
    plt.close(fig)


def test_dist_offset_stat_on_a_non_offset_graph_raises(offset_stat_df):
    with pytest.raises(ValueError, match="offset_stat only applies"):
        p.dist("hist", offset_stat_df, x="val", cols="grp", offset_stat="mean", show=False)


def test_dist_offset_stat_invalid_value_raises(offset_stat_df):
    with pytest.raises(ValueError, match="Invalid offset_stat"):
        p.dist("half_offset", offset_stat_df, x="val", cols="grp", offset_stat="median", show=False)


def test_dist_offset_stat_invalid_errorbar_raises(offset_stat_df):
    with pytest.raises(ValueError, match="Invalid offset_stat_errorbar"):
        p.dist("half_offset", offset_stat_df, x="val", cols="grp", offset_stat="mean",
               offset_stat_reps="rep", offset_stat_errorbar="ci", show=False)


def test_dist_offset_stat_reps_draws_an_error_bar_per_bar(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_graph="bar", offset_stat_reps="rep",
                       show=False)
    stat_ax = fig.axes[1]
    assert len(stat_ax.patches) == 2
    # barh is called per baseline, so each bar brings its own container plus an error bar container.
    kinds = [type(container).__name__ for container in stat_ax.containers]
    assert kinds.count("BarContainer") == 2
    assert kinds.count("ErrorbarContainer") == 2
    plt.close(fig)


def test_dist_offset_stat_reps_error_bar_spans_the_replicate_sd(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_graph="bar", offset_stat_reps="rep",
                       show=False)
    rep_means = offset_stat_df[offset_stat_df["grp"] == "a"].groupby("rep")["val"].mean()
    # The first cols value is the top baseline; its error bar segment runs mean +/- sd.
    segments = [line.get_segments()[0] for line in fig.axes[1].collections]
    top = max(segments, key=lambda seg: seg[0][1])
    assert top[0][0] == pytest.approx(rep_means.mean() - rep_means.std(ddof=1))
    assert top[1][0] == pytest.approx(rep_means.mean() + rep_means.std(ddof=1))
    assert top[0][1] == pytest.approx(0.625) # Rides with its bar, 0.25 spacing above the baseline
    plt.close(fig)


def test_dist_offset_stat_reps_bar_summarizes_replicates_not_pooled_values(offset_stat_df):
    # Unequal replicate sizes make the mean of the replicate means differ from the pooled mean.
    df = offset_stat_df[~((offset_stat_df["grp"] == "a") & (offset_stat_df["rep"].str.endswith("r3"))
                          & (offset_stat_df.index % 2 == 0))]
    fig, axes = p.dist("half_offset", df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_graph="bar", offset_stat_reps="rep",
                       show=False)
    values = [value for _, value in _bars(fig.axes[1])]
    rep_means = df[df["grp"] == "a"].groupby("rep")["val"].mean()
    assert values[0] == pytest.approx(rep_means.mean())
    assert values[0] != pytest.approx(df.loc[df["grp"] == "a", "val"].mean())
    plt.close(fig)


def test_dist_offset_stat_reps_errorbar_none_draws_no_error_bars(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp",
                       offset_stat="mean", offset_stat_graph="bar", offset_stat_reps="rep",
                       offset_stat_errorbar=None, show=False)
    assert len(fig.axes[1].collections) == 0
    plt.close(fig)


def test_dist_offset_stat_axis_leaves_room_for_the_error_bars(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", offset_stat="mean",
                       offset_stat_graph="bar", offset_stat_reps="rep", offset_stat_dims=(0, 0),
                       show=False)
    stat_ax = fig.axes[1]
    upper = max(seg[0][1][0] for seg in (line.get_segments() for line in stat_ax.collections))
    assert stat_ax.get_xlim()[1] >= upper
    plt.close(fig)


def test_grid_combines_panels_and_saves(tmp_path, dist_df, cat_df):
    fig, panel_axes = p.grid(
        panels=[
            {"func": p.dist, "graph": "half_offset", "df": dist_df, "x": "val", "cols": "grp"},
            {"func": p.cat, "graph": "bar", "df": cat_df, "x": "grp", "y": "val"},
        ],
        width_ratios=[2, 1], file=str(tmp_path / "grid.png"), show=False,
    )
    assert isinstance(fig, plt.Figure)
    assert len(panel_axes) == 2
    assert (tmp_path / "grid.png").exists()
    plt.close(fig)


def test_grid_panel_keeps_its_own_faceting(dist_df):
    df = dist_df.assign(rep=["r1", "r2"] * (len(dist_df) // 2))
    fig, panel_axes = p.grid(
        panels=[
            {"func": p.dist, "graph": "kde", "df": df, "x": "val", "cols": "grp", "facetx": "rep"},
            {"func": p.cat, "graph": "box", "df": df, "x": "grp", "y": "val"},
        ],
        show=False,
    )
    assert len(panel_axes[0]) == 2  # One axes per facet value inside the panel
    assert len(panel_axes[1]) == 1
    plt.close(fig)


def test_grid_figure_legend_folds_panel_legends(dist_df, cat_df):
    fig, _ = p.grid(
        panels=[
            {"func": p.dist, "graph": "kde", "df": dist_df, "x": "val", "cols": "grp"},
            {"func": p.cat, "graph": "bar", "df": cat_df, "x": "grp", "y": "val"},
        ],
        legend_mode="figure", show=False,
    )
    assert len(fig.legends) == 1
    assert all(len(subfig.legends) == 0 for subfig in fig.subfigs)
    plt.close(fig)


def test_grid_leaves_unused_cells_blank(dist_df):
    fig, panel_axes = p.grid(
        panels=[{"func": p.dist, "graph": "hist", "df": dist_df, "x": "val"}] * 3,
        nrows=2, ncols=2, show=False,
    )
    assert len(panel_axes) == 3
    assert fig.subfigs[-1].get_visible() is False
    plt.close(fig)


def test_graph_method_given_a_figure_does_not_save(tmp_path, dist_df):
    # A method drawing into someone else's figure leaves saving & showing to them.
    fig = plt.figure()
    subfig = fig.subfigures(1, 1)
    p.dist("hist", dist_df, x="val", fig=subfig, file=str(tmp_path / "panel.png"), show=False)
    assert not (tmp_path / "panel.png").exists()
    assert len(subfig.axes) == 1
    plt.close(fig)


def test_grid_rejects_bad_panels_and_interactive_files(dist_df):
    panel = {"func": p.dist, "graph": "hist", "df": dist_df, "x": "val"}
    with pytest.raises(ValueError):
        p.grid(panels=[], show=False)
    with pytest.raises(ValueError):
        p.grid(panels=[{"graph": "hist", "df": dist_df, "x": "val"}], show=False)
    with pytest.raises(ValueError):
        p.grid(panels=[panel] * 3, nrows=1, ncols=2, show=False)
    with pytest.raises(ValueError):
        p.grid(panels=[panel], file="grid.html", show=False)
    with pytest.raises(ValueError):
        p.grid(panels=[panel], legend_mode="nope", show=False)
    plt.close("all")


# Offset distribution statistic panels: cat() graph types & the % difference panel


@pytest.fixture
def offset_ref_df():
    """2 groups sharing the same 3 replicate labels, so a baseline can be paired with a reference."""
    rng = np.random.default_rng(11)
    rows = []
    for grp, center in [("a", 100.0), ("ref", 50.0)]:
        for rep in ["r1", "r2", "r3"]:
            for value in rng.lognormal(mean=math.log(center), sigma=0.2, size=30):
                rows.append({"val": value, "grp": grp, "rep": rep})
    return pd.DataFrame(rows)


@pytest.fixture
def offset_tied_df():
    """3 replicates holding identical values, so their statistics land on top of each other."""
    rows = []
    for grp in ["a", "b"]:
        for rep in ["r1", "r2", "r3"]:
            for value in [10.0, 20.0, 30.0]:
                rows.append({"val": value, "grp": grp, "rep": rep})
    return pd.DataFrame(rows)


def _points(ax):
    """(x, y) of every scatter point, in drawing order."""
    return [tuple(point) for coll in ax.collections
            if type(coll).__name__ == "PathCollection"
            for point in coll.get_offsets()]


def _hand_stats(**reps):
    """A _offset_stats() style dictionary built by hand, so the arithmetic is exact."""
    return {cat: {"center": float(np.mean(list(values.values()))), "lower": None, "upper": None,
                  "n": len(values), "reps": values}
            for cat, values in reps.items()}


def test_offset_stats_hands_back_the_replicate_statistics(offset_stat_df):
    stats = p._offset_stats(offset_stat_df, x="val", cols="grp", stat="mean", reps="rep")
    expected = offset_stat_df[offset_stat_df["grp"] == "a"].groupby("rep", sort=False)["val"].mean()
    assert stats["a"]["reps"] == pytest.approx(expected.to_dict())
    assert stats["a"]["n"] == len(expected)


def test_offset_stats_without_reps_holds_the_pooled_statistic(offset_stat_df):
    stats = p._offset_stats(offset_stat_df, x="val", cols="grp", stat="mean")
    assert list(stats["a"]["reps"]) == [None]
    assert stats["a"]["reps"][None] == pytest.approx(stats["a"]["center"])


def test_offset_stats_relative_puts_the_reference_at_zero():
    stats = _hand_stats(a={"r1": 120.0, "r2": 150.0, "r3": 180.0},
                        ref={"r1": 100.0, "r2": 100.0, "r3": 100.0})
    relative = p._offset_stats_relative(stats, ["ref"], stat="mean")
    assert relative["ref"]["center"] == pytest.approx(0)
    assert relative["a"]["reps"] == pytest.approx({"r1": 20.0, "r2": 50.0, "r3": 80.0})
    assert relative["a"]["center"] == pytest.approx(50.0)


def test_offset_stats_relative_pairs_each_replicate_with_its_own_reference():
    # A reference that moves replicate to replicate cancels out, so every % difference is the same.
    stats = _hand_stats(a={"r1": 200.0, "r2": 400.0}, ref={"r1": 100.0, "r2": 200.0})
    relative = p._offset_stats_relative(stats, ["ref"], stat="mean")
    assert relative["a"]["reps"] == pytest.approx({"r1": 100.0, "r2": 100.0})
    assert relative["a"]["lower"] == pytest.approx(relative["a"]["upper"]) # No spread left


def test_offset_stats_relative_spread_is_the_sd_of_the_percent_differences():
    stats = _hand_stats(a={"r1": 120.0, "r2": 150.0, "r3": 180.0},
                        ref={"r1": 100.0, "r2": 100.0, "r3": 100.0})
    relative = p._offset_stats_relative(stats, ["ref"], stat="mean", errorbar="sd")
    spread = float(np.std([20.0, 50.0, 80.0], ddof=1))
    assert relative["a"]["lower"] == pytest.approx(50.0 - spread)
    assert relative["a"]["upper"] == pytest.approx(50.0 + spread)


def test_offset_stats_relative_sem_divides_the_sd_by_root_n():
    stats = _hand_stats(a={"r1": 120.0, "r2": 150.0, "r3": 180.0},
                        ref={"r1": 100.0, "r2": 100.0, "r3": 100.0})
    sd = p._offset_stats_relative(stats, ["ref"], stat="mean", errorbar="sd")["a"]
    sem = p._offset_stats_relative(stats, ["ref"], stat="mean", errorbar="sem")["a"]
    assert (sem["upper"] - sem["center"]) == pytest.approx((sd["upper"] - sd["center"]) / math.sqrt(3))


def test_offset_stats_relative_combines_several_references_geometrically():
    stats = _hand_stats(a={"r1": 400.0}, lo={"r1": 100.0}, hi={"r1": 400.0})
    relative = p._offset_stats_relative(stats, ["lo", "hi"], stat="gmean")
    assert relative["a"]["center"] == pytest.approx(100.0) # 400 against the gmean of 100 & 400, which is 200


def test_offset_stats_relative_combines_several_references_arithmetically():
    stats = _hand_stats(a={"r1": 400.0}, lo={"r1": 100.0}, hi={"r1": 300.0})
    relative = p._offset_stats_relative(stats, ["lo", "hi"], stat="mean")
    assert relative["a"]["center"] == pytest.approx(100.0) # 400 against the mean of 100 & 300, which is 200


def test_offset_stats_relative_falls_back_to_the_pooled_reference():
    # 'a' has a replicate the reference never saw, so that replicate reads against the whole reference.
    stats = _hand_stats(a={"r1": 200.0, "r2": 300.0}, ref={"r1": 100.0})
    relative = p._offset_stats_relative(stats, ["ref"], stat="mean")
    assert relative["a"]["reps"] == pytest.approx({"r1": 100.0, "r2": 200.0})


def test_offset_stats_relative_errorbar_none_leaves_out_the_bounds():
    stats = _hand_stats(a={"r1": 120.0, "r2": 180.0}, ref={"r1": 100.0, "r2": 100.0})
    relative = p._offset_stats_relative(stats, ["ref"], stat="mean", errorbar=None)
    assert relative["a"]["lower"] is None and relative["a"]["upper"] is None


def test_offset_stats_relative_missing_reference_is_all_none():
    stats = _hand_stats(a={"r1": 120.0})
    assert p._offset_stats_relative(stats, ["ref"], stat="mean") == {"a": None}


def test_offset_stats_relative_invalid_errorbar_raises():
    with pytest.raises(ValueError, match="Invalid errorbar"):
        p._offset_stats_relative(_hand_stats(a={"r1": 1.0}), ["a"], errorbar="ci")


def test_dist_offset_stat_defaults_to_a_swarm(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_reps="rep", show=False)
    stat_ax = fig.axes[1]
    assert len(stat_ax.patches) == 0 # Points, not bars
    reps = offset_stat_df.groupby(["grp", "rep"]).ngroups
    assert len(_points(stat_ax)) == reps # 1 point per replicate per baseline
    plt.close(fig)


def test_dist_offset_stat_swarm_draws_the_replicate_statistics(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_reps="rep", show=False)
    drawn = sorted(x for x, _ in _points(fig.axes[1]))
    expected = sorted(offset_stat_df.groupby(["grp", "rep"], sort=False)["val"].mean())
    assert drawn == pytest.approx(expected)
    plt.close(fig)


def test_dist_offset_stat_swarm_spreads_overlapping_points_off_the_baseline(offset_tied_df):
    fig, axes = p.dist("half_offset", offset_tied_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_reps="rep", show=False)
    top = sorted(point for point in _points(fig.axes[1]) if point[1] > 0.375)
    assert len({round(y, 6) for _, y in top}) == 3 # The 3 identical values step apart
    assert len({round(x, 6) for x, _ in top}) == 1 # & keep their value
    plt.close(fig)


def test_dist_offset_stat_swarm_leaves_a_lone_point_on_its_baseline(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", show=False) # No replicates, so 1 point per baseline
    ys = sorted(y for _, y in _points(fig.axes[1]))
    # half_offset spaces baselines 0.5 apart, & the default shift is a quarter of that spacing.
    assert ys == pytest.approx([0.125, 0.625])
    plt.close(fig)


def test_dist_offset_stat_strip_jitters_within_the_band(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_graph="strip", offset_stat_reps="rep",
                       offset_stat_height=0.5, show=False)
    top = [y for _, y in _points(fig.axes[1]) if y > 0.25]
    assert len(top) == offset_stat_df[offset_stat_df["grp"] == "a"]["rep"].nunique()
    # The band is offset_spacing (0.5) x offset_stat_height (0.5), centered on the shifted baseline.
    assert all(abs(y - 0.625) <= 0.25 / 2 for y in top)
    plt.close(fig)


def test_dist_offset_stat_point_draws_one_marker_with_an_error_bar(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_graph="point", offset_stat_reps="rep",
                       show=False)
    stat_ax = fig.axes[1]
    kinds = [type(container).__name__ for container in stat_ax.containers]
    assert kinds.count("ErrorbarContainer") == 2
    assert len(stat_ax.patches) == 0
    plt.close(fig)


def test_dist_offset_stat_box_draws_a_box_per_baseline(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_graph="box", offset_stat_reps="rep",
                       show=False)
    assert len(fig.axes[1].patches) == 2
    plt.close(fig)


def test_dist_offset_stat_violin_draws_a_filled_shape_per_baseline(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_graph="violin", offset_stat_reps="rep",
                       show=False)
    polys = [coll for coll in fig.axes[1].collections if "Poly" in type(coll).__name__]
    assert len(polys) == 2
    plt.close(fig)


def test_dist_offset_stat_composite_graph_draws_both_layers(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_graph="bar_swarm", offset_stat_reps="rep",
                       show=False)
    stat_ax = fig.axes[1]
    assert len(stat_ax.patches) == 2 # The bars
    assert len(_points(stat_ax)) == offset_stat_df.groupby(["grp", "rep"]).ngroups # & every replicate
    plt.close(fig)


def test_dist_offset_stat_box_falls_back_to_a_point_without_replicates(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", cols_order=["a", "b"],
                       offset_stat="mean", offset_stat_graph="box", show=False)
    stat_ax = fig.axes[1]
    assert len(stat_ax.patches) == 0 # 1 replicate has no quartiles to draw
    kinds = [type(container).__name__ for container in stat_ax.containers]
    assert kinds.count("ErrorbarContainer") == 2
    plt.close(fig)


def test_dist_offset_stat_invalid_graph_raises(offset_stat_df):
    with pytest.raises(ValueError, match="Invalid offset_stat_graph"):
        p.dist("half_offset", offset_stat_df, x="val", cols="grp", offset_stat="mean",
               offset_stat_graph="pie", show=False)


def test_dist_offset_stat_size_sets_the_marker_size(offset_stat_df):
    sizes = []
    for size in (4, 9):
        fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", offset_stat="mean",
                           offset_stat_reps="rep", offset_stat_size=size, show=False)
        sizes.append(fig.axes[1].collections[0].get_sizes()[0])
        plt.close(fig)
    assert sizes == pytest.approx([16, 81]) # s is the marker area, so size squared


def test_dist_offset_stat_dims_default_to_autoscale(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", offset_stat="mean",
                       offset_stat_reps="rep", show=False)
    dist_ax, stat_ax = fig.axes
    # Autoscaled over the statistics, which span far less than the events the traces are drawn from.
    assert stat_ax.get_xlim()[1] < dist_ax.get_xlim()[1]
    assert stat_ax.get_xlim()[1] >= max(x for x, _ in _points(stat_ax))
    plt.close(fig)


def test_dist_offset_stat_ref_adds_a_third_panel(offset_ref_df):
    fig, axes = p.dist("half_offset", offset_ref_df, x="val", cols="grp", cols_order=["a", "ref"],
                       offset_stat="mean", offset_stat_reps="rep", offset_stat_ref="ref",
                       show=False)
    assert len(fig.axes) == 3
    plt.close(fig)


def test_dist_offset_stat_ref_puts_the_reference_at_zero(offset_ref_df):
    fig, axes = p.dist("half_offset", offset_ref_df, x="val", cols="grp", cols_order=["a", "ref"],
                       offset_stat="mean", offset_stat_reps="rep", offset_stat_ref="ref",
                       show=False)
    ref_ax = fig.axes[2]
    bottom = [x for x, y in _points(ref_ax) if y < 0.375] # 'ref' is the lower baseline
    assert bottom == pytest.approx([0, 0, 0])
    expected = (offset_ref_df[offset_ref_df["grp"] == "a"].groupby("rep")["val"].mean()
                / offset_ref_df[offset_ref_df["grp"] == "ref"].groupby("rep")["val"].mean() - 1) * 100
    assert sorted(x for x, y in _points(ref_ax) if y > 0.375) == pytest.approx(sorted(expected))
    plt.close(fig)


def test_dist_offset_stat_ref_panel_marks_zero(offset_ref_df):
    fig, axes = p.dist("half_offset", offset_ref_df, x="val", cols="grp", offset_stat="mean",
                       offset_stat_reps="rep", offset_stat_ref="ref", show=False)
    assert _vlines(fig.axes[2]) == pytest.approx([0])
    plt.close(fig)


def test_dist_offset_stat_ref_panel_is_linear_on_a_log_scale(offset_ref_df):
    fig, axes = p.dist("half_offset", offset_ref_df, x="val", cols="grp", x_axis_scale="log",
                       offset_stat="auto", offset_stat_reps="rep", offset_stat_ref="ref", show=False)
    assert fig.axes[1].get_xscale() == "log" # The statistic keeps the x-axis scale
    assert fig.axes[2].get_xscale() == "linear" # A % difference does not
    plt.close(fig)


def test_dist_offset_stat_ref_axis_defaults_to_percent_difference(offset_ref_df):
    fig, axes = p.dist("half_offset", offset_ref_df, x="val", cols="grp", offset_stat="mean",
                       offset_stat_reps="rep", offset_stat_ref="ref", show=False)
    assert fig.axes[2].get_xlabel() == "% difference"
    plt.close(fig)
    fig, axes = p.dist("half_offset", offset_ref_df, x="val", cols="grp", offset_stat="mean",
                       offset_stat_reps="rep", offset_stat_ref="ref",
                       offset_stat_ref_axis="% vs NIC", show=False)
    assert fig.axes[2].get_xlabel() == "% vs NIC"
    plt.close(fig)


def test_dist_offset_stat_ref_turns_the_statistic_on(offset_ref_df):
    fig, axes = p.dist("half_offset", offset_ref_df, x="val", cols="grp", offset_stat_ref="ref",
                       offset_stat_reps="rep", show=False)
    assert len(fig.axes) == 3
    assert fig.axes[1].get_xlabel() == "Mean" # 'auto' on a linear x-axis
    plt.close(fig)


def test_dist_offset_stat_ref_accepts_several_values(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", offset_stat="mean",
                       offset_stat_ref=["a", "b"], show=False)
    assert len(fig.axes) == 3
    plt.close(fig)


def test_dist_offset_stat_ref_unknown_value_raises(offset_stat_df):
    with pytest.raises(ValueError, match="offset_stat_ref is not among the traces"):
        p.dist("half_offset", offset_stat_df, x="val", cols="grp", offset_stat="mean",
               offset_stat_ref="nope", show=False)


def test_dist_offset_stat_ref_dims_match_raises(offset_stat_df):
    with pytest.raises(ValueError, match="offset_stat_ref_dims cannot match"):
        p.dist("half_offset", offset_stat_df, x="val", cols="grp", offset_stat="mean",
               offset_stat_ref="a", offset_stat_ref_dims="match", show=False)


def test_dist_offset_stat_ref_hides_every_panel_of_an_empty_facet(offset_stat_df):
    df = offset_stat_df[~((offset_stat_df["line"] == "L2") & (offset_stat_df["block"] == "top"))]
    fig, axes = p.dist("half_offset", df, x="val", cols="grp", facetx="line", facety="block",
                       offset_stat="mean", offset_stat_reps="rep", offset_stat_ref="a", show=False)
    hidden = [ax for ax in fig.axes if not ax.get_visible()]
    assert len(hidden) == 3 # The distribution panel & both statistic panels
    plt.close(fig)


def test_dist_offset_stat_ref_width_applies_to_both_panels(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", offset_stat="mean",
                       offset_stat_ref="a", offset_stat_width=0.5, show=False)
    widths = [ax.get_position().width for ax in fig.axes]
    assert widths[1] == pytest.approx(widths[0] / 2, rel=0.05)
    assert widths[2] == pytest.approx(widths[0] / 2, rel=0.05)
    plt.close(fig)


def test_dist_offset_stat_swarm_closes_in_on_the_values(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", x_axis_scale="log",
                       offset_stat="gmean", offset_stat_reps="rep", show=False)
    stat_ax = fig.axes[1]
    values = [x for x, _ in _points(stat_ax)]
    lower, upper = stat_ax.get_xlim()
    assert lower < min(values) and upper > max(values) # Every point is in view
    # & the window hugs them rather than running out to the whole decades around them.
    assert lower > 10 ** math.floor(math.log10(min(values)))
    assert upper < 10 ** math.ceil(math.log10(max(values)))
    plt.close(fig)


def test_dist_offset_stat_swarm_does_not_pin_a_linear_panel_to_zero(offset_tied_df):
    # Every statistic is 20, so a window that closed in on them cannot reach 0 on its own.
    fig, axes = p.dist("half_offset", offset_tied_df, x="val", cols="grp",
                       offset_stat="mean", offset_stat_reps="rep", show=False)
    assert fig.axes[1].get_xlim()[0] > 0 # No bar to run from the axis, so 0 is not held in view
    plt.close(fig)
    fig, axes = p.dist("half_offset", offset_tied_df, x="val", cols="grp", offset_stat="mean",
                       offset_stat_graph="bar", offset_stat_reps="rep", show=False)
    assert fig.axes[1].get_xlim()[0] == 0 # A bar has to run from somewhere
    plt.close(fig)


def test_dist_offset_stat_bar_keeps_zero_in_view(offset_stat_df):
    fig, axes = p.dist("half_offset", offset_stat_df, x="val", cols="grp", offset_stat="mean",
                       offset_stat_graph="bar_swarm", offset_stat_reps="rep", show=False)
    assert fig.axes[1].get_xlim()[0] == 0 # A bar layer anchors the window again
    plt.close(fig)


def test_dist_offset_stat_ref_panel_keeps_zero_in_view(offset_ref_df):
    fig, axes = p.dist("half_offset", offset_ref_df, x="val", cols="grp", offset_stat="mean",
                       offset_stat_reps="rep", offset_stat_ref="ref", show=False)
    lower, upper = fig.axes[2].get_xlim()
    assert lower <= 0 <= upper # The reference it is read against
    plt.close(fig)


# --------------------------------------------------------------------------- #
# log_limits() & the minor tick labels that go with it
# --------------------------------------------------------------------------- #
def _decades(dims):
    """How many decades a window spans."""
    return math.log10(dims[1]) - math.log10(dims[0])


def test_log_limits_closes_in_on_values_inside_a_decade():
    """The pair of rounding helpers opens a lopsided window; this is what that cost."""
    values = [8730, 12000, 17614]
    old = (p.round_down_pow_10(min(values)), p.round_up_pow_10(max(values)))
    assert old == (8000, 100000)                      # 1.1 decades for 0.3 decades of data
    dims = p.log_limits(values)
    assert dims[0] < min(values) and dims[1] > max(values)
    assert _decades(dims) == pytest.approx(0.30 + 2 * 0.06, abs=0.01)
    assert 0.3 / _decades(dims) > 0.7                 # the data fills the axis it is drawn on


def test_log_limits_keeps_the_decade_ends_for_wider_values():
    """Past a decade the rounded ends land the majors on both ends of the axis."""
    values = [120, 50000]
    assert p.log_limits(values) == (p.round_down_pow_10(120), p.round_up_pow_10(50000))


def test_log_limits_pad_and_decades_are_adjustable():
    assert _decades(p.log_limits([100, 200], pad=0.5)) == pytest.approx(math.log10(2) + 1, abs=1e-9)
    # A lower threshold hands the same values back to the decade ends
    assert p.log_limits([100, 200], decades=0.2) == (100, 1000)


def test_log_limits_leaves_out_what_a_log_axis_cannot_draw():
    """log10 of 0 or of a negative value is undefined, so they cannot set a window."""
    assert p.log_limits([-5, 0, 100, 300]) == p.log_limits([100, 300])
    assert p.log_limits([float('nan'), 100, 300]) == p.log_limits([100, 300])


def test_log_limits_with_nothing_positive_leaves_the_axis_alone():
    assert p.log_limits([0, -1]) == (0, 0)
    assert p.log_limits([]) == (0, 0)


def test_log_limits_gives_one_value_a_window_of_its_own():
    dims = p.log_limits([1000, 1000])
    assert dims[0] < 1000 < dims[1]
    assert _decades(dims) == pytest.approx(1 + 2 * 0.06, abs=1e-9)


def test_autoscale_xy_takes_its_log_window_from_log_limits():
    df = pd.DataFrame({'x': [1, 2, 3], 'y': [8730, 12000, 17614]})
    _, y_dims = p.autoscale_xy(df, 'x', 'y', (0, 0), (0, 0), 'linear', 'log')
    assert y_dims == p.log_limits(df['y'])


@pytest.mark.parametrize('dims, labeled', [
    ((5e4, 3e5), True),     # 1e5 is the only major in it, so the minors carry the axis
    ((1e2, 1e5), False),    # 4 majors; labeled minors would only crowd them
])
def test_log_ticks_labels_the_minors_only_inside_a_decade(dims, labeled):
    fig, ax = plt.subplots()
    ax.set_xscale('log')
    ax.set_xlim(*dims)
    p.log_ticks(ax, 'x')
    formatter = ax.xaxis.get_minor_formatter()
    assert isinstance(formatter, matplotlib.ticker.LogFormatterSciNotation) == labeled
    assert bool([t for t in ax.get_xminorticklabels() if t.get_text()]) == labeled
    plt.close(fig)


def test_log_ticks_minor_labels_can_be_asked_for_and_sized():
    fig, ax = plt.subplots()
    ax.set_xscale('log')
    ax.set_xlim(1e2, 1e5)
    p.log_ticks(ax, 'x', label_minors=True, minor_size=7)
    assert [t.get_text() for t in ax.get_xminorticklabels() if t.get_text()]
    assert {t.get_fontsize() for t in ax.get_xminorticklabels()} == {7.0}
    plt.close(fig)


def test_log_ticks_leaves_a_linear_axis_alone():
    fig, ax = plt.subplots()
    ax.set_xlim(0, 0.5)   # under a decade wide, but not a log axis
    p.log_ticks(ax, 'x')
    assert not isinstance(ax.xaxis.get_minor_formatter(),
                          matplotlib.ticker.LogFormatterSciNotation)
    plt.close(fig)


def test_scat_log_axis_closes_in_on_its_points():
    """What a gMFI scatter needs: every point inside 1 decade, & ticks to read them against."""
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [5.1e4, 9e4, 1.5e5, 2.4e5]})
    fig, axes = p.scat(graph='scat', df=df, x='x', y='y', y_axis_scale='log', show=False)
    ax = np.ravel(axes)[0]
    assert ax.get_ylim() == pytest.approx(p.log_limits(df['y']))
    assert [t.get_text() for t in ax.get_yminorticklabels() if t.get_text()]
    plt.close(fig)


# --------------------------------------------------------------------------- #
# stack_layers() & shade_color()
# --------------------------------------------------------------------------- #
@pytest.fixture
def layered():
    """2 facets x 3 categories x 2 replicates of 3 measures taken on the same rows."""
    return pd.DataFrame([
        {'Line': line, 'Mutant': mutant, 'Replicate': rep,
         'Precise': base + rep, 'Indels': 2.0, 'Non-indels': 1.0}
        for line in ('A', 'B')
        for mutant, base in (('m1', 10), ('m2', 20), ('m3', 30))
        for rep in (0, 1)
    ])


def test_shade_color_blends_toward_white():
    assert p.shade_color((1, 0, 0), 1) == (1, 0, 0)
    assert p.shade_color((1, 0, 0), 0.5) == (1, 0.5, 0.5)
    assert p.shade_color('red', 0) == (1, 1, 1)


def test_layer_spec_forms_agree():
    for layer in ('Precise', ('Precise',), {'y': 'Precise'}):
        assert p._layer_spec(layer)['y'] == 'Precise'
    assert p._layer_spec(('Precise', 'Precise edit'))['label'] == 'Precise edit'
    assert p._layer_spec('Precise')['label'] == 'Precise'  # the column names itself


def test_layer_spec_spreads_the_shades_and_cycles_the_hatches():
    specs = [p._layer_spec(name, i, 3) for i, name in enumerate(('a', 'b', 'c'))]
    assert [spec['shade'] for spec in specs] == [1.0, 0.625, 0.25]
    assert [spec['hatch'] for spec in specs] == p._LAYER_HATCHES[:3]
    assert p._layer_spec('a', 0, 1)['shade'] == 1.0  # a lone layer keeps its whole color


def test_layer_spec_keeps_a_hatch_that_was_asked_for():
    assert p._layer_spec({'y': 'a', 'hatch': None}, 1, 3)['hatch'] is None
    assert p._layer_spec({'y': 'a', 'hatch': 'xx'}, 0, 3)['hatch'] == 'xx'
    assert p._layer_spec({'y': 'a', 'shade': 0.4}, 0, 3)['shade'] == 0.4


@pytest.mark.parametrize('layer, match', [
    ({'y': 'a', 'colour': 'red'}, r"Invalid layer setting\(s\) \['colour'\]"),
    ({'label': 'a'}, "needs a 'y' column"),
    (('a', 'b', 'c'), 'Invalid layer'),
    (42, 'Invalid layer'),
])
def test_layer_spec_rejects_bad_input(layer, match):
    with pytest.raises(ValueError, match=match):
        p._layer_spec(layer)


def test_stack_layers_draws_one_bar_per_layer_and_category(layered):
    fig, axes = p.stack_layers(df=layered, x='Mutant', layers=['Precise', 'Indels', 'Non-indels'],
                               facetx='Line', show=False)
    ax = axes[0, 0]
    assert len(ax.patches) == 3 * 3                       # 3 categories x 3 layers
    heights = [patch.get_height() for patch in ax.patches]
    assert heights[:3] == [10.5, 20.5, 30.5]              # layer 1: the mean over replicates
    assert heights[3:6] == [2.0, 2.0, 2.0]
    bottoms = [patch.get_y() for patch in ax.patches]
    assert bottoms[3:6] == heights[:3]                    # layer 2 starts on top of layer 1
    assert bottoms[6:9] == [12.5, 22.5, 32.5]
    plt.close(fig)


def test_stack_layers_keeps_every_layer_in_its_category_color(layered):
    """A 2nd palette would cost the category the color it carries in every other plot."""
    fig, axes = p.stack_layers(df=layered, x='Mutant', layers=['Precise', 'Indels'],
                               palette_or_cmap={'m1': (1, 0, 0), 'm2': (0, 1, 0), 'm3': (0, 0, 1)},
                               show=False)
    patches = axes[0, 0].patches
    assert patches[0].get_facecolor()[:3] == (1, 0, 0)                  # full color on the bottom
    assert patches[2].get_facecolor()[:3] == pytest.approx((0, 0, 1))
    assert patches[3].get_facecolor()[:3] == pytest.approx(p.shade_color((1, 0, 0), 0.25))
    assert [patch.get_hatch() for patch in patches[:2]] == [None, None]
    assert patches[3].get_hatch() == p._LAYER_HATCHES[1]
    plt.close(fig)


def test_stack_layers_points_mark_each_replicate_total(layered):
    fig, axes = p.stack_layers(df=layered, x='Mutant', layers=['Precise', 'Indels', 'Non-indels'],
                               facetx='Line', show=False)
    totals = sorted(np.concatenate([c.get_offsets()[:, 1] for c in axes[0, 0].collections]))
    assert totals == [13.0, 14.0, 23.0, 24.0, 33.0, 34.0]  # Precise + 2 + 1, per replicate
    plt.close(fig)


def test_stack_layers_total_column_can_be_named(layered):
    layered['Edited'] = layered['Precise']
    fig, axes = p.stack_layers(df=layered, x='Mutant', layers=['Precise', 'Indels'],
                               total='Edited', facetx='Line', show=False)
    totals = sorted(np.concatenate([c.get_offsets()[:, 1] for c in axes[0, 0].collections]))
    assert totals == [10.0, 11.0, 20.0, 21.0, 30.0, 31.0]
    assert axes[0, 0].get_ylabel() == 'Edited'
    plt.close(fig)


def test_stack_layers_marks_a_category_it_has_no_rows_for(layered):
    """A blank bar reads as a measured 0, which is not what a missing sample is."""
    thinned = layered[~((layered['Line'] == 'B') & (layered['Mutant'] == 'm3'))]
    fig, axes = p.stack_layers(df=thinned, x='Mutant', layers=['Precise', 'Indels'],
                               facetx='Line', show=False)
    assert [text.get_text() for text in axes[0, 1].texts] == ['n.d.']
    assert len(axes[0, 1].patches) == 3 * 2   # the bar is still there, at 0
    assert axes[0, 1].patches[2].get_height() == 0
    assert [text.get_text() for text in axes[0, 0].texts] == []
    plt.close(fig)


def test_stack_layers_legend_names_the_layers_in_grey(layered):
    """The layers are told apart by shade & hatch, so their keys cannot carry a category color."""
    fig, _ = p.stack_layers(df=layered, x='Mutant', layers=[('Precise', 'Precise edit'), 'Indels'],
                            show=False)
    assert len(fig.legends) == 1
    assert [t.get_text() for t in fig.legends[0].get_texts()] == [
        'Precise edit', 'Indels', 'Replicate (total)']
    face = fig.legends[0].legend_handles[0].get_facecolor()[:3]
    assert face[0] == face[1] == face[2]
    plt.close(fig)


def test_stack_layers_without_points_drops_their_legend_entry(layered):
    fig, axes = p.stack_layers(df=layered, x='Mutant', layers=['Precise'], points=False, show=False)
    assert [t.get_text() for t in fig.legends[0].get_texts()] == ['Precise']
    assert len(axes[0, 0].collections) == 0
    plt.close(fig)


def test_stack_layers_writes_a_file(layered, tmp_path):
    out = tmp_path / 'layers.png'
    fig, _ = p.stack_layers(df=layered, x='Mutant', layers=['Precise'], file=str(out), show=False)
    assert out.exists() and out.stat().st_size > 0
    plt.close(fig)


def test_stack_layers_leaves_a_supplied_figure_to_its_caller(layered, tmp_path):
    out = tmp_path / 'not_written.png'
    fig = plt.figure(figsize=(3, 2))
    returned, _ = p.stack_layers(df=layered, x='Mutant', layers=['Precise'], fig=fig,
                                 file=str(out), figsize=(9, 9), show=False)
    assert returned is fig
    assert tuple(fig.get_size_inches()) == (3.0, 2.0)
    assert not out.exists()
    plt.close(fig)


@pytest.mark.parametrize('kwargs, match', [
    (dict(x='Nope', layers=['Precise']), "x 'Nope' is not a column"),
    (dict(x='Mutant', layers=['Nope']), r"Layer column\(s\) \['Nope'\] are not columns"),
    (dict(x='Mutant', layers=[]), 'needs at least 1 layer'),
    (dict(x='Mutant', layers=['Precise'], total='Nope'), "total 'Nope' is not a column"),
])
def test_stack_layers_checks_its_columns(layered, kwargs, match):
    with pytest.raises(ValueError, match=match):
        p.stack_layers(df=layered, show=False, **kwargs)


# --------------------------------------------------------------------------- #
# a deliberately blank axis name
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize('value, blank', [
    (None, True), ('', True), ([], True), ([''], True),
    (' ', False), ('  ', False), ([' '], False), ('Reads (%)', False),
])
def test_is_blank_label_counts_only_none_and_empty(value, blank):
    """' ' is a name the caller asked for; only None & '' mean they set nothing."""
    assert p._is_blank_label(value) is blank


def test_stack_layers_keeps_a_blank_axis_name_blank(layered):
    """The tick labels already name the categories, so ' ' drops the axis name under them."""
    _, axes = p.stack_layers(df=layered, x='Mutant', layers=['Precise'], facetx='Line',
                             x_axis=' ', show=False)
    assert {ax.get_xlabel() for ax in np.ravel(axes)} == {' '}
    plt.close('all')


def test_stack_layers_unset_axis_name_still_falls_back_to_the_column(layered):
    _, axes = p.stack_layers(df=layered, x='Mutant', layers=['Precise'], facetx='Line', show=False)
    assert np.ravel(axes)[0].get_xlabel() == 'Mutant'
    plt.close('all')


@pytest.mark.parametrize('graph_fn, kwargs', [
    (p.scat, dict(graph='scat', x='x', y='y')),
    (p.cat, dict(graph='swarm', x='cat', y='y')),
])
def test_faceted_graphs_keep_a_blank_axis_name_blank(graph_fn, kwargs):
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'cat': ['a', 'b', 'a', 'b'],
                       'y': [1.0, 2.0, 3.0, 4.0], 'Line': ['A', 'A', 'B', 'B']})
    _, axes = graph_fn(df=df, facetx='Line', x_axis=' ', show=False, **kwargs)
    assert {ax.get_xlabel() for ax in np.ravel(axes)} == {' '}
    plt.close('all')
