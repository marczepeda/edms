"""
Tests for edms.bio.fcs

Covers offset_summary(), the half-offset distribution + statistic + % difference figure:
- the panel arithmetic it owns (how many statistic panels a stat/ref combination asks for,
  and the figure size that keeps every distribution panel the same width as the grid grows)
- the axis names it fills in (gMFI for a geometric mean of a log scaled channel; the
  '% diff. vs <ref>' label for one reference versus several)
- flooring non-positive events before a log scale, and leaving them alone on a linear one
- the caller-owns-the-figure contract it inherits from plot.dist()

The expected figure sizes are the hand-tuned sizes from the figures this function was
factored out of: a 2 facet column figure is 13 x 7 inches with no statistic panel,
19.5 x 7 with one and 26 x 7 with two.
"""
import numpy as np
import pandas as pd
import pytest
import matplotlib.pyplot as plt

from edms.bio import fcs


@pytest.fixture
def events():
    """3 replicates x 3 samples x 2 cell lines of a lognormal channel."""
    rng = np.random.default_rng(0)
    frames = [
        pd.DataFrame({'PE-A': 10 ** rng.normal(mu, 0.25, 200), 'Line': line,
                      'Sample': sample, 'Replicate': rep})
        for line in ('A', 'B')
        for sample, mu in (('ctrl', 4.0), ('s1', 4.3), ('s2', 4.6))
        for rep in '123'
    ]
    return pd.concat(frames, ignore_index=True)


def _summary(events, **kwargs):
    kwargs.setdefault('facetx', 'Line')
    fig, axes = fcs.offset_summary(df=events, x='PE-A', cols='Sample',
                                   cols_order=['ctrl', 's1', 's2'], reps='Replicate',
                                   show=False, **kwargs)
    return fig, axes


# --------------------------------------------------------------------------- #
# Panel arithmetic
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize('kwargs, panels_per_facet, figsize', [
    (dict(stat=None), 1, (13.0, 7.0)),            # distribution only
    (dict(), 2, (19.5, 7.0)),                     # + the statistic
    (dict(ref='ctrl'), 3, (26.0, 7.0)),           # + the % difference
])
def test_panel_count_and_figsize_per_stat_combination(events, kwargs, panels_per_facet, figsize):
    fig, axes = _summary(events, **kwargs)
    assert len(fig.axes) == panels_per_facet * 2  # 2 facet columns
    assert axes.shape == (1, 2)                   # dist() hands back the distribution panels
    assert tuple(fig.get_size_inches()) == figsize
    plt.close(fig)


def test_ref_turns_the_statistic_panel_on_by_itself(events):
    """A % difference is of a statistic, so ref alone still buys both panels."""
    fig, _ = _summary(events, stat=None, ref='ctrl')
    assert len(fig.axes) == 3 * 2
    assert tuple(fig.get_size_inches()) == (26.0, 7.0)
    plt.close(fig)


def test_facety_grows_the_figure_down_instead_of_across(events):
    fig, axes = _summary(events, ref='ctrl', facetx=None, facety='Line')
    assert axes.shape == (2, 1)
    assert tuple(fig.get_size_inches()) == (13.0, 14.0)  # 1 column of panels, 2 rows
    plt.close(fig)


def test_stat_width_widens_only_the_statistic_panels(events):
    """The distribution panel keeps its width; a wider statistic panel grows the figure."""
    fig, _ = _summary(events, ref='ctrl', stat_width=1.0)
    assert tuple(fig.get_size_inches()) == (6.5 * 3 * 2, 7.0)
    plt.close(fig)


def test_panel_sets_the_distribution_panel_size(events):
    fig, _ = _summary(events, ref='ctrl', panel=(4.0, 5.0))
    assert tuple(fig.get_size_inches()) == (4.0 * 2 * 2, 5.0)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Axis names
# --------------------------------------------------------------------------- #
def test_log_geometric_mean_panel_is_labeled_gMFI(events):
    fig, _ = _summary(events, ref='ctrl')
    assert [ax.get_xlabel() for ax in fig.axes[:3]] == ['PE-A', 'gMFI', '% diff. vs ctrl']
    plt.close(fig)


def test_linear_scale_does_not_claim_a_gMFI(events):
    """gMFI is a geometric mean of a log scaled channel; an arithmetic mean is not one."""
    fig, _ = _summary(events, ref='ctrl', x_axis_scale='linear')
    assert fig.axes[1].get_xlabel().lower() == 'mean'
    plt.close(fig)


def test_several_references_are_named_as_a_group(events):
    fig, _ = _summary(events, ref=['ctrl', 's1'])
    assert fig.axes[2].get_xlabel() == '% diff. vs references'
    plt.close(fig)


def test_axis_names_can_be_overridden(events):
    fig, _ = _summary(events, ref='ctrl', stat_axis='gMFI (PE)',
                      ref_axis='% diff. vs controls', x_axis='APC- (Singlet) PE-A')
    assert [ax.get_xlabel() for ax in fig.axes[:3]] == [
        'APC- (Singlet) PE-A', 'gMFI (PE)', '% diff. vs controls']
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Flooring
# --------------------------------------------------------------------------- #
def test_non_positive_events_are_floored_not_dropped(events, capsys):
    """A compensated channel holds negatives; dropping them would thin the lowest baselines."""
    events.loc[events.index[:5], 'PE-A'] = -3.0
    before = len(events)
    fig, _ = _summary(events, ref='ctrl', floor=1)
    assert 'Floored 5' in capsys.readouterr().out
    assert len(events) == before          # the caller's frame is untouched
    assert (events['PE-A'] < 0).sum() == 5
    plt.close(fig)


def test_linear_scale_leaves_negative_events_alone(events, capsys):
    events.loc[events.index[:5], 'PE-A'] = -3.0
    fig, _ = _summary(events, x_axis_scale='linear', floor=1)
    assert 'Floored' not in capsys.readouterr().out
    plt.close(fig)


def test_floor_none_skips_flooring(events, capsys):
    """Opting out hands the negatives to log10, which is the footgun floor exists for."""
    events.loc[events.index[:5], 'PE-A'] = -3.0
    with pytest.warns(RuntimeWarning, match='invalid value encountered in log10'):
        fig, _ = _summary(events, floor=None)
    assert 'Floored' not in capsys.readouterr().out
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Figure ownership & validation
# --------------------------------------------------------------------------- #
def test_owning_the_figure_writes_the_file(events, tmp_path):
    out = tmp_path / 'summary.png'
    fig, _ = _summary(events, ref='ctrl', file=str(out))
    assert out.exists() and out.stat().st_size > 0
    plt.close(fig)


def test_a_supplied_figure_is_drawn_into_and_left_to_the_caller(events, tmp_path):
    """Matching plot.dist(): a caller-supplied figure is not sized, laid out or saved here."""
    out = tmp_path / 'not_written.png'
    fig = plt.figure(figsize=(3, 2))
    returned, _ = _summary(events, ref='ctrl', fig=fig, file=str(out), panel=(99, 99))
    assert returned is fig
    assert tuple(fig.get_size_inches()) == (3.0, 2.0)  # panel ignored
    assert len(fig.axes) == 6                          # but the panels were drawn
    assert not out.exists()
    plt.close(fig)


def test_unknown_cols_column_is_rejected_by_name(events):
    with pytest.raises(ValueError, match="cols 'Nope' is not a column"):
        fcs.offset_summary(df=events, x='PE-A', cols='Nope', show=False)


# --------------------------------------------------------------------------- #
# Metric panels
# --------------------------------------------------------------------------- #
@pytest.fixture
def measures():
    """1 row per sample & replicate of 2 per-sample metrics on very different scales."""
    return pd.DataFrame([
        {'Line': line, 'Sample': sample, 'Replicate': rep,
         'Edit %': edit + i, 'Cut %': (edit + i) / 25}
        for line in ('A', 'B')
        for sample, edit in (('ctrl', 0.4), ('s1', 31.0), ('s2', 62.0))
        for i, rep in enumerate('123')
    ])


def _rows(fig):
    """Visible panels grouped into facet rows, each ordered left to right."""
    rows = {}
    for ax in fig.axes:
        if ax.get_visible():
            rows.setdefault(round(ax.get_position().y0, 6), []).append(ax)
    return [sorted(axs, key=lambda ax: ax.get_position().x0)
            for _, axs in sorted(rows.items(), reverse=True)]


def test_metric_panels_are_added_after_the_statistic_panels(events, measures):
    fig, _ = _summary(events, ref='ctrl', metrics=['Edit %', 'Cut %'], metrics_df=measures)
    assert len(fig.axes) == 5 * 2                                  # 2 facet columns of 5 panels
    assert tuple(fig.get_size_inches()) == (6.5 * (1 + 0.5 * 4) * 2, 7.0)
    first, second = _rows(fig)[0][:5], _rows(fig)[0][5:]
    assert [ax.get_xlabel() for ax in first] == ['PE-A', 'gMFI', '% diff. vs ctrl', 'Edit %', 'Cut %']
    assert [ax.get_xlabel() for ax in second] == [ax.get_xlabel() for ax in first]
    plt.close(fig)


def test_each_metric_autoscales_over_its_own_zero_anchored_window(events, measures):
    """2 metrics on 1 scale would hide whichever is smaller, & a bar has to start at the axis."""
    fig, _ = _summary(events, ref='ctrl', metrics=['Edit %', 'Cut %'], metrics_df=measures)
    edit, cut = _rows(fig)[0][3], _rows(fig)[0][4]
    assert edit.get_xlim() == (0, 65)   # 64 at the top of the error bar, rounded to a 5 tick
    assert cut.get_xlim() == (0, 3)     # 2.56, rounded to a 0.5 tick
    assert _rows(fig)[0][8].get_xlim() == edit.get_xlim()  # & the same window in the other facet
    plt.close(fig)


def test_metric_window_can_be_set_outright(events, measures):
    fig, _ = _summary(events, ref='ctrl', metrics_df=measures,
                      metrics=[{'x': 'Edit %', 'dims': (0, 100)}])
    assert _rows(fig)[0][3].get_xlim() == (0, 100)
    plt.close(fig)


def test_metric_spec_forms_agree(events, measures):
    """A column name, a (column, axis) pair & a dictionary are the same panel."""
    for spec, axis in [('Edit %', 'Edit %'), (('Edit %', 'Editing (%)'), 'Editing (%)'),
                       ({'x': 'Edit %', 'axis': 'Editing (%)'}, 'Editing (%)')]:
        fig, _ = _summary(events, ref='ctrl', metrics=[spec], metrics_df=measures)
        panel = _rows(fig)[0][3]
        assert panel.get_xlabel() == axis
        assert panel.get_xlim() == (0, 65)
        plt.close(fig)


def test_only_the_bottom_facet_row_names_the_metric_axis(events, measures):
    fig, _ = _summary(events, ref='ctrl', facetx=None, facety='Line',
                      metrics=['Edit %'], metrics_df=measures)
    top, bottom = _rows(fig)
    assert top[3].get_xlabel() == ''
    assert bottom[3].get_xlabel() == 'Edit %'
    # dist() shares x down each column, so the tick labels are the bottom row's; the metric panel
    # is outside that grid & follows whatever its distribution panel does.
    def labeled(ax):
        return bool(ax.xaxis.get_tick_params(which='major').get('labelbottom', True))
    assert [labeled(top[3]), labeled(bottom[3])] == [labeled(top[0]), labeled(bottom[0])]
    assert [labeled(top[3]), labeled(bottom[3])] == [False, True]
    plt.close(fig)


def test_metric_panel_shares_the_distribution_panel_baselines(events, measures):
    """The marks sit on the ridges they belong to, so a row reads across the whole figure."""
    fig, _ = _summary(events, ref='ctrl', metrics=['Edit %'], metrics_df=measures)
    dist_ax, metric_ax = _rows(fig)[0][0], _rows(fig)[0][3]
    assert metric_ax.get_ylim() == dist_ax.get_ylim()
    assert len(metric_ax.patches) == 3  # 1 bar per baseline
    assert sorted(round(patch.get_y() + patch.get_height() / 2, 6) for patch in metric_ax.patches) == [
        0.125, 0.625, 1.125]  # (n-1-k) * 0.5 spacing, raised by 0.25 of it
    plt.close(fig)


def test_a_baseline_the_metric_has_nothing_for_is_marked(events, measures):
    """A blank would read as 0, which is a measurement rather than a missing one."""
    fig, _ = _summary(events, ref='ctrl', metrics=['Edit %'],
                      metrics_df=measures[~((measures['Line'] == 'A') & (measures['Sample'] == 's2'))])
    missing, drawn = _rows(fig)[0][3], _rows(fig)[0][7]
    assert [text.get_text() for text in missing.texts] == ['n.d.']
    assert len(missing.patches) == 2 and len(drawn.patches) == 3
    plt.close(fig)


def test_metric_width_widens_only_that_panel(events, measures):
    fig, _ = _summary(events, ref='ctrl', metrics_df=measures,
                      metrics=[{'x': 'Edit %', 'width': 1.0}])
    assert tuple(fig.get_size_inches()) == (6.5 * (1 + 0.5 + 0.5 + 1.0) * 2, 7.0)
    plt.close(fig)


def test_every_panel_lands_on_one_grid(events, measures):
    """tight_layout cannot see the added panels, so margins & gaps place all of them in inches."""
    fig, _ = _summary(events, ref='ctrl', facetx=None, facety='Line',
                      metrics=['Edit %', 'Cut %'], metrics_df=measures)
    width, height = fig.get_size_inches()
    rows = _rows(fig)
    assert [len(row) for row in rows] == [5, 5]
    for row in rows:
        assert len({round(ax.get_position().y0, 6) for ax in row}) == 1       # 1 baseline per row
        assert len({round(ax.get_position().height, 6) for ax in row}) == 1
        for left, right in zip(row, row[1:]):                                 # no overlap, 1 gap
            assert round((right.get_position().x0 - left.get_position().x1) * width, 6) == 0.29
    assert round(rows[0][0].get_position().x0 * width, 6) == 2.2              # the left margin
    assert round((1 - rows[0][-1].get_position().x1) * width, 6) == 0.13
    assert round((rows[0][0].get_position().y0 - rows[1][0].get_position().y1) * height, 6) == 0.75
    plt.close(fig)


def test_metrics_are_drawn_into_a_supplied_figure_but_not_saved(events, measures, tmp_path):
    out = tmp_path / 'not_written.png'
    fig = plt.figure(figsize=(12, 6))
    returned, _ = _summary(events, ref='ctrl', metrics=['Edit %'], metrics_df=measures,
                           fig=fig, file=str(out))
    assert returned is fig
    assert len(fig.axes) == 4 * 2
    assert tuple(fig.get_size_inches()) == (12.0, 6.0)
    assert not out.exists()
    plt.close(fig)


@pytest.mark.parametrize('metrics, kwargs, match', [
    (['Nope'], {}, "Metric panel x 'Nope' is not a column"),
    ([{'x': 'Edit %', 'colour': 'red'}], {}, r"Invalid metric panel setting\(s\) \['colour'\]"),
    ([('Edit %', 'Editing (%)', 'extra')], {}, 'Invalid metric panel'),
    ([42], {}, 'Invalid metric panel'),
    (['Edit %'], {'metrics_df': None}, "has no dataframe"),
])
def test_metric_specs_are_checked(events, measures, metrics, kwargs, match):
    kwargs = {'metrics_df': measures, **kwargs}
    with pytest.raises(ValueError, match=match):
        _summary(events, metrics=metrics, **kwargs)


def test_metric_dataframe_must_carry_the_baselines_and_facets(events, measures):
    with pytest.raises(ValueError, match=r"is missing column\(s\) \['Line'\]"):
        _summary(events, metrics=['Edit %'], metrics_df=measures.drop(columns='Line'))


def test_metric_dims_round_out_to_a_readable_tick():
    """_metric_dims() keeps 0 in view & rounds to a half decade, on either side of it."""
    def stats(*values, reps=()):
        return {(0, 0): {str(i): {'center': value, 'lower': None, 'upper': None,
                                  'reps': {str(k): rep for k, rep in enumerate(reps)}}
                         for i, value in enumerate(values)}}
    assert fcs._metric_dims(stats(79.7)) == (0, 80)
    # A swarm draws every replicate, so one outside its own error bar still has to fit
    assert fcs._metric_dims(stats(5.1, reps=(5.1, 5.6))) == (0, 6)
    assert fcs._metric_dims(stats(5.1)) == (0, 5.5)
    assert fcs._metric_dims(stats(-12.0, 4.0)) == (-15, 5)
    assert fcs._metric_dims(stats(0.0)) == (0, 1)        # nothing to scale over
    assert fcs._metric_dims(stats(240.0, 3e4), scale='log') == (100, 1e5)  # whole decades
