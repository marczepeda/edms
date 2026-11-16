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
