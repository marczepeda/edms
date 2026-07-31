"""
Tests for edms.gen.stat

Hand-computed expected values on small inputs where feasible, cross-checked
against scipy/numpy for the statistical-test wrappers.
"""
import numpy as np
import pandas as pd
import pytest
from scipy.stats import (
    ttest_ind, ttest_rel, f_oneway, mannwhitneyu, wilcoxon, fisher_exact,
)

from edms.gen import stat as st


# --------------------------------------------------------------------------- #
# describe()
# --------------------------------------------------------------------------- #

def test_describe_basic_stats():
    df = pd.DataFrame({"a": [1.0, 2.0, 3.0, 4.0]})
    out = st.describe(df, cols=["a"])
    row = out.iloc[0]
    assert row["mean"] == pytest.approx(2.5)
    assert row["median"] == pytest.approx(2.5)
    assert row["min"] == 1.0
    assert row["max"] == 4.0
    assert row["range"] == pytest.approx(3.0)
    assert row["count"] == 4
    assert row["sum"] == pytest.approx(10.0)


def test_describe_from_file(tmp_path):
    df = pd.DataFrame({"a": [1.0, 2.0, 3.0]})
    pt = tmp_path / "data.csv"
    df.to_csv(pt, index=False)
    out = st.describe(str(pt), cols=["a"])
    assert out.iloc[0]["mean"] == pytest.approx(2.0)


def test_describe_save(tmp_path):
    df = pd.DataFrame({"a": [1.0, 2.0, 3.0]})
    st.describe(df, cols=["a"], dir=str(tmp_path), file="out.csv")
    assert (tmp_path / "out.csv").exists()


# --------------------------------------------------------------------------- #
# difference()
# --------------------------------------------------------------------------- #

def test_difference_ttest_ind_two_groups():
    df = pd.DataFrame({
        "val": [1, 2, 3, 10, 11, 12],
        "grp": ["a", "a", "a", "b", "b", "b"],
    })
    out = st.difference(df, data_col="val", compare_col="grp", compare=["a", "b"],
                         same=False, para=True)
    expected_t, expected_p = ttest_ind(df[df.grp == "a"]["val"], df[df.grp == "b"]["val"])
    assert out.iloc[0]["p_value"] == pytest.approx(expected_p)
    assert out.iloc[0]["test"] == "Student's T-test"


def test_difference_ttest_rel_two_groups():
    df = pd.DataFrame({
        "val": [1, 2, 3, 2, 3, 4],
        "grp": ["a", "a", "a", "b", "b", "b"],
    })
    out = st.difference(df, data_col="val", compare_col="grp", compare=["a", "b"],
                         same=True, para=True)
    expected_t, expected_p = ttest_rel(df[df.grp == "a"]["val"], df[df.grp == "b"]["val"])
    assert out.iloc[0]["p_value"] == pytest.approx(expected_p)
    assert out.iloc[0]["test"] == "Paired Student's T-test"


def test_difference_anova_three_groups():
    df = pd.DataFrame({
        "val": [1, 2, 3, 10, 11, 12, 20, 21, 22],
        "grp": ["a"] * 3 + ["b"] * 3 + ["c"] * 3,
    })
    out = st.difference(df, data_col="val", compare_col="grp", compare=["a", "b", "c"],
                         same=False, para=True)
    expected_f, expected_p = f_oneway(
        df[df.grp == "a"]["val"], df[df.grp == "b"]["val"], df[df.grp == "c"]["val"]
    )
    anova_row = out[out["test"] == "1-way Anova Test"].iloc[0]
    assert anova_row["p_value"] == pytest.approx(expected_p)


def test_difference_mannwhitney():
    df = pd.DataFrame({
        "val": [1, 2, 3, 10, 11, 12],
        "grp": ["a", "a", "a", "b", "b", "b"],
    })
    out = st.difference(df, data_col="val", compare_col="grp", compare=["a", "b"],
                         same=False, para=False)
    expected_u, expected_p = mannwhitneyu(
        df[df.grp == "a"]["val"], df[df.grp == "b"]["val"], alternative="two-sided"
    )
    assert out.iloc[0]["p_value"] == pytest.approx(expected_p)


def test_difference_wilcoxon():
    df = pd.DataFrame({
        "val": [1, 2, 3, 4, 2, 3, 5, 6],
        "grp": ["a"] * 4 + ["b"] * 4,
    })
    out = st.difference(df, data_col="val", compare_col="grp", compare=["a", "b"],
                         same=True, para=False)
    expected_w, expected_p = wilcoxon(
        df[df.grp == "a"]["val"], df[df.grp == "b"]["val"]
    )
    assert out.iloc[0]["p_value"] == pytest.approx(expected_p)


def test_difference_invalid_compare_length_raises():
    # `compare` with fewer than 2 entries (para=True, same=False) has no
    # valid statistical test to run, so difference() should raise a clear,
    # catchable ValueError rather than crash with UnboundLocalError.
    df = pd.DataFrame({"val": [1, 2], "grp": ["a", "a"]})
    with pytest.raises(ValueError, match="Invalid compare"):
        st.difference(df, data_col="val", compare_col="grp", compare=["a"], same=False, para=True)


def test_difference_invalid_compare_length_raises_same_true():
    # Same bug/fix applies to the same=True, para=True branch.
    df = pd.DataFrame({"val": [1, 2], "grp": ["a", "a"]})
    with pytest.raises(ValueError, match="Invalid compare"):
        st.difference(df, data_col="val", compare_col="grp", compare=["a"], same=True, para=True)


# --------------------------------------------------------------------------- #
# correlation()
# --------------------------------------------------------------------------- #

def test_correlation_non_tidy_matches_pandas_corr():
    df = pd.DataFrame({"x": [1, 2, 3, 4], "y": [2, 4, 6, 9]})
    out = st.correlation(df, value_cols=["x", "y"], plot=False)
    expected = df[["x", "y"]].corr(method="pearson")
    pd.testing.assert_frame_equal(out, expected)


def test_correlation_tidy_pivot():
    df = pd.DataFrame({
        "sample": ["s1", "s1", "s2", "s2"],
        "gene": ["g1", "g2", "g1", "g2"],
        "val": [1.0, 2.0, 3.0, 4.0],
    })
    out = st.correlation(df, var_cols=["sample", "gene"], value_cols=["val"], plot=False)
    # pivot uses index=sample, columns=gene -> corr() then correlates the
    # resulting *columns*, i.e. the gene columns (g1, g2), each with one
    # value per sample.
    assert set(out.columns) == {"g1", "g2"}


def test_correlation_save(tmp_path):
    df = pd.DataFrame({"x": [1, 2, 3], "y": [1, 2, 3]})
    st.correlation(df, value_cols=["x", "y"], plot=False, dir=str(tmp_path), file_data="corr.csv")
    assert (tmp_path / "corr.csv").exists()


# --------------------------------------------------------------------------- #
# weighted_correlation()
# --------------------------------------------------------------------------- #

def test_weighted_correlation_unweighted_matches_numpy_pearson():
    df = pd.DataFrame({"x": [1.0, 2.0, 3.0, 4.0], "y": [2.0, 1.0, 4.0, 3.0]})
    out = st.weighted_correlation(df, x="x", y="y", method="pearson")
    expected = np.corrcoef(df["x"], df["y"])[0, 1]
    assert out == pytest.approx(expected)


def test_weighted_correlation_all_equal_weights_matches_unweighted():
    df = pd.DataFrame({
        "x": [1.0, 2.0, 3.0, 4.0],
        "y": [2.0, 1.0, 4.0, 3.0],
        "w": [2.0, 2.0, 2.0, 2.0],
    })
    weighted = st.weighted_correlation(df, x="x", y="y", weight="w", method="pearson")
    unweighted = st.weighted_correlation(df, x="x", y="y", method="pearson")
    assert weighted == pytest.approx(unweighted)


def test_weighted_correlation_extreme_weight_dominates():
    # Point 0 given overwhelming weight should pull the correlation toward
    # whatever a 2-point correlation involving mostly that point implies;
    # here we just check the weighted correlation differs from unweighted
    # when weights are highly skewed and non-trivial.
    df = pd.DataFrame({
        "x": [1.0, 2.0, 3.0, -10.0],
        "y": [1.0, 2.0, 3.0, 10.0],
        "w": [1.0, 1.0, 1.0, 1000.0],
    })
    weighted = st.weighted_correlation(df, x="x", y="y", weight="w", method="pearson")
    unweighted = st.weighted_correlation(df, x="x", y="y", method="pearson")
    assert weighted != pytest.approx(unweighted)


def test_weighted_correlation_invalid_method_raises():
    df = pd.DataFrame({"x": [1.0, 2.0], "y": [1.0, 2.0]})
    with pytest.raises(ValueError):
        st.weighted_correlation(df, x="x", y="y", method="bogus")


# --------------------------------------------------------------------------- #
# corr_line() / weighted_corr_line()
# --------------------------------------------------------------------------- #

def test_corr_line_matches_polyfit():
    df = pd.DataFrame({"x": [1.0, 2.0, 3.0, 4.0], "y": [2.0, 4.0, 6.0, 8.0]})
    a, b = st.corr_line(df, x="x", y="y")
    expected_b, expected_a = np.polyfit(df["x"], df["y"], 1)
    assert a == pytest.approx(expected_a)
    assert b == pytest.approx(expected_b)
    # y = 2x exactly
    assert b == pytest.approx(2.0)
    assert a == pytest.approx(0.0, abs=1e-9)


def test_corr_line_draws_on_provided_axis():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    df = pd.DataFrame({"x": [1.0, 2.0, 3.0], "y": [1.0, 2.0, 3.0]})
    n_lines_before = len(ax.lines)
    st.corr_line(df, x="x", y="y", ax=ax)
    assert len(ax.lines) == n_lines_before + 1
    plt.close(fig)


def test_weighted_corr_line_equal_weights_matches_corr_line():
    df = pd.DataFrame({"x": [1.0, 2.0, 3.0, 4.0], "y": [1.0, 3.0, 2.0, 5.0], "w": [1.0] * 4})
    a1, b1 = st.corr_line(df, x="x", y="y")
    a2, b2 = st.weighted_corr_line(df, x="x", y="y", weight="w")
    assert a1 == pytest.approx(a2)
    assert b1 == pytest.approx(b2)


# --------------------------------------------------------------------------- #
# compare()
# --------------------------------------------------------------------------- #

def _compare_input_df():
    # 2 samples per condition, 2 variables. Condition 'ctrl' is the baseline.
    rows = []
    for cond, samples in [("ctrl", ["s1", "s2"]), ("treat", ["s3", "s4"])]:
        for s in samples:
            for var, count in [("v1", 100), ("v2", 50)]:
                rows.append({"sample": s, "cond": cond, "var": var, "count": count})
    return pd.DataFrame(rows)


def test_compare_aggregate_fc_is_one_when_counts_identical():
    df = _compare_input_df()
    out = st.compare(df, sample="sample", cond="cond", cond_comp="ctrl", var="var",
                      count="count", column_prefix="raw")
    # counts are identical between conditions -> FC should be ~1 for every variable
    assert out["raw_FC"].apply(lambda v: v == pytest.approx(1.0, abs=1e-6)).all()
    assert set(out["var"]) == {"v1", "v2"}
    assert (out["raw_compare"] == "ctrl").all()


def test_compare_replicate_mode_runs_and_has_fc():
    rows = []
    for cond in ["ctrl", "treat"]:
        for rep in ["r1", "r2"]:
            for var, count in [("v1", 100), ("v2", 50)]:
                rows.append({"sample": f"{cond}_{rep}", "cond": cond, "replicate": rep,
                             "var": var, "count": count})
    df = pd.DataFrame(rows)
    out = st.compare(df, sample="sample", cond="cond", cond_comp="ctrl", var="var",
                      count="count", replicate="replicate", column_prefix="raw")
    # The replicate-mode summary step names its FC/variability columns
    # without the column_prefix (unlike the aggregate-mode branch), so the
    # actual output column is 'FC_mean', not 'raw_FC'.
    assert "FC_mean" in out.columns
    assert (out["cond"] == "treat").all()
    # counts are identical across replicates/conditions -> FC should be ~1
    assert out["FC_mean"].apply(lambda v: v == pytest.approx(1.0, abs=1e-6)).all()


def test_compare_column_prefix_conflict_raises():
    df = _compare_input_df()
    df["pre_count"] = 1
    with pytest.raises(ValueError):
        st.compare(df, sample="sample", cond="cond", cond_comp="ctrl", var="var",
                   count="count", column_prefix="pre")


def test_compare_negative_pseudocount_raises():
    df = _compare_input_df()
    with pytest.raises(ValueError):
        st.compare(df, sample="sample", cond="cond", cond_comp="ctrl", var="var",
                   count="count", column_prefix="raw", pseudocount=-1)


def test_compare_default_column_prefix_works():
    # Using the documented default column_prefix='' should work for
    # ordinary, valid input (previously always raised ValueError because
    # the rename-collision check compared the count column's own name
    # against itself).
    df = _compare_input_df()
    out = st.compare(df, sample="sample", cond="cond", cond_comp="ctrl", var="var", count="count")
    assert out["FC"].apply(lambda v: v == pytest.approx(1.0, abs=1e-6)).all()
    assert set(out["var"]) == {"v1", "v2"}
    assert (out["compare"] == "ctrl").all()


# --------------------------------------------------------------------------- #
# odds_ratio()
# --------------------------------------------------------------------------- #

def test_odds_ratio_matches_fisher_exact_manual():
    # Simple 2-condition, 2-variable (mut vs WT) table.
    df = pd.DataFrame({
        "cond": ["ctrl", "ctrl", "treat", "treat"],
        "var": ["WT", "mut", "WT", "mut"],
        "count": [100, 10, 80, 40],
    })
    out = st.odds_ratio(df, cond="cond", cond_comp="ctrl", var="var", var_comp="WT",
                         count="count", pseudocount=1)
    row = out[out["var"] == "mut"].iloc[0]

    pc = 1
    table = [[40 + pc, 80 + pc], [10 + pc, 100 + pc]]
    expected_or, expected_p = fisher_exact(table)
    assert row["fisher_exact_odds_ratio"] == pytest.approx(expected_or)
    assert row["fisher_exact_pval"] == pytest.approx(expected_p)
    assert row["log2_fisher_exact_odds_ratio"] == pytest.approx(np.log2(expected_or))


def test_odds_ratio_missing_var_comp_raises():
    df = pd.DataFrame({
        "cond": ["ctrl", "treat"],
        "var": ["mut", "mut"],
        "count": [10, 20],
    })
    with pytest.raises(ValueError):
        st.odds_ratio(df, cond="cond", cond_comp="ctrl", var="var", var_comp="WT", count="count")


# --------------------------------------------------------------------------- #
# zscore()
# --------------------------------------------------------------------------- #

def test_zscore_manual_computation():
    df = pd.DataFrame({
        "cond": ["a", "a", "a", "a", "b", "b", "b", "b"],
        "kind": ["base", "base", "other", "other", "base", "base", "other", "other"],
        "val": [1.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0, 40.0],
    })
    out = st.zscore(df, val="val", cond_col="cond", var_col="kind", var="base")

    base_a = df[(df.cond == "a") & (df.kind == "base")]["val"]
    mu_a, sigma_a = base_a.mean(), base_a.std(ddof=1)
    expected_a_other = (5.0 - mu_a) / sigma_a
    got_a_other = out[(out.cond == "a") & (out.kind == "other")]["val_z"].iloc[0]
    assert got_a_other == pytest.approx(expected_a_other)

    base_b = df[(df.cond == "b") & (df.kind == "base")]["val"]
    mu_b, sigma_b = base_b.mean(), base_b.std(ddof=1)
    expected_b_other = (30.0 - mu_b) / sigma_b
    got_b_other = out[(out.cond == "b") & (out.kind == "other")]["val_z"].iloc[0]
    assert got_b_other == pytest.approx(expected_b_other)


def test_zscore_no_baseline_rows_gives_nan():
    df = pd.DataFrame({
        "cond": ["a", "a"],
        "kind": ["other", "other"],
        "val": [1.0, 2.0],
    })
    out = st.zscore(df, val="val", cond_col="cond", var_col="kind", var="base")
    assert out["val_z"].isna().all()


def test_zscore_missing_columns_raises_keyerror():
    df = pd.DataFrame({"cond": ["a"], "val": [1.0]})
    with pytest.raises(KeyError):
        st.zscore(df, val="val", cond_col="cond", var_col="missing_col", var="base")


def test_zscore_custom_out_col_and_save(tmp_path):
    df = pd.DataFrame({
        "cond": ["a", "a"],
        "kind": ["base", "other"],
        "val": [1.0, 3.0],
    })
    out = st.zscore(df, val="val", cond_col="cond", var_col="kind", var="base",
                     out_col="z_custom", dir=str(tmp_path), file="z.csv")
    assert "z_custom" in out.columns
    assert (tmp_path / "z.csv").exists()
