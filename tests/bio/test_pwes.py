"""
Tests for edms.bio.pwes — Proximity-Weighted Enrichment Score (PWES) clustering
and its plotting functions.

Covers the pure numeric helpers (hex/rgb conversion, gauss, hill,
calculate_pw_score, calculate_pwes, cluster_pwes, shuffle_pwes, fillgaps,
make_cluster_color_dict, _safe_slug, pymol_script, process_pdb/get_pairwise_dist)
with hand-verified expected values, and the matplotlib plotting functions
(hist, cat, torn, heatmap, heatmap_cbar, clustermap) with small synthetic
DataFrames, asserting they run without error.

Special attention is given to torn()'s secondary-structure ("DSSP track")
code path, since a previous bug there referenced an undefined name `FC`
instead of the `scores_col` parameter at the line computing `ss_y`
(`ss_y = min(df[f'log2({FC})'])-2*0.5` -> fixed to
`ss_y = min(df[scores_col])-2*0.5`). A regression test exercises this exact
path (individual=True, DSSP + chain_id supplied) so a reintroduction of that
class of bug would be caught.
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy.cluster.hierarchy import linkage

from edms.bio import pwes


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


# --------------------------------------------------------------------------- #
# HEX <-> RGB conversion
# --------------------------------------------------------------------------- #

def test_hex_to_rgb():
    assert pwes.hex_to_rgb("#FF0000") == (1.0, 0.0, 0.0)
    assert pwes.hex_to_rgb("00FF00") == (0.0, 1.0, 0.0)
    assert pwes.hex_to_rgb("0000FF") == (0.0, 0.0, 1.0)


def test_hex_list_to_rgb_list():
    assert pwes.hex_list_to_rgb_list(["#FF0000", "#0000FF"]) == [(1.0, 0.0, 0.0), (0.0, 0.0, 1.0)]


def test_tuple_to_hex_roundtrip():
    assert pwes.tuple_to_hex((1.0, 0.0, 0.0)) == "FF0000"
    assert pwes.tuple_to_hex((0.0, 1.0, 0.0)) == "00FF00"


def test_tuple_to_hex_rejects_non_float_tuple():
    # ints (not floats) in [0,1] -> guard clause returns None
    assert pwes.tuple_to_hex((0, 0, 0)) is None


def test_tuple_to_hex_rejects_non_tuple():
    assert pwes.tuple_to_hex("not a tuple") is None


# --------------------------------------------------------------------------- #
# fillgaps
# --------------------------------------------------------------------------- #

def test_fillgaps_inserts_missing_positions_as_nan():
    idx = ["A0001", "A0003"]
    df = pd.DataFrame([[1, 2], [3, 4]], index=idx, columns=idx)
    out = pwes.fillgaps(df, {"A": (1, 3)})
    assert out.index.tolist() == ["A0001", "A0002", "A0003"]
    assert out.loc["A0001", "A0001"] == 1
    assert out.loc["A0001", "A0003"] == 2
    assert out.loc["A0003", "A0001"] == 3
    assert out.loc["A0003", "A0003"] == 4
    assert math.isnan(out.loc["A0002", "A0002"])
    assert math.isnan(out.loc["A0001", "A0002"])


# --------------------------------------------------------------------------- #
# process_pdb / get_pairwise_dist (distance factor)
# --------------------------------------------------------------------------- #

_MINI_PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00  0.00           C
ATOM      4  N   GLY A   2       0.000   0.000   0.000  1.00  0.00           N
ATOM      5  CA  GLY A   2       4.000   0.000   0.000  1.00  0.00           C
ATOM      6  N   ALA B   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      7  CA  ALA B   1       1.000   3.000   0.000  1.00  0.00           C
END
"""


@pytest.fixture()
def mini_pdb(tmp_path):
    pt = tmp_path / "mini.pdb"
    pt.write_text(_MINI_PDB)
    return str(pt)


def test_process_pdb_and_get_pairwise_dist(mini_pdb):
    df_centroids = pwes.process_pdb(mini_pdb, chains=["A"])
    assert sorted(df_centroids["aa_num"].astype(int).tolist()) == [1, 2]
    assert set(df_centroids["label"]) == {"A0001", "A0002"}

    dist = pwes.get_pairwise_dist(df_centroids, chains=["A"])
    # CA of residue 1 at (1,0,0), CA of residue 2 at (4,0,0) -> euclidean distance 3.0
    assert dist.loc["A0001", "A0002"] == pytest.approx(3.0)
    assert dist.loc["A0001", "A0001"] == pytest.approx(0.0)


# --------------------------------------------------------------------------- #
# gauss / hill
# --------------------------------------------------------------------------- #

def test_gauss_matches_hand_formula():
    assert pwes.gauss(0, 1) == pytest.approx(1.0)
    assert pwes.gauss(2, 1) == pytest.approx(math.exp(-(2 * 2) / (2 * 1 * 1)))
    assert pwes.gauss(3, 2) == pytest.approx(math.exp(-(3 * 3) / (2 * 2 * 2)))


def test_hill_matches_hand_formula():
    lfc, m, theta = 2, 3, 2
    expected = (lfc ** m) / ((lfc ** m) + (theta ** m))
    assert pwes.hill(lfc, m, theta) == pytest.approx(expected)
    assert pwes.hill(theta, m, theta) == pytest.approx(0.5)  # at lfc==theta, value is always 0.5


# --------------------------------------------------------------------------- #
# calculate_pw_score
# --------------------------------------------------------------------------- #

def test_calculate_pw_score_matches_hand_computed_tanh():
    df_score = pd.DataFrame({"epegRNA_ID": ["a", "b", "c"], "score": [1.0, 2.0, 3.0]})
    out = pwes.calculate_pw_score(df_score, "score", tanh_a=1.0)

    # pairwise sums: upper-triangle {a,b}=3 {a,c}=4 {b,c}=5 -> mean=4, std(ddof=1)=1
    # df_pws[i,j] = tanh(sum-4)
    expected = {
        ("a", "a"): math.tanh(2 - 4),
        ("a", "b"): math.tanh(3 - 4),
        ("a", "c"): math.tanh(4 - 4),
        ("b", "b"): math.tanh(4 - 4),
        ("b", "c"): math.tanh(5 - 4),
        ("c", "c"): math.tanh(6 - 4),
    }
    for (i, j), val in expected.items():
        assert out.loc[i, j] == pytest.approx(val)
        assert out.loc[j, i] == pytest.approx(val)  # symmetric


# --------------------------------------------------------------------------- #
# calculate_pwes
# --------------------------------------------------------------------------- #

def test_calculate_pwes_multiplies_and_zeroes_negatives():
    list_aas = [1, 2, 3]
    df_gauss = pd.DataFrame(
        np.array([[1, 0.5, 0.1], [0.5, 1, 0.5], [0.1, 0.5, 1]]), index=list_aas, columns=list_aas
    )
    df_pws = pd.DataFrame(np.array([[1.0, -2.0, 3.0], [-2.0, 1.0, -1.0], [3.0, -1.0, 1.0]]))

    sorted_out, unsorted_out = pwes.calculate_pwes(df_gauss, df_pws.copy(), list_aas, pos_only=True)
    expected = np.array([[1.0, 0.0, 0.3], [0.0, 1.0, 0.0], [0.3, 0.0, 1.0]])
    assert np.allclose(unsorted_out.values, expected)
    assert np.allclose(sorted_out.values, expected)  # already sorted in this example


def test_calculate_pwes_pos_only_false_keeps_negatives():
    list_aas = [1, 2]
    df_gauss = pd.DataFrame(np.array([[1, 1], [1, 1]]), index=list_aas, columns=list_aas)
    df_pws = pd.DataFrame(np.array([[1.0, -2.0], [-2.0, 1.0]]))
    _, unsorted_out = pwes.calculate_pwes(df_gauss, df_pws.copy(), list_aas, pos_only=False)
    assert unsorted_out.loc[1, 2] == pytest.approx(-2.0)


# --------------------------------------------------------------------------- #
# cluster_pwes / get_clus_aa
# --------------------------------------------------------------------------- #

def test_cluster_pwes_and_get_clus_aa():
    list_aas = [1, 2, 3]
    df_score = pd.DataFrame({"AA Number": [1, 2, 3], "log2(FC)": [1.0, -2.0, 0.5]})
    # 1 & 2 close together, 3 far away -> 2 clusters at a moderate threshold
    df_pws = pd.DataFrame(np.array([[0.0, 1.0, 4.0], [1.0, 0.0, 2.0], [4.0, 2.0, 0.0]]))

    df_clus, link = pwes.cluster_pwes(df_pws, df_score, list_aas, t=5, x_col="AA Number")
    assert "cl_new" in df_clus.columns
    assert len(df_clus) == 3
    assert link.shape == (2, 4)  # scipy linkage matrix for 3 observations: n-1 merges, 4 cols

    aas_dict = pwes.get_clus_aa(df_clus, "AA Number")
    all_aas = sorted(sum(aas_dict.values(), []))
    assert all_aas == [1, 2, 3]


def test_get_clus_aa_missing_column_raises():
    df_clus = pd.DataFrame({"cl_new": [1, 2]})
    with pytest.raises(Exception):
        pwes.get_clus_aa(df_clus, "AA Number")


# --------------------------------------------------------------------------- #
# shuffle_pwes
# --------------------------------------------------------------------------- #

def test_shuffle_pwes_identity_matrix_stays_identity():
    list_aas = [1, 2, 3]
    df_gauss = pd.DataFrame(
        np.array([[1, 0.5, 0.1], [0.5, 1, 0.5], [0.1, 0.5, 1]]), index=list_aas, columns=list_aas
    )
    df_pws = pd.DataFrame(np.eye(3))
    sorted_out, _ = pwes.shuffle_pwes(df_gauss, df_pws.copy(), list_aas, nrand=5)
    # permuting an identity matrix's rows/cols together always yields the identity matrix,
    # and the diagonal of df_gauss is always 1, so the average over shuffles stays identity.
    assert np.allclose(sorted_out.values, np.eye(3))


# --------------------------------------------------------------------------- #
# _safe_slug
# --------------------------------------------------------------------------- #

def test_safe_slug_replaces_separators_and_whitespace():
    assert pwes._safe_slug("some / weird\tname??") == "some_-_weird_name-"


def test_safe_slug_empty_string_becomes_EMPTY():
    assert pwes._safe_slug("   ") == "EMPTY"


def test_safe_slug_truncates_to_max_len():
    long_name = "A" * 200
    assert len(pwes._safe_slug(long_name, max_len=10)) == 10


# --------------------------------------------------------------------------- #
# make_cluster_color_dict
# --------------------------------------------------------------------------- #

def test_make_cluster_color_dict_sorted_and_rgb_tuples():
    out = pwes.make_cluster_color_dict([3, 1, 2])
    assert list(out.keys()) == [1, 2, 3]
    for rgb in out.values():
        assert len(rgb) == 3
        assert all(0.0 <= c <= 1.0 for c in rgb)


# --------------------------------------------------------------------------- #
# pymol_script
# --------------------------------------------------------------------------- #

def test_pymol_script_writes_expected_commands(tmp_path):
    residue_dict = {1: [168, 172], 2: ["O0010"]}
    colors = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    pwes.pymol_script(
        dir=str(tmp_path), file="script.pml", pdb_filename="/some/path/model.pdb",
        residue_dict=residue_dict, colors=colors,
    )
    text = (tmp_path / "script.pml").read_text()
    assert "load /some/path/model.pdb" in text
    assert "resi 168 and name CA or resi 172 and name CA" in text
    assert "resi 10 and chain O and name CA" in text
    assert "color 0xFF0000, Cluster1" in text
    assert "color 0x00FF00, Cluster2" in text


def test_pymol_script_skips_empty_residue_lists(tmp_path):
    pwes.pymol_script(
        dir=str(tmp_path), file="script.pml", pdb_filename="model.pdb",
        residue_dict={1: []}, colors=[(1.0, 0.0, 0.0)],
    )
    text = (tmp_path / "script.pml").read_text()
    assert "Cluster1" not in text


# --------------------------------------------------------------------------- #
# hist() / cat() - simple categorical plots
# --------------------------------------------------------------------------- #

def test_hist_runs_without_error():
    df_clus = pd.DataFrame({"cl_new": [1, 1, 2, 2, 2, 3]})
    pwes.hist(df_clus, show=False)


def test_hist_missing_cluster_col_raises():
    df_clus = pd.DataFrame({"other_col": [1, 2]})
    with pytest.raises(AssertionError):
        pwes.hist(df_clus, show=False)


def test_cat_runs_without_error():
    df_clus = pd.DataFrame({"cl_new": [1, 1, 2, 2], "log2(FC)": [0.5, -0.5, 1.0, -1.0]})
    pwes.cat(df_clus, show=False)


def test_cat_missing_scores_col_raises():
    df_clus = pd.DataFrame({"cl_new": [1, 2]})
    with pytest.raises(AssertionError):
        pwes.cat(df_clus, scores_col="not_a_col", show=False)


# --------------------------------------------------------------------------- #
# heatmap / heatmap_cbar / clustermap
# --------------------------------------------------------------------------- #

def test_heatmap_runs_without_error():
    labels = ["A1", "A2", "A3", "X1"]
    data = np.array([[1, 0.5, 0.2, 0], [0.5, 1, 0.3, 0], [0.2, 0.3, 1, 0], [0, 0, 0, 1]])
    df_scaled = pd.DataFrame(data, index=labels, columns=labels)
    pwes.heatmap(df_scaled, show=False)


def test_heatmap_cbar_runs_without_error():
    pwes.heatmap_cbar(show=False)
    pwes.heatmap_cbar(orientation="horizontal", show=False)


def test_clustermap_runs_without_error():
    labels = ["A1", "A2", "A3", "A4"]
    data = np.array(
        [[1, 0.9, 0.1, 0.05], [0.9, 1, 0.15, 0.1], [0.1, 0.15, 1, 0.8], [0.05, 0.1, 0.8, 1]]
    )
    df_scaled = pd.DataFrame(data, index=labels, columns=labels)
    link = linkage(data, method="ward")
    df_clusters = pd.DataFrame({"cl_new": [1, 1, 2, 2]}, index=labels)
    pwes.clustermap(df_scaled, link, df_clusters, show=False)


# --------------------------------------------------------------------------- #
# torn() - core tornado plot, including the secondary-structure/DSSP path
# --------------------------------------------------------------------------- #

def _write_mini_dssp(tmp_path, chain="A", n=10, ss_codes=None):
    """Build a minimal legacy-format DSSP file parseable by edms.data.dssp.parse_segments."""
    if ss_codes is None:
        ss_codes = ["H"] * 5 + [" "] * 5
    lines = [
        "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA \n"
    ]
    for i in range(n):
        resnum = i + 1
        line = list(" " * 80)
        for j, c in enumerate(str(resnum).rjust(5)):
            line[5 + j] = c
        line[11] = chain
        line[13] = "A"
        line[16] = ss_codes[i]
        lines.append("".join(line) + "\n")
    pt = tmp_path / "mini.dssp"
    pt.write_text("".join(lines))
    return str(pt)


@pytest.fixture()
def torn_inputs():
    df = pd.DataFrame({
        "Edit": ["A2C", "A3D", "A4E", "A5F", "A7G", "A8H", "A9I", "A10K"],
        "log2(FC)": [-2.0, -1.0, 0.5, 1.5, -0.8, 0.2, 1.1, -1.4],
        "weight": [1, 2, 3, 4, 5, 6, 7, 8],
    })
    df_clus = pd.DataFrame({
        "AA Number": [2, 3, 4, 5, 7, 8, 9, 10],
        "log2(FC)": [-2.0, -1.0, 0.5, 1.5, -0.8, 0.2, 1.1, -1.4],
        "cl_new": [1, 1, 1, 1, 2, 2, 2, 2],
        "Change": ["Polar", "Acidic", "Acidic", "Nonpolar", "Polar", "Basic", "Nonpolar", "Basic"],
        "weight": [1, 2, 3, 4, 5, 6, 7, 8],
    })
    return df, df_clus


def test_torn_basic_runs_and_returns_dataframe(torn_inputs):
    df, df_clus = torn_inputs
    out = pwes.torn(df=df, df_clus=df_clus, individual=True, show=False,
                     display_labels=False, size="weight", display_legend=True)
    assert isinstance(out, pd.DataFrame)
    assert "AA Number" in out.columns


def test_torn_individual_false_runs(torn_inputs):
    df, df_clus = torn_inputs
    # individual=False path re-derives 'AA Number'/'Change' from 'Edit' via add_label_info,
    # so df_clus only needs the raw 'Edit' column (not the precomputed 'AA Number'/'Change').
    df_clus_raw = df_clus.drop(columns=["AA Number", "Change"]).assign(Edit=df["Edit"].values)
    out = pwes.torn(df=df, df_clus=df_clus_raw, individual=False, show=False,
                     display_labels=False, size="weight", display_legend=True)
    assert isinstance(out, pd.DataFrame)


def test_torn_secondary_structure_dssp_path_regression(tmp_path, torn_inputs):
    """
    Regression test for the fixed bug at pwes.py ~line 858:
    `ss_y = min(df[f'log2({FC})'])-2*0.5` (undefined name `FC`) was fixed to
    `ss_y = min(df[scores_col])-2*0.5`. This exercises exactly that DSSP
    secondary-structure-track code path within torn(individual=True), which
    previously crashed with NameError before the fix.
    """
    df, df_clus = torn_inputs
    dssp_pt = _write_mini_dssp(tmp_path)

    out = pwes.torn(
        df=df, df_clus=df_clus, individual=True, DSSP=dssp_pt, chain_id="A",
        show=False, display_labels=False, size="weight", display_legend=True,
    )
    assert isinstance(out, pd.DataFrame)


def test_torn_secondary_structure_dssp_path_custom_scores_col(tmp_path):
    """Same regression as above, but with a non-default scores_col name to make sure
    the fix genuinely reads the parameter rather than a hardcoded column name."""
    df = pd.DataFrame({
        "Edit": ["A2C", "A3D", "A4E", "A7G", "A8H"],
        "my_score": [-2.0, -1.0, 0.5, -0.8, 0.2],
        "weight": [1, 2, 3, 4, 5],
    })
    df_clus = pd.DataFrame({
        "AA Number": [2, 3, 4, 7, 8],
        "my_score": [-2.0, -1.0, 0.5, -0.8, 0.2],
        "cl_new": [1, 1, 1, 2, 2],
        "Change": ["Polar", "Acidic", "Acidic", "Polar", "Basic"],
        "weight": [1, 2, 3, 4, 5],
    })
    dssp_pt = _write_mini_dssp(tmp_path)
    out = pwes.torn(
        df=df, df_clus=df_clus, scores_col="my_score", individual=True,
        DSSP=dssp_pt, chain_id="A", show=False, display_labels=False,
        size="weight", display_legend=True,
    )
    assert isinstance(out, pd.DataFrame)
    assert "my_score" in out.columns


# --------------------------------------------------------------------------- #
# torn() with size left at its default (None), or explicitly disabled
# (size=False), must not crash even though no size column is being used
# (regression test: size_norm now defaults to None up front, instead of only
# being assigned inside the `size is not None and size in df.columns` branch).
# --------------------------------------------------------------------------- #

def test_torn_default_size_param_runs_without_error(torn_inputs):
    df, df_clus = torn_inputs
    out = pwes.torn(df=df, df_clus=df_clus, individual=True, show=False, display_labels=False)
    assert isinstance(out, pd.DataFrame)


def test_torn_size_false_runs_without_error(torn_inputs):
    df, df_clus = torn_inputs
    out = pwes.torn(df=df, df_clus=df_clus, individual=True, show=False,
                     display_labels=False, size=False)
    assert isinstance(out, pd.DataFrame)
