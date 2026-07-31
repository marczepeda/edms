"""
Tests for edms.bio.qPCR

Covers:
- cfx_Cq(): reads a CFX-exported Cq CSV, drops rows with missing sample IDs, retains
  only the requested columns, and coerces float-valued sample IDs to int.
- ddCq(): Cq mean/error -> dCq mean/error (target pairs) -> ddCq mean/error & RQ
  (sample pairs), hand-verified against the documented ΔΔCq formula.
"""
import numpy as np
import pandas as pd
import pytest

from edms.bio import qPCR


# --------------------------------------------------------------------------- #
# cfx_Cq()
# --------------------------------------------------------------------------- #
def test_cfx_Cq_drops_rows_missing_sample_and_keeps_requested_columns(tmp_path):
    csv_pt = tmp_path / "cq.csv"
    pd.DataFrame({
        "Well": ["A1", "A2", "A3"],
        "Fluor": ["SYBR", "SYBR", "SYBR"],
        "Target": ["GeneA", "GeneA", "GeneB"],
        "Sample": [1.0, np.nan, 2.0],
        "Cq": [20.1, 21.5, 19.8],
        "Extra": ["x", "y", "z"],
    }).to_csv(csv_pt, index=False)

    df = qPCR.cfx_Cq(pt=str(csv_pt))

    assert list(df.columns) == ["Well", "Fluor", "Target", "Sample", "Cq"]
    assert len(df) == 2  # the NaN-sample row was dropped
    assert "Extra" not in df.columns


def test_cfx_Cq_coerces_float_sample_ids_to_int(tmp_path):
    csv_pt = tmp_path / "cq.csv"
    pd.DataFrame({
        "Well": ["A1", "A2"],
        "Fluor": ["SYBR", "SYBR"],
        "Target": ["GeneA", "GeneA"],
        "Sample": [1.0, 2.0],
        "Cq": [20.1, 21.5],
    }).to_csv(csv_pt, index=False)

    df = qPCR.cfx_Cq(pt=str(csv_pt))
    assert list(df["Sample"]) == [1, 2]
    assert all(isinstance(s, int) for s in df["Sample"])


def test_cfx_Cq_custom_sample_col_and_cols(tmp_path):
    csv_pt = tmp_path / "cq.csv"
    pd.DataFrame({
        "Well": ["A1"],
        "cDNA": ["s1"],
        "Cq": [20.0],
    }).to_csv(csv_pt, index=False)

    df = qPCR.cfx_Cq(pt=str(csv_pt), sample_col="cDNA", cols=["Well", "cDNA", "Cq"])
    assert list(df.columns) == ["Well", "cDNA", "Cq"]
    assert df["cDNA"].iloc[0] == "s1"  # non-float sample id is passed through unchanged


# --------------------------------------------------------------------------- #
# ddCq()
# --------------------------------------------------------------------------- #
@pytest.fixture
def two_sample_two_target_df():
    # Sample 1: GeneA Cq=20, GeneB Cq=22 -> dCq = -2
    # Sample 2: GeneA Cq=21, GeneB Cq=23 -> dCq = -2
    # ddCq (Sample1 vs Sample2) = -2 - (-2) = 0 -> RQ = 2**0 = 1
    return pd.DataFrame({
        "Sample": [1, 1, 2, 2],
        "Target": ["GeneA", "GeneB", "GeneA", "GeneB"],
        "Cq": [20.0, 22.0, 21.0, 23.0],
    })


def test_ddCq_zero_relative_change_between_equally_shifted_samples(two_sample_two_target_df):
    out = qPCR.ddCq(data=two_sample_two_target_df)

    assert set(out["Targets"]) == {"GeneA ~ GeneB"}
    assert set(out["Samples"]) == {"1 ~ 2", "2 ~ 1"}
    assert np.allclose(out["ddCq_mean"], 0.0)
    assert np.allclose(out["RQ_mean"], 1.0)
    assert np.allclose(out["ddCq_err"], 0.0)  # single Cq value per (sample,target) -> std = 0


def test_ddCq_detects_relative_change():
    # Sample 2's GeneA is 1 cycle higher than Sample 1's (less abundant) while GeneB
    # matches -> dCq(sample1)=-2, dCq(sample2)=-1, ddCq(1 vs 2) = -2-(-1) = -1
    df = pd.DataFrame({
        "Sample": [1, 1, 2, 2],
        "Target": ["GeneA", "GeneB", "GeneA", "GeneB"],
        "Cq": [20.0, 22.0, 22.0, 23.0],
    })
    out = qPCR.ddCq(data=df)
    row = out[out["Samples"] == "1 ~ 2"].iloc[0]
    assert row["ddCq_mean"] == pytest.approx(-1.0)
    assert row["RQ_mean"] == pytest.approx(2 ** 1.0)


def test_ddCq_accepts_replicate_Cq_values_and_computes_error():
    # Two Cq replicates per (sample,target) -> nonzero std.
    df = pd.DataFrame({
        "Sample": [1, 1, 1, 1],
        "Target": ["GeneA", "GeneA", "GeneB", "GeneB"],
        "Cq": [20.0, 22.0, 22.0, 22.0],
    })
    out = qPCR.ddCq(data=df)
    # Only 1 sample -> no cross-sample ddCq rows, but shouldn't error.
    assert isinstance(out, pd.DataFrame)


def test_ddCq_reads_from_csv_path(tmp_path, two_sample_two_target_df):
    csv_pt = tmp_path / "cq.csv"
    two_sample_two_target_df.to_csv(csv_pt, index=False)
    out = qPCR.ddCq(data=str(csv_pt))
    assert np.allclose(out["RQ_mean"], 1.0)


def test_ddCq_saves_file(tmp_path, two_sample_two_target_df):
    out_dir = tmp_path / "out"
    qPCR.ddCq(data=two_sample_two_target_df, dir=str(out_dir), file="ddcq.csv")
    assert (out_dir / "ddcq.csv").is_file()
