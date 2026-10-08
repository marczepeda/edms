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
    qPCR.ddCq(data=two_sample_two_target_df, file=str(out_dir / "ddcq.csv"))
    assert (out_dir / "ddcq.csv").is_file()


# --------------------------------------------------------------------------- #
# cfx_amp() & amp()
# --------------------------------------------------------------------------- #
@pytest.fixture
def amp_csv(tmp_path):
    # CFX amplification export: unnamed index column, Cycle, then one column per well
    cycles = np.arange(1, 11)
    sig = lambda mid: 1000 / (1 + np.exp(-(cycles - mid)))
    csv_pt = tmp_path / "amp.csv"
    pd.DataFrame({"Unnamed": "", "Cycle": cycles, "A1": sig(5), "A2": sig(7), "B1": np.zeros(10)}
                 ).rename(columns={"Unnamed": ""}).to_csv(csv_pt, index=False)
    return csv_pt


@pytest.fixture
def cq_csv(tmp_path):
    csv_pt = tmp_path / "cq.csv"
    pd.DataFrame({"Well": ["A01", "A02", "B01"], "Target": ["GeneA", "GeneA", np.nan],
                  "Sample": [np.nan] * 3, "Cq": [5.0, 7.0, np.nan]}).to_csv(csv_pt)
    return csv_pt


def test_cfx_amp_tidy_format(amp_csv):
    df = qPCR.cfx_amp(pt=str(amp_csv))
    assert {"Well", "Row", "Column", "Cycle", "RFU"} <= set(df.columns)
    assert not any(str(c).startswith("Unnamed") for c in df.columns)
    assert len(df) == 30
    assert set(df["Row"]) == {"A", "B"} and set(df["Column"]) == {1, 2}


def test_cfx_amp_annot_merges_padded_wells_and_drops_empty_cols(amp_csv, cq_csv):
    df = qPCR.cfx_amp(pt=str(amp_csv), annot=str(cq_csv))
    assert "Sample" not in df.columns  # all-NaN column dropped
    assert df.loc[df["Well"] == "A2", "Cq"].iloc[0] == 7.0  # A02 -> A2


def test_cfx_amp_filters_wells_and_drop_no_Cq(amp_csv, cq_csv):
    assert set(qPCR.cfx_amp(pt=str(amp_csv), wells=["A01", "B1"])["Well"]) == {"A1", "B1"}
    assert set(qPCR.cfx_amp(pt=str(amp_csv), wells_exclude=["A1"])["Well"]) == {"A2", "B1"}
    assert set(qPCR.cfx_amp(pt=str(amp_csv), annot=str(cq_csv), drop_no_Cq=True)["Well"]) == {"A1", "A2"}
    with pytest.raises(ValueError):
        qPCR.cfx_amp(pt=str(amp_csv), drop_no_Cq=True)


def test_cfx_amp_baseline_subtraction(amp_csv):
    df = qPCR.cfx_amp(pt=str(amp_csv), baseline=(1, 3))
    base = df[(df["Well"] == "A2") & (df["Cycle"] <= 3)]["RFU"]
    assert base.mean() == pytest.approx(0.0)


def test_amp_one_line_per_well_and_saves(amp_csv, tmp_path):
    out = tmp_path / "out" / "amp.png"
    fig, axes = qPCR.amp(df=str(amp_csv), threshold=100, file=str(out), dpi=50, show=False)
    assert out.is_file()
    ax = axes.flat[0]
    assert len([l for l in ax.lines if l.get_linestyle() == "-"]) == 3  # wells not aggregated
    assert any(l.get_linestyle() == "--" for l in ax.lines)  # threshold


def test_find_cfx_single_multiple_missing(tmp_path):
    assert qPCR.find_cfx("Quantification Cq Results", dir=str(tmp_path), required=False) is None
    with pytest.raises(FileNotFoundError):
        qPCR.find_cfx("Quantification Cq Results", dir=str(tmp_path))
    (tmp_path / "run -  Quantification Cq Results_0.csv").write_text("Well\n")
    assert qPCR.find_cfx("Quantification Cq Results", dir=str(tmp_path)).endswith("Results_0.csv")
    (tmp_path / "run2 -  Quantification Cq Results_0.csv").write_text("Well\n")
    with pytest.raises(ValueError):
        qPCR.find_cfx("Quantification Cq Results", dir=str(tmp_path))


def test_amp_auto_detects_amp_and_cq_files(tmp_path, monkeypatch, amp_csv, cq_csv):
    amp_csv.rename(tmp_path / "run -  Quantification Amplification Results_SYBR.csv")
    cq_csv.rename(tmp_path / "run -  Quantification Cq Results_0.csv")
    monkeypatch.chdir(tmp_path)
    fig, axes = qPCR.amp(cols="Target", drop_no_Cq=True, show=False)  # Cq & Target come from auto-detected annot
    assert len(axes.flat[0].lines) >= 2


def test_ddCq_auto_detects_cq_file(tmp_path, monkeypatch, two_sample_two_target_df):
    two_sample_two_target_df.to_csv(tmp_path / "run -  Quantification Cq Results_0.csv", index=False)
    monkeypatch.chdir(tmp_path)
    out = qPCR.ddCq()
    assert np.allclose(out["RQ_mean"], 1.0)
