"""
Tests for edms.bio.plate

Covers:
- Helpers: _is_blank_row, _looks_like_plate_header, _normalize_cell, _infer_plate_size
- parse_csv_to_tidy(): single & multi-block CSVs, well-label construction, replicate
  counting (both replicate_requires_value=True/False), default plate naming, and the
  "no header row detected" error path.
- make(): pivoting tidy data back into a plate-shaped DataFrame, single & combined
  value columns.
"""
import os

import pandas as pd
import pytest

from edms.bio import plate


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def test_is_blank_row_true_for_all_empty():
    assert plate._is_blank_row(["", " ", ""]) is True


def test_is_blank_row_false_when_any_nonempty():
    assert plate._is_blank_row(["", "x", ""]) is False


def test_looks_like_plate_header_true_for_numeric_row():
    assert plate._looks_like_plate_header(["96-well 1", "1", "2", "3"]) is True


def test_looks_like_plate_header_true_for_well_in_first_cell():
    assert plate._looks_like_plate_header(["Well ID", "a", "b"]) is True


def test_looks_like_plate_header_false_for_empty_first_cell():
    assert plate._looks_like_plate_header(["", "1", "2"]) is False


def test_looks_like_plate_header_false_for_data_row():
    assert plate._looks_like_plate_header(["A", "s1", "s2"]) is False


def test_normalize_cell_strips_whitespace():
    assert plate._normalize_cell("  hi  ") == "hi"
    assert plate._normalize_cell(None) == ""


@pytest.mark.parametrize("rows,cols,expected", [
    (2, 3, 6),
    (3, 4, 12),
    (4, 6, 24),
    (6, 8, 48),
    (8, 12, 96),
    (5, 5, None),
])
def test_infer_plate_size(rows, cols, expected):
    assert plate._infer_plate_size(rows, cols) == expected


# --------------------------------------------------------------------------- #
# parse_csv_to_tidy()
# --------------------------------------------------------------------------- #
def _write_csv(tmp_path, text, name="plate.csv"):
    pt = tmp_path / name
    pt.write_text(text)
    return pt


def test_parse_csv_to_tidy_single_block(tmp_path):
    csv_text = "96-well 1,1,2\nA,s1,s2\nB,s3,s1\n"
    pt = _write_csv(tmp_path, csv_text)

    df = plate.parse_csv_to_tidy(pt)

    assert list(df.columns) == ["plate", "col", "row", "well", "value", "replicate"]
    assert len(df) == 4
    # Well labels are row label + numeric column label.
    wells = set(df["well"])
    assert wells == {"A1", "A2", "B1", "B2"}
    # 's1' appears twice across the block -> replicate 1 then 2 (stable row/col order).
    s1_rows = df[df["value"] == "s1"].sort_values("well")
    assert list(s1_rows["replicate"]) == [1, 2]


def test_parse_csv_to_tidy_skips_empty_cells_by_default(tmp_path):
    csv_text = "96-well 1,1,2\nA,s1,\nB,,s2\n"
    pt = _write_csv(tmp_path, csv_text)
    df = plate.parse_csv_to_tidy(pt)
    assert len(df) == 2
    assert set(df["value"]) == {"s1", "s2"}


def test_parse_csv_to_tidy_multi_block(tmp_path):
    csv_text = (
        "96-well 1,1,2\n"
        "A,s1,s2\n"
        "\n"
        "96-well 2,1,2\n"
        "A,s3,s4\n"
    )
    pt = _write_csv(tmp_path, csv_text)
    df = plate.parse_csv_to_tidy(pt)
    assert set(df["plate"]) == {"96-well 1", "96-well 2"}
    assert len(df) == 4


def test_parse_csv_to_tidy_no_header_raises(tmp_path):
    csv_text = "just,some,data\nwithout,a,header\n"
    pt = _write_csv(tmp_path, csv_text)
    with pytest.raises(ValueError, match="No plate header rows detected"):
        plate.parse_csv_to_tidy(pt)


def test_parse_csv_to_tidy_replicate_requires_value_false_counts_blanks(tmp_path):
    csv_text = "96-well 1,1,2\nA,s1,\nB,,s1\n"
    pt = _write_csv(tmp_path, csv_text, name="p2.csv")
    df = plate.parse_csv_to_tidy(pt, skip_empty_cells=False, replicate_requires_value=False)
    # groupby('value').cumcount()+1 groups blank cells together as one 'value' too.
    s1_rows = df[df["value"] == "s1"]
    assert sorted(s1_rows["replicate"]) == [1, 2]


def test_parse_csv_to_tidy_saves_file(tmp_path):
    csv_text = "96-well 1,1,2\nA,s1,s2\n"
    pt = _write_csv(tmp_path, csv_text)
    out_dir = tmp_path / "out"
    plate.parse_csv_to_tidy(pt, dir=str(out_dir), file="tidy.csv")
    assert os.path.isfile(out_dir / "tidy.csv")


# --------------------------------------------------------------------------- #
# make()
# --------------------------------------------------------------------------- #
def test_make_pivots_tidy_dataframe_to_plate_shape():
    df = pd.DataFrame({
        "plate": ["p1", "p1", "p1", "p1"],
        "row": ["A", "A", "B", "B"],
        "col": [1, 2, 1, 2],
        "value": ["s1", "s2", "s3", "s4"],
    })
    pivot = plate.make(df=df, values="value")
    assert pivot.loc[("p1", "A"), 1] == "s1"
    assert pivot.loc[("p1", "A"), 2] == "s2"
    assert pivot.loc[("p1", "B"), 1] == "s3"
    assert pivot.loc[("p1", "B"), 2] == "s4"


def test_make_combines_multiple_value_columns():
    df = pd.DataFrame({
        "plate": ["p1", "p1"],
        "row": ["A", "A"],
        "col": [1, 2],
        "gene": ["gA", "gB"],
        "rep": ["1", "2"],
    })
    pivot = plate.make(df=df, values=["gene", "rep"])
    assert pivot.loc[("p1", "A"), 1] == "gA_1"
    assert pivot.loc[("p1", "A"), 2] == "gB_2"


def test_make_return_df_false_returns_none():
    df = pd.DataFrame({
        "plate": ["p1"],
        "row": ["A"],
        "col": [1],
        "value": ["s1"],
    })
    result = plate.make(df=df, values="value", return_df=False)
    assert result is None


def test_make_from_csv_path(tmp_path):
    csv_pt = tmp_path / "tidy.csv"
    pd.DataFrame({
        "plate": ["p1", "p1"],
        "row": ["A", "A"],
        "col": [1, 2],
        "value": ["s1", "s2"],
    }).to_csv(csv_pt, index=False)
    pivot = plate.make(df=str(csv_pt), values="value")
    assert pivot.loc[("p1", "A"), 1] == "s1"
