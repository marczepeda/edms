"""
Tests for edms.gen.io

I/O functions - all reads/writes are confined to tmp_path fixtures; nothing
touches the real filesystem outside pytest's temp area.
"""
import os

import pandas as pd
import pytest

from edms.gen import io


# --------------------------------------------------------------------------- #
# Parsing Python literals
# --------------------------------------------------------------------------- #

def test_df_try_parse_converts_literal_strings():
    df = pd.DataFrame({"a": ["[1, 2, 3]", "hello"], "b": ["{'x': 1}", "5"]})
    out = io.df_try_parse(df)
    assert out["a"].iloc[0] == [1, 2, 3]
    assert out["a"].iloc[1] == "hello"
    assert out["b"].iloc[0] == {"x": 1}
    assert out["b"].iloc[1] == 5


def test_recursive_json_decode_nested():
    obj = {"a": '{"b": [1, 2, "x"]}', "c": "not json"}
    out = io.recursive_json_decode(obj)
    # "x" is not valid JSON so it stays a string; the outer JSON string is
    # decoded into a real dict/list structure.
    assert out == {"a": {"b": [1, 2, "x"]}, "c": "not json"}


def test_recursive_json_decode_list_of_dicts():
    obj = ['{"x": 1}', '{"y": 2}']
    out = io.recursive_json_decode(obj)
    assert out == [{"x": 1}, {"y": 2}]


# --------------------------------------------------------------------------- #
# get() / get_dir()
# --------------------------------------------------------------------------- #

def test_get_csv_round_trip(tmp_path):
    df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
    pt = tmp_path / "data.csv"
    df.to_csv(pt, index=False)
    out = io.get(str(pt))
    pd.testing.assert_frame_equal(out, df)


def test_get_tsv_round_trip(tmp_path):
    df = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
    pt = tmp_path / "data.tsv"
    df.to_csv(pt, index=False, sep="\t")
    out = io.get(str(pt))
    pd.testing.assert_frame_equal(out, df)


def test_get_csv_literal_eval(tmp_path):
    df = pd.DataFrame({"a": ["[1, 2]", "[3, 4]"]})
    pt = tmp_path / "data.csv"
    df.to_csv(pt, index=False)
    out = io.get(str(pt), literal_eval=True)
    assert out["a"].iloc[0] == [1, 2]
    assert out["a"].iloc[1] == [3, 4]


def test_get_xlsx_round_trip(tmp_path):
    df1 = pd.DataFrame({"a": [1, 2]})
    df2 = pd.DataFrame({"b": [3, 4]})
    pt = tmp_path / "data.xlsx"
    with pd.ExcelWriter(pt) as writer:
        df1.to_excel(writer, sheet_name="sheet1", index=False)
        df2.to_excel(writer, sheet_name="sheet2", index=False)
    out = io.get(str(pt))
    assert set(out.keys()) == {"sheet1", "sheet2"}
    pd.testing.assert_frame_equal(out["sheet1"], df1)
    pd.testing.assert_frame_equal(out["sheet2"], df2)


def test_get_unknown_suffix_falls_back_to_csv(tmp_path):
    df = pd.DataFrame({"a": [1, 2]})
    pt = tmp_path / "data.txt"
    df.to_csv(pt, index=False)
    out = io.get(str(pt))
    pd.testing.assert_frame_equal(out, df)


def test_get_dir(tmp_path):
    pd.DataFrame({"a": [1]}).to_csv(tmp_path / "one.csv", index=False)
    pd.DataFrame({"b": [2]}).to_csv(tmp_path / "two.csv", index=False)
    (tmp_path / "ignore.txt").write_text("not a csv")
    dc = io.get_dir(str(tmp_path))
    assert set(dc.keys()) == {"one", "two"}
    assert dc["one"]["a"].iloc[0] == 1
    assert dc["two"]["b"].iloc[0] == 2


# --------------------------------------------------------------------------- #
# save() / save_dir()
# --------------------------------------------------------------------------- #

def test_save_dataframe_csv(tmp_path):
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    io.save(obj=df, file="out.csv", dir=str(tmp_path))
    out = pd.read_csv(tmp_path / "out.csv")
    pd.testing.assert_frame_equal(out, df)


def test_save_dataframe_isolate_cols(tmp_path):
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4], "c": [5, 6]})
    io.save(obj=df, file="out.csv", dir=str(tmp_path), cols=["a", "c"])
    out = pd.read_csv(tmp_path / "out.csv")
    assert list(out.columns) == ["a", "c"]


def test_save_dataframe_invalid_cols_type_raises(tmp_path):
    df = pd.DataFrame({"a": [1]})
    with pytest.raises(ValueError):
        io.save(obj=df, file="out.csv", dir=str(tmp_path), cols=[1, 2])


def test_save_dataframe_tsv(tmp_path):
    df = pd.DataFrame({"a": [1], "b": [2]})
    io.save(obj=df, file="out.tsv", dir=str(tmp_path))
    out = pd.read_csv(tmp_path / "out.tsv", sep="\t")
    pd.testing.assert_frame_equal(out, df)


def test_save_dataframe_xlsx(tmp_path):
    df = pd.DataFrame({"a": [1], "b": [2]})
    io.save(obj=df, file="out.xlsx", dir=str(tmp_path))
    out = pd.read_excel(tmp_path / "out.xlsx")
    pd.testing.assert_frame_equal(out, df)


def test_save_dataframe_with_index(tmp_path):
    df = pd.DataFrame({"a": [1, 2]}, index=["r1", "r2"])
    io.save(obj=df, file="out.csv", dir=str(tmp_path), id=True)
    out = pd.read_csv(tmp_path / "out.csv", index_col=0)
    assert list(out.index) == ["r1", "r2"]


def test_save_list_sorted(tmp_path):
    io.save(obj=[3, 1, 2], file="out.csv", dir=str(tmp_path))
    content = (tmp_path / "out.csv").read_text().strip()
    assert content == "1,2,3"


def test_save_list_unsorted(tmp_path):
    io.save(obj=[3, 1, 2], file="out.csv", dir=str(tmp_path), sort=False)
    content = (tmp_path / "out.csv").read_text().strip()
    assert content == "3,1,2"


def test_save_set(tmp_path):
    io.save(obj={2, 1, 3}, file="out.csv", dir=str(tmp_path))
    content = (tmp_path / "out.csv").read_text().strip()
    assert content == "1,2,3"


def test_save_dict_to_xlsx(tmp_path):
    dc = {"a": pd.DataFrame({"x": [1]}), "b": pd.DataFrame({"x": [2]})}
    io.save(obj=dc, file="out.xlsx", dir=str(tmp_path))
    result = pd.read_excel(tmp_path / "out.xlsx", sheet_name=None)
    assert set(result.keys()) == {"a", "b"}
    assert result["a"]["x"].iloc[0] == 1
    assert result["b"]["x"].iloc[0] == 2


def test_save_unsupported_type_raises(tmp_path):
    with pytest.raises(ValueError):
        io.save(obj=42, file="out.csv", dir=str(tmp_path))


def test_save_dir(tmp_path):
    dc = {"one": pd.DataFrame({"a": [1]}), "two": pd.DataFrame({"a": [2]})}
    io.save_dir(dir=str(tmp_path), dc=dc)
    assert (tmp_path / "one.csv").exists()
    assert (tmp_path / "two.csv").exists()


def test_save_file_none_returns_without_error(tmp_path, capsys):
    # file=None -> check_outpath returns (None, None); save() should warn
    # and return without raising.
    io.save(obj=pd.DataFrame({"a": [1]}), file=None, dir=str(tmp_path))
    captured = capsys.readouterr()
    assert "not saved" in captured.out


# --------------------------------------------------------------------------- #
# excel_csvs()
# --------------------------------------------------------------------------- #

def test_excel_csvs(tmp_path):
    df1 = pd.DataFrame({"a": [1, 2]})
    df2 = pd.DataFrame({"b": [3, 4]})
    xlsx_pt = tmp_path / "book.xlsx"
    with pd.ExcelWriter(xlsx_pt) as writer:
        df1.to_excel(writer, sheet_name="sheet1", index=False)
        df2.to_excel(writer, sheet_name="sheet2", index=False)

    out_dir = tmp_path / "csvs"
    io.excel_csvs(str(xlsx_pt), dir=str(out_dir))

    out1 = pd.read_csv(out_dir / "sheet1.csv")
    out2 = pd.read_csv(out_dir / "sheet2.csv")
    pd.testing.assert_frame_equal(out1, df1)
    pd.testing.assert_frame_equal(out2, df2)


# --------------------------------------------------------------------------- #
# df_to_dc_txt() / dc_txt_to_df()
# --------------------------------------------------------------------------- #

def test_df_to_dc_txt_and_back_roundtrip():
    df = pd.DataFrame({"col1": [1, 2], "col2": ["a", "b"]})
    txt = io.df_to_dc_txt(df)
    assert "col1" in txt and "col2" in txt
    # dc_txt_to_df's default transpose=True is what makes the round trip
    # symmetric: df_to_dc_txt() emits {index: {col: val, ...}, ...}, and
    # transposing after ast.literal_eval puts columns back on the columns
    # axis and the original index back on the index axis.
    back = io.dc_txt_to_df(txt)
    assert list(back["col1"]) == [1, 2]
    assert list(back["col2"]) == ["a", "b"]


# --------------------------------------------------------------------------- #
# in_subs() / out_subs()
# --------------------------------------------------------------------------- #

def test_in_subs_basename(tmp_path):
    (tmp_path / "sample1.txt").write_text("a")
    (tmp_path / "sample2.txt").write_text("b")
    io.in_subs(str(tmp_path), suf=".txt")
    assert (tmp_path / "sample1" / "sample1.txt").exists()
    assert (tmp_path / "sample2" / "sample2.txt").exists()


def test_in_subs_prefix_sep(tmp_path):
    (tmp_path / "sampleA_R1.txt").write_text("a")
    (tmp_path / "sampleA_R2.txt").write_text("b")
    io.in_subs(str(tmp_path), suf=".txt", group_by="prefix", prefix_sep="_")
    assert (tmp_path / "sampleA" / "sampleA_R1.txt").exists()
    assert (tmp_path / "sampleA" / "sampleA_R2.txt").exists()


def test_in_subs_prefix_len(tmp_path):
    (tmp_path / "abcdef.txt").write_text("a")
    io.in_subs(str(tmp_path), suf=".txt", group_by="prefix", prefix_len=3)
    assert (tmp_path / "abc" / "abcdef.txt").exists()


def test_in_subs_prefix_missing_sep_and_len_raises(tmp_path):
    (tmp_path / "a.txt").write_text("x")
    with pytest.raises(ValueError):
        io.in_subs(str(tmp_path), suf=".txt", group_by="prefix", prefix_sep=None, prefix_len=None)


def test_in_subs_invalid_dir_raises():
    with pytest.raises(ValueError):
        io.in_subs("/nonexistent/dir/path/xyz", suf=".txt")


def test_out_subs(tmp_path):
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "file1.txt").write_text("a")
    io.out_subs(str(tmp_path))
    assert (tmp_path / "file1.txt").exists()
    assert not sub.exists()


def test_out_subs_resolves_name_conflicts(tmp_path):
    (tmp_path / "file1.txt").write_text("top-level")
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "file1.txt").write_text("nested")
    io.out_subs(str(tmp_path))
    assert (tmp_path / "file1.txt").read_text() == "top-level"
    assert (tmp_path / "file1_1.txt").read_text() == "nested"


def test_out_subs_invalid_dir_raises():
    with pytest.raises(ValueError):
        io.out_subs("/nonexistent/dir/path/xyz")


# --------------------------------------------------------------------------- #
# create_sh()
# --------------------------------------------------------------------------- #

def test_create_sh_writes_expected_content(tmp_path):
    io.create_sh(dir=str(tmp_path), file="job.sh", cores=2, mem=2000, email="test@example.com")
    content = (tmp_path / "job.sh").read_text()
    assert "#!/bin/bash" in content
    assert "-n 2" in content
    assert "--mem=2000" in content
    assert "test@example.com" in content


def test_create_sh_requires_sh_extension(tmp_path):
    with pytest.raises(ValueError):
        io.create_sh(dir=str(tmp_path), file="job.txt", email="test@example.com")


# --------------------------------------------------------------------------- #
# combine()
# --------------------------------------------------------------------------- #

def test_combine_concatenates_with_headers(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    (in_dir / "b.txt").write_text("second\n")
    (in_dir / "a.txt").write_text("first\n")
    out_dir = tmp_path / "out"

    out_path = io.combine(in_dir=str(in_dir), out_dir=str(out_dir), out_file="combined.txt")
    content = out_path.read_text()

    assert "### SOURCE_FILE: a.txt" in content
    assert "### SOURCE_FILE: b.txt" in content
    # natural-sorted: a.txt before b.txt
    assert content.index("a.txt") < content.index("b.txt")
    assert "first" in content and "second" in content


def test_combine_no_matching_files_raises(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    (in_dir / "irrelevant.csv").write_text("x")
    with pytest.raises(FileNotFoundError):
        io.combine(in_dir=str(in_dir), out_dir=str(tmp_path / "out"), out_file="combined.txt")


def test_combine_empty_suffixes_raises(tmp_path):
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    (in_dir / "a.txt").write_text("x")
    with pytest.raises(ValueError):
        io.combine(in_dir=str(in_dir), out_dir=str(tmp_path / "out"), out_file="combined.txt", suffixes=[])


def test_combine_recursive(tmp_path):
    in_dir = tmp_path / "in"
    (in_dir / "sub").mkdir(parents=True)
    (in_dir / "top.txt").write_text("top\n")
    (in_dir / "sub" / "nested.txt").write_text("nested\n")
    out_path = io.combine(in_dir=str(in_dir), out_dir=str(tmp_path / "out"), out_file="combined.txt", recursive=True)
    content = out_path.read_text()
    assert "top" in content and "nested" in content


# --------------------------------------------------------------------------- #
# split_R1_R2()
# --------------------------------------------------------------------------- #

def test_split_r1_r2(tmp_path):
    (tmp_path / "sample_R1_001.fastq").write_text("r1")
    (tmp_path / "sample_R2_001.fastq").write_text("r2")
    (tmp_path / "other.fastq").write_text("other")
    io.split_R1_R2(str(tmp_path))
    assert (tmp_path / "R1" / "sample_R1_001.fastq").exists()
    assert (tmp_path / "R2" / "sample_R2_001.fastq").exists()
    # non-matching file untouched
    assert (tmp_path / "other.fastq").exists()


# --------------------------------------------------------------------------- #
# Directory methods
# --------------------------------------------------------------------------- #

def test_relative_paths(tmp_path):
    (tmp_path / "a.txt").write_text("x")
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "b.txt").write_text("y")
    paths = io.relative_paths(str(tmp_path))
    assert set(paths) == {"a.txt", os.path.join("sub", "b.txt")}


def test_sorted_file_names(tmp_path):
    (tmp_path / "b.csv").write_text("x")
    (tmp_path / "a.csv").write_text("y")
    (tmp_path / "c.txt").write_text("z")
    names = io.sorted_file_names(str(tmp_path), suf=".csv")
    assert names == ["a.csv", "b.csv"]
