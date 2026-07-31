"""
Tests for edms.gen.tidy

Pure dataframe/dict/string transform utilities - tests use hand-computed
expected values on small inputs rather than "doesn't crash" checks.
"""
import pandas as pd
import pytest

from edms.gen import tidy as t


# --------------------------------------------------------------------------- #
# Dataframe methods
# --------------------------------------------------------------------------- #

def test_reorder_cols_keep_true():
    df = pd.DataFrame({"a": [1], "b": [2], "c": [3]})
    out = t.reorder_cols(df, cols=["c", "a"], keep=True)
    assert list(out.columns) == ["c", "a", "b"]


def test_reorder_cols_keep_false_drops_others():
    df = pd.DataFrame({"a": [1], "b": [2], "c": [3]})
    out = t.reorder_cols(df, cols=["b"], keep=False)
    assert list(out.columns) == ["b"]


def test_zip_cols():
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    result = list(t.zip_cols(df, ["a", "b"]))
    assert result == [(1, 3), (2, 4)]


def test_missing_cols():
    df = pd.DataFrame({"a": [1], "b": [2]})
    assert t.missing_cols(df, ["a", "b", "c", "d"]) == ["c", "d"]
    assert t.missing_cols(df, ["a", "b"]) == []


def test_merge_with_string_id():
    data = pd.DataFrame({"id": [1, 2, 3]})
    meta = pd.DataFrame({"id": [1, 2, 3], "name": ["x", "y", "z"], "val": [10, 20, 30]})
    out = t.merge(data, meta, id="id", cols=["name", "val"])
    assert list(out["name"]) == ["x", "y", "z"]
    assert list(out["val"]) == [10, 20, 30]


def test_merge_with_list_id():
    data = pd.DataFrame({"data_id": [1, 2]})
    meta = pd.DataFrame({"meta_id": [1, 2], "name": ["x", "y"]})
    out = t.merge(data, meta, id=["data_id", "meta_id"], cols=["name"])
    assert list(out["name"]) == ["x", "y"]


def test_merge_invalid_id_prints_error(capsys):
    data = pd.DataFrame({"id": [1]})
    meta = pd.DataFrame({"id": [1], "name": ["x"]})
    out = t.merge(data, meta, id=123, cols=["name"])
    captured = capsys.readouterr()
    assert "Error" in captured.out
    # dataframe returned unchanged (no 'name' column added)
    assert "name" not in out.columns


def test_shared_group():
    df = pd.DataFrame({
        "shared": ["g1", "g1", "g2"],
        "other": ["a", "b", "c"],
    })
    out = t.shared_group(df, shared="shared", group=["other"])
    # each row gets a list of unique 'other' values sharing its 'shared' value
    row0 = out[out["other"] == "a"].iloc[0]
    assert set(row0["other_list"]) == {"a", "b"}
    row2 = out[out["other"] == "c"].iloc[0]
    assert row2["other_list"] == ["c"]


def test_unique_tuples_str_col():
    df = pd.DataFrame({"a": [1, 2, 1, 3]})
    assert t.unique_tuples(df, "a") == [(1,), (2,), (3,)]


def test_unique_tuples_multi_col_order_preserved():
    df = pd.DataFrame({"a": [1, 2, 1], "b": ["x", "y", "x"]})
    assert t.unique_tuples(df, ["a", "b"]) == [(1, "x"), (2, "y")]


def test_vcs_ordered():
    df = pd.DataFrame({"a": [1, 1, 2, 3, 3, 3]})
    vcs = t.vcs_ordered(df, "a")
    assert list(vcs.index) == [(1,), (2,), (3,)]
    assert list(vcs.values) == [2, 1, 3]


def test_unpack_singleton_iterables():
    df = pd.DataFrame({
        "a": [[1], (2,), {3}, [1, 2], [], "text"],
    })
    out = t.unpack_singleton_iterables(df)
    assert out["a"].tolist() == [1, 2, 3, [1, 2], [], "text"]


def test_unpack_singleton_iterables_not_inplace_by_default():
    df = pd.DataFrame({"a": [[1]]})
    out = t.unpack_singleton_iterables(df, inplace=False)
    assert out["a"].tolist() == [1]
    # original untouched
    assert df["a"].tolist() == [[1]]


def test_drop_duplicates_n_keep_first_2():
    df = pd.DataFrame({"a": [1, 1, 1, 2, 2]})
    out = t.drop_duplicates_n(df, subset="a", n=2)
    assert out["a"].tolist() == [1, 1, 2, 2]


def test_drop_duplicates_n_default_n1():
    df = pd.DataFrame({"a": [1, 1, 2]})
    out = t.drop_duplicates_n(df)
    assert out["a"].tolist() == [1, 2]


# --------------------------------------------------------------------------- #
# Dictionary methods
# --------------------------------------------------------------------------- #

def test_filter_kwargs():
    out = t.filter_kwargs(["a", "b"], a=1, b=None, c=3)
    # b is filtered out because its value is None
    assert out == {"a": 1}


def test_comb_dcs():
    ls = [{"a": 1, "b": 2}, {"a": 3, "b": 4}]
    out = t.comb_dcs(ls)
    assert out == {"a": [1, 3], "b": [2, 4]}


# --------------------------------------------------------------------------- #
# Methods for dictionary containing dataframes
# --------------------------------------------------------------------------- #

def test_split_by():
    series = ["a, b", "c", 5, "d, e, f"]
    out = t.split_by(series, by=", ")
    assert out == ["a", "b", "c", "d", "e", "f"]


def test_isolate_want_exact_list():
    dc = {"k1": pd.DataFrame({"col": ["a", "b", "c"]})}
    out = t.isolate(dc, col="col", get=["a", "c"], want=True, exact=True)
    assert out["k1"]["col"].tolist() == ["a", "c"]


def test_isolate_not_want_exact_list():
    dc = {"k1": pd.DataFrame({"col": ["a", "b", "c"]})}
    out = t.isolate(dc, col="col", get=["a"], want=False, exact=True)
    assert out["k1"]["col"].tolist() == ["b", "c"]


def test_isolate_scalar_get():
    dc = {"k1": pd.DataFrame({"col": ["a", "b", "a"]})}
    out = t.isolate(dc, col="col", get="a", want=True)
    assert out["k1"]["col"].tolist() == ["a", "a"]


def test_isolate_get_none_want_true_returns_null_rows():
    dc = {"k1": pd.DataFrame({"col": ["a", None, "c"]})}
    out = t.isolate(dc, col="col", get=None, want=True)
    assert out["k1"].shape[0] == 1
    assert out["k1"]["col"].isnull().all()


def test_modify_with_static_value():
    dc = {"k1": pd.DataFrame({"a": [1, 2]})}
    out = t.modify(dc, col="b", val=99)
    assert out["k1"]["b"].tolist() == [99, 99]
    # original untouched
    assert "b" not in dc["k1"].columns


def test_modify_with_callable():
    dc = {"k1": pd.DataFrame({"a": [1, 2, 3]})}
    out = t.modify(dc, col="b", val=lambda row: row["a"] * 2, axis=1)
    assert out["k1"]["b"].tolist() == [2, 4, 6]


def test_melt_with_id_vars_only():
    dc = {"k1": pd.DataFrame({"id": [1, 2], "x": [10, 20], "y": [30, 40]})}
    out = t.melt(dc, id_vars=["id"])
    df = out["k1"]
    assert set(df.columns) == {"id", "variable", "value"}
    assert sorted(df["variable"].unique().tolist()) == ["x", "y"]


def test_melt_missing_id_vars_raises():
    dc = {"k1": pd.DataFrame({"x": [1]})}
    with pytest.raises(TypeError):
        t.melt(dc, id_vars=["missing"])


def test_melt_no_vars_raises():
    dc = {"k1": pd.DataFrame({"x": [1]})}
    with pytest.raises(TypeError):
        t.melt(dc)


def test_join():
    dc = {"g1": pd.DataFrame({"a": [1, 2]}), "g2": pd.DataFrame({"a": [3]})}
    out = t.join(dc, col="key")
    assert out.shape[0] == 3
    assert out[out["key"] == "g1"]["a"].tolist() == [1, 2]
    assert out[out["key"] == "g2"]["a"].tolist() == [3]


def test_split():
    df = pd.DataFrame({"grp": ["a", "b", "a"], "val": [1, 2, 3]})
    out = t.split(df, key="grp")
    assert set(out.keys()) == {"a", "b"}
    assert out["a"]["val"].tolist() == [1, 3]
    assert out["b"]["val"].tolist() == [2]


def test_join_split_roundtrip():
    dc = {"a": pd.DataFrame({"val": [1, 3]}), "b": pd.DataFrame({"val": [2]})}
    joined = t.join(dc, col="grp")
    split_back = t.split(joined, key="grp")
    assert split_back["a"]["val"].tolist() == [1, 3]
    assert split_back["b"]["val"].tolist() == [2]


# --------------------------------------------------------------------------- #
# Dictionary <-> list interconversion
# --------------------------------------------------------------------------- #

def test_dc_to_ls_simple():
    dc = {"a": 1, "b": 2}
    out = t.dc_to_ls(dc, sep=".")
    assert set(out) == {"a.1", "b.2"}


def test_dc_to_ls_nested():
    dc = {"a": {"b": 1}}
    out = t.dc_to_ls(dc, sep=".")
    assert out == ["a.b.1"]


def test_ls_to_dc_simple():
    out = t.ls_to_dc(["a.1", "b.2"], sep=".")
    assert out == {"a": "1", "b": "2"}


def test_ls_to_dc_nested():
    out = t.ls_to_dc(["a.b.1"], sep=".")
    assert out == {"a": {"b": "1"}}


def test_dc_to_ls_ls_to_dc_roundtrip_nested():
    dc = {"a": {"b": "1", "c": "2"}, "d": "3"}
    ls = t.dc_to_ls(dc, sep=".")
    back = t.ls_to_dc(ls, sep=".")
    assert back == dc


# --------------------------------------------------------------------------- #
# String methods
# --------------------------------------------------------------------------- #

def test_find_all():
    assert t.find_all("abcabcabc", "abc") == [0, 3, 6]
    assert t.find_all("abc", "xyz") == []


def test_find_all_overlapping():
    assert t.find_all("aaaa", "aa") == [0, 1, 2]


def test_split_nth():
    assert t.split_nth("a-b-c-d-e", "-", 2) == ["a-b", "c-d", "e"]


def test_split_nth_n_larger_than_pieces():
    assert t.split_nth("a-b", "-", 5) == ["a-b"]


def test_natural_key_sorts_numbers_in_human_order():
    names = ["file10.txt", "file2.txt", "file1.txt"]
    out = sorted(names, key=t.natural_key)
    assert out == ["file1.txt", "file2.txt", "file10.txt"]


def test_natural_key_case_insensitive():
    assert t.natural_key("ABC") == t.natural_key("abc")
