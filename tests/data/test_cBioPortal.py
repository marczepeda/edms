"""Tests for edms.data.cBioPortal.

cBioPortal.mutations() does not perform any HTTP calls itself -- it processes
an already-downloaded cBioPortal mutation TSV/DataFrame (GENIE cohort export)
and translates cBioPortal protein-change notation into EDMS edit notation.
No network access is exercised anywhere in this module, so there is nothing
to mock at the HTTP layer; instead we focus on:
  - correct notation translation for every branch of mutations()
  - the WT-mismatch sanity checks (ValueError paths)
  - the str-path input branch (mocking io.get, not the network)
  - that config=True writes only within a redirected HOME (never real ~/.config)
"""
import os

import pandas as pd
import pytest

from edms.data import cBioPortal as cbp

# 1-indexed reference protein used across tests:
# pos:  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
# aa :  M  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
WT = "MACDEFGHIKLMNPQRSTVWY"


def _row(gene="GENE1", cancer_type="Lung", protein_change="K10Q", mutation_type="Missense_Mutation"):
    return {
        "Gene": gene,
        "Cancer Type": cancer_type,
        "Protein Change": protein_change,
        "Mutation Type": mutation_type,
    }


def test_missense_kept_as_is_and_counts_aggregated():
    df = pd.DataFrame(
        [
            _row(cancer_type="Lung", protein_change="K10Q"),
            _row(cancer_type="Lung", protein_change="K10Q"),
            _row(cancer_type="Breast", protein_change="K10Q"),
        ]
    )
    out = cbp.mutations(df=df, wt=WT, config=False)

    row = out[out["Protein Change"] == "K10Q"].iloc[0]
    assert row["Edit Change"] == "K10Q"
    assert row["counts"] == 3
    # Cancer Types ordered by descending frequency: Lung (2) before Breast (1)
    assert row["Cancer Types"] == "Lung, Breast"


def test_fusions_are_excluded():
    df = pd.DataFrame(
        [
            _row(protein_change="K10Q", mutation_type="Missense_Mutation"),
            _row(protein_change="FUSION1", mutation_type="fusion"),
        ]
    )
    out = cbp.mutations(df=df, wt=WT, config=False)
    assert "FUSION1" not in out["Protein Change"].values
    assert set(out["Mutation Type"]) == {"Missense_Mutation"}


def test_rows_missing_required_fields_are_dropped():
    df = pd.DataFrame(
        [
            _row(protein_change="K10Q"),
            {"Gene": "GENE1", "Cancer Type": "Lung", "Protein Change": None, "Mutation Type": "Missense_Mutation"},
        ]
    )
    out = cbp.mutations(df=df, wt=WT, config=False)
    assert len(out) == 1


@pytest.mark.parametrize(
    "protein_change,mutation_type,expected",
    [
        # In_Frame_Del: deletion + insertion (delins), 2 residue numbers
        ("K10_N12delinsQR", "In_Frame_Del", "KLMN10QRN"),
        # In_Frame_Del: single amino acid deletion
        ("K10del", "In_Frame_Del", "KL10L"),
        # In_Frame_Del: multiple amino acid deletion
        ("K10_N12del", "In_Frame_Del", "KLMN10N"),
        # In_Frame_Ins: deletion + insertion (delins), 1 residue number
        ("K10delinsQR", "In_Frame_Ins", "K10QR"),
        # In_Frame_Ins: single amino acid duplication
        ("K10dup", "In_Frame_Ins", "K10KK"),
        # In_Frame_Ins: multiple amino acid duplication
        ("K10_M12dup", "In_Frame_Ins", "M12MKLM"),
        # In_Frame_Ins: amino acid insertion between consecutive residues
        ("K10_L11insXY", "In_Frame_Ins", "K10KXY"),
    ],
)
def test_edit_notation_translation_branches(protein_change, mutation_type, expected):
    df = pd.DataFrame([_row(protein_change=protein_change, mutation_type=mutation_type)])
    out = cbp.mutations(df=df, wt=WT, config=False)
    assert out.iloc[0]["Edit Change"] == expected


def test_in_frame_ins_unknown_pattern_keeps_original_notation():
    # 2 numbers, contains neither 'delins', 'dup', nor 'ins' -> falls through to else-branch
    df = pd.DataFrame([_row(protein_change="K10_N12weird", mutation_type="In_Frame_Ins")])
    out = cbp.mutations(df=df, wt=WT, config=False)
    assert out.iloc[0]["Edit Change"] == "K10_N12weird"


@pytest.mark.parametrize(
    "protein_change,mutation_type",
    [
        ("Z10del", "In_Frame_Del"),  # single deletion: wt[9] is 'K', not 'Z'
        ("Z10_N12del", "In_Frame_Del"),  # multiple deletion
        ("Z10_N12delinsQR", "In_Frame_Del"),  # delins
        ("Z10delinsQR", "In_Frame_Ins"),  # ins delins
        ("Z10dup", "In_Frame_Ins"),  # single duplication
    ],
)
def test_wt_mismatch_sanity_check_raises_value_error(protein_change, mutation_type):
    df = pd.DataFrame([_row(protein_change=protein_change, mutation_type=mutation_type)])
    with pytest.raises(ValueError, match="mismatch with WT protein"):
        cbp.mutations(df=df, wt=WT, config=False)


def test_in_frame_ins_insertion_non_consecutive_positions_raises_value_error():
    # num2 (12) != num1(10) + 1 -> should raise the "non-consecutive" ValueError
    df = pd.DataFrame([_row(protein_change="K10_M12insXY", mutation_type="In_Frame_Ins")])
    with pytest.raises(ValueError, match="non-consecutive"):
        cbp.mutations(df=df, wt=WT, config=False)


def test_multiple_genes_produce_concatenated_dataframe():
    df = pd.DataFrame(
        [
            _row(gene="GENE1", protein_change="K10Q"),
            _row(gene="GENE2", protein_change="K10Q"),
        ]
    )
    out = cbp.mutations(df=df, wt=WT, config=False)
    assert set(out["Gene"]) == {"GENE1", "GENE2"}
    assert len(out) == 2


def test_string_path_input_loads_via_io_get(mocker):
    df = pd.DataFrame([_row(protein_change="K10Q")])
    mock_get = mocker.patch.object(cbp.io, "get", return_value=df)

    out = cbp.mutations(df="some/path/mutations.tsv", wt=WT, config=False)

    mock_get.assert_called_once_with("some/path/mutations.tsv")
    assert len(out) == 1


def test_dir_file_saves_output_csv_within_tmp_path(tmp_path):
    df = pd.DataFrame([_row(protein_change="K10Q")])
    out = cbp.mutations(df=df, wt=WT, config=False, dir=str(tmp_path), file="out.csv")

    out_file = tmp_path / "out.csv"
    assert out_file.exists()
    saved = pd.read_csv(out_file)
    assert saved.iloc[0]["Protein Change"] == "K10Q"


def test_config_true_writes_only_within_redirected_home(home_dir):
    """config=True (the default) saves a per-gene CSV to
    ~/.config/edms/cBioPortal_mutations/{gene}.csv. Verify it lands inside
    the redirected HOME (tmp_path) and never touches the real home dir.
    """
    df = pd.DataFrame([_row(gene="GENE1", protein_change="K10Q")])
    cbp.mutations(df=df, wt=WT, config=True)

    expected = home_dir / ".config" / "edms" / "cBioPortal_mutations" / "GENE1.csv"
    assert expected.exists()
