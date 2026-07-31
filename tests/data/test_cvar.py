"""Tests for edms.data.cvar (ClinVar).

Like cBioPortal.py and cosmic.py, this module performs no HTTP calls -- it
processes an already-downloaded ClinVar export DataFrame/TSV. There is
nothing to mock at the network layer.
"""
import pandas as pd
import pytest

from edms.data import cvar


# ---------------------------------------------------------------------------
# mutations() / prevalence()
# ---------------------------------------------------------------------------

def test_mutations_filters_gene_and_parses_position():
    df = pd.DataFrame(
        {
            "Gene(s)": ["GENE1", "GENE1", "OTHERGENE"],
            "Protein change": ["R123H", None, "G45D"],
        }
    )
    out = cvar.mutations(df=df, gene_name="GENE1")

    # Only GENE1 rows with non-null 'Protein change' survive
    assert len(out) == 1
    row = out.iloc[0]
    assert row["AA_position"] == 123
    assert row["AA_before"] == "R"
    assert row["AA_after"] == "H"


def test_mutations_saves_to_tmp_path(tmp_path):
    df = pd.DataFrame({"Gene(s)": ["GENE1"], "Protein change": ["R123H"]})
    cvar.mutations(df=df, gene_name="GENE1", dir=str(tmp_path), file="out.csv")
    assert (tmp_path / "out.csv").exists()


def test_mutations_accepts_string_path(mocker):
    df = pd.DataFrame({"Gene(s)": ["GENE1"], "Protein change": ["R123H"]})
    mock_get = mocker.patch.object(cvar.io, "get", return_value=df)

    out = cvar.mutations(df="some/path.tsv", gene_name="GENE1")

    mock_get.assert_called_once_with(pt="some/path.tsv")
    assert len(out) == 1


def test_prevalence_sorted_by_frequency():
    df = pd.DataFrame({"Protein change": ["R123H", "R123H", "G45D"]})
    assert cvar.prevalence(df) == ["R123H", "G45D"]


# ---------------------------------------------------------------------------
# priority_muts()
# ---------------------------------------------------------------------------

@pytest.fixture
def clinvar_priority_setup():
    df_clinvar = pd.DataFrame(
        {"Gene(s)": ["GENE1"] * 3, "Protein change": ["K10Q", "K10Q", "L20M"]}
    )
    pegRNAs_shared = pd.DataFrame(
        {
            "Target_name": ["GENE1", "GENE1"],
            "Edits": [["K10Q"], ["L20M"]],
            "Spacer_sequence": ["AAA", "CCC"],
            "PBS_sequence": ["TTT", "GGG"],
        }
    )
    pegRNAs = pd.DataFrame(
        {
            "Edit": ["K10Q", "L20M"],
            "Spacer_sequence": ["AAA", "CCC"],
            "PBS_sequence": ["TTT", "GGG"],
        }
    )
    return df_clinvar, pegRNAs_shared, pegRNAs


def test_priority_muts_picks_most_prevalent_available_edit(clinvar_priority_setup):
    df_clinvar, pegRNAs_shared, _ = clinvar_priority_setup
    out = cvar.priority_muts(pegRNAs_shared=pegRNAs_shared, df_clinvar=df_clinvar)
    assert list(out["Priority_mut"]) == ["K10Q", "L20M"]


def test_priority_muts_handles_string_edits_column():
    df_clinvar = pd.DataFrame({"Gene(s)": ["GENE1"] * 2, "Protein change": ["K10Q", "L20M"]})
    pegRNAs_shared = pd.DataFrame(
        {
            "Target_name": ["GENE1"],
            "Edits": ["['K10Q', 'L20M']"],
            "Spacer_sequence": ["AAA"],
            "PBS_sequence": ["TTT"],
        }
    )
    out = cvar.priority_muts(pegRNAs_shared=pegRNAs_shared, df_clinvar=df_clinvar)
    assert out.iloc[0]["Priority_mut"] in {"K10Q", "L20M"}


def test_priority_muts_reprocesses_raw_clinvar_dataframe_with_nans():
    # df_clinvar still containing NaNs in 'Protein change' -> mutations() re-run
    df_clinvar_raw = pd.DataFrame(
        {"Gene(s)": ["GENE1", "GENE1"], "Protein change": ["K10Q", None]}
    )
    pegRNAs_shared = pd.DataFrame(
        {
            "Target_name": ["GENE1"],
            "Edits": [["K10Q"]],
            "Spacer_sequence": ["AAA"],
            "PBS_sequence": ["TTT"],
        }
    )
    out = cvar.priority_muts(pegRNAs_shared=pegRNAs_shared, df_clinvar=df_clinvar_raw)
    assert out.iloc[0]["Priority_mut"] == "K10Q"


# ---------------------------------------------------------------------------
# priority_edits()
# ---------------------------------------------------------------------------

def test_priority_edits_counts_clinvar_occurrences(clinvar_priority_setup):
    df_clinvar, pegRNAs_shared, pegRNAs = clinvar_priority_setup
    shared_out = cvar.priority_muts(pegRNAs_shared=pegRNAs_shared, df_clinvar=df_clinvar)

    out = cvar.priority_edits(pegRNAs=pegRNAs, pegRNAs_shared=shared_out, df_clinvar=df_clinvar)

    counts = dict(zip(out["Edit"], out["ClinVar_count"]))
    assert counts == {"K10Q": 2, "L20M": 1}
