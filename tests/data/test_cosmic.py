"""Tests for edms.data.cosmic.

Like cBioPortal.py, this module never performs any HTTP calls: it processes
an already-downloaded COSMIC export DataFrame/TSV. There is nothing to mock
at the network layer. We mock only the plotting entry points
(edms.gen.plot.dist / edms.gen.plot.stack) so tests stay fast and hermetic,
and use tmp_path for any on-disk output.
"""
import pandas as pd
import pytest

from edms.data import cosmic


# ---------------------------------------------------------------------------
# mutations() / prevalence()
# ---------------------------------------------------------------------------

def test_mutations_parses_aa_position_before_after():
    df = pd.DataFrame(
        {
            "AA Mutation": ["p.R123H", "p.R123H", "p.G45_A46del", "p.Unknown0"],
            "Type": [
                "Substitution - Missense",
                "Substitution - Missense",
                "Deletion - In frame",
                "Unknown",
            ],
        }
    )
    out = cosmic.mutations(df=df)

    # 'Unknown' type row dropped entirely
    assert "Unknown" not in out["Type"].values
    assert len(out) == 3

    row = out[out["AA Mutation"] == "p.R123H"].iloc[0]
    assert row["AA_position"] == 123
    assert row["AA_before"] == "R"
    assert row["AA_after"] == "H"
    assert row["AA_mut"] == "R123H"


def test_mutations_accepts_string_path(mocker):
    df = pd.DataFrame({"AA Mutation": ["p.R123H"], "Type": ["Substitution - Missense"]})
    mock_get = mocker.patch.object(cosmic.io, "get", return_value=df)

    out = cosmic.mutations(df="some/path.tsv")

    mock_get.assert_called_once_with(pt="some/path.tsv")
    assert len(out) == 1


def test_mutations_saves_to_tmp_path(tmp_path):
    df = pd.DataFrame({"AA Mutation": ["p.R123H"], "Type": ["Substitution - Missense"]})
    cosmic.mutations(df=df, dir=str(tmp_path), file="out.csv")
    assert (tmp_path / "out.csv").exists()


def test_prevalence_sorted_by_frequency():
    df = pd.DataFrame({"AA_mut": ["R123H", "R123H", "G45D", "R123H", "G45D"]})
    assert cosmic.prevalence(df) == ["R123H", "G45D"]


# ---------------------------------------------------------------------------
# cds_group()
# ---------------------------------------------------------------------------

def test_cds_group_assigns_region_and_calls_plot_dist(mocker):
    mock_dist = mocker.patch("edms.data.cosmic.p.dist")

    df_cosmic = pd.DataFrame(
        {
            "AA_position": [5, 15, 25, 35, 45, 55, 65],
            "Type": [
                "Substitution - Missense",
                "Substitution - Nonsense",
                "Insertion - Frameshift",
                "Insertion - In frame",
                "Deletion - Frameshift",
                "Deletion - In frame",
                "Substitution - coding silent",  # dropped by cds_group
            ],
            "AA_mut": ["A5B", "C15D", "E25F", "G35H", "I45J", "K55L", "M65N"],
        }
    )
    df_cds = pd.DataFrame(
        {"CDS": ["exon1", "exon2"], "start": [1, 31], "end": [30, 70], "gene": ["GENE1", "GENE1"]}
    )

    cosmic.cds_group(df_cosmic=df_cosmic, df_cds=df_cds, out_dir=None)

    mock_dist.assert_called_once()
    kwargs = mock_dist.call_args.kwargs
    # Silent mutation (position 65) must have been dropped before grouping
    assert 65 not in kwargs["df"]["AA_position"].values
    regions = dict(zip(kwargs["df"]["AA_position"], kwargs["df"]["CDS"]))
    assert regions == {5: "exon1", 15: "exon1", 25: "exon1", 35: "exon2", 45: "exon2", 55: "exon2"}
    assert kwargs["bins"] == 70
    assert kwargs["y_axis"] == "GENE1"
    assert kwargs["title"] == "COSMIC Mutations"


def test_cds_group_position_outside_all_regions_gets_zero(mocker):
    mock_dist = mocker.patch("edms.data.cosmic.p.dist")

    df_cosmic = pd.DataFrame(
        {"AA_position": [1000], "Type": ["Substitution - Missense"], "AA_mut": ["Z1000W"]}
    )
    df_cds = pd.DataFrame({"CDS": ["exon1"], "start": [1], "end": [30], "gene": ["GENE1"]})

    cosmic.cds_group(df_cosmic=df_cosmic, df_cds=df_cds, out_dir=None)

    kwargs = mock_dist.call_args.kwargs
    assert kwargs["df"].iloc[0]["CDS"] == 0


def test_cds_group_accepts_string_paths(mocker):
    mocker.patch("edms.data.cosmic.p.dist")
    df_cosmic = pd.DataFrame(
        {"AA_position": [5], "Type": ["Substitution - Missense"], "AA_mut": ["A5B"]}
    )
    df_cds = pd.DataFrame({"CDS": ["exon1"], "start": [1], "end": [30], "gene": ["GENE1"]})
    mock_get = mocker.patch.object(cosmic.io, "get", side_effect=[df_cosmic, df_cds])

    cosmic.cds_group(df_cosmic="cosmic.csv", df_cds="cds.csv", out_dir=None)

    assert mock_get.call_count == 2


# ---------------------------------------------------------------------------
# priority_muts() / priority_edits() (prime editing helpers)
# ---------------------------------------------------------------------------

@pytest.fixture
def cosmic_priority_setup():
    df_cosmic = pd.DataFrame(
        {"AA_mut": ["K10Q", "K10Q", "L20M"], "Type": ["Substitution - Missense"] * 3}
    )
    pegRNAs_shared = pd.DataFrame(
        {
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
    return df_cosmic, pegRNAs_shared, pegRNAs


def test_priority_muts_picks_most_prevalent_available_edit(cosmic_priority_setup):
    df_cosmic, pegRNAs_shared, _ = cosmic_priority_setup
    out = cosmic.priority_muts(pegRNAs_shared=pegRNAs_shared, df_cosmic=df_cosmic)
    assert list(out["Priority_mut"]) == ["K10Q", "L20M"]


def test_priority_muts_handles_string_edits_column():
    df_cosmic = pd.DataFrame({"AA_mut": ["K10Q", "L20M", "L20M"], "Type": ["Substitution - Missense"] * 3})
    pegRNAs_shared = pd.DataFrame(
        {"Edits": ["['K10Q', 'L20M']"], "Spacer_sequence": ["AAA"], "PBS_sequence": ["TTT"]}
    )
    out = cosmic.priority_muts(pegRNAs_shared=pegRNAs_shared, df_cosmic=df_cosmic)
    # L20M (count 2) is more prevalent than K10Q (count 1) -> chosen first
    assert out.iloc[0]["Priority_mut"] == "L20M"


def test_priority_muts_reprocesses_raw_cosmic_dataframe():
    # df_cosmic without AA_mut column -> should be run through mutations() first
    df_cosmic_raw = pd.DataFrame(
        {"AA Mutation": ["p.K10Q"], "Type": ["Substitution - Missense"]}
    )
    pegRNAs_shared = pd.DataFrame(
        {"Edits": [["K10Q"]], "Spacer_sequence": ["AAA"], "PBS_sequence": ["TTT"]}
    )
    out = cosmic.priority_muts(pegRNAs_shared=pegRNAs_shared, df_cosmic=df_cosmic_raw)
    assert out.iloc[0]["Priority_mut"] == "K10Q"


def test_priority_edits_counts_cosmic_occurrences(cosmic_priority_setup):
    df_cosmic, pegRNAs_shared, pegRNAs = cosmic_priority_setup
    shared_out = cosmic.priority_muts(pegRNAs_shared=pegRNAs_shared, df_cosmic=df_cosmic)
    out = cosmic.priority_edits(pegRNAs=pegRNAs, pegRNAs_shared=shared_out, df_cosmic=df_cosmic)

    counts = dict(zip(out["Edit"], out["COSMIC_count"]))
    assert counts == {"K10Q": 2, "L20M": 1}


# ---------------------------------------------------------------------------
# editor_mutations()
# ---------------------------------------------------------------------------

def test_editor_mutations_classifies_and_saves(tmp_path, mocker):
    mock_stack = mocker.patch("edms.data.cosmic.p.stack")

    df_bescan = pd.DataFrame(
        {
            "gene": ["GENE1"],
            "AtoG_mutations": ["A5B;C15D"],
            "CtoT_mutations": ["E25F/G35H"],
        }
    )
    df_cosmic = pd.DataFrame(
        {
            "AA_mut": [
                "A5B", "C15D", "E25F", "G35H", "I45J", "K55L", "M65N", "P75Q", "R85S",
            ],
            "Type": [
                "Substitution - Missense",  # A5B -> BE
                "Substitution - Nonsense",  # C15D -> BE
                "Substitution - Missense",  # E25F -> BE
                "Substitution - Nonsense",  # G35H -> BE
                "Insertion - Frameshift",   # PE indel
                "Deletion - In frame",      # PE indel
                "Substitution - coding silent",  # complex/other
                "Insertion - In frame",     # PE indel
                "Deletion - Frameshift",    # PE indel
            ],
        }
    )

    cosmic.editor_mutations(df_cosmic=df_cosmic, df_bescan=df_bescan, out_dir=str(tmp_path))

    assert (tmp_path / "GENE1_BE_COSMIC.csv").exists()
    assert (tmp_path / "GENE1_COSMIC_types.csv").exists()
    editor_csv = tmp_path / "GENE1_editor_COSMIC.csv"
    assert editor_csv.exists()

    df_editor = pd.read_csv(editor_csv)
    counts = dict(zip(df_editor["mutation"], df_editor["COSMIC"]))
    assert counts["Substitution (ABE/CBE)"] == 4  # A5B, C15D, E25F, G35H
    assert counts["Indel (PE)"] == 4  # I45J, K55L, P75Q, R85S
    assert counts["Substitution (PE)"] == 0
    assert counts["Complex"] == 0
    mock_stack.assert_called_once()


def test_editor_mutations_no_out_dir_skips_saving(mocker):
    mocker.patch("edms.data.cosmic.p.stack")
    mock_save = mocker.patch.object(cosmic.io, "save")

    df_bescan = pd.DataFrame({"gene": ["GENE1"], "AtoG_mutations": [""], "CtoT_mutations": [""]})
    df_cosmic = pd.DataFrame(
        {
            "AA_mut": ["A5B", "C15D", "I45J", "K55L", "P75Q", "R85S"],
            "Type": [
                "Substitution - Missense",
                "Substitution - Nonsense",
                "Insertion - Frameshift",
                "Deletion - In frame",
                "Insertion - In frame",
                "Deletion - Frameshift",
            ],
        }
    )

    cosmic.editor_mutations(df_cosmic=df_cosmic, df_bescan=df_bescan, out_dir=None)

    mock_save.assert_not_called()
