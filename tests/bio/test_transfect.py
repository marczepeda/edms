"""
Tests for edms.bio.transfect

Covers:
- PE3(): epegRNA/ngRNA/PE plasmid merge, uL calculations (hand-computed against the
  documented defaults), tube-B master-mix scale-up, and 96-well plate placement.
- virus(): VSVG/GagPol/transfer categorization, ng & uL calculations, and the
  Optimem/P3000/L3000 bookend rows.

All numeric expectations below were derived from the documented formulas in
transfect.py and cross-checked by running the functions directly (see PR description);
they are not just re-derivations of the source code.
"""
import pandas as pd
import pytest

from edms.bio import transfect


# --------------------------------------------------------------------------- #
# PE3()
# --------------------------------------------------------------------------- #
@pytest.fixture
def plasmids_df():
    return pd.DataFrame({
        "Plasmid": ["epeg1", "ng1", "pMUZ86.7"],
        "Description": ["epeg desc", "ng desc", "PE desc"],
        "Colony": ["c1", "c2", "c3"],
        "ng/uL": [100.0, 200.0, 300.0],
    })


@pytest.fixture
def epegRNAs_df():
    return pd.DataFrame({"pegRNA_number": [1], "Name": ["epeg1"]})


@pytest.fixture
def ngRNAs_df():
    return pd.DataFrame({"pegRNA_number": [1], "Name": ["ng1"]})


def test_PE3_default_uL_calculations(plasmids_df, epegRNAs_df, ngRNAs_df):
    pivots = transfect.PE3(plasmids=plasmids_df, epegRNAs=epegRNAs_df, ngRNAs=ngRNAs_df)
    tube_A = pivots["Transfection"].iloc[0]

    # mm_x=1.1, reps=3: epegRNA_ng=66/100 ng_uL -> 2.178 -> round 2.18
    assert tube_A["epegRNA uL"] == pytest.approx(2.18)
    # ngRNA_ng=22/200 ng_uL -> 0.363 -> round 0.36
    assert tube_A["ngRNA uL"] == pytest.approx(0.36)
    # PE_ng=200/300 ng_uL -> 2.2
    assert tube_A["PE uL"] == pytest.approx(2.2)
    # Total = mm_x*reps*well_uL/2 = 1.1*3*10/2 = 16.5
    assert tube_A["Total uL"] == pytest.approx(16.5)
    # Optimem = Total - epeg - ng - PE
    assert tube_A["Optimem uL"] == pytest.approx(16.5 - 2.18 - 0.36 - 2.2)
    assert tube_A["Tube"] == "A1"


def test_PE3_merges_epeg_ngRNA_and_PE_names(plasmids_df, epegRNAs_df, ngRNAs_df):
    pivots = transfect.PE3(plasmids=plasmids_df, epegRNAs=epegRNAs_df, ngRNAs=ngRNAs_df)
    tube_A = pivots["Transfection"].iloc[0]
    assert tube_A["epegRNA"] == "epeg1"
    assert tube_A["ngRNA"] == "ng1"
    assert tube_A["PE"] == "pMUZ86.7"
    assert tube_A["PE Colony"] == "c3"
    assert tube_A["PE ng/uL"] == 300.0


def test_PE3_tube_B_scales_with_number_of_conditions(plasmids_df, epegRNAs_df, ngRNAs_df):
    pivots = transfect.PE3(plasmids=plasmids_df, epegRNAs=epegRNAs_df, ngRNAs=ngRNAs_df)
    transfection = pivots["Transfection"]
    tube_B = transfection[transfection["Tube"] == "B"].iloc[0]
    n_conditions = 1  # one epegRNA/ngRNA pair in the fixture
    mm_x, reps, well_uL = 1.1, 3, 10
    assert tube_B["Optimem uL"] == pytest.approx(round(9 / 10 * mm_x * reps * mm_x * n_conditions, 2))
    assert tube_B["L2000 uL"] == pytest.approx(round(1 / 10 * well_uL / 2 * mm_x * reps * mm_x * n_conditions, 2))


def test_PE3_places_conditions_on_96_well_plate(plasmids_df, epegRNAs_df, ngRNAs_df):
    pivots = transfect.PE3(plasmids=plasmids_df, epegRNAs=epegRNAs_df, ngRNAs=ngRNAs_df)
    plate_pivot = pivots["96-well plates"]
    # Placement starts at row B, column 2 by design (row/col indices start at 1,
    # deliberately skipping the outermost row/column of the 96-well plate -- see
    # the "excluding outer wells" comment in transfect.PE3()).
    assert plate_pivot.loc[(1, "B"), 2] == "A1"


def test_PE3_multiple_ngRNAs_for_same_pegRNA(plasmids_df, epegRNAs_df):
    ngRNAs = pd.DataFrame({"pegRNA_number": [1, 1], "Name": ["ng1", "ng2"]})
    plasmids = pd.concat([plasmids_df, pd.DataFrame({
        "Plasmid": ["ng2"], "Description": ["ng2 desc"], "Colony": ["c4"], "ng/uL": [150.0],
    })], ignore_index=True)
    pivots = transfect.PE3(plasmids=plasmids, epegRNAs=epegRNAs_df, ngRNAs=ngRNAs)
    assert len(pivots["Transfection"]) == 3  # 2 tube A rows + 1 tube B row
    assert set(pivots["Transfection"]["ngRNA"].dropna()) == {"ng1", "ng2"}


# --------------------------------------------------------------------------- #
# virus()
# --------------------------------------------------------------------------- #
@pytest.fixture
def virus_plasmids_df():
    return pd.DataFrame({
        "Plasmid": ["pMUZ26.6", "pMUZ26.7", "transfer1.2"],
        "Description": ["VSVG", "GagPol", "Transfer"],
        "Colony": ["c1", "c2", "c3"],
        "ng/uL": [100.0, 200.0, 50.0],
    })


def test_virus_default_calculations(virus_plasmids_df):
    df_virus = transfect.virus(plasmids=virus_plasmids_df.copy())

    optimem_row = df_virus[(df_virus["Description"] == "Optimem") & (df_virus["Tube"] == "A")].iloc[0]
    p3000_row = df_virus[df_virus["Description"] == "P3000 Enhancer"].iloc[0]
    vsvg_row = df_virus[df_virus["Plasmid"] == "pMUZ26.6"].iloc[0]
    gagpol_row = df_virus[df_virus["Plasmid"] == "pMUZ26.7"].iloc[0]
    transfer_row = df_virus[df_virus["Plasmid"] == "transfer1.2"].iloc[0]
    l3000_row = df_virus[df_virus["Description"] == "L3000"].iloc[0]

    # transfers_num = 1 (only the non-VSVG/GagPol plasmid)
    assert optimem_row["uL"] == pytest.approx(275.0)  # 1*1.1*1*500/2
    assert p3000_row["uL"] == pytest.approx(6.6)       # 1*1.1*1*6
    assert l3000_row["uL"] == pytest.approx(7.7)       # 1*1.1*1*7

    assert vsvg_row["ng"] == pytest.approx(825.0)      # 1*1.1*1*750
    assert vsvg_row["uL"] == pytest.approx(8.25)       # 825/100
    assert gagpol_row["ng"] == pytest.approx(1650.0)   # 1*1.1*1*1500
    assert gagpol_row["uL"] == pytest.approx(8.25)     # 1650/200

    assert transfer_row["ng"] == pytest.approx(750.0)  # 750*1
    assert transfer_row["uL"] == pytest.approx(15.0)   # 750/50
    assert transfer_row["Tube"] == "1.2"               # digits parsed from plasmid name


def test_virus_row_order_places_A_before_transfer_before_B(virus_plasmids_df):
    df_virus = transfect.virus(plasmids=virus_plasmids_df.copy())
    orders = df_virus["Order"].tolist()
    assert orders == sorted(orders)
    assert df_virus.iloc[0]["Tube"] == "A"
    assert df_virus.iloc[-1]["Tube"] == "B"


def test_virus_multiple_transfer_plasmids_scale_transfers_num():
    plasmids = pd.DataFrame({
        "Plasmid": ["pMUZ26.6", "pMUZ26.7", "transfer1.1", "transfer1.2"],
        "Description": ["VSVG", "GagPol", "Transfer1", "Transfer2"],
        "Colony": ["c1", "c2", "c3", "c4"],
        "ng/uL": [100.0, 200.0, 50.0, 50.0],
    })
    df_virus = transfect.virus(plasmids=plasmids.copy())
    optimem_row = df_virus[(df_virus["Description"] == "Optimem") & (df_virus["Tube"] == "A")].iloc[0]
    # transfers_num = 2 now.
    assert optimem_row["uL"] == pytest.approx(2 * 1.1 * 1 * 500 / 2)
