"""
Tests for edms.bio.sanger

Covers pure-Python / pandas logic:
- group_boundaries(): consecutive-integer run grouping
- min_sec(): decimal-minutes -> (minutes, seconds)
- thermocycler(): per-primer-pair thermocycler step tables (n must be the int 1, and
  the amplicon-length column is 'PCR1 bp', unlike ngs.py's PCR2-based n=1/n=2 setup)
- pcr_mm() / pcr_mm_ultra(): NEB Q5 / NEBNext Ultra II master-mix uL calculations

Note: group_boundaries(), min_sec(), pcr_mm(), and pcr_mm_ultra() are byte-for-byte
duplicates of the ngs.py versions (this module predates a ngs/sanger split), so the
numeric expectations here mirror tests/bio/test_ngs.py by design.
"""
import pandas as pd
import pytest

from edms.bio import sanger


# --------------------------------------------------------------------------- #
# group_boundaries()
# --------------------------------------------------------------------------- #
def test_group_boundaries_empty():
    assert sanger.group_boundaries([]) == []


def test_group_boundaries_multiple_runs():
    assert sanger.group_boundaries([1, 2, 3, 5, 6, 9]) == [(1, 3), (5, 6), (9, 9)]


# --------------------------------------------------------------------------- #
# min_sec()
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("decimal_minutes,expected", [
    (2.5, (2, 30)),
    (0.25, (0, 15)),
    (3.0, (3, 0)),
])
def test_min_sec(decimal_minutes, expected):
    assert sanger.min_sec(decimal_minutes) == expected


# --------------------------------------------------------------------------- #
# thermocycler()
# --------------------------------------------------------------------------- #
@pytest.fixture
def pcr1_df():
    return pd.DataFrame({
        "ID": ["s1", "s2", "s3"],
        "PCR1 FWD": ["F1", "F1", "F1"],
        "PCR1 REV": ["R1", "R1", "R1"],
        "PCR1 ID": ["p1", "p2", "p3"],
        "PCR1 Tm": [65, 65, 65],
        "PCR1 bp": [300, 300, 300],
    })


def test_thermocycler_step_table_structure(pcr1_df):
    dc = sanger.thermocycler(df=pcr1_df, n=1, cycles=30)
    assert list(dc.keys()) == ["F1_R1_65°C"]

    table = dc["F1_R1_65°C"]
    assert list(table["Temperature"]) == ["98°C", "98°C", "65°C", "72°C", "72°C", "4°C", ""]
    assert list(table["Repeat"]) == ["", "30 cycles", "30 cycles", "30 cycles", "", "", ""]


def test_thermocycler_anneal_time_from_amplicon_length(pcr1_df):
    # bp=300 -> floor(300/500)/2+0.5 = 0.5 -> min_sec(0.5) = (0,30) -> "30s"
    dc = sanger.thermocycler(df=pcr1_df, n=1, cycles=30)
    assert dc["F1_R1_65°C"]["Time"].iloc[3] == "30s"


def test_thermocycler_title_ranges_consecutive_ids(pcr1_df):
    dc = sanger.thermocycler(df=pcr1_df, n=1, cycles=30)
    assert dc["F1_R1_65°C"].index[-1] == "F1_R1: p1 -> p3"


def test_thermocycler_requires_PCR1_bp_column():
    # sanger.thermocycler() always reads 'PCR1 bp' (unlike ngs.thermocycler() which
    # always reads 'PCR2 bp' regardless of n); a dataframe without it should raise.
    df = pd.DataFrame({
        "ID": ["s1"], "PCR1 FWD": ["F1"], "PCR1 REV": ["R1"],
        "PCR1 ID": ["p1"], "PCR1 Tm": [65],
    })
    with pytest.raises(KeyError):
        sanger.thermocycler(df=df, n=1)


def test_thermocycler_invalid_n_raises(pcr1_df):
    with pytest.raises(ValueError, match="n must be 1"):
        sanger.thermocycler(df=pcr1_df, n=2)


# --------------------------------------------------------------------------- #
# pcr_mm() / pcr_mm_ultra()
# --------------------------------------------------------------------------- #
def test_pcr_mm_default_uL_calculations():
    primers = pd.Series([2], index=pd.MultiIndex.from_tuples([("F1", "R1")]))
    mm = sanger.pcr_mm(primers=primers, template="gDNA", template_uL=5)
    table = mm[("F1", "R1")]

    def uL_for(component):
        return table.loc[table["Component"] == component, "uL"].iloc[0]

    assert uL_for("5x Q5 Reaction Buffer") == pytest.approx(5.0)
    assert uL_for("dNTPs") == pytest.approx(0.5)
    assert uL_for("F1") == pytest.approx(1.25)
    assert uL_for("gDNA") == pytest.approx(5.0)
    assert uL_for("Total") == pytest.approx(25.0)


def test_pcr_mm_ultra_default_uL_calculations():
    primers = pd.Series([2], index=pd.MultiIndex.from_tuples([("F1", "R1")]))
    mm = sanger.pcr_mm_ultra(primers=primers, template="gDNA", template_uL=5)
    table = mm[("F1", "R1")]

    def uL_for(component):
        return table.loc[table["Component"] == component, "uL"].iloc[0]

    assert uL_for("NEBNext Ultra II Q5 2x MM") == pytest.approx(10.0)
    assert uL_for("F1") == pytest.approx(1.0)
    assert uL_for("Total") == pytest.approx(20.0)


# --------------------------------------------------------------------------- #
# pcrs() -- light integration test of plate placement + master mixes + thermocycler
# --------------------------------------------------------------------------- #
def test_pcrs_end_to_end_plate_and_mm(pcr1_df, tmp_path):
    pivots, pcr1_mms, pcr1_thermo = sanger.pcrs(df=pcr1_df)

    # 3 samples sharing one primer pair -> one 96-well pivot block with 3 filled wells.
    well_pivot = pivots["96-well_ID"]
    filled = well_pivot.stack(future_stack=True).dropna()
    assert set(filled.values) == {"s1", "s2", "s3"}

    assert list(pcr1_mms.keys()) == [("F1", "R1")]
    assert list(pcr1_thermo.keys()) == ["F1_R1_65°C"]
