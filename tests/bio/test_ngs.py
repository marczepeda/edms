"""
Tests for edms.bio.ngs

Covers pure-Python / pandas logic:
- group_boundaries(): consecutive-integer run grouping
- min_sec(): decimal-minutes -> (minutes, seconds)
- thermocycler(): per-primer-pair thermocycler step tables, including the annealing
  time derived from amplicon length and the ID-range title for grouped vs.
  non-consecutive reads
- pcr_mm() / pcr_mm_ultra(): NEB Q5 / NEBNext Ultra II master-mix uL calculations
- umis(): gDNA/molecule/read calculations for UMI-based genotyping
- hamming_distance() / hamming_distance_matrix(): pairwise sequence distance

All numeric expectations were hand-derived from the documented formulas and
cross-checked by running the functions directly.
"""
import math

import numpy as np
import pandas as pd
import pytest

from edms.bio import ngs


# --------------------------------------------------------------------------- #
# group_boundaries()
# --------------------------------------------------------------------------- #
def test_group_boundaries_empty():
    assert ngs.group_boundaries([]) == []


def test_group_boundaries_single_run():
    assert ngs.group_boundaries([1, 2, 3]) == [(1, 3)]


def test_group_boundaries_multiple_runs():
    assert ngs.group_boundaries([1, 2, 3, 5, 6, 9]) == [(1, 3), (5, 6), (9, 9)]


def test_group_boundaries_deduplicates_and_sorts():
    assert ngs.group_boundaries([3, 1, 2, 2, 1]) == [(1, 3)]


# --------------------------------------------------------------------------- #
# min_sec()
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("decimal_minutes,expected", [
    (2.5, (2, 30)),
    (0.25, (0, 15)),
    (0.5, (0, 30)),
    (3.0, (3, 0)),
])
def test_min_sec(decimal_minutes, expected):
    assert ngs.min_sec(decimal_minutes) == expected


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
        "PCR2 bp": [300, 300, 300],
    })


def test_thermocycler_step_table_structure(pcr1_df):
    dc = ngs.thermocycler(df=pcr1_df, n="1", cycles=30)
    assert list(dc.keys()) == ["F1_R1_65°C"]

    table = dc["F1_R1_65°C"]
    assert list(table["Temperature"]) == ["98°C", "98°C", "65°C", "72°C", "72°C", "4°C", ""]
    assert list(table["Repeat"]) == ["", "30 cycles", "30 cycles", "30 cycles", "", "", ""]


def test_thermocycler_anneal_time_derived_from_amplicon_length(pcr1_df):
    # bp=300 -> floor(300/500)/2+0.5 = 0.5 -> min_sec(0.5) = (0, 30) -> "30s"
    dc = ngs.thermocycler(df=pcr1_df, n="1", cycles=30)
    table = dc["F1_R1_65°C"]
    assert table["Time"].iloc[3] == "30s"


def test_thermocycler_anneal_time_minutes_and_seconds():
    df = pd.DataFrame({
        "ID": ["s1"], "PCR1 FWD": ["F1"], "PCR1 REV": ["R1"],
        "PCR1 ID": ["p1"], "PCR1 Tm": [65], "PCR2 bp": [2300],
    })
    # floor(2300/500)/2+0.5 = floor(4.6)/2+0.5 = 4/2+0.5 = 2.5 -> (2, 30) -> "2min 30s"
    dc = ngs.thermocycler(df=df, n="1", cycles=30)
    table = dc["F1_R1_65°C"]
    assert table["Time"].iloc[3] == "2min 30s"


def test_thermocycler_title_ranges_consecutive_ids(pcr1_df):
    dc = ngs.thermocycler(df=pcr1_df, n="1", cycles=30)
    table = dc["F1_R1_65°C"]
    # 3 consecutive rows sharing the same primers collapse to a "start -> end" range.
    assert table.index[-1] == "F1_R1: p1 -> p3"


def test_thermocycler_title_lists_non_consecutive_ids():
    df = pd.DataFrame({
        "ID": ["s1", "s2", "s3"],
        "PCR1 FWD": ["F1", "F2", "F1"],
        "PCR1 REV": ["R1", "R2", "R1"],
        "PCR1 ID": ["p1", "p2", "p3"],
        "PCR1 Tm": [65, 60, 65],
        "PCR2 bp": [300, 300, 300],
    })
    dc = ngs.thermocycler(df=df, n="1", cycles=30)
    # s1 (row0) and s3 (row2) share F1/R1 but are not adjacent -> comma-separated, not a range.
    assert dc["F1_R1_65°C"].index[-1] == "F1_R1: p1, p3"
    assert dc["F2_R2_60°C"].index[-1] == "F2_R2: p2"


def test_thermocycler_pcr2_defaults_and_fwd_override(pcr1_df):
    df = pcr1_df.rename(columns={"PCR1 FWD": "PCR2 FWD_orig"})
    df["PCR2 ID"] = df["PCR1 ID"]
    df["PCR2 Tm"] = df["PCR1 Tm"]
    df["PCR2 REV"] = "R1"
    dc = ngs.thermocycler(df=df, n="2")
    # thermocycler() forces PCR2 FWD to the literal string 'PCR2-FWD' for n='2'.
    assert list(dc.keys()) == ["PCR2-FWD_R1_65°C"]


def test_thermocycler_invalid_n_raises(pcr1_df):
    with pytest.raises(ValueError):
        ngs.thermocycler(df=pcr1_df, n="3")


# --------------------------------------------------------------------------- #
# pcr_mm() / pcr_mm_ultra()
# --------------------------------------------------------------------------- #
def test_pcr_mm_default_uL_calculations():
    primers = pd.Series([2], index=pd.MultiIndex.from_tuples([("F1", "R1")]))
    mm = ngs.pcr_mm(primers=primers, template="gDNA", template_uL=5)
    table = mm[("F1", "R1")]

    def uL_for(component):
        return table.loc[table["Component"] == component, "uL"].iloc[0]

    assert uL_for("5x Q5 Reaction Buffer") == pytest.approx(5.0)   # 1/5*25
    assert uL_for("dNTPs") == pytest.approx(0.5)                    # 0.2/10*25
    assert uL_for("F1") == pytest.approx(1.25)                      # 0.5/10*25
    assert uL_for("gDNA") == pytest.approx(5.0)                     # template_uL passthrough
    assert uL_for("Q5 Polymerase") == pytest.approx(0.25)           # 0.02/2*25
    assert uL_for("Total") == pytest.approx(25.0)


def test_pcr_mm_scales_by_reaction_count_and_mm_x():
    primers = pd.Series([3], index=pd.MultiIndex.from_tuples([("F1", "R1")]))
    mm = ngs.pcr_mm(primers=primers, template="gDNA", template_uL=5, mm_x=1.1)
    table = mm[("F1", "R1")]
    total_row = table[table["Component"] == "Total"].iloc[0]
    assert total_row["uL MM"] == pytest.approx(round(25.0 * 3 * 1.1, 2))


def test_pcr_mm_ultra_default_uL_calculations():
    primers = pd.Series([2], index=pd.MultiIndex.from_tuples([("F1", "R1")]))
    mm = ngs.pcr_mm_ultra(primers=primers, template="gDNA", template_uL=5)
    table = mm[("F1", "R1")]

    def uL_for(component):
        return table.loc[table["Component"] == component, "uL"].iloc[0]

    assert uL_for("NEBNext Ultra II Q5 2x MM") == pytest.approx(10.0)  # 1/2*20
    assert uL_for("F1") == pytest.approx(1.0)                           # 0.5/10*20
    assert uL_for("gDNA") == pytest.approx(5.0)
    assert uL_for("Total") == pytest.approx(20.0)


# --------------------------------------------------------------------------- #
# umis()
# --------------------------------------------------------------------------- #
def test_umis_default_calculation():
    ug, molecules, reads_needed, reads_needed_samples = ngs.umis(genotypes=10, samples=2)
    # ug = genotypes * cell_coverage * ug_gDNA_per_cell = 10*1000*6e-6
    assert ug == pytest.approx(10 * 1000 * 6e-6)
    # molecules = genotypes * cell_coverage * ploidy_per_cell
    assert molecules == 10 * 1000 * 2
    # reads_needed = molecules * umi_coverage
    assert reads_needed == 20000 * 5
    assert reads_needed_samples == reads_needed * 2


def test_umis_scales_with_custom_parameters():
    ug, molecules, reads_needed, reads_needed_samples = ngs.umis(
        genotypes=5, samples=1, cell_coverage=500, ploidy_per_cell=1, umi_coverage=10,
    )
    assert molecules == 5 * 500 * 1
    assert reads_needed == molecules * 10
    assert reads_needed_samples == reads_needed


# --------------------------------------------------------------------------- #
# hamming_distance() / hamming_distance_matrix()
# --------------------------------------------------------------------------- #
def test_hamming_distance_counts_mismatches():
    assert ngs.hamming_distance("ACGT", "ACGA") == 1
    assert ngs.hamming_distance("AAAA", "TTTT") == 4
    assert ngs.hamming_distance("ACGT", "ACGT") == 0


def test_hamming_distance_requires_equal_length():
    with pytest.raises(ValueError, match="equal length"):
        ngs.hamming_distance("ACG", "ACGT")


def test_hamming_distance_matrix_is_symmetric_with_zero_diagonal():
    df = pd.DataFrame({"ID": ["a", "b", "c"], "seq": ["ACGT", "ACGA", "TTTT"]})
    dm = ngs.hamming_distance_matrix(df=df, id="ID", seqs="seq")

    values = dm.to_numpy()
    assert np.array_equal(values, values.T)
    assert np.all(np.diag(values) == 0)
    assert dm.loc["ACGT", "b"] == 1
    assert dm.loc["ACGT", "c"] == 3


def test_hamming_distance_matrix_saves_file(tmp_path):
    df = pd.DataFrame({"ID": ["a", "b"], "seq": ["ACGT", "ACGA"]})
    out_dir = tmp_path / "out"
    ngs.hamming_distance_matrix(df=df, id="ID", seqs="seq", dir=str(out_dir), file="dm.csv")
    assert (out_dir / "dm.csv").is_file()
