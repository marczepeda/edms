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


def test_pcr_mm_syber_adds_stain_from_water():
    primers = pd.Series([2], index=pd.MultiIndex.from_tuples([("F1", "R1")]))
    plain = ngs.pcr_mm(primers=primers, template="gDNA", template_uL=5)[("F1", "R1")]
    table = ngs.pcr_mm(primers=primers, template="gDNA", template_uL=5, syber=True)[("F1", "R1")]

    def uL_for(t, component):
        return t.loc[t["Component"] == component, "uL"].iloc[0]

    assert list(table["Component"]).index("10x SYBR Green DNA Stain") == list(table["Component"]).index("gDNA") - 1
    assert uL_for(table, "10x SYBR Green DNA Stain") == pytest.approx(2.5)   # 1/10*25
    assert uL_for(table, "Nuclease-free H2O") == pytest.approx(uL_for(plain, "Nuclease-free H2O") - 2.5)
    assert table["uL"].iloc[:-1].sum() == pytest.approx(25.0)
    assert list(table.index) == list(range(1, 10))


def test_pcr_mm_ultra_syber_adds_stain_from_water():
    primers = pd.Series([2], index=pd.MultiIndex.from_tuples([("F1", "R1")]))
    table = ngs.pcr_mm_ultra(primers=primers, template="gDNA", template_uL=5, syber=True)[("F1", "R1")]
    assert table.loc[table["Component"] == "10x SYBR Green DNA Stain", "uL"].iloc[0] == pytest.approx(2.0)  # 1/10*20
    assert table["uL"].iloc[:-1].sum() == pytest.approx(20.0)


# --------------------------------------------------------------------------- #
# pcrs()
# --------------------------------------------------------------------------- #
def _pcrs_df(n=3):
    return pd.DataFrame({
        "ID": [f"g{i}" for i in range(n)],
        "PCR1 ID": [f"p1_{i}" for i in range(n)],
        "PCR1 FWD": ["F1"] * n, "PCR1 REV": ["R1"] * n,
        "PCR1 Tm": [65] * n, "PCR1 bp": [200] * n,
        "PCR2 ID": [f"p2_{i}" for i in range(n)],
        "PCR2 FWD": ["P5"] * n, "PCR2 REV": [f"P7_{i}" for i in range(n)],
        "PCR2 Tm": [65] * n, "PCR2 bp": [300] * n,
    })


@pytest.mark.parametrize("ultra", [False, True])
def test_pcrs_syber_only_in_pcr1(ultra):
    out = ngs.pcrs(df=_pcrs_df(), syber=True, ultra=ultra)
    pcr1_mms, pcr2_mms = out[1], out[2]
    assert all("10x SYBR Green DNA Stain" in list(mm["Component"]) for mm in pcr1_mms.values())
    assert not any("10x SYBR Green DNA Stain" in list(mm["Component"]) for mm in pcr2_mms.values())
    for mm in pcr1_mms.values():
        assert (mm["uL"] >= 0).all()


def _rows_cols(pivot):
    return set(pivot.index.get_level_values("row")), set(pivot.columns)


@pytest.mark.parametrize("split", [True, False])
def test_pcrs_uses_outer_wells_by_default(split):
    pivots = ngs.pcrs(df=_pcrs_df(n=70), split_pcr1_primers=split)[0]
    for key in (["96-well_PCR1 ID", "96-well_PCR2 ID"] if split else ["PCR1 ID", "PCR2 ID"]):
        rows, cols = _rows_cols(pivots[key])
        assert "A" in rows and 1 in cols


@pytest.mark.parametrize("split", [True, False])
@pytest.mark.parametrize("inner1,inner2", [(True, False), (False, True), (True, True)])
def test_pcrs_inner_flags_exclude_outer_wells(split, inner1, inner2):
    pivots = ngs.pcrs(df=_pcrs_df(n=70), split_pcr1_primers=split, inner1=inner1, inner2=inner2)[0]
    keys = ["96-well_PCR1 ID", "96-well_PCR2 ID"] if split else ["PCR1 ID", "PCR2 ID"]
    for key, inner in zip(keys, [inner1, inner2]):
        rows, cols = _rows_cols(pivots[key])
        if inner:
            assert rows <= set("BCDEFG") and cols <= set(range(2, 12))
        else:
            assert "A" in rows and 1 in cols
        assert pivots[key].notna().sum().sum() == 70   # 60 inner wells -> spills onto plate 2


def _wells(pivot):
    """{(row, column): value} for filled wells on the first plate"""
    first = pivot.loc[pivot.index.get_level_values(0)[0]]
    return {(r, c): v for r, row in first.iterrows() for c, v in row.items() if pd.notna(v)}


@pytest.mark.parametrize("total_uL,expected", [
    (100, [50, 50]),
    (240, [50, 50, 50, 50, 40]),
])
def test_pcrs_splits_large_volumes_across_wells(total_uL, expected):
    pivots = ngs.pcrs(df=_pcrs_df(n=1), pcr1_total_uL=total_uL)[0]
    ids, uLs = _wells(pivots["96-well_PCR1 ID"]), _wells(pivots["96-well_PCR1 uL"])
    wells = [("A", c) for c in range(1, len(expected) + 1)]
    assert sorted(ids) == wells and set(ids.values()) == {"p1_0"}
    assert [uLs[w] for w in wells] == expected
    assert "96-well_PCR2 uL" not in pivots   # PCR2 stays at 20 uL -> 1 well each


def test_pcrs_split_wells_wrap_to_next_row():
    pivots = ngs.pcrs(df=_pcrs_df(n=3), pcr1_total_uL=240, split_pcr1_primers=False)[0]
    ids, uLs = _wells(pivots["PCR1 ID"]), _wells(pivots["PCR1 uL"])
    third = [("A", 11), ("A", 12), ("B", 1), ("B", 2), ("B", 3)]
    assert [ids[w] for w in third] == ["p1_2"] * 5
    assert [uLs[w] for w in third] == [50, 50, 50, 50, 40]


def test_pcrs_split_wells_respect_inner():
    pivots = ngs.pcrs(df=_pcrs_df(n=4), pcr2_total_uL=120, inner2=True)[0]   # 3 wells each; 10 inner wells per row
    ids = _wells(pivots["96-well_PCR2 ID"])
    assert [ids[w] for w in [("B", 11), ("C", 2), ("C", 3)]] == ["p2_3"] * 3
    assert not any(r in ("A", "H") or c in (1, 12) for r, c in ids)


def test_pcrs_no_uL_pivots_by_default():
    pivots = ngs.pcrs(df=_pcrs_df())[0]
    assert not any(k.endswith("uL") for k in pivots)


def _pcr1_only_df(n=4):
    df = _pcrs_df(n=n)
    return df.drop(columns=[c for c in df.columns if c.startswith("PCR2")])


@pytest.mark.parametrize("split", [True, False])
@pytest.mark.parametrize("ultra", [False, True])
def test_pcrs_without_pcr2_skips_pcr2_and_needs_no_pcr2_columns(tmp_path, split, ultra):
    out_file = tmp_path / "plan.xlsx"
    pivots, pcr1_mms, pcr2_mms, pcr1_thermo, pcr2_thermo = ngs.pcrs(
        df=_pcr1_only_df(), pcr2=False, split_pcr1_primers=split, ultra=ultra, file=str(out_file), pcr2_total_uL=100,
    )
    assert pcr1_mms and pcr1_thermo
    assert pcr2_mms == {} and pcr2_thermo == {}
    assert not any("PCR2" in k for k in pivots)
    assert not any("PCR2" in s for s in pd.ExcelFile(out_file).sheet_names)
    # Extension time falls back to PCR1 bp (200 bp -> 30s)
    assert list(pcr1_thermo.values())[0].loc["2", "Time"].tolist()[-1] == "30s"


def test_pcrs_without_pcr2_umi():
    df = _pcr1_only_df()
    df["UMI"] = [True, False, True, False]
    out = ngs.pcrs(df=df, pcr2=False)
    assert out[2] == {} and out[-1] == {}
    assert out[3] and out[4] and out[5]   # UMI PCR1, PCR1.5, & non-UMI PCR1 thermocyclers


# --------------------------------------------------------------------------- #
# excel_colors() / pcrs() Excel styling
# --------------------------------------------------------------------------- #
def test_excel_colors_cycles_set3_and_reuses_colors():
    colors = ngs.excel_colors([f"v{i}" for i in range(13)] + ["v0"], "Set3")
    assert len(colors) == 13
    assert colors["v0"] == "#8dd3c7" and colors["v12"] == colors["v0"]   # Set3 has 12 colors -> cycles


def test_excel_colors_accepts_color_list():
    assert ngs.excel_colors(["x", "y", "z"], ["#ff0000", "blue"]) == {"x": "#ff0000", "y": "#0000ff", "z": "#ff0000"}


def test_pcrs_excel_alternating_rows_and_value_colors(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    df = _pcrs_df(n=4)
    df.loc[3, "ID"] = "g0"   # same gDNA ID in two samples
    out = tmp_path / "plan.xlsx"
    ngs.pcrs(df=df, file=str(out), pcr1_total_uL=100)
    wb = openpyxl.load_workbook(out)

    def rgb(cell):
        return cell.fill.start_color.rgb[-6:].upper() if cell.fill.fill_type else None

    ws = wb["96-well_ID"]   # header row 1; A row 2: g0 g0 g1 g1 g2 g2 g0 g0
    assert [ws.cell(row=2, column=c).value for c in range(3, 11)] == ["g0", "g0", "g1", "g1", "g2", "g2", "g0", "g0"]
    fills = [rgb(ws.cell(row=2, column=c)) for c in range(3, 11)]
    assert fills[0] == fills[1] == fills[6] == fills[7] == "8DD3C7"
    assert len({fills[0], fills[2], fills[4]}) == 3

    ws = wb["F1_R1"]   # master mix: alternating white / light gray rows
    assert [rgb(ws.cell(row=r, column=2)) for r in range(2, 6)] == ["FFFFFF", "EDEDED", "FFFFFF", "EDEDED"]


def test_pcrs_excel_autofit_and_left_aligned_thermocycler_reactions(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    df = _pcrs_df(n=4)
    df.loc[0, "ID"] = "MUZ360-701-long-gDNA-name"
    out = tmp_path / "plan.xlsx"
    ngs.pcrs(df=df, file=str(out))
    wb = openpyxl.load_workbook(out)

    ws = wb["96-well_ID"]   # columns: plate label, row, 1, 2, ...
    assert ws.column_dimensions["C"].width >= len("MUZ360-701-long-gDNA-name")
    assert ws.column_dimensions["A"].width >= len("96-well plate (PCR1)")
    assert ws.column_dimensions["D"].width < ws.column_dimensions["C"].width   # fits per column

    plan = wb["NGS Plan"]   # shared sheet keeps the widest fit across tables
    assert plan.column_dimensions["C"].width >= len("MUZ360-701-long-gDNA-name")

    ws = wb["F1_R1_65°C"]   # thermocycler: no autofit; reactions row left aligned
    assert "A" not in ws.column_dimensions or not ws.column_dimensions["A"].customWidth
    last = ws.cell(row=ws.max_row, column=1)
    assert last.value == "F1_R1: p1_0 -> p1_3"
    assert last.alignment.horizontal == "left"


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
    ngs.hamming_distance_matrix(df=df, id="ID", seqs="seq", file=str(out_dir / "dm.csv"))
    assert (out_dir / "dm.csv").is_file()
