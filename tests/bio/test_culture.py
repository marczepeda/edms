"""
Tests for edms.bio.culture

The numeric expectations are the master mixes from five Benchling entries
(MUZ361, MUZ378, MUZ382, MUZ386, MUZ387), which between them cover polybrene
at transduction, puromycin at selection, and DHT and ARV-110 at treatment, in
both ug/mL and nM. One formula has to reproduce all of them, because the only
thing that changes between those steps is the units and the volumes.

MUZ361 also seeds fifteen clones from confluency rather than a cell count,
and splits each harvest between a transduction plate and a competitive growth
plate; its seed table is the expectation for seed_plan().
"""
import pandas as pd
import pytest

from edms.bio import culture


def _value(df, component, unit=None):
    rows = df[df["Component"] == component]
    if unit is not None:
        rows = rows[rows["Unit"] == unit]
    return rows["Value"].iloc[0]


# --------------------------------------------------------------------------- #
# drug_mm(): checked against the recorded master mixes
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("drug,stock,desired,drug_mL,cell_mL,wells,mm_x,conc,stock_mL,media_mL", [
    # MUZ378 Day 2: puromycin selection
    ("Puro", 10000, 2, 1.7, 0.3, 45, 1.1, 2.3529, 0.0198, 84.1302),
    # MUZ378 Day 0: polybrene at transduction
    ("Polybrene", 10000, 8, 1.0, 1.0, 18, 1.2, 16.0, 0.0346, 21.5654),
    # MUZ386 Day 5: DHT, in nM
    ("DHT", 10000, 1, 0.5, 0.5, 168, 1.1, 2.0, 0.0185, 92.3815),
    # MUZ386 Day 5: ARV-110, in nM
    ("ARV-110", 100000, 100, 0.5, 0.5, 21, 1.1, 200.0, 0.0231, 11.5269),
])
def test_drug_mm_reproduces_the_recorded_mixes(drug, stock, desired, drug_mL,
                                               cell_mL, wells, mm_x, conc,
                                               stock_mL, media_mL):
    df = culture.drug_mm(stock=stock, desired=desired, wells=wells,
                         drug_media_mL=drug_mL, cell_media_mL=cell_mL,
                         mm_x=mm_x, drug=drug)
    assert _value(df, f"{drug}, MM") == pytest.approx(conc, abs=1e-4)
    assert _value(df, f"{drug} stock, MM") == pytest.approx(stock_mL, abs=1e-4)
    assert _value(df, "Media, MM") == pytest.approx(media_mL, abs=1e-4)


def test_drug_mm_concentrates_by_the_volume_ratio():
    """The whole point: the drug covers the cell volume added on top of it."""
    df = culture.drug_mm(stock=1000, desired=10, wells=1,
                         drug_media_mL=1.0, cell_media_mL=1.0)
    assert _value(df, "Drug, MM") == pytest.approx(20.0)


def test_drug_mm_is_unchanged_when_cells_add_no_volume():
    df = culture.drug_mm(stock=1000, desired=10, wells=1,
                         drug_media_mL=1.0, cell_media_mL=0.0)
    assert _value(df, "Drug, MM") == pytest.approx(10.0)


def test_drug_mm_media_and_stock_sum_to_the_total():
    df = culture.drug_mm(stock=10000, desired=2, wells=45,
                         drug_media_mL=1.7, cell_media_mL=0.3, mm_x=1.1)
    total = 1.7 * 45 * 1.1
    assert (_value(df, "Media, MM") + _value(df, "Drug stock, MM")
            == pytest.approx(total, abs=1e-4))


def test_drug_mm_carries_the_unit_through():
    df = culture.drug_mm(stock=10000, desired=1, wells=10, unit="nM", drug="DHT")
    assert _value(df, "DHT, stock", "nM") == 10000
    assert _value(df, "DHT, MM", "nM") == pytest.approx(2.0)


def test_drug_mm_rejects_a_desired_above_the_stock():
    with pytest.raises(ValueError, match="above the stock"):
        culture.drug_mm(stock=10, desired=100, wells=4)


def test_drug_mm_rejects_a_zero_drug_volume():
    with pytest.raises(ValueError, match="drug_media_mL"):
        culture.drug_mm(stock=1000, desired=1, wells=4, drug_media_mL=0)


def test_drug_mm_rejects_a_zero_stock():
    with pytest.raises(ValueError, match="stock"):
        culture.drug_mm(stock=0, desired=1, wells=4)


# --------------------------------------------------------------------------- #
# cells_mm()
# --------------------------------------------------------------------------- #
def test_cells_mm_reproduces_the_recorded_mix():
    """MUZ378: 2.5e5 cells/mL counted, 9.6e5 wanted per well, 18 wells."""
    df = culture.cells_mm(flow_per_mL=2.5e5, cells_per_well=9.6e5, wells=18,
                          cell_media_mL=1.0, mm_x=1.1)
    assert _value(df, "Culture to spin, MM") == pytest.approx(76.032, abs=1e-3)
    assert _value(df, "Resuspend in, MM") == pytest.approx(19.8, abs=1e-3)


def test_cells_mm_resuspension_delivers_the_wanted_density():
    df = culture.cells_mm(flow_per_mL=2.5e5, cells_per_well=9.6e5, wells=18)
    per_mL = _value(df, "Resuspension, MM")
    assert per_mL * 1.0 == pytest.approx(9.6e5)  # 1.0 mL transferred per well


def test_cells_mm_rejects_a_zero_count():
    with pytest.raises(ValueError, match="flow_per_mL"):
        culture.cells_mm(flow_per_mL=0, cells_per_well=1e5, wells=4)


# --------------------------------------------------------------------------- #
# transduction()
# --------------------------------------------------------------------------- #
def test_transduction_returns_both_mixes():
    tables = culture.transduction(wells=18, flow_per_mL=2.5e5, cells_per_well=9.6e5)
    assert set(tables) == {"Polybrene", "Cells"}
    assert all(isinstance(df, pd.DataFrame) for df in tables.values())


def test_transduction_matches_the_recorded_entry():
    """MUZ378 Day 0, both columns of the Transduction MM1 table."""
    tables = culture.transduction(wells=18, flow_per_mL=2.5e5,
                                  cells_per_well=9.6e5, virus_mL_per_well=1.0,
                                  cell_media_mL=1.0)
    pb, cells = tables["Polybrene"], tables["Cells"]
    assert _value(pb, "Polybrene, MM") == pytest.approx(16.0)
    assert _value(pb, "Polybrene stock, MM") == pytest.approx(0.0346, abs=1e-4)
    assert _value(pb, "Media, MM") == pytest.approx(21.5654, abs=1e-4)
    assert _value(cells, "Culture to spin, MM") == pytest.approx(76.032, abs=1e-3)
    assert _value(cells, "Resuspend in, MM") == pytest.approx(19.8, abs=1e-3)


def test_transduction_polybrene_covers_the_virus_volume():
    """Polybrene is made up in the cell media but must dose the whole well."""
    tables = culture.transduction(wells=6, flow_per_mL=1e6, cells_per_well=1e5,
                                  virus_mL_per_well=3.0, cell_media_mL=1.0,
                                  polybrene_desired=8)
    assert _value(tables["Polybrene"], "Polybrene, MM") == pytest.approx(32.0)


def test_transduction_uses_separate_multipliers():
    """Polybrene media is cheap to overmake; cells are not."""
    tables = culture.transduction(wells=10, flow_per_mL=1e6, cells_per_well=1e5,
                                  mm_x_polybrene=1.5, mm_x_cells=1.0)
    assert _value(tables["Polybrene"], "MM") == 1.5
    assert _value(tables["Cells"], "MM") == 1.0


# --------------------------------------------------------------------------- #
# MUZ361: polybrene at three scales, including a 0.25 mL well and a 1.4x mix
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("wells,mm_x,media_mL,stock_mL,leftover_mL", [
    (106, 1.4, 0.5, 0.1187, 74.0813),   # MM8
    (12, 1.2, 0.5, 0.0115, 7.1885),     # MM1
    (64, 1.4, 0.25, 0.0358, 22.3642),   # MM2, half-volume wells
])
def test_muz361_polybrene_mixes(wells, mm_x, media_mL, stock_mL, leftover_mL):
    """Same formula, three scales, two well volumes, two multipliers."""
    df = culture.drug_mm(stock=10000, desired=8, wells=wells,
                         drug_media_mL=media_mL, cell_media_mL=media_mL,
                         mm_x=mm_x, drug="Polybrene")
    assert _value(df, "Polybrene, MM") == 16.0
    assert _value(df, "Polybrene stock, MM") == stock_mL
    assert _value(df, "Media, MM") == leftover_mL


def test_muz361_cells_at_quarter_volume():
    """Halving the well volume halves the resuspension but not the spin."""
    full = culture.cells_mm(flow_per_mL=9.5e5, cells_per_well=1.8e5, wells=9,
                            cell_media_mL=0.5, mm_x=1.1)
    half = culture.cells_mm(flow_per_mL=9.5e5, cells_per_well=1.8e5, wells=9,
                            cell_media_mL=0.25, mm_x=1.1)
    assert _value(full, "Culture to spin, MM") == 1.8758
    assert _value(half, "Culture to spin, MM") == 1.8758
    assert _value(full, "Resuspend in, MM") == 4.95
    assert _value(half, "Resuspend in, MM") == 2.475


# --------------------------------------------------------------------------- #
# per_mL()
# --------------------------------------------------------------------------- #
def test_per_mL_scales_linearly_with_confluency():
    assert culture.per_mL("15cm", 100) == culture.CONFLUENT_PER_ML["15-cm"]
    assert culture.per_mL("15cm", 50) == culture.CONFLUENT_PER_ML["15-cm"] / 2


def test_per_mL_accepts_vessel_aliases():
    assert culture.per_mL("15-cm", 80) == culture.per_mL("15cm", 80)
    assert culture.per_mL("6 well", 80) == culture.per_mL("6-well", 80)


@pytest.mark.parametrize("confluency", [0, -10, 101])
def test_per_mL_rejects_impossible_confluency(confluency):
    with pytest.raises(ValueError, match="confluency"):
        culture.per_mL("15cm", confluency)


def test_per_mL_rejects_vessel_without_a_calibration():
    """A 10 cm dish is a real vessel; nobody has told us what it yields."""
    with pytest.raises(KeyError, match="No confluent density"):
        culture.per_mL("10cm", 100)


def test_per_mL_override_adds_a_vessel():
    assert culture.per_mL("10cm", 50, confluent_per_mL={"10cm": 1e6}) == 5e5


# --------------------------------------------------------------------------- #
# seed_plan(): checked against the MUZ361 seed table
# --------------------------------------------------------------------------- #
# sample, confluency, vessel, transduction spin, CG spin, sum
MUZ361_SEED = [
    ("L1", 85, "15cm", 0.3812, 0.2541, 0.6353),
    ("L102", 100, "15cm", 0.3240, 0.2160, 0.5400),
    ("L201", 75, "6-well", 0.7200, 0.4800, 1.2000),
    ("L202", 85, "6-well", 0.6353, 0.4235, 1.0588),
    ("L203", 50, "6-well", 1.0800, 0.7200, 1.8000),
    ("L207", 60, "6-well", 0.9000, 0.6000, 1.5000),
    ("L208", 45, "6-well", 1.2000, 0.8000, 2.0000),
    ("L209", 90, "6-well", 0.6000, 0.4000, 1.0000),
    ("L211", 80, "6-well", 0.6750, 0.4500, 1.1250),
    ("L217", 85, "6-well", 0.6353, 0.4235, 1.0588),
    ("L224", 90, "6-well", 0.6000, 0.4000, 1.0000),
]


@pytest.fixture
def muz361_samples():
    return pd.DataFrame([{"Sample": s, "Confluency %": c, "Vessel": v}
                         for s, c, v, *_ in MUZ361_SEED])


def test_seed_plan_reproduces_muz361(muz361_samples):
    """Eleven clones, two vessels, confluencies from 45 to 100 percent."""
    plan = culture.seed_plan(muz361_samples, cells_per_well=1.8e5, wells=3,
                             arm="Transduce")
    plan = culture.seed_plan(plan, cells_per_well=1.2e5, wells=3, arm="CG")
    for (sample, _, _, transduce, cg, total) in MUZ361_SEED:
        row = plan[plan["Sample"] == sample].iloc[0]
        assert row["Transduce spin, mL"] == transduce, sample
        assert row["CG spin, mL"] == cg, sample
        assert row["Sum spin, mL"] == total, sample


def test_seed_plan_resuspension_is_the_same_for_every_sample(muz361_samples):
    """Only the spin varies with confluency; every sample ends up alike."""
    plan = culture.seed_plan(muz361_samples, cells_per_well=1.8e5, wells=3)
    assert set(plan["Seed resuspend, mL"]) == {0.9}


def test_seed_plan_leaves_one_arm_without_a_total(muz361_samples):
    """A sum column with one term in it would just be noise."""
    plan = culture.seed_plan(muz361_samples, cells_per_well=1.8e5, wells=3)
    assert "Sum spin, mL" not in plan.columns


def test_seed_plan_total_is_last(muz361_samples):
    plan = culture.seed_plan(muz361_samples, cells_per_well=1.8e5, wells=3,
                             arm="Transduce")
    plan = culture.seed_plan(plan, cells_per_well=1.2e5, wells=3, arm="CG")
    assert list(plan.columns)[-1] == "Sum spin, mL"


def test_seed_plan_does_not_touch_the_input(muz361_samples):
    before = muz361_samples.copy()
    culture.seed_plan(muz361_samples, cells_per_well=1.8e5, wells=3)
    pd.testing.assert_frame_equal(muz361_samples, before)


def test_seed_plan_halves_the_spin_when_twice_as_confluent():
    df = pd.DataFrame([{"Sample": "a", "Confluency %": 45, "Vessel": "6-well"},
                       {"Sample": "b", "Confluency %": 90, "Vessel": "6-well"}])
    plan = culture.seed_plan(df, cells_per_well=1.8e5, wells=3)
    spins = list(plan["Seed spin, mL"])
    assert spins[0] == pytest.approx(spins[1] * 2)


def test_seed_plan_prefers_a_measured_concentration(muz361_samples):
    """A line that was counted should not be estimated."""
    df = muz361_samples.copy()
    df["Cells, flow"] = [None] * len(df)
    df.loc[0, "Cells, flow"] = 1e6
    plan = culture.seed_plan(df, cells_per_well=1.8e5, wells=3)
    assert plan.loc[0, "Seed spin, mL"] == round(1.8e5 * 3 * 1.2 / 1e6, 4)
    # everything else still comes from confluency
    assert plan.loc[1, "Seed spin, mL"] == 0.3240


def test_seed_plan_scale_column_seeds_a_sample_heavier(muz361_samples):
    """MUZ361 put 2x the intended cells on one line's transduction plate."""
    df = muz361_samples.copy()
    df["Heavy"] = [2.0] + [1.0] * (len(df) - 1)
    plan = culture.seed_plan(df, cells_per_well=1.8e5, wells=3, scale_col="Heavy")
    assert plan.loc[0, "Seed spin, mL"] == 2 * 0.3812
    assert plan.loc[0, "Seed resuspend, mL"] == 1.8


def test_seed_plan_requires_its_columns():
    with pytest.raises(KeyError, match="Confluency %"):
        culture.seed_plan(pd.DataFrame({"Sample": ["a"]}), cells_per_well=1e5,
                          wells=3)


def test_seed_plan_rejects_no_wells(muz361_samples):
    with pytest.raises(ValueError, match="wells"):
        culture.seed_plan(muz361_samples, cells_per_well=1e5, wells=0)


def test_seed_plan_rejects_a_measured_zero(muz361_samples):
    df = muz361_samples.copy()
    df["Cells, flow"] = 0.0
    with pytest.raises(ValueError, match="greater than 0"):
        culture.seed_plan(df, cells_per_well=1e5, wells=3)


def test_seed_plan_writes_a_file(muz361_samples, tmp_path):
    out = tmp_path / "plan.csv"
    culture.seed_plan(muz361_samples, cells_per_well=1e5, wells=3, file=str(out))
    assert out.exists()


# --------------------------------------------------------------------------- #
# competitive_growth()
# --------------------------------------------------------------------------- #
def test_competitive_growth_returns_both_arms_and_the_competitor(muz361_samples):
    tables = culture.competitive_growth(muz361_samples)
    assert set(tables) == {"Seed plan", "LCG"}
    plan = tables["Seed plan"]
    assert {"Transduce spin, mL", "CG spin, mL", "Sum spin, mL"} <= set(plan.columns)


def test_competitive_growth_matches_muz361_seed_table(muz361_samples):
    plan = culture.competitive_growth(muz361_samples)["Seed plan"]
    for (sample, _, _, transduce, cg, total) in MUZ361_SEED:
        row = plan[plan["Sample"] == sample].iloc[0]
        assert (row["Transduce spin, mL"], row["CG spin, mL"],
                row["Sum spin, mL"]) == (transduce, cg, total), sample


def test_competitive_growth_competitor_matches_muz361():
    """15 clones plus one day 0 sample, at a 1.2x multiplier of its own."""
    df = pd.DataFrame([{"Sample": f"c{i}", "Confluency %": 80, "Vessel": "6-well"}
                       for i in range(15)])
    lcg = culture.competitive_growth(df)["LCG"]
    assert _value(lcg, "Samples") == 16
    assert _value(lcg, "LCG to spin, MM") == 4.1472
    assert _value(lcg, "Resuspend in, MM") == 17.28


def test_competitive_growth_skips_the_transduction_arm(muz361_samples):
    plan = culture.competitive_growth(muz361_samples, transduce_wells=0)["Seed plan"]
    assert "Transduce spin, mL" not in plan.columns
    # and with one arm there is nothing to total
    assert "Sum spin, mL" not in plan.columns


def test_competitive_growth_competitor_can_be_counted(muz361_samples):
    lcg = culture.competitive_growth(muz361_samples,
                                     competitor_per_mL=5e5)["LCG"]
    assert _value(lcg, "LCG, source") == "measured"
    assert _value(lcg, "LCG, cells") == 500000.0


def test_competitive_growth_competitor_is_named(muz361_samples):
    tables = culture.competitive_growth(muz361_samples, competitor="LNCaP-GFP")
    assert "LNCaP-GFP" in tables
    assert _value(tables["LNCaP-GFP"], "LNCaP-GFP, source") == "100% of a 15-cm"


def test_competitive_growth_rejects_an_empty_table():
    with pytest.raises(ValueError, match="at least one sample"):
        culture.competitive_growth(pd.DataFrame(columns=["Sample", "Confluency %",
                                                         "Vessel"]))


def test_competitive_growth_rejects_a_measured_zero(muz361_samples):
    with pytest.raises(ValueError, match="greater than 0"):
        culture.competitive_growth(muz361_samples, competitor_per_mL=0)


def test_competitive_growth_writes_a_file(muz361_samples, tmp_path):
    out = tmp_path / "cg.xlsx"
    culture.competitive_growth(muz361_samples, file=str(out))
    assert out.exists()
