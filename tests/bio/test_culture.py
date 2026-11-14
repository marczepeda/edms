"""
Tests for edms.bio.culture

The numeric expectations are the master mixes from four Benchling entries
(MUZ378, MUZ382, MUZ386, MUZ387), which between them cover polybrene at
transduction, puromycin at selection, and DHT and ARV-110 at treatment, in
both ug/mL and nM. One formula has to reproduce all of them, because the only
thing that changes between those steps is the units and the volumes.
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
