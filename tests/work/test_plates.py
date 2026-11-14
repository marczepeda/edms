"""
Tests for edms.work.plates

Covers the vessel table, the spellings people actually type, growth-area
scaling, and the two transfection wrappers. The property that matters most is
that a run in a calculator's anchor format (96-well for PE3, 6-well for virus)
comes out with exactly the numbers edms.bio.transfect produced before the
plate layer existed, so every anchor case is checked against a direct call.
"""
import pandas as pd
import pytest

from edms.bio import transfect
from edms.work import plates


@pytest.fixture
def packaging():
    """VSVG + GagPol + two transfer plasmids."""
    return pd.DataFrame({
        "Plasmid": ["pMUZ26.6", "pMUZ26.7", "pMUZ310", "pMUZ311"],
        "Description": ["VSVG", "GagPol", "transfer A", "transfer B"],
        "Colony": ["c1", "c1", "c2", "c1"],
        "ng/uL": [412.0, 388.5, 265.3, 301.8],
    })


@pytest.fixture
def pe3_tables():
    """plasmids / epegRNAs / ngRNAs in the shape transfect.PE3 expects."""
    plasmids = pd.DataFrame({
        "Plasmid": ["pMUZ101", "pMUZ102", "pMUZ201", "pMUZ202", "pMUZ86.7"],
        "Description": ["epeg A", "epeg B", "ng A", "ng B", "PEmax"],
        "Colony": ["c1", "c2", "c1", "c2", "c1"],
        "ng/uL": [128.4, 96.2, 208.9, 176.3, 305.1],
    })
    epegRNAs = pd.DataFrame({"pegRNA_number": [1, 2],
                             "Name": ["pMUZ101", "pMUZ102"]})
    ngRNAs = pd.DataFrame({"pegRNA_number": [1, 2],
                           "Name": ["pMUZ201", "pMUZ202"]})
    return plasmids, epegRNAs, ngRNAs


# --------------------------------------------------------------------------- #
# Formats
# --------------------------------------------------------------------------- #
def test_names_are_ordered_smallest_area_first():
    assert plates.names() == ["96-well", "48-well", "24-well", "12-well",
                              "6-well", "10-cm", "15-cm"]


def test_every_requested_format_is_present():
    for fmt in ["96-well", "48-well", "24-well", "12-well", "6-well",
                "10-cm", "15-cm"]:
        assert plates.spec(fmt)["plate_format"] == fmt


@pytest.mark.parametrize("spelling,expected", [
    ("96-well", "96-well"),
    ("96", "96-well"),
    ("96 well", "96-well"),
    ("96well", "96-well"),
    ("  96-WELL  ", "96-well"),
    ("6-wells", "6-well"),
    ("10-cm", "10-cm"),
    ("10cm", "10-cm"),
    ("10 cm", "10-cm"),
    ("10cm dish", "10-cm"),
    ("15 cm dish", "15-cm"),
])
def test_canonical_accepts_the_usual_spellings(spelling, expected):
    assert plates.canonical(spelling) == expected


def test_canonical_rejects_an_unknown_format():
    with pytest.raises(KeyError, match="Unknown plate format"):
        plates.canonical("384-well")


def test_canonical_rejects_a_non_string():
    with pytest.raises(TypeError, match="must be a string"):
        plates.canonical(96)


def test_spec_reports_geometry_and_seeding():
    s = plates.spec("6-well")
    assert s == {"plate_format": "6-well", "wells": 6, "area_cm2": 9.6,
                 "media_mL": 2.0, "seed_HEK293T": 450000}


def test_seeding_density_is_constant_across_formats():
    # Every format is the 96-well density, so confluency is format-independent
    for name in plates.names():
        s = plates.spec(name)
        exact = s["area_cm2"] * plates.SEED_PER_CM2
        assert s["seed_HEK293T"] == pytest.approx(exact, rel=0.02)


def test_seeding_keeps_the_labs_own_96_well_number():
    assert plates.spec("96-well")["seed_HEK293T"] == 15000


def test_scale_is_the_growth_area_ratio():
    assert plates.scale("6-well", "96-well") == pytest.approx(9.6 / 0.32)
    assert plates.scale("96-well", "96-well") == 1.0


def test_table_has_one_row_per_format():
    df = plates.table()
    assert list(df["plate_format"]) == plates.names()
    assert set(df.columns) == {"plate_format", "wells", "area_cm2",
                               "media_mL", "seed_HEK293T"}


def test_table_does_not_write_a_file_without_a_path(tmp_path):
    plates.table()
    assert list(tmp_path.iterdir()) == []


# --------------------------------------------------------------------------- #
# geometry()
# --------------------------------------------------------------------------- #
def test_geometry_without_conditions_reports_per_well_only():
    df = plates.geometry("24-well", reps=2)
    assert len(df) == 1
    assert "Wells needed" not in df.columns
    assert df.iloc[0]["Format"] == "24-well"
    assert df.iloc[0]["Replicates"] == 2


def test_geometry_with_conditions_totals_wells_media_and_cells():
    df = plates.geometry("12-well", reps=2, conditions=2)
    row = df.iloc[0]
    assert row["Wells needed"] == 4
    assert row["Plates needed"] == pytest.approx(4 / 12, abs=0.01)
    assert row["Media needed (mL)"] == pytest.approx(4 * 1.0)
    assert row["Cells needed"] == 4 * plates.spec("12-well")["seed_HEK293T"]


def test_geometry_rejects_an_unknown_format():
    with pytest.raises(KeyError, match="Unknown plate format"):
        plates.geometry("384-well")


# --------------------------------------------------------------------------- #
# PE3()
# --------------------------------------------------------------------------- #
def test_PE3_at_96_well_matches_a_direct_transfect_call(pe3_tables):
    plasmids, epegRNAs, ngRNAs = pe3_tables
    direct = transfect.PE3(plasmids.copy(), epegRNAs.copy(), ngRNAs.copy())
    scaled = plates.PE3(plasmids.copy(), epegRNAs.copy(), ngRNAs.copy(),
                        plate_format="96-well")
    for key in direct:
        pd.testing.assert_frame_equal(direct[key], scaled[key])


def test_PE3_scales_dna_by_growth_area(pe3_tables):
    plasmids, epegRNAs, ngRNAs = pe3_tables
    small = plates.PE3(plasmids.copy(), epegRNAs.copy(), ngRNAs.copy(),
                       plate_format="96-well")
    big = plates.PE3(plasmids.copy(), epegRNAs.copy(), ngRNAs.copy(),
                     plate_format="24-well")
    factor = plates.scale("24-well", "96-well")
    assert (big["Transfection"]["epegRNA uL"].iloc[0]
            == pytest.approx(small["Transfection"]["epegRNA uL"].iloc[0] * factor,
                             rel=0.01))


def test_PE3_keeps_the_plate_map_at_96_well(pe3_tables):
    plasmids, epegRNAs, ngRNAs = pe3_tables
    tables = plates.PE3(plasmids.copy(), epegRNAs.copy(), ngRNAs.copy(),
                        plate_format="96-well")
    assert "96-well plates" in tables


@pytest.mark.parametrize("fmt", ["24-well", "6-well", "10-cm"])
def test_PE3_drops_the_96_well_map_for_other_formats(pe3_tables, fmt):
    """An 8x12 grid does not describe a 24-well plate; a note is safer."""
    plasmids, epegRNAs, ngRNAs = pe3_tables
    tables = plates.PE3(plasmids.copy(), epegRNAs.copy(), ngRNAs.copy(),
                        plate_format=fmt)
    assert "96-well plates" not in tables
    assert f"{plates.LAYOUT_NOTE} ({fmt})" in tables
    assert "Transfection" in tables


def test_PE3_rejects_an_unknown_format(pe3_tables):
    plasmids, epegRNAs, ngRNAs = pe3_tables
    with pytest.raises(KeyError, match="Unknown plate format"):
        plates.PE3(plasmids, epegRNAs, ngRNAs, plate_format="384-well")


# --------------------------------------------------------------------------- #
# virus()
# --------------------------------------------------------------------------- #
def test_virus_at_6_well_matches_a_direct_transfect_call(packaging):
    direct = transfect.virus(packaging.copy())
    scaled = plates.virus(packaging.copy(), plate_format="6-well")
    pd.testing.assert_frame_equal(direct, scaled)


def test_virus_scales_dna_by_growth_area(packaging):
    six = plates.virus(packaging.copy(), plate_format="6-well")
    twelve = plates.virus(packaging.copy(), plate_format="12-well")
    factor = plates.scale("12-well", "6-well")
    vsvg_6 = six[six["Description"] == "VSVG"]["ng"].iloc[0]
    vsvg_12 = twelve[twelve["Description"] == "VSVG"]["ng"].iloc[0]
    assert vsvg_12 == pytest.approx(vsvg_6 * factor, rel=0.01)


def test_virus_scales_the_lipofectamine_volumes(packaging):
    """P3000/L3000 were fixed per-well constants; scaling them is the point."""
    six = plates.virus(packaging.copy(), plate_format="6-well")
    ten = plates.virus(packaging.copy(), plate_format="10-cm")
    factor = plates.scale("10-cm", "6-well")
    for reagent in ("P3000 Enhancer", "L3000"):
        small = six[six["Description"] == reagent]["uL"].iloc[0]
        big = ten[ten["Description"] == reagent]["uL"].iloc[0]
        assert big == pytest.approx(small * factor, rel=0.01)


def test_virus_does_not_mutate_the_input_frame(packaging):
    before = packaging.copy()
    plates.virus(packaging, plate_format="24-well")
    pd.testing.assert_frame_equal(packaging, before)


def test_virus_rejects_an_unknown_format(packaging):
    with pytest.raises(KeyError, match="Unknown plate format"):
        plates.virus(packaging, plate_format="384-well")
