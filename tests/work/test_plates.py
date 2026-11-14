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
                              "6-well", "6-cm", "10-cm", "15-cm"]


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
        # Seeding numbers are rounded to two significant figures so they are
        # worth pipetting, which can move one by up to ~5%
        assert s["seed_HEK293T"] == pytest.approx(exact, rel=0.05)


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


# Thermo Fisher's Lipofectamine 3000 lentiviral protocol, Tables 1 and 2:
# per-vessel P3000, Lipofectamine 3000 and total DNA. Their own scale-up is
# growth-area proportional off the 6-well anchor, which is what plates.virus
# does, so these are an outside check on the scaling rather than a
# restatement of it.
THERMO_LENTI = [
    # vessel, P3000 uL, L3000 uL, total DNA ug
    ("6-well", 6, 7, 3.0),
    ("6-cm", 13, 15, 6.7),
    ("10-cm", 35, 41, 17.3),
]


@pytest.mark.parametrize("vessel,p3000,l3000,dna_ug", THERMO_LENTI)
def test_virus_scaling_matches_thermos_own_table(vessel, p3000, l3000, dna_ug):
    """Growth-area scaling off the 6-well lands on Thermo's published rows."""
    factor = plates.scale(vessel, "6-well")
    assert 6 * factor == pytest.approx(p3000, rel=0.05)
    assert 7 * factor == pytest.approx(l3000, rel=0.05)
    assert 3.0 * factor == pytest.approx(dna_ug, rel=0.05)


def test_virus_plan_scales_by_the_same_factor_thermo_uses(packaging):
    """Whatever the plan's absolute volumes are, the ratio is the area ratio."""
    six = plates.virus(packaging.copy(), plate_format="6-well")
    ten = plates.virus(packaging.copy(), plate_format="10-cm")
    factor = plates.scale("10-cm", "6-well")
    for reagent in ("P3000 Enhancer", "L3000"):
        small = six[six["Description"] == reagent]["uL"].iloc[0]
        big = ten[ten["Description"] == reagent]["uL"].iloc[0]
        assert big / small == pytest.approx(factor, rel=0.01)
    # and Thermo scales 35/6 and 41/7 the same way
    assert 35 / 6 == pytest.approx(factor, rel=0.05)
    assert 41 / 7 == pytest.approx(factor, rel=0.05)


def test_virus_total_dna_matches_thermos_three_micrograms():
    """750 + 1500 + 750 ng is Thermo's 3 ug at 3:1 packaging to transfer."""
    from edms.bio import transfect
    import inspect

    defaults = inspect.signature(transfect.virus).parameters
    vsvg = defaults["VSVG_ng"].default
    gagpol = defaults["GagPol_ng"].default
    transfer = defaults["transfer_ng"].default
    assert vsvg + gagpol + transfer == 3000
    assert (vsvg + gagpol) / transfer == 3.0
    assert defaults["P3000_uL_per_well"].default == 6
    assert defaults["L3000_uL_per_well"].default == 7
    assert defaults["well_uL"].default == 500


def test_virus_seeding_density_is_thermos_lentiviral_number():
    """1.2e6 per 6-well is what the protocol specifies, not the PE3 density."""
    plan = plates.geometry("6-well", seed_per_cm2=plates.VIRUS_SEED_PER_CM2)
    assert plan["Seed HEK293T per well"].iloc[0] == 1200000
    assert plates.VIRUS_SEED_PER_CM2 > plates.SEED_PER_CM2


def test_geometry_defaults_to_the_transfection_density():
    assert (plates.geometry("6-well")["Seed HEK293T per well"].iloc[0]
            != plates.geometry("6-well",
                               seed_per_cm2=plates.VIRUS_SEED_PER_CM2
                               )["Seed HEK293T per well"].iloc[0])


def test_vessel_areas_are_thermos_not_cornings():
    """The bench uses Thermo plates; MUZ361 names them. The two differ."""
    assert plates.spec("48-well")["area_cm2"] == 1.1     # Corning: 0.95
    assert plates.spec("12-well")["area_cm2"] == 3.5     # Corning: 3.8
    assert plates.spec("15-cm")["area_cm2"] == 145.0     # Corning: 152


@pytest.mark.parametrize("spelling,expected", [
    ("60mm", "6-cm"), ("100 mm", "10-cm"), ("150mm", "15-cm"),
    ("6", "6-well"), ("6cm", "6-cm"),
])
def test_canonical_resolves_thermos_millimetre_spellings(spelling, expected):
    """'6' is a 6-well plate to anyone who types it, not a 60 mm dish."""
    assert plates.canonical(spelling) == expected


def test_virus_does_not_mutate_the_input_frame(packaging):
    before = packaging.copy()
    plates.virus(packaging, plate_format="24-well")
    pd.testing.assert_frame_equal(packaging, before)


def test_virus_rejects_an_unknown_format(packaging):
    with pytest.raises(KeyError, match="Unknown plate format"):
        plates.virus(packaging, plate_format="384-well")


# --------------------------------------------------------------------------- #
# expansion()
#
# MUZ361 expanded sorted clones from 96-well up to a 15 cm dish, and recorded
# why it skipped 12-well: cells clumped into the middle of the well and died
# there twice, so it went 24-well straight to 6-well.
# --------------------------------------------------------------------------- #
def test_expansion_skips_12_well_by_default():
    ladder = plates.expansion("96-well", "15-cm")
    assert "12-well" not in set(ladder["From"]) | set(ladder["To"])
    step = ladder[ladder["From"] == "24-well"].iloc[0]
    assert step["To"] == "6-well"
    assert "skipping 12-well" in step["Note"]


def test_expansion_can_put_12_well_back():
    ladder = plates.expansion("96-well", "15-cm", skip=())
    assert "12-well" in set(ladder["To"])
    assert list(ladder[ladder["From"] == "24-well"]["To"]) == ["12-well"]


def test_expansion_split_is_the_growth_area_ratio():
    ladder = plates.expansion("96-well", "24-well")
    for _, row in ladder.iterrows():
        assert row["Area ratio"] == round(
            plates.scale(to=row["To"], frm=row["From"]), 2)


def test_expansion_flags_a_wide_split():
    """Skipping a rung widens the step either side of it."""
    narrow = plates.expansion("96-well", "6-well", skip=())
    wide = plates.expansion("96-well", "6-well")
    assert all(not n for n in narrow["Note"])
    assert any("wide split" in n for n in wide["Note"])


def test_expansion_covers_every_rung_between_the_ends():
    ladder = plates.expansion("96-well", "15-cm", skip=())
    assert list(ladder["From"]) == ["96-well", "48-well", "24-well",
                                    "12-well", "6-well", "6-cm", "10-cm"]
    assert ladder["To"].iloc[-1] == "15-cm"


def test_expansion_accepts_aliases():
    assert plates.expansion("96", "15cm").equals(
        plates.expansion("96-well", "15-cm"))


def test_expansion_refuses_to_go_down():
    with pytest.raises(ValueError, match="not larger than"):
        plates.expansion("15-cm", "96-well")


def test_expansion_refuses_a_zero_length_ladder():
    with pytest.raises(ValueError, match="not larger than"):
        plates.expansion("6-well", "6-well")


def test_expansion_refuses_to_skip_its_own_endpoint():
    with pytest.raises(ValueError, match="cannot be left out of its own ladder"):
        plates.expansion("12-well", "15-cm")


def test_expansion_writes_a_file(tmp_path):
    out = tmp_path / "ladder.csv"
    plates.expansion("96-well", "6-well", file=str(out))
    assert out.exists()
