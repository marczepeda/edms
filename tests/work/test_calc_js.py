"""
Tests for edms/work/static/calc.js

calc.js is a second implementation of calculators that already exist in
Python, which is the one genuinely dangerous thing about the interactive
protocols: a browser number that has drifted from the package looks exactly
as authoritative as one that has not.

So every ported target is run under node and compared to the Python, cell for
cell, over the values the bundled templates actually use and over randomised
ones. A port that drifts fails here rather than at the bench.

Skipped, not passed, when node is unavailable.
"""
import json
import random
import shutil
import subprocess

import pandas as pd
import pytest

from edms.work import interactive, plates

NODE = shutil.which("node")
pytestmark = pytest.mark.skipif(NODE is None,
                                reason="node is needed to run calc.js")

DRIVER = """
const path = process.argv[2];
const EDMS = require(path);
EDMS.consts = JSON.parse(process.argv[3]);
const target = process.argv[4];
const args = JSON.parse(process.argv[5]);
let out;
try {
  out = {ok: true, value: EDMS.calc[target](args)};
} catch (e) {
  out = {ok: false, error: e.message};
}
process.stdout.write(JSON.stringify(out));
"""


@pytest.fixture(scope="module")
def driver(tmp_path_factory):
    p = tmp_path_factory.mktemp("js") / "driver.js"
    p.write_text(DRIVER, encoding="utf-8")
    return str(p)


@pytest.fixture(scope="module")
def consts():
    return json.dumps(interactive.constants())


def run_js(driver, consts, target, args):
    """Call one calc.js target under node and return its table.

    __floats__ is what resolveArgs() in app.js attaches: the arguments Python
    would print with a trailing .0, which JSON cannot carry on its own. Adding
    it here means the parity is between the Python and what the page actually
    calls, rather than between the Python and a stripped-down version of it.
    """
    args = dict(args)
    args["__floats__"] = sorted(k for k, v in args.items()
                                if isinstance(v, float))
    calc = str(interactive.STATIC / "calc.js")
    proc = subprocess.run(
        [NODE, driver, calc, consts, target, json.dumps(args)],
        capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, proc.stderr
    out = json.loads(proc.stdout)
    assert out["ok"], out.get("error")
    return out["value"]


def same(js, df, where=""):
    """Compare a calc.js table to a DataFrame, cell for cell."""
    want = interactive.table_data(df)
    assert js["columns"] == want["columns"], f"columns differ {where}"
    assert len(js["rows"]) == len(want["rows"]), f"row count differs {where}"
    for i, (got, exp) in enumerate(zip(js["rows"], want["rows"])):
        for col in want["columns"]:
            g, e = got.get(col), exp.get(col)
            if isinstance(e, float) and isinstance(g, (int, float)):
                assert g == pytest.approx(e, rel=1e-9, abs=1e-9), \
                    f"{where} row {i} col {col}: {g} != {e}"
            else:
                assert str(g) == str(e), f"{where} row {i} col {col}: {g!r} != {e!r}"


# --------------------------------------------------------------------------- #
# The file itself
# --------------------------------------------------------------------------- #
def test_calc_js_parses():
    proc = subprocess.run([NODE, "--check", str(interactive.STATIC / "calc.js")],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr


def test_live_targets_are_all_callable(driver, consts):
    """live_targets() is read out of the file, so it cannot name a target the
    file does not actually define."""
    for target in interactive.live_targets():
        proc = subprocess.run(
            [NODE, "-e",
             f'const E=require({str(interactive.STATIC / "calc.js")!r});'
             f'process.exit(typeof E.calc[{target!r}] === "function" ? 0 : 1)'],
            capture_output=True, text=True)
        assert proc.returncode == 0, f"{target} is not defined in calc.js"


def test_every_live_target_exists_in_python():
    """A target calc.js implements that Python does not would be a number with
    no authority behind it."""
    import importlib

    for target in interactive.live_targets():
        module, _, func = target.partition(":")
        assert callable(getattr(importlib.import_module(module), func)), target


# --------------------------------------------------------------------------- #
# Rounding
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("value,digits", [
    (0.125, 2), (0.135, 2), (2.675, 2), (1.005, 2), (0.5, 0), (1.5, 0),
    (2.5, 0), (-0.125, 2), (16.666666, 2), (0.2333333, 4), (1234.5, 0),
])
def test_pyround_matches_python(driver, consts, value, digits):
    """Python rounds half to even; JavaScript's Math.round does not."""
    calc = str(interactive.STATIC / "calc.js")
    proc = subprocess.run(
        [NODE, "-e",
         f"const E=require({calc!r});"
         f"process.stdout.write(String(E.calc.pyRound({value},{digits})))"],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert float(proc.stdout) == round(value, digits)


# --------------------------------------------------------------------------- #
# geometry
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("fmt", plates.names())
@pytest.mark.parametrize("density", [None, plates.VIRUS_SEED_PER_CM2])
def test_geometry_matches(driver, consts, fmt, density):
    args = {"plate_format": fmt, "reps": 2, "conditions": 3}
    if density is not None:
        args["seed_per_cm2"] = density
    js = run_js(driver, consts, "edms.work.plates:geometry", args)
    same(js, plates.geometry(**args), where=f"geometry {fmt}")


@pytest.mark.parametrize("fmt", plates.names())
def test_geometry_without_conditions_matches(driver, consts, fmt):
    args = {"plate_format": fmt, "reps": 1}
    js = run_js(driver, consts, "edms.work.plates:geometry", args)
    same(js, plates.geometry(**args), where=f"geometry bare {fmt}")


# --------------------------------------------------------------------------- #
# virus
# --------------------------------------------------------------------------- #
def plasmid_frame(concs):
    names = ["pMUZ26.6", "pMUZ26.7"] + [f"pMUZ3{i}.1" for i in range(len(concs) - 2)]
    return pd.DataFrame({
        "Plasmid": names,
        "Description": ["VSVG", "GagPol"] + ["transfer"] * (len(concs) - 2),
        "Colony": ["c1"] * len(concs),
        "ng/uL": list(concs)})


@pytest.mark.parametrize("fmt", plates.names())
def test_virus_matches(driver, consts, fmt):
    df = plasmid_frame([412.0, 388.5, 301.2, 520.4])
    args = {"plasmids": interactive.table_data(df), "plate_format": fmt,
            "reps": 2, "mm_x": 1.1}
    js = run_js(driver, consts, "edms.work.plates:virus", args)
    same(js, plates.virus(df.copy(), plate_format=fmt, reps=2, mm_x=1.1),
         where=f"virus {fmt}")


def test_virus_matches_on_random_inputs(driver, consts):
    """Rounding is where two implementations drift, so push odd values at it."""
    rng = random.Random(20261114)
    for _ in range(12):
        n = rng.randint(3, 6)
        concs = [round(rng.uniform(8.0, 1400.0), 1) for _ in range(n)]
        df = plasmid_frame(concs)
        fmt = rng.choice(plates.names())
        reps = rng.randint(1, 4)
        mm_x = rng.choice([1.0, 1.05, 1.1, 1.2, 1.35])
        args = {"plasmids": interactive.table_data(df), "plate_format": fmt,
                "reps": reps, "mm_x": mm_x}
        js = run_js(driver, consts, "edms.work.plates:virus", args)
        same(js, plates.virus(df.copy(), plate_format=fmt, reps=reps, mm_x=mm_x),
             where=f"virus {fmt} reps={reps} mm_x={mm_x} concs={concs}")


@pytest.mark.parametrize("fmt", plates.names())
def test_virus_per_well_matches(driver, consts, fmt):
    js = run_js(driver, consts, "edms.work.plates:virus_per_well",
                {"plate_format": fmt})
    same(js, plates.virus_per_well(fmt), where=f"per_well {fmt}")


@pytest.mark.parametrize("fmt", ["96-well", "48-well", "24-well", "12-well",
                                 "6-well"])
def test_virus_dilutions_matches(driver, consts, fmt):
    df = plasmid_frame([412.0, 388.5, 301.2, 520.4, 188.7])
    args = {"plasmids": interactive.table_data(df), "plate_format": fmt,
            "reps": 2, "mm_x": 1.1}
    js = run_js(driver, consts, "edms.work.plates:virus_dilutions", args)
    same(js, plates.virus_dilutions(df.copy(), plate_format=fmt, reps=2,
                                    mm_x=1.1),
         where=f"dilutions {fmt}")


def test_virus_dilutions_matches_on_random_inputs(driver, consts):
    rng = random.Random(778)
    for _ in range(8):
        concs = [round(rng.uniform(15.0, 900.0), 2) for _ in range(5)]
        df = plasmid_frame(concs)
        fmt = rng.choice(["96-well", "48-well", "24-well", "12-well"])
        args = {"plasmids": interactive.table_data(df), "plate_format": fmt,
                "reps": rng.randint(1, 3), "mm_x": 1.1}
        js = run_js(driver, consts, "edms.work.plates:virus_dilutions", args)
        same(js, plates.virus_dilutions(df.copy(), plate_format=fmt,
                                        reps=args["reps"], mm_x=1.1),
             where=f"dilutions {fmt} {concs}")


# --------------------------------------------------------------------------- #
# expansion
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("start,end", [
    ("96-well", "6-well"), ("96-well", "15-cm"), ("24-well", "10-cm"),
    ("48-well", "6-cm"), ("6-well", "15-cm"),
])
def test_expansion_matches(driver, consts, start, end):
    js = run_js(driver, consts, "edms.work.plates:expansion",
                {"start": start, "end": end})
    same(js, plates.expansion(start, end), where=f"expansion {start}->{end}")


def test_expansion_matches_with_nothing_skipped(driver, consts):
    """12-well is in the default skip, so the ladder that keeps it is the one
    that exercises every rung."""
    js = run_js(driver, consts, "edms.work.plates:expansion",
                {"start": "96-well", "end": "15-cm", "skip": []})
    same(js, plates.expansion("96-well", "15-cm", skip=[]),
         where="expansion full ladder")


# --------------------------------------------------------------------------- #
# Multi-table calculators
# --------------------------------------------------------------------------- #
def same_tables(js, result, where=""):
    """Compare a calc.js result to whatever the Python returned, after
    edms.work.plan has normalised it into named tables. The names matter as
    much as the numbers: a rendered protocol looks a table up by the name the
    Python gave it, so a port that computes the right figures under the wrong
    key paints nothing."""
    from edms.work import plan

    want = (plan._tables(result, ["result"])
            if not isinstance(result, pd.DataFrame) else {"_": result})
    got = js.get("__tables__") if isinstance(js, dict) else None
    if got is None:
        got = {"_": js}
    assert sorted(got) == sorted(want), f"table names differ {where}"
    for name in want:
        same(got[name], want[name], where=f"{where}/{name}")


def frame_data(df):
    """A DataFrame as the page would receive it, NaN included."""
    return interactive.table_data(df)


# --------------------------------------------------------------------------- #
# Number formatting primitives
# --------------------------------------------------------------------------- #
def js_eval(expr):
    calc = str(interactive.STATIC / "calc.js")
    proc = subprocess.run(
        [NODE, "-e", f"const E=require({calc!r});"
                     f"process.stdout.write(String({expr}))"],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return proc.stdout


@pytest.mark.parametrize("value", [
    0, 1, 100, 100.0, 0.5, 1.5, 0.1, 1e-5, 6e-6, 123456.0, 1234567.0,
    0.000123456, 99999.5, 1e6, -42.5, 0.75, 20.0, 1.19,
])
def test_pyg_matches_python_g_format(value):
    """The ':g' in a note cell: '100 ng at 412 ng/uL', never '100.0 ng'."""
    assert js_eval(f"E.calc.pyG({value!r})") == format(value, "g")


@pytest.mark.parametrize("vals", [
    ["P5-FWD", "P7-REV"], ["F1", "R1"], ["a'b", "c"], ["plain"],
    ["with space", "x"],
])
def test_pytuple_matches_python_repr(vals):
    """plan._tables names a master mix after str((fwd, rev)), so the page has
    to spell a Python tuple exactly."""
    assert js_eval(f"E.calc.pyTuple({json.dumps(vals)})") == str(tuple(vals))


# --------------------------------------------------------------------------- #
# clone
# --------------------------------------------------------------------------- #
from edms.bio import clone, culture, ngs  # noqa: E402


PCR_COMPONENTS = pd.DataFrame([
    {"Component": "5x Q5 buffer", "Stock": 5, "Desired": 1, "Unit": "x"},
    {"Component": "dNTPs", "Stock": 10, "Desired": 0.2, "Unit": "mM"},
    {"Component": "FWD primer", "Stock": 10, "Desired": 0.5, "Unit": "uM"},
    {"Component": "REV primer", "Stock": 10, "Desired": 0.5, "Unit": "uM"},
    {"Component": "Template", "uL": 1.0, "In MM": False},
    {"Component": "Q5 polymerase", "Stock": 2, "Desired": 0.02,
     "Unit": "U/uL", "Note": "add last"},
])


@pytest.mark.parametrize("rxn_uL,samples,mm_x", [
    (50, 3, 1.1), (20, 1, 1.05), (25, 8, 1.2), (10, 2, 1.0),
])
def test_rxn_mm_matches(driver, consts, rxn_uL, samples, mm_x):
    js = run_js(driver, consts, "edms.bio.clone:rxn_mm",
                {"df": frame_data(PCR_COMPONENTS), "rxn_uL": rxn_uL,
                 "samples": samples, "mm_x": mm_x})
    same(js, clone.rxn_mm(PCR_COMPONENTS, rxn_uL=rxn_uL, samples=samples,
                          mm_x=mm_x),
         where=f"rxn_mm {rxn_uL}/{samples}/{mm_x}")


def test_rxn_mm_drops_the_note_column_the_same_way(driver, consts):
    """The Note column is only there when a component has one, and the page
    has to agree about that or the header row shifts."""
    bare = PCR_COMPONENTS.drop(columns="Note")
    js = run_js(driver, consts, "edms.bio.clone:rxn_mm",
                {"df": frame_data(bare), "rxn_uL": 50, "samples": 2})
    assert "Note" not in js["columns"]
    same(js, clone.rxn_mm(bare, rxn_uL=50, samples=2), where="rxn_mm bare")


def test_rxn_mm_overfill_is_an_error_on_both_sides(driver, consts):
    """Water is the balance, so too much of anything is an error rather than
    a silently overfilled tube -- in the page as well."""
    calc = str(interactive.STATIC / "calc.js")
    # The template is a fixed 1 uL, so a 1 uL reaction cannot hold it and
    # the 0.33 uL of everything else beside it
    args = {"df": frame_data(PCR_COMPONENTS), "rxn_uL": 1, "samples": 1}
    proc = subprocess.run(
        [NODE, driver, calc, consts, "edms.bio.clone:rxn_mm", json.dumps(args)],
        capture_output=True, text=True, timeout=60)
    assert json.loads(proc.stdout)["ok"] is False
    with pytest.raises(ValueError):
        clone.rxn_mm(PCR_COMPONENTS, rxn_uL=1, samples=1)


@pytest.mark.parametrize("kwargs", [
    {}, {"oligos": 6}, {"oligos": 3, "rxn_uL": 20, "mm_x": 1.1},
    {"oligos": 12, "pnk_uL": 1.0, "buffer_stock": 10, "buffer_desired": 2},
])
def test_anneal_mm_matches(driver, consts, kwargs):
    js = run_js(driver, consts, "edms.bio.clone:anneal_mm", kwargs)
    same_tables(js, clone.anneal_mm(**kwargs), where=f"anneal_mm {kwargs}")


@pytest.mark.parametrize("kwargs", [
    {}, {"inserts": 1, "samples": 4, "insert_uL": 2.0},
    {"inserts": 5, "samples": 2, "shared_inserts": True, "cycles": 30},
    {"inserts": 2, "vector_ng": 150, "vector_stock": 412.5},
    {"inserts": 4, "samples": 3, "rxn_uL": 20, "mm_x": 1.2,
     "enzyme": "FastDigest BsmBI", "vector": "pMUZ1"},
])
def test_golden_gate_matches(driver, consts, kwargs):
    js = run_js(driver, consts, "edms.bio.clone:golden_gate", kwargs)
    same_tables(js, clone.golden_gate(**kwargs), where=f"golden_gate {kwargs}")


FRAGMENTS = pd.DataFrame([
    {"Fragment": "backbone", "kb": 5.2, "ng/uL": 88.4},
    {"Fragment": "insert A", "kb": 0.35, "ng/uL": 22.1},
    {"Fragment": "insert B", "kb": 1.8, "ng/uL": 140.0},
])


@pytest.mark.parametrize("kwargs", [
    {}, {"rxn_uL": 20, "mm_uL": 10},
    {"pmol": 0.03, "pmol_small": 0.08, "small_kb": 2.0},
])
def test_gibson_matches(driver, consts, kwargs):
    js = run_js(driver, consts, "edms.bio.clone:gibson",
                dict(df=frame_data(FRAGMENTS), **kwargs))
    same(js, clone.gibson(FRAGMENTS, **kwargs), where=f"gibson {kwargs}")


def test_gibson_matches_on_random_inputs(driver, consts):
    rng = random.Random(4242)
    for _ in range(10):
        n = rng.randint(2, 5)
        df = pd.DataFrame([
            {"Fragment": f"f{i}", "kb": round(rng.uniform(0.15, 8.0), 2),
             "ng/uL": round(rng.uniform(12.0, 400.0), 1)} for i in range(n)])
        kwargs = {"rxn_uL": rng.choice([10, 15, 20]),
                  "pmol": rng.choice([0.02, 0.025, 0.04])}
        js = run_js(driver, consts, "edms.bio.clone:gibson",
                    dict(df=frame_data(df), **kwargs))
        same(js, clone.gibson(df, **kwargs), where=f"gibson random {kwargs}")


# --------------------------------------------------------------------------- #
# culture
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("kwargs", [
    {"stock": 1000, "desired": 1, "wells": 24},
    {"stock": 10000, "desired": 8, "wells": 6, "drug_media_mL": 0.5,
     "cell_media_mL": 0.25, "mm_x": 1.2, "drug": "Polybrene"},
    {"stock": 50, "desired": 0.5, "wells": 96, "drug_media_mL": 0.1,
     "cell_media_mL": 0.1, "mm_x": 1.05, "drug": "DHT", "unit": "nM"},
])
def test_drug_mm_matches(driver, consts, kwargs):
    js = run_js(driver, consts, "edms.bio.culture:drug_mm", kwargs)
    same(js, culture.drug_mm(**kwargs), where=f"drug_mm {kwargs}")


@pytest.mark.parametrize("kwargs", [
    {"wells": 12, "flow_per_mL": 8.4e5, "cells_per_well": 1.8e5},
    {"wells": 3, "flow_per_mL": 1.2e6, "cells_per_well": 1.2e5,
     "virus_mL_per_well": 0.5, "cell_media_mL": 0.25, "polybrene_desired": 10,
     "mm_x_polybrene": 1.15, "mm_x_cells": 1.2},
])
def test_transduction_matches(driver, consts, kwargs):
    js = run_js(driver, consts, "edms.bio.culture:transduction", kwargs)
    same_tables(js, culture.transduction(**kwargs),
                where=f"transduction {kwargs}")


@pytest.mark.parametrize("kwargs", [
    {"target_clones": 24},
    {"target_clones": 96, "plating_efficiency": 0.05, "plate_format": "48-well"},
    {"target_clones": 10, "plating_efficiency": 0.3,
     "conditioned_fraction": 0.25, "media_mL": 0.15},
    {"target_clones": 1, "plating_efficiency": 1.0, "conditioned_fraction": 0},
])
def test_sort_plan_matches(driver, consts, kwargs):
    js = run_js(driver, consts, "edms.bio.culture:sort_plan", kwargs)
    same(js, culture.sort_plan(**kwargs), where=f"sort_plan {kwargs}")


@pytest.mark.parametrize("kwargs", [
    {"harvest_mL": 10, "vessel": "15-cm", "confluency": 80},
    {"harvest_mL": 5, "vessel": "10-cm", "confluency": 75, "recovery": 0.75},
    {"harvest_mL": 4, "cells_per_mL": 3.2e6, "cells_per_vial": 2e6,
     "dmso_pct": 5, "vial_mL": 1.5},
])
def test_freeze_down_matches(driver, consts, kwargs):
    js = run_js(driver, consts, "edms.bio.culture:freeze_down", kwargs)
    same(js, culture.freeze_down(**kwargs), where=f"freeze_down {kwargs}")


SAMPLES = pd.DataFrame([
    {"Sample": "A", "Vessel": "15-cm", "Confluency %": 80},
    {"Sample": "B", "Vessel": "6-well", "Confluency %": 60},
    {"Sample": "C", "Vessel": "15-cm", "Confluency %": 100,
     "Cells, flow": 1.4e6},
])


@pytest.mark.parametrize("kwargs", [
    {}, {"transduce_wells": 0},
    {"cg_wells": 6, "extra_samples": 2, "competitor": "LCG-GFP",
     "competitor_vessel": "6-well", "competitor_confluency": 85},
])
def test_competitive_growth_matches(driver, consts, kwargs):
    """Two arms off one harvest, including the 'Sum spin, mL' column that only
    appears once there is more than one."""
    js = run_js(driver, consts, "edms.bio.culture:competitive_growth",
                dict(df=frame_data(SAMPLES), **kwargs))
    same_tables(js, culture.competitive_growth(SAMPLES, **kwargs),
                where=f"competitive_growth {kwargs}")


# --------------------------------------------------------------------------- #
# ngs
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("kwargs", [
    {"genotypes": 1000},
    {"genotypes": 4, "samples": 6, "cell_coverage": 500, "umi_coverage": 10},
    {"genotypes": 12000, "samples": 2, "ng_per_pcr": 150, "ploidy_per_cell": 4},
])
def test_coverage_matches(driver, consts, kwargs):
    js = run_js(driver, consts, "edms.bio.ngs:coverage", kwargs)
    same(js, ngs.coverage(**kwargs), where=f"coverage {kwargs}")


def ngs_sheet(n=6, umi=None):
    rows = []
    for i in range(n):
        row = {"ID": f"S{i + 1}", "PCR1 ID": f"P1_{i + 1}",
               "PCR1 FWD": f"F{i % 2 + 1}", "PCR1 REV": f"R{i % 2 + 1}",
               "PCR1 Tm": 65 + (i % 2), "PCR1 bp": 420 + 100 * (i % 3),
               "PCR2 ID": f"P2_{i + 1}", "PCR2 FWD": "PCR2-FWD",
               "PCR2 REV": f"idx{i + 1}", "PCR2 Tm": 60, "PCR2 bp": 300}
        if umi is not None:
            row["UMI"] = umi[i]
        rows.append(row)
    return pd.DataFrame(rows)


@pytest.mark.parametrize("label,n,umi,kwargs", [
    ("plain", 6, None, {}),
    ("no pcr2", 6, None, {"pcr2": False}),
    ("unsplit", 6, None, {"split_pcr1_primers": False}),
    ("ultra", 6, None, {"ultra": True}),
    ("syber", 6, None, {"syber": True}),
    ("split wells", 6, None, {"pcr1_total_uL": 120, "pcr2_total_uL": 60}),
    ("cycles", 6, None, {"pcr1_cycles": 35, "pcr2_cycles": 10}),
    ("inner wells", 14, None, {"inner1": True, "inner2": True}),
    ("two plates", 40, None, {}),
    ("umi", 6, [True, True, False, False, True, False], {}),
    ("umi ultra", 6, [True, False, True, False, True, False], {"ultra": True}),
])
def test_pcrs_matches(driver, consts, label, n, umi, kwargs):
    """The whole NGS plan: plate layouts, master mixes and thermocycler
    programs, named the way edms.work.plan names them."""
    df = ngs_sheet(n, umi)
    js = run_js(driver, consts, "edms.bio.ngs:pcrs",
                dict(df=frame_data(df), **kwargs))
    same_tables(js, ngs.pcrs(df.copy(), **kwargs), where=f"pcrs {label}")


# --------------------------------------------------------------------------- #
# PE3
# --------------------------------------------------------------------------- #
PE3_PLASMIDS = pd.DataFrame({
    "Plasmid": ["pMUZ86.7", "pMUZ100.1", "pMUZ100.2", "pMUZ100.3",
                "pMUZ200.1", "pMUZ200.2", "pMUZ200.3"],
    "Description": ["PE", "epeg 1", "epeg 2", "epeg 3", "ng 1", "ng 2", "ng 3"],
    "Colony": ["c1", "c2", "c3", "c4", "c5", "c6", "c7"],
    "ng/uL": [512.0, 288.1, 301.4, 150.5, 190.2, 265.0, 410.8]})
PE3_EPEG = pd.DataFrame({"pegRNA_number": [1, 2, 3],
                         "Name": ["pMUZ100.1", "pMUZ100.2", "pMUZ100.3"]})
PE3_NG = pd.DataFrame({"pegRNA_number": [1, 2, 3, 1],
                       "Name": ["pMUZ200.1", "pMUZ200.2", "pMUZ200.3",
                                "pMUZ200.2"]})


@pytest.mark.parametrize("kwargs", [
    {}, {"plate_format": "24-well"}, {"reps": 2, "mm_x": 1.2, "well_uL": 12},
    {"plate_format": "48-well", "reps": 1, "epegRNA_ng": 80, "ngRNA_ng": 30,
     "PE_ng": 250},
    {"plate_format": "6-well"},
])
def test_pe3_matches(driver, consts, kwargs):
    """Away from 96-well the plate map is replaced by a note, because the
    8x12 grid stops describing the plate."""
    args = dict(plasmids=frame_data(PE3_PLASMIDS),
                epegRNAs=frame_data(PE3_EPEG), ngRNAs=frame_data(PE3_NG),
                **kwargs)
    js = run_js(driver, consts, "edms.work.plates:PE3", args)
    same_tables(js, plates.PE3(plasmids=PE3_PLASMIDS, epegRNAs=PE3_EPEG,
                               ngRNAs=PE3_NG, **kwargs),
                where=f"PE3 {kwargs}")


def test_pe3_wraps_onto_a_second_plate_the_same_way(driver, consts):
    epeg = pd.DataFrame({"pegRNA_number": list(range(1, 16)),
                         "Name": [f"pMUZ300.{i}" for i in range(1, 16)]})
    ngr = pd.DataFrame({"pegRNA_number": list(range(1, 16)),
                        "Name": [f"pMUZ400.{i}" for i in range(1, 16)]})
    plasmids = pd.DataFrame({
        "Plasmid": ["pMUZ86.7"] + list(epeg["Name"]) + list(ngr["Name"]),
        "Description": ["PE"] + ["e"] * 15 + ["n"] * 15,
        "Colony": ["c"] * 31,
        "ng/uL": [512.0] + [200.0 + i for i in range(30)]})
    for reps in (2, 3):
        args = dict(plasmids=frame_data(plasmids), epegRNAs=frame_data(epeg),
                    ngRNAs=frame_data(ngr), reps=reps)
        js = run_js(driver, consts, "edms.work.plates:PE3", args)
        same_tables(js, plates.PE3(plasmids=plasmids, epegRNAs=epeg,
                                   ngRNAs=ngr, reps=reps),
                    where=f"PE3 wrap reps={reps}")


# --------------------------------------------------------------------------- #
# Nothing is left behind
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# plates:layout
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("items,fmt,per_item,skip_outer", [
    (352, "96-well", 1, False),      # the ngRNA oligo set, four plates
    (352, "48-well", 2, False),      # two colonies each, fifteen plates
    (96, "96-well", 1, False),       # exactly full: four plates must not be five
    (97, "96-well", 1, False),       # one over, so the off-by-one shows
    (0, "96-well", 1, False),        # nothing to place
    (1, "6-well", 1, False),
    (60, "96-well", 1, True),        # exactly fills the inner wells
    (61, "96-well", 1, True),
    (10, "24-well", 3, False),
])
def test_layout_matches_python(driver, consts, items, fmt, per_item, skip_outer):
    args = {"items": items, "plate_format": fmt, "per_item": per_item,
            "skip_outer": skip_outer}
    js = run_js(driver, consts, "edms.work.plates:layout", args)
    same_tables(js, plates.layout(**args), where=f"{items} in {fmt}")


def test_layout_of_a_real_table_matches_python(driver, consts):
    """The oligos arrive as an input table, not a count, and the page has to
    carry their own columns through into the map."""
    items = pd.DataFrame({
        "Oligo_ID": [f"nMUZ380_{i}" for i in range(1, 15)],
        "Edit": ["F266L"] * 7 + ["R233H"] * 7,
    })
    args = {"plate_format": "6-well", "per_item": 1}
    js = run_js(driver, consts,
                "edms.work.plates:layout",
                dict(args, items=interactive.table_data(items)))
    same_tables(js, plates.layout(items, **args), where="oligo table")


def test_layout_breaks_plates_on_a_group_the_same_way(driver, consts):
    items = pd.DataFrame({"Pool": ["p1"] * 5 + ["p2"] * 5 + ["p3"] * 2})
    args = {"plate_format": "6-well", "by": "Pool"}
    js = run_js(driver, consts, "edms.work.plates:layout",
                dict(args, items=interactive.table_data(items)))
    same_tables(js, plates.layout(items, **args), where="grouped")


@pytest.mark.parametrize("args", [
    {"items": 4, "plate_format": "10-cm"},          # a dish has no grid
    {"items": 4, "plate_format": "6-well", "skip_outer": True},
    {"items": 4, "plate_format": "96-well", "per_item": 0},
])
def test_layout_fails_on_both_sides(driver, consts, args):
    """An impossible layout is an error in the page too, not a blank table."""
    with pytest.raises((KeyError, ValueError)):
        plates.layout(**args)
    calc = str(interactive.STATIC / "calc.js")
    sent = dict(args)
    sent["__floats__"] = []
    proc = subprocess.run(
        [NODE, driver, calc, consts, "edms.work.plates:layout", json.dumps(sent)],
        capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, proc.stderr
    assert json.loads(proc.stdout)["ok"] is False


# --------------------------------------------------------------------------- #
# clone:dilute
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("args", [
    {"samples": 352, "sample_uL": 10, "fold": 10, "mm_x": 1.05},
    {"samples": 352, "sample_uL": 20, "fold": 10, "mm_x": 1.05},
    {"samples": 4, "sample_uL": 2.0, "fold": 10, "in_place": False},
    {"samples": 1, "sample_uL": 7.5, "fold": 2.5, "mm_x": 1.0},
    {"samples": 0, "sample_uL": 10, "fold": 10},
    {"samples": 96, "sample_uL": 10, "fold": 1},
])
def test_dilute_matches_python(driver, consts, args):
    from edms.bio import clone

    js = run_js(driver, consts, "edms.bio.clone:dilute", args)
    same(js, clone.dilute(**args), where=str(args))


def test_dilute_fold_label_matches_python_g_format(driver, consts):
    """Python writes the heading with format(fold, 'g'), so 2.5 keeps its
    decimal and 10.0 loses it."""
    from edms.bio import clone

    for fold in (10, 10.0, 2.5, 1.25):
        args = {"samples": 2, "sample_uL": 10, "fold": fold}
        js = run_js(driver, consts, "edms.bio.clone:dilute", args)
        same(js, clone.dilute(**args), where=f"fold={fold!r}")


# --------------------------------------------------------------------------- #
# plates:layout, label_col
# --------------------------------------------------------------------------- #
def test_layout_label_col_moves_the_column_and_is_reported(driver, consts):
    """Whatever draws the plate labels wells from this, so the page and the
    package have to agree on which column it is."""
    items = pd.DataFrame({
        "Insert": ["GCaa", "GCtt", "GCgg"],
        "Oligo_ID": ["n1", "n2", "n3"],
    })
    args = {"plate_format": "6-well", "label_col": "Oligo_ID"}
    js = run_js(driver, consts, "edms.work.plates:layout",
                dict(args, items=interactive.table_data(items)))
    want = plates.layout(items, **args)
    same_tables(js, want, where="label_col")
    assert list(want["Plate map"].columns)[:3] == ["Plate", "Well", "Oligo_ID"]
    assert want["Plates"].iloc[0]["Labelled by"] == "Oligo_ID"


def test_layout_without_label_col_falls_back_the_same_way(driver, consts):
    items = pd.DataFrame({"Insert": ["a", "b"], "Oligo_ID": ["n1", "n2"]})
    js = run_js(driver, consts, "edms.work.plates:layout",
                {"items": interactive.table_data(items), "plate_format": "6-well"})
    want = plates.layout(items, "6-well")
    same_tables(js, want, where="no label_col")
    assert want["Plates"].iloc[0]["Labelled by"] == "Insert"


def test_layout_rejects_an_unknown_label_col_on_both_sides(driver, consts):
    items = pd.DataFrame({"Insert": ["a"]})
    with pytest.raises(KeyError, match="to label the wells with"):
        plates.layout(items, "6-well", label_col="nope")
    calc = str(interactive.STATIC / "calc.js")
    sent = {"items": interactive.table_data(items), "plate_format": "6-well",
            "label_col": "nope", "__floats__": []}
    proc = subprocess.run(
        [NODE, driver, calc, consts, "edms.work.plates:layout", json.dumps(sent)],
        capture_output=True, text=True, timeout=60)
    assert json.loads(proc.stdout)["ok"] is False


# --------------------------------------------------------------------------- #
# plate.js, drawing the wells a layout actually uses
# --------------------------------------------------------------------------- #
PLATE_DRIVER = """
const vm = require('vm');
const fs = require('fs');
const ctx = {module: {exports: {}}, console};
ctx.global = ctx;
vm.createContext(ctx);
vm.runInContext(fs.readFileSync(process.argv[2], 'utf8'), ctx);
vm.runInContext(fs.readFileSync(process.argv[3], 'utf8'), ctx);
ctx.EDMS.consts = JSON.parse(process.argv[4]);
const opts = JSON.parse(process.argv[5]);
const svg = ctx.EDMS.plate.svg(opts);
const on = (svg.match(/pv-well pv-on/g) || []).length;
const titles = (svg.match(/<title>([^<]*)<\\/title>/g) || [])
  .map(t => t.replace(/<\\/?title>/g, ''));
process.stdout.write(JSON.stringify({on, titles, caption: /pv-cap-left">([^<]*)</.exec(svg)[1]}));
"""


@pytest.fixture(scope="module")
def plate_driver(tmp_path_factory):
    p = tmp_path_factory.mktemp("js") / "plate_driver.js"
    p.write_text(PLATE_DRIVER, encoding="utf-8")
    return str(p)


def run_plate(plate_driver, consts, opts):
    proc = subprocess.run(
        [NODE, plate_driver,
         str(interactive.STATIC / "calc.js"), str(interactive.STATIC / "plate.js"),
         consts, json.dumps(opts)],
        capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, proc.stderr
    return json.loads(proc.stdout)


def test_plate_fills_from_a1_when_given_only_a_count(plate_driver, consts):
    out = run_plate(plate_driver, consts, {"format": "96-well", "used": 3})
    assert out["on"] == 3
    assert out["titles"][:3] == ["A1", "A2", "A3"]


def test_plate_fills_the_named_wells_when_given_them(plate_driver, consts):
    """A plan that skips the edge wells does not start at A1, and a drawing
    that counted from A1 would put every sample in the wrong well."""
    wells = plates.wells("96-well", skip_outer=True)
    out = run_plate(plate_driver, consts,
                    {"format": "96-well", "used": len(wells), "wells": wells})
    assert out["on"] == 60
    assert "A1" in out["titles"]          # drawn, but empty
    assert "60 of 96 wells" in out["caption"]
    filled = [t for t in out["titles"] if " - " in t]
    assert not filled  # no labels given


def test_plate_labels_follow_the_named_wells(plate_driver, consts):
    out = run_plate(plate_driver, consts, {
        "format": "96-well", "used": 2,
        "wells": ["B2", "C3"], "labels": ["first", "second"]})
    assert out["on"] == 2
    assert "B2 - first" in out["titles"]
    assert "C3 - second" in out["titles"]


def test_plate_ignores_a_well_that_is_not_on_the_plate(plate_driver, consts):
    out = run_plate(plate_driver, consts,
                    {"format": "6-well", "used": 3, "wells": ["A1", "Z9", "B3"]})
    assert out["on"] == 2


def test_every_template_compute_target_is_live():
    """The point of the exercise: no bundled protocol has a table that falls
    back to static. A new template with a new calculator fails here until it
    is ported."""
    from edms.work import template

    targets = set()
    for name in template.list_templates():
        tpl = template.load(name)
        for step in tpl.steps:
            if step.compute:
                targets.add(step.compute)
    missing = sorted(targets - interactive.live_targets())
    assert not missing, f"not ported to calc.js: {', '.join(missing)}"
