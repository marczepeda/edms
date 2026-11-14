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
    """Call one calc.js target under node and return its table."""
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
