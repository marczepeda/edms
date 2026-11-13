"""
Tests for edms.work.plan

Covers parameter substitution, the @input/$param arg sigils, normalising
whatever a compute function returns, skipping a step that needs an optional
input, and the strict/forgiving behaviour when a calculation raises.

The PE3 integration test at the bottom checks that a real edms.bio function
can be driven from a template.
"""
import pandas as pd
import pytest

from edms.work import plan as plan_mod
from edms.work import run as run_mod
from edms.work import spec as spec_mod
from edms.work import template as template_mod


# --------------------------------------------------------------------------- #
# fmt()
# --------------------------------------------------------------------------- #
def test_fmt_substitutes_params():
    assert plan_mod.fmt("Use {scale} mL", {"scale": 5}) == "Use 5 mL"


def test_fmt_leaves_unknown_keys_alone():
    assert plan_mod.fmt("Use {ghost} mL", {"scale": 5}) == "Use {ghost} mL"


def test_fmt_survives_a_stray_brace():
    text = "Spin at 300 g { then resuspend"
    assert plan_mod.fmt(text, {}) == text


def test_fmt_handles_empty_text():
    assert plan_mod.fmt("", {"a": 1}) == ""


# --------------------------------------------------------------------------- #
# Argument resolution
# --------------------------------------------------------------------------- #
def test_input_refs_finds_nested_references():
    args = {"a": "@one", "b": ["@two", 3], "c": {"d": "@three"}, "e": "$param"}
    assert plan_mod.input_refs(args) == {"one", "two", "three"}


def test_input_refs_ignores_escaped_sigils():
    assert plan_mod.input_refs({"a": "@@literal"}) == set()


def test_resolve_arg_maps_sigils_to_inputs_and_params():
    inputs = {"samples": pd.DataFrame({"id": [1]})}
    params = {"scale": 5}
    assert plan_mod._resolve_arg("$scale", inputs, params, "s", "k") == 5
    assert plan_mod._resolve_arg("@samples", inputs, params, "s", "k") is inputs["samples"]
    assert plan_mod._resolve_arg("plain", inputs, params, "s", "k") == "plain"


def test_resolve_arg_unescapes_doubled_sigils():
    assert plan_mod._resolve_arg("@@name", {}, {}, "s", "k") == "@name"
    assert plan_mod._resolve_arg("$$name", {}, {}, "s", "k") == "$name"


def test_resolve_arg_reports_a_missing_reference_usefully():
    with pytest.raises(KeyError, match="wants parameter 'ghost'"):
        plan_mod._resolve_arg("$ghost", {}, {"scale": 1}, "mystep", "myarg")


# --------------------------------------------------------------------------- #
# Output normalisation
# --------------------------------------------------------------------------- #
def test_tables_names_a_bare_dataframe_from_the_template():
    df = pd.DataFrame({"a": [1]})
    assert plan_mod._tables(df, ["built"]) == {"built": df}


def test_tables_keeps_dict_keys():
    out = plan_mod._tables({"first": pd.DataFrame({"a": [1]})}, [])
    assert list(out) == ["first"]


def test_tables_flattens_a_tuple_of_dicts():
    # edms.bio.ngs.pcrs() returns a tuple of dicts of DataFrames
    result = ({"pcr1": pd.DataFrame({"a": [1]})}, {"pcr2": pd.DataFrame({"b": [2]})})
    assert list(plan_mod._tables(result, [])) == ["pcr1", "pcr2"]


def test_tables_disambiguates_colliding_names():
    result = ({"same": pd.DataFrame({"a": [1]})}, {"same": pd.DataFrame({"b": [2]})})
    assert list(plan_mod._tables(result, [])) == ["same", "same_2"]


def test_tables_ignores_non_dataframes():
    out = plan_mod._tables({"a": "not a table", "b": pd.DataFrame({"x": [1]})}, [])
    assert list(out) == ["b"]


def test_sheet_name_truncates_and_strips_illegal_characters():
    used = set()
    assert plan_mod._sheet_name("a/b:c" + "x" * 40, used) == ("a_b_c" + "x" * 40)[:31]


def test_sheet_name_disambiguates_within_a_workbook():
    used = set()
    first = plan_mod._sheet_name("sheet", used)
    second = plan_mod._sheet_name("sheet", used)
    assert (first, second) == ("sheet", "sheet_2")


# --------------------------------------------------------------------------- #
# resolve()
# --------------------------------------------------------------------------- #
def test_resolve_fills_prose_from_the_runs_params(started_run):
    run_dir, _ = started_run
    protocol = plan_mod.resolve(run_dir)
    assert protocol.steps[0].text == "Prepare 5 x widget."


def test_resolve_picks_up_a_deviation_without_re_initialising(started_run):
    run_dir, _ = started_run
    run_mod.deviate(run_dir, "prep", param="scale", to=9)
    assert plan_mod.resolve(run_dir).steps[0].text == "Prepare 9 x widget."


def test_resolve_runs_compute_steps_and_writes_them(started_run):
    run_dir, _ = started_run
    protocol = plan_mod.resolve(run_dir)
    build = [s for s in protocol.steps if s.id == "build"][0]
    assert list(build.tables) == ["built"]
    assert (run_dir / "outputs" / "build.xlsx").exists()


def test_resolve_records_outputs_on_the_ledger(started_run):
    run_dir, _ = started_run
    plan_mod.resolve(run_dir)
    entry = run_mod._find(run_mod.load(run_dir), "build")
    assert entry["outputs"] == ["outputs/build.xlsx"]


def test_resolve_carries_the_status_through_from_the_ledger(started_run):
    run_dir, _ = started_run
    run_mod.done(run_dir, "prep", note="went fine")
    prep = plan_mod.resolve(run_dir).steps[0]
    assert prep.status == "done"
    assert prep.notes[-1]["text"] == "went fine"


def test_resolve_skips_a_step_whose_optional_input_is_absent(started_run):
    run_dir, _ = started_run
    step = [s for s in plan_mod.resolve(run_dir).steps if s.id == "optional_step"][0]
    assert step.compute_skipped and "extra" in step.compute_skipped
    assert step.error is None       # absent-by-design is not a failure
    assert step.tables == {}


def test_resolve_with_compute_off_runs_nothing(started_run):
    run_dir, _ = started_run
    protocol = plan_mod.resolve(run_dir, compute=False)
    assert all(s.tables == {} for s in protocol.steps)
    assert not (run_dir / "outputs").exists()


def test_resolve_records_a_failing_calculation_without_aborting(started_run, tmp_path):
    run_dir, state = started_run
    # Point the build step at a function that will raise on these arguments
    state["template"]["source"] = str(_broken_template(tmp_path))
    run_mod.save(state, run_dir)
    protocol = plan_mod.resolve(run_dir)
    assert protocol.errors and protocol.errors[0][0] == "build"
    assert protocol.steps[0].text == "Prepare 5 x widget."   # the rest still renders


def test_resolve_strict_raises_on_a_failing_calculation(started_run, tmp_path):
    run_dir, state = started_run
    state["template"]["source"] = str(_broken_template(tmp_path))
    run_mod.save(state, run_dir)
    with pytest.raises(Exception):
        plan_mod.resolve(run_dir, strict=True)


def _broken_template(tmp_path):
    p = tmp_path / "broken.toml"
    p.write_text("""
name = "minimal"
version = 1
[params]
scale = 2
label = "widget"
[[inputs]]
name = "samples"
[[steps]]
id = "prep"
title = "Prepare"
text = "Prepare {scale} x {label}."
[[steps]]
id = "build"
title = "Build"
compute = "edms.work.plan:call"
[steps.args]
compute = "not_a_module:nope"
args = {}
""", encoding="utf-8")
    return p


def test_call_reports_an_unimportable_module():
    with pytest.raises(ImportError, match="Cannot import"):
        plan_mod.call("no_such_module_xyz:f", {})


def test_call_reports_a_missing_function():
    with pytest.raises(AttributeError, match="has no function 'nope'"):
        plan_mod.call("pandas:nope", {})


# --------------------------------------------------------------------------- #
# load_inputs()
# --------------------------------------------------------------------------- #
def test_load_inputs_reads_the_frozen_copies(started_run):
    run_dir, state = started_run
    inputs = plan_mod.load_inputs(state, run_dir)
    assert list(inputs["samples"]["id"]) == ["s1", "s2"]


def test_load_inputs_reports_a_missing_frozen_file(started_run):
    run_dir, state = started_run
    (run_dir / state["inputs"]["samples"]["frozen"]).unlink()
    with pytest.raises(FileNotFoundError, match="is missing"):
        plan_mod.load_inputs(state, run_dir)


def test_load_inputs_requires_a_sheet_for_a_multi_tab_workbook(started_run, tmp_path):
    run_dir, state = started_run
    book = run_dir / "inputs" / "samples.xlsx"
    with pd.ExcelWriter(book) as writer:
        pd.DataFrame({"id": [1]}).to_excel(writer, sheet_name="one", index=False)
        pd.DataFrame({"id": [2]}).to_excel(writer, sheet_name="two", index=False)
    state["inputs"]["samples"]["frozen"] = "inputs/samples.xlsx"
    with pytest.raises(ValueError, match="several tabs"):
        plan_mod.load_inputs(state, run_dir)


# --------------------------------------------------------------------------- #
# Integration: driving a real edms.bio function from a template
# --------------------------------------------------------------------------- #
def test_pe3_template_computes_the_transfection_plan(tmp_path):
    plasmids = tmp_path / "plasmids.csv"
    pd.DataFrame({
        "Plasmid": ["epeg1", "ng1", "pMUZ86.7"],
        "Description": ["epeg", "ng", "PE"],
        "Colony": ["c1", "c2", "c3"],
        "ng/uL": [100.0, 200.0, 300.0],
    }).to_csv(plasmids, index=False)

    epegs = tmp_path / "epegRNAs.csv"
    pd.DataFrame({"pegRNA_number": [1], "Name": ["epeg1"]}).to_csv(epegs, index=False)
    ngs = tmp_path / "ngRNAs.csv"
    pd.DataFrame({"pegRNA_number": [1], "Name": ["ng1"]}).to_csv(ngs, index=False)

    tpl = template_mod.load("pe3_transfection")
    run_spec = spec_mod.RunSpec(
        template="pe3_transfection",
        run_id="pe3_it",
        inputs={"plasmids": str(plasmids), "epegRNAs": str(epegs), "ngRNAs": str(ngs)},
        source=tmp_path / "run.toml",
    )
    run_dir, _ = run_mod.init(run_spec, tpl=tpl, dir=tmp_path / "runs")
    protocol = plan_mod.resolve(run_dir)

    assert protocol.errors == []
    step = [s for s in protocol.steps if s.id == "transfection_plan"][0]
    tube_A = step.tables["Transfection"].iloc[0]
    # Defaults: mm_x=1.1, reps=3, epegRNA_ng=66 at 100 ng/uL -> 2.18
    assert tube_A["epegRNA uL"] == pytest.approx(2.18)

    # The NGS sheet was not supplied, so that step waits rather than failing
    ngs_step = [s for s in protocol.steps if s.id == "ngs_plan"][0]
    assert ngs_step.compute_skipped and ngs_step.error is None


def test_pe3_template_params_reach_the_bio_function(tmp_path):
    plasmids = tmp_path / "plasmids.csv"
    pd.DataFrame({
        "Plasmid": ["epeg1", "ng1", "pMUZ86.7"],
        "Description": ["epeg", "ng", "PE"],
        "Colony": ["c1", "c2", "c3"],
        "ng/uL": [100.0, 200.0, 300.0],
    }).to_csv(plasmids, index=False)
    epegs = tmp_path / "epegRNAs.csv"
    pd.DataFrame({"pegRNA_number": [1], "Name": ["epeg1"]}).to_csv(epegs, index=False)
    ngs = tmp_path / "ngRNAs.csv"
    pd.DataFrame({"pegRNA_number": [1], "Name": ["ng1"]}).to_csv(ngs, index=False)

    run_spec = spec_mod.RunSpec(
        template="pe3_transfection", run_id="pe3_reps",
        inputs={"plasmids": str(plasmids), "epegRNAs": str(epegs), "ngRNAs": str(ngs)},
        params={"reps": 6}, source=tmp_path / "run.toml")
    run_dir, _ = run_mod.init(run_spec, dir=tmp_path / "runs")
    protocol = plan_mod.resolve(run_dir)

    step = [s for s in protocol.steps if s.id == "transfection_plan"][0]
    # reps=6: 1.1 * 6 * 66 / 100 = 4.356 -> 4.36
    assert step.tables["Transfection"].iloc[0]["epegRNA uL"] == pytest.approx(4.36)
