"""
Tests for edms.work.run

Covers run creation (frozen inputs, resolved params), every status transition,
the dependency guard, and deviations - the part of edms.work that the one-shot
planners elsewhere in the package do not have.
"""
import json

import pytest

from edms.work import run as run_mod
from edms.work import spec as spec_mod


# --------------------------------------------------------------------------- #
# init()
# --------------------------------------------------------------------------- #
def test_init_creates_a_ledger_with_every_step_pending(started_run):
    run_dir, state = started_run
    assert (run_dir / "run.json").exists()
    assert [s["id"] for s in state["steps"]] == ["prep", "build", "optional_step"]
    assert {s["status"] for s in state["steps"]} == {"pending"}


def test_init_resolves_params_template_defaults_then_spec(started_run):
    _, state = started_run
    assert state["params"] == {"scale": 5, "label": "widget"}


def test_init_freezes_input_tables_into_the_run(started_run, samples_csv):
    run_dir, state = started_run
    frozen = run_dir / state["inputs"]["samples"]["frozen"]
    assert frozen.exists()
    assert frozen.read_text() == samples_csv.read_text()
    # The original path is recorded, but the copy is what the run uses
    assert state["inputs"]["samples"]["source"] == str(samples_csv)


def test_init_is_unaffected_by_later_edits_to_the_source_table(started_run, samples_csv):
    run_dir, state = started_run
    samples_csv.write_text("id\nchanged\n")
    frozen = run_dir / state["inputs"]["samples"]["frozen"]
    assert "changed" not in frozen.read_text()


def test_init_names_the_run_directory_by_date_and_label(started_run):
    run_dir, state = started_run
    assert run_dir.name == state["run_id"]
    assert state["run_id"].endswith("_unit")


def test_init_rejects_a_missing_required_input(tpl, template_file, tmp_path):
    spec = spec_mod.RunSpec(template=str(template_file), inputs={}, source=tmp_path / "r.toml")
    with pytest.raises(ValueError, match="needs input\\(s\\).*samples"):
        run_mod.init(spec, tpl=tpl, dir=tmp_path / "runs")


def test_init_cleans_up_after_rejecting_a_missing_input(tpl, template_file, tmp_path):
    spec = spec_mod.RunSpec(template=str(template_file), run_id="x",
                            inputs={}, source=tmp_path / "r.toml")
    with pytest.raises(ValueError):
        run_mod.init(spec, tpl=tpl, dir=tmp_path / "runs")
    # The half-built run directory must not be left behind
    assert list((tmp_path / "runs").glob("*")) == []


def test_init_rejects_an_input_file_that_does_not_exist(tpl, template_file, tmp_path):
    spec = spec_mod.RunSpec(template=str(template_file),
                            inputs={"samples": str(tmp_path / "ghost.csv")},
                            source=tmp_path / "r.toml")
    with pytest.raises(FileNotFoundError, match="not found"):
        run_mod.init(spec, tpl=tpl, dir=tmp_path / "runs")


def test_init_rejects_an_undeclared_parameter(tpl, template_file, samples_csv, tmp_path):
    spec = spec_mod.RunSpec(template=str(template_file),
                            inputs={"samples": str(samples_csv)},
                            params={"ghost": 1}, source=tmp_path / "r.toml")
    with pytest.raises(ValueError, match="no parameter\\(s\\): ghost"):
        run_mod.init(spec, tpl=tpl, dir=tmp_path / "runs")


def test_init_rejects_an_undeclared_input(tpl, template_file, samples_csv, tmp_path):
    spec = spec_mod.RunSpec(template=str(template_file),
                            inputs={"samples": str(samples_csv), "ghost": str(samples_csv)},
                            source=tmp_path / "r.toml")
    with pytest.raises(ValueError, match="does not declare input\\(s\\): ghost"):
        run_mod.init(spec, tpl=tpl, dir=tmp_path / "runs")


def test_init_refuses_to_clobber_an_existing_run(run_spec, tpl, tmp_path):
    run_mod.init(run_spec, tpl=tpl, dir=tmp_path / "runs")
    with pytest.raises(FileExistsError, match="already exists"):
        run_mod.init(run_spec, tpl=tpl, dir=tmp_path / "runs")


def test_init_overwrite_replaces_an_existing_run(run_spec, tpl, tmp_path):
    run_dir, _ = run_mod.init(run_spec, tpl=tpl, dir=tmp_path / "runs")
    (run_dir / "stray.txt").write_text("x")
    run_mod.init(run_spec, tpl=tpl, dir=tmp_path / "runs", overwrite=True)
    assert not (run_dir / "stray.txt").exists()


# --------------------------------------------------------------------------- #
# Status transitions
# --------------------------------------------------------------------------- #
def test_done_marks_a_step_finished_and_timestamps_it(started_run):
    run_dir, _ = started_run
    state = run_mod.done(run_dir, "prep")
    entry = run_mod._find(state, "prep")
    assert entry["status"] == "done"
    assert entry["started"] and entry["finished"]


def test_done_on_a_pending_step_starts_it_too(started_run):
    run_dir, _ = started_run
    entry = run_mod._find(run_mod.done(run_dir, "prep"), "prep")
    assert entry["started"] is not None


def test_done_refuses_a_step_with_unfinished_dependencies(started_run):
    run_dir, _ = started_run
    with pytest.raises(ValueError, match="depends on unfinished step\\(s\\): prep"):
        run_mod.done(run_dir, "build")


def test_force_overrides_the_dependency_guard(started_run):
    run_dir, _ = started_run
    state = run_mod.done(run_dir, "build", force=True)
    assert run_mod._find(state, "build")["status"] == "done"


def test_a_skipped_dependency_still_unblocks_the_next_step(started_run):
    run_dir, _ = started_run
    run_mod.skip(run_dir, "prep", why="not needed")
    state = run_mod.done(run_dir, "build")
    assert run_mod._find(state, "build")["status"] == "done"


def test_start_refuses_a_step_that_is_already_done(started_run):
    run_dir, _ = started_run
    run_mod.done(run_dir, "prep")
    with pytest.raises(ValueError, match="already done"):
        run_mod.start(run_dir, "prep")


def test_reopen_returns_a_step_to_pending(started_run):
    run_dir, _ = started_run
    run_mod.done(run_dir, "prep")
    entry = run_mod._find(run_mod.reopen(run_dir, "prep", why="repeat"), "prep")
    assert entry["status"] == "pending"
    assert entry["started"] is None and entry["finished"] is None
    assert "repeat" in entry["notes"][-1]["text"]


def test_fail_keeps_the_step_open(started_run):
    run_dir, _ = started_run
    state = run_mod.fail(run_dir, "prep", why="contamination")
    assert run_mod._find(state, "prep")["status"] == "failed"
    assert "prep" in [e["id"] for e in run_mod.next_steps(state)]


def test_note_appends_without_changing_status(started_run):
    run_dir, _ = started_run
    state = run_mod.note(run_dir, "prep", text="looked fine")
    entry = run_mod._find(state, "prep")
    assert entry["status"] == "pending"
    assert entry["notes"][-1]["text"] == "looked fine"


def test_unknown_step_is_rejected(started_run):
    run_dir, _ = started_run
    with pytest.raises(KeyError, match="no step 'ghost'"):
        run_mod.note(run_dir, "ghost", text="x")


# --------------------------------------------------------------------------- #
# Deviations
# --------------------------------------------------------------------------- #
def test_deviate_updates_the_parameter_and_records_the_old_value(started_run):
    run_dir, _ = started_run
    state = run_mod.deviate(run_dir, "prep", param="scale", to="3", why="short on cells")
    assert state["params"]["scale"] == 3           # parsed to an int, not "3"
    dev = run_mod._find(state, "prep")["deviations"][-1]
    assert (dev["from"], dev["to"], dev["why"]) == (5, 3, "short on cells")


def test_deviate_is_structured_rather_than_prose(started_run):
    run_dir, _ = started_run
    run_mod.deviate(run_dir, "prep", param="scale", to=3)
    dev = run_mod.load(run_dir)["steps"][0]["deviations"][-1]
    assert set(dev) == {"ts", "param", "from", "to", "why"}


def test_deviate_rejects_an_unknown_parameter(started_run):
    run_dir, _ = started_run
    with pytest.raises(KeyError, match="no parameter 'ghost'"):
        run_mod.deviate(run_dir, "prep", param="ghost", to=1)


# --------------------------------------------------------------------------- #
# Queries and persistence
# --------------------------------------------------------------------------- #
def test_next_steps_only_offers_unblocked_work(started_run):
    run_dir, state = started_run
    assert [e["id"] for e in run_mod.next_steps(state)] == ["prep"]
    state = run_mod.done(run_dir, "prep")
    assert [e["id"] for e in run_mod.next_steps(state)] == ["build"]


def test_progress_counts_done_and_skipped_as_closed(started_run):
    run_dir, _ = started_run
    run_mod.done(run_dir, "prep")
    run_mod.skip(run_dir, "build", why="n/a")
    assert run_mod.progress(run_mod.load(run_dir)) == (2, 3)


def test_events_record_an_audit_trail(started_run):
    run_dir, _ = started_run
    run_mod.done(run_dir, "prep", note="ok")
    run_mod.deviate(run_dir, "prep", param="scale", to=1)
    actions = [e["action"] for e in run_mod.load(run_dir)["events"]]
    assert actions == ["init", "done", "deviate"]


def test_summary_shows_status_marks_and_what_is_next(started_run):
    run_dir, _ = started_run
    run_mod.done(run_dir, "prep")
    text = run_mod.summary(run_mod.load(run_dir))
    assert "[x] prep" in text
    assert "[ ] build" in text
    assert "next: build" in text


def test_load_rejects_a_ledger_from_a_future_edms(started_run):
    run_dir, state = started_run
    state["state_version"] = 99
    (run_dir / "run.json").write_text(json.dumps(state))
    with pytest.raises(ValueError, match="state_version"):
        run_mod.load(run_dir)


def test_load_without_a_run_file_is_rejected(tmp_path):
    with pytest.raises(FileNotFoundError, match="No run.json"):
        run_mod.load(tmp_path)
