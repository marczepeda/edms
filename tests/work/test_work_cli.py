"""
Tests for edms.work.cli

Argument-parsing wiring, plus an end-to-end pass through the commands a run
actually goes through: init -> done -> deviate -> plan. The CLI functions are
called the way edms/main.py calls them, as func(**vars(args)).
"""
import argparse

import pandas as pd
import pytest

from edms.work import cli
from edms.work import run as run_mod


def _build_parser():
    parser = argparse.ArgumentParser(prog="edms")
    subparsers = parser.add_subparsers(dest="command")
    # As in edms/main.py: a concrete formatter, since the None default is
    # unusable for --help.
    cli.add_subparser(subparsers, argparse.HelpFormatter)
    return parser


def _dispatch(parser, argv):
    """Parse and call, exactly as edms/main.py does."""
    args = vars(parser.parse_args(argv))
    func = args.pop("func")
    args.pop("command", None)
    return func(**args)


# --------------------------------------------------------------------------- #
# Wiring
# --------------------------------------------------------------------------- #
def test_add_subparser_builds_without_error():
    assert _build_parser() is not None


def test_work_help_does_not_crash(capsys):
    parser = _build_parser()
    with pytest.raises(SystemExit) as e:
        parser.parse_args(["work", "--help"])
    assert e.value.code == 0
    assert "Plan experiments" in capsys.readouterr().out


@pytest.mark.parametrize("name", [
    "templates", "show", "where", "example", "init", "plan", "status",
    "start", "done", "skip", "fail", "reopen", "note", "deviate",
])
def test_every_subcommand_is_attached(name, capsys):
    parser = _build_parser()
    with pytest.raises(SystemExit) as e:
        parser.parse_args(["work", name, "--help"])
    assert e.value.code == 0


def test_init_parses_repeatable_inputs_and_overrides():
    parser = _build_parser()
    args = parser.parse_args([
        "work", "init", "-t", "minimal",
        "-I", "samples=s.csv", "-I", "extra=e.csv",
        "-S", "scale=5", "-d", "runs/",
    ])
    assert args.input == ["samples=s.csv", "extra=e.csv"]
    assert args.set == ["scale=5"]
    assert args.dir == "runs/"
    assert args.plan is True


def test_no_plan_flag_turns_rendering_off():
    parser = _build_parser()
    assert parser.parse_args(["work", "init", "-t", "m", "-n"]).plan is False


def test_plan_no_compute_flag_turns_calculation_off():
    parser = _build_parser()
    assert parser.parse_args(["work", "plan", "r", "-n"]).compute is False


# --------------------------------------------------------------------------- #
# Errors arrive as one line, not a traceback
# --------------------------------------------------------------------------- #
def test_unknown_template_exits_cleanly():
    parser = _build_parser()
    with pytest.raises(SystemExit, match="Unknown template"):
        _dispatch(parser, ["work", "show", "no_such_template"])


def test_init_without_template_or_spec_exits_cleanly():
    parser = _build_parser()
    with pytest.raises(SystemExit, match="Pass a run spec"):
        _dispatch(parser, ["work", "init"])


def test_missing_run_directory_exits_cleanly(tmp_path):
    parser = _build_parser()
    with pytest.raises(SystemExit, match="No run.json"):
        _dispatch(parser, ["work", "status", str(tmp_path)])


# --------------------------------------------------------------------------- #
# End to end
# --------------------------------------------------------------------------- #
def test_full_cli_round_trip(template_file, samples_csv, tmp_path, capsys):
    parser = _build_parser()
    runs = tmp_path / "runs"

    _dispatch(parser, [
        "work", "init", "-t", str(template_file),
        "-I", f"samples={samples_csv}", "-S", "scale=4",
        "-r", "cli", "-d", str(runs),
    ])
    run_dir = next(runs.glob("*_cli"))
    assert (run_dir / "run.json").exists()
    assert (run_dir / "protocol.html").exists()
    assert (run_dir / "outputs" / "build.xlsx").exists()

    _dispatch(parser, ["work", "done", str(run_dir), "-s", "prep", "-n", "ok"])
    assert run_mod._find(run_mod.load(run_dir), "prep")["status"] == "done"

    _dispatch(parser, ["work", "deviate", str(run_dir),
                       "-s", "prep", "-p", "scale", "-v", "8", "-w", "ran short"])
    assert run_mod.load(run_dir)["params"]["scale"] == 8

    _dispatch(parser, ["work", "plan", str(run_dir), "-f", "md"])
    md = (run_dir / "protocol.md").read_text()
    assert "Prepare 8 x widget." in md          # deviation reached the prose
    assert "Deviation: scale 4 -> 8" in md

    _dispatch(parser, ["work", "status", str(run_dir)])
    assert "[x] prep" in capsys.readouterr().out


def test_init_from_a_spec_file(template_file, samples_csv, tmp_path):
    spec = tmp_path / "run.toml"
    spec.write_text(f"""
template = "{template_file}"
run_id = "fromspec"

[inputs]
samples = "{samples_csv}"

[params]
scale = 3
""", encoding="utf-8")

    parser = _build_parser()
    _dispatch(parser, ["work", "init", "-i", str(spec), "-d", str(tmp_path / "runs")])
    run_dir = next((tmp_path / "runs").glob("*_fromspec"))
    assert run_mod.load(run_dir)["params"]["scale"] == 3


def test_cli_overrides_beat_the_spec_file(template_file, samples_csv, tmp_path):
    spec = tmp_path / "run.toml"
    spec.write_text(f"""
template = "{template_file}"
run_id = "override"
[inputs]
samples = "{samples_csv}"
[params]
scale = 3
""", encoding="utf-8")

    parser = _build_parser()
    _dispatch(parser, ["work", "init", "-i", str(spec), "-S", "scale=11",
                       "-d", str(tmp_path / "runs")])
    run_dir = next((tmp_path / "runs").glob("*_override"))
    assert run_mod.load(run_dir)["params"]["scale"] == 11


def test_example_writes_a_spec(tmp_path, capsys):
    parser = _build_parser()
    _dispatch(parser, ["work", "example", "-d", str(tmp_path)])
    assert (tmp_path / "run.toml").exists()


def test_templates_lists_the_bundled_template(capsys):
    parser = _build_parser()
    _dispatch(parser, ["work", "templates"])
    assert "pe3_transfection" in capsys.readouterr().out


def test_where_prints_a_real_path(capsys):
    parser = _build_parser()
    _dispatch(parser, ["work", "where", "pe3_transfection"])
    out = capsys.readouterr().out.strip()
    assert out.endswith("pe3_transfection.toml")


# --------------------------------------------------------------------------- #
# Adding an input part-way through a run
# --------------------------------------------------------------------------- #
def test_input_command_adds_a_table_and_re_renders(template_file, samples_csv, tmp_path):
    parser = _build_parser()
    runs = tmp_path / "runs"
    _dispatch(parser, ["work", "init", "-t", str(template_file),
                       "-I", f"samples={samples_csv}", "-r", "late", "-d", str(runs)])
    run_dir = next(runs.glob("*_late"))
    assert "Not calculated" in (run_dir / "protocol.html").read_text()

    extra = tmp_path / "extra.csv"
    pd.DataFrame({"id": [1, 2]}).to_csv(extra, index=False)
    _dispatch(parser, ["work", "input", str(run_dir), "-I", f"extra={extra}"])

    assert "extra" in run_mod.load(run_dir)["inputs"]
    assert "Not calculated" not in (run_dir / "protocol.html").read_text()


def test_input_command_requires_a_pair():
    parser = _build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["work", "input", "somewhere"])


def test_input_command_refuses_to_swap_without_replace(template_file, samples_csv, tmp_path):
    parser = _build_parser()
    runs = tmp_path / "runs"
    _dispatch(parser, ["work", "init", "-t", str(template_file),
                       "-I", f"samples={samples_csv}", "-r", "swap", "-d", str(runs)])
    run_dir = next(runs.glob("*_swap"))
    with pytest.raises(SystemExit, match="already froze"):
        _dispatch(parser, ["work", "input", str(run_dir),
                           "-I", f"samples={samples_csv}", "-n"])
