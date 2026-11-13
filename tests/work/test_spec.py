"""
Tests for edms.work.spec

Covers reading a spec file, building one from CLI arguments, layering the two,
and resolving input paths relative to the spec file.
"""
import pytest

from edms.work import spec as spec_mod


def _write(tmp_path, body, name="run.toml"):
    p = tmp_path / name
    p.write_text(body, encoding="utf-8")
    return p


def test_load_reads_template_inputs_and_params(tmp_path):
    p = _write(tmp_path, """
template = "minimal"
run_id = "abc"

[inputs]
samples = "samples.csv"

[params]
scale = 7
""")
    s = spec_mod.load(p)
    assert s.template == "minimal"
    assert s.run_id == "abc"
    assert s.inputs == {"samples": "samples.csv"}
    assert s.params == {"scale": 7}


def test_load_without_template_is_rejected(tmp_path):
    p = _write(tmp_path, '[inputs]\nsamples = "s.csv"\n')
    with pytest.raises(ValueError, match='needs a top-level "template"'):
        spec_mod.load(p)


def test_load_missing_file_is_rejected(tmp_path):
    with pytest.raises(FileNotFoundError):
        spec_mod.load(tmp_path / "absent.toml")


def test_relative_input_paths_resolve_against_the_spec_file(tmp_path):
    sub = tmp_path / "proj"
    sub.mkdir()
    p = _write(sub, 'template = "m"\n[inputs]\nsamples = "data/s.csv"\n')
    s = spec_mod.load(p)
    assert s.input_path("samples") == (sub / "data" / "s.csv").resolve()


def test_absolute_input_paths_are_left_alone(tmp_path):
    target = tmp_path / "elsewhere" / "s.csv"
    p = _write(tmp_path, f'template = "m"\n[inputs]\nsamples = "{target}"\n')
    assert spec_mod.load(p).input_path("samples") == target


def test_from_args_parses_pairs_and_types():
    s = spec_mod.from_args(
        template="minimal",
        input=["samples=s.csv"],
        set=["scale=5", "label=widget", "ratio=1.5"],
    )
    assert s.inputs == {"samples": "s.csv"}
    # params go through try_parse(), so numbers arrive as numbers
    assert s.params == {"scale": 5, "label": "widget", "ratio": 1.5}


def test_from_args_rejects_a_pair_without_equals():
    with pytest.raises(ValueError, match="Expected 'name=value'"):
        spec_mod.from_args(template="m", set=["scale"])


def test_from_args_rejects_an_empty_key():
    with pytest.raises(ValueError, match="Expected 'name=value'"):
        spec_mod.from_args(template="m", set=["=5"])


def test_from_args_keeps_paths_with_equals_signs_intact():
    s = spec_mod.from_args(template="m", input=["samples=/tmp/a=b.csv"])
    assert s.inputs == {"samples": "/tmp/a=b.csv"}


def test_merge_layers_cli_over_file(tmp_path):
    base = spec_mod.load(_write(tmp_path, """
template = "minimal"
[inputs]
samples = "s.csv"
extra = "e.csv"
[params]
scale = 2
label = "from_file"
"""))
    override = spec_mod.from_args(template="", input=["extra=other.csv"], set=["scale=9"])
    merged = spec_mod.merge(base, override)

    assert merged.template == "minimal"          # empty override does not win
    assert merged.inputs["samples"] == "s.csv"   # untouched key survives
    assert merged.inputs["extra"] == "other.csv" # overridden
    assert merged.params == {"scale": 9, "label": "from_file"}
    assert merged.source == base.source          # paths still resolve from the file


def test_example_writes_a_loadable_spec(tmp_path):
    out = spec_mod.example(dir=tmp_path)
    assert out.exists()
    s = spec_mod.load(out)
    assert s.template == "pe3_transfection"
    assert "plasmids" in s.inputs
