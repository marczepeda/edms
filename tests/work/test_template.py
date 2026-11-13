"""
Tests for edms.work.template

Covers loading, the bundled template, and every validation rule: duplicate
step ids, forward and unknown dependencies, malformed compute targets, and
step args that reference inputs or params the template never declares.
"""
import textwrap

import pytest

from edms.work import template as template_mod


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #
def test_load_reads_params_inputs_and_steps(tpl):
    assert tpl.name == "minimal"
    assert tpl.params == {"scale": 2, "label": "widget"}
    assert [i.name for i in tpl.inputs] == ["samples", "extra"]
    assert [s.id for s in tpl.steps] == ["prep", "build", "optional_step"]


def test_load_marks_optional_inputs(tpl):
    assert tpl.inputs[0].optional is False
    assert tpl.inputs[1].optional is True


def test_load_keeps_compute_args(tpl):
    build = tpl.step("build")
    assert build.compute == "pandas:DataFrame"
    assert build.args == {"data": {"n": [1, 2, 3]}}
    assert build.depends_on == ["prep"]


def test_step_lookup_raises_for_unknown_id(tpl):
    with pytest.raises(KeyError, match="no step 'nope'"):
        tpl.step("nope")


def test_bundled_pe3_template_loads():
    assert "pe3_transfection" in template_mod.list_templates()
    tpl = template_mod.load("pe3_transfection")
    assert tpl.step("transfection_plan").compute == "edms.bio.transfect:PE3"
    assert tpl.step("ngs_plan").compute == "edms.bio.ngs:pcrs"
    # The NGS sheet is not known at the start of a run, so it must be optional
    assert [i for i in tpl.inputs if i.name == "ngs_samples"][0].optional is True


def test_resolve_path_rejects_unknown_name():
    with pytest.raises(FileNotFoundError, match="Unknown template"):
        template_mod.resolve_path("not_a_real_template")


def test_resolve_path_rejects_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        template_mod.resolve_path(str(tmp_path / "nope.toml"))


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #
def _write(tmp_path, body):
    p = tmp_path / "t.toml"
    p.write_text(textwrap.dedent(body), encoding="utf-8")
    return str(p)


def test_template_without_name_is_rejected(tmp_path):
    with pytest.raises(ValueError, match='needs a top-level "name"'):
        template_mod.load(_write(tmp_path, """
            version = 1
            [[steps]]
            id = "a"
        """))


def test_template_without_steps_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="has no steps"):
        template_mod.load(_write(tmp_path, 'name = "x"\n'))


def test_duplicate_step_ids_are_rejected(tmp_path):
    with pytest.raises(ValueError, match="duplicate step id 'a'"):
        template_mod.load(_write(tmp_path, """
            name = "x"
            [[steps]]
            id = "a"
            [[steps]]
            id = "a"
        """))


def test_unknown_dependency_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="depends on unknown step 'ghost'"):
        template_mod.load(_write(tmp_path, """
            name = "x"
            [[steps]]
            id = "a"
            depends_on = ["ghost"]
        """))


def test_forward_dependency_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="which comes later"):
        template_mod.load(_write(tmp_path, """
            name = "x"
            [[steps]]
            id = "a"
            depends_on = ["b"]
            [[steps]]
            id = "b"
        """))


def test_malformed_compute_target_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="module.path:function"):
        template_mod.load(_write(tmp_path, """
            name = "x"
            [[steps]]
            id = "a"
            compute = "pandas.DataFrame"
        """))


def test_arg_referencing_undeclared_input_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="references input 'ghost'"):
        template_mod.load(_write(tmp_path, """
            name = "x"
            [[steps]]
            id = "a"
            compute = "pandas:DataFrame"
            [steps.args]
            data = "@ghost"
        """))


def test_arg_referencing_undeclared_param_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="references param 'ghost'"):
        template_mod.load(_write(tmp_path, """
            name = "x"
            [[steps]]
            id = "a"
            compute = "pandas:DataFrame"
            [steps.args]
            data = "$ghost"
        """))
