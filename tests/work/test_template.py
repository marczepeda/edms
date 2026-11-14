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
    assert tpl.step("transfection_plan").compute == "edms.work.plates:PE3"
    assert tpl.step("ngs_plan").compute == "edms.bio.ngs:pcrs"
    # The NGS sheet is not known at the start of a run, so it must be optional
    assert [i for i in tpl.inputs if i.name == "ngs_samples"][0].optional is True


def test_bundled_virus_template_loads():
    assert "virus_production" in template_mod.list_templates()
    tpl = template_mod.load("virus_production")
    assert tpl.step("transfection_plan").compute == "edms.work.plates:virus"
    assert [s.id for s in tpl.steps][:2] == ["plate_plan", "seed"]


def test_bundled_ngs_prep_template_loads():
    assert "ngs_prep" in template_mod.list_templates()
    tpl = template_mod.load("ngs_prep")
    assert tpl.step("coverage").compute == "edms.bio.ngs:coverage"
    assert tpl.step("pcr_plan").compute == "edms.bio.ngs:pcrs"
    # gDNA concentrations arrive after the samples are quantified
    assert [i for i in tpl.inputs if i.name == "gdna"][0].optional is True


def test_bundled_transduction_template_loads():
    assert "transduction" in template_mod.list_templates()
    tpl = template_mod.load("transduction")
    assert tpl.step("transduce").compute == "edms.bio.culture:transduction"
    # Selection and treatment are the same calculation with different units
    assert tpl.step("selection").compute == "edms.bio.culture:drug_mm"
    assert tpl.step("treatment").compute == "edms.bio.culture:drug_mm"


def test_bundled_competitive_growth_template_loads():
    assert "competitive_growth" in template_mod.list_templates()
    tpl = template_mod.load("competitive_growth")
    assert tpl.step("seed_plan").compute == "edms.bio.culture:competitive_growth"
    # Confluency has to be read before anything is calculated from it
    assert [s.id for s in tpl.steps][:2] == ["confluency", "seed_plan"]


def test_bundled_golden_gate_template_loads():
    assert "golden_gate" in template_mod.list_templates()
    tpl = template_mod.load("golden_gate")
    assert tpl.step("anneal").compute == "edms.bio.clone:anneal_mm"
    assert tpl.step("assembly_plan").compute == "edms.bio.clone:golden_gate"


def test_bundled_gibson_template_loads():
    assert "gibson_assembly" in template_mod.list_templates()
    tpl = template_mod.load("gibson_assembly")
    assert tpl.step("pcr_plan").compute == "edms.bio.clone:rxn_mm"
    assert tpl.step("assembly_plan").compute == "edms.bio.clone:gibson"
    # Fragments cannot be quantified until the gel extraction is done
    assert [i for i in tpl.inputs if i.name == "fragments"][0].optional is True


def test_bundled_clone_isolation_template_loads():
    assert "clone_isolation" in template_mod.list_templates()
    tpl = template_mod.load("clone_isolation")
    assert tpl.step("sort_plan").compute == "edms.bio.culture:sort_plan"
    assert tpl.step("expansion_plan").compute == "edms.work.plates:expansion"
    assert tpl.step("freeze_plan").compute == "edms.bio.culture:freeze_down"


def test_clone_isolation_sort_waits_on_the_conditioned_media():
    """The media takes a week to make and the sort cannot run without it."""
    tpl = template_mod.load("clone_isolation")
    assert "conditioned_media" in tpl.step("sort").depends_on
    assert [s.id for s in tpl.steps][0] == "conditioned_media"


@pytest.mark.parametrize("name", template_mod.list_templates())
def test_every_bundled_template_validates(name):
    """Each shipped template must load and pass every validation rule."""
    tpl = template_mod.load(name)
    assert tpl.steps
    assert tpl.title


@pytest.mark.parametrize("name", template_mod.list_templates())
def test_every_bundled_template_resolves_its_compute_targets(name):
    """A compute target that cannot be imported would only fail at render."""
    import importlib

    for step in template_mod.load(name).steps:
        if not step.compute:
            continue
        module, _, func = step.compute.partition(":")
        assert callable(getattr(importlib.import_module(module), func))


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
