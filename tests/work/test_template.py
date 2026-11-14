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


def test_virus_template_plans_the_small_format_dilutions():
    """Scaling below a 12-well needs working stocks, so the template has a
    step for them rather than leaving sub-microlitre volumes in the plan."""
    tpl = template_mod.load("virus_production")
    step = tpl.step("dilutions")
    assert step.compute == "edms.work.plates:virus_dilutions"
    assert "transfection_plan" in step.depends_on
    assert tpl.params["min_uL"] == 1.0


def test_virus_template_tells_you_the_per_well_complex_volume():
    """The master mix table never says what one well takes, and at 96-well
    the 6-well number would overflow it."""
    tpl = template_mod.load("virus_production")
    step = tpl.step("transfect")
    assert step.compute == "edms.work.plates:virus_per_well"
    assert any("95-99%" in c for c in step.checks)


# --------------------------------------------------------------------------- #
# Components
# --------------------------------------------------------------------------- #
def _component(tmp_path, body, name="c.toml"):
    p = tmp_path / name
    p.write_text(textwrap.dedent(body), encoding="utf-8")
    return p


def test_bundled_components_are_listed():
    assert "cloning_tail" in template_mod.list_components()


def test_resolve_component_rejects_an_unknown_name():
    with pytest.raises(FileNotFoundError, match="Unknown component"):
        template_mod.resolve_component("not_a_component")


def test_include_expands_in_the_order_it_is_written(tmp_path):
    comp = _component(tmp_path, """
        name = "mid"
        [[steps]]
        id = "b"
        [[steps]]
        id = "c"
    """)
    tpl = template_mod.load(_write(tmp_path, f"""
        name = "x"
        [[steps]]
        id = "a"
        [[steps]]
        include = "{comp}"
        [[steps]]
        id = "d"
        depends_on = ["c"]
    """))
    assert [s.id for s in tpl.steps] == ["a", "b", "c", "d"]


def test_included_steps_remember_where_they_came_from(tmp_path):
    comp = _component(tmp_path, """
        name = "mid"
        [[steps]]
        id = "b"
    """)
    tpl = template_mod.load(_write(tmp_path, f"""
        name = "x"
        [[steps]]
        id = "a"
        [[steps]]
        include = "{comp}"
    """))
    assert tpl.step("a").component is None
    assert tpl.step("b").component == "mid"
    assert tpl.components == ["mid"]


def test_include_can_say_what_it_follows(tmp_path):
    """A component cannot name the step it comes after, because it does not
    know which protocol it is in."""
    comp = _component(tmp_path, """
        name = "tail"
        [[steps]]
        id = "b"
        [[steps]]
        id = "c"
        depends_on = ["b"]
    """)
    tpl = template_mod.load(_write(tmp_path, f"""
        name = "x"
        [[steps]]
        id = "a"
        [[steps]]
        include = "{comp}"
        depends_on = ["a"]
    """))
    assert tpl.step("b").depends_on == ["a"]      # entry step attaches
    assert tpl.step("c").depends_on == ["b"]      # internal one is untouched


def test_template_params_win_over_a_components_defaults(tmp_path):
    """Including a component must never change a number the protocol set."""
    comp = _component(tmp_path, """
        name = "mid"
        [params]
        temp = 30
        other = 1
        [[steps]]
        id = "b"
    """)
    tpl = template_mod.load(_write(tmp_path, f"""
        name = "x"
        [params]
        temp = 37
        [[steps]]
        include = "{comp}"
    """))
    assert tpl.params["temp"] == 37
    assert tpl.params["other"] == 1


def test_an_include_can_override_a_component_param(tmp_path):
    comp = _component(tmp_path, """
        name = "mid"
        [params]
        temp = 30
        [[steps]]
        id = "b"
    """)
    tpl = template_mod.load(_write(tmp_path, f"""
        name = "x"
        [[steps]]
        include = "{comp}"
        [steps.params]
        temp = 25
    """))
    assert tpl.params["temp"] == 25


def test_an_include_cannot_set_a_param_the_component_does_not_have(tmp_path):
    """Otherwise a typo is a parameter that silently does nothing."""
    comp = _component(tmp_path, """
        name = "mid"
        [params]
        temp = 30
        [[steps]]
        id = "b"
    """)
    with pytest.raises(ValueError, match="does not declare it"):
        template_mod.load(_write(tmp_path, f"""
            name = "x"
            [[steps]]
            include = "{comp}"
            [steps.params]
            tempe = 25
        """))


def test_a_prefix_renames_steps_params_and_their_references(tmp_path):
    """The same component twice in one protocol has to stay two things."""
    comp = _component(tmp_path, """
        name = "round"
        [params]
        temp = 30
        [[steps]]
        id = "grow"
        title = "Grow at {temp}"
        text = "Hold at {temp} C."
        checks = ["Was {temp} C"]
        [[steps]]
        id = "check"
        depends_on = ["grow"]
    """)
    tpl = template_mod.load(_write(tmp_path, f"""
        name = "x"
        [[steps]]
        include = "{comp}"
        prefix = "one_"
        [[steps]]
        include = "{comp}"
        prefix = "two_"
        [steps.params]
        temp = 37
    """))
    assert [s.id for s in tpl.steps] == ["one_grow", "one_check",
                                         "two_grow", "two_check"]
    assert tpl.params["one_temp"] == 30
    assert tpl.params["two_temp"] == 37
    assert tpl.step("one_grow").text == "Hold at {one_temp} C."
    assert tpl.step("two_grow").text == "Hold at {two_temp} C."
    assert tpl.step("two_grow").title == "Grow at {two_temp}"
    assert tpl.step("two_grow").checks == ["Was {two_temp} C"]
    assert tpl.step("two_check").depends_on == ["two_grow"]


def test_a_prefix_repoints_compute_args(tmp_path):
    comp = _component(tmp_path, """
        name = "r"
        [params]
        n = 3
        [[steps]]
        id = "calc"
        compute = "pandas:DataFrame"
        [steps.args]
        data = "$n"
    """)
    tpl = template_mod.load(_write(tmp_path, f"""
        name = "x"
        [[steps]]
        include = "{comp}"
        prefix = "a_"
    """))
    assert tpl.step("a_calc").args == {"data": "$a_n"}


def test_an_unprefixed_include_leaves_everything_alone(tmp_path):
    comp = _component(tmp_path, """
        name = "r"
        [params]
        n = 3
        [[steps]]
        id = "calc"
        text = "{n} of them"
    """)
    tpl = template_mod.load(_write(tmp_path, f"""
        name = "x"
        [[steps]]
        include = "{comp}"
    """))
    assert tpl.step("calc").text == "{n} of them"
    assert tpl.params["n"] == 3


def test_a_step_cannot_be_both_an_include_and_a_step(tmp_path):
    with pytest.raises(ValueError, match='both "include" and "id"'):
        template_mod.load(_write(tmp_path, """
            name = "x"
            [[steps]]
            include = "cloning_tail"
            id = "a"
        """))


def test_a_component_with_no_steps_is_rejected(tmp_path):
    comp = _component(tmp_path, 'name = "empty"\n')
    with pytest.raises(ValueError, match="has no steps"):
        template_mod.load(_write(tmp_path, f"""
            name = "x"
            [[steps]]
            include = "{comp}"
        """))


def test_a_component_colliding_with_a_template_step_is_rejected(tmp_path):
    """Validation runs after expansion, so the usual rules still apply."""
    comp = _component(tmp_path, """
        name = "mid"
        [[steps]]
        id = "a"
    """)
    with pytest.raises(ValueError, match="duplicate step id 'a'"):
        template_mod.load(_write(tmp_path, f"""
            name = "x"
            [[steps]]
            id = "a"
            [[steps]]
            include = "{comp}"
        """))


def test_a_component_brings_its_inputs(tmp_path):
    comp = _component(tmp_path, """
        name = "mid"
        [[inputs]]
        name = "sheet"
        optional = true
        [[steps]]
        id = "b"
    """)
    tpl = template_mod.load(_write(tmp_path, f"""
        name = "x"
        [[steps]]
        include = "{comp}"
    """))
    assert [i.name for i in tpl.inputs] == ["sheet"]
    assert tpl.inputs[0].optional is True


# --------------------------------------------------------------------------- #
# The bundled cloning component
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name", ["golden_gate", "gibson_assembly"])
def test_both_cloning_templates_share_the_tail(name):
    tpl = template_mod.load(name)
    assert tpl.components == ["cloning_tail"]
    assert [s.id for s in tpl.steps][-4:] == ["transform", "inoculate",
                                              "miniprep", "sequence"]
    assert tpl.step("transform").depends_on == ["assemble"]


def test_gibson_now_sequences_what_it_assembled():
    """It ended at the miniprep before the tail was shared, so an assembly
    was never verified."""
    tpl = template_mod.load("gibson_assembly")
    assert tpl.step("sequence").component == "cloning_tail"
    assert "junction" in tpl.step("sequence").text


@pytest.mark.parametrize("name,mL,temp,hours", [
    ("golden_gate", 3, 37, 20),
    ("gibson_assembly", 4, 30, 22),
])
def test_sharing_the_tail_kept_each_protocols_own_numbers(name, mL, temp, hours):
    """The component's defaults are the old values, and each template's own
    params still win: nobody's protocol changed by being refactored."""
    tpl = template_mod.load(name)
    assert tpl.params["culture_mL"] == mL
    assert tpl.params["culture_temp"] == temp
    assert tpl.params["culture_hours"] == hours


@pytest.mark.parametrize("name", template_mod.list_components())
def test_every_bundled_component_loads_on_its_own(name, tmp_path):
    """A component is not a protocol, but it still has to parse and validate
    when something includes it."""
    tpl = template_mod.load(_write(tmp_path, f"""
        name = "probe"
        [[steps]]
        include = "{name}"
    """))
    assert tpl.steps
    assert all(s.component == name for s in tpl.steps)


def test_an_include_param_can_alias_a_template_param(tmp_path):
    """Copying the value in would go stale the first time somebody deviated
    the original, so '$name' repoints instead."""
    comp = _component(tmp_path, """
        name = "mid"
        [params]
        n = 1
        [[steps]]
        id = "b"
        text = "{n} of them"
        compute = "pandas:DataFrame"
        [steps.args]
        data = "$n"
    """)
    tpl = template_mod.load(_write(tmp_path, f"""
        name = "x"
        [params]
        outer = 7
        [[steps]]
        include = "{comp}"
        [steps.params]
        n = "$outer"
    """))
    assert tpl.step("b").args == {"data": "$outer"}
    assert tpl.step("b").text == "{outer} of them"
    assert "n" not in tpl.params          # the template's own supplies it
    assert tpl.params["outer"] == 7


def test_an_alias_survives_a_prefix(tmp_path):
    """An alias points out of the component, so it is never prefixed."""
    comp = _component(tmp_path, """
        name = "mid"
        [params]
        n = 1
        keep = 2
        [[steps]]
        id = "b"
        text = "{n} and {keep}"
    """)
    tpl = template_mod.load(_write(tmp_path, f"""
        name = "x"
        [params]
        outer = 7
        [[steps]]
        include = "{comp}"
        prefix = "p_"
        [steps.params]
        n = "$outer"
    """))
    assert tpl.step("p_b").text == "{outer} and {p_keep}"
    assert tpl.params["p_keep"] == 2
    assert "p_n" not in tpl.params


# --------------------------------------------------------------------------- #
# The other bundled components
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name,follows,source", [
    ("transduction", "pellet", "the gDNA plate"),
    ("competitive_growth", "flow", "the kept aliquots"),
])
def test_both_screen_arms_share_the_gdna_handoff(name, follows, source):
    tpl = template_mod.load(name)
    assert "gdna_handoff" in tpl.components
    assert tpl.step("extract").depends_on == [follows]
    assert tpl.params["extract_from"] == source
    # both were already running these, and sharing did not change them
    assert tpl.params["gdna_reagent"] == "Lucigen QuickExtract"
    assert tpl.params["gdna_extract_uL"] == 50


def test_competitive_growth_kept_its_flow_export_check():
    """The handoff component drops it, so it moved to the step it is about."""
    tpl = template_mod.load("competitive_growth")
    assert any("day 0 reference" in c for c in tpl.step("flow").checks)


@pytest.mark.parametrize("name", ["pe3_transfection", "virus_production"])
def test_both_transfections_share_the_head(name):
    tpl = template_mod.load(name)
    assert "transfection_head" in tpl.components
    assert [s.id for s in tpl.steps][:2] == ["plate_plan", "seed"]
    assert tpl.step("seed").depends_on == ["plate_plan"]


def test_the_shared_head_keeps_each_confluency_target():
    """70-90% for a transfection and 95-99% for packaging is the whole of
    what differs, so it had better survive being shared."""
    assert template_mod.load("pe3_transfection").params["confluency"] == "70-90%"
    assert template_mod.load("virus_production").params["confluency"] == "95-99%"


def test_the_shared_head_keeps_each_seeding_density():
    tpl = template_mod.load("virus_production")
    assert tpl.params["seed_per_cm2"] == 125000
    pe3 = template_mod.load("pe3_transfection")
    assert pe3.params["seed_per_cm2"] != tpl.params["seed_per_cm2"]


def test_the_shared_head_keeps_the_protocol_specific_warning():
    """A shared step flattens the reasoning first; seed_note is where it
    survives."""
    virus = template_mod.load("virus_production")
    assert "{seed_note}" in virus.step("seed").text   # substituted at resolve
    assert "two-fold loss in titer" in virus.params["seed_note"]
    pe3 = template_mod.load("pe3_transfection")
    assert "Lipofection toxicity" in pe3.params["seed_note"]


def test_virus_counts_its_wells_from_the_transfer_plasmids():
    """conditions aliases transfers, so a run that changes one changes both."""
    tpl = template_mod.load("virus_production")
    assert tpl.step("plate_plan").args["conditions"] == "$transfers"


def test_pe3_leaves_its_condition_count_unknown():
    """It cannot be known until the epegRNA sheet is read, and 0 says so."""
    tpl = template_mod.load("pe3_transfection")
    assert tpl.params["conditions"] == 0
