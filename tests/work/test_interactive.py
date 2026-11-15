"""
Tests for edms.work.interactive

The payload is what the page recalculates from, so the things worth pinning
are that every number in it came out of the Python objects rather than being
written down a second time, that a compute target calc.js has not implemented
is reported as static instead of silently stale, and that nothing here can
write to a run.
"""
import json

import pandas as pd
import pytest

from edms.work import interactive, plates, render


# --------------------------------------------------------------------------- #
# Assets
# --------------------------------------------------------------------------- #
def test_bundled_assets_are_present():
    for name in ("calc.js", "plate.js", "app.js"):
        assert interactive.asset(name).strip()


def test_asset_rejects_an_unknown_name():
    with pytest.raises(FileNotFoundError, match="No bundled asset"):
        interactive.asset("nope.js")


def test_asset_cannot_be_used_to_read_outside_the_package():
    with pytest.raises(FileNotFoundError):
        interactive.asset("../plates.py")


def test_copying_a_deviate_command_does_not_assume_a_secure_context():
    """A protocol is read on a phone, where the clipboard API is not there.

    navigator.clipboard is undefined rather than failing outside a secure
    context, so reaching straight for .writeText throws on the property read
    and takes any rejection handler with it. localhost counts as secure, which
    is why this was invisible on the machine that rendered the file and broken
    on every route to a phone: file:// and plain http from the laptop both.
    """
    src = interactive.asset("app.js")
    guard = "navigator.clipboard && navigator.clipboard.writeText"
    # The API is called once, and only after something checked it is there
    assert src.count("navigator.clipboard.writeText(") == 1
    assert 0 <= src.index(guard) < src.index("navigator.clipboard.writeText(")
    # ...and there is a way through that needs no API at all, because on a
    # phone there will not be one
    assert "document.execCommand('copy')" in src
    assert "edms-copy-text" in src


# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #
def test_constants_come_from_the_python_objects():
    """No vessel table is written down twice; the browser gets this one."""
    c = interactive.constants()
    assert c["vessels"] == {k: dict(v) for k, v in plates._VESSELS.items()}
    assert c["SEED_PER_CM2"] == plates.SEED_PER_CM2
    assert c["VIRUS_SEED_PER_CM2"] == plates.VIRUS_SEED_PER_CM2
    assert c["DILUTION_STEPS"] == list(plates.DILUTION_STEPS)
    assert c["VIRUS_ANCHOR"] == plates.VIRUS_ANCHOR


def _code_only(js):
    """calc.js without its comments, which discuss numbers by name."""
    import re

    js = re.sub(r"/\*.*?\*/", "", js, flags=re.S)
    return re.sub(r"(?m)^\s*//.*$", "", js)


def test_calc_js_has_no_vessel_table_of_its_own():
    """The areas are what would drift silently, so they are read from
    EDMS.consts and never written here."""
    js = interactive.asset("calc.js")
    assert "EDMS.consts.vessels" in js
    code = _code_only(js)
    # the distinctive areas; 1.1 and 3.5 are left out because they collide
    # with ordinary defaults like mm_x
    for area in (0.32, 1.9, 9.6, 21.5, 56.7, 145.0):
        assert str(area) not in code, f"{area} is written into calc.js"


def test_constants_are_json_serialisable():
    json.dumps(interactive.constants(), allow_nan=False)


# --------------------------------------------------------------------------- #
# live_targets
# --------------------------------------------------------------------------- #
def test_live_targets_is_read_from_the_file():
    """The set is parsed out of calc.js rather than listed anywhere, so it
    cannot claim a target the file does not define. Every bundled calculator
    is ported now, which makes the absent case the one worth pinning: a name
    nothing implements must not appear."""
    live = interactive.live_targets()
    assert "edms.work.plates:virus" in live
    assert "edms.work.plates:geometry" in live
    assert "edms.bio.ngs:pcrs" in live
    assert "edms.bio.clone:invented" not in live
    for target in live:
        assert "C['%s']" % target in interactive.asset("calc.js")


def test_every_live_target_is_a_real_python_function():
    import importlib

    for target in interactive.live_targets():
        module, _, func = target.partition(":")
        assert callable(getattr(importlib.import_module(module), func))


# --------------------------------------------------------------------------- #
# Payload
# --------------------------------------------------------------------------- #
def test_table_data_keeps_column_order():
    df = pd.DataFrame({"b": [1], "a": [2]})
    assert interactive.table_data(df)["columns"] == ["b", "a"]


def test_table_data_turns_nan_into_null():
    df = pd.DataFrame({"x": [float("nan")]})
    assert interactive.table_data(df)["rows"][0]["x"] is None


def test_payload_is_json_serialisable(virus_run):
    json.dumps(interactive.payload(*virus_run), allow_nan=False)


def test_payload_marks_unported_steps_static(virus_run):
    protocol, inputs = virus_run
    data = interactive.payload(protocol, inputs)
    for step in data["steps"]:
        assert step["live"] == (step["compute"] in interactive.live_targets())


def test_payload_carries_the_unresolved_args(virus_run):
    """The page resolves $param itself, so a change repoints the same way
    plan.resolve would."""
    protocol, inputs = virus_run
    data = interactive.payload(protocol, inputs)
    plan_step = [s for s in data["steps"] if s["id"] == "transfection_plan"][0]
    assert plan_step["args"]["plate_format"] == "$plate_format"
    assert plan_step["args"]["plasmids"] == "@plasmids"


def test_payload_includes_the_inputs_so_they_can_be_edited(virus_run):
    protocol, inputs = virus_run
    data = interactive.payload(protocol, inputs)
    assert "plasmids" in data["inputs"]
    assert "ng/uL" in data["inputs"]["plasmids"]["columns"]


def test_payload_carries_the_row_cap_the_static_render_uses(virus_run):
    """Otherwise the page agrees with itself only until you touch something.

    A 352-well plate map is written out cut to MAX_ROWS with a note pointing
    at the .xlsx; before this, the first recompute replaced it with every row
    and dropped the note, so one table said two different things about how
    much of it you were looking at.
    """
    assert interactive.payload(*virus_run)["max_rows"] == render.MAX_ROWS
    assert "D.max_rows" in interactive.asset("app.js")


def test_payload_carries_the_plate_drawing_cap(virus_run):
    """Fifteen plates of 48-well is a wall of pictures, not a glance."""
    assert interactive.payload(*virus_run)["max_plates"] == render.MAX_PLATES
    assert "D.max_plates" in interactive.asset("app.js")


def test_a_layout_step_gets_somewhere_to_draw_its_plates():
    """render.py leaves the host empty and app.js fills it, so a step whose
    plates are never drawn is a missing div rather than a visible error."""
    import inspect

    assert "drawLayout" in interactive.asset("app.js")
    assert ":layout" in inspect.getsource(render.to_html)


def test_payload_offers_every_format(virus_run):
    assert interactive.payload(*virus_run)["formats"] == plates.names()


# --------------------------------------------------------------------------- #
# Rendering
# --------------------------------------------------------------------------- #
def test_html_embeds_the_scripts_and_data(virus_run):
    html = render.to_html(virus_run[0])
    assert 'id="edms-data"' in html
    assert "EDMS.calc" in html
    assert 'id="edms-controls"' in html


def test_html_can_be_rendered_flat(virus_run):
    """The record does not have to carry a calculator."""
    html = render.to_html(virus_run[0], interactive=False)
    assert 'id="edms-data"' not in html
    assert 'id="edms-controls"' not in html
    assert "EDMS.calc" not in html
    assert virus_run[0].run_id in html


def test_html_tables_are_addressable_by_step(virus_run):
    html = render.to_html(virus_run[0])
    assert 'data-step="transfection_plan"' in html
    assert 'data-table="transfection_mix"' in html


def test_html_marks_a_geometry_step_for_a_plate_map(virus_run):
    assert 'data-plate="plate_plan"' in render.to_html(virus_run[0])


def test_html_stripes_its_tables(virus_run):
    assert "tbody tr:nth-child(even)" in render.to_html(virus_run[0])


def test_flat_html_still_stripes_its_tables(virus_run):
    """Striping is presentation, not interactivity, so it is in both."""
    assert "tbody tr:nth-child(even)" in render.to_html(virus_run[0],
                                                        interactive=False)


def test_html_closes_the_data_script_safely(virus_run):
    """A '</script>' inside the JSON would end the block early."""
    html = render.to_html(virus_run[0])
    start = html.index('id="edms-data"')
    block = html[start:html.index("</script>", start)]
    assert "</" not in block.split(">", 1)[1]


def test_rendering_does_not_write_to_the_run(virus_run, tmp_path):
    """A protocol opened in a browser is a record; the page is a scratchpad."""
    protocol, _ = virus_run
    before = (protocol.run_dir / "run.json").read_text()
    render.to_html(protocol)
    assert (protocol.run_dir / "run.json").read_text() == before
