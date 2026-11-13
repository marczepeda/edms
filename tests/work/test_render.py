"""
Tests for edms.work.render

Covers both output formats: that status, notes, deviations and computed tables
all reach the page, that large tables are truncated rather than dumped whole,
and that the HTML is self-contained and theme-aware.
"""
import pandas as pd
import pytest

from edms.work import plan as plan_mod
from edms.work import render as render_mod
from edms.work import run as run_mod


@pytest.fixture
def protocol(started_run):
    run_dir, _ = started_run
    run_mod.done(run_dir, "prep", note="went fine")
    run_mod.deviate(run_dir, "prep", param="scale", to=9, why="short on cells")
    return plan_mod.resolve(run_dir)


# --------------------------------------------------------------------------- #
# HTML
# --------------------------------------------------------------------------- #
def test_html_is_a_standalone_document(protocol):
    html = render_mod.to_html(protocol)
    assert html.startswith("<!doctype html>")
    assert "<style>" in html          # CSS inlined, no external stylesheet
    assert "<link" not in html
    assert "<script" not in html


def test_html_includes_every_step(protocol):
    html = render_mod.to_html(protocol)
    for step in protocol.steps:
        assert step.title in html


def test_html_shows_status_notes_and_deviations(protocol):
    html = render_mod.to_html(protocol)
    assert 'class="step done"' in html
    assert "went fine" in html
    assert "short on cells" in html
    assert "scale" in html


def test_html_renders_computed_tables(protocol):
    html = render_mod.to_html(protocol)
    assert "<table>" in html
    assert "built" in html            # the output table's name


def test_html_escapes_values_from_the_ledger(started_run):
    run_dir, _ = started_run
    run_mod.note(run_dir, "prep", text="<script>alert(1)</script>")
    html = render_mod.to_html(plan_mod.resolve(run_dir))
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;" in html


def test_html_supports_both_colour_schemes(protocol):
    html = render_mod.to_html(protocol)
    assert "prefers-color-scheme: dark" in html
    assert '[data-theme="dark"]' in html


def test_html_theme_can_be_pinned(protocol):
    assert 'data-theme="dark"' in render_mod.to_html(protocol, theme="dark")
    assert "data-theme" not in render_mod.to_html(protocol).split("<head>")[0]


def test_html_truncates_a_long_table(protocol):
    big = pd.DataFrame({"n": range(render_mod.MAX_ROWS + 25)})
    protocol.steps[1].tables = {"big": big}
    html = render_mod.to_html(protocol)
    assert "25 more row(s)" in html
    assert html.count("<tr>") < len(big)


def test_html_reports_a_failed_calculation(protocol):
    protocol.steps[1].error = "ValueError: bad input"
    assert "Could not compute" in render_mod.to_html(protocol)


def test_html_reports_a_step_waiting_on_an_optional_input(protocol):
    html = render_mod.to_html(protocol)
    assert "Not calculated" in html


# --------------------------------------------------------------------------- #
# Markdown
# --------------------------------------------------------------------------- #
def test_markdown_has_a_heading_per_step(protocol):
    md = render_mod.to_markdown(protocol)
    for step in protocol.steps:
        assert step.title in md
    assert md.startswith("# ")


def test_markdown_marks_status_and_checks(protocol):
    md = render_mod.to_markdown(protocol)
    assert "## [x] Prepare" in md
    assert "- [ ] Bench is clear" in md


def test_markdown_renders_tables_without_tabulate(protocol):
    md = render_mod.to_markdown(protocol)
    assert "| n |" in md
    assert "| --- |" in md


def test_markdown_escapes_pipes_inside_cells(protocol):
    protocol.steps[1].tables = {"t": pd.DataFrame({"a": ["x|y"]})}
    assert r"x\|y" in render_mod.to_markdown(protocol)


def test_markdown_truncates_a_long_table(protocol):
    protocol.steps[1].tables = {"t": pd.DataFrame({"n": range(render_mod.MAX_ROWS + 5)})}
    assert "5 more row(s)" in render_mod.to_markdown(protocol)


def test_markdown_carries_notes_and_deviations(protocol):
    md = render_mod.to_markdown(protocol)
    assert "went fine" in md
    assert "Deviation: scale 5 -> 9" in md


def test_markdown_uses_no_html_entities(protocol):
    # It gets pasted into a notebook entry, so it has to read as plain text
    md = render_mod.to_markdown(protocol)
    assert "&middot;" not in md
    assert "·" in md


# --------------------------------------------------------------------------- #
# write()
# --------------------------------------------------------------------------- #
def test_write_defaults_into_the_run_directory(protocol):
    out = render_mod.write(protocol)
    assert out == protocol.run_dir / "protocol.html"
    assert out.read_text().startswith("<!doctype html>")


def test_write_markdown_uses_the_md_extension(protocol):
    assert render_mod.write(protocol, format="md").name == "protocol.md"


def test_write_accepts_an_explicit_path(protocol, tmp_path):
    out = render_mod.write(protocol, file=tmp_path / "nested" / "p.html")
    assert out.exists()


def test_write_rejects_an_unknown_format(protocol):
    with pytest.raises(ValueError, match="'html' or 'md'"):
        render_mod.write(protocol, format="pdf")
