"""
Smoke tests for edms.gen.cli argument-parsing wiring.

Per task scope, we don't test argparse wiring in depth here - just confirm
the gen subparsers can be attached to a top-level parser without crashing,
and that a couple of representative subcommands parse valid arguments.
"""
import argparse
import sys

import pytest

from edms.gen import cli


def _build_parser():
    parser = argparse.ArgumentParser(prog="edms")
    subparsers = parser.add_subparsers(dest="command")
    # NOTE: add_subparser()'s own default is formatter_class=None, but that
    # default is unusable (argparse crashes formatting --help without a
    # real formatter class) - production code (edms/main.py) always passes
    # a concrete formatter (a RichHelpFormatter subclass), so we do the
    # same here rather than exercise the unusable default.
    cli.add_subparser(subparsers, argparse.HelpFormatter)
    return parser


def test_add_subparser_builds_without_error():
    parser = _build_parser()
    assert parser is not None


def test_plot_subcommand_help_does_not_crash(capsys):
    parser = _build_parser()
    # add_subparser() attaches 'plot', 'stat', 'io', 'com', 'html'
    # subcommands directly (there is no top-level 'gen' subcommand).
    # argparse's --help exits via SystemExit(0); just confirm it doesn't
    # raise anything else and prints usage info.
    try:
        parser.parse_args(["plot", "--help"])
    except SystemExit as e:
        assert e.code == 0
    captured = capsys.readouterr()
    assert "usage" in captured.out.lower()


def test_pipeline_help_documents_every_reserved_column(monkeypatch, capsys):
    """The `edms io pipeline -h` column reference must not drift from the parser.

    add_subparser() prints the reference and exits when it sees this exact
    command in sys.argv (the same trick edms pe designer uses, because
    RichHelpFormatter reflows hand-laid-out block text).
    """
    from edms.gen import io

    monkeypatch.setattr(sys, "argv", ["edms", "io", "pipeline", "-h"])
    with pytest.raises(SystemExit):
        _build_parser()

    out = capsys.readouterr().out
    for column in io.PIPELINE_RESERVED:
        assert column in out, f"'{column}' is reserved but absent from `edms io pipeline -h`"
    # The '#' conventions and the fall-through flag rule are documented too
    assert "#notes" in out
    assert "--fastq_dir" in out


def test_pipeline_help_not_triggered_by_other_commands(monkeypatch, capsys):
    # A path that merely contains "pipeline" must not hijack another command's help
    monkeypatch.setattr(sys, "argv", ["edms", "fastq", "trim", "-i", "pipeline.csv", "-h"])
    _build_parser()  # must return normally, not SystemExit
    assert "column reference" not in capsys.readouterr().out
