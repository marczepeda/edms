"""Tests for edms.data.dssp.

dssp.retrieve() shells out to the `mkdssp` binary rather than making HTTP
calls; subprocess.run is always mocked so tests never require mkdssp to be
installed. parse_segments()/ssa_str()/dssp_pdb_id() are pure file parsers
tested against small hand-built legacy-DSSP-format text (no real dssp
output needed). All tests that exercise the ~/.config/edms/* default cache
paths use the `home_dir` fixture to redirect HOME into tmp_path.
"""
import os

import pandas as pd
import pytest

from edms.data import dssp


def _dssp_line(resnum, chain, aa, ss):
    """Build one legacy-DSSP-format residue line with correct fixed columns.

    resnum -> line[5:10], chain -> line[11], aa -> line[13], ss -> line[16]
    """
    line = list(" " * 40)
    s = str(resnum)
    line[10 - len(s):10] = list(s)
    line[11] = chain
    line[13] = aa
    line[16] = ss
    return "".join(line)


HEADER_LINE = "HEADER    PLANT SEED PROTEIN                      30-APR-81   1CRN"
DATA_HEADER = "  #  RESIDUE AA STRUCTURE"


def _write_dssp_file(tmp_path, lines, name="test.dssp"):
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return str(path)


# ---------------------------------------------------------------------------
# retrieve() -- shells out to mkdssp via subprocess.run
# ---------------------------------------------------------------------------

def test_retrieve_from_filepath_builds_expected_mkdssp_command(tmp_path, home_dir, mocker):
    mock_run = mocker.patch("edms.data.dssp.subprocess.run")

    cif_path = tmp_path / "mystruct.cif"
    cif_path.write_text("dummy cif content")

    dssp.retrieve(str(cif_path), dir=str(tmp_path))

    mock_run.assert_called_once()
    command = mock_run.call_args.args[0]
    assert command.startswith("mkdssp ")
    assert str(cif_path) in command
    assert str(tmp_path / "mystruct.dssp") in command
    assert mock_run.call_args.kwargs.get("shell") is True


def test_retrieve_from_config_pdb_dir_finds_matching_file(home_dir, mocker):
    mock_run = mocker.patch("edms.data.dssp.subprocess.run")

    pdb_dir = home_dir / ".config" / "edms" / "PDB"
    pdb_dir.mkdir(parents=True)
    (pdb_dir / "1abc.cif").write_text("dummy")

    dssp.retrieve("1ABC")

    mock_run.assert_called_once()
    command = mock_run.call_args.args[0]
    assert "1abc.cif" in command
    assert "1ABC.dssp" in command


def test_retrieve_missing_config_dir_raises_file_not_found(home_dir, mocker):
    mock_run = mocker.patch("edms.data.dssp.subprocess.run")
    # No ~/.config/edms/PDB directory exists at all -> os.listdir raises -> caught -> re-raised
    with pytest.raises(FileNotFoundError):
        dssp.retrieve("1ABC")
    mock_run.assert_not_called()


def test_retrieve_no_matching_file_in_existing_config_dir_raises_file_not_found(home_dir, mocker):
    mock_run = mocker.patch("edms.data.dssp.subprocess.run")

    pdb_dir = home_dir / ".config" / "edms" / "PDB"
    pdb_dir.mkdir(parents=True)
    (pdb_dir / "unrelated.cif").write_text("dummy")  # exists, but does not match "zzzz"

    with pytest.raises(FileNotFoundError):
        dssp.retrieve("zzzz")

    assert not mock_run.called


# ---------------------------------------------------------------------------
# parse_segments()
# ---------------------------------------------------------------------------

def test_parse_segments_groups_contiguous_runs_by_code(tmp_path):
    lines = [
        HEADER_LINE,
        DATA_HEADER,
        _dssp_line(1, "A", "M", "H"),
        _dssp_line(2, "A", "A", "H"),
        _dssp_line(3, "A", "C", "H"),
        _dssp_line(4, "A", "D", "E"),
        _dssp_line(5, "A", "E", "E"),
        _dssp_line(6, "A", "F", " "),
        _dssp_line(1, "B", "X", "H"),  # different chain, must be excluded
    ]
    path = _write_dssp_file(tmp_path, lines)

    df = dssp.parse_segments(path, chain_id="A")

    assert list(df["ss_code"]) == ["H", "E", " "]
    assert list(df["start"]) == [1, 4, 6]
    assert list(df["end"]) == [3, 5, 6]
    assert list(df["ss_description"]) == ["α-helix", "β-strand", "coil"]
    assert list(df["ss_color"]) == ["turquoise", "gold", "darkgray"]
    assert (df["chain"] == "A").all()


def test_parse_segments_chain_break_splits_segment(tmp_path):
    # Non-consecutive resnum with the same ss_code must start a new segment
    lines = [
        DATA_HEADER,
        _dssp_line(1, "A", "M", "H"),
        _dssp_line(2, "A", "A", "H"),
        _dssp_line(5, "A", "C", "H"),  # gap: 3,4 missing -> new segment
    ]
    path = _write_dssp_file(tmp_path, lines)

    df = dssp.parse_segments(path, chain_id="A")

    assert len(df) == 2
    assert list(df["start"]) == [1, 5]
    assert list(df["end"]) == [2, 5]


def test_parse_segments_unknown_code_uses_unknown_color(tmp_path):
    lines = [DATA_HEADER, _dssp_line(1, "A", "M", "Z")]
    path = _write_dssp_file(tmp_path, lines)

    df = dssp.parse_segments(path, chain_id="A", unknown_color="pink")

    assert df.iloc[0]["ss_color"] == "pink"
    assert df.iloc[0]["ss_description"] == "Z"
    assert df.iloc[0]["ss_name"] == "Z #1"


def test_parse_segments_no_records_raises_value_error(tmp_path):
    path = _write_dssp_file(tmp_path, ["some unrelated header line"])
    with pytest.raises(ValueError, match="No DSSP residue records"):
        dssp.parse_segments(path, chain_id="A")


def test_parse_segments_missing_chain_raises_value_error(tmp_path):
    lines = [DATA_HEADER, _dssp_line(1, "A", "M", "H")]
    path = _write_dssp_file(tmp_path, lines)
    with pytest.raises(ValueError, match="No residues found for chain"):
        dssp.parse_segments(path, chain_id="Z")


# ---------------------------------------------------------------------------
# ssa_str()
# ---------------------------------------------------------------------------

def test_ssa_str_extracts_codes_for_chain(tmp_path):
    lines = [
        DATA_HEADER,
        _dssp_line(1, "A", "M", "H"),
        _dssp_line(2, "A", "A", "H"),
        _dssp_line(3, "A", "D", "E"),
        _dssp_line(1, "B", "X", "T"),
    ]
    path = _write_dssp_file(tmp_path, lines)
    assert dssp.ssa_str(path, "A") == "HHE"
    assert dssp.ssa_str(path, "B") == "T"


def test_ssa_str_missing_header_raises_value_error(tmp_path):
    path = _write_dssp_file(tmp_path, ["no header here"])
    with pytest.raises(ValueError, match="Could not find DSSP residue table"):
        dssp.ssa_str(path, "A")


# ---------------------------------------------------------------------------
# dssp_pdb_id()
# ---------------------------------------------------------------------------

def test_dssp_pdb_id_extracts_from_fixed_width_header(tmp_path):
    path = _write_dssp_file(tmp_path, [HEADER_LINE, DATA_HEADER])
    assert dssp.dssp_pdb_id(path) == "1CRN"


def test_dssp_pdb_id_returns_none_when_no_header(tmp_path):
    path = _write_dssp_file(tmp_path, [DATA_HEADER])
    assert dssp.dssp_pdb_id(path) is None


# ---------------------------------------------------------------------------
# pymol_color_defs_from_ss_map()
# ---------------------------------------------------------------------------

def test_pymol_color_defs_from_ss_map_converts_all_entries():
    lines = dssp.pymol_color_defs_from_ss_map(dssp.ss_map)
    assert len(lines) == len(dssp.ss_map)
    assert any(line.startswith("set_color dssp_coil,") for line in lines)
    assert any(line.startswith("set_color dssp_H,") for line in lines)


# ---------------------------------------------------------------------------
# pymol_ssa()
# ---------------------------------------------------------------------------

def test_pymol_ssa_writes_script_with_fetch_for_pdb_id(tmp_path, mocker):
    mock_run = mocker.patch("edms.data.dssp.subprocess.run")

    lines = [
        DATA_HEADER,
        _dssp_line(1, "A", "M", "H"),
        _dssp_line(2, "A", "A", "H"),
        _dssp_line(3, "A", "D", "E"),
    ]
    dssp_path = _write_dssp_file(tmp_path, lines)

    dssp.pymol_ssa(
        dssp_file=dssp_path,
        chain_id="A",
        dir=str(tmp_path),
        file="script.pml",
        pdb_id_or_filename="1ABC",
        execute=False,
    )

    out_file = tmp_path / "script.pml"
    assert out_file.exists()
    content = out_file.read_text()
    assert "fetch 1ABC, prot" in content
    assert "color dssp_H, (prot and chain A and resi 1-2)" in content
    assert "color dssp_E, (prot and chain A and resi 3-3)" in content
    mock_run.assert_not_called()  # execute=False


def test_pymol_ssa_execute_true_invokes_subprocess(tmp_path, mocker):
    mock_run = mocker.patch("edms.data.dssp.subprocess.run")
    lines = [DATA_HEADER, _dssp_line(1, "A", "M", "H")]
    dssp_path = _write_dssp_file(tmp_path, lines)

    dssp.pymol_ssa(
        dssp_file=dssp_path,
        chain_id="A",
        dir=str(tmp_path),
        file="script.pml",
        pdb_id_or_filename="1ABC",
        execute=True,
    )

    mock_run.assert_called_once()
    assert "pymol" in mock_run.call_args.args[0]
    assert str(tmp_path / "script.pml") in mock_run.call_args.args[0]


def test_pymol_ssa_missing_chain_raises_value_error(tmp_path, mocker):
    # parse_segments() itself raises for a missing chain before pymol_ssa's
    # own (dead-code) "No DSSP segments found" empty-check is ever reached.
    mocker.patch("edms.data.dssp.subprocess.run")
    lines = [DATA_HEADER, _dssp_line(1, "A", "M", "H")]
    dssp_path = _write_dssp_file(tmp_path, lines)

    with pytest.raises(ValueError, match="No residues found for chain"):
        dssp.pymol_ssa(
            dssp_file=dssp_path,
            chain_id="Z",
            dir=str(tmp_path),
            file="script.pml",
        )
