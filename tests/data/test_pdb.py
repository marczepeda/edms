"""Tests for edms.data.pdb.

retrieve() downloads structure files from the RCSB PDB REST API via
urllib.request (NOT the `requests` library) -- every test mocks
urllib.request.urlopen and never touches the network. compute_residue_neighbors()
and compute_residue_contacts() are pure local structural-geometry computations
(via MDAnalysis) that never touch the network either; they're exercised
against tiny hand-built PDB files written to tmp_path.
"""
import urllib.error

import pandas as pd
import pytest

from edms.data import pdb


# ---------------------------------------------------------------------------
# retrieve()
# ---------------------------------------------------------------------------

def test_retrieve_downloads_expected_url_and_saves_file(tmp_path, mocker):
    fake_resp = mocker.MagicMock()
    fake_resp.read.return_value = b"ATOM some pdb content"
    fake_resp.__enter__.return_value = fake_resp
    mock_urlopen = mocker.patch("edms.data.pdb.urllib.request.urlopen", return_value=fake_resp)

    pdb.retrieve("1abc", suf=".pdb", dir=str(tmp_path), config=False)

    mock_urlopen.assert_called_once_with("https://files.rcsb.org/download/1abc.pdb")
    out_file = tmp_path / "1abc.pdb"
    assert out_file.exists()
    assert out_file.read_text() == "ATOM some pdb content"


def test_retrieve_lowercases_id_in_url_and_filename(tmp_path, mocker):
    fake_resp = mocker.MagicMock()
    fake_resp.read.return_value = b"content"
    fake_resp.__enter__.return_value = fake_resp
    mock_urlopen = mocker.patch("edms.data.pdb.urllib.request.urlopen", return_value=fake_resp)

    pdb.retrieve("1ABC", suf=".cif", dir=str(tmp_path), config=False)

    mock_urlopen.assert_called_once_with("https://files.rcsb.org/download/1abc.cif")
    assert (tmp_path / "1abc.cif").exists()


def test_retrieve_multiple_suffixes_downloads_each(tmp_path, mocker):
    fake_resp = mocker.MagicMock()
    fake_resp.read.return_value = b"content"
    fake_resp.__enter__.return_value = fake_resp
    mock_urlopen = mocker.patch("edms.data.pdb.urllib.request.urlopen", return_value=fake_resp)

    pdb.retrieve("1abc", suf=[".pdb", ".cif"], dir=str(tmp_path), config=False)

    urls = [c.args[0] for c in mock_urlopen.call_args_list]
    assert urls == [
        "https://files.rcsb.org/download/1abc.pdb",
        "https://files.rcsb.org/download/1abc.cif",
    ]
    assert (tmp_path / "1abc.pdb").exists()
    assert (tmp_path / "1abc.cif").exists()


def test_retrieve_respects_custom_base_url(tmp_path, mocker):
    fake_resp = mocker.MagicMock()
    fake_resp.read.return_value = b"content"
    fake_resp.__enter__.return_value = fake_resp
    mock_urlopen = mocker.patch("edms.data.pdb.urllib.request.urlopen", return_value=fake_resp)

    pdb.retrieve(
        "1abc", suf=".pdb", dir=str(tmp_path), config=False, base_url="https://example.org/api/"
    )

    mock_urlopen.assert_called_once_with("https://example.org/api/1abc.pdb")


@pytest.mark.parametrize("bad_suf", [".txt", ".fasta"])
def test_retrieve_invalid_suffix_raises_value_error(bad_suf):
    with pytest.raises(ValueError, match="must end with"):
        pdb.retrieve("1abc", suf=bad_suf)


def test_retrieve_http_error_raises_runtime_error_with_context(tmp_path, mocker):
    err = urllib.error.HTTPError(url="x", code=404, msg="Not Found", hdrs=None, fp=None)
    mocker.patch("edms.data.pdb.urllib.request.urlopen", side_effect=err)

    with pytest.raises(RuntimeError, match=r"HTTP 404"):
        pdb.retrieve("9xyz", suf=".pdb", dir=str(tmp_path), config=False)


def test_retrieve_url_error_raises_runtime_error_with_context(tmp_path, mocker):
    err = urllib.error.URLError("connection refused")
    mocker.patch("edms.data.pdb.urllib.request.urlopen", side_effect=err)

    with pytest.raises(RuntimeError, match="Failed to reach RCSB PDB"):
        pdb.retrieve("1abc", suf=".pdb", dir=str(tmp_path), config=False)


def test_retrieve_config_true_writes_only_within_redirected_home(home_dir, mocker):
    fake_resp = mocker.MagicMock()
    fake_resp.read.return_value = b"content"
    fake_resp.__enter__.return_value = fake_resp
    mocker.patch("edms.data.pdb.urllib.request.urlopen", return_value=fake_resp)

    pdb.retrieve("1abc", suf=".pdb", dir=None, config=True)

    expected = home_dir / ".config" / "edms" / "PDB" / "1abc.pdb"
    assert expected.exists()


# ---------------------------------------------------------------------------
# compute_residue_neighbors() / compute_residue_contacts()
# (pure local geometry computation -- no network involved)
# ---------------------------------------------------------------------------

TWO_RESIDUE_PDB = """\
ATOM      1  N   ALA A   1      -1.000   0.000   0.000  1.00 50.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00 50.00           C
ATOM      3  C   ALA A   1       1.000   0.000   0.000  1.00 50.00           C
ATOM      4  O   ALA A   1       1.500   1.000   0.000  1.00 50.00           O
ATOM      5  N   GLY A   2       3.000   0.000   0.000  1.00 80.00           N
ATOM      6  CA  GLY A   2       4.000   0.000   0.000  1.00 80.00           C
ATOM      7  C   GLY A   2       5.000   0.000   0.000  1.00 80.00           C
ATOM      8  O   GLY A   2       5.500   1.000   0.000  1.00 80.00           O
END
"""

THREE_RESIDUE_CONTACT_PDB = """\
ATOM      1  N   ALA A   1      -1.000   0.000   0.000  1.00 50.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00 50.00           C
ATOM      3  C   ALA A   1       1.000   0.000   0.000  1.00 50.00           C
ATOM      4  O   ALA A   1       1.500   1.000   0.000  1.00 50.00           O
ATOM      5  CB  ALA A   1       0.000   1.500   0.000  1.00 50.00           C
ATOM      6  N   GLY A   2     100.000   0.000   0.000  1.00 80.00           N
ATOM      7  CA  GLY A   2     101.000   0.000   0.000  1.00 80.00           C
ATOM      8  C   GLY A   2     102.000   0.000   0.000  1.00 80.00           C
ATOM      9  O   GLY A   2     102.500   1.000   0.000  1.00 80.00           O
ATOM     10  N   SER A   3       2.000   0.000   0.000  1.00 60.00           N
ATOM     11  CA  SER A   3       3.000   0.000   0.000  1.00 60.00           C
ATOM     12  C   SER A   3       4.000   0.000   0.000  1.00 60.00           C
ATOM     13  O   SER A   3       4.500   1.000   0.000  1.00 60.00           O
ATOM     14  CB  SER A   3       3.000   1.700   0.000  1.00 60.00           C
END
"""


@pytest.fixture
def two_residue_pdb(tmp_path):
    p = tmp_path / "two_res.pdb"
    p.write_text(TWO_RESIDUE_PDB)
    return str(p)


@pytest.fixture
def three_residue_pdb(tmp_path):
    p = tmp_path / "three_res.pdb"
    p.write_text(THREE_RESIDUE_CONTACT_PDB)
    return str(p)


def test_compute_residue_neighbors_finds_pairs_within_cutoff(two_residue_pdb):
    df = pdb.compute_residue_neighbors(two_residue_pdb, cutoff=10.0)
    assert len(df) == 2  # (res1, res2) and (res2, res1)
    assert set(df["query_resid"]) == {1, 2}
    assert (df["distance"] == 4.0).all()  # CA-CA centers of geometry are 4 A apart


def test_compute_residue_neighbors_beyond_cutoff_returns_empty(two_residue_pdb):
    df = pdb.compute_residue_neighbors(two_residue_pdb, cutoff=1.0)
    assert df.empty


def test_compute_residue_neighbors_confidence_weighting_uses_bfactor(two_residue_pdb):
    # B-factors: residue1 atoms=50, residue2 atoms=80 -> conf = 0.5, 0.8
    # effective_distance = distance / (conf_i * conf_j) = 4.0 / (0.5*0.8) = 10.0
    df = pdb.compute_residue_neighbors(two_residue_pdb, cutoff=10.0, weigh_by_confidence=True)
    assert df["effective_distance"].round(2).tolist() == [10.0, 10.0]


def test_compute_residue_neighbors_empty_selection_warns_and_returns_empty_df(two_residue_pdb, capsys):
    df = pdb.compute_residue_neighbors(two_residue_pdb, selection="nucleic")
    assert df.empty
    assert list(df.columns) == [
        "query_index", "query_resid", "query_resname", "query_chain",
        "partner_index", "partner_resid", "partner_resname", "partner_chain",
        "distance", "effective_distance",
    ]
    captured = capsys.readouterr()
    assert "Warning" in captured.out


def test_compute_residue_neighbors_rejects_unsupported_format(tmp_path):
    bad_file = tmp_path / "structure.xyz"
    bad_file.write_text("not a structure")
    with pytest.raises(ValueError, match="Unsupported file format"):
        pdb.compute_residue_neighbors(str(bad_file))


def test_compute_residue_contacts_finds_min_heavy_atom_distance(three_residue_pdb):
    df = pdb.compute_residue_contacts(three_residue_pdb, cutoff=6.0)
    assert not df.empty
    # ALA(1) CB and SER(3) CB are the closest non-backbone atoms (~1.66 A)
    pair = df[(df["query_resid"] == 1) & (df["partner_resid"] == 3)].iloc[0]
    assert pair["min_atom_distance"] == pytest.approx(1.6553, abs=1e-3)


def test_compute_residue_contacts_skips_adjacent_residues_same_chain(three_residue_pdb):
    df = pdb.compute_residue_contacts(three_residue_pdb, cutoff=200.0)
    # GLY(2) is adjacent to both ALA(1) and SER(3)... but resid 1 & 3 are NOT
    # adjacent (diff=2) and ARE reported; resid pairs with |diff|==1 must never appear.
    diffs = (df["query_resid"] - df["partner_resid"]).abs()
    assert not (diffs == 1).any()


def test_compute_residue_contacts_backbone_exclusion_toggle(three_residue_pdb):
    # With backbone-backbone pairs allowed, the very close backbone O/N atoms
    # of ALA and SER should be found (and be at least as close as heavy CB-CB).
    df_excluded = pdb.compute_residue_contacts(
        three_residue_pdb, cutoff=6.0, exclude_backbone_backbone=True
    )
    df_included = pdb.compute_residue_contacts(
        three_residue_pdb, cutoff=6.0, exclude_backbone_backbone=False
    )
    assert not df_excluded.empty
    assert not df_included.empty
    excluded_min = df_excluded["min_atom_distance"].min()
    included_min = df_included["min_atom_distance"].min()
    assert included_min <= excluded_min
