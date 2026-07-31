"""Tests for edms.data.uniprot.

Only retrieve() performs an HTTP call (via urllib.request against the
UniProtKB REST flat-file endpoint); it is fully mocked in every test below.
The remainder of this module is pure parsing (flat-file / Feature-viewer
JSON) and matplotlib plotting; parsing functions are tested against small
hand-built inputs, and plotting helpers get light smoke tests under the
Agg backend configured in tests/conftest.py.
"""
import json
import urllib.error

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from edms.data import uniprot as up


# ---------------------------------------------------------------------------
# retrieve()
# ---------------------------------------------------------------------------

def test_retrieve_downloads_expected_url_and_saves_flat_file(tmp_path, mocker):
    fake_resp = mocker.MagicMock()
    fake_resp.read.return_value = b"ID   FOXA1_HUMAN  Reviewed; 472 AA."
    fake_resp.__enter__.return_value = fake_resp
    mock_urlopen = mocker.patch("edms.data.uniprot.urllib.request.urlopen", return_value=fake_resp)

    up.retrieve("P55317", dir=str(tmp_path), config=False)

    mock_urlopen.assert_called_once_with("https://rest.uniprot.org/uniprotkb/P55317.txt")
    out_file = tmp_path / "P55317.txt"
    assert out_file.exists()
    assert out_file.read_text() == "ID   FOXA1_HUMAN  Reviewed; 472 AA."


def test_retrieve_respects_custom_base_url(tmp_path, mocker):
    fake_resp = mocker.MagicMock()
    fake_resp.read.return_value = b"content"
    fake_resp.__enter__.return_value = fake_resp
    mock_urlopen = mocker.patch("edms.data.uniprot.urllib.request.urlopen", return_value=fake_resp)

    up.retrieve("P55317", dir=str(tmp_path), config=False, base_url="https://example.org/api/")

    mock_urlopen.assert_called_once_with("https://example.org/api/P55317.txt")


def test_retrieve_http_error_raises_runtime_error_with_context(tmp_path, mocker):
    err = urllib.error.HTTPError(url="x", code=404, msg="Not Found", hdrs=None, fp=None)
    mocker.patch("edms.data.uniprot.urllib.request.urlopen", side_effect=err)

    with pytest.raises(RuntimeError, match=r"HTTP 404"):
        up.retrieve("BADACC", dir=str(tmp_path), config=False)


def test_retrieve_url_error_raises_runtime_error_with_context(tmp_path, mocker):
    err = urllib.error.URLError("connection refused")
    mocker.patch("edms.data.uniprot.urllib.request.urlopen", side_effect=err)

    with pytest.raises(RuntimeError, match="Failed to reach UniProt REST API"):
        up.retrieve("P55317", dir=str(tmp_path), config=False)


def test_retrieve_config_true_writes_only_within_redirected_home(home_dir, mocker):
    fake_resp = mocker.MagicMock()
    fake_resp.read.return_value = b"content"
    fake_resp.__enter__.return_value = fake_resp
    mocker.patch("edms.data.uniprot.urllib.request.urlopen", return_value=fake_resp)

    up.retrieve("P55317", dir=None, config=True)

    expected = home_dir / ".config" / "edms" / "UniProt" / "P55317.txt"
    assert expected.exists()


# ---------------------------------------------------------------------------
# Flat file parsing: load_flat_file() / _parse_flat_text()
# ---------------------------------------------------------------------------

FLAT_FILE_TEXT = """\
ID   FOXA1_HUMAN             Reviewed;         472 AA.
AC   P55317; B2R9H6; B7ZAP5; Q9H2A0;
FT   HELIX           172..191
FT                   /evidence="ECO:0007829"
FT   STRAND          200..205
FT   MOD_RES         269
FT                   /note="Phosphoserine"
FT                   /evidence="ECO:0000269|PubMed:12345678"
FT   CARBOHYD        100
FT                   /note="N-linked (GlcNAc...) asparagine"
SQ   SEQUENCE   472 AA;  51100 MW;  ABCDEF12345 CRC64;
     MKKQ
//
"""


def test_load_flat_file_parses_accession_and_entry_type():
    entry = up.load_flat_file(FLAT_FILE_TEXT)
    assert entry.accession == "P55317"
    assert entry.entry_type == "Reviewed"
    assert len(entry.features) == 4


def test_load_flat_file_parses_feature_locations_and_qualifiers():
    entry = up.load_flat_file(FLAT_FILE_TEXT)
    by_type = {f.type: f for f in entry.features}

    helix = by_type["HELIX"]
    assert (helix.start, helix.end) == (172, 191)
    assert helix.raw["evidence_str"] == "ECO:0007829"

    strand = by_type["STRAND"]
    assert (strand.start, strand.end) == (200, 205)

    mod_res = by_type["MOD_RES"]
    assert (mod_res.start, mod_res.end) == (269, 269)
    assert mod_res.description == "Phosphoserine"
    assert mod_res.raw["evidence_str"] == "ECO:0000269|PubMed:12345678"


def test_load_flat_file_from_path(tmp_path):
    path = tmp_path / "P55317.txt"
    path.write_text(FLAT_FILE_TEXT)
    entry = up.load_flat_file(str(path))
    assert entry.accession == "P55317"


def test_load_flat_file_missing_accession_raises_value_error():
    with pytest.raises(ValueError, match="does not contain an 'AC' line"):
        up.load_flat_file("ID   FOO_HUMAN  Reviewed; 10 AA.\n")


def test_load_flat_file_rejects_non_string_non_path():
    with pytest.raises(TypeError):
        up.load_flat_file(12345)


def test_load_flat_file_unreviewed_entry_type():
    text = "ID   FOO_HUMAN  Unreviewed;  10 AA.\nAC   Q99999;\n"
    entry = up.load_flat_file(text)
    assert entry.entry_type == "Unreviewed"


# ---------------------------------------------------------------------------
# secondary_structure_from_flat_file()
# ---------------------------------------------------------------------------

def test_secondary_structure_from_flat_file_extracts_helix_and_strand():
    df = up.secondary_structure_from_flat_file(FLAT_FILE_TEXT)
    assert list(df["type"]) == ["α-helix", "β-strand"]
    assert list(df["start"]) == [172, 200]
    assert list(df["end"]) == [191, 205]
    assert list(df["name"]) == ["α-helix #1", "β-strand #1"]
    assert (df["accession"] == "P55317").all()


def test_secondary_structure_from_flat_file_numbers_multiple_helices():
    text = (
        "ID   FOO_HUMAN  Reviewed;  10 AA.\n"
        "AC   Q99999;\n"
        "FT   HELIX           1..5\n"
        "FT   HELIX           10..15\n"
    )
    df = up.secondary_structure_from_flat_file(text)
    assert list(df["name"]) == ["α-helix #1", "α-helix #2"]


# ---------------------------------------------------------------------------
# normalize_ptm_description() / ptms_from_flat_file()
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "desc,expected",
    [
        ("Phosphoserine", "Phosphorylation"),
        ("N6-acetyllysine", "Acetylation"),
        ("N6,N6-dimethyllysine", "Methylation"),
        ("Interchain (with C-256)", "interchain (with c-256)"),  # no keyword match -> lowercased passthrough
    ],
)
def test_normalize_ptm_description_keyword_matches(desc, expected):
    assert up.normalize_ptm_description(desc) == expected


def test_normalize_ptm_description_none_returns_none():
    assert up.normalize_ptm_description(None) is None


@pytest.mark.parametrize(
    "desc",
    [
        "N-linked (GlcNAc...) asparagine",
        "O-linked (GalNAc...) threonine",
        "C-linked (Man) tryptophan",
    ],
)
def test_normalize_ptm_description_classifies_real_glycosylation_text(desc):
    assert up.normalize_ptm_description(desc) == "Glycosylation"


def test_ptms_from_flat_file_extracts_and_normalizes():
    df = up.ptms_from_flat_file(FLAT_FILE_TEXT)
    assert set(df["type"]) == {"MOD_RES", "CARBOHYD"}
    mod_res_row = df[df["type"] == "MOD_RES"].iloc[0]
    assert mod_res_row["description"] == "Phosphoserine"
    assert mod_res_row["normalized_description"] == "Phosphorylation"

    carbohyd_row = df[df["type"] == "CARBOHYD"].iloc[0]
    assert carbohyd_row["description"] == "N-linked (GlcNAc...) asparagine"
    assert carbohyd_row["normalized_description"] == "Glycosylation"


def test_ptms_from_flat_file_excludes_non_ptm_features():
    df = up.ptms_from_flat_file(FLAT_FILE_TEXT)
    assert "HELIX" not in df["type"].values
    assert "STRAND" not in df["type"].values


# ---------------------------------------------------------------------------
# Feature-viewer JSON parsing: load_uniprot_feature_viewer_json()
# ---------------------------------------------------------------------------

def _feature_viewer_dict():
    return {
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "primaryAccession": "P55317",
        "features": [
            {
                "type": "Natural variant",
                "location": {
                    "start": {"value": 10, "modifier": "EXACT"},
                    "end": {"value": 10, "modifier": "EXACT"},
                },
                "description": "In cancer.",
                "featureId": "VAR_001",
                "evidences": [{"evidenceCode": "ECO:0000269", "source": "PubMed", "id": "12345"}],
                "featureCrossReferences": [{"database": "dbSNP", "id": "rs123"}],
                "alternativeSequence": {"originalSequence": "R", "alternativeSequences": ["H"]},
            }
        ],
        "extraAttributes": {"countByFeatureType": {"Natural variant": 1}},
    }


def test_load_uniprot_feature_viewer_json_from_dict():
    entry = up.load_uniprot_feature_viewer_json(_feature_viewer_dict())
    assert entry.accession == "P55317"
    assert entry.entry_type == "UniProtKB reviewed (Swiss-Prot)"
    assert len(entry.features) == 1

    feat = entry.features[0]
    assert feat.type == "Natural variant"
    assert (feat.start, feat.end) == (10, 10)
    assert feat.evidences[0].evidence_code == "ECO:0000269"
    assert feat.cross_references[0].database == "dbSNP"
    assert feat.alternative_sequence.original == "R"
    assert feat.alternative_sequence.alternatives == ["H"]


def test_load_uniprot_feature_viewer_json_from_string():
    small = {"entryType": "x", "primaryAccession": "P1", "features": []}
    entry = up.load_uniprot_feature_viewer_json(json.dumps(small))
    assert entry.accession == "P1"


def test_load_uniprot_feature_viewer_json_from_long_string():
    # A realistically-sized raw JSON string (with several features) easily
    # exceeds typical OS filename length limits (~255 chars); this must be
    # parsed as JSON text rather than attempted as a file path.
    raw_json = json.dumps(_feature_viewer_dict())
    assert len(raw_json) > 255
    entry = up.load_uniprot_feature_viewer_json(raw_json)
    assert entry.accession == "P55317"


def test_load_uniprot_feature_viewer_json_from_file(tmp_path):
    path = tmp_path / "P55317.json"
    path.write_text(json.dumps(_feature_viewer_dict()))
    entry = up.load_uniprot_feature_viewer_json(str(path))
    assert entry.accession == "P55317"


def test_load_uniprot_feature_viewer_json_missing_accession_raises():
    with pytest.raises(ValueError, match="primaryAccession"):
        up.load_uniprot_feature_viewer_json({"entryType": "x"})


def test_load_uniprot_feature_viewer_json_rejects_bad_type():
    with pytest.raises(TypeError):
        up.load_uniprot_feature_viewer_json(12345)


def test_feature_to_dataframe_flattens_evidence_and_xrefs():
    entry = up.load_uniprot_feature_viewer_json(_feature_viewer_dict())
    df = entry.to_dataframe()
    row = df.iloc[0]
    assert row["accession"] == "P55317"
    assert row["feature_type"] == "Natural variant"
    assert row["start"] == 10
    assert row["length"] == 1
    assert bool(row["has_alternative_sequence"]) is True
    assert row["original_sequence"] == "R"
    assert row["alternative_sequences"] == "H"
    assert row["evidence_codes"] == "ECO:0000269"
    assert row["xref_databases"] == "dbSNP"


def test_filter_features_matches_type():
    entry = up.load_uniprot_feature_viewer_json(_feature_viewer_dict())
    matches = entry.filter_features("Natural variant")
    assert len(matches) == 1
    assert entry.filter_features("Nonexistent") == []


def test_load_many_from_files_and_entries_to_dataframe(tmp_path):
    data1 = _feature_viewer_dict()
    data2 = _feature_viewer_dict()
    data2["primaryAccession"] = "Q99999"

    path1 = tmp_path / "a.json"
    path2 = tmp_path / "b.json"
    path1.write_text(json.dumps(data1))
    path2.write_text(json.dumps(data2))

    entries = up.load_many_from_files([path1, path2])
    assert [e.accession for e in entries] == ["P55317", "Q99999"]

    df = up.entries_to_dataframe(entries)
    assert set(df["accession"]) == {"P55317", "Q99999"}
    assert len(df) == 2


def test_entries_to_dataframe_empty_list_returns_empty_dataframe():
    df = up.entries_to_dataframe([])
    assert df.empty


# ---------------------------------------------------------------------------
# Plotting smoke tests (Agg backend set in tests/conftest.py; dir=None,
# file=None everywhere so nothing is written to disk).
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def test_draw_helix_adds_patches_to_axis():
    fig, ax = plt.subplots()
    n_before = len(ax.patches)
    up.draw_helix(ax, start=0, end=10, height=1.0)
    assert len(ax.patches) > n_before


def test_draw_helix_noop_when_end_before_start():
    fig, ax = plt.subplots()
    up.draw_helix(ax, start=10, end=5)
    assert len(ax.patches) == 0


def test_draw_strand_adds_arrow_patch():
    fig, ax = plt.subplots()
    up.draw_strand(ax, start=0, end=10)
    assert len(ax.patches) == 1


def test_draw_strand_noop_for_zero_length():
    fig, ax = plt.subplots()
    up.draw_strand(ax, start=5, end=5)
    assert len(ax.patches) == 0


def test_draw_loop_plots_sine_wave_line():
    fig, ax = plt.subplots()
    n_before = len(ax.lines)
    up.draw_loop(ax, start=0, end=10)
    assert len(ax.lines) == n_before + 1


def test_draw_ss_track_smoke_runs_without_saving(mocker):
    mocker.patch("edms.data.uniprot.p.save_fig")  # never touch disk
    df = pd.DataFrame(
        {
            "accession": ["P55317", "P55317"],
            "type": ["α-helix", "β-strand"],
            "start": [10, 30],
            "end": [20, 35],
            "name": ["α-helix #1", "β-strand #1"],
        }
    )
    up.draw_ss_track(df, dir=None, file=None, show=False)
    # No exception -> smoke test passed; a figure was created
    assert len(plt.get_fignums()) >= 1


def test_draw_ss_track_string_accession_not_found_raises_file_not_found(home_dir):
    with pytest.raises(FileNotFoundError):
        up.draw_ss_track("NOTANACCESSION", show=False)
