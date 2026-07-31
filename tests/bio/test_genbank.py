"""
Tests for edms.bio.genbank

Covers:
- viewer(): parses a synthetic GenBank record (built in-test with Bio.SeqRecord /
  Bio.SeqFeature, written to a temp file via Bio.SeqIO), and verifies:
  - features are extracted with correct start/end coordinates
  - label priority order (label > gene > product > note > 'Unnamed')
  - default vs. user-overridden feature colors
  - `exclude` filters out feature types
  - `region` filters to only features overlapping the requested window
No plot windows are opened (matplotlib Agg backend is set in tests/conftest.py) and
`show=False` is always passed so pyplot.show() is a no-op under Agg.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from edms.bio import genbank


def _write_genbank(pt, features, seq_len=60, record_id="test_construct"):
    seq = Seq("ATGC" * (seq_len // 4 + 1))[:seq_len]
    rec = SeqRecord(seq, id=record_id, name=record_id, description="synthetic test construct")
    rec.annotations["molecule_type"] = "DNA"
    rec.features = features
    SeqIO.write(rec, str(pt), "genbank")
    return pt


@pytest.fixture
def gb_path(tmp_path):
    features = [
        SeqFeature(FeatureLocation(0, 60), type="source", qualifiers={"organism": ["synthetic"]}),
        SeqFeature(FeatureLocation(0, 30), type="gene", qualifiers={"gene": ["testGene"]}),
        SeqFeature(FeatureLocation(0, 30), type="CDS", qualifiers={"gene": ["testGene"], "product": ["test protein"]}),
        SeqFeature(FeatureLocation(40, 55), type="misc_feature", qualifiers={"note": ["a note only"]}),
    ]
    pt = tmp_path / "construct.gb"
    _write_genbank(pt, features)
    yield str(pt)
    plt.close("all")


def test_viewer_returns_features_with_correct_coordinates(gb_path):
    features = genbank.viewer(pt=gb_path, show=False, return_features=True)
    by_type = {typ: gf for typ, gf in features}

    assert set(by_type.keys()) == {"source", "gene", "CDS", "misc_feature"}
    assert by_type["source"].start == 0 and by_type["source"].end == 60
    assert by_type["gene"].start == 0 and by_type["gene"].end == 30
    assert by_type["CDS"].start == 0 and by_type["CDS"].end == 30
    assert by_type["misc_feature"].start == 40 and by_type["misc_feature"].end == 55


def test_viewer_label_priority_gene_over_product(gb_path):
    features = genbank.viewer(pt=gb_path, show=False, return_features=True)
    by_type = {typ: gf for typ, gf in features}

    # CDS has both 'gene' and 'product' qualifiers; 'gene' should win per priority order.
    assert by_type["CDS"].label == "testGene"
    assert by_type["gene"].label == "testGene"


def test_viewer_label_falls_back_to_note(gb_path):
    features = genbank.viewer(pt=gb_path, show=False, return_features=True)
    by_type = {typ: gf for typ, gf in features}
    assert by_type["misc_feature"].label == "a note only"


def test_viewer_label_falls_back_to_unnamed(gb_path):
    features = genbank.viewer(pt=gb_path, show=False, return_features=True)
    by_type = {typ: gf for typ, gf in features}
    # 'source' only has an 'organism' qualifier, which is not in the label priority list.
    assert by_type["source"].label == "Unnamed"


def test_viewer_default_color_scheme(gb_path):
    features = genbank.viewer(pt=gb_path, show=False, return_features=True)
    by_type = {typ: gf for typ, gf in features}
    assert by_type["gene"].color == "#00B200"
    assert by_type["CDS"].color == "#FFFF00"
    assert by_type["source"].color == "#0000FF"
    # Type not in the built-in color_scheme dict falls back to gray.
    assert by_type["misc_feature"].color == "#AAAAAA"


def test_viewer_feature_colors_override(gb_path):
    features = genbank.viewer(
        pt=gb_path, show=False, return_features=True,
        feature_colors={"gene": "#123456"},
    )
    by_type = {typ: gf for typ, gf in features}
    assert by_type["gene"].color == "#123456"
    # Non-overridden types keep their default color.
    assert by_type["CDS"].color == "#FFFF00"


def test_viewer_exclude_removes_feature_type(gb_path):
    features = genbank.viewer(pt=gb_path, show=False, return_features=True, exclude=["source", "misc_feature"])
    types = {typ for typ, _ in features}
    assert types == {"gene", "CDS"}


def test_viewer_region_filters_to_overlapping_features(gb_path):
    # Region [35, 60] should only catch the misc_feature (40-55) and source (0-60,
    # since its end (60) falls within the region), not gene/CDS (0-30).
    features = genbank.viewer(pt=gb_path, show=False, return_features=True, region=(35, 60))
    types = {typ for typ, _ in features}
    assert "misc_feature" in types
    assert "gene" not in types
    assert "CDS" not in types


def test_viewer_return_features_false_returns_none(gb_path):
    result = genbank.viewer(pt=gb_path, show=False, return_features=False)
    assert result is None


def test_viewer_saves_plot_file(gb_path, tmp_path):
    out_dir = tmp_path / "plots"
    genbank.viewer(pt=gb_path, show=False, dir=str(out_dir), file="construct.png", return_features=False)
    assert os.path.isfile(out_dir / "construct.png")
