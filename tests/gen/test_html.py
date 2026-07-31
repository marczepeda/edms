"""
Tests for edms.gen.html
"""
import time

import pytest

from edms.gen import html as ht


def _write_html(pt, title, body="body"):
    pt.write_text(f"<html><head><title>{title}</title></head><body>{body}</body></html>")


def test_make_html_index_basic(tmp_path):
    _write_html(tmp_path / "a.html", "Alpha Page")
    _write_html(tmp_path / "b.html", "Beta Page")

    out_path = ht.make_html_index(dir=str(tmp_path), file="index.html", preview=False)
    assert out_path.exists()
    content = out_path.read_text()
    assert "Alpha Page" in content
    assert "Beta Page" in content
    # icon copied alongside
    assert (tmp_path / "index" / "python.svg").exists()


def test_make_html_index_excludes_self_and_named_excludes(tmp_path):
    _write_html(tmp_path / "index.html", "Should Be Excluded")
    _write_html(tmp_path / "keep.html", "Keep Me")
    out_path = ht.make_html_index(dir=str(tmp_path), file="index.html", preview=False)
    content = out_path.read_text()
    assert "Keep Me" in content


def test_make_html_index_requires_html_extension(tmp_path):
    _write_html(tmp_path / "a.html", "A")
    with pytest.raises(ValueError):
        ht.make_html_index(dir=str(tmp_path), file="index.txt")


def test_make_html_index_no_files_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        ht.make_html_index(dir=str(tmp_path), file="index.html")


def test_make_html_index_falls_back_to_stem_when_no_title(tmp_path):
    (tmp_path / "notitle.html").write_text("<html><body>no title here</body></html>")
    out_path = ht.make_html_index(dir=str(tmp_path), file="index.html", preview=False)
    content = out_path.read_text()
    assert "notitle" in content


def test_make_html_index_label_name(tmp_path):
    _write_html(tmp_path / "a.html", "Alpha Page")
    out_path = ht.make_html_index(dir=str(tmp_path), file="index.html", label="name", preview=False)
    content = out_path.read_text()
    assert "a.html" in content


def test_make_html_index_sort_name(tmp_path):
    _write_html(tmp_path / "b.html", "Zeta")
    _write_html(tmp_path / "a.html", "Alpha")
    out_path = ht.make_html_index(dir=str(tmp_path), file="index.html", sort="name", label="name", preview=False)
    content = out_path.read_text()
    # a.html should appear before b.html in the markup when sorted by name
    assert content.index("a.html") < content.index("b.html")


def test_make_html_index_sort_mtime(tmp_path):
    _write_html(tmp_path / "old.html", "Old")
    time.sleep(0.02)
    _write_html(tmp_path / "new.html", "New")
    out_path = ht.make_html_index(dir=str(tmp_path), file="index.html", sort="mtime", label="name", preview=False)
    content = out_path.read_text()
    # mtime sort is reverse=True -> newest first
    assert content.index("new.html") < content.index("old.html")


def test_make_html_index_recursive_grouping(tmp_path):
    sub = tmp_path / "subdir"
    sub.mkdir()
    _write_html(tmp_path / "top.html", "Top Page")
    _write_html(sub / "nested.html", "Nested Page")
    out_path = ht.make_html_index(dir=str(tmp_path), file="index.html", recursive=True, preview=False)
    content = out_path.read_text()
    assert "Top Page" in content
    assert "Nested Page" in content
    assert "root" in content  # group label for top-level files
    assert "subdir" in content  # group label for nested files


def test_make_html_index_image_types_included(tmp_path):
    _write_html(tmp_path / "a.html", "A Page")
    (tmp_path / "pic.png").write_bytes(b"\x89PNG\r\n\x1a\n" + b"0" * 20)
    out_path = ht.make_html_index(
        dir=str(tmp_path), file="index.html", image_types=[".png"], label="name", preview=False
    )
    content = out_path.read_text()
    assert "pic.png" in content


def test_make_html_index_with_preview_iframe(tmp_path):
    _write_html(tmp_path / "a.html", "Alpha")
    out_path = ht.make_html_index(dir=str(tmp_path), file="index.html", preview=True)
    content = out_path.read_text()
    assert '<iframe id="viewer"' in content


def test_make_html_index_grid_cols_applied(tmp_path):
    _write_html(tmp_path / "a.html", "Alpha")
    out_path = ht.make_html_index(dir=str(tmp_path), file="index.html", grid_cols=7, preview=False)
    content = out_path.read_text()
    assert "repeat(7, minmax(240px, 1fr))" in content


def test_make_html_index_output_file_relative_to_dir(tmp_path):
    _write_html(tmp_path / "a.html", "Alpha")
    out_path = ht.make_html_index(dir=str(tmp_path), file="custom_index.html", preview=False)
    assert out_path == tmp_path / "custom_index.html"
    assert out_path.exists()
