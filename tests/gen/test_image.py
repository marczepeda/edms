"""
Tests for edms.gen.image
"""
from PIL import Image

from edms.gen import image as img


def _make_png(pt, size=(20, 10), color=(255, 0, 0)):
    Image.new("RGB", size, color).save(pt)


# --------------------------------------------------------------------------- #
# crop()
# --------------------------------------------------------------------------- #

def test_crop_with_box_dims(tmp_path):
    in_dir = tmp_path / "in"
    out_dir = tmp_path / "out"
    in_dir.mkdir()
    _make_png(in_dir / "a.png", size=(20, 10))

    img.crop(str(in_dir), str(out_dir), box_dims=(0, 0, 10, 5))

    out_img = Image.open(out_dir / "a.png")
    assert out_img.size == (10, 5)


def test_crop_with_box_fracs(tmp_path):
    in_dir = tmp_path / "in"
    out_dir = tmp_path / "out"
    in_dir.mkdir()
    _make_png(in_dir / "a.png", size=(100, 50))

    img.crop(str(in_dir), str(out_dir), box_fracs=(0.0, 0.0, 0.5, 0.5))

    out_img = Image.open(out_dir / "a.png")
    assert out_img.size == (50, 25)


def test_crop_ignores_non_image_files(tmp_path):
    in_dir = tmp_path / "in"
    out_dir = tmp_path / "out"
    in_dir.mkdir()
    _make_png(in_dir / "a.png", size=(20, 10))
    (in_dir / "notes.txt").write_text("hello")

    img.crop(str(in_dir), str(out_dir), box_dims=(0, 0, 10, 5))

    assert (out_dir / "a.png").exists()
    assert not (out_dir / "notes.txt").exists()


# --------------------------------------------------------------------------- #
# convert()
# --------------------------------------------------------------------------- #

def test_convert_png_to_jpeg(tmp_path):
    in_dir = tmp_path / "in"
    out_dir = tmp_path / "out"
    in_dir.mkdir()
    _make_png(in_dir / "a.png")

    img.convert(str(in_dir), str(out_dir), suffix=".jpeg")

    out_file = out_dir / "a.jpeg"
    assert out_file.exists()
    out_img = Image.open(out_file)
    assert out_img.format == "JPEG"


# --------------------------------------------------------------------------- #
# combine()
# --------------------------------------------------------------------------- #

def test_combine_images_into_pdf(tmp_path):
    in_dir = tmp_path / "in"
    out_dir = tmp_path / "out"
    in_dir.mkdir()
    _make_png(in_dir / "a.png")
    _make_png(in_dir / "b.png", color=(0, 255, 0))

    img.combine(str(in_dir), str(out_dir), out_file="combined.pdf")

    out_pdf = out_dir / "combined.pdf"
    assert out_pdf.exists()
    assert out_pdf.stat().st_size > 0


def test_combine_no_images_prints_message(tmp_path, capsys):
    in_dir = tmp_path / "in"
    out_dir = tmp_path / "out"
    in_dir.mkdir()
    (in_dir / "notes.txt").write_text("hello")

    img.combine(str(in_dir), str(out_dir), out_file="combined.pdf")

    captured = capsys.readouterr()
    assert "No valid images found" in captured.out
    assert not (out_dir / "combined.pdf").exists()


# --------------------------------------------------------------------------- #
# info()
# --------------------------------------------------------------------------- #

def test_info_returns_dataframe_with_expected_columns(tmp_path):
    _make_png(tmp_path / "a.png", size=(30, 15))

    df = img.info(str(tmp_path))

    assert list(df.columns) == ["Filename", "Width", "Height", "Mode", "Format", "Size (bytes)"]
    row = df[df["Filename"] == "a.png"].iloc[0]
    assert row["Width"] == 30
    assert row["Height"] == 15
    assert row["Format"] == "PNG"


def test_info_skips_non_image_files(tmp_path, capsys):
    (tmp_path / "notes.txt").write_text("hello")
    df = img.info(str(tmp_path))
    assert df.empty
