"""
Tests for edms.gen.com

com.py touches the shell/OS directly (subprocess, and reading/writing shell
config files), so extra care is taken to keep everything inside tmp_path and
to monkeypatch com.shell_configs so tests never touch the user's real
~/.bashrc or ~/.zshrc.
"""
import gzip
import os

import pytest

from edms.gen import com


# --------------------------------------------------------------------------- #
# access()
# --------------------------------------------------------------------------- #

def test_access_runs_chmod_without_error(tmp_path, capsys):
    (tmp_path / "file.txt").write_text("hello")
    # Should complete without raising, regardless of whether the underlying
    # chmod produces platform-specific warnings on stderr.
    com.access(str(tmp_path))
    captured = capsys.readouterr()
    assert "terminal:" in captured.out


# --------------------------------------------------------------------------- #
# smaller_fastq()
# --------------------------------------------------------------------------- #

def test_smaller_fastq_gzipped_limits_reads(tmp_path):
    lines = "\n".join(f"line{i}" for i in range(1, 21)) + "\n"  # 20 lines = 5 reads
    gz_path = tmp_path / "sample.fastq.gz"
    with gzip.open(gz_path, "wt") as f:
        f.write(lines)

    com.smaller_fastq(str(tmp_path), reads=2, suf=".fastq.gz")

    out_dir = tmp_path / "2_reads"
    out_file = out_dir / "sample.fastq"
    assert out_file.exists()
    out_lines = out_file.read_text().splitlines()
    assert len(out_lines) == 8  # 4 lines/read * 2 reads
    assert out_lines[0] == "line1"


def test_smaller_fastq_unzipped_limits_reads(tmp_path):
    lines = "\n".join(f"line{i}" for i in range(1, 13)) + "\n"  # 12 lines = 3 reads
    fq_path = tmp_path / "sample.fastq"
    fq_path.write_text(lines)

    com.smaller_fastq(str(tmp_path), reads=1, suf=".fastq")

    out_file = tmp_path / "1_reads" / "sample.fastq"
    assert out_file.exists()
    out_lines = out_file.read_text().splitlines()
    assert len(out_lines) == 4
    assert out_lines == ["line1", "line2", "line3", "line4"]


# --------------------------------------------------------------------------- #
# detect_shell()
# --------------------------------------------------------------------------- #

def test_detect_shell_zsh(monkeypatch):
    monkeypatch.setenv("SHELL", "/bin/zsh")
    assert com.detect_shell() == "zsh"


def test_detect_shell_bash(monkeypatch):
    monkeypatch.setenv("SHELL", "/bin/bash")
    assert com.detect_shell() == "bash"


def test_detect_shell_unsupported_raises(monkeypatch):
    monkeypatch.setenv("SHELL", "/bin/fish")
    with pytest.raises(ValueError):
        com.detect_shell()


# --------------------------------------------------------------------------- #
# create_export_var() / view_export_vars()
#
# shell_configs is monkeypatched to point at files inside tmp_path so these
# tests can never touch the user's real ~/.bashrc or ~/.zshrc.
# --------------------------------------------------------------------------- #

@pytest.fixture
def fake_shell_config(tmp_path, monkeypatch):
    fake_rc = tmp_path / "fake_zshrc"
    monkeypatch.setitem(com.shell_configs, "zsh", fake_rc)
    return fake_rc


def test_create_export_var_appends_line(fake_shell_config):
    com.create_export_var("MYVAR", "/some/path", shell="zsh")
    content = fake_shell_config.read_text()
    assert 'export MYVAR="/some/path"' in content


def test_create_export_var_uppercases_name(fake_shell_config):
    com.create_export_var("myvar", "/some/path", shell="zsh")
    content = fake_shell_config.read_text()
    assert 'export MYVAR=' in content


def test_create_export_var_no_duplicate(fake_shell_config, capsys):
    com.create_export_var("MYVAR", "/some/path", shell="zsh")
    com.create_export_var("MYVAR", "/some/path", shell="zsh")
    content = fake_shell_config.read_text()
    assert content.count('export MYVAR="/some/path"') == 1
    captured = capsys.readouterr()
    assert "already exists" in captured.out


def test_create_export_var_unsupported_shell_raises(tmp_path, monkeypatch):
    with pytest.raises(ValueError):
        com.create_export_var("MYVAR", "/some/path", shell="fish")


def test_view_export_vars_prints_exports(fake_shell_config, capsys):
    fake_shell_config.write_text("export FOO=\"bar\"\nsome other line\nexport BAZ=\"qux\"\n")
    com.view_export_vars(shell="zsh")
    captured = capsys.readouterr()
    assert 'export FOO="bar"' in captured.out
    assert 'export BAZ="qux"' in captured.out
    assert "some other line" not in captured.out


def test_view_export_vars_missing_file_reports_not_found(tmp_path, monkeypatch, capsys):
    missing = tmp_path / "does_not_exist_rc"
    monkeypatch.setitem(com.shell_configs, "zsh", missing)
    com.view_export_vars(shell="zsh")
    captured = capsys.readouterr()
    assert "not found" in captured.out
