"""Shared fixtures for tests/data (edms.data.* external-database wrappers).

These modules default to writing "cache" files under ~/.config/edms/... .
To guarantee that no test ever touches the real home directory, every test
that exercises a code path with such a default must use the `home_dir`
fixture below, which repoints HOME (and USERPROFILE, for completeness) at a
pytest tmp_path for the duration of the test.
"""
import os

import pytest


@pytest.fixture
def home_dir(tmp_path, monkeypatch):
    """Redirect the user's home directory to an isolated tmp_path.

    Any code that calls os.path.expanduser("~/...") during the test will
    resolve inside tmp_path instead of the real user home directory.
    """
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))  # harmless on POSIX, safe on Windows
    return tmp_path


def make_fake_response(text: str = "", content: bytes | None = None):
    """Build a MagicMock that behaves like the object returned by
    urllib.request.urlopen(...) when used as a context manager.
    """
    from unittest.mock import MagicMock

    if content is None:
        content = text.encode("utf-8")

    resp = MagicMock()
    resp.read.return_value = content
    resp.__enter__.return_value = resp
    resp.__exit__.return_value = False
    return resp
