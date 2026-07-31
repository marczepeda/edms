import os
import re
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

REPO_ROOT = Path(__file__).resolve().parent.parent
TESTS_OUT = REPO_ROOT / "tests_out"


@pytest.fixture(autouse=True)
def _run_in_tests_out(request):
    """Safety net: run every test with cwd inside tests_out/, so any function that
    falls back to writing relative to the current directory (e.g. a missing `dir=`
    argument) lands in a gitignored scratch folder instead of the repo root.
    Tests that use tmp_path directly are unaffected; this only catches the rest.
    """
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", request.node.name)
    out_dir = TESTS_OUT / request.node.module.__name__.replace(".", "/") / safe_name
    out_dir.mkdir(parents=True, exist_ok=True)
    prev_cwd = os.getcwd()
    os.chdir(out_dir)
    try:
        yield out_dir
    finally:
        os.chdir(prev_cwd)
