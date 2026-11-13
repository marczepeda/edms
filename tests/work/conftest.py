"""
Shared fixtures for edms.work tests.

Everything is built inside tmp_path: templates, input tables, and run
directories. Nothing reads or writes the user's real run directories.
"""
import textwrap

import pandas as pd
import pytest

from edms.work import run as run_mod
from edms.work import spec as spec_mod
from edms.work import template as template_mod


MINIMAL = """
name = "minimal"
version = 1
title = "Minimal protocol"
description = "Two steps, one of which computes."

[params]
scale = 2
label = "widget"

[[inputs]]
name = "samples"
required_columns = ["id"]

[[inputs]]
name = "extra"
optional = true

[[steps]]
id = "prep"
title = "Prepare"
day = 0
text = "Prepare {scale} x {label}."
checks = ["Bench is clear"]

[[steps]]
id = "build"
title = "Build a table"
day = 1
depends_on = ["prep"]
compute = "pandas:DataFrame"
outputs = ["built"]

[steps.args.data]
n = [1, 2, 3]

[[steps]]
id = "optional_step"
title = "Needs the optional input"
depends_on = ["build"]
compute = "pandas:DataFrame"

[steps.args]
data = "@extra"
"""


@pytest.fixture
def template_file(tmp_path):
    """A minimal template on disk."""
    p = tmp_path / "minimal.toml"
    p.write_text(textwrap.dedent(MINIMAL), encoding="utf-8")
    return p


@pytest.fixture
def tpl(template_file):
    """The minimal template, loaded."""
    return template_mod.load(str(template_file))


@pytest.fixture
def samples_csv(tmp_path):
    """A small input table."""
    p = tmp_path / "samples.csv"
    pd.DataFrame({"id": ["s1", "s2"], "value": [10, 20]}).to_csv(p, index=False)
    return p


@pytest.fixture
def run_spec(template_file, samples_csv, tmp_path):
    """A run spec pointing at the minimal template."""
    return spec_mod.RunSpec(
        template=str(template_file),
        run_id="unit",
        inputs={"samples": str(samples_csv)},
        params={"scale": 5},
        source=tmp_path / "run.toml",
    )


@pytest.fixture
def started_run(run_spec, tpl, tmp_path):
    """An initialised run directory; returns (run_dir, state)."""
    return run_mod.init(run_spec, tpl=tpl, dir=tmp_path / "runs")
