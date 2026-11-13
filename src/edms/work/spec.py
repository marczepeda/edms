'''
Module: spec.py
Author: Marc Zepeda
Created: 2026-11-12
Description: Run specs

A run spec is the simple input file: which template to run, which tables to
feed it, and which parameters to override. Everything else comes from the
template's defaults.

Usage:
[Run spec]
- RunSpec: template, inputs, and parameter overrides for one run
- load(): read a run spec from a .toml file
- from_args(): build a run spec from CLI arguments
- example(): write an annotated starting spec
'''
from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

from ..utils import try_parse

EXAMPLE_SPEC = '''# edms work run spec
# Fill this in, then:  edms work init -i run.toml -d runs/

template = "pe3_transfection"   # bundled template name, or a path to your own .toml

# Optional: a short label for the run directory. Defaults to the template name.
# run_id = "PRNP_tiling"

[inputs]
# Each key matches an input the template declares; each value is a file path.
plasmids  = "plasmids.csv"
epegRNAs  = "epegRNAs.csv"
ngRNAs    = "ngRNAs.csv"

[params]
# Anything you leave out keeps the template default.
# Check them with:  edms work show pe3_transfection
reps = 3

# [sheets]
# For .xlsx inputs, name the tab to read.
# plasmids = "Plasmids"
'''

@dataclass
class RunSpec:
    '''
    RunSpec: what one run of a protocol needs.

    Attributes:
    template (str): bundled template name or path to a .toml template
    run_id (str | None): short label for the run directory
    inputs (dict): input name -> file path
    params (dict): parameter overrides for the template defaults
    sheets (dict): input name -> sheet name, for .xlsx inputs
    source (Path | None): file the spec was read from; relative input paths
        resolve against its directory
    '''
    template: str
    run_id: str | None = None
    inputs: dict = field(default_factory=dict)
    params: dict = field(default_factory=dict)
    sheets: dict = field(default_factory=dict)
    source: Path | None = None

    def input_path(self, name: str) -> Path:
        '''
        input_path(): absolute path for one input, resolved against the spec file.

        Parameters:
        name (str): input name
        '''
        raw = Path(self.inputs[name]).expanduser()
        if raw.is_absolute() or self.source is None:
            return raw
        return (self.source.parent / raw).resolve()

def load(pt: str | Path) -> RunSpec:
    '''
    load(): read a run spec from a .toml file.

    Parameters:
    pt (str | Path): path to the run spec

    Dependencies: tomllib
    '''
    path = Path(pt).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f'No run spec at {path}')
    with open(path, 'rb') as fh:
        raw = tomllib.load(fh)

    if 'template' not in raw:
        raise ValueError(f'{path}: run spec needs a top-level "template"')

    return RunSpec(
        template=raw['template'],
        run_id=raw.get('run_id'),
        inputs=dict(raw.get('inputs', {})),
        params=dict(raw.get('params', {})),
        sheets=dict(raw.get('sheets', {})),
        source=path,
    )

def from_args(template: str, input: list[str] = None, set: list[str] = None,
              sheet: list[str] = None, run_id: str = None) -> RunSpec:
    '''
    from_args(): build a run spec from CLI arguments instead of a file.

    Parameters:
    template (str): bundled template name or path to a .toml template
    input (list, optional): 'name=path' pairs
    set (list, optional): 'param=value' pairs
    sheet (list, optional): 'name=sheet' pairs for .xlsx inputs
    run_id (str, optional): short label for the run directory

    Dependencies: utils.try_parse()
    '''
    return RunSpec(
        template=template,
        run_id=run_id,
        inputs=_pairs(input, parse=False),
        params=_pairs(set, parse=True),
        sheets=_pairs(sheet, parse=False),
        source=Path.cwd() / '_cli_',
    )

def merge(base: RunSpec, override: RunSpec) -> RunSpec:
    '''
    merge(): layer CLI arguments over a spec file.

    Parameters:
    base (RunSpec): spec read from a file
    override (RunSpec): spec built from CLI arguments; set fields win
    '''
    return RunSpec(
        template=override.template or base.template,
        run_id=override.run_id or base.run_id,
        inputs={**base.inputs, **override.inputs},
        params={**base.params, **override.params},
        sheets={**base.sheets, **override.sheets},
        source=base.source,
    )

def example(dir: str | Path = '.', file: str = 'run.toml') -> Path:
    '''
    example(): write an annotated starting spec and return its path.

    Parameters:
    dir (str | Path, optional): directory to write into (Default: '.')
    file (str, optional): file name (Default: 'run.toml')
    '''
    out = Path(dir).expanduser() / file
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(EXAMPLE_SPEC, encoding='utf-8')
    return out

def _pairs(items, parse: bool) -> dict:
    '''
    _pairs(): turn ['a=1','b=2'] into {'a':1,'b':2}.

    Parameters:
    items (list | None): 'key=value' strings
    parse (bool): run values through try_parse() to get real ints/floats/lists
    '''
    out = {}
    for item in items or []:
        if '=' not in item:
            raise ValueError(f"Expected 'name=value', got '{item}'")
        key, _, value = item.partition('=')
        key = key.strip()
        if not key:
            raise ValueError(f"Expected 'name=value', got '{item}'")
        out[key] = try_parse(value.strip()) if parse else value.strip()
    return out
