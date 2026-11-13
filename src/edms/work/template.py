'''
Module: template.py
Author: Marc Zepeda
Created: 2026-11-12
Description: Protocol templates

A template is declarative TOML: parameters with defaults, the input tables the
protocol needs, and an ordered list of steps. A step is prose, or it binds to a
function that already exists elsewhere in edms (e.g. edms.bio.transfect:PE3)
and lets that function do the arithmetic.

Usage:
[Paths]
- bundled_dir(): directory holding the templates shipped with edms
- list_templates(): names of the bundled protocol templates
- resolve_path(): turn a template name or file path into a path

[Loading]
- load(): read and validate a protocol template
'''
from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

# Status values a step may hold in a run ledger
STATUSES = ('pending', 'in_progress', 'done', 'skipped', 'failed')

@dataclass
class Input:
    '''
    Input: an input table a protocol template expects.

    Attributes:
    name (str): key the run spec and step args refer to
    description (str): human-readable note about the table
    required_columns (list): columns that must be present
    optional (bool): whether a run may omit the table
    '''
    name: str
    description: str = ''
    required_columns: list[str] = field(default_factory=list)
    optional: bool = False

@dataclass
class Step:
    '''
    Step: one stage of a protocol.

    Attributes:
    id (str): short unique identifier used by the CLI and depends_on
    title (str): heading shown in the rendered protocol
    day (float | None): day relative to the start of the run
    text (str): protocol prose; {param} is substituted from the run params
    depends_on (list): ids of steps that must finish first
    compute (str | None): 'module.path:function' to call for this step
    args (dict): keyword arguments for compute; '@name' is an input table and
        '$name' a parameter, anything else a literal
    outputs (list): names for the tables compute returns (optional; used for
        file naming when the function returns a bare DataFrame or a tuple)
    checks (list): things to confirm before marking the step done
    duration (str | None): rough hands-on or elapsed time
    '''
    id: str
    title: str
    day: float | None = None
    text: str = ''
    depends_on: list[str] = field(default_factory=list)
    compute: str | None = None
    args: dict = field(default_factory=dict)
    outputs: list[str] = field(default_factory=list)
    checks: list[str] = field(default_factory=list)
    duration: str | None = None

@dataclass
class Template:
    '''
    Template: a protocol template.

    Attributes:
    name (str): template name; also the bundled file stem
    version (int): bump when a change would alter an existing run's meaning
    title (str): human-readable title
    description (str): what the protocol is for
    params (dict): parameter defaults a run spec may override
    inputs (list): Input entries
    steps (list): Step entries, in order
    source (Path | None): file the template was read from
    '''
    name: str
    version: int = 1
    title: str = ''
    description: str = ''
    params: dict = field(default_factory=dict)
    inputs: list[Input] = field(default_factory=list)
    steps: list[Step] = field(default_factory=list)
    source: Path | None = None

    def step(self, id: str) -> Step:
        '''
        step(): return the step with this id.

        Parameters:
        id (str): step identifier
        '''
        for s in self.steps:
            if s.id == id:
                return s
        raise KeyError(f"Template '{self.name}' has no step '{id}'. "
                       f"Steps: {', '.join(s.id for s in self.steps)}")

# Paths
def bundled_dir() -> Path:
    '''
    bundled_dir(): directory holding the templates shipped with edms.
    '''
    return Path(__file__).resolve().parent / 'templates'

def list_templates() -> list[str]:
    '''
    list_templates(): names of the bundled protocol templates.
    '''
    d = bundled_dir()
    if not d.is_dir():
        return []
    return sorted(p.stem for p in d.glob('*.toml') if not p.stem.startswith('_'))

def resolve_path(name_or_path: str) -> Path:
    '''
    resolve_path(): turn a template name or file path into a path.

    Parameters:
    name_or_path (str): bundled template name (e.g. 'pe3_transfection') or a
        path to a .toml template file

    Dependencies: bundled_dir(), list_templates()
    '''
    p = Path(name_or_path)
    if p.suffix.lower() == '.toml' or p.exists():
        if not p.exists():
            raise FileNotFoundError(f'No template file at {p}')
        return p
    bundled = bundled_dir() / f'{name_or_path}.toml'
    if bundled.exists():
        return bundled
    available = ', '.join(list_templates()) or '(none)'
    raise FileNotFoundError(
        f"Unknown template '{name_or_path}'. Bundled templates: {available}. "
        f'Pass a path to a .toml file to use your own.')

# Loading
def load(name_or_path: str) -> Template:
    '''
    load(): read and validate a protocol template.

    Parameters:
    name_or_path (str): bundled template name or path to a .toml template

    Dependencies: tomllib, resolve_path()
    '''
    path = resolve_path(name_or_path)
    with open(path, 'rb') as fh:
        raw = tomllib.load(fh)

    if 'name' not in raw:
        raise ValueError(f'{path}: template needs a top-level "name"')

    inputs = []
    for i, entry in enumerate(raw.get('inputs', [])):
        if 'name' not in entry:
            raise ValueError(f'{path}: inputs[{i}] needs a "name"')
        inputs.append(Input(
            name=entry['name'],
            description=entry.get('description', ''),
            required_columns=list(entry.get('required_columns', [])),
            optional=bool(entry.get('optional', False)),
        ))

    steps = []
    for i, entry in enumerate(raw.get('steps', [])):
        if 'id' not in entry:
            raise ValueError(f'{path}: steps[{i}] needs an "id"')
        steps.append(Step(
            id=entry['id'],
            title=entry.get('title', entry['id']),
            day=entry.get('day'),
            text=entry.get('text', ''),
            depends_on=list(entry.get('depends_on', [])),
            compute=entry.get('compute'),
            args=dict(entry.get('args', {})),
            outputs=list(entry.get('outputs', [])),
            checks=list(entry.get('checks', [])),
            duration=entry.get('duration'),
        ))

    tpl = Template(
        name=raw['name'],
        version=int(raw.get('version', 1)),
        title=raw.get('title', raw['name']),
        description=raw.get('description', ''),
        params=dict(raw.get('params', {})),
        inputs=inputs,
        steps=steps,
        source=path,
    )
    validate(tpl)
    return tpl

def validate(tpl: Template) -> Template:
    '''
    validate(): check a template for duplicate ids, bad dependencies, bad
    compute targets, and step args that reference things the template does not
    declare. Raises ValueError on the first problem found.

    Parameters:
    tpl (Template): template to check
    '''
    where = str(tpl.source) if tpl.source else tpl.name

    if not tpl.steps:
        raise ValueError(f'{where}: template has no steps')

    seen = set()
    for s in tpl.steps:
        if s.id in seen:
            raise ValueError(f"{where}: duplicate step id '{s.id}'")
        seen.add(s.id)

    input_names = {i.name for i in tpl.inputs}
    done = set()
    for s in tpl.steps:
        for dep in s.depends_on:
            if dep not in seen:
                raise ValueError(f"{where}: step '{s.id}' depends on unknown step '{dep}'")
            if dep not in done:
                raise ValueError(f"{where}: step '{s.id}' depends on '{dep}', which comes later")
        if s.compute is not None and ':' not in s.compute:
            raise ValueError(f"{where}: step '{s.id}' compute must look like "
                             f"'module.path:function', got '{s.compute}'")
        for key, value in s.args.items():
            for ref, kind, known in ((value, '@', input_names), (value, '$', set(tpl.params))):
                if isinstance(ref, str) and ref.startswith(kind) and not ref.startswith(kind * 2):
                    target = ref[1:]
                    if target not in known:
                        label = 'input' if kind == '@' else 'param'
                        raise ValueError(
                            f"{where}: step '{s.id}' arg '{key}' references "
                            f"{label} '{target}', which the template does not declare")
        done.add(s.id)

    return tpl
