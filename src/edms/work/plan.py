'''
Module: plan.py
Author: Marc Zepeda
Created: 2026-11-12
Description: Resolve a template and a run into a concrete protocol

Substitutes the run's parameters into each step's prose and, where a step
binds to a function, calls that function so the protocol arrives with the
arithmetic already done. The functions live in edms.bio; this module only
wires them up.

Usage:
[Resolution]
- ResolvedStep: one step with its text filled in and its tables computed
- ResolvedProtocol: the whole protocol for one run
- resolve(): template + run ledger -> ResolvedProtocol
- load_inputs(): read a run's frozen input tables

[Helpers]
- fmt(): substitute {param} without failing on unknown keys
- call(): import and call a step's compute target
'''
from __future__ import annotations

import importlib
import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from ..gen import io
from . import run as run_mod
from . import template as template_mod

@dataclass
class ResolvedStep:
    '''
    ResolvedStep: one step of a protocol, ready to follow.

    Attributes:
    id (str): step identifier
    title (str): heading
    day (float | None): day relative to the start of the run
    text (str): prose with {param} substituted
    duration (str | None): rough time
    checks (list): things to confirm before marking it done
    depends_on (list): ids of steps that must finish first
    status (str): status from the run ledger
    notes (list): notes recorded against the step
    deviations (list): parameter changes recorded against the step
    compute (str | None): the function this step called
    tables (dict): name -> DataFrame produced by compute
    error (str | None): why compute failed, if it did
    compute_skipped (str | None): why compute was not attempted, when an
        optional input the step needs was not supplied
    '''
    id: str
    title: str
    day: float | None = None
    text: str = ''
    duration: str | None = None
    checks: list = field(default_factory=list)
    depends_on: list = field(default_factory=list)
    status: str = 'pending'
    notes: list = field(default_factory=list)
    deviations: list = field(default_factory=list)
    compute: str | None = None
    tables: dict = field(default_factory=dict)
    error: str | None = None
    compute_skipped: str | None = None

@dataclass
class ResolvedProtocol:
    '''
    ResolvedProtocol: a protocol resolved against one run.

    Attributes:
    run_id (str): run identifier
    title (str): template title
    description (str): template description
    template (dict): name, version, title, source
    params (dict): resolved parameters, including any recorded deviations
    inputs (dict): input name -> DataFrame
    steps (list): ResolvedStep entries
    run_dir (Path): the run directory
    errors (list): (step id, message) for every compute that failed
    '''
    run_id: str
    title: str
    description: str
    template: dict
    params: dict
    inputs: dict
    steps: list
    run_dir: Path
    errors: list = field(default_factory=list)

# Helpers
class _SafeDict(dict):
    '''
    _SafeDict: leaves {unknown} alone instead of raising KeyError.
    '''
    def __missing__(self, key):
        return '{' + key + '}'

def fmt(text: str, params: dict) -> str:
    '''
    fmt(): substitute {param} into protocol prose without failing on unknown
    keys or on stray braces.

    Parameters:
    text (str): prose from a template step
    params (dict): resolved run parameters
    '''
    if not text:
        return ''
    try:
        return text.format_map(_SafeDict(params))
    except (ValueError, IndexError):
        return text

def _resolve_arg(value, inputs: dict, params: dict, step_id: str, key: str):
    '''
    _resolve_arg(): turn one template arg into a real Python value.

    '@name' is an input table, '$name' a parameter, '@@'/'$$' an escaped
    literal, and anything else passes through untouched.

    Parameters:
    value: the raw arg from the template
    inputs (dict): input name -> DataFrame
    params (dict): resolved run parameters
    step_id (str): step the arg belongs to, for error messages
    key (str): arg name, for error messages
    '''
    if isinstance(value, str):
        for sigil, source, label in (('@', inputs, 'input'), ('$', params, 'parameter')):
            if value.startswith(sigil * 2):
                return value[1:]
            if value.startswith(sigil):
                name = value[1:]
                if name not in source:
                    raise KeyError(
                        f"Step '{step_id}' arg '{key}' wants {label} '{name}', "
                        f"which this run does not have. Available: "
                        f"{', '.join(sorted(source)) or '(none)'}")
                return source[name]
        return value
    if isinstance(value, list):
        return [_resolve_arg(v, inputs, params, step_id, key) for v in value]
    if isinstance(value, dict):
        return {k: _resolve_arg(v, inputs, params, step_id, key) for k, v in value.items()}
    return value

def input_refs(args) -> set:
    '''
    input_refs(): the input names a step's args reference with '@'.

    Parameters:
    args: a step's args dict, or any value nested inside it
    '''
    found = set()
    if isinstance(args, str):
        if args.startswith('@') and not args.startswith('@@'):
            found.add(args[1:])
    elif isinstance(args, dict):
        for v in args.values():
            found |= input_refs(v)
    elif isinstance(args, (list, tuple)):
        for v in args:
            found |= input_refs(v)
    return found

def call(compute: str, args: dict):
    '''
    call(): import and call a step's compute target.

    Parameters:
    compute (str): 'module.path:function'
    args (dict): already-resolved keyword arguments

    Dependencies: importlib
    '''
    module_path, _, func_name = compute.partition(':')
    try:
        module = importlib.import_module(module_path)
    except ImportError as e:
        raise ImportError(f"Cannot import '{module_path}' for compute '{compute}': {e}") from e
    func = getattr(module, func_name, None)
    if func is None:
        raise AttributeError(f"'{module_path}' has no function '{func_name}'")
    return func(**args)

def _tables(result, names: list) -> dict:
    '''
    _tables(): normalise whatever a compute function returned into a dict of
    DataFrames. Handles a bare DataFrame, a dict of them, and the tuples of
    dicts that some edms planners return.

    Parameters:
    result: the compute function's return value
    names (list): preferred output names from the template
    '''
    out = {}

    def add(key, value):
        key = str(key)
        if key in out:
            i = 2
            while f'{key}_{i}' in out:
                i += 1
            key = f'{key}_{i}'
        out[key] = value

    def walk(obj, hint):
        if isinstance(obj, pd.DataFrame):
            add(hint, obj)
        elif isinstance(obj, dict):
            for k, v in obj.items():
                walk(v, k)
        elif isinstance(obj, (list, tuple)):
            for i, v in enumerate(obj):
                walk(v, names[i] if i < len(names) else f'{hint}_{i + 1}')

    walk(result, names[0] if names else 'result')
    return out

def _sheet_name(name: str, used: set) -> str:
    '''
    _sheet_name(): make an Excel-safe, unique sheet name.

    Parameters:
    name (str): desired name
    used (set): names already taken in this workbook
    '''
    safe = re.sub(r'[\\/*?:\[\]]', '_', str(name))[:31] or 'Sheet'
    base = safe
    i = 2
    while safe in used:
        suffix = f'_{i}'
        safe = base[:31 - len(suffix)] + suffix
        i += 1
    used.add(safe)
    return safe

def load_inputs(state: dict, run_dir: str | Path) -> dict:
    '''
    load_inputs(): read a run's frozen input tables.

    Parameters:
    state (dict): run ledger
    run_dir (str | Path): run directory

    Dependencies: gen.io.get()
    '''
    run_dir = Path(run_dir).expanduser()
    out = {}
    for name, meta in state.get('inputs', {}).items():
        path = run_dir / meta['frozen']
        if not path.exists():
            raise FileNotFoundError(f"Frozen input '{name}' is missing from {path}")
        obj = io.get(pt=str(path))
        if isinstance(obj, dict):
            sheet = meta.get('sheet')
            if sheet is None:
                if len(obj) == 1:
                    obj = next(iter(obj.values()))
                else:
                    raise ValueError(
                        f"Input '{name}' has several tabs ({', '.join(obj)}). "
                        f'Name one with a [sheets] entry in the run spec.')
            elif sheet not in obj:
                raise ValueError(
                    f"Input '{name}' has no tab '{sheet}'. Tabs: {', '.join(obj)}")
            else:
                obj = obj[sheet]
        out[name] = obj
    return out

# Resolution
def resolve(run_dir: str | Path, tpl=None, strict: bool = False,
            compute: bool = True, write: bool = True) -> ResolvedProtocol:
    '''
    resolve(): turn a run ledger plus its template into a concrete protocol.

    A compute step that raises is recorded on the step and in protocol.errors
    rather than aborting the whole render, so a protocol still comes out when
    one calculation is not ready yet. Pass strict=True to raise instead.

    Parameters:
    run_dir (str | Path): run directory
    tpl (Template, optional): loaded template (Default: None -> from the ledger)
    strict (bool, optional): raise on the first compute failure (Default: False)
    compute (bool, optional): run compute steps at all (Default: True)
    write (bool, optional): write computed tables to run_dir/outputs (Default: True)

    Dependencies: pandas, run.load(), template.load(), load_inputs(), fmt(), call()
    '''
    run_dir = Path(run_dir).expanduser()
    state = run_mod.load(run_dir)

    if tpl is None:
        ref = state['template']
        source = ref.get('source')
        tpl = template_mod.load(source if source and Path(source).exists() else ref['name'])

    params = dict(state['params'])
    inputs = load_inputs(state, run_dir) if compute else {}
    by_id = {e['id']: e for e in state['steps']}
    optional_absent = {i.name for i in tpl.inputs if i.optional} - set(inputs)

    steps = []
    errors = []
    out_dir = run_dir / 'outputs'

    for s in tpl.steps:
        entry = by_id.get(s.id, {})
        resolved = ResolvedStep(
            id=s.id,
            title=s.title,
            day=s.day,
            text=fmt(s.text, params),
            duration=s.duration,
            checks=[fmt(c, params) for c in s.checks],
            depends_on=list(s.depends_on),
            status=entry.get('status', 'pending'),
            notes=list(entry.get('notes', [])),
            deviations=list(entry.get('deviations', [])),
            compute=s.compute,
        )

        if s.compute and compute:
            absent = sorted(input_refs(s.args) & optional_absent)
            if absent:
                resolved.compute_skipped = (
                    f"needs optional input(s) {', '.join(absent)}, which this run "
                    f'does not supply')
                steps.append(resolved)
                continue
            try:
                args = {k: _resolve_arg(v, inputs, params, s.id, k) for k, v in s.args.items()}
                resolved.tables = _tables(call(s.compute, args), s.outputs)
            except Exception as e:
                message = f'{type(e).__name__}: {e}'
                resolved.error = message
                errors.append((s.id, message))
                if strict:
                    raise
            else:
                if write and resolved.tables:
                    out_dir.mkdir(parents=True, exist_ok=True)
                    path = out_dir / f'{s.id}.xlsx'
                    used = set()
                    with pd.ExcelWriter(path, engine='openpyxl') as writer:
                        for name, table in resolved.tables.items():
                            table.to_excel(writer, sheet_name=_sheet_name(name, used),
                                           index=False)
                    run_mod.record_outputs(run_dir, s.id,
                                           [str(path.relative_to(run_dir))])

        steps.append(resolved)

    return ResolvedProtocol(
        run_id=state['run_id'],
        title=tpl.title or tpl.name,
        description=tpl.description,
        template=state['template'],
        params=params,
        inputs=inputs,
        steps=steps,
        run_dir=run_dir,
        errors=errors,
    )
