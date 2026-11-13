'''
Module: run.py
Author: Marc Zepeda
Created: 2026-11-12
Description: Run state

A run is a directory holding run.json: the frozen template reference, the
resolved parameters, copies of the input tables, and the status of every step.
Generators elsewhere in edms are one-shot; this is the part that remembers.

Usage:
[Run directory]
- new_run_id(): build a run id from the date and a label
- init(): create a run directory and its ledger
- load(): read a run ledger
- save(): write a run ledger
- add_input(): freeze a further input table into a run that is already going

[Status]
- start(): mark a step in progress
- done(): mark a step finished
- skip(): mark a step deliberately not done
- fail(): mark a step failed
- reopen(): return a step to pending
- note(): attach a free-text note to a step
- deviate(): record that a parameter was changed at the bench
- next_steps(): steps whose dependencies are satisfied
- summary(): printable status table
'''
from __future__ import annotations

import datetime
import json
import re
import shutil
from pathlib import Path

from . import template as template_mod

RUN_FILE = 'run.json'
STATE_VERSION = 1
OPEN_STATUSES = ('pending', 'in_progress', 'failed')
CLOSED_STATUSES = ('done', 'skipped')

# Run directory
def _now() -> str:
    '''
    _now(): current local time as an ISO-8601 string to the second.
    '''
    return datetime.datetime.now().isoformat(timespec='seconds')

def _slug(text: str) -> str:
    '''
    _slug(): make a string safe for a directory name.

    Parameters:
    text (str): raw label
    '''
    out = re.sub(r'[^A-Za-z0-9_.-]+', '_', str(text)).strip('_')
    return out or 'run'

def new_run_id(template: str, run_id: str = None, date: str = None) -> str:
    '''
    new_run_id(): build a run id from the date and a label.

    Parameters:
    template (str): template name, used when no label is given
    run_id (str, optional): label for the run (Default: None -> template name)
    date (str, optional): YYYY-MM-DD (Default: None -> today)

    Dependencies: datetime, _slug()
    '''
    day = date or datetime.date.today().isoformat()
    label = _slug(run_id or Path(template).stem)
    if label.startswith(day):
        return label
    return f'{day}_{label}'

def init(spec, tpl=None, dir: str | Path = '.', overwrite: bool = False) -> tuple[Path, dict]:
    '''
    init(): create a run directory and its ledger.

    Freezes the resolved parameters and copies every input table into the run
    directory, so the run still makes sense after the template or the source
    spreadsheets move on.

    Parameters:
    spec (RunSpec): the run spec
    tpl (Template, optional): loaded template (Default: None -> load spec.template)
    dir (str | Path, optional): parent directory for the run (Default: '.')
    overwrite (bool, optional): replace an existing run directory (Default: False)

    Dependencies: shutil, template.load(), new_run_id(), save()
    '''
    if tpl is None:
        tpl = template_mod.load(spec.template)

    run_id = new_run_id(tpl.name, spec.run_id)
    run_dir = Path(dir).expanduser() / run_id

    if run_dir.exists():
        if not overwrite:
            raise FileExistsError(
                f'{run_dir} already exists. Pass overwrite=True to replace it, '
                f'or give a different run_id.')
        shutil.rmtree(run_dir)
    (run_dir / 'inputs').mkdir(parents=True)

    # Resolved params: template defaults, then spec overrides
    unknown = set(spec.params) - set(tpl.params)
    if unknown:
        raise ValueError(
            f"Template '{tpl.name}' has no parameter(s): {', '.join(sorted(unknown))}. "
            f"Known parameters: {', '.join(sorted(tpl.params)) or '(none)'}")
    params = {**tpl.params, **spec.params}

    # Freeze inputs into the run directory
    inputs = {}
    missing = []
    for decl in tpl.inputs:
        if decl.name not in spec.inputs:
            if not decl.optional:
                missing.append(decl.name)
            continue
        src = spec.input_path(decl.name)
        if not src.exists():
            raise FileNotFoundError(f"Input '{decl.name}' not found at {src}")
        frozen = run_dir / 'inputs' / f'{decl.name}{src.suffix}'
        shutil.copy2(src, frozen)
        inputs[decl.name] = {
            'source': str(src),
            'frozen': str(frozen.relative_to(run_dir)),
            'sheet': spec.sheets.get(decl.name),
        }
    if missing:
        shutil.rmtree(run_dir)
        raise ValueError(
            f"Template '{tpl.name}' needs input(s) the spec does not provide: "
            f"{', '.join(missing)}")

    extra = set(spec.inputs) - {i.name for i in tpl.inputs}
    if extra:
        raise ValueError(
            f"Template '{tpl.name}' does not declare input(s): {', '.join(sorted(extra))}")

    state = {
        'state_version': STATE_VERSION,
        'run_id': run_id,
        'created': _now(),
        'updated': _now(),
        'template': {
            'name': tpl.name,
            'version': tpl.version,
            'title': tpl.title,
            'source': str(tpl.source) if tpl.source else None,
        },
        'params': params,
        'inputs': inputs,
        'steps': [
            {
                'id': s.id,
                'title': s.title,
                'day': s.day,
                'depends_on': list(s.depends_on),
                'status': 'pending',
                'started': None,
                'finished': None,
                'notes': [],
                'deviations': [],
                'outputs': [],
            }
            for s in tpl.steps
        ],
        'events': [{'ts': _now(), 'step': None, 'action': 'init',
                    'detail': f'{tpl.name}@{tpl.version}'}],
    }
    save(state, run_dir)
    return run_dir, state

def run_path(run_dir: str | Path) -> Path:
    '''
    run_path(): path to the ledger inside a run directory.

    Parameters:
    run_dir (str | Path): run directory, or the run.json itself
    '''
    p = Path(run_dir).expanduser()
    return p if p.name == RUN_FILE else p / RUN_FILE

def load(run_dir: str | Path) -> dict:
    '''
    load(): read a run ledger.

    Parameters:
    run_dir (str | Path): run directory, or the run.json itself

    Dependencies: json, run_path()
    '''
    path = run_path(run_dir)
    if not path.exists():
        raise FileNotFoundError(
            f'No {RUN_FILE} at {path}. Start a run with: edms work init')
    with open(path, encoding='utf-8') as fh:
        state = json.load(fh)
    version = state.get('state_version')
    if version != STATE_VERSION:
        raise ValueError(
            f'{path}: ledger state_version {version!r} does not match '
            f'{STATE_VERSION} supported by this edms')
    return state

def save(state: dict, run_dir: str | Path) -> Path:
    '''
    save(): write a run ledger.

    Parameters:
    state (dict): ledger
    run_dir (str | Path): run directory, or the run.json itself

    Dependencies: json, run_path()
    '''
    path = run_path(run_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    state['updated'] = _now()
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(state, fh, indent=2)
        fh.write('\n')
    return path

def _template_for(state: dict):
    '''
    _template_for(): load the template a run was started from.

    Falls back to the bundled name when the file it was read from has moved.

    Parameters:
    state (dict): ledger
    '''
    ref = state['template']
    source = ref.get('source')
    return template_mod.load(source if source and Path(source).exists() else ref['name'])

def add_input(run_dir, name: str, pt: str | Path, sheet: str = None,
              replace: bool = False, tpl=None) -> dict:
    '''
    add_input(): freeze a further input table into a run that is already going.

    Inputs are frozen at init so that a finished run still describes what was
    actually used, but some tables are not known on day one - an NGS sample
    sheet, say, which depends on primers chosen after the harvest. This adds
    one, copies it in beside the others, and records where it came from.

    Parameters:
    run_dir (str | Path): run directory
    name (str): input name, as the template declares it
    pt (str | Path): path to the table
    sheet (str, optional): tab to read, for .xlsx inputs
    replace (bool, optional): overwrite an input the run already froze
        (Default: False)
    tpl (Template, optional): loaded template (Default: None -> from the ledger)

    Dependencies: shutil, load(), save(), _template_for()
    '''
    run_dir = run_path(run_dir).parent
    state = load(run_dir)
    if tpl is None:
        tpl = _template_for(state)

    declared = {i.name for i in tpl.inputs}
    if name not in declared:
        raise ValueError(
            f"Template '{tpl.name}' does not declare an input '{name}'. "
            f"Inputs: {', '.join(sorted(declared)) or '(none)'}")
    if name in state['inputs'] and not replace:
        raise ValueError(
            f"Run already froze an input '{name}'. Pass replace=True to swap it, "
            f'which is worth doing only when the first copy was wrong.')

    src = Path(pt).expanduser().resolve()
    if not src.exists():
        raise FileNotFoundError(f"Input '{name}' not found at {src}")

    frozen = run_dir / 'inputs' / f'{name}{src.suffix}'
    frozen.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, frozen)

    action = 'replace_input' if name in state['inputs'] else 'add_input'
    state['inputs'][name] = {
        'source': str(src),
        'frozen': str(frozen.relative_to(run_dir)),
        'sheet': sheet,
    }
    _event(state, None, action, f'{name} <- {src}')
    save(state, run_dir)
    return state

# Status
def _find(state: dict, step_id: str) -> dict:
    '''
    _find(): return the ledger entry for a step.

    Parameters:
    state (dict): ledger
    step_id (str): step identifier
    '''
    for entry in state['steps']:
        if entry['id'] == step_id:
            return entry
    raise KeyError(f"Run has no step '{step_id}'. "
                   f"Steps: {', '.join(e['id'] for e in state['steps'])}")

def _blockers(state: dict, entry: dict) -> list[str]:
    '''
    _blockers(): dependencies of a step that are not finished yet.

    Parameters:
    state (dict): ledger
    entry (dict): step entry
    '''
    return [dep for dep in entry['depends_on']
            if _find(state, dep)['status'] not in CLOSED_STATUSES]

def _event(state: dict, step_id, action: str, detail=None):
    '''
    _event(): append to the run's audit trail.

    Parameters:
    state (dict): ledger
    step_id (str | None): step the event concerns
    action (str): what happened
    detail (optional): free-form detail
    '''
    state['events'].append({'ts': _now(), 'step': step_id,
                            'action': action, 'detail': detail})

def start(run_dir, step: str, force: bool = False) -> dict:
    '''
    start(): mark a step in progress.

    Parameters:
    run_dir (str | Path): run directory
    step (str): step identifier
    force (bool, optional): start even with unfinished dependencies (Default: False)

    Dependencies: load(), save(), _find(), _blockers()
    '''
    state = load(run_dir)
    entry = _find(state, step)
    if entry['status'] == 'done':
        raise ValueError(f"Step '{step}' is already done. Use reopen() to redo it.")
    blockers = _blockers(state, entry)
    if blockers and not force:
        raise ValueError(
            f"Step '{step}' depends on unfinished step(s): {', '.join(blockers)}. "
            f'Pass force=True to start anyway.')
    entry['status'] = 'in_progress'
    entry['started'] = entry['started'] or _now()
    _event(state, step, 'start', 'forced' if blockers else None)
    save(state, run_dir)
    return state

def done(run_dir, step: str, note: str = None, force: bool = False) -> dict:
    '''
    done(): mark a step finished, starting it first if it was still pending.

    Parameters:
    run_dir (str | Path): run directory
    step (str): step identifier
    note (str, optional): note to attach
    force (bool, optional): finish even with unfinished dependencies (Default: False)

    Dependencies: load(), save(), _find(), _blockers()
    '''
    state = load(run_dir)
    entry = _find(state, step)
    if entry['status'] == 'done':
        raise ValueError(
            f"Step '{step}' is already done. Use note() to add to it, or "
            f'reopen() to redo it.')
    blockers = _blockers(state, entry)
    if blockers and not force:
        raise ValueError(
            f"Step '{step}' depends on unfinished step(s): {', '.join(blockers)}. "
            f'Pass force=True to finish anyway.')
    entry['started'] = entry['started'] or _now()
    entry['finished'] = _now()
    entry['status'] = 'done'
    if note:
        entry['notes'].append({'ts': _now(), 'text': note})
    _event(state, step, 'done', note)
    save(state, run_dir)
    return state

def skip(run_dir, step: str, why: str = None) -> dict:
    '''
    skip(): mark a step deliberately not done.

    Parameters:
    run_dir (str | Path): run directory
    step (str): step identifier
    why (str, optional): reason, worth recording

    Dependencies: load(), save(), _find()
    '''
    state = load(run_dir)
    entry = _find(state, step)
    entry['status'] = 'skipped'
    entry['finished'] = _now()
    if why:
        entry['notes'].append({'ts': _now(), 'text': f'skipped: {why}'})
    _event(state, step, 'skip', why)
    save(state, run_dir)
    return state

def fail(run_dir, step: str, why: str = None) -> dict:
    '''
    fail(): mark a step failed, so it stays visible as unfinished work.

    Parameters:
    run_dir (str | Path): run directory
    step (str): step identifier
    why (str, optional): what went wrong

    Dependencies: load(), save(), _find()
    '''
    state = load(run_dir)
    entry = _find(state, step)
    entry['status'] = 'failed'
    if why:
        entry['notes'].append({'ts': _now(), 'text': f'failed: {why}'})
    _event(state, step, 'fail', why)
    save(state, run_dir)
    return state

def reopen(run_dir, step: str, why: str = None) -> dict:
    '''
    reopen(): return a step to pending, e.g. to repeat it.

    Parameters:
    run_dir (str | Path): run directory
    step (str): step identifier
    why (str, optional): reason

    Dependencies: load(), save(), _find()
    '''
    state = load(run_dir)
    entry = _find(state, step)
    entry['status'] = 'pending'
    entry['started'] = None
    entry['finished'] = None
    if why:
        entry['notes'].append({'ts': _now(), 'text': f'reopened: {why}'})
    _event(state, step, 'reopen', why)
    save(state, run_dir)
    return state

def note(run_dir, step: str, text: str) -> dict:
    '''
    note(): attach a free-text note to a step.

    Parameters:
    run_dir (str | Path): run directory
    step (str): step identifier
    text (str): the note

    Dependencies: load(), save(), _find()
    '''
    state = load(run_dir)
    entry = _find(state, step)
    entry['notes'].append({'ts': _now(), 'text': text})
    _event(state, step, 'note', text)
    save(state, run_dir)
    return state

def deviate(run_dir, step: str, param: str, to, why: str = None) -> dict:
    '''
    deviate(): record that a parameter was changed at the bench.

    Kept as structured data rather than prose so that a question like 'every
    run where we dropped cell density' stays answerable later.

    Parameters:
    run_dir (str | Path): run directory
    step (str): step the change applies to
    param (str): parameter name, as the template declares it
    to: the value actually used
    why (str, optional): reason for the change

    Dependencies: load(), save(), _find(), utils.try_parse()
    '''
    from ..utils import try_parse

    state = load(run_dir)
    entry = _find(state, step)
    if param not in state['params']:
        raise KeyError(
            f"Run has no parameter '{param}'. "
            f"Parameters: {', '.join(sorted(state['params'])) or '(none)'}")
    value = try_parse(to) if isinstance(to, str) else to
    record = {'ts': _now(), 'param': param, 'from': state['params'][param],
              'to': value, 'why': why}
    entry['deviations'].append(record)
    state['params'][param] = value
    _event(state, step, 'deviate', f'{param}: {record["from"]} -> {value}')
    save(state, run_dir)
    return state

def record_outputs(run_dir, step: str, paths: list) -> dict:
    '''
    record_outputs(): note the files a step produced.

    Parameters:
    run_dir (str | Path): run directory
    step (str): step identifier
    paths (list): file paths, relative to the run directory where possible

    Dependencies: load(), save(), _find()
    '''
    state = load(run_dir)
    entry = _find(state, step)
    for p in paths:
        if p not in entry['outputs']:
            entry['outputs'].append(str(p))
    _event(state, step, 'outputs', ', '.join(str(p) for p in paths))
    save(state, run_dir)
    return state

def next_steps(state: dict) -> list[dict]:
    '''
    next_steps(): steps that are not finished and whose dependencies are.

    Parameters:
    state (dict): ledger
    '''
    return [e for e in state['steps']
            if e['status'] in OPEN_STATUSES and not _blockers(state, e)]

def progress(state: dict) -> tuple[int, int]:
    '''
    progress(): (closed steps, total steps).

    Parameters:
    state (dict): ledger
    '''
    closed = sum(1 for e in state['steps'] if e['status'] in CLOSED_STATUSES)
    return closed, len(state['steps'])

MARKS = {'pending': '[ ]', 'in_progress': '[~]', 'done': '[x]',
         'skipped': '[-]', 'failed': '[!]'}

def summary(state: dict) -> str:
    '''
    summary(): printable status table for the CLI.

    Parameters:
    state (dict): ledger

    Dependencies: progress(), next_steps()
    '''
    closed, total = progress(state)
    tpl = state['template']
    lines = [
        f"{state['run_id']}  ({tpl['name']}@{tpl['version']})",
        f'  {closed}/{total} steps closed, updated {state["updated"]}',
        '',
    ]
    for e in state['steps']:
        day = '' if e['day'] is None else f'  day {e["day"]:g}'
        line = f'  {MARKS.get(e["status"], "[?]")} {e["id"]:<22} {e["title"]}{day}'
        lines.append(line)
        for dev in e['deviations']:
            lines.append(f'        ~ {dev["param"]}: {dev["from"]} -> {dev["to"]}'
                         + (f' ({dev["why"]})' if dev['why'] else ''))
        for n in e['notes']:
            lines.append(f'        . {n["text"]}')
    upcoming = next_steps(state)
    lines.append('')
    if upcoming:
        lines.append('  next: ' + ', '.join(e['id'] for e in upcoming))
    else:
        lines.append('  next: nothing open - run complete')
    return '\n'.join(lines)
