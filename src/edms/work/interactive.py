'''
Module: interactive.py
Author: Marc Zepeda
Created: 2026-11-14
Description: Data and assets that make a rendered protocol recalculate

A protocol written to HTML is a record, and a record you cannot poke at is
half of what you want at the bench. The other half is 'what would this look
like in a 24-well', asked at the moment it matters rather than back at a
terminal.

The arithmetic therefore exists twice: once in Python, where it is
authoritative, and once in static/calc.js, where it runs in the page. That is
a real cost and it is paid deliberately, because the alternative is a server
the protocol cannot be emailed away from. Three things keep the copy honest:

- Every constant is serialised out of the Python objects here, so no number
  is written down twice. calc.js has no vessel table of its own.
- live_targets() reads back which compute targets calc.js actually implements,
  so a table whose calculation was never ported renders static and says so
  rather than quietly showing a stale number.
- tests/work/test_calc_js.py runs calc.js under node against the Python.

Nothing here writes to a run. A protocol edited in the browser is a scratch
calculation; keeping a change means recording it with edms work deviate, and
the page builds that command for you.

Usage:
[Assets]
- asset(): read one of the bundled static files
- live_targets(): compute targets calc.js implements

[Payload]
- constants(): the numbers calc.js needs, out of the Python objects
- payload(): a resolved protocol as the data the page recalculates from
'''
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

STATIC = Path(__file__).parent / 'static'

def asset(name: str) -> str:
    '''
    asset(): read one of the bundled static files.

    Parameters:
    name (str): file name within edms/work/static
    '''
    if '/' in name or '\\' in name or name in ('', '.', '..'):
        raise FileNotFoundError(f'No bundled asset {name!r}: a name, not a path')
    path = (STATIC / name).resolve()
    if path.parent != STATIC.resolve() or not path.is_file():
        raise FileNotFoundError(f'No bundled asset {name!r} in {STATIC}')
    return path.read_text(encoding='utf-8')

def live_targets() -> set[str]:
    '''
    live_targets(): the compute targets calc.js implements.

    Read out of the file rather than listed here, so the two cannot disagree
    about what recalculates. A target absent from this set renders as a static
    table with a note saying so.

    Dependencies: re,asset()
    '''
    return set(re.findall(r"C\['([a-zA-Z0-9_.]+:[a-zA-Z0-9_]+)'\]", asset('calc.js')))

def constants() -> dict:
    '''
    constants(): the numbers calc.js needs, taken from the Python objects.

    Dependencies: edms.work.plates,edms.bio.clone,edms.bio.culture
    '''
    from . import plates
    from ..bio import clone, culture

    return {
        'CONFLUENT_CELLS': dict(culture.CONFLUENT_CELLS),
        'CONFLUENT_PER_ML': dict(culture.CONFLUENT_PER_ML),
        'vessels': {name: dict(v) for name, v in plates._VESSELS.items()},
        'aliases': dict(plates._ALIASES),
        # The well grid, so calc.js can name a well without a table of its own
        'grid': {name: list(shape) for name, shape in plates._GRID.items()},
        'ROW_LETTERS': plates._ROW_LETTERS,
        'SEED_PER_CM2': plates.SEED_PER_CM2,
        'VIRUS_SEED_PER_CM2': plates.VIRUS_SEED_PER_CM2,
        'LAYOUT_NOTE': plates.LAYOUT_NOTE,
        'PE3_ANCHOR': plates.PE3_ANCHOR,
        'VIRUS_ANCHOR': plates.VIRUS_ANCHOR,
        'MIN_PIPETTE_UL': plates.MIN_PIPETTE_UL,
        'DILUTION_STEPS': list(plates.DILUTION_STEPS),
        'WIDE_SPLIT': plates.WIDE_SPLIT,
        'EXPANSION_SKIP_DEFAULT': list(plates.EXPANSION_SKIP_DEFAULT),
        'NG_PER_PMOL_KB': clone.NG_PER_PMOL_KB,
    }

def _jsonable(value):
    '''
    _jsonable(): a cell or parameter as something json can carry.

    Parameters:
    value: any cell from a table or params block
    '''
    if value is None or isinstance(value, (bool, str, int)):
        return value
    if isinstance(value, float):
        # NaN and the infinities have no JSON spelling, and a cell that reads
        # 'nan' in a browser is worse than an empty one
        if value != value or abs(value) == float('inf'):
            return None
        return value
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    item = getattr(value, 'item', None)
    if callable(item):
        try:
            return _jsonable(item())
        except (TypeError, ValueError):
            pass
    return str(value)

def table_data(df: pd.DataFrame) -> dict:
    '''
    table_data(): a DataFrame as columns and row objects.

    Parameters:
    df (DataFrame): any table from a resolved step or an input

    Dependencies: pandas,_jsonable()
    '''
    columns = [str(c) for c in df.columns]
    rows = [{str(c): _jsonable(row[c]) for c in df.columns}
            for _, row in df.iterrows()]
    return {'columns': columns, 'rows': rows}

def payload(protocol, inputs: dict=None) -> dict:
    '''
    payload(): a resolved protocol as the data the page recalculates from.

    Parameters:
    protocol (ResolvedProtocol): from plan.resolve()
    inputs (dict, optional): name -> DataFrame, the run's loaded inputs
        (Default: None, meaning the page can still edit params)

    Dependencies: pandas,constants(),live_targets(),table_data()
    '''
    from . import plates

    live = live_targets()
    steps = []
    for s in protocol.steps:
        if not s.compute or not s.tables:
            continue
        args = getattr(s, 'args', {}) or {}
        steps.append({'id': s.id,
                      'compute': s.compute,
                      'args': _jsonable(args),
                      # A float literal in the step's own args, for the same
                      # reason as 'floats' below
                      'float_args': sorted(k for k, v in args.items()
                                           if isinstance(v, float)),
                      'tables': [str(n) for n in s.tables],
                      'live': s.compute in live})

    # The row cap the static render uses. Without it the page agrees with
    # itself only until you touch something: a 352-row plate map arrives cut
    # to MAX_ROWS with a note pointing at the .xlsx, and the first recompute
    # would quietly replace it with all 352 and drop the note.
    from .render import MAX_ROWS

    return {'run_id': protocol.run_id,
            'max_rows': MAX_ROWS,
            'consts': constants(),
            'formats': plates.names(),
            'params': {k: _jsonable(v) for k, v in protocol.params.items()},
            # Which parameters are Python floats. JSON cannot tell 1.0 from 1
            # and JavaScript has one number type, so a master mix headed
            # 'uL, MM (3 x 1.1)' would come back as '... x 1)' if the template
            # happened to write 1.0. The page needs the distinction only for
            # the handful of values that get printed into a label.
            'floats': sorted(k for k, v in protocol.params.items()
                             if isinstance(v, float)),
            'inputs': {name: table_data(df)
                       for name, df in (inputs or {}).items()},
            'steps': steps}

def scripts(protocol, inputs: dict=None) -> str:
    '''
    scripts(): the script tags that make a rendered protocol recalculate.

    Parameters:
    protocol (ResolvedProtocol): from plan.resolve()
    inputs (dict, optional): name -> DataFrame, the run's loaded inputs

    Dependencies: json,asset(),payload()
    '''
    data = json.dumps(payload(protocol, inputs), allow_nan=False)
    data = data.replace('</', r'<\/')
    return ('<script id="edms-data" type="application/json">' + data + '</script>\n'
            '<script>\n' + asset('calc.js') + '\n</script>\n'
            '<script>\n' + asset('plate.js') + '\n</script>\n'
            '<script>\n' + asset('app.js') + '\n</script>')
