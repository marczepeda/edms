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

import re
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
    component (str | None): the component this step was included from, or
        None when the template declares it directly
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
    component: str | None = None

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
    components (list): names of the components this template includes, in
        the order they were included
    source (Path | None): file the template was read from
    '''
    name: str
    version: int = 1
    title: str = ''
    description: str = ''
    params: dict = field(default_factory=dict)
    inputs: list[Input] = field(default_factory=list)
    steps: list[Step] = field(default_factory=list)
    components: list[str] = field(default_factory=list)
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

def components_dir() -> Path:
    '''
    components_dir(): directory holding the reusable protocol components.
    '''
    return bundled_dir() / 'components'

def list_components() -> list[str]:
    '''
    list_components(): names of the bundled components.

    Dependencies: components_dir()
    '''
    d = components_dir()
    if not d.is_dir():
        return []
    return sorted(p.stem for p in d.glob('*.toml') if not p.stem.startswith('_'))

def resolve_component(name_or_path: str, relative_to: Path=None) -> Path:
    '''
    resolve_component(): turn a component name or file path into a path.

    A name resolves against the bundled components; a path resolves against
    the including template's own directory first, so somebody's copied
    template can sit next to its own components.

    Parameters:
    name_or_path (str): bundled component name or a path to a .toml file
    relative_to (Path, optional): the including template's path

    Dependencies: components_dir(),list_components()
    '''
    p = Path(name_or_path)
    if p.suffix.lower() == '.toml':
        if relative_to is not None and not p.is_absolute():
            beside = Path(relative_to).resolve().parent / p
            if beside.exists():
                return beside
        if p.exists():
            return p
        raise FileNotFoundError(f'No component file at {p}')
    bundled = components_dir() / f'{name_or_path}.toml'
    if bundled.exists():
        return bundled
    available = ', '.join(list_components()) or '(none)'
    raise FileNotFoundError(
        f"Unknown component '{name_or_path}'. Bundled components: "
        f'{available}. Pass a path to a .toml file to use your own.')

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
def _step_from(entry: dict, where: str, index: int,
               component: str=None) -> Step:
    '''
    _step_from(): one Step out of a parsed TOML entry.

    Parameters:
    entry (dict): the [[steps]] table
    where (str): file the entry came from, for error messages
    index (int): position in the array, for error messages
    component (str, optional): component the step came from (Default: None)
    '''
    if 'id' not in entry:
        raise ValueError(f'{where}: steps[{index}] needs an "id"')
    return Step(
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
        component=component,
    )

def _rename(text: str, mapping: dict) -> str:
    '''
    _rename(): rewrite {param} references after a prefix was applied.

    Parameters:
    text (str): step prose or a check
    mapping (dict): old parameter name -> new name

    Dependencies: re
    '''
    if not text or not mapping:
        return text
    def sub(m):
        name = m.group(1)
        return '{' + mapping.get(name, name) + '}'
    return re.sub(r'\{([A-Za-z_][A-Za-z0-9_]*)\}', sub, text)

def _expand(entry: dict, where: str, index: int, source: Path) -> tuple:
    '''
    _expand(): read one [[steps]] include and return what it contributes.

    An include is a step entry carrying "include" in place of "id": it stands
    where it is written, so one ordered array still describes the protocol.
    A prefix makes a component safe to include twice, by renaming its steps
    and its parameters together and rewriting the {param} references in its
    own prose to match.

    A depends_on entry is prefixed only when it names a step inside the same
    component. Anything else is a reference out to the including template,
    which is how a component attaches to what comes before it.

    The include may also carry its own depends_on, which attaches to the
    component's entry steps -- the ones with nothing inside the component to
    wait for. A component cannot name the step it follows, because it does
    not know which protocol it is in; this is where that is said.

    Parameters:
    entry (dict): the [[steps]] table carrying "include"
    where (str): including file, for error messages
    index (int): position in the array, for error messages
    source (Path): the including template's path

    Dependencies: tomllib,resolve_component(),_step_from(),_rename()
    '''
    name = entry['include']
    prefix = entry.get('prefix', '')
    if not isinstance(prefix, str):
        raise ValueError(f'{where}: steps[{index}] prefix must be a string, '
                         f'got {type(prefix).__name__}')

    path = resolve_component(name, relative_to=source)
    with open(path, 'rb') as fh:
        raw = tomllib.load(fh)
    if 'name' not in raw:
        raise ValueError(f'{path}: component needs a top-level "name"')
    if not raw.get('steps'):
        raise ValueError(f'{path}: component "{raw["name"]}" has no steps')

    own = [_step_from(e, str(path), i, component=raw['name'])
           for i, e in enumerate(raw['steps'])]
    inside = {st.id for st in own}

    params = dict(raw.get('params', {}))
    overrides = dict(entry.get('params', {}))
    unknown = sorted(set(overrides) - set(params))
    if unknown:
        raise ValueError(
            f'{where}: steps[{index}] sets {", ".join(unknown)} for component '
            f'"{raw["name"]}", which does not declare '
            f'{"them" if len(unknown) > 1 else "it"}. Its parameters are: '
            f'{", ".join(sorted(params)) or "(none)"}')
    params.update(overrides)

    renames = {k: prefix + k for k in params} if prefix else {}
    params = {renames.get(k, k): v for k, v in params.items()}

    attach = list(entry.get('depends_on', []))
    for st in own:
        internal = [d for d in st.depends_on if d in inside]
        if attach and not internal:
            st.depends_on = attach + [d for d in st.depends_on if d not in attach]
        st.id = prefix + st.id
        st.depends_on = [prefix + d if d in inside else d
                         for d in st.depends_on]
        if prefix:
            st.text = _rename(st.text, renames)
            st.title = _rename(st.title, renames)
            st.checks = [_rename(c, renames) for c in st.checks]
            st.args = {k: ('$' + renames[v[1:]]
                           if isinstance(v, str) and v.startswith('$')
                           and v[1:] in renames else v)
                       for k, v in st.args.items()}

    inputs = []
    for i, e in enumerate(raw.get('inputs', [])):
        if 'name' not in e:
            raise ValueError(f'{path}: inputs[{i}] needs a "name"')
        inputs.append(Input(
            name=e['name'],
            description=e.get('description', ''),
            required_columns=list(e.get('required_columns', [])),
            optional=bool(e.get('optional', False)),
        ))

    return raw['name'], own, params, inputs

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
    components = []
    from_components = {}
    for i, entry in enumerate(raw.get('steps', [])):
        if 'include' in entry:
            if 'id' in entry:
                raise ValueError(
                    f'{path}: steps[{i}] has both "include" and "id"; an '
                    f'include stands for the component\'s own steps and '
                    f'cannot also be a step')
            cname, own, cparams, cinputs = _expand(entry, str(path), i, path)
            components.append(cname)
            for key, value in cparams.items():
                from_components.setdefault(key, value)
            inputs.extend(c for c in cinputs
                          if c.name not in {x.name for x in inputs})
            steps.extend(own)
            continue
        steps.append(_step_from(entry, str(path), i))

    tpl = Template(
        name=raw['name'],
        version=int(raw.get('version', 1)),
        title=raw.get('title', raw['name']),
        description=raw.get('description', ''),
        # A component's defaults fill in behind the template's own, so
        # including one never silently changes a value the template set.
        params={**from_components, **dict(raw.get('params', {}))},
        inputs=inputs,
        steps=steps,
        components=components,
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
