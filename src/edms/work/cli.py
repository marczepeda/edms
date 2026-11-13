'''
src/edms/work/cli.py            Command Line Interface for EDMS Work module
├── __init__.py                 Initializer
├── template.py                 Protocol templates
├── spec.py                     Run specs
├── plan.py                     Resolve a template and a run into a protocol
├── run.py                      Run state
└── render.py                   Render a resolved protocol

Usage:
[Template commands]
- templates(): list the bundled protocol templates
- show(): print a template's parameters, inputs, and steps
- where(): print the path to a template file

[Run commands]
- example(): write an annotated starting run spec
- init(): start a run from a spec
- input(): add an input table to a run that is already going
- plan(): resolve and render the protocol for a run
- status(): print a run's status

[Status commands]
- start(), done(), skip(), fail(), reopen(), note(), deviate()

[Main subparser method]
- add_subparser(): Attach all work-related subparsers to the top-level CLI.
'''
import argparse
import functools
from pathlib import Path

from rich import print as rprint

from . import plan as plan_mod
from . import render as render_mod
from . import run as run_mod
from . import spec as spec_mod
from . import template as template_mod

def _clean(func):
    '''
    _clean(): turn the errors a user can actually cause into a one-line message
    instead of a traceback. Anything unexpected still raises in full.

    Parameters:
    func (callable): CLI command function
    '''
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except (ValueError, KeyError, FileNotFoundError, FileExistsError) as e:
            # KeyError's str() adds quotes around the message; drop them
            message = e.args[0] if isinstance(e, KeyError) and e.args else e
            raise SystemExit(f'edms work {func.__name__}: {message}') from None
    return wrapper

# Template commands
def templates():
    '''
    templates(): list the bundled protocol templates
    '''
    names = template_mod.list_templates()
    if not names:
        rprint('[yellow]No bundled templates found.[/yellow]')
        return
    rprint('[green]Bundled protocol templates:[/green]')
    for name in names:
        tpl = template_mod.load(name)
        rprint(f'  [cyan]{name}[/cyan]  v{tpl.version}  {tpl.title}')
    rprint('\nInspect one with: [cyan]edms work show <name>[/cyan]')

def show(template: str):
    '''
    show(): print a template's parameters, inputs, and steps

    Parameters:
    template (str): bundled template name or path to a .toml template
    '''
    tpl = template_mod.load(template)
    rprint(f'[green]{tpl.title}[/green]  ([cyan]{tpl.name}[/cyan] v{tpl.version})')
    if tpl.description:
        rprint(tpl.description.strip())
    rprint(f'\n[green]Source:[/green] {tpl.source}')

    rprint('\n[green]Inputs[/green]')
    for i in tpl.inputs:
        tag = ' [dim](optional)[/dim]' if i.optional else ''
        rprint(f'  [cyan]{i.name}[/cyan]{tag}  {i.description}')
        if i.required_columns:
            rprint(f'    columns: {", ".join(i.required_columns)}')

    rprint('\n[green]Parameters[/green]')
    for k, v in tpl.params.items():
        rprint(f'  [cyan]{k}[/cyan] = {v!r}')

    rprint('\n[green]Steps[/green]')
    for s in tpl.steps:
        day = '' if s.day is None else f'  day {s.day:g}'
        rprint(f'  [cyan]{s.id}[/cyan]{day}  {s.title}')
        if s.compute:
            rprint(f'    computes: {s.compute}')
        if s.depends_on:
            rprint(f'    after: {", ".join(s.depends_on)}')

def where(template: str):
    '''
    where(): print the path to a template file, for copying it

    Parameters:
    template (str): bundled template name or path to a .toml template
    '''
    print(template_mod.resolve_path(template))

# Run commands
def example(dir: str = '.', file: str = 'run.toml'):
    '''
    example(): write an annotated starting run spec

    Parameters:
    dir (str, optional): directory to write into (Default: '.')
    file (str, optional): file name (Default: 'run.toml')
    '''
    out = spec_mod.example(dir=dir, file=file)
    rprint(f'[green]Wrote example run spec to[/green] {out}')
    rprint(f'Edit it, then: [cyan]edms work init -i {out.name}[/cyan]')

def init(spec: str = None, template: str = None, input: list = None, set: list = None,
         sheet: list = None, run_id: str = None, dir: str = '.',
         overwrite: bool = False, plan: bool = True):
    '''
    init(): start a run from a spec

    Parameters:
    spec (str, optional): path to a run spec .toml
    template (str, optional): template name or path; required without a spec file
    input (list, optional): 'name=path' pairs, layered over the spec file
    set (list, optional): 'param=value' pairs, layered over the spec file
    sheet (list, optional): 'name=sheet' pairs for .xlsx inputs
    run_id (str, optional): short label for the run directory
    dir (str, optional): parent directory for the run (Default: '.')
    overwrite (bool, optional): replace an existing run directory (Default: False)
    plan (bool, optional): render the protocol straight away (Default: True)

    Dependencies: spec.load(), spec.from_args(), spec.merge(), run.init()
    '''
    if spec is None and template is None:
        raise SystemExit('Pass a run spec with -i, or a template with -t. '
                         'Write a starting spec with: edms work example')

    overrides = spec_mod.from_args(template=template or '', input=input, set=set,
                                   sheet=sheet, run_id=run_id)
    run_spec = spec_mod.merge(spec_mod.load(spec), overrides) if spec else overrides

    tpl = template_mod.load(run_spec.template)
    run_dir, state = run_mod.init(run_spec, tpl=tpl, dir=dir, overwrite=overwrite)
    rprint(f'[green]Started run[/green] {run_dir}')
    rprint(run_mod.summary(state))

    if plan:
        _plan(run_dir, tpl=tpl)

def input(run_dir: str, input: list = None, sheet: list = None,
          replace: bool = False, plan: bool = True):
    '''
    input(): add an input table to a run that is already going

    For tables that are not known on day one, such as an NGS sample sheet that
    depends on primers chosen after the harvest.

    Parameters:
    run_dir (str): run directory
    input (list, optional): 'name=path' pairs; repeatable
    sheet (list, optional): 'name=sheet' pairs for .xlsx inputs
    replace (bool, optional): overwrite an input the run already froze
    plan (bool, optional): re-render the protocol afterwards (Default: True)

    Dependencies: spec._pairs(), run.add_input()
    '''
    pairs = spec_mod._pairs(input, parse=False)
    if not pairs:
        raise SystemExit('Pass at least one -I NAME=PATH')
    sheets = spec_mod._pairs(sheet, parse=False)

    for name, pt in pairs.items():
        state = run_mod.add_input(run_dir, name=name, pt=pt,
                                  sheet=sheets.get(name), replace=replace)
        rprint(f'[green]Added input[/green] {name} [dim]<-[/dim] {pt}')

    if plan:
        _plan(run_dir)

def _plan(run_dir, tpl=None, format: str = 'html', file: str = None,
          strict: bool = False, compute: bool = True, theme: str = None):
    '''
    _plan(): resolve and render, shared by init() and plan()

    Parameters:
    run_dir (str | Path): run directory
    tpl (Template, optional): already-loaded template
    format (str, optional): 'html' or 'md' (Default: 'html')
    file (str, optional): output path
    strict (bool, optional): raise on the first compute failure (Default: False)
    compute (bool, optional): run compute steps (Default: True)
    theme (str, optional): 'light' or 'dark' to pin an HTML theme
    '''
    protocol = plan_mod.resolve(run_dir, tpl=tpl, strict=strict, compute=compute)
    out = render_mod.write(protocol, file=file, format=format, theme=theme)
    for step_id, message in protocol.errors:
        rprint(f'[yellow]Could not compute {step_id}:[/yellow] {message}')
    skipped = [s.id for s in protocol.steps if s.compute_skipped]
    if skipped:
        rprint(f'[yellow]Not calculated yet:[/yellow] {", ".join(skipped)} '
               f'(waiting on optional inputs)')
    rprint(f'[green]Wrote protocol to[/green] {out}')
    return out

def plan(run_dir: str, format: str = 'html', file: str = None, strict: bool = False,
         compute: bool = True, theme: str = None):
    '''
    plan(): resolve and render the protocol for a run

    Re-runnable at any point: it picks up parameter deviations and any inputs
    added since the run started.

    Parameters:
    run_dir (str): run directory
    format (str, optional): 'html' or 'md' (Default: 'html')
    file (str, optional): output path (Default: None -> protocol.html in the run)
    strict (bool, optional): raise on the first compute failure (Default: False)
    compute (bool, optional): run compute steps (Default: True)
    theme (str, optional): 'light' or 'dark' to pin an HTML theme

    Dependencies: plan.resolve(), render.write()
    '''
    _plan(run_dir, format=format, file=file, strict=strict, compute=compute, theme=theme)

def status(run_dir: str):
    '''
    status(): print a run's status

    Parameters:
    run_dir (str): run directory

    Dependencies: run.load(), run.summary()
    '''
    print(run_mod.summary(run_mod.load(run_dir)))

# Status commands
def start(run_dir: str, step: str, force: bool = False):
    '''
    start(): mark a step in progress

    Parameters:
    run_dir (str): run directory
    step (str): step identifier
    force (bool, optional): start even with unfinished dependencies
    '''
    print(run_mod.summary(run_mod.start(run_dir, step, force=force)))

def done(run_dir: str, step: str, note: str = None, force: bool = False):
    '''
    done(): mark a step finished

    Parameters:
    run_dir (str): run directory
    step (str): step identifier
    note (str, optional): note to attach
    force (bool, optional): finish even with unfinished dependencies
    '''
    print(run_mod.summary(run_mod.done(run_dir, step, note=note, force=force)))

def skip(run_dir: str, step: str, why: str = None):
    '''
    skip(): mark a step deliberately not done

    Parameters:
    run_dir (str): run directory
    step (str): step identifier
    why (str, optional): reason
    '''
    print(run_mod.summary(run_mod.skip(run_dir, step, why=why)))

def fail(run_dir: str, step: str, why: str = None):
    '''
    fail(): mark a step failed

    Parameters:
    run_dir (str): run directory
    step (str): step identifier
    why (str, optional): what went wrong
    '''
    print(run_mod.summary(run_mod.fail(run_dir, step, why=why)))

def reopen(run_dir: str, step: str, why: str = None):
    '''
    reopen(): return a step to pending, e.g. to repeat it

    Parameters:
    run_dir (str): run directory
    step (str): step identifier
    why (str, optional): reason
    '''
    print(run_mod.summary(run_mod.reopen(run_dir, step, why=why)))

def note(run_dir: str, step: str, note: str):
    '''
    note(): attach a free-text note to a step

    Parameters:
    run_dir (str): run directory
    step (str): step identifier
    note (str): the note
    '''
    print(run_mod.summary(run_mod.note(run_dir, step, text=note)))

def deviate(run_dir: str, step: str, param: str, value: str, why: str = None):
    '''
    deviate(): record that a parameter was changed at the bench

    Parameters:
    run_dir (str): run directory
    step (str): step the change applies to
    param (str): parameter name, as the template declares it
    value (str): the value actually used
    why (str, optional): reason for the change
    '''
    state = run_mod.deviate(run_dir, step, param=param, to=value, why=why)
    print(run_mod.summary(state))
    rprint('[yellow]Re-render the protocol to pick this up:[/yellow] '
           f'[cyan]edms work plan {run_dir}[/cyan]')

# Main subparser method
def add_subparser(subparsers, formatter_class=None):
    """
    add_subparser(): Attach all work-related subparsers to the top-level CLI.

    Parameters:
    subparsers (argparse._SubParsersAction): The subparsers object to attach the work subparsers to.
    formatter_class (type, optional): The formatter class to use for the subparsers.
    """
    if formatter_class is None:
        formatter_class = argparse.HelpFormatter

    parser_work = subparsers.add_parser(
        "work", help="Plan experiments and track them",
        description="Plan experiments from protocol templates and track them as they run",
        formatter_class=formatter_class)
    sub = parser_work.add_subparsers()

    def add(name, help, func):
        p = sub.add_parser(name, help=help, description=help, formatter_class=formatter_class)
        p.set_defaults(func=_clean(func))
        return p

    # Template commands
    add("templates", "List the bundled protocol templates", templates)

    p_show = add("show", "Print a template's parameters, inputs, and steps", show)
    p_show.add_argument("template", type=str, help="Template name or path to a .toml template")

    p_where = add("where", "Print the path to a template file", where)
    p_where.add_argument("template", type=str, help="Template name or path to a .toml template")

    # Run commands
    p_example = add("example", "Write an annotated starting run spec", example)
    p_example.add_argument("-d", "--dir", type=str, default='.', help="Directory to write into (Default: .)")
    p_example.add_argument("-f", "--file", type=str, default='run.toml', help="File name (Default: run.toml)")

    p_init = add("init", "Start a run from a spec", init)
    p_init.add_argument("-i", "--spec", type=str, help="Run spec .toml (see: edms work example)")
    p_init.add_argument("-t", "--template", type=str, help="Template name or path; required without -i")
    p_init.add_argument("-I", "--input", action="append", metavar="NAME=PATH", help="Input table; repeatable")
    p_init.add_argument("-S", "--set", action="append", metavar="PARAM=VALUE", help="Parameter override; repeatable")
    p_init.add_argument("-H", "--sheet", action="append", metavar="NAME=SHEET", help="Tab to read for an .xlsx input; repeatable")
    p_init.add_argument("-r", "--run_id", type=str, help="Short label for the run directory")
    p_init.add_argument("-d", "--dir", type=str, default='.', help="Parent directory for the run (Default: .)")
    p_init.add_argument("-o", "--overwrite", action="store_true", help="Replace an existing run directory")
    p_init.add_argument("-n", "--no_plan", dest="plan", action="store_false", default=True, help="Don't render the protocol straight away")

    p_input = add("input", "Add an input table to a run that is already going", input)
    p_input.add_argument("run_dir", type=str, help="Run directory")
    p_input.add_argument("-I", "--input", action="append", metavar="NAME=PATH", required=True, help="Input table; repeatable")
    p_input.add_argument("-H", "--sheet", action="append", metavar="NAME=SHEET", help="Tab to read for an .xlsx input; repeatable")
    p_input.add_argument("-R", "--replace", action="store_true", help="Overwrite an input the run already froze")
    p_input.add_argument("-n", "--no_plan", dest="plan", action="store_false", default=True, help="Don't re-render the protocol afterwards")

    p_plan = add("plan", "Resolve and render the protocol for a run", plan)
    p_plan.add_argument("run_dir", type=str, help="Run directory")
    p_plan.add_argument("-f", "--format", type=str, choices=['html', 'md'], default='html', help="Output format (Default: html)")
    p_plan.add_argument("-o", "--file", type=str, help="Output path (Default: protocol.html in the run directory)")
    p_plan.add_argument("-s", "--strict", action="store_true", help="Stop at the first calculation that fails")
    p_plan.add_argument("-n", "--no_compute", dest="compute", action="store_false", default=True, help="Render prose only, skipping every calculation")
    p_plan.add_argument("-T", "--theme", type=str, choices=['light', 'dark'], help="Pin the HTML theme instead of following the viewer")

    p_status = add("status", "Print a run's status", status)
    p_status.add_argument("run_dir", type=str, help="Run directory")

    # Status commands
    p_start = add("start", "Mark a step in progress", start)
    p_start.add_argument("run_dir", type=str, help="Run directory")
    p_start.add_argument("-s", "--step", type=str, required=True, help="Step identifier")
    p_start.add_argument("-F", "--force", action="store_true", help="Start even with unfinished dependencies")

    p_done = add("done", "Mark a step finished", done)
    p_done.add_argument("run_dir", type=str, help="Run directory")
    p_done.add_argument("-s", "--step", type=str, required=True, help="Step identifier")
    p_done.add_argument("-n", "--note", type=str, help="Note to attach")
    p_done.add_argument("-F", "--force", action="store_true", help="Finish even with unfinished dependencies")

    for name, help_text, func in (
        ("skip", "Mark a step deliberately not done", skip),
        ("fail", "Mark a step failed", fail),
        ("reopen", "Return a step to pending, e.g. to repeat it", reopen),
    ):
        p = add(name, help_text, func)
        p.add_argument("run_dir", type=str, help="Run directory")
        p.add_argument("-s", "--step", type=str, required=True, help="Step identifier")
        p.add_argument("-w", "--why", type=str, help="Reason, worth recording")

    p_note = add("note", "Attach a free-text note to a step", note)
    p_note.add_argument("run_dir", type=str, help="Run directory")
    p_note.add_argument("-s", "--step", type=str, required=True, help="Step identifier")
    p_note.add_argument("-n", "--note", type=str, required=True, help="The note")

    p_dev = add("deviate", "Record that a parameter was changed at the bench", deviate)
    p_dev.add_argument("run_dir", type=str, help="Run directory")
    p_dev.add_argument("-s", "--step", type=str, required=True, help="Step the change applies to")
    p_dev.add_argument("-p", "--param", type=str, required=True, help="Parameter name, as the template declares it")
    p_dev.add_argument("-v", "--value", type=str, required=True, help="The value actually used")
    p_dev.add_argument("-w", "--why", type=str, help="Reason for the change")
