'''
Module: render.py
Author: Marc Zepeda
Created: 2026-11-12
Description: Render a resolved protocol

Markdown for pasting into a notebook entry, HTML for following at the bench.
Both carry the current status, the notes, and any recorded deviations, so the
rendered protocol and the ledger never disagree.

Usage:
[Rendering]
- to_markdown(): protocol -> Markdown
- to_html(): protocol -> standalone HTML
- write(): render to a file and return its path
'''
from __future__ import annotations

import datetime
import html as html_mod
from pathlib import Path

import pandas as pd

MAX_ROWS = 60

MARKS = {'pending': '&#9744;', 'in_progress': '&#9744;', 'done': '&#9745;',
         'skipped': '&#9647;', 'failed': '&#9888;'}
MD_MARKS = {'pending': '[ ]', 'in_progress': '[~]', 'done': '[x]',
            'skipped': '[-]', 'failed': '[!]'}

CSS = '''
:root {
  --bg: #ffffff; --fg: #1a1a1a; --muted: #666; --line: #e3e3e3;
  --accent: #2b6cb0; --warn: #b7791f; --err: #c53030; --done: #2f855a;
  --card: #fafafa;
}
:root:not([data-theme="light"]) {
  @media (prefers-color-scheme: dark) {
    --bg: #15171a; --fg: #e8e8e8; --muted: #9aa0a6; --line: #2d3137;
    --accent: #7eb3e8; --warn: #e0b35c; --err: #f08a8a; --done: #7fc79b;
    --card: #1c1f23;
  }
}
:root[data-theme="dark"] {
  --bg: #15171a; --fg: #e8e8e8; --muted: #9aa0a6; --line: #2d3137;
  --accent: #7eb3e8; --warn: #e0b35c; --err: #f08a8a; --done: #7fc79b;
  --card: #1c1f23;
}
* { box-sizing: border-box; }
body {
  background: var(--bg); color: var(--fg); margin: 0;
  font: 15px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
}
.wrap { max-width: 980px; margin: 0 auto; padding: 32px 16px 80px; }
h1 { font-size: 1.6rem; margin: 0 0 4px; }
h2 { font-size: 1.1rem; margin: 0; }
.sub { color: var(--muted); font-size: .9rem; margin-bottom: 24px; }
.bar { height: 6px; background: var(--line); border-radius: 3px; overflow: hidden; margin: 16px 0 28px; }
.bar > i { display: block; height: 100%; background: var(--done); }
.meta { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 8px 20px; margin-bottom: 28px; }
.meta div { font-size: .85rem; }
.meta b { color: var(--muted); font-weight: 500; display: block; }
.step { border: 1px solid var(--line); border-radius: 8px; padding: 16px 18px; margin-bottom: 14px; background: var(--card); }
.step.done { opacity: .72; }
.step.failed { border-color: var(--err); }
.hd { display: flex; align-items: baseline; gap: 10px; flex-wrap: wrap; }
.mark { font-size: 1.15rem; line-height: 1; }
.day { color: var(--muted); font-size: .8rem; border: 1px solid var(--line); border-radius: 999px; padding: 1px 8px; }
.id { color: var(--muted); font-size: .78rem; font-family: ui-monospace, SFMono-Regular, Menlo, monospace; }
.text { white-space: pre-wrap; margin: 12px 0 0; }
.checks { margin: 12px 0 0; padding-left: 20px; color: var(--muted); font-size: .9rem; }
.note { font-size: .85rem; color: var(--muted); border-left: 2px solid var(--line); padding-left: 10px; margin-top: 8px; }
.dev { font-size: .85rem; color: var(--warn); border-left: 2px solid var(--warn); padding-left: 10px; margin-top: 8px; }
.err { font-size: .85rem; color: var(--err); margin-top: 10px; }
table { border-collapse: collapse; width: 100%; font-size: .82rem; }
.scroll { overflow-x: auto; margin-top: 12px; border: 1px solid var(--line); border-radius: 6px; }
th, td { text-align: left; padding: 5px 9px; border-bottom: 1px solid var(--line); white-space: nowrap; }
th { color: var(--muted); font-weight: 600; position: sticky; top: 0; background: var(--card); }
tr:last-child td { border-bottom: 0; }
.tname { font-size: .8rem; color: var(--accent); margin-top: 14px; font-weight: 600; }
.more { font-size: .78rem; color: var(--muted); padding: 6px 9px; }
@media print {
  .step { break-inside: avoid; border-color: #ccc; background: #fff; }
  body { background: #fff; color: #000; }
}
'''

def _esc(x) -> str:
    '''
    _esc(): HTML-escape a value.

    Parameters:
    x: value to escape
    '''
    return html_mod.escape('' if x is None else str(x))

def _table_html(df: pd.DataFrame) -> str:
    '''
    _table_html(): render a DataFrame as an HTML table, truncated.

    Parameters:
    df (pd.DataFrame): table to render
    '''
    shown = df.head(MAX_ROWS)
    head = ''.join(f'<th>{_esc(c)}</th>' for c in shown.columns)
    rows = ''.join(
        '<tr>' + ''.join(f'<td>{_esc(v)}</td>' for v in row) + '</tr>'
        for row in shown.itertuples(index=False, name=None))
    more = ''
    if len(df) > MAX_ROWS:
        more = (f'<div class="more">{len(df) - MAX_ROWS} more row(s) - '
                f'see the .xlsx in outputs/</div>')
    return f'<div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{rows}</tbody></table>{more}</div>'

def _md_table(df: pd.DataFrame) -> list[str]:
    '''
    _md_table(): render a DataFrame as GitHub-flavoured Markdown rows.

    Written out here rather than using DataFrame.to_markdown(), which needs
    tabulate; edms does not depend on it.

    Parameters:
    df (pd.DataFrame): table to render
    '''
    def cell(v):
        return '' if v is None else str(v).replace('|', r'\|').replace('\n', ' ')

    cols = [cell(c) for c in df.columns]
    if not cols:
        return []
    out = ['| ' + ' | '.join(cols) + ' |',
           '| ' + ' | '.join('---' for _ in cols) + ' |']
    for row in df.itertuples(index=False, name=None):
        out.append('| ' + ' | '.join(cell(v) for v in row) + ' |')
    return out

def to_html(protocol, theme: str = None) -> str:
    '''
    to_html(): render a resolved protocol as standalone HTML.

    Parameters:
    protocol (ResolvedProtocol): from plan.resolve()
    theme (str, optional): 'light' or 'dark' to pin the theme (Default: None -> follow the viewer)

    Dependencies: datetime, html, pandas
    '''
    closed = sum(1 for s in protocol.steps if s.status in ('done', 'skipped'))
    total = len(protocol.steps) or 1
    pct = round(100 * closed / total)

    meta = [
        ('Run', protocol.run_id),
        ('Template', f"{protocol.template['name']}@{protocol.template['version']}"),
        ('Progress', f'{closed}/{len(protocol.steps)} steps'),
        ('Generated', datetime.datetime.now().strftime('%Y-%m-%d %H:%M')),
    ]
    meta_html = ''.join(f'<div><b>{_esc(k)}</b>{_esc(v)}</div>' for k, v in meta)

    params_html = ''
    if protocol.params:
        rows = ''.join(f'<tr><td>{_esc(k)}</td><td>{_esc(v)}</td></tr>'
                       for k, v in protocol.params.items())
        params_html = (f'<div class="step"><div class="hd"><h2>Parameters</h2></div>'
                       f'<div class="scroll"><table><thead><tr><th>Parameter</th>'
                       f'<th>Value</th></tr></thead><tbody>{rows}</tbody></table></div></div>')

    steps_html = []
    for s in protocol.steps:
        parts = [f'<div class="hd">'
                 f'<span class="mark">{MARKS.get(s.status, "&#9744;")}</span>'
                 f'<h2>{_esc(s.title)}</h2>'
                 f'<span class="id">{_esc(s.id)}</span>']
        if s.day is not None:
            parts.append(f'<span class="day">day {s.day:g}</span>')
        if s.duration:
            parts.append(f'<span class="day">{_esc(s.duration)}</span>')
        parts.append('</div>')

        if s.text:
            parts.append(f'<p class="text">{_esc(s.text)}</p>')
        if s.checks:
            items = ''.join(f'<li>{_esc(c)}</li>' for c in s.checks)
            parts.append(f'<ul class="checks">{items}</ul>')
        for name, table in s.tables.items():
            parts.append(f'<div class="tname">{_esc(name)}</div>')
            parts.append(_table_html(table))
        if s.compute_skipped:
            parts.append(f'<div class="note">Not calculated: {_esc(s.compute_skipped)}</div>')
        if s.error:
            parts.append(f'<div class="err">Could not compute: {_esc(s.error)}</div>')
        for dev in s.deviations:
            why = f' &mdash; {_esc(dev["why"])}' if dev.get('why') else ''
            parts.append(f'<div class="dev">{_esc(dev["param"])}: '
                         f'{_esc(dev["from"])} &rarr; {_esc(dev["to"])}{why}</div>')
        for n in s.notes:
            parts.append(f'<div class="note">{_esc(n["text"])} '
                         f'<span class="id">{_esc(n["ts"])}</span></div>')

        steps_html.append(f'<div class="step {s.status}">' + ''.join(parts) + '</div>')

    attr = f' data-theme="{theme}"' if theme in ('light', 'dark') else ''
    desc = f'<p class="sub">{_esc(protocol.description)}</p>' if protocol.description else ''

    return f'''<!doctype html>
<html lang="en"{attr}>
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{_esc(protocol.run_id)}</title>
<style>{CSS}</style>
</head>
<body>
<div class="wrap">
<h1>{_esc(protocol.title)}</h1>
{desc}
<div class="meta">{meta_html}</div>
<div class="bar"><i style="width:{pct}%"></i></div>
{params_html}
{''.join(steps_html)}
</div>
</body>
</html>
'''

def to_markdown(protocol) -> str:
    '''
    to_markdown(): render a resolved protocol as Markdown, for pasting into a
    notebook entry.

    Parameters:
    protocol (ResolvedProtocol): from plan.resolve()

    Dependencies: datetime, pandas
    '''
    closed = sum(1 for s in protocol.steps if s.status in ('done', 'skipped'))
    lines = [
        f'# {protocol.title}',
        '',
        f"**Run:** `{protocol.run_id}`  ",
        f"**Template:** `{protocol.template['name']}@{protocol.template['version']}`  ",
        f'**Progress:** {closed}/{len(protocol.steps)} steps  ',
        f"**Generated:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
        '',
    ]
    if protocol.description:
        lines += [protocol.description, '']

    if protocol.params:
        lines += ['## Parameters', '', '| Parameter | Value |', '| --- | --- |']
        lines += [f'| {k} | {v} |' for k, v in protocol.params.items()]
        lines += ['']

    for s in protocol.steps:
        head = f'## {MD_MARKS.get(s.status, "[ ]")} {s.title}'
        bits = [f'`{s.id}`']
        if s.day is not None:
            bits.append(f'day {s.day:g}')
        if s.duration:
            bits.append(s.duration)
        lines += [head, '', ' \u00b7 '.join(bits), '']
        if s.text:
            lines += [s.text.strip(), '']
        if s.checks:
            lines += [f'- [ ] {c}' for c in s.checks] + ['']
        for name, table in s.tables.items():
            lines += [f'**{name}**', '']
            lines += _md_table(table.head(MAX_ROWS))
            if len(table) > MAX_ROWS:
                lines += ['', f'_{len(table) - MAX_ROWS} more row(s) in outputs/{s.id}.xlsx_']
            lines += ['']
        if s.compute_skipped:
            lines += [f'> Not calculated: {s.compute_skipped}', '']
        if s.error:
            lines += [f'> Could not compute: {s.error}', '']
        for dev in s.deviations:
            why = f' - {dev["why"]}' if dev.get('why') else ''
            lines += [f'> Deviation: {dev["param"]} {dev["from"]} -> {dev["to"]}{why}', '']
        for n in s.notes:
            lines += [f'> {n["text"]} _({n["ts"]})_', '']

    return '\n'.join(lines).rstrip() + '\n'

def write(protocol, file: str | Path = None, format: str = 'html',
          theme: str = None) -> Path:
    '''
    write(): render a protocol to a file and return its path.

    Parameters:
    protocol (ResolvedProtocol): from plan.resolve()
    file (str | Path, optional): output path (Default: None -> protocol.html or
        protocol.md inside the run directory)
    format (str, optional): 'html' or 'md' (Default: 'html')
    theme (str, optional): 'light' or 'dark' to pin an HTML theme

    Dependencies: to_html(), to_markdown()
    '''
    if format not in ('html', 'md'):
        raise ValueError(f"format must be 'html' or 'md', got '{format}'")
    out = Path(file).expanduser() if file else protocol.run_dir / f'protocol.{format}'
    out.parent.mkdir(parents=True, exist_ok=True)
    text = to_html(protocol, theme=theme) if format == 'html' else to_markdown(protocol)
    out.write_text(text, encoding='utf-8')
    return out
