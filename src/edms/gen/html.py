from __future__ import annotations

from pathlib import Path
from itertools import groupby
import html
import re
from typing import Iterable, Optional
import shutil
import importlib.resources as pkg_resources
import edms.resources.icon as icon_pkg
from edms.gen.tidy import natural_key
from edms.gen.plot import re_un_cap

_TITLE_RE = re.compile(r"<title[^>]*>(.*?)</title>", re.IGNORECASE | re.DOTALL)

def _extract_html_title(path: Path, max_bytes: int = 200_000) -> Optional[str]:
    """Best-effort: read first max_bytes and extract <title>...</title>."""
    try:
        data = path.read_bytes()[:max_bytes]
        text = data.decode("utf-8", errors="ignore")
        m = _TITLE_RE.search(text)
        if not m:
            return None
        title = re.sub(r"\s+", " ", m.group(1)).strip()
        return title or None
    except Exception:
        return None

def make_html_index(
    dir: str | Path,
    file: str | Path = "index.html",
    recursive: bool = False,
    exclude: Iterable[str] = ("index.html", "index.pdf"),
    sort: str = "title",  # "title" | "name" | "mtime"
    preview: bool = True,
    grid_cols: int = 3,
    image_types: list[str] | None = None,
    preview_height_px: int = 900,
    icon: str = "python",
) -> Path:
    """
    make_html_index(): Create an index HTML that links to other HTML files in `dir`.
    - Uses <title> from each HTML file if available; falls back to stem/filename.
    - Makes titles into buttons.
    - Optionally embeds a preview iframe that updates when you click a button.
    - Displays plot links in a responsive grid; `grid_cols` controls the default column count.
    - Grouping: if recursive=True, groups by subdirectory; otherwise all plots are shown in a single grid.

    Parameters:
    dir (str | Path): Directory to scan for .html files.
    file (str | Path, optional): Output HTML file path (absolute or relative to `dir`).
    recursive (bool, optional): Whether to search subdirectories.
    exclude (Iterable[str], optional): Filenames to exclude (case insensitive).
    sort (str): Sort by "title", "name", or "mtime" (modification time).
    preview (bool, optional): Whether to include an iframe preview panel.
    grid_cols (int, optional): Number of grid columns (default 3).
    image_types (list[str] | None, optional): List of image file extensions to include (e.g. ['.png','.jpg','.gif']). If None, only .html files are included.
    preview_height_px (int, optional): Height of the preview iframe in pixels.
    icon (str, optional): Name of the SVG icon file (without .svg) to use as favicon.
    """
    dir = Path(dir).expanduser().resolve()
    out_path = (dir / file) if not Path(file).is_absolute() else Path(file)

    with pkg_resources.path(icon_pkg, f"{icon}.svg") as svg_path:
        shutil.copy(svg_path, Path(dir) / f"{icon}.svg")

    exts = {".html"}
    if image_types:
        exts |= {ext.lower() if ext.startswith('.') else f'.{ext.lower()}' for ext in image_types}

    pattern = "**/*" if recursive else "*"
    files = [p for p in dir.glob(pattern) if p.is_file() and p.suffix.lower() in exts]
    
    excl = {e.lower() for e in exclude}
    files = [p for p in files if p.name.lower() not in excl]

    items = []
    for p in files:
        if p.suffix.lower() == ".html":
            title = _extract_html_title(p) or p.stem
        else:  # image file
            title = p.stem

        # grouping key: subdirectory if recursive; unused otherwise
        if recursive:
            try:
                rel_parent = p.parent.relative_to(dir)
                group = rel_parent.as_posix() if rel_parent.as_posix() != "." else "root"
            except Exception:
                group = "root"
        else:
            group = None

        rel = p.relative_to(out_path.parent).as_posix()
        items.append(
            {
                "path": p,
                "rel": rel,
                "title": title,
                "name": p.name,
                "mtime": p.stat().st_mtime,
                "group": group,
            }
        )

    if sort == "mtime":
        items.sort(key=lambda d: d["mtime"], reverse=True)
    elif sort == "name":
        items.sort(key=lambda d: natural_key(d["name"]))
    else:  # title (default)
        items.sort(key=lambda d: natural_key(d["title"]))

    if not items:
        raise FileNotFoundError(f"No .html files found in: {dir}")

    grid_cols = max(1, int(grid_cols))

    # Build HTML
    def _card_html(it):
        title_esc = html.escape(it["title"])
        rel_esc = html.escape(it["rel"])
        return f"""
        <div class="card">
        <div class="title">{title_esc}</div>
        <div class="actions">
            <button class="btn" data-href="{rel_esc}" onclick="openPlot(this)">Open</button>
            <a class="link" href="{rel_esc}" target="_blank" rel="noopener">New tab</a>
        </div>
        <div class="meta">{html.escape(it["name"])}</div>
        </div>
        """.strip()

    # groups_html entries include their own `.grid` blocks (recursive) or a single `.grid` (non-recursive)
    groups_html = []

    if recursive:
        # group items by subdirectory
        items_by_group = sorted(items, key=lambda d: (natural_key(d["group"]), natural_key(d["title"])))
        for group_name, grp in groupby(items_by_group, key=lambda d: d["group"]):
            cards = "\n".join(_card_html(it) for it in grp)
            groups_html.append(
                f"""
                <div class="group">
                <div class="group-title">{html.escape(group_name)}</div>
                <div class="grid">
                    {cards}
                </div>
                </div>
                """.strip()
            )
    else:
        cards = "\n".join(_card_html(it) for it in items)
        groups_html.append(f'<div class="grid">{cards}</div>')

    first_rel = html.escape(items[0]["rel"])

    preview_html = ""
    js_preview = ""
    if preview:
        preview_html = f"""
        <div class="preview">
          <div class="previewbar">
            <div id="currentTitle" class="current">Showing: {html.escape(items[0]["title"])}</div>
            <div class="spacer"></div>
            <a id="currentLink" class="link" href="{first_rel}" target="_blank" rel="noopener">Open in new tab</a>
          </div>
          <iframe id="viewer" src="{first_rel}" style="height:{preview_height_px}px;"></iframe>
        </div>
        """.strip()

        js_preview = """
        function openPlot(btn) {
          const href = btn.getAttribute('data-href');
          const title = btn.closest('.card').querySelector('.title').textContent;
          const iframe = document.getElementById('viewer');
          const currentTitle = document.getElementById('currentTitle');
          const currentLink = document.getElementById('currentLink');

          iframe.src = href;
          currentTitle.textContent = 'Showing: ' + title;
          currentLink.href = href;

          // highlight selected
          document.querySelectorAll('.card').forEach(c => c.classList.remove('active'));
          btn.closest('.card').classList.add('active');
        }

        // default highlight first card
        window.addEventListener('DOMContentLoaded', () => {
          const first = document.querySelector('.card');
          if (first) first.classList.add('active');
        });
        """.strip()

    page = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width,initial-scale=1" />
  <title>EDMS {re_un_cap('.'.join(file.split('.')[:-1]))}</title>
  <link rel="icon" type="image/svg+xml" href="{dir}/{icon}.svg">
  <style>
    :root {{
      --bg: #0b0c10;
      --panel: #111318;
      --card: #151823;
      --text: #e8e8ea;
      --muted: #a8abb5;
      --accent: #6ea8fe;
      --border: #2a2f3a;
    }}
    body {{
      margin: 0; padding: 20px;
      background: var(--bg); color: var(--text);
      font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Arial, "Helvetica Neue", Helvetica;
    }}
    h1 {{ margin: 0 0 12px; font-size: 20px; }}

    .wrap {{ display: grid; grid-template-columns: 1fr; gap: 16px; }}

    .topbar {{
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 12px;
    }}

    .grid {{
      display: grid;
      grid-template-columns: repeat({grid_cols}, minmax(240px, 1fr));
      gap: 10px;
      align-items: start;
    }}

    .card {{
      background: var(--card);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 10px;
    }}
    .card.active {{ border-color: var(--accent); box-shadow: 0 0 0 2px rgba(110,168,254,0.15) inset; }}

    .group {{ margin-top: 10px; }}
    .group:first-child {{ margin-top: 0; }}

    .group-title {{
    font-size: 12px;
    color: var(--muted);
    margin: 6px 4px 10px;
    letter-spacing: 0.06em;
    text-transform: uppercase;
    }}

    .title {{ font-size: 13px; font-weight: 650; line-height: 1.2; word-break: break-word; }}
    .meta {{ margin-top: 6px; font-size: 11px; color: var(--muted); word-break: break-word; }}
    .actions {{ margin-top: 8px; display: flex; gap: 10px; align-items: center; }}

    .btn {{
      background: var(--accent);
      color: #07101f;
      border: 0;
      border-radius: 10px;
      padding: 7px 10px;
      font-weight: 650;
      cursor: pointer;
    }}
    .btn:hover {{ filter: brightness(1.05); }}
    .link {{ color: var(--muted); text-decoration: none; font-size: 12px; }}
    .link:hover {{ color: var(--text); text-decoration: underline; }}

    .preview {{
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 12px;
      overflow: hidden;
    }}
    .previewbar {{
      display: flex; gap: 10px; align-items: center;
      padding: 10px 12px;
      border-bottom: 1px solid var(--border);
    }}
    .current {{ font-size: 13px; color: var(--text); }}
    .spacer {{ flex: 1; }}
    iframe {{
      width: 100%;
      border: 0;
      background: #fff;
    }}

    @media (max-width: 1100px) {{
      .grid {{ grid-template-columns: 1fr; }}
    }}
  </style>
</head>
<body>
  <h1>EDMS {re_un_cap('.'.join(file.split('.')[:-1]))}</h1>
  <div class="wrap">
    <div class="topbar">
      {"".join(groups_html)}
    </div>
    {preview_html}
  </div>

  <script>
    {js_preview}
  </script>
</body>
</html>
"""

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(page, encoding="utf-8")
    return out_path