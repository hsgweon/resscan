#!/usr/bin/env python3

# resscan/report_common.py
"""
Pieces shared by the HTML reports (index, HomScan and VarScan pages): the colour
theme, the page header, sortable tables and the "download as TSV" helper.

Every report is a single self-contained file. Nothing is fetched from the
network, so the pages open offline and can be archived with the results.
"""
import html
import json
import re


def sanitize_filename(name):
    """Converts a family name into the file name used for its report page."""
    name = name.replace(' ', '_').replace(';', '_').replace('/', '_')
    name = re.sub(r'[^\w\-.]', '', name)
    return name.lower() + ".html"


# Colours are CSS variables so that one switch changes the whole page. The dark
# theme is the default on screen; printing always uses the light theme.
_LIGHT_VARS = """
            --bg: #ffffff; --panel: #f6f7f9; --panel-border: #d5d9e0; --plot-bg: #ffffff; --plot-border: #c9ced6;
            --text: #1c2026; --text-strong: #000000; --text-muted: #5b6470; --rule: #d5d9e0;
            --ref: #0a7a2f; --que: #8a5a00; --marker: #c81e2c; --mismatch-bg: #c81e2c; --mismatch-fg: #ffffff;
            --pid-full: #1f9d3a; --pid-partial: #d99a00; --unique: #7a5a00; --depth: #4c78a8; --axis: #7b8490;
            --link: #1557b0; --row-alt: #eef0f3; --highlight: #dbe9ff; --shadow: rgba(0,0,0,0.08);
"""

THEME_CSS = """
        :root {
            --bg: #121212; --panel: #1e1e1e; --panel-border: #333333; --plot-bg: #000000; --plot-border: #222222;
            --text: #e0e0e0; --text-strong: #ffffff; --text-muted: #9a9a9a; --rule: #333333;
            --ref: #00ff41; --que: #ffcc00; --marker: #ff4757; --mismatch-bg: #ff4757; --mismatch-fg: #000000;
            --pid-full: #00ff41; --pid-partial: #ffcc00; --unique: #ffd700; --depth: #6ea8dc; --axis: #666666;
            --link: #7ab8ff; --row-alt: #232323; --highlight: #2c3e50; --shadow: rgba(0,0,0,0.3);
        }
        :root[data-theme="light"] {""" + _LIGHT_VARS + """        }
        @media print {
            :root {""" + _LIGHT_VARS + """            }
            .no-print { display: none !important; }
            .seq-container, .chart-container { overflow: visible !important; }
        }
        *, *::before, *::after { box-sizing: border-box; }
        body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 2em; background-color: var(--bg); color: var(--text); }
        .container { max-width: 1200px; margin: auto; }
        a { color: var(--link); }
        h1 { color: var(--text-strong); border-bottom: 2px solid var(--rule); padding-bottom: 0.5em; margin-top: 0.4em; }
        h2 { color: var(--text-strong); margin-top: 2em; border-bottom: 1px solid var(--rule); padding-bottom: 0.3em; }
        .page-bar { display: flex; justify-content: space-between; align-items: center; gap: 1em; flex-wrap: wrap; font-size: 0.9em; color: var(--text-muted); }
        .page-bar .sample { color: var(--text-strong); font-weight: 600; }
        .toolbar { display: flex; gap: 0.6em; flex-wrap: wrap; }
        button.tool { background-color: var(--panel); color: var(--text); border: 1px solid var(--panel-border); border-radius: 6px; padding: 0.35em 0.8em; font-size: 0.9em; cursor: pointer; }
        button.tool:hover { border-color: var(--text-muted); }
        .legend { background-color: var(--panel); padding: 1em; border-radius: 8px; margin: 1.5em 0; display: flex; gap: 2em; align-items: center; flex-wrap: wrap; border: 1px solid var(--panel-border); }
        .legend h3 { margin: 0; padding-right: 1em; border-right: 1px solid var(--rule); color: var(--text-strong); font-size: 1em; }
        .legend-item { display: flex; align-items: center; gap: 0.5em; font-size: 0.9em; }
        .note { color: var(--text-muted); font-size: 0.9em; line-height: 1.5; margin: 0.6em 0; }
        .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 1em; margin: 1em 0; }
        .stat-card { background-color: var(--panel); border-radius: 8px; padding: 1em 1.2em; box-shadow: 0 4px 6px var(--shadow); border: 1px solid var(--panel-border); }
        .stat-card h3 { margin: 0 0 0.4em 0; color: var(--text-muted); font-size: 0.9em; font-weight: 600; }
        .stat-card p { margin: 0; color: var(--text-strong); font-size: 1.5em; font-weight: 700; }
        .stat-card small { display: block; margin-top: 0.3em; color: var(--text-muted); font-size: 0.8em; font-weight: 400; }
        table.data { width: 100%; border-collapse: collapse; margin: 1em 0; font-size: 0.9em; }
        table.data th, table.data td { border: 1px solid var(--panel-border); padding: 6px 10px; text-align: left; vertical-align: top; }
        table.data th { background-color: var(--panel); color: var(--text-strong); font-weight: 600; }
        .table-wrap { overflow-x: auto; }
        table.data td.num, table.data th.num { text-align: right; font-variant-numeric: tabular-nums; }
        table.data tbody tr:nth-child(even) { background-color: var(--row-alt); }
        table.sortable th { cursor: pointer; user-select: none; }
        table.sortable th[data-dir="asc"]::after { content: " \\25B2"; font-size: 0.7em; }
        table.sortable th[data-dir="desc"]::after { content: " \\25BC"; font-size: 0.7em; }
        input.filter { background-color: var(--panel); color: var(--text); border: 1px solid var(--panel-border); border-radius: 6px; padding: 0.4em 0.7em; font-size: 0.9em; min-width: 260px; }
"""

THEME_JS = """
        (function () {
            // Theme: remembered between pages of the same report.
            try { const saved = localStorage.getItem('resscan-theme'); if (saved) document.documentElement.setAttribute('data-theme', saved); } catch (e) {}
        })();
        function toggleTheme() {
            const root = document.documentElement;
            const next = root.getAttribute('data-theme') === 'light' ? 'dark' : 'light';
            root.setAttribute('data-theme', next);
            try { localStorage.setItem('resscan-theme', next); } catch (e) {}
        }
        function downloadTSV(filename, rows) {
            // rows: array of arrays; the first row is the header.
            const clean = v => String(v === null || v === undefined ? '' : v).replace(/[\\t\\n\\r]+/g, ' ');
            const text = rows.map(r => r.map(clean).join('\\t')).join('\\n') + '\\n';
            const link = document.createElement('a');
            link.href = URL.createObjectURL(new Blob([text], { type: 'text/tab-separated-values' }));
            link.download = filename;
            document.body.appendChild(link); link.click(); link.remove();
        }
        function makeSortable(table) {
            table.querySelectorAll('thead th').forEach((th, col) => {
                th.addEventListener('click', () => {
                    const dir = th.getAttribute('data-dir') === 'desc' ? 'asc' : 'desc';
                    table.querySelectorAll('thead th').forEach(o => o.removeAttribute('data-dir'));
                    th.setAttribute('data-dir', dir);
                    const body = table.tBodies[0];
                    const key = tr => { const c = tr.cells[col]; return c.hasAttribute('data-sort') ? c.getAttribute('data-sort') : c.textContent.trim(); };
                    const rows = Array.from(body.rows);
                    const numeric = rows.every(r => key(r) === '' || !isNaN(parseFloat(key(r))));
                    rows.sort((a, b) => {
                        const x = key(a), y = key(b);
                        const cmp = numeric ? (parseFloat(x) || 0) - (parseFloat(y) || 0) : x.localeCompare(y);
                        return dir === 'asc' ? cmp : -cmp;
                    });
                    rows.forEach(r => body.appendChild(r));
                });
            });
        }
        function attachFilter(input, table) {
            input.addEventListener('input', () => {
                const needle = input.value.toLowerCase();
                Array.from(table.tBodies[0].rows).forEach(r => { r.style.display = r.textContent.toLowerCase().includes(needle) ? '' : 'none'; });
            });
        }
        document.addEventListener('DOMContentLoaded', () => { document.querySelectorAll('table.sortable').forEach(makeSortable); });
"""


def page_bar(sample_name=None, index_href=None, kind=None, extra_buttons=""):
    """The strip at the top of every page: sample, link back to the index, and tools."""
    left = []
    if index_href:
        left.append(f'<a href="{html.escape(index_href, quote=True)}">&larr; Sample overview</a>')
    if sample_name:
        left.append(f'Sample <span class="sample">{html.escape(sample_name)}</span>')
    if kind:
        left.append(html.escape(kind))
    return f"""
        <div class="page-bar">
            <div>{' &nbsp;|&nbsp; '.join(left)}</div>
            <div class="toolbar no-print">
                {extra_buttons}
                <button class="tool" onclick="toggleTheme()">Light / dark</button>
                <button class="tool" onclick="window.print()">Print</button>
            </div>
        </div>"""


def title_text(page_title, sample_name=None):
    return f"{sample_name}: {page_title}" if sample_name else page_title


def js_data(value):
    """JSON for embedding inside a <script> element."""
    return json.dumps(value).replace("</", "<\\/")
