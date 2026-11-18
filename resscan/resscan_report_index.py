#!/usr/bin/env python3

# resscan/resscan_report_index.py
"""
Writes <prefix>_index.html: one overview page per sample that lists every HomScan
and VarScan result, links to the per-family pages, and records how the results
were produced. It only reads the final tables already written by the pipeline.
"""
import argparse
import csv
import html
import os
import sys
try:
    from resscan.utils import BColors
    from resscan.report_common import THEME_CSS, THEME_JS, page_bar, title_text, js_data, sanitize_filename
except ImportError:
    from utils import BColors
    from report_common import THEME_CSS, THEME_JS, page_bar, title_text, js_data, sanitize_filename


def read_table(path):
    """Rows of a ResScan TSV as dictionaries; an absent or empty file gives no rows."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames:
            return []
        reader.fieldnames = [name.lstrip('#') for name in reader.fieldnames]
        return list(reader)

def read_key_values(path):
    pairs = []
    if os.path.exists(path):
        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#') or '\t' not in line:
                    continue
                key, value = line.rstrip('\n').split('\t', 1)
                pairs.append((key, value))
    return pairs

def number(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default

def fmt(value, decimals):
    try:
        return f"{float(value):,.{decimals}f}"
    except (TypeError, ValueError):
        return html.escape(str(value)) if value not in (None, '') else ''

def group_members(row):
    """References accounted for by a HomScan row, most likely first."""
    members = []
    for item in (row.get('Allocation_Proportions') or '').split(','):
        if ':' in item:
            reference, share = item.rsplit(':', 1)
            members.append((reference.strip(), share.strip()))
    return members

def short_name(reference_id):
    parts = reference_id.split('__')
    return parts[1] if len(parts) > 1 else reference_id

def load_reference_names(metadata_path):
    """{'ARO_3002601': (gene name, gene family)} from the database metadata, for naming members of group rows."""
    names = {}
    if not metadata_path or not os.path.exists(metadata_path):
        return names
    with open(metadata_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            key = row.get('Sequence_ID', '').split('__')[0]
            if key and key not in names:
                names[key] = (row.get('ARO_Name') or row.get('Card_Short_Name') or key, row.get('AMR_Gene_Family', ''))
    return names

def tidy_classes(text):
    """Drug classes of a row, each listed once."""
    seen = []
    for item in text.replace('|', ';').split(';'):
        item = item.strip()
        if item and item not in seen:
            seen.append(item)
    return "; ".join(seen)

def page_link(family, html_dir, dir_name, label=None):
    """Link to a family page if that page exists."""
    file_name = sanitize_filename(family)
    if html_dir and os.path.exists(os.path.join(html_dir, file_name)):
        return f'<a href="{html.escape(dir_name + "/" + file_name, quote=True)}">{html.escape(label or family)}</a>'
    return html.escape(label or family)

def build_index(sample_name, homscan_rows, varscan_rows, uscg_rows, provenance, totals,
                homscan_html_dir, varscan_html_dir, homscan_dir_name, varscan_dir_name, reference_names=None):
    reference_names = reference_names or {}
    hom_body, hom_tsv = [], [["gene_family", "gene", "references_in_row", "drug_class", "resistance_mechanism", "reads",
                              "fragments", "lateral_coverage_%", "copies_per_cell", "reads_per_kb_per_Gb"]]
    for row in sorted(homscan_rows, key=lambda r: -number(r.get('RPKPC'))):
        members = group_members(row)
        family = row.get('AMR_Gene_Family', '')
        family_cell = page_link(family, homscan_html_dir, homscan_dir_name)
        if len(members) > 1:
            top_id = (row.get('Top_ARO') or members[0][0]).split('__')[0]
            top, top_family = reference_names.get(top_id, (short_name(top_id), ''))
            member_names = row.get('ARO_Name') or ", ".join(m for m, _ in members)
            gene = (f'{html.escape(top)} <span class="tag" title="{html.escape(member_names, quote=True)}">'
                    f'group of {len(members)}</span>')
            gene_text = f"{top} (group of {len(members)})"
            if family == 'multiple' and top_family:
                # The group spans several families: link to the page of its most likely member.
                family_cell = 'several: ' + page_link(top_family, homscan_html_dir, homscan_dir_name)
                family = f"several ({top_family} most likely)"
        else:
            gene_text = row.get('ARO_Name') or short_name(row.get('ARO', ''))
            gene = html.escape(gene_text)
        hom_body.append(
            f'<tr><td>{family_cell}</td><td>{gene}</td>'
            f'<td>{html.escape(tidy_classes(row.get("Drug_Class", "")))}</td>'
            f'<td class="num">{fmt(row.get("Read_Count"), 0)}</td><td class="num">{fmt(row.get("Fragment_Count"), 0)}</td>'
            f'<td class="num">{fmt(row.get("Lateral_Coverage_%"), 1)}</td>'
            f'<td class="num" data-sort="{number(row.get("RPKPC"))}">{fmt(row.get("RPKPC"), 4)}</td>'
            f'<td class="num" data-sort="{number(row.get("RPKG"))}">{fmt(row.get("RPKG"), 2)}</td></tr>')
        hom_tsv.append([family, gene_text, len(members) or 1, tidy_classes(row.get('Drug_Class', '')), row.get('Resistance_Mechanism', ''),
                        row.get('Read_Count', ''), row.get('Fragment_Count', ''), row.get('Lateral_Coverage_%', ''),
                        row.get('RPKPC', ''), row.get('RPKG', '')])

    var_body, var_tsv = [], [["gene_family", "entry", "mutation", "drug_class", "reads", "fragments", "copies_per_cell",
                              "reads_per_kb_per_Gb"]]
    for row in sorted(varscan_rows, key=lambda r: -number(r.get('Read_Count'))):
        family = row.get('AMR_Gene_Family', '')
        reference = row.get('ARO', '')
        mutation = reference.split('__')[-1].replace('_', ', ') if reference.count('__') >= 3 else ''
        var_body.append(
            f'<tr><td>{page_link(family, varscan_html_dir, varscan_dir_name)}</td>'
            f'<td>{html.escape(row.get("ARO_Name") or short_name(reference))}</td><td>{html.escape(mutation)}</td>'
            f'<td>{html.escape(row.get("Drug_Class", "").replace(";", "; "))}</td>'
            f'<td class="num">{fmt(row.get("Read_Count"), 0)}</td><td class="num">{fmt(row.get("Fragment_Count"), 0)}</td>'
            f'<td class="num" data-sort="{number(row.get("RPKPC"))}">{fmt(row.get("RPKPC"), 4)}</td>'
            f'<td class="num" data-sort="{number(row.get("RPKG"))}">{fmt(row.get("RPKG"), 2)}</td></tr>')
        var_tsv.append([family, row.get('ARO_Name') or short_name(reference), mutation, row.get('Drug_Class', ''),
                        row.get('Read_Count', ''), row.get('Fragment_Count', ''), row.get('RPKPC', ''), row.get('RPKG', '')])

    total_copies = sum(number(r.get('RPKPC')) for r in homscan_rows)
    marker_reads = sum(number(r.get('Read_Count')) for r in uscg_rows)
    markers_detected = sum(1 for r in uscg_rows if number(r.get('Read_Count')) > 0)
    cards = []
    if totals.get('reads'):
        cards.append(f'<div class="stat-card"><h3>Reads analysed</h3><p>{totals["reads"]:,.0f}</p><small>{totals.get("bases", 0) / 1e9:,.2f} Gb</small></div>')
    cards.append(f'<div class="stat-card"><h3>Acquired resistance genes</h3><p>{len(homscan_rows):,}</p><small>rows in the HomScan table</small></div>')
    cards.append(f'<div class="stat-card"><h3>Total copies per cell</h3><p>{total_copies:,.3f}</p><small>sum over the HomScan table</small></div>')
    cards.append(f'<div class="stat-card"><h3>Resistance mutations</h3><p>{len(varscan_rows):,}</p><small>{sum(number(r.get("Read_Count")) for r in varscan_rows):,.0f} supporting reads</small></div>')
    if uscg_rows:
        cards.append(f'<div class="stat-card"><h3>Marker genes detected</h3><p>{markers_detected} of {len(uscg_rows)}</p><small>{marker_reads:,.0f} reads</small></div>')

    marker_rows = "".join(
        f'<tr><td>{html.escape(r.get("USCG_ID", ""))}</td><td class="num">{fmt(r.get("Gene_Length_AA"), 0)}</td>'
        f'<td class="num">{fmt(r.get("Read_Count"), 0)}</td><td class="num">{fmt(r.get("Fragment_Count"), 0)}</td>'
        f'<td class="num">{fmt(r.get("RPK"), 2)}</td></tr>' for r in uscg_rows)
    provenance_rows = "".join(f'<tr><td>{html.escape(k.replace("_", " "))}</td><td>{html.escape(v)}</td></tr>' for k, v in provenance)
    version = dict(provenance).get('ResScan_Version', '')
    card_version = dict(provenance).get('Database_CARD_Version', '')
    subtitle = " &nbsp;|&nbsp; ".join(x for x in [f"ResScan {html.escape(version)}" if version else "",
                                                 f"CARD {html.escape(card_version)}" if card_version else "",
                                                 html.escape(dict(provenance).get('Run_Date', ''))] if x)
    safe = sanitize_filename(sample_name).replace('.html', '')
    buttons = (f'<button class="tool" onclick="downloadTSV({html.escape(js_data(safe + "_homscan_overview.tsv"), quote=True)}, HOMSCAN_TABLE)">HomScan TSV</button>'
               f'<button class="tool" onclick="downloadTSV({html.escape(js_data(safe + "_varscan_overview.tsv"), quote=True)}, VARSCAN_TABLE)">VarScan TSV</button>')

    return f"""<!DOCTYPE html>
<html lang="en-GB">
<head>
    <meta charset="UTF-8">
    <title>{html.escape(title_text('ResScan report', sample_name))}</title>
    <style>{THEME_CSS}
        .tag {{ display: inline-block; font-size: 0.8em; color: var(--text-muted); border: 1px solid var(--panel-border); border-radius: 10px; padding: 0 0.6em; margin-left: 0.4em; cursor: help; white-space: nowrap; }}
        .section-head {{ display: flex; justify-content: space-between; align-items: flex-end; gap: 1em; flex-wrap: wrap; }}
        details {{ margin: 1.5em 0; }}
        details summary {{ cursor: pointer; color: var(--text-strong); font-weight: 600; font-size: 1.1em; }}
        @media print {{ details {{ display: none; }} }}
    </style>
    <script>{THEME_JS}</script>
</head>
<body>
    <div class="container">
        {page_bar(sample_name, None, "Overview", buttons)}
        <h1>ResScan report: {html.escape(sample_name)}</h1>
        <p class="note">{subtitle}</p>
        <div class="stats-grid">{''.join(cards)}</div>

        <div class="section-head"><h2 style="flex: 1;">HomScan: acquired resistance genes</h2>
            <input class="filter no-print" id="hom-filter" type="search" placeholder="Filter by family, gene or drug class"></div>
        <p class="note">One row per gene, or per group of genes the reads cannot tell apart (hover over &ldquo;group of n&rdquo; for the members;
        the gene named is the most likely one). Copies per cell is the gene's read density divided by the mean density of the
        single-copy marker genes. Click a family for the read-level view; click a column heading to sort.</p>
        {'<table class="data sortable" id="hom-table"><thead><tr><th>Gene family</th><th>Gene</th><th>Drug class</th><th class="num">Reads</th><th class="num">Fragments</th><th class="num">Lateral coverage %</th><th class="num">Copies per cell</th><th class="num">Reads per kb per Gb</th></tr></thead><tbody>' + ''.join(hom_body) + '</tbody></table>' if hom_body else '<p class="note">No acquired resistance genes were detected.</p>'}

        <div class="section-head"><h2 style="flex: 1;">VarScan: resistance-conferring mutations</h2>
            <input class="filter no-print" id="var-filter" type="search" placeholder="Filter by family, entry or drug class"></div>
        <p class="note">One row per database entry with at least one read carrying every mutation of that entry.
        Click a family to see the supporting reads aligned to the reference.</p>
        {'<table class="data sortable" id="var-table"><thead><tr><th>Gene family</th><th>Entry</th><th>Mutation(s)</th><th>Drug class</th><th class="num">Reads</th><th class="num">Fragments</th><th class="num">Copies per cell</th><th class="num">Reads per kb per Gb</th></tr></thead><tbody>' + ''.join(var_body) + '</tbody></table>' if var_body else '<p class="note">No resistance-conferring mutations were confirmed.</p>'}

        <details><summary>Single-copy marker genes ({markers_detected} of {len(uscg_rows)} detected)</summary>
            <p class="note">The per-cell denominator is the mean read density over these genes. Even counts across the genes indicate a reliable estimate.</p>
            <table class="data sortable"><thead><tr><th>Marker gene</th><th class="num">Length (aa)</th><th class="num">Reads</th><th class="num">Fragments</th><th class="num">Reads per kb</th></tr></thead><tbody>{marker_rows}</tbody></table>
        </details>
        <details><summary>How these results were produced</summary>
            <table class="data"><tbody>{provenance_rows}</tbody></table>
        </details>
    </div>
    <script>
        const HOMSCAN_TABLE = {js_data(hom_tsv)};
        const VARSCAN_TABLE = {js_data(var_tsv)};
        document.addEventListener('DOMContentLoaded', () => {{
            [['hom-filter', 'hom-table'], ['var-filter', 'var-table']].forEach(([i, t]) => {{
                const input = document.getElementById(i), table = document.getElementById(t);
                if (input && table) attachFilter(input, table); else if (input) input.style.display = 'none';
            }});
        }});
    </script>
</body>
</html>
"""

def main():
    parser = argparse.ArgumentParser(description="Writes the per-sample overview page that links the HomScan and VarScan reports.")
    parser.add_argument("--output-prefix", required=True, help="The run's output prefix (the results directory), as given to resscan -o.")
    parser.add_argument("--sample-name", default=None, help="Sample name for the heading. Default: the prefix name.")
    parser.add_argument("--metadata", default=None, help="Database metadata TSV, used to name the most likely member of grouped rows.")
    args = parser.parse_args()

    results_dir = os.path.abspath(args.output_prefix)
    prefix = os.path.basename(results_dir.rstrip('/'))
    sample_name = args.sample_name or prefix
    path = lambda suffix: os.path.join(results_dir, f"{prefix}_{suffix}")

    totals = {}
    total_bases_file = os.path.join(results_dir, "tmp", "scgscan", f"{prefix}_total_bases.txt")
    if os.path.exists(total_bases_file):
        values = open(total_bases_file).read().split()
        if len(values) >= 2:
            totals = {'bases': number(values[0]), 'reads': number(values[1])}

    homscan_dir_name, varscan_dir_name = f"{prefix}_homscan_html", f"{prefix}_varscan_html"
    page = build_index(
        sample_name, read_table(path("homscan.tsv")), read_table(path("varscan.tsv")), read_table(path("uscg.tsv")),
        read_key_values(path("provenance.txt")), totals,
        os.path.join(results_dir, homscan_dir_name), os.path.join(results_dir, varscan_dir_name),
        homscan_dir_name, varscan_dir_name, load_reference_names(args.metadata))
    output_path = path("index.html")
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(page)
    print(BColors.green(f"--- Wrote sample overview: {output_path} ---"))

if __name__ == "__main__":
    main()
