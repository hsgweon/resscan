#!/usr/bin/env python3

# resscan/homscan_visualise.py
import argparse
import os
import sys
import csv
import html
import math
import json
import re
from collections import defaultdict
try:
    from resscan.utils import BColors
    from resscan.report_common import THEME_CSS, THEME_JS, page_bar, title_text, js_data, sanitize_filename
except ImportError:
    from utils import BColors
    from report_common import THEME_CSS, THEME_JS, page_bar, title_text, js_data, sanitize_filename

# Reads drawn individually per reference. Above this the pile is an evenly spaced
# sample; the depth track and every statistic are still computed from all reads.
DEFAULT_MAX_READS_DRAWN = 500
# A family with hundreds of near-identical references (KPC, TEM) would still draw
# hundreds of thousands of reads, so the per-reference limit is lowered until the
# page as a whole stays within this budget (never below MIN_READS_DRAWN each).
PAGE_READ_BUDGET = 20000
MIN_READS_DRAWN = 50
# References listed when a read is clicked. A read from a large allele family
# matches every member, so the list is limited to its best matches to keep pages small.
MAX_HITS_LISTED_PER_READ = 20

def reads_drawn_limit(max_reads_drawn, reference_count):
    """Per-reference limit on drawn reads for a page showing reference_count references."""
    if max_reads_drawn <= 0:
        return 0
    return min(max_reads_drawn, max(MIN_READS_DRAWN, PAGE_READ_BUDGET // max(1, reference_count)))

def compact_read_hits(query_ids, hits_by_query_id):
    """
    The data behind the click-on-a-read table, in compact form:
    references = [[reference, family], ...] and, per read,
    [number of references matched, [reference index, nucleotide identity, protein identity,
    uniqueness, position, length], ...] holding its MAX_HITS_LISTED_PER_READ best matches.
    """
    references, index_of, reads = [], {}, {}
    for query_id in sorted(query_ids):
        hits = hits_by_query_id.get(query_id, [])
        best = sorted(hits, key=lambda h: (-float(h['protein_pid']), -float(h['nucleotide_pid']), h['ARO_matched']))
        rows = []
        for hit in best[:MAX_HITS_LISTED_PER_READ]:
            key = hit['ARO_matched']
            if key not in index_of:
                index_of[key] = len(references)
                references.append([key, hit['AMR_Gene_Family']])
            rows.append([index_of[key], hit['nucleotide_pid'], hit['protein_pid'], hit['uniqueness'],
                         hit['position_on_ref'], hit['nucleotide_denominator']])
        reads[query_id] = [len(hits)] + rows
    return {'references': references, 'reads': reads}

def load_metadata(metadata_path):
    """Loads metadata to map AROs to AMR_Gene_Families."""
    print(BColors.cyan(f"--- Loading metadata from: {metadata_path} ---"))
    if not os.path.exists(metadata_path):
        print(BColors.red(f"Error: Metadata file not found at '{metadata_path}'"), file=sys.stderr)
        sys.exit(1)
    
    aro_to_family = {}
    family_to_aros = defaultdict(list)
    try:
        with open(metadata_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                aro_id = row['Sequence_ID'] 
                family = row['AMR_Gene_Family'] 
                aro_to_family[aro_id] = family
                # Treat the entire family string as a single unique key for consistency
                family_to_aros[family].append(aro_id)
    except (Exception, KeyError) as e:
        print(BColors.red(f"Error reading or parsing metadata file: {e}"), file=sys.stderr)
        sys.exit(1)
    print(BColors.green(f"--- Loaded metadata for {len(family_to_aros)} unique family strings and {len(aro_to_family)} AROs."))
    return aro_to_family, family_to_aros

def load_reference_sequences(db_path, aro_ids):
    """Loads sequences for a specific list of AROs."""
    print(BColors.cyan(f"--- Loading reference sequences for {len(aro_ids)} AROs ---"))
    aro_sequences = {}
    aro_ids_set = set(aro_ids)
    
    with open(db_path, 'r') as f:
        current_id = None
        seq_parts = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id in aro_ids_set:
                    aro_sequences[current_id] = "".join(seq_parts)
                current_id = line[1:] 
                seq_parts = []
            elif current_id in aro_ids_set:
                seq_parts.append(line)
    
    if current_id in aro_ids_set:
        aro_sequences[current_id] = "".join(seq_parts)
        
    print(BColors.green(f"--- Found sequences for {len(aro_sequences)} of the requested AROs."))
    return aro_sequences

def load_all_hits(input_files, pid_cutoff, pid_type, aro_to_family):
    """
    Loads all hits from multiple files, filters by the chosen PID type, and enriches with family info.
    """
    pid_column = f"{pid_type}_pid"
    print(BColors.cyan(f"--- Loading all hits with {pid_column} >= {pid_cutoff*100:.2f}% ---"))
    all_hits_raw = []
    pid_cutoff_percent = pid_cutoff * 100.0
    for file_path in input_files:
        try:
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    try:
                        if float(row[pid_column]) >= pid_cutoff_percent:
                            aro_val = row.get('ARO_matched') or row.get('ARO')
                            row['AMR_Gene_Family'] = aro_to_family.get(aro_val, 'Unknown_Family')
                            all_hits_raw.append(row)
                    except (ValueError, KeyError):
                        continue
        except FileNotFoundError:
            print(BColors.yellow(f"Warning: Input file not found, skipping: {file_path}"), file=sys.stderr)
    
    hits_by_query_id = defaultdict(list)
    for hit in all_hits_raw:
        aro_val = hit.get('ARO_matched') or hit.get('ARO')
        simplified_hit = {
            'ARO_matched': aro_val,
            'AMR_Gene_Family': hit['AMR_Gene_Family'],
            'nucleotide_pid': hit['nucleotide_pid'],
            'protein_pid': hit['protein_pid'],
            'uniqueness': hit.get('uniqueness', 'UNKNOWN'),
            'position_on_ref': hit['position_on_ref'],
            'nucleotide_denominator': hit['nucleotide_denominator'],
        }
        hits_by_query_id[hit['query_id']].append(simplified_hit)

    print(BColors.green(f"--- Loaded {len(all_hits_raw)} hits from {len(hits_by_query_id)} unique reads passing the filter."))
    return all_hits_raw, hits_by_query_id

class _LanePacker:
    """
    Packs read rectangles into display lanes.

    Each read is placed in the lowest-numbered lane whose previous read has
    already ended, which is what determines the stacked appearance of the plot.
    Searching that lane by scanning the list is quadratic in the number of reads,
    and a deeply covered gene can carry tens of thousands, so the lane ends are
    held in a segment tree over their minimum and the search descends it instead.
    The lane chosen is identical to the one a left-to-right scan would pick.

    Lane ends are initialised to 0, so an unused lane always reads as free and is
    taken only once every lower-numbered lane is occupied.
    """

    def __init__(self, capacity):
        size = 1
        while size < max(1, capacity):
            size *= 2
        self._size = size
        self._min_end = [0] * (2 * size)
        self._lane_count = 0

    def assign(self, start, end):
        node = 1
        while node < self._size:
            node *= 2
            if self._min_end[node] > start:      # left subtree full, go right
                node += 1
        lane = node - self._size

        if lane >= self._lane_count:
            self._lane_count = lane + 1

        node = lane + self._size
        self._min_end[node] = end
        node //= 2
        while node:
            self._min_end[node] = min(self._min_end[2 * node], self._min_end[2 * node + 1])
            node //= 2
        return lane

    def lane_count(self):
        return self._lane_count if self._lane_count else 1


def compute_depth(ref_seq_len, hits):
    """Per-base depth over the reference (index 0 unused) from every hit, not only those drawn."""
    change = [0] * (ref_seq_len + 2)
    for hit in hits:
        start = max(1, int(hit['position_on_ref']))
        end = min(ref_seq_len + 1, int(hit['position_on_ref']) + int(hit['nucleotide_denominator']))
        if end > start:
            change[start] += 1
            change[end] -= 1
    depth, running = [0] * (ref_seq_len + 1), 0
    for position in range(1, ref_seq_len + 1):
        running += change[position]
        depth[position] = running
    return depth

def select_reads_to_draw(hits, max_reads_drawn):
    """All hits if they fit, otherwise an evenly spaced sample along the reference."""
    if max_reads_drawn <= 0 or len(hits) <= max_reads_drawn:
        return list(hits)
    ordered = sorted(hits, key=lambda h: (int(h['position_on_ref']), h['query_id']))
    step = len(ordered) / max_reads_drawn
    return [ordered[int(i * step)] for i in range(max_reads_drawn)]

def depth_track_svg(depth, ref_seq_len, top, height):
    """A filled step profile of depth along the reference, with its scale."""
    max_depth = max(depth) if depth else 0
    if max_depth == 0:
        return ""
    base = top + height
    y_of = lambda d: base - (d / max_depth) * height
    points = [f"1,{base:.1f}"]
    previous = 0
    for position in range(1, ref_seq_len + 1):
        if depth[position] != previous:
            points.append(f"{position},{y_of(previous):.1f}")
            points.append(f"{position},{y_of(depth[position]):.1f}")
            previous = depth[position]
    points.append(f"{ref_seq_len + 1},{y_of(previous):.1f}")
    points.append(f"{ref_seq_len + 1},{base:.1f}")
    return (f'<polygon class="depth-area" points="{" ".join(points)}"></polygon>'
            f'<line class="depth-rule" x1="0" y1="{top}" x2="{ref_seq_len}" y2="{top}"></line>'
            f'<text class="depth-label" x="4" y="{top - 4}">Depth (all reads), maximum {max_depth:,}&times;</text>')

def reference_number(reference_id):
    """'ARO_3002601__aadA__H' -> 'ARO_3002601', the form used for members of a group row."""
    return reference_id.split('__')[0]

def load_summary_table(summary_path):
    """
    Reads the final HomScan table and returns how each reference is reported:
    'own' maps a reference to its own row (reads that point to it alone), and
    'groups' maps a reference number to the group rows it belongs to (reads that
    several references share, with Allocation_Proportions giving estimated shares).
    """
    summary = {'own': {}, 'groups': defaultdict(list)}
    if not summary_path or not os.path.exists(summary_path):
        return summary
    try:
        with open(summary_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            reader.fieldnames = [name.lstrip('#') for name in reader.fieldnames]
            for row in reader:
                members = []
                for item in (row.get('Allocation_Proportions') or '').split(','):
                    if ':' in item:
                        reference, share = item.rsplit(':', 1)
                        members.append((reference.strip(), share.strip()))
                entry = {'reads': row.get('Read_Count', ''), 'copies_per_cell': row.get('RPKPC', '')}
                if len(members) <= 1:
                    summary['own'][row.get('ARO', '')] = entry
                else:
                    for reference, share in members:
                        summary['groups'][reference_number(reference)].append(dict(entry, share=share, group_size=len(members)))
    except Exception as e:
        print(BColors.yellow(f"Warning: could not read summary table '{summary_path}': {e}"), file=sys.stderr)
    return summary

def final_rows_for(summary, aro_id):
    """(own row or None, list of group rows) for one reference."""
    return summary['own'].get(aro_id), summary['groups'].get(reference_number(aro_id), [])

def short_name(reference_id):
    parts = reference_id.split('__')
    return parts[1] if len(parts) > 1 else reference_id

def generate_single_aro_plot_html(aro_id, ref_seq_len, hits, original_read_count, pid_type,
                                  max_reads_drawn=DEFAULT_MAX_READS_DRAWN, final_rows=(None, [])):
    """
    Generates the section for one reference: statistics, a depth profile from all
    reads, and the read pile (capped at max_reads_drawn reads).
    Returns (html, statistics, reads drawn).
    """
    depth = compute_depth(ref_seq_len, hits)
    covered = sum(1 for d in depth[1:] if d > 0)
    coverage = (covered / ref_seq_len * 100.0) if ref_seq_len else 0.0
    mean_depth = (sum(depth[1:]) / ref_seq_len) if ref_seq_len else 0.0
    unique_reads = sum(1 for h in hits if h.get('uniqueness') == 'UNIQUE')

    drawn = select_reads_to_draw(hits, max_reads_drawn)
    drawn_for_layout = sorted(drawn, key=lambda x: float(x['nucleotide_pid']), reverse=True)
    rect_height, gap = 8, 2
    packer = _LanePacker(len(drawn_for_layout))
    lanes = []
    for hit in drawn_for_layout:
        start = int(hit['position_on_ref'])
        lanes.append(packer.assign(start, start + int(hit['nucleotide_denominator'])))
    max_lanes = packer.lane_count()

    depth_top, depth_height = 18, 60
    pile_top = depth_top + depth_height + 26
    pid_column_for_color = f"{pid_type}_pid"
    svg_elements = [depth_track_svg(depth, ref_seq_len, depth_top, depth_height)]
    pile_label = (f"Reads: {len(drawn):,} of {original_read_count:,} shown, evenly spaced along the reference"
                  if len(drawn) < original_read_count else f"Reads: all {original_read_count:,} shown")
    svg_elements.append(f'<text class="depth-label" x="4" y="{pile_top - 8}">{pile_label}</text>')
    for hit, lane in zip(drawn_for_layout, lanes):
        y_pos = pile_top + (max_lanes - lane - 1) * (rect_height + gap)
        classes = "r " + ("full" if float(hit[pid_column_for_color]) >= 100.0 else "part")
        if hit.get('uniqueness', 'UNKNOWN') == 'UNIQUE':
            classes += " uniq"
        svg_elements.append(
            f'<rect class="{classes}" x="{hit["position_on_ref"]}" y="{y_pos}" '
            f'width="{hit["nucleotide_denominator"]}" height="{rect_height}" rx="1" '
            f'data-q="{html.escape(hit["query_id"], quote=True)}" data-n="{hit["nucleotide_pid"]}" '
            f'data-p="{hit["protein_pid"]}"></rect>')

    axis_elements = []
    axis_y_pos = pile_top + max_lanes * (rect_height + gap) + 10
    tick_height, tick_interval = 5, max(1, int(10 ** math.floor(math.log10(max(2, ref_seq_len) / 2))))
    if ref_seq_len / tick_interval > 10: tick_interval *= 2
    tick_interval = max(1, int(tick_interval))

    axis_elements.append(f'<line class="axis-line" x1="0" y1="{axis_y_pos}" x2="{ref_seq_len}" y2="{axis_y_pos}"></line>')
    for i in range(0, ref_seq_len + 1, tick_interval):
        anchor = ' style="text-anchor: start;"' if i == 0 else ''
        axis_elements.append(f'<line class="tick" x1="{i}" y1="{axis_y_pos}" x2="{i}" y2="{axis_y_pos + tick_height}"></line>')
        axis_elements.append(f'<text class="tick-label" x="{i}" y="{axis_y_pos + tick_height + 12}"{anchor}>{i}</text>')

    axis_label_y = axis_y_pos + tick_height + 28
    axis_elements.append(f'<text class="axis-label" x="{ref_seq_len / 2}" y="{axis_label_y}">Reference position (bp)</text>')
    svg_width, svg_height = max(ref_seq_len + 2, 420), axis_label_y + 10

    own_row, group_rows = final_rows
    group_reads = sum(int(float(g['reads'])) for g in group_rows if str(g['reads']).replace('.', '', 1).isdigit())
    group_text = ""
    if group_rows:
        largest = max(g['group_size'] for g in group_rows)
        group_text = (f'Also a member of {len(group_rows)} group row{"s" if len(group_rows) != 1 else ""} '
                      f'({group_reads:,} reads shared with up to {largest - 1} other references the reads cannot tell apart)')
    if own_row is not None:
        final_card = (f'<p>{html.escape(str(own_row["reads"]))} read{"" if str(own_row["reads"]) == "1" else "s"}</p>'
                      f'<small>{html.escape(str(own_row["copies_per_cell"]))} copies per cell in its own row. {group_text}</small>')
    elif group_rows:
        final_card = f'<p>shared only</p><small>No row of its own. {group_text}</small>'
    else:
        final_card = '<p>&ndash;</p><small>no row in the final table</small>'

    stats = {
        'reference': aro_id, 'length': ref_seq_len, 'reads': original_read_count, 'unique_reads': unique_reads,
        'lateral_coverage_pct': round(coverage, 2), 'mean_depth': round(mean_depth, 2), 'max_depth': max(depth) if depth else 0,
        'reads_drawn': len(drawn),
        'final_reads': own_row['reads'] if own_row else '', 'final_copies_per_cell': own_row['copies_per_cell'] if own_row else '',
        'group_rows': len(group_rows), 'group_reads': group_reads,
    }
    section = f"""
    <section class="aro-plot" id="plot-{html.escape(aro_id, quote=True)}" data-aro-id="{html.escape(aro_id, quote=True)}">
        <h2>{html.escape(aro_id)}</h2>
        <div class="stats-grid">
            <div class="stat-card"><h3>Reference length</h3><p>{ref_seq_len:,} bp</p></div>
            <div class="stat-card"><h3>Reads mapped</h3><p class="reads-mapped-count">{original_read_count:,}</p><small>{unique_reads:,} match only this reference</small></div>
            <div class="stat-card"><h3>Lateral coverage</h3><p class="lateral-coverage-value">{coverage:.2f}%</p><small>mean depth {mean_depth:,.1f}&times;</small></div>
            <div class="stat-card"><h3>In the final table</h3>{final_card}</div>
        </div>
        <div class="chart-container">
            <svg width="{svg_width}" height="{svg_height}" data-ref-len="{ref_seq_len}">
                {''.join(svg_elements)}
                {''.join(axis_elements)}
            </svg>
        </div>
    </section>
    """
    return section, stats, drawn

def generate_main_html_page(family_name, plot_htmls, hits_by_query_id_json, pid_type,
                            sample_name=None, index_href=None, table_rows=None, max_reads_drawn=DEFAULT_MAX_READS_DRAWN):
    """Generates the page for one gene family."""
    table_rows = table_rows or []
    families = [f.strip() for f in family_name.split(';') if f.strip()]
    family_note = ""
    if len(families) > 1:
        family_note = (f'<p class="note">The references on this page are annotated in CARD with {len(families)} gene families at once '
                       f'({html.escape("; ".join(families))}), so they are reported under the combined name.</p>')
    tsv_name = sanitize_filename(f"{sample_name}_{family_name}" if sample_name else family_name).replace('.html', '_homscan_references.tsv')
    buttons = f'<button class="tool" onclick="downloadTSV({html.escape(js_data(tsv_name), quote=True)}, REFERENCE_TABLE)">Download TSV</button>'
    legend_html = f"""
    <div class="legend">
        <h3>Read colour: {pid_type} identity</h3>
        <div class="legend-item"><span class="color-box" style="background-color: var(--pid-full);"></span>100%</div>
        <div class="legend-item"><span class="color-box" style="background-color: var(--pid-partial);"></span>below 100%</div>
        <div class="legend-item"><span class="color-box" style="background-color: transparent; border: 2px solid var(--unique);"></span>Outlined: matches only this reference</div>
        <div class="legend-item"><span class="color-box" style="background-color: var(--depth);"></span>Depth</div>
    </div>
    <p class="note">Each bar is one read placed where it aligns on the reference; hover for its identity and click to list every
    reference it matched. An outlined read matched this reference and no other, so it is direct evidence for it. The depth profile
    and all statistics use every mapped read; when a reference has more than {max_reads_drawn:,} reads, the pile shows an evenly spaced
    sample of them (the label above each pile says how many). &ldquo;Reads mapped&rdquo; counts every read that matches the reference.
    In the final table a reference has its own row for the reads that point to it alone, while reads shared by several references
    are reported in group rows; the last card gives both.</p>
    {family_note}
    """

    script_js = f"""
    <script>
        const allReadHits = {hits_by_query_id_json};
        const REFERENCE_TABLE = {js_data(table_rows)};
        const tooltip = document.getElementById('tooltip');
        const readDetailsOverlay = document.getElementById('read-details-overlay');
        const detailsReadIdSpan = document.getElementById('details-read-id');
        const detailsContentDiv = document.getElementById('details-content');
        const escapeHtml = s => String(s).replace(/[&<>"]/g, c => ({{'&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;'}}[c]));

        function showTooltip(evt) {{
            const rect = evt.target;
            if (!rect.classList || !rect.classList.contains('r')) {{ hideTooltip(); return; }}
            const entry = allReadHits.reads[rect.getAttribute('data-q')];
            const matched = entry ? entry[0] : 1;
            tooltip.style.display = 'block';
            tooltip.innerHTML = 'Read: ' + escapeHtml(rect.getAttribute('data-q')) + '<br>Nucleotide identity: '
                + escapeHtml(rect.getAttribute('data-n')) + '%<br>Protein identity: ' + escapeHtml(rect.getAttribute('data-p'))
                + '%<br>' + (matched === 1 ? 'Matches only this reference' : 'Matches ' + matched + ' references') + ' (click to list)';
            let x = evt.clientX + 15, y = evt.clientY + 15;
            if (x + tooltip.offsetWidth > window.innerWidth) x = evt.clientX - tooltip.offsetWidth - 15;
            if (y + tooltip.offsetHeight > window.innerHeight) y = evt.clientY - tooltip.offsetHeight - 15;
            tooltip.style.left = x + 'px'; tooltip.style.top = y + 'px';
        }}
        function hideTooltip() {{ tooltip.style.display = 'none'; }}

        function showReadDetails(event, queryId) {{
            event.stopPropagation();
            const entry = allReadHits.reads[queryId];
            if (!entry) return;
            const rect = event.target;
            const currentPlotAroId = rect.closest('.aro-plot').getAttribute('data-aro-id');
            const rows = entry.slice(1).map(h => ({{
                family: allReadHits.references[h[0]][1], reference: allReadHits.references[h[0]][0],
                nuc: h[1], prot: h[2], uniqueness: h[3], position: h[4], length: h[5] }}));
            if (!rows.some(r => r.reference === currentPlotAroId)) {{
                // The reference being viewed is not among the best matches listed: show it first.
                rows.unshift({{ family: '', reference: currentPlotAroId, nuc: rect.getAttribute('data-n'), prot: rect.getAttribute('data-p'),
                               uniqueness: '', position: rect.getAttribute('x'), length: rect.getAttribute('width') }});
            }}
            detailsReadIdSpan.textContent = queryId;
            let contentHtml = '<p>This read matched ' + entry[0] + ' reference' + (entry[0] === 1 ? '' : 's') + '.'
                + (entry[0] > rows.length ? ' The ' + rows.length + ' listed are its best matches and the one being viewed.' : '') + '</p>'
                + '<table class="data"><thead><tr><th>Gene family</th><th>Reference</th><th>Nucleotide identity</th><th>Protein identity</th><th>Uniqueness</th><th>Position</th><th>Length</th></tr></thead><tbody>';
            rows.forEach(r => {{
                const rowClassAttr = (r.reference === currentPlotAroId) ? 'highlight-row' : '';
                contentHtml += `<tr class="${{rowClassAttr}}"><td>${{escapeHtml(r.family)}}</td><td>${{escapeHtml(r.reference)}}</td><td>${{escapeHtml(r.nuc)}}%</td><td>${{escapeHtml(r.prot)}}%</td><td>${{escapeHtml(r.uniqueness)}}</td><td>${{escapeHtml(r.position)}}</td><td>${{escapeHtml(r.length)}}</td></tr>`;
            }});
            contentHtml += '</tbody></table>';
            detailsContentDiv.innerHTML = contentHtml;
            readDetailsOverlay.style.visibility = 'visible'; readDetailsOverlay.style.opacity = 1;
        }}

        function hideReadDetails() {{
            readDetailsOverlay.style.opacity = 0;
            setTimeout(() => {{ readDetailsOverlay.style.visibility = 'hidden'; }}, 300);
        }}

        // One set of listeners for the whole page rather than three per read.
        const plots = document.getElementById('plots-container');
        plots.addEventListener('mousemove', showTooltip);
        plots.addEventListener('mouseleave', hideTooltip);
        plots.addEventListener('click', evt => {{
            if (evt.target.classList && evt.target.classList.contains('r')) showReadDetails(evt, evt.target.getAttribute('data-q'));
        }});
        readDetailsOverlay.addEventListener('click', evt => {{ if (evt.target === readDetailsOverlay) hideReadDetails(); }});
    </script>
    """

    overview_rows = "".join(
        f'<tr><td><a href="#plot-{html.escape(str(r[0]), quote=True)}">{html.escape(short_name(str(r[0])))}</a></td>'
        f'<td class="num">{r[1]:,}</td><td class="num">{r[2]:,}</td><td class="num">{r[3]:,}</td>'
        f'<td class="num">{r[4]:.2f}</td><td class="num">{r[5]:.2f}</td>'
        f'<td class="num">{html.escape(str(r[8]))}</td><td class="num">{html.escape(str(r[9]))}</td><td class="num">{r[10]}</td></tr>'
        for r in table_rows[1:])
    overview = f"""
        <div class="table-wrap"><table class="data sortable">
            <thead><tr><th>Reference</th><th class="num">Length (bp)</th><th class="num">Reads mapped</th><th class="num">Only this reference</th>
            <th class="num">Lateral coverage %</th><th class="num">Mean depth</th><th class="num">Reads in own row</th><th class="num">Copies per cell</th><th class="num">Group rows</th></tr></thead>
            <tbody>{overview_rows}</tbody>
        </table></div>""" if len(table_rows) > 2 else ""

    return f"""<!DOCTYPE html>
<html lang="en-GB">
<head>
    <meta charset="UTF-8">
    <title>{html.escape(title_text('HomScan: ' + family_name, sample_name))}</title>
    <style>{THEME_CSS}
        .chart-container {{ position: relative; width: 100%; margin-top: 1em; overflow-x: auto; padding: 1em; background-color: var(--plot-bg); border-radius: 8px; box-shadow: 0 4px 6px var(--shadow); border: 1px solid var(--plot-border); }}
        #tooltip {{ position: fixed; display: none; padding: 8px; background-color: var(--panel); color: var(--text); border-radius: 4px; font-size: 12px; pointer-events: none; white-space: pre; z-index: 1000; border: 1px solid var(--panel-border); }}
        rect.r {{ cursor: pointer; }}
        rect.full {{ fill: var(--pid-full); }}
        rect.part {{ fill: var(--pid-partial); }}
        rect.uniq {{ stroke: var(--unique); stroke-width: 1.5; }}
        .depth-area {{ fill: var(--depth); opacity: 0.85; }}
        .depth-rule {{ stroke: var(--axis); stroke-width: 1; stroke-dasharray: 3 3; }}
        .depth-label {{ font-size: 11px; fill: var(--text-muted); }}
        .axis-line, .tick {{ stroke: var(--axis); stroke-width: 2; }}
        .tick-label {{ font-size: 10px; fill: var(--text-muted); text-anchor: middle; }}
        .axis-label {{ font-size: 12px; font-weight: 600; fill: var(--text-muted); text-anchor: middle; }}
        .aro-plot {{ margin-top: 3em; }}
        .aro-plot h2 {{ word-break: break-all; }}
        .color-box {{ width: 15px; height: 15px; border-radius: 3px; border: 1px solid var(--panel-border); display: inline-block; }}
        #read-details-overlay {{ position: fixed; top: 0; left: 0; width: 100%; height: 100%; background-color: rgba(0, 0, 0, 0.6); display: flex; justify-content: center; align-items: center; z-index: 2000; visibility: hidden; opacity: 0; transition: opacity 0.3s ease; }}
        #read-details-box {{ background-color: var(--panel); padding: 25px; border-radius: 10px; box-shadow: 0 10px 30px rgba(0, 0, 0, 0.5); max-width: 1200px; max-height: 80vh; overflow-y: auto; position: relative; color: var(--text); border: 1px solid var(--panel-border); }}
        #read-details-box h3 {{ margin-top: 0; color: var(--text-strong); border-bottom: 1px solid var(--rule); padding-bottom: 10px; font-size: 1.3em; }}
        #read-details-box .close-btn {{ position: absolute; top: 10px; right: 15px; font-size: 1.5em; cursor: pointer; color: var(--text-muted); }}
        #read-details-box .highlight-row {{ background-color: var(--highlight) !important; font-weight: bold; }}
        @media print {{ .aro-plot {{ break-inside: avoid; }} }}
    </style>
    <script>{THEME_JS}</script>
</head>
<body>
    <div class="container">
        {page_bar(sample_name, index_href, "HomScan", buttons)}
        <h1>Coverage for gene family: {html.escape(family_name)}</h1>
        {legend_html}
        {overview}
        <div id="plots-container">{''.join(plot_htmls)}</div>
    </div>
    <div id="tooltip"></div>
    <div id="read-details-overlay"><div id="read-details-box"><span class="close-btn" onclick="hideReadDetails()">×</span><h3>Read details: <span id="details-read-id"></span></h3><div id="details-content"></div></div></div>
    {script_js}
</body>
</html>
"""

def main():
    parser = argparse.ArgumentParser(description="Generates interactive HTML visualisations.")
    parser.add_argument("-i", "--input-files", required=True, help="Comma-delimited list of input `_hits.tsv` files.")
    parser.add_argument("-d", "--db", required=True, help="Path to reference FASTA.")
    parser.add_argument("--metadata", required=True, help="Path to metadata TSV.")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory.")
    parser.add_argument("--pid-cutoff", type=float, default=0.9)
    parser.add_argument("--pid-type", choices=['protein', 'nucleotide'], default='protein')
    parser.add_argument("--summary-table", default=None, help="Final HomScan table, used to show how each reference is reported after ambiguous reads are resolved.")
    parser.add_argument("--max-reads-drawn", type=int, default=DEFAULT_MAX_READS_DRAWN, help=f"Most reads drawn individually per reference (fewer on pages with very many references); 0 draws all. Depth and statistics always use every read. Default: {DEFAULT_MAX_READS_DRAWN}")
    parser.add_argument("--sample-name", default=None, help="Sample name shown in the page title and heading.")
    parser.add_argument("--index-href", default=None, help="Relative link back to the sample overview page.")
    args = parser.parse_args()

    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]
    os.makedirs(args.output_dir, exist_ok=True)

    aro_to_family, family_to_aros = load_metadata(args.metadata)
    all_hits_raw, all_hits_by_query_id = load_all_hits(input_files_list, args.pid_cutoff, args.pid_type, aro_to_family)
    row_of_reference = load_summary_table(args.summary_table)

    hits_by_aro = defaultdict(list)
    for hit in all_hits_raw:
        hits_by_aro[hit['ARO_matched']].append(hit)
    ref_sequences = load_reference_sequences(args.db, list(hits_by_aro.keys()))

    families_with_hits = {h['AMR_Gene_Family'] for h in all_hits_raw}
    table_header = ["reference", "length_bp", "reads_mapped", "reads_matching_only_this_reference", "lateral_coverage_%",
                    "mean_depth", "max_depth", "reads_drawn", "final_table_reads", "final_table_copies_per_cell",
                    "group_rows_including_this_reference", "reads_in_those_group_rows"]

    for family_name in sorted(list(families_with_hits)):
        print(BColors.cyan(f"--- Processing family: {family_name} ---"))
        aros_to_visualise = [a for a in family_to_aros.get(family_name, []) if a in hits_by_aro]
        if not aros_to_visualise: continue

        plot_htmls, table_rows, page_query_ids = [], [table_header], set()
        drawn_limit = reads_drawn_limit(args.max_reads_drawn, len(aros_to_visualise))
        for aro_id in sorted(aros_to_visualise):
            if aro_id not in ref_sequences: continue
            aro_hits = hits_by_aro[aro_id]
            section, stats, drawn = generate_single_aro_plot_html(
                aro_id, len(ref_sequences[aro_id]), aro_hits, len(aro_hits), args.pid_type,
                drawn_limit, final_rows_for(row_of_reference, aro_id))
            plot_htmls.append(section)
            table_rows.append([stats['reference'], stats['length'], stats['reads'], stats['unique_reads'],
                               stats['lateral_coverage_pct'], stats['mean_depth'], stats['max_depth'], stats['reads_drawn'],
                               stats['final_reads'], stats['final_copies_per_cell'], stats['group_rows'], stats['group_reads']])
            page_query_ids.update(h['query_id'] for h in drawn)

        # Embed the click-on-a-read details only for reads drawn on this page, and in
        # compact form, so that page size follows what the page shows.
        page_hits = compact_read_hits(page_query_ids, all_hits_by_query_id)
        final_html = generate_main_html_page(family_name, plot_htmls, js_data(page_hits), args.pid_type,
                                             args.sample_name, args.index_href, table_rows, drawn_limit)
        output_path = os.path.join(args.output_dir, sanitize_filename(family_name))
        with open(output_path, 'w', encoding='utf-8') as f: f.write(final_html)

if __name__ == "__main__":
    main()
