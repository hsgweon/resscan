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

class BColors:
    """A helper class to add color to terminal output."""
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    enabled = sys.stdout.isatty()

    @classmethod
    def _colorize(cls, color_code, text):
        return f"{color_code}{text}{cls.ENDC}" if cls.enabled else text
    @classmethod
    def cyan(cls, text): return cls._colorize(cls.OKCYAN, text)
    @classmethod
    def green(cls, text): return cls._colorize(cls.OKGREEN, text)
    @classmethod
    def red(cls, text): return cls._colorize(cls.FAIL, text)
    @classmethod
    def yellow(cls, text): return cls._colorize(cls.WARNING, text)

def sanitize_filename(name):
    """Converts a string into a safe filename, consistent with varscan_visualise."""
    # Replace spaces and semicolons with underscores
    name = name.replace(' ', '_').replace(';', '_').replace('/', '_')
    # Remove any other characters that aren't alphanumeric, underscores, or hyphens
    name = re.sub(r'[^\w\-.]', '', name)
    return name.lower() + ".html"

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

def get_color_for_pid(pid):
    if pid >= 100.0:
        return "#00ff41" # Matrix Green for 100%
    else:
        return "#ffcc00" # Goldenrod for < 100%

def generate_single_aro_plot_html(aro_id, ref_seq_len, hits, original_read_count, pid_type):
    """Generates the HTML for a single ARO's plot with dark theme."""
    stats = {
        "ref_len": f"{ref_seq_len:,}",
        "num_reads_all": f"{original_read_count:,}",
        "coverage": "N/A"
    }

    hits_for_initial_layout = sorted(hits, key=lambda x: float(x['nucleotide_pid']), reverse=True)
    rect_height, gap = 8, 2
    lanes = [0]
    
    js_hits_data = []
    for hit in hits_for_initial_layout:
        start = int(hit['position_on_ref'])
        length = int(hit['nucleotide_denominator'])
        assigned_lane = -1
        for i, end_pos in enumerate(lanes):
            if start >= end_pos:
                assigned_lane = i
                break
        if assigned_lane == -1:
            assigned_lane = len(lanes)
            lanes.append(0)
        lanes[assigned_lane] = start + length
        hit_copy = hit.copy()
        hit_copy['draw_lane'] = assigned_lane
        js_hits_data.append(hit_copy)

    max_lanes = len(lanes) if lanes else 1
    svg_elements, axis_elements = [], []
    svg_filter_defs = """<defs><filter id="goldGlow" x="-50%" y="-50%" width="200%" height="200%"><feDropShadow dx="0" dy="0" stdDeviation="3" flood-color="#FFD700" flood-opacity="0.8"/></filter></defs>"""

    pid_column_for_color = f"{pid_type}_pid"

    for hit in js_hits_data:
        pid_for_color = float(hit[pid_column_for_color])
        color = get_color_for_pid(pid_for_color)
        y_pos = (max_lanes - hit['draw_lane'] - 1) * (rect_height + gap)
        uniqueness = hit.get('uniqueness', 'UNKNOWN')
        stroke_attr = 'stroke="#B8860B" stroke-width="1.5"' if uniqueness == 'UNIQUE' else ''
        filter_attr = 'filter="url(#goldGlow)"' if uniqueness == 'UNIQUE' else ''
        
        tooltip_text_short = (
            f"Read: {html.escape(hit['query_id'])} ({uniqueness})<br>"
            f"Nuc PID: {hit['nucleotide_pid']}%<br>"
            f"Prot PID: {hit['protein_pid']}%"
        )
        
        svg_elements.append(
            f'<rect x="{hit["position_on_ref"]}" y="{y_pos}" '
            f'width="{hit["nucleotide_denominator"]}" height="{rect_height}" '
            f'fill="{color}" rx="1" {stroke_attr} {filter_attr} '
            f'data-query-id="{html.escape(hit["query_id"])}"'
            f'data-aro-matched="{html.escape(hit["ARO_matched"])}"'
            f'data-nuc-pid="{hit["nucleotide_pid"]}"'
            f'data-prot-pid="{hit["protein_pid"]}"'
            f'data-pos-on-ref="{hit["position_on_ref"]}"'
            f'data-aln-len="{hit["nucleotide_denominator"]}"'
            f'data-tooltip-short="{tooltip_text_short}" '
            f'onmousemove="showTooltip(event)" onmouseout="hideTooltip()" '
            f'onclick="showReadDetails(event, \'{html.escape(hit["query_id"])}\')"></rect>'
        )

    axis_y_pos = max_lanes * (rect_height + gap) + 10
    tick_height, tick_interval = 5, max(1, int(10 ** math.floor(math.log10(ref_seq_len / 2))))
    if ref_seq_len / tick_interval > 10: tick_interval *= 2
    tick_interval = max(1, int(tick_interval))

    axis_elements.append(f'<line class="axis-line" x1="0" y1="{axis_y_pos}" x2="{ref_seq_len}" y2="{axis_y_pos}"></line>')
    for i in range(0, ref_seq_len + 1, tick_interval):
        axis_elements.append(f'<line class="tick" x1="{i}" y1="{axis_y_pos}" x2="{i}" y2="{axis_y_pos + tick_height}"></line>')
        axis_elements.append(f'<text class="tick-label" x="{i}" y="{axis_y_pos + tick_height + 12}">{i}</text>')
    
    axis_label_y = axis_y_pos + tick_height + 28
    axis_elements.append(f'<text class="axis-label" x="{ref_seq_len / 2}" y="{axis_label_y}">Reference Position (bp)</text>')

    svg_width, svg_height = ref_seq_len, axis_label_y + 10

    return f"""
    <section class="aro-plot" id="plot-{html.escape(aro_id)}" data-aro-id="{html.escape(aro_id)}">
        <h2>{html.escape(aro_id)}</h2>
        <div class="stats-grid">
            <div class="stat-card"><h3>Reference Length</h3><p>{stats['ref_len']} bp</p></div>
            <div class="stat-card"><h3>Reads Mapped</h3><p class="reads-mapped-count">{stats['num_reads_all']}</p></div>
            <div class="stat-card"><h3>Lateral Coverage</h3><p class="lateral-coverage-value">{stats['coverage']}</p></div>
        </div>
        <div class="chart-container">
            <svg width="{svg_width}" height="{svg_height}" data-ref-len="{ref_seq_len}">
                {svg_filter_defs}
                {''.join(svg_elements)}
                {''.join(axis_elements)}
            </svg>
        </div>
    </section>
    """

def generate_main_html_page(family_name, plot_htmls, hits_by_query_id_json, pid_type):
    """Generates the final HTML page with full dark theme."""
    legend_html = f"""
    <div class="legend">
        <h3>PID: {pid_type.capitalize()}</h3>
        <div class="legend-item"><span class="color-box" style="background-color: #00ff41;"></span>100.0%</div>
        <div class="legend-item"><span class="color-box" style="background-color: #ffcc00;"></span>< 100.0%</div>
        <div class="legend-item" style="margin-left: 2em;"><span class="color-box" style="background-color: transparent; border: 1.5px solid #B8860B; box-shadow: 0 0 5px #FFD700;"></span>Unique Read</div>
    </div>
    """

    script_js = f"""
    <script>
        const allReadHits = {hits_by_query_id_json};
        const tooltip = document.getElementById('tooltip');
        const readDetailsOverlay = document.getElementById('read-details-overlay');
        const detailsReadIdSpan = document.getElementById('details-read-id');
        const detailsContentDiv = document.getElementById('details-content');

        function showTooltip(evt) {{
            const target = evt.target;
            const tooltipText = target.getAttribute('data-tooltip-short');
            if (!tooltipText) return;
            tooltip.style.display = 'block';
            tooltip.innerHTML = tooltipText.replace(/\\n/g, '<br/>');
            let x = evt.clientX + 15, y = evt.clientY + 15;
            if (x + tooltip.offsetWidth > window.innerWidth) x = evt.clientX - tooltip.offsetWidth - 15;
            if (y + tooltip.offsetHeight > window.innerHeight) y = evt.clientY - tooltip.offsetHeight - 15;
            tooltip.style.left = x + 'px'; tooltip.style.top = y + 'px';
        }}
        function hideTooltip() {{ tooltip.style.display = 'none'; }}

        function showReadDetails(event, queryId) {{
            event.stopPropagation();
            const hitsForThisRead = allReadHits[queryId];
            if (!hitsForThisRead) return;
            detailsReadIdSpan.textContent = queryId;
            let contentHtml = '<p>This read mapped to the following AROs:</p><table><thead><tr><th>Gene Family</th><th>ARO</th><th>Nuc PID</th><th>Prot PID</th><th>Uniqueness</th><th>Pos</th><th>Len</th></tr></thead><tbody>';
            const currentPlotAroId = event.target.closest('.aro-plot').getAttribute('data-aro-id');
            hitsForThisRead.forEach(hit => {{
                let rowClassAttr = (hit.ARO_matched === currentPlotAroId) ? 'highlight-row' : '';
                contentHtml += `<tr class="${{rowClassAttr}}"><td>${{hit.AMR_Gene_Family}}</td><td>${{hit.ARO_matched}}</td><td>${{hit.nucleotide_pid}}%</td><td>${{hit.protein_pid}}%</td><td>${{hit.uniqueness}}</td><td>${{hit.position_on_ref}}</td><td>${{hit.nucleotide_denominator}}</td></tr>`;
            }});
            contentHtml += '</tbody></table>';
            detailsContentDiv.innerHTML = contentHtml;
            readDetailsOverlay.style.visibility = 'visible'; readDetailsOverlay.style.opacity = 1;
        }}

        function hideReadDetails() {{
            readDetailsOverlay.style.opacity = 0;
            setTimeout(() => {{ readDetailsOverlay.style.visibility = 'hidden'; }}, 300);
        }}

        function calculateLateralCoverage(plot) {{
            const svg = plot.querySelector('svg');
            const refLength = parseInt(svg.getAttribute('data-ref-len'));
            const rects = svg.querySelectorAll('rect');
            if (refLength === 0) return '0.00%';
            const coverageSet = new Set();
            rects.forEach(rect => {{
                const start = parseInt(rect.getAttribute('data-pos-on-ref'));
                const length = parseInt(rect.getAttribute('data-aln-len'));
                for (let i = start; i < start + length; i++) coverageSet.add(i);
            }});
            return ((coverageSet.size / refLength) * 100).toFixed(2) + '%';
        }}

        document.addEventListener('DOMContentLoaded', () => {{
            document.querySelectorAll('.aro-plot').forEach(plot => {{
                plot.querySelector('.lateral-coverage-value').textContent = calculateLateralCoverage(plot);
            }});
        }});
    </script>
    """

    return f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>ResScan Visualisation: {html.escape(family_name)}</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 2em; background-color: #121212; color: #e0e0e0; }}
        .container {{ max-width: 1200px; margin: auto; }}
        h1 {{ color: #ffffff; border-bottom: 2px solid #333; padding-bottom: 0.5em; }}
        h2 {{ color: #ffffff; margin-top: 2em; border-bottom: 1px solid #333; padding-bottom: 0.3em; display: flex; align-items: center; gap: 0.5em;}}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1em; margin: 1em 0; }}
        .stat-card {{ background-color: #1e1e1e; border-radius: 8px; padding: 1.2em; box-shadow: 0 4px 6px rgba(0,0,0,0.3); border: 1px solid #333; }}
        .stat-card h3 {{ margin: 0 0 0.5em 0; color: #aaa; font-size: 1em; font-weight: 600; }}
        .stat-card p {{ margin: 0; color: #ffffff; font-size: 1.6em; font-weight: 700; }}
        .chart-container {{ position: relative; width: 100%; margin-top: 1em; overflow-x: auto; padding: 1em; background-color: #000000; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.3); border: 1px solid #222; }}
        #tooltip {{ position: fixed; display: none; padding: 8px; background-color: rgba(0, 0, 0, 0.9); color: white; border-radius: 4px; font-size: 12px; pointer-events: none; white-space: pre; z-index: 1000; border: 1px solid #444; }}
        .axis-line, .tick {{ stroke: #555; stroke-width: 2; }}
        .tick-label {{ font-size: 10px; fill: #888; text-anchor: middle; }}
        .axis-label {{ font-size: 12px; font-weight: 600; fill: #aaa; text-anchor: middle; }}
        .aro-plot {{ margin-top: 3em; }}
        .legend {{ background-color: #1e1e1e; padding: 1em; border-radius: 8px; margin: 2em 0; display: flex; gap: 2em; align-items: center; flex-wrap: wrap; border: 1px solid #333; }}
        .legend h3 {{ margin: 0; padding-right: 1em; border-right: 1px solid #333; color: #fff; }}
        .legend-item {{ display: flex; align-items: center; gap: 0.5em; font-size: 0.9em; color: #ddd; }}
        .color-box {{ width: 15px; height: 15px; border-radius: 3px; border: 1px solid #444; }}
        #read-details-overlay {{ position: fixed; top: 0; left: 0; width: 100%; height: 100%; background-color: rgba(0, 0, 0, 0.7); display: flex; justify-content: center; align-items: center; z-index: 2000; visibility: hidden; opacity: 0; transition: opacity 0.3s ease; }}
        #read-details-box {{ background-color: #1e1e1e; padding: 25px; border-radius: 10px; box-shadow: 0 10px 30px rgba(0, 0, 0, 0.5); max-width: 1200px; max-height: 80vh; overflow-y: auto; position: relative; color: #eee; border: 1px solid #444; }}
        #read-details-box h3 {{ margin-top: 0; color: #fff; border-bottom: 1px solid #333; padding-bottom: 10px; font-size: 1.5em; }}
        #read-details-box .close-btn {{ position: absolute; top: 10px; right: 15px; font-size: 1.5em; cursor: pointer; color: #aaa; }}
        #read-details-box table {{ width: 100%; border-collapse: collapse; margin-top: 15px; }}
        #read-details-box th, #read-details-box td {{ border: 1px solid #333; padding: 8px 12px; text-align: left; }}
        #read-details-box th {{ background-color: #252525; font-weight: 600; color: #fff; }}
        #read-details-box tr:nth-child(even) {{ background-color: #222; }}
        #read-details-box .highlight-row {{ background-color: #2c3e50 !important; font-weight: bold; color: #00ff41; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Coverage for Gene Family: {html.escape(family_name)}</h1>
        {legend_html}
        <div id="plots-container">{''.join(plot_htmls)}</div>
    </div>
    <div id="tooltip"></div>
    <div id="read-details-overlay"><div id="read-details-box"><span class="close-btn" onclick="hideReadDetails()">×</span><h3>Read Details: <span id="details-read-id"></span></h3><div id="details-content"></div></div></div>
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
    args = parser.parse_args()

    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]
    os.makedirs(args.output_dir, exist_ok=True)

    aro_to_family, family_to_aros = load_metadata(args.metadata)
    all_hits_raw, all_hits_by_query_id = load_all_hits(input_files_list, args.pid_cutoff, args.pid_type, aro_to_family)
    
    families_with_hits = {h['AMR_Gene_Family'] for h in all_hits_raw}

    for family_name in sorted(list(families_with_hits)):
        print(BColors.cyan(f"--- Processing family: {family_name} ---"))
        target_aros = set(family_to_aros.get(family_name, []))
        family_hits_raw = [h for h in all_hits_raw if h['ARO_matched'] in target_aros]
        
        hits_by_aro = defaultdict(list)
        for hit in family_hits_raw:
            hits_by_aro[hit['ARO_matched']].append(hit)

        aros_to_visualise = list(hits_by_aro.keys())
        if not aros_to_visualise: continue

        ref_sequences = load_reference_sequences(args.db, aros_to_visualise)
        plot_htmls = []
        for aro_id in sorted(aros_to_visualise):
            if aro_id not in ref_sequences: continue
            aro_hits = hits_by_aro[aro_id]
            plot_htmls.append(generate_single_aro_plot_html(aro_id, len(ref_sequences[aro_id]), aro_hits, len(aro_hits), args.pid_type))

        final_html = generate_main_html_page(family_name, plot_htmls, json.dumps(all_hits_by_query_id), args.pid_type)
        output_path = os.path.join(args.output_dir, sanitize_filename(family_name))
        with open(output_path, 'w', encoding='utf-8') as f: f.write(final_html)

if __name__ == "__main__":
    main()