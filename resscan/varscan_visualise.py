#!/usr/bin/env python3

# resscan/varscan_visualise.py
import argparse
import os
import sys
import csv
import html
import re
from collections import Counter, defaultdict
try:
    from resscan.utils import BColors
    from resscan.report_common import THEME_CSS, THEME_JS, page_bar, title_text, js_data, sanitize_filename
except ImportError:
    from utils import BColors
    from report_common import THEME_CSS, THEME_JS, page_bar, title_text, js_data, sanitize_filename

# Bases shown either side of the read in the default (windowed) view.
WINDOW_FLANK = 10

_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TO_AA = {a + b + c: _AMINO[16 * i + 4 * j + k]
               for i, a in enumerate(_BASES) for j, b in enumerate(_BASES) for k, c in enumerate(_BASES)}


def load_metadata(metadata_path):
    """Loads metadata to map AROs to AMR_Gene_Families."""
    aro_to_family = {}
    if not os.path.exists(metadata_path):
        print(BColors.red(f"Error: Metadata file not found at '{metadata_path}'"), file=sys.stderr)
        sys.exit(1)

    try:
        with open(metadata_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                aro_to_family[row['Sequence_ID']] = row.get('AMR_Gene_Family', 'Unknown_Family')
    except Exception as e:
        print(BColors.red(f"Error reading metadata: {e}"), file=sys.stderr)
        sys.exit(1)
    return aro_to_family

def load_hits(hits_path, cutoff):
    """
    Returns {query_id: hit row} for the reads that pass the identity cutoff, or
    None if the hits file is unavailable (in which case nothing is filtered).
    """
    if not hits_path or not os.path.exists(hits_path):
        return None

    passing = {}
    try:
        with open(hits_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            cutoff_percent = cutoff * 100.0
            for row in reader:
                try:
                    if float(row['nucleotide_pid']) >= cutoff_percent:
                        passing[row['query_id']] = row
                except (ValueError, KeyError):
                    continue
    except Exception as e:
        print(BColors.red(f"Warning: Could not read hits file for PID filtering: {e}"), file=sys.stderr)
        return None
    return passing

def load_pid_whitelist(hits_path, cutoff):
    """Returns a set of query_ids that pass the PID cutoff."""
    hits = load_hits(hits_path, cutoff)
    return None if hits is None else set(hits)

def parse_mutations(header):
    """Returns (variant type, [(original, position, new), ...]) for an alignment header."""
    type_match = re.search(r'__([RVO])__', header)
    var_type = type_match.group(1) if type_match else 'R'

    if 'mutation:' in header:
        mutation_field = header.rsplit('mutation:', 1)[1]
    else:
        # Older alignment files: fall back to the mutation suffix of the entry name.
        mutation_field = re.split(r'__[RVO]__', header)[-1]
    mutations = [(o, int(p), n) for o, p, n in re.findall(r'([a-zA-Z*])(\d+)([a-zA-Z*])', mutation_field)]
    return var_type, mutations

def parse_header_info(header):
    """
    Determines if the variant is Nucleotide (R) or Protein (V/O) and its positions (plural).

    Positions are read from the mutation field only. Entry names can themselves
    look like a substitution (the "P2x" in PBP2x, "P1a" in PBP1a, "E2A" in GE2A),
    and reading the whole header would mark those as mutation sites.
    """
    var_type, mutations = parse_mutations(header)
    return var_type, [position for _, position, _ in mutations]

def parse_alignments(file_path):
    """Parses the multiline variant alignment format."""
    alignments = []
    current_block = None

    with open(file_path, 'r') as f:
        for line in f:
            clean = line.strip()
            if not clean: continue

            if clean.startswith('>'):
                if current_block: alignments.append(current_block)
                current_block = {
                    'header': clean,
                    'refseq': "",
                    'queseq': "",
                    'raw_ref_line': "",
                    'aro': clean.split('|')[2].strip() if '|' in clean else "Unknown"
                }
            elif clean.startswith('REFSEQ:'):
                if not current_block['raw_ref_line']:
                    current_block['raw_ref_line'] = clean
                current_block['refseq'] += re.sub(r'^REFSEQ(\s*\(\d+\))?:\s*', '', clean)
            elif clean.startswith('QUESEQ:'):
                current_block['queseq'] += clean.replace('QUESEQ:', '').strip()

        if current_block: alignments.append(current_block)
    return alignments

def generate_marker_line(ref_seq, var_type, positions, start_pos):
    """Generates the track mapping multiple non-gapped coordinates to gapped alignment indices."""
    dots = ["." for _ in range(len(ref_seq))]
    if not positions: return "".join(dots)

    all_target_nucs = set()
    for pos in positions:
        if var_type == 'R':
            all_target_nucs.add(pos)
        else:
            nuc_base = ((pos - 1) * 3) + 1
            all_target_nucs.update([nuc_base, nuc_base + 1, nuc_base + 2])

    current_nuc_coord = start_pos
    for i, char in enumerate(ref_seq):
        if char != '-':
            if current_nuc_coord in all_target_nucs:
                dots[i] = "^"
            current_nuc_coord += 1

    return "".join(dots)

def highlight_mismatches(ref_seq, que_seq):
    """
    Returns the read line as HTML with every aligned position that differs from
    the reference (substitution, insertion or deletion) wrapped in a highlight.
    Unaligned positions (dots) are left as they are.
    """
    parts = []
    run, run_is_mismatch = [], False
    for ref_char, que_char in zip(ref_seq, que_seq):
        is_mismatch = que_char != '.' and que_char.upper() != ref_char.upper()
        if run and is_mismatch != run_is_mismatch:
            text = html.escape("".join(run))
            parts.append(f'<span class="mismatch">{text}</span>' if run_is_mismatch else text)
            run = []
        run.append(que_char)
        run_is_mismatch = is_mismatch
    if run:
        text = html.escape("".join(run))
        parts.append(f'<span class="mismatch">{text}</span>' if run_is_mismatch else text)
    return "".join(parts)

def describe_sites(ref_seq, que_seq, var_type, mutations, start_coord=1):
    """
    For each mutation of the entry, what the reference and the read carry at the
    site: the base for rRNA entries, the codon and its amino acid for proteins.
    """
    column_of = {}
    coord = start_coord
    for i, char in enumerate(ref_seq):
        if char != '-':
            column_of[coord] = i
            coord += 1

    sites = []
    for original, position, new in mutations:
        label = f"{original}{position}{new}"
        if var_type == 'R':
            col = column_of.get(position)
            if col is None: continue
            sites.append({'label': label, 'reference': ref_seq[col].upper(), 'read': que_seq[col].upper(),
                          'reference_aa': '', 'read_aa': ''})
        else:
            cols = [column_of.get((position - 1) * 3 + k) for k in (1, 2, 3)]
            if None in cols: continue
            ref_codon = "".join(ref_seq[c] for c in cols).upper()
            read_codon = "".join(que_seq[c] for c in cols).upper()
            sites.append({'label': label, 'reference': ref_codon, 'read': read_codon,
                          'reference_aa': CODON_TO_AA.get(ref_codon, '?'), 'read_aa': CODON_TO_AA.get(read_codon, '?')})
    return sites

def build_block_view(block, flank=WINDOW_FLANK):
    """
    Works out what to display for one read: the part of the alignment around the
    read (the default view), where that window sits on the reference, and the
    evidence at each mutation site.
    """
    ref_seq, que_seq = block['refseq'], block['queseq']
    var_type, mutations = parse_mutations(block['header'])
    coord_match = re.search(r'\((\d+)\)', block['raw_ref_line'])
    start_coord = int(coord_match.group(1)) if coord_match else 1
    marker = generate_marker_line(ref_seq, var_type, [p for _, p, _ in mutations], start_coord)

    aligned_cols = [i for i, char in enumerate(que_seq) if char != '.']
    first, last = (aligned_cols[0], aligned_cols[-1]) if aligned_cols else (0, len(que_seq) - 1)
    win_start, win_end = max(0, first - flank), min(len(ref_seq), last + 1 + flank)

    bases_before = lambda col: sum(1 for char in ref_seq[:col] if char != '-')
    aligned_pairs = [(r, q) for r, q in zip(ref_seq[first:last + 1], que_seq[first:last + 1])]
    return {
        'query_id': block['header'].split('|')[0].strip().lstrip('>').strip(),
        'aro': block['aro'],
        'var_type': var_type,
        'mutation': "_".join(f"{o}{p}{n}" for o, p, n in mutations),
        'ref_window': ref_seq[win_start:win_end],
        'que_window': que_seq[win_start:win_end],
        'marker_window': marker[win_start:win_end],
        # Reference bases before and up to the end of the window (0-based, gap-free),
        # so the page can rebuild the full-length view around it.
        'window_ref_from': bases_before(win_start),
        'window_ref_to': bases_before(win_end),
        'ref_first_coord': start_coord,
        'aligned_from': start_coord + bases_before(first),
        'aligned_to': start_coord + bases_before(last + 1) - 1,
        'mismatches': sum(1 for r, q in aligned_pairs if r.upper() != q.upper()),
        'aligned_columns': len(aligned_pairs),
        'sites': describe_sites(ref_seq, que_seq, var_type, mutations, start_coord),
        'full_reference': ref_seq.replace('-', ''),
    }

def _site_text(site):
    if site['reference_aa']:
        return (f"{site['label']}: reference {site['reference']} ({site['reference_aa']}), "
                f"read {site['read']} ({site['read_aa']})")
    return f"{site['label']}: reference {site['reference']}, read {site['read']}"

def _fragment_id(query_id):
    return re.sub(r'_\d+$', '', query_id)

def generate_family_page(family_name, blocks, hits=None, sample_name=None, index_href=None):
    """Generates the page for one gene family: a summary per entry, then every supporting read."""
    hits = hits or {}
    views_by_entry = defaultdict(list)
    for b in blocks:
        view = build_block_view(b)
        hit = hits.get(view['query_id'])
        view['identity'] = hit.get('nucleotide_pid', '') if hit else ''
        views_by_entry[view['aro']].append(view)

    summary_rows, sections, references, tsv_rows = [], [], {}, []
    tsv_rows.append(["read", "entry", "mutation", "nucleotide_identity_%", "reference_from", "reference_to",
                     "site", "reference_state", "read_state", "reference_residue", "read_residue"])
    for entry_number, aro in enumerate(sorted(views_by_entry), 1):
        views = views_by_entry[aro]
        anchor = f"entry-{entry_number}"
        references[aro] = views[0]['full_reference']
        identities = [float(v['identity']) for v in views if v['identity']]
        mean_identity = f"{sum(identities) / len(identities):.2f}" if identities else ""
        fragments = len({_fragment_id(v['query_id']) for v in views})

        # What the reads carry at each site, e.g. "E350Q: CAA (Q) x 300, CAG (Q) x 14".
        observed = defaultdict(Counter)
        for v in views:
            for site in v['sites']:
                state = f"{site['read']} ({site['read_aa']})" if site['read_aa'] else site['read']
                observed[site['label']][state] += 1
        evidence = "; ".join(
            f"{label}: " + ", ".join(f"{state} &times; {n:,}" for state, n in counts.most_common())
            for label, counts in observed.items())

        entry_label = aro.split('__', 1)[1] if '__' in aro else aro
        summary_rows.append(
            f'<tr><td><a href="#{anchor}">{html.escape(entry_label)}</a></td>'
            f'<td>{html.escape(views[0]["mutation"].replace("_", ", "))}</td>'
            f'<td class="num">{len(views):,}</td><td class="num">{fragments:,}</td>'
            f'<td class="num">{mean_identity}</td><td>{evidence}</td></tr>')

        block_html = []
        for v in views:
            sites_text = " &nbsp;|&nbsp; ".join(html.escape(_site_text(s)) for s in v['sites'])
            identity_text = f"identity {html.escape(v['identity'])}% &nbsp;|&nbsp; " if v['identity'] else ""
            marker_html = html.escape(v['marker_window']).replace('^', '<span class="marker">^</span>')
            lines = "\n".join([
                '<span class="pos-text">POS: <span class="pos-line"></span></span>',
                f'<span class="ref-text">REF: <span class="ref-line">{html.escape(v["ref_window"])}</span></span>',
                f'<span class="que-text">QUE: <span class="que-line">{highlight_mismatches(v["ref_window"], v["que_window"])}</span></span>',
                f'<span class="mrk-text">MRK: <span class="mrk-line">{marker_html}</span></span>',
            ])
            block_html.append(
                f'<div class="block">'
                f'<span class="header-label"><b>{html.escape(v["query_id"])}</b> &nbsp;|&nbsp; {identity_text}'
                f'reference positions {v["aligned_from"]:,}&ndash;{v["aligned_to"]:,}</span>'
                f'<span class="site-label">{sites_text}</span>'
                f'<div class="seq-container" data-aro="{html.escape(aro, quote=True)}" '
                f'data-from="{v["window_ref_from"]}" data-to="{v["window_ref_to"]}" '
                f'data-first-coord="{v["ref_first_coord"]}">{lines}</div></div>')
            for s in v['sites'] or [{'label': '', 'reference': '', 'read': '', 'reference_aa': '', 'read_aa': ''}]:
                tsv_rows.append([v['query_id'], aro, v['mutation'], v['identity'], v['aligned_from'], v['aligned_to'],
                                 s['label'], s['reference'], s['read'], s['reference_aa'], s['read_aa']])

        sections.append(
            f'<h2 id="{anchor}">{html.escape(entry_label)}</h2>'
            f'<p class="note">{len(views):,} supporting read{"s" if len(views) != 1 else ""} '
            f'from {fragments:,} fragment{"s" if fragments != 1 else ""}. '
            f'Reference length {len(views[0]["full_reference"]):,} bases.</p>' + "".join(block_html))

    total_fragments = len({_fragment_id(v['query_id']) for views in views_by_entry.values() for v in views})
    tsv_name = sanitize_filename(f"{sample_name}_{family_name}" if sample_name else family_name).replace('.html', '_varscan_reads.tsv')
    buttons = ('<button class="tool" id="view-toggle" onclick="toggleFullReference()">Show full reference</button>'
               f'<button class="tool" onclick="downloadTSV({html.escape(js_data(tsv_name), quote=True)}, READ_TABLE)">Download TSV</button>')

    return f"""<!DOCTYPE html>
<html lang="en-GB">
<head>
    <meta charset="UTF-8">
    <title>{html.escape(title_text('VarScan: ' + family_name, sample_name))}</title>
    <style>{THEME_CSS}
        .block {{ background-color: var(--panel); border-radius: 8px; padding: 1em 1.2em; margin-bottom: 1.2em; box-shadow: 0 4px 6px var(--shadow); border: 1px solid var(--panel-border); border-left: 5px solid var(--axis); }}
        .header-label {{ display: block; font-family: monospace; font-size: 0.9em; color: var(--text-muted); }}
        .header-label b {{ color: var(--text); }}
        .site-label {{ display: block; font-family: monospace; font-size: 0.9em; color: var(--text); margin: 0.3em 0 0.7em 0; padding-bottom: 0.5em; border-bottom: 1px solid var(--rule); }}
        .seq-container {{ position: relative; font-family: 'Courier New', Courier, monospace; white-space: pre; background-color: var(--plot-bg); color: var(--text); padding: 0.8em 1em; border-radius: 5px; overflow-x: auto; line-height: 1.5; border: 1px solid var(--plot-border); }}
        .pos-text {{ color: var(--text-muted); }}
        .ref-text {{ color: var(--ref); }}
        .que-text {{ color: var(--que); }}
        .mrk-text {{ color: var(--text-muted); }}
        .mismatch {{ color: var(--mismatch-fg); background-color: var(--mismatch-bg); }}
        .marker {{ color: var(--marker); font-weight: bold; }}
        @media print {{ .block {{ break-inside: avoid; }} }}
    </style>
    <script>{THEME_JS}</script>
</head>
<body>
    <div class="container">
        {page_bar(sample_name, index_href, "VarScan", buttons)}
        <h1>Variant confirmation: {html.escape(family_name)}</h1>
        <div class="legend">
            <h3>Key</h3>
            <span class="legend-item"><b style="color: var(--ref);">REF</b> Reference carrying the mutation</span>
            <span class="legend-item"><b style="color: var(--que);">QUE</b> Read</span>
            <span class="legend-item"><b style="color: var(--marker);">^</b> Mutation site</span>
            <span class="legend-item"><b class="mismatch">&nbsp;A&nbsp;</b> Read differs from reference</span>
        </div>
        <p class="note">Each read is shown with {WINDOW_FLANK} reference bases either side; use &ldquo;Show full reference&rdquo; to see it
        in the context of the whole reference. POS gives the reference coordinate of the base above each number's first digit.
        Only reads that carry every mutation of an entry are listed. Reads that cover a site without the mutation are
        neither shown nor counted here.</p>
        <div class="stats-grid">
            <div class="stat-card"><h3>Supporting reads</h3><p>{len(blocks):,}</p><small>after the identity filter</small></div>
            <div class="stat-card"><h3>Fragments</h3><p>{total_fragments:,}</p><small>read pairs counted once</small></div>
            <div class="stat-card"><h3>Entries with support</h3><p>{len(views_by_entry):,}</p></div>
        </div>
        <table class="data sortable">
            <thead><tr><th>Entry</th><th>Mutation(s)</th><th class="num">Reads</th><th class="num">Fragments</th><th class="num">Mean identity %</th><th>Observed in the reads at each site</th></tr></thead>
            <tbody>{''.join(summary_rows)}</tbody>
        </table>
        {''.join(sections)}
    </div>
    <script>
        const REFERENCES = {js_data(references)};
        const READ_TABLE = {js_data(tsv_rows)};
        let fullView = false;

        // Coordinate ruler: a number is written starting above every tenth reference base.
        function ruler(gappedReference, firstCoord) {{
            const out = new Array(gappedReference.length).fill(' ');
            let coord = firstCoord, freeFrom = 0;
            for (let i = 0; i < gappedReference.length; i++) {{
                if (gappedReference[i] === '-') continue;
                if (coord % 10 === 0 && i >= freeFrom) {{
                    const label = String(coord);
                    if (i + label.length <= out.length) {{
                        for (let k = 0; k < label.length; k++) out[i + k] = label[k];
                        freeFrom = i + label.length + 1;
                    }}
                }}
                coord++;
            }}
            return out.join('');
        }}
        function centreOnMarker(box) {{
            const marker = box.querySelector('.marker');
            if (marker) box.scrollLeft = Math.max(0, marker.offsetLeft - box.clientWidth / 2);
        }}
        function render(box) {{
            const ref = box.querySelector('.ref-line'), que = box.querySelector('.que-line'), mrk = box.querySelector('.mrk-line');
            if (!box.dataset.refWindow) {{                       // remember the windowed view once
                box.dataset.refWindow = ref.textContent; box.dataset.queWindow = que.innerHTML; box.dataset.mrkWindow = mrk.innerHTML;
            }}
            const from = parseInt(box.dataset.from), to = parseInt(box.dataset.to), first = parseInt(box.dataset.firstCoord);
            const full = REFERENCES[box.dataset.aro] || '';
            if (fullView && full) {{
                const left = '.'.repeat(from), right = '.'.repeat(Math.max(0, full.length - to));
                ref.textContent = full.slice(0, from) + box.dataset.refWindow + full.slice(to);
                que.innerHTML = left + box.dataset.queWindow + right;
                mrk.innerHTML = left + box.dataset.mrkWindow + right;
                box.querySelector('.pos-line').textContent = ruler(ref.textContent, first);
            }} else {{
                ref.textContent = box.dataset.refWindow; que.innerHTML = box.dataset.queWindow; mrk.innerHTML = box.dataset.mrkWindow;
                box.querySelector('.pos-line').textContent = ruler(ref.textContent, first + from);
            }}
            centreOnMarker(box);
        }}
        function toggleFullReference() {{
            fullView = !fullView;
            document.getElementById('view-toggle').textContent = fullView ? 'Show read window' : 'Show full reference';
            document.querySelectorAll('.seq-container').forEach(render);
        }}
        window.addEventListener('load', () => document.querySelectorAll('.seq-container').forEach(render));
    </script>
</body>
</html>
"""

def main():
    parser = argparse.ArgumentParser(description="Generates interactive VarScan reports.")
    parser.add_argument("--variant-hits", help="Path to variant hits TSV (required for PID filtering).")
    parser.add_argument("--variant-alignments", required=True, help="Path to alignments TXT.")
    parser.add_argument("--metadata", required=True, help="Path to metadata TSV.")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory.")
    parser.add_argument("--pid-cutoff", type=float, default=0.95, help="Minimum PID to display in report (0.0-1.0). Default: 0.95")
    parser.add_argument("--sample-name", default=None, help="Sample name shown in the page title and heading.")
    parser.add_argument("--index-href", default=None, help="Relative link back to the sample overview page.")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(BColors.cyan(f"--- Loading data for {os.path.basename(args.variant_alignments)} ---"))
    aro_to_family = load_metadata(args.metadata)
    all_blocks = parse_alignments(args.variant_alignments)

    hits = load_hits(args.variant_hits, args.pid_cutoff)
    if hits is not None:
        print(BColors.cyan(f"--- Filtering alignments by PID >= {args.pid_cutoff*100.0}% ---"))

    family_groups = defaultdict(list)
    filtered_count = 0

    for b in all_blocks:
        q_id = b['header'].split('|')[0].strip().replace('> ', '')

        if hits is not None and q_id not in hits:
            continue

        filtered_count += 1
        family_name = aro_to_family.get(b['aro'], "Uncategorised_Variants")
        family_groups[family_name].append(b)

    print(BColors.green(f"--- Processing {filtered_count} blocks after filtering. ---"))

    for family, blocks in family_groups.items():
        print(BColors.green(f"Writing report: {family}"))
        html_out = generate_family_page(family, blocks, hits, args.sample_name, args.index_href)
        with open(os.path.join(args.output_dir, sanitize_filename(family)), 'w', encoding='utf-8') as f:
            f.write(html_out)

    print(BColors.cyan(f"--- Done. Reports saved in {args.output_dir} ---"))

if __name__ == "__main__":
    main()
