#!/usr/bin/env python3

# resscan/varscan_visualise.py
import argparse
import os
import sys
import csv
import html
import re
from collections import defaultdict

class BColors:
    """A helper class to add colour to terminal output."""
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

def sanitize_filename(name):
    """Converts a string into a safe filename, replacing complex chars with underscores."""
    name = name.replace(' ', '_').replace(';', '_').replace('/', '_')
    name = re.sub(r'[^\w\-.]', '', name)
    return name.lower() + ".html"

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

def load_pid_whitelist(hits_path, cutoff):
    """Returns a set of query_ids that pass the PID cutoff."""
    whitelist = set()
    if not hits_path or not os.path.exists(hits_path):
        return None
    
    try:
        with open(hits_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            cutoff_percent = cutoff * 100.0
            for row in reader:
                try:
                    if float(row['nucleotide_pid']) >= cutoff_percent:
                        whitelist.add(row['query_id'])
                except (ValueError, KeyError):
                    continue
    except Exception as e:
        print(BColors.red(f"Warning: Could not read hits file for PID filtering: {e}"), file=sys.stderr)
        return None
    return whitelist

def parse_header_info(header):
    """Determines if the variant is Nucleotide (R) or Protein (V/O) and its positions (plural)."""
    type_match = re.search(r'__([RVO])__', header)
    var_type = type_match.group(1) if type_match else 'R'
    
    pos_matches = re.findall(r'[a-zA-Z*](\d+)[a-zA-Z*]', header)
    positions = [int(p) for p in pos_matches]
    
    return var_type, positions

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

def generate_family_page(family_name, blocks):
    """Generates HTML with a full dark theme. All blue accents removed."""
    
    alignment_html_list = []
    for b in blocks:
        var_type, positions = parse_header_info(b['header'])
        coord_match = re.search(r'\((\d+)\)', b['raw_ref_line'])
        start_coord = int(coord_match.group(1)) if coord_match else 1
        
        marker = generate_marker_line(b['refseq'], var_type, positions, start_coord)
        marker_html = html.escape(marker).replace('^', '<span class="marker">^</span>')

        alignment_html_list.append(f"""
        <div class="block">
            <span class="header-label">{html.escape(b['header'])}</span>
            <div class="seq-container">
<span class="ref-text">REF: {html.escape(b['refseq'])}</span>
<span class="que-text">QUE: {html.escape(b['queseq'])}</span>
<span class="mrk-text">MRK: {marker_html}</span>
            </div>
        </div>
        """)

    return f"""
<!DOCTYPE html>
<html lang="en-GB">
<head>
    <meta charset="UTF-8">
    <title>VarScan: {html.escape(family_name)}</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 2em; background-color: #121212; color: #e0e0e0; }}
        .container {{ max-width: 1200px; margin: auto; }}
        h1 {{ color: #ffffff; border-bottom: 2px solid #333; padding-bottom: 0.5em; }}
        .legend {{ background-color: #1e1e1e; padding: 1em; border-radius: 8px; margin: 2em 0; display: flex; gap: 2em; align-items: center; box-shadow: 0 4px 6px rgba(0,0,0,0.3); border: 1px solid #333; }}
        /* Alignment block border changed from blue to dark grey */
        .block {{ background-color: #1e1e1e; border-radius: 8px; padding: 1.5em; margin-bottom: 2em; box-shadow: 0 4px 6px rgba(0,0,0,0.3); border-left: 5px solid #444; border-top: 1px solid #333; border-right: 1px solid #333; border-bottom: 1px solid #333; }}
        .header-label {{ display: block; font-family: monospace; font-size: 0.9em; color: #888; margin-bottom: 1em; border-bottom: 1px solid #333; padding-bottom: 5px; }}
        .seq-container {{ font-family: 'Courier New', Courier, monospace; white-space: pre; background-color: #000000; color: #dcdcdc; padding: 1.2em; border-radius: 5px; overflow-x: auto; line-height: 1.6; border: 1px solid #222; }}
        .ref-text {{ color: #00ff41; }}
        .que-text {{ color: #ffcc00; }}
        .mrk-text {{ color: #ffffff; }}
        .marker {{ color: #ff4757; font-weight: bold; text-shadow: 0 0 8px rgba(255, 71, 87, 0.9); }}
        .stat-card {{ background-color: #1e1e1e; border-radius: 8px; padding: 1.2em; box-shadow: 0 4px 6px rgba(0,0,0,0.3); margin-bottom: 1em; border: 1px solid #333; }}
        .stat-card h3 {{ margin: 0 0 0.5em 0; color: #aaa; font-size: 1em; }}
        .stat-card p {{ margin: 0; color: #ffffff; font-size: 1.6em; font-weight: 700; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>Variant Confirmation: {html.escape(family_name)}</h1>
        <div class="legend">
            <strong>Key:</strong>
            <span><b style="color: #00ff41;">REF</b> Reference</span>
            <span><b style="color: #ffcc00;">QUE</b> Query Read</span>
            <span><b style="color: #ff4757;">^</b> Mutation Site(s)</span>
        </div>
        <div class="stat-card">
            <h3>Total Supporting Reads (Filtered)</h3>
            <p>{len(blocks)}</p>
        </div>
        {''.join(alignment_html_list)}
    </div>
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
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    
    print(BColors.cyan(f"--- Loading data for {os.path.basename(args.variant_alignments)} ---"))
    aro_to_family = load_metadata(args.metadata)
    all_blocks = parse_alignments(args.variant_alignments)
    
    whitelist = load_pid_whitelist(args.variant_hits, args.pid_cutoff)
    if whitelist is not None:
        print(BColors.cyan(f"--- Filtering alignments by PID >= {args.pid_cutoff*100.0}% ---"))
    
    family_groups = defaultdict(list)
    filtered_count = 0

    for b in all_blocks:
        q_id = b['header'].split('|')[0].strip().replace('> ', '')
        
        if whitelist is not None and q_id not in whitelist:
            continue
            
        filtered_count += 1
        family_name = aro_to_family.get(b['aro'], "Uncategorised_Variants")
        family_groups[family_name].append(b)

    print(BColors.green(f"--- Processing {filtered_count} blocks after filtering. ---"))

    for family, blocks in family_groups.items():
        print(BColors.green(f"Writing report: {family}"))
        html_out = generate_family_page(family, blocks)
        with open(os.path.join(args.output_dir, sanitize_filename(family)), 'w', encoding='utf-8') as f:
            f.write(html_out)

    print(BColors.cyan(f"--- Done. Reports saved in {args.output_dir} ---"))

if __name__ == "__main__":
    main()