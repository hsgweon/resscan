#!/usr/bin/env python3

# resscan/homscan_tabulate_and_normalise.py
import argparse
import os
import sys
import csv
import re
from collections import defaultdict, Counter
try:
    from resscan.utils import BColors, load_metadata, load_uscg_metrics, load_sequencing_metrics
except ImportError:
    from utils import BColors, load_metadata, load_uscg_metrics, load_sequencing_metrics

ARO_PATTERN = re.compile(r'ARO_\d+')

def load_priors(filepath):
    """Loads MAP priors from a tab-separated file."""
    if not filepath:
        return {}
    print(BColors.cyan(f"--- Loading MAP priors from: {filepath} ---"))
    priors = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 1:
                    priors[parts[0]] = float(parts[1]) if len(parts) > 1 else 1.0
    except FileNotFoundError:
        print(BColors.yellow(f"Warning: Priors file not found at '{filepath}'."), file=sys.stderr)
        return {}
    print(BColors.green(f"--- Loaded {len(priors)} prior entries."))
    return priors

def load_and_filter_hits(input_files, pid_cutoff_fraction, pid_type):
    """Reads all input TSV files, filters by the chosen PID type, and returns hits grouped by query_id."""
    hits_by_query = defaultdict(list)
    pid_cutoff_percent = pid_cutoff_fraction * 100.0
    pid_column = f"{pid_type}_pid"
    print(BColors.cyan(f"--- Processing {len(input_files)} input file(s) with {pid_column} >= {pid_cutoff_percent:.2f}%... ---"))

    for file_path in input_files:
        try:
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    try:
                        if float(row[pid_column]) >= pid_cutoff_percent:
                            hits_by_query[row['query_id']].append(row)
                    except (ValueError, KeyError):
                        continue
        except FileNotFoundError:
            print(BColors.yellow(f"Warning: Input file not found, skipping: {file_path}"), file=sys.stderr)
    return hits_by_query

def group_hits_by_final_category(hits_by_query, aro_to_family, consensus_threshold):
    """Resolves each read into its final category and groups the actual hit data by that category."""
    hits_by_final_category = defaultdict(list)
    for query_id, hits in hits_by_query.items():
        unique_aros_hit = {h['ARO_matched'] for h in hits}

        if len(unique_aros_hit) == 1:
            aro = list(unique_aros_hit)[0]
            family = aro_to_family.get(aro, "Unknown_Family")
            category_key = (family, aro)
        else:
            family_aro_pairs = [f"{aro_to_family.get(aro, 'Unknown_Family')};{aro}" for aro in unique_aros_hit]
            sorted_pairs_str = "|".join(sorted(family_aro_pairs))
            aro_description = f"multiple[{sorted_pairs_str}]"

            family_counts = Counter(aro_to_family.get(aro, "Unknown_Family") for aro in unique_aros_hit)
            most_common_family, count = family_counts.most_common(1)[0]
            family_description = most_common_family if count / len(unique_aros_hit) >= consensus_threshold else "multiple"
            category_key = (family_description, aro_description)

        hits_by_final_category[category_key].extend(hits)

    return hits_by_final_category

def calculate_coverage_for_hit_group(list_of_hits, aro_to_length):
    """Calculates the lateral coverage for a specific group of hits."""
    hits_by_aro = defaultdict(list)
    for hit in list_of_hits:
        hits_by_aro[hit['ARO_matched']].append(hit)

    coverages = []
    for aro, hits in hits_by_aro.items():
        coverage_set = set()
        for hit in hits:
            try:
                start = int(hit['position_on_ref'])
                length = int(hit['nucleotide_denominator'])
                coverage_set.update(range(start, start + length))
            except (ValueError, KeyError):
                continue

        total_len = aro_to_length.get(aro)
        if total_len and total_len > 0:
            coverages.append((len(coverage_set) / total_len) * 100.0)

    return sum(coverages) / len(coverages) if coverages else 0.0

def get_annotation(aro, list_of_hits, aro_to_annotation):
    """Returns (ARO_Name, Drug_Class, Resistance_Mechanism) for a single or multiple ARO entry."""
    if aro.startswith("multiple["):
        constituent_aros = ARO_PATTERN.findall(aro)
    elif aro == "multiple":
        constituent_aros = list({h['ARO_matched'] for h in list_of_hits})
    else:
        constituent_aros = [aro]

    annotations = [aro_to_annotation.get(a, {}) for a in constituent_aros]

    def _join(field):
        vals = [a.get(field, 'N/A') for a in annotations if a.get(field, 'N/A') != 'N/A']
        return " | ".join(dict.fromkeys(vals)) if vals else "N/A"

    return _join("ARO_Name"), _join("Drug_Class"), _join("Resistance_Mechanism")

def calculate_normalised_metrics(hits_by_final_category, aro_to_length, aro_to_annotation,
                                  uscg_rpk, uscg_fpk, total_bases, total_reads=None, total_fragments=None):
    """Calculates all metrics for each AMR entry, including inline annotation."""
    print(BColors.cyan("--- Calculating metrics... ---"))
    results = []
    total_bases_in_gb = (total_bases / 1e9) if total_bases > 0 else 0
    total_reads_in_millions = (total_reads / 1e6) if (total_reads and total_reads > 0) else None
    total_fragments_in_millions = (total_fragments / 1e6) if (total_fragments and total_fragments > 0) else None

    for category_key, list_of_hits in hits_by_final_category.items():
        family, aro = category_key

        all_read_ids = [h['query_id'] for h in list_of_hits]
        read_count = len(set(all_read_ids))
        fragment_ids = {re.sub(r'_\d+$', '', r_id) for r_id in all_read_ids}
        fragment_count = len(fragment_ids)

        coverage = calculate_coverage_for_hit_group(list_of_hits, aro_to_length)

        gene_length_bp = 0
        if aro.startswith("multiple[") or aro == "multiple":
            all_aros_in_group = {h['ARO_matched'] for h in list_of_hits}
            lengths = [aro_to_length.get(a, 0) for a in all_aros_in_group]
            valid_lengths = [l for l in lengths if l > 0]
            if valid_lengths:
                gene_length_bp = sum(valid_lengths) / len(valid_lengths)
        else:
            gene_length_bp = aro_to_length.get(aro, 0)

        rpk = (read_count / (gene_length_bp / 1000.0)) if gene_length_bp > 0 else 0.0
        rpkg = (rpk / total_bases_in_gb) if total_bases_in_gb > 0 else 0.0
        rpkm_val = f"{rpk / total_reads_in_millions:.4f}" if total_reads_in_millions else "NA"
        rpkpc_val, rpkpmc_val, rpkpgc_val = "NA", "NA", "NA"
        if uscg_rpk is not None and uscg_rpk > 1e-9:
            rpkpc = rpk / uscg_rpk
            rpkpc_val = f"{rpkpc:.4f}"
            rpkpmc_val = f"{rpkpc * 1_000_000:.2f}"
            rpkpgc_val = f"{rpkpc * 1_000_000_000:.2f}"
        elif uscg_rpk is not None and rpk == 0.0:
            rpkpc_val, rpkpmc_val, rpkpgc_val = "0.0000", "0.00", "0.00"

        fpk = (fragment_count / (gene_length_bp / 1000.0)) if gene_length_bp > 0 else 0.0
        fpkg = (fpk / total_bases_in_gb) if total_bases_in_gb > 0 else 0.0
        fpkm_val = f"{fpk / total_fragments_in_millions:.4f}" if total_fragments_in_millions else "NA"
        fpkpc_val, fpkpmc_val, fpkpgc_val = "NA", "NA", "NA"
        if uscg_fpk is not None and uscg_fpk > 1e-9:
            fpkpc = fpk / uscg_fpk
            fpkpc_val = f"{fpkpc:.4f}"
            fpkpmc_val = f"{fpkpc * 1_000_000:.2f}"
            fpkpgc_val = f"{fpkpc * 1_000_000_000:.2f}"
        elif uscg_fpk is not None and fpk == 0.0:
            fpkpc_val, fpkpmc_val, fpkpgc_val = "0.0000", "0.00", "0.00"

        aro_name, drug_class, resistance_mechanism = get_annotation(aro, list_of_hits, aro_to_annotation)

        results.append({
            "family": family,
            "aro": aro,
            "ARO_Name": aro_name,
            "Drug_Class": drug_class,
            "Resistance_Mechanism": resistance_mechanism,
            "Read_Count": read_count,
            "Fragment_Count": fragment_count,
            "Lateral_Coverage_%": f"{coverage:.2f}",
            "Gene_Length_bp": f"{gene_length_bp:.2f}" if isinstance(gene_length_bp, float) else str(gene_length_bp),
            "RPK": f"{rpk:.4f}", "FPK": f"{fpk:.4f}",
            "RPKG": f"{rpkg:.4f}", "FPKG": f"{fpkg:.4f}",
            "RPKM": rpkm_val, "FPKM": fpkm_val,
            "RPKPC": rpkpc_val, "FPKPC": fpkpc_val,
            "RPKPMC": rpkpmc_val, "FPKPMC": fpkpmc_val,
            "RPKPGC": rpkpgc_val, "FPKPGC": fpkpgc_val,
        })

    results.sort(key=lambda x: (x['family'], x['aro']))
    return results

def aggregate_for_clean_summary(hits_by_final_category):
    """Aggregates detailed hit groups into clean groups."""
    clean_hit_groups = defaultdict(list)
    for (family, aro), list_of_hits in hits_by_final_category.items():
        if aro.startswith("multiple["):
            clean_key = (family, "multiple")
        else:
            clean_key = (family, aro)
        clean_hit_groups[clean_key].extend(list_of_hits)
    return clean_hit_groups

def run_map_solver(detailed_results, metric_col, priors, base_prior, prior_strength):
    """Runs an EM-like iterative solver to distribute ambiguous read abundance across AROs."""
    print(BColors.cyan(f"--- Running MAP solver using metric: {metric_col} ---"))

    def get_val(row, col):
        v = row.get(col, 0)
        try:
            f = float(v)
            return f if f > 0 else 0.0
        except (ValueError, TypeError):
            return 0.0

    actual_col = metric_col
    for fallback in [metric_col, 'RPKG', 'RPK', 'Read_Count']:
        if any(get_val(r, fallback) > 0 for r in detailed_results):
            actual_col = fallback
            break

    if actual_col != metric_col:
        print(BColors.yellow(f"Warning: '{metric_col}' is all zero/NA; falling back to '{actual_col}' for MAP solver."), file=sys.stderr)

    all_aros = set()
    for r in detailed_results:
        all_aros.update(ARO_PATTERN.findall(r['aro']))
    all_aros_list = sorted(all_aros)

    if not all_aros_list:
        return {}

    unique_abun = defaultdict(float)
    for r in detailed_results:
        if not r['aro'].startswith('multiple'):
            aros = ARO_PATTERN.findall(r['aro'])
            if aros:
                unique_abun[aros[0]] += get_val(r, actual_col)

    static_abun = {a: unique_abun.get(a, 0.0) for a in all_aros_list}
    for aro, hist_count in priors.items():
        if aro in static_abun:
            static_abun[aro] += hist_count * prior_strength * 0.1
    for a in all_aros_list:
        static_abun[a] += base_prior * 1e-6

    ambiguous = []
    for r in detailed_results:
        if r['aro'].startswith('multiple['):
            aros = sorted(set(ARO_PATTERN.findall(r['aro'])))
            val = get_val(r, actual_col)
            if aros:
                ambiguous.append((aros, val))

    if not ambiguous:
        return dict(static_abun)

    final_abun = dict(static_abun)
    for _ in range(50):
        prev = dict(final_abun)
        ambig_contrib = defaultdict(float)
        for aros, val in ambiguous:
            weights = [prev.get(a, 0.0) for a in aros]
            total_w = sum(weights)
            props = [w / total_w for w in weights] if total_w > 1e-9 else [1.0 / len(aros)] * len(aros)
            for a, p in zip(aros, props):
                ambig_contrib[a] += val * p

        new_abun = {a: static_abun[a] + ambig_contrib.get(a, 0.0) for a in all_aros_list}
        if all(abs(new_abun.get(a, 0) - prev.get(a, 0)) < 1e-6 for a in all_aros_list):
            final_abun = new_abun
            break
        final_abun = new_abun

    return final_abun

def get_map_allocation(family, detailed_results, final_abun):
    """Returns (Top_ARO, Allocation_Proportions) for a family's multiple row."""
    family_rows = [r for r in detailed_results if r['family'] == family]
    all_family_aros = sorted(set(
        a for r in family_rows for a in ARO_PATTERN.findall(r['aro'])
    ))

    if not all_family_aros:
        return "None", "No Evidence"

    fam_weights = {a: final_abun.get(a, 0.0) for a in all_family_aros}
    total_w = sum(fam_weights.values())

    if total_w > 1e-9:
        props = {a: w / total_w for a, w in fam_weights.items()}
        top_aro = max(props, key=props.get)
        if list(props.values()).count(props[top_aro]) > 1:
            top_aro = 'N/A'
        sorted_p = sorted(props.items(), key=lambda x: x[1], reverse=True)
        alloc = ",".join([f"{a}:{p * 100:.1f}%" for a, p in sorted_p])
    else:
        top_aro = "None"
        alloc = "No Evidence"

    return top_aro, alloc

def write_summary_file(results, output_path, include_gene_length=False, include_map_cols=False):
    """Writes a normalised summary file."""
    print(BColors.cyan(f"--- Writing normalised summary to: {output_path} ---"))
    try:
        with open(output_path, 'w', newline='') as f:
            header = ["AMR_Gene_Family", "ARO", "ARO_Name", "Drug_Class", "Resistance_Mechanism",
                      "Read_Count", "Fragment_Count", "Lateral_Coverage_%"]
            if include_gene_length:
                header.append("Gene_Length_bp")
            header += ["RPK", "FPK", "RPKG", "FPKG", "RPKM", "FPKM",
                       "RPKPC", "FPKPC", "RPKPMC", "FPKPMC", "RPKPGC", "FPKPGC"]
            if include_map_cols:
                header += ["Top_ARO", "Allocation_Proportions"]

            writer = csv.DictWriter(f, fieldnames=header, delimiter='\t', extrasaction='ignore')
            writer.writeheader()

            for row_data in results:
                row_to_write = {
                    "AMR_Gene_Family": row_data["family"],
                    "ARO": row_data["aro"],
                    "ARO_Name": row_data["ARO_Name"],
                    "Drug_Class": row_data["Drug_Class"],
                    "Resistance_Mechanism": row_data["Resistance_Mechanism"],
                    "Read_Count": row_data["Read_Count"],
                    "Fragment_Count": row_data["Fragment_Count"],
                    "Lateral_Coverage_%": row_data["Lateral_Coverage_%"],
                    "Gene_Length_bp": row_data["Gene_Length_bp"],
                    "RPK": row_data["RPK"], "FPK": row_data["FPK"],
                    "RPKG": row_data["RPKG"], "FPKG": row_data["FPKG"],
                    "RPKM": row_data["RPKM"], "FPKM": row_data["FPKM"],
                    "RPKPC": row_data["RPKPC"], "FPKPC": row_data["FPKPC"],
                    "RPKPMC": row_data["RPKPMC"], "FPKPMC": row_data["FPKPMC"],
                    "RPKPGC": row_data["RPKPGC"], "FPKPGC": row_data["FPKPGC"],
                    "Top_ARO": row_data.get("Top_ARO", ""),
                    "Allocation_Proportions": row_data.get("Allocation_Proportions", ""),
                }
                writer.writerow(row_to_write)
        print(BColors.green(f"--- Successfully wrote report: {output_path} ---"))
    except IOError as e:
        print(BColors.red(f"Error writing report file: {e}"), file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description="Combines, tabulates, and normalizes AMR homology hits, producing both a detailed and a clean summary."
    )
    parser.add_argument("-i", "--input-files", required=True, help="Comma-delimited list of input `_hits.tsv` files.")
    parser.add_argument("--metadata", required=True, help="Path to the AMR database metadata file.")
    parser.add_argument("--uscg-report", required=True, help="Path to the USCG quantification report.")
    parser.add_argument("--total-bases-file", required=True, help="Path to a file containing the total number of bases.")
    parser.add_argument("--tmp-dir", default=".", help="Directory to store the output report files.")
    parser.add_argument("--output-prefix", required=True, help="Prefix for the output summary files.")
    parser.add_argument("--pid-cutoff", type=float, default=0.9, help="Minimum percent identity to consider a hit (0.0-1.0 scale).")
    parser.add_argument("--pid-type", choices=['protein', 'nucleotide'], default='protein', help="PID type to use for filtering. Default: protein")
    parser.add_argument("--consensus", type=float, default=0.9, help="Minimum fraction of non-unique hits that must map to the same gene family to reach a consensus.")
    parser.add_argument("--map-priors-file", default=None, help="Path to a tab-separated file of priors for the MAP solver.")
    parser.add_argument("--map-metric-column", default='FPKPMC', help="The numeric column to use for MAP abundance resolution. Default: FPKPMC")
    parser.add_argument("--map-base-prior", type=float, default=1.0, help="Baseline prior pseudo-count. Default: 1.0")
    parser.add_argument("--map-prior-strength", type=float, default=1.0, help="Multiplier for the influence of the priors file. Default: 1.0")
    args = parser.parse_args()

    if not 0.0 <= args.pid_cutoff <= 1.0 or not 0.0 <= args.consensus <= 1.0:
        print(BColors.red("Error: --pid-cutoff and --consensus must be fractions between 0.0 and 1.0."), file=sys.stderr)
        sys.exit(1)

    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]
    if not input_files_list:
        print(BColors.red("Error: No input files provided."), file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.tmp_dir, exist_ok=True)

    aro_to_family, aro_to_length, aro_to_annotation = load_metadata(args.metadata)
    uscg_rpk, uscg_fpk = load_uscg_metrics(args.uscg_report)
    total_bases, total_reads, total_fragments = load_sequencing_metrics(args.total_bases_file)
    priors = load_priors(args.map_priors_file)

    hits_by_query = load_and_filter_hits(input_files_list, args.pid_cutoff, args.pid_type)
    print(BColors.green(f"--- Filtered and processed {len(hits_by_query)} unique reads from all files. ---"))

    hits_by_final_category = group_hits_by_final_category(hits_by_query, aro_to_family, args.consensus)
    print(BColors.green(f"--- Aggregated reads into {len(hits_by_final_category)} detailed categories. ---"))

    detailed_normalised_results = calculate_normalised_metrics(
        hits_by_final_category, aro_to_length, aro_to_annotation,
        uscg_rpk, uscg_fpk, total_bases, total_reads, total_fragments)
    detailed_output_path = os.path.join(args.tmp_dir, f"{args.output_prefix}_homscan_detailed.tsv")
    write_summary_file(detailed_normalised_results, detailed_output_path, include_gene_length=True, include_map_cols=False)

    clean_hit_groups = aggregate_for_clean_summary(hits_by_final_category)
    print(BColors.green(f"--- Aggregated detailed results into {len(clean_hit_groups)} clean categories. ---"))

    clean_normalised_results = calculate_normalised_metrics(
        clean_hit_groups, aro_to_length, aro_to_annotation,
        uscg_rpk, uscg_fpk, total_bases, total_reads, total_fragments)

    final_abun = run_map_solver(detailed_normalised_results, args.map_metric_column, priors, args.map_base_prior, args.map_prior_strength)

    for row in clean_normalised_results:
        if row['aro'] == 'multiple':
            top_aro, alloc = get_map_allocation(row['family'], detailed_normalised_results, final_abun)
            row['Top_ARO'] = top_aro
            row['Allocation_Proportions'] = alloc
        else:
            row['Top_ARO'] = row['aro']
            row['Allocation_Proportions'] = f"{row['aro']}:100.0%"

    clean_output_path = os.path.join(args.tmp_dir, f"{args.output_prefix}_homscan.tsv")
    write_summary_file(clean_normalised_results, clean_output_path, include_gene_length=False, include_map_cols=True)

    print(BColors.green("\n\n--- All Processing Complete ---"))
    print("The following output files were generated:")
    print(f"- {detailed_output_path}")
    print(f"- {clean_output_path}")

if __name__ == "__main__":
    main()
