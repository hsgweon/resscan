#!/usr/bin/env python3

# resscan/varscan_tabulate_and_normalise.py
import argparse
import os
import sys
import csv
import re
from collections import Counter, defaultdict
try:
    from resscan.utils import BColors, load_metadata, load_uscg_metrics, load_sequencing_metrics
except ImportError:
    from utils import BColors, load_metadata, load_uscg_metrics, load_sequencing_metrics

def count_variant_hits_and_fragments(input_files, pid_cutoff_fraction):
    """
    Reads all input _variant_hits.tsv files, filters by PID, and counts confirmed
    reads and unique fragments for each variant ARO.
    """
    read_counts = Counter()
    fragment_sets = defaultdict(set)
    pid_cutoff_percent = pid_cutoff_fraction * 100.0

    print(BColors.cyan(f"--- Processing {len(input_files)} variant hits file(s) with PID >= {pid_cutoff_percent:.2f}%... ---"))
    for file_path in input_files:
        try:
            with open(file_path, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    try:
                        if float(row['nucleotide_pid']) >= pid_cutoff_percent:
                            aro = row['ARO_matched']
                            query_id = row['query_id']
                            read_counts[aro] += 1
                            fragment_id = re.sub(r'_\d+$', '', query_id)
                            fragment_sets[aro].add(fragment_id)
                    except (ValueError, KeyError):
                        print(BColors.yellow(f"Warning: Skipping row in {file_path} due to missing or invalid columns."), file=sys.stderr)
                        continue
        except FileNotFoundError:
            print(BColors.yellow(f"Warning: Input file not found, skipping: {file_path}"), file=sys.stderr)

    return read_counts, fragment_sets

def calculate_and_write_summary(read_counts, fragment_sets, aro_to_family, aro_to_length, aro_to_annotation,
                                uscg_rpk, uscg_fpk, total_bases, total_reads, total_fragments, output_path):
    """Calculates all metrics and writes the final normalised variant summary."""
    print(BColors.cyan(f"--- Writing normalised variant summary to: {output_path} ---"))
    total_bases_in_gb = (total_bases / 1e9) if total_bases > 0 else 0
    total_reads_in_millions = (total_reads / 1e6) if (total_reads and total_reads > 0) else None
    total_fragments_in_millions = (total_fragments / 1e6) if (total_fragments and total_fragments > 0) else None

    results = []
    for aro, r_count in read_counts.items():
        family = aro_to_family.get(aro, "Unknown_Family")
        gene_length_bp = aro_to_length.get(aro, 0)
        f_count = len(fragment_sets.get(aro, set()))
        ann = aro_to_annotation.get(aro, {})

        rpk = (r_count / (gene_length_bp / 1000.0)) if gene_length_bp > 0 else 0.0
        rpkg = (rpk / total_bases_in_gb) if total_bases_in_gb > 0 else 0.0
        rpkm_val = f"{rpk / total_reads_in_millions:.4f}" if total_reads_in_millions else "NA"
        rpkpc_val, rpkpmc_val, rpkpgc_val = "NA", "NA", "NA"
        if uscg_rpk is not None and uscg_rpk > 1e-9:
            rpkpc = rpk / uscg_rpk
            rpkpc_val = f"{rpkpc:.4f}"
            rpkpmc_val = f"{rpkpc * 1_000_000:.2f}"
            rpkpgc_val = f"{rpkpc * 1_000_000_000:.2f}"
        elif uscg_rpk is not None and uscg_rpk > 1e-9 and rpk == 0.0:
            rpkpc_val, rpkpmc_val, rpkpgc_val = "0.0000", "0.00", "0.00"

        fpk = (f_count / (gene_length_bp / 1000.0)) if gene_length_bp > 0 else 0.0
        fpkg = (fpk / total_bases_in_gb) if total_bases_in_gb > 0 else 0.0
        fpkm_val = f"{fpk / total_fragments_in_millions:.4f}" if total_fragments_in_millions else "NA"
        fpkpc_val, fpkpmc_val, fpkpgc_val = "NA", "NA", "NA"
        if uscg_fpk is not None and uscg_fpk > 1e-9:
            fpkpc = fpk / uscg_fpk
            fpkpc_val = f"{fpkpc:.4f}"
            fpkpmc_val = f"{fpkpc * 1_000_000:.2f}"
            fpkpgc_val = f"{fpkpc * 1_000_000_000:.2f}"
        elif uscg_fpk is not None and uscg_fpk > 1e-9 and fpk == 0.0:
            fpkpc_val, fpkpmc_val, fpkpgc_val = "0.0000", "0.00", "0.00"

        results.append({
            "family": family,
            "aro": aro,
            "ARO_Name": ann.get("ARO_Name", "N/A"),
            "Drug_Class": ann.get("Drug_Class", "N/A"),
            "Resistance_Mechanism": ann.get("Resistance_Mechanism", "N/A"),
            "Read_Count": r_count,
            "Fragment_Count": f_count,
            "RPK": f"{rpk:.4f}", "FPK": f"{fpk:.4f}",
            "RPKG": f"{rpkg:.4f}", "FPKG": f"{fpkg:.4f}",
            "RPKM": rpkm_val, "FPKM": fpkm_val,
            "RPKPC": rpkpc_val, "FPKPC": fpkpc_val,
            "RPKPMC": rpkpmc_val, "FPKPMC": fpkpmc_val,
            "RPKPGC": rpkpgc_val, "FPKPGC": fpkpgc_val,
        })

    results.sort(key=lambda x: (x['family'], x['aro']))

    try:
        with open(output_path, 'w', newline='') as f:
            header = [
                "AMR_Gene_Family", "ARO", "ARO_Name", "Drug_Class", "Resistance_Mechanism",
                "Read_Count", "Fragment_Count",
                "RPK", "FPK", "RPKG", "FPKG", "RPKM", "FPKM",
                "RPKPC", "FPKPC", "RPKPMC", "FPKPMC", "RPKPGC", "FPKPGC"
            ]
            writer = csv.DictWriter(f, fieldnames=header, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for row_data in results:
                writer.writerow({
                    "AMR_Gene_Family": row_data["family"],
                    "ARO": row_data["aro"],
                    "ARO_Name": row_data["ARO_Name"],
                    "Drug_Class": row_data["Drug_Class"],
                    "Resistance_Mechanism": row_data["Resistance_Mechanism"],
                    "Read_Count": row_data["Read_Count"],
                    "Fragment_Count": row_data["Fragment_Count"],
                    "RPK": row_data["RPK"], "FPK": row_data["FPK"],
                    "RPKG": row_data["RPKG"], "FPKG": row_data["FPKG"],
                    "RPKM": row_data["RPKM"], "FPKM": row_data["FPKM"],
                    "RPKPC": row_data["RPKPC"], "FPKPC": row_data["FPKPC"],
                    "RPKPMC": row_data["RPKPMC"], "FPKPMC": row_data["FPKPMC"],
                    "RPKPGC": row_data["RPKPGC"], "FPKPGC": row_data["FPKPGC"],
                })
        print(BColors.green(f"--- Successfully wrote report: {output_path} ---"))
    except IOError as e:
        print(BColors.red(f"Error writing report file: {e}"), file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(
        description="Tabulates and normalises confirmed variant hits into a summary report."
    )
    parser.add_argument("-i", "--input-files", required=True, help="Comma-delimited list of input `_variant_hits.tsv` files.")
    parser.add_argument("--metadata", required=True, help="Path to the AMR database metadata file (TSV format, must include 'SeqNucLength' column).")
    parser.add_argument("--uscg-report", required=True, help="Path to the USCG quantification report containing the overall RPK and FPK.")
    parser.add_argument("--total-bases-file", required=True, help="Path to a file containing the total number of bases for RPKG calculation.")
    parser.add_argument("--tmp-dir", default=".", help="Directory to store the output report file. Defaults to the current directory.")
    parser.add_argument("--output-prefix", required=True, help="Prefix for the output summary file (e.g., 'MyProject').")
    parser.add_argument("--pid-cutoff", type=float, default=0.95, help="Minimum nucleotide PID to consider a hit (0.0-1.0 scale). Default: 0.95")
    args = parser.parse_args()

    if not 0.0 <= args.pid_cutoff <= 1.0:
        print(BColors.red("Error: --pid-cutoff must be a fraction between 0.0 and 1.0."), file=sys.stderr)
        sys.exit(1)

    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]
    if not input_files_list:
        print(BColors.red("Error: No input files provided."), file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.tmp_dir, exist_ok=True)

    aro_to_family, aro_to_length, aro_to_annotation = load_metadata(args.metadata)
    uscg_rpk, uscg_fpk = load_uscg_metrics(args.uscg_report)
    total_bases, total_reads, total_fragments = load_sequencing_metrics(args.total_bases_file)

    read_counts, fragment_sets = count_variant_hits_and_fragments(input_files_list, args.pid_cutoff)
    print(BColors.green(f"--- Aggregated results for {len(read_counts)} unique variants passing filters. ---"))

    summary_output_path = os.path.join(args.tmp_dir, f"{args.output_prefix}_varscan.tsv")
    calculate_and_write_summary(read_counts, fragment_sets, aro_to_family, aro_to_length, aro_to_annotation,
                                uscg_rpk, uscg_fpk, total_bases, total_reads, total_fragments, summary_output_path)

    print(BColors.green("\n\n--- Variant Tabulation and Normalization Complete ---"))

if __name__ == "__main__":
    main()
