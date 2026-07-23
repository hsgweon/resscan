#!/usr/bin/env python3

# resscan/resscan_consolidate_all.py
import argparse
import os
import sys
from pathlib import Path
import pandas as pd

try:
    from resscan.utils import BColors
except ImportError:
    from utils import BColors

FINAL_HEADERS = {
    "homscan": [
        "AMR_Gene_Family", "ARO", "ARO_Name", "Drug_Class", "Resistance_Mechanism",
        "Read_Count", "Fragment_Count", "Lateral_Coverage_%",
        "RPK", "FPK", "RPKG", "FPKG", "RPKM", "FPKM",
        "RPKPC", "FPKPC", "RPKPMC", "FPKPMC", "RPKPGC", "FPKPGC",
        "Top_ARO", "Allocation_Proportions"
    ],
    "homscan_detailed": [
        "AMR_Gene_Family", "ARO", "ARO_Name", "Drug_Class", "Resistance_Mechanism",
        "Read_Count", "Fragment_Count", "Lateral_Coverage_%", "Gene_Length_bp",
        "RPK", "FPK", "RPKG", "FPKG", "RPKM", "FPKM",
        "RPKPC", "FPKPC", "RPKPMC", "FPKPMC", "RPKPGC", "FPKPGC"
    ],
    "varscan": [
        "AMR_Gene_Family", "ARO", "ARO_Name", "Drug_Class", "Resistance_Mechanism",
        "Read_Count", "Fragment_Count",
        "RPK", "FPK", "RPKG", "FPKG", "RPKM", "FPKM",
        "RPKPC", "FPKPC", "RPKPMC", "FPKPMC", "RPKPGC", "FPKPGC"
    ],
}

def get_report_type(filename):
    if "homscan_detailed.tsv" in filename: return "homscan_detailed"
    if "homscan.tsv" in filename: return "homscan"
    if "varscan.tsv" in filename: return "varscan"
    return None

def reformat_and_write_summary(source_path, destination_path, report_type):
    header = FINAL_HEADERS.get(report_type)

    is_empty = not source_path.exists() or os.path.getsize(source_path) == 0

    try:
        if is_empty:
            with open(destination_path, 'w', newline='') as f_out:
                header_with_hash = list(header)
                header_with_hash[0] = "#" + header_with_hash[0]
                f_out.write('\t'.join(header_with_hash) + '\n')
            print(BColors.cyan(f"Created empty report (header only): {destination_path.name}"))
            return

        df = pd.read_csv(source_path, sep='\t')

        if df.empty:
            with open(destination_path, 'w', newline='') as f_out:
                header_with_hash = list(header)
                header_with_hash[0] = "#" + header_with_hash[0]
                f_out.write('\t'.join(header_with_hash) + '\n')
            return

        final_df = pd.DataFrame()
        for col in header:
            if col in df.columns:
                final_df[col] = df[col]
            else:
                final_df[col] = "N/A"

        with open(destination_path, 'w', newline='') as f_out:
            header_with_hash = list(final_df.columns)
            header_with_hash[0] = "#" + header_with_hash[0]
            f_out.write('\t'.join(header_with_hash) + '\n')
            final_df.to_csv(f_out, sep='\t', index=False, header=False, float_format='%.6f', na_rep='N/A')

        print(BColors.green(f"Successfully reformatted and copied: {destination_path.name}"))

    except Exception as e:
        print(BColors.red(f"Error processing file {source_path.name}: {e}"), file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description="Consolidates final summary reports.")
    parser.add_argument("--output-prefix", required=True, help="Full path to output prefix/directory")
    args = parser.parse_args()

    main_output_dir = Path(args.output_prefix)
    filename_prefix = main_output_dir.name

    print(BColors.cyan(f"--- Consolidating and reformatting final reports ---"))

    tasks = [
        (main_output_dir / 'tmp' / 'homscan' / f"{filename_prefix}_homscan.tsv", "homscan"),
        (main_output_dir / 'tmp' / 'homscan' / f"{filename_prefix}_homscan_detailed.tsv", "homscan_detailed"),
        (main_output_dir / 'tmp' / 'varscan' / f"{filename_prefix}_varscan.tsv", "varscan"),
    ]

    for source_path, report_type in tasks:
        destination_path = main_output_dir / source_path.name
        reformat_and_write_summary(source_path, destination_path, report_type)

    print(BColors.green("\n--- Consolidation Complete ---"))

if __name__ == "__main__":
    main()
