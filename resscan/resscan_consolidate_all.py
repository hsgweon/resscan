#!/usr/bin/env python3

# resscan/resscan_consolidate_all.py
import argparse
import os
import sys
import re
from pathlib import Path
import pandas as pd

class BColors:
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    enabled = sys.stdout.isatty()
    @classmethod
    def _colorize(cls, c, t): return f"{c}{t}{cls.ENDC}" if cls.enabled else t
    @classmethod
    def cyan(cls, t): return cls._colorize(cls.OKCYAN, t)
    @classmethod
    def green(cls, t): return cls._colorize(cls.OKGREEN, t)
    @classmethod
    def red(cls, t): return cls._colorize(cls.FAIL, t)

# Updated Headers to match your new Tidy format
FINAL_HEADERS = {
    "homscan": [
        "AMR_Gene_Family", "ARO", "Read_Count", "Fragment_Count", 
        "Lateral_Coverage_%", "Gene_Length_bp", "RPK", "FPK", "RPKG", "FPKG", 
        "RPKPC", "FPKPC", "RPKPMC", "FPKPMC"
    ],
    "varscan": [
        "AMR_Gene_Family", "ARO", "Read_Count", "Fragment_Count", 
        "Gene_Length_bp", "RPK", "FPK", "RPKG", "FPKG", "RPKPC", "FPKPC", "RPKPMC", "FPKPMC"
    ],
    "homscan_MAP": [
        "AMR_Gene_Family", "ARO", "Read_Count", "Fragment_Count", 
        "Lateral_Coverage_%", "Gene_Length_bp", "RPK", "FPK", "RPKG", "FPKG", 
        "RPKPC", "FPKPC", "RPKPMC", "FPKPMC", "Top_ARO", "Allocation_Proportions"
    ]
}

def get_report_type(filename):
    if "homscan_MAP.tsv" in filename: return "homscan_MAP"
    if "homscan.tsv" in filename: return "homscan"
    if "varscan.tsv" in filename: return "varscan"
    return None

def reformat_and_write_summary(source_path, destination_path, report_type):
    header = FINAL_HEADERS.get(report_type)
    
    # Check if we should produce a header-only file due to missing/empty source
    is_empty = not source_path.exists() or os.path.getsize(source_path) == 0

    try:
        if is_empty:
            # Create a clean file with just the header
            with open(destination_path, 'w', newline='') as f_out:
                header_with_hash = list(header)
                header_with_hash[0] = "#" + header_with_hash[0]
                f_out.write('\t'.join(header_with_hash) + '\n')
            print(BColors.cyan(f"Created empty report (header only): {destination_path.name}"))
            return

        # If data exists, process normally
        df = pd.read_csv(source_path, sep='\t')
        
        if df.empty:
            with open(destination_path, 'w', newline='') as f_out:
                header_with_hash = list(header)
                header_with_hash[0] = "#" + header_with_hash[0]
                f_out.write('\t'.join(header_with_hash) + '\n')
            return

        # Create a clean final dataframe based on requested headers
        final_df = pd.DataFrame()
        for col in header:
            if col in df.columns:
                final_df[col] = df[col]
            else:
                final_df[col] = "N/A"
        
        # Write out with a # in the header for compatibility
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

    # Mapping source directories to their specific filenames
    tasks = [
        (main_output_dir / 'tmp' / 'homscan' / f"{filename_prefix}_homscan.tsv", "homscan"),
        (main_output_dir / 'tmp' / 'homscan' / f"{filename_prefix}_homscan_MAP.tsv", "homscan_MAP"),
        (main_output_dir / 'tmp' / 'varscan' / f"{filename_prefix}_varscan.tsv", "varscan")
    ]
    
    for source_path, report_type in tasks:
        destination_path = main_output_dir / source_path.name
        reformat_and_write_summary(source_path, destination_path, report_type)

    print(BColors.green("\n--- Consolidation Complete ---"))

if __name__ == "__main__":
    main()