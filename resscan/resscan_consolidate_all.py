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

FINAL_HEADERS = {
    "homscan": [
        "AMR_Gene_Family", "ARO", "CARD_Short_Name_Misc", "Read_Count", "Fragment_Count", 
        "Lateral_Coverage_%", "Gene_Length_bp", "RPK", "FPK", "RPKG", "FPKG", 
        "RPKPC", "FPKPC", "RPKPMC", "FPKPMC", "AMR_Gene_Family;ARO"
    ],
    "varscan": [
        "AMR_Gene_Family", "ARO", "CARD_Short_Name_Misc", "Read_Count", "Fragment_Count", 
        "Gene_Length_bp", "RPK", "FPK", "RPKG", "FPKG", "RPKPC", "FPKPC", "RPKPMC", "FPKPMC",
        "AMR_Gene_Family;ARO"
    ],
    "homscan_MAP": [
        "AMR_Gene_Family", "ARO", "CARD_Short_Name_Misc", "Read_Count", "Fragment_Count", 
        "Lateral_Coverage_%", "Gene_Length_bp", "RPK", "FPK", "RPKG", "FPKG", 
        "RPKPC", "FPKPC", "RPKPMC", "FPKPMC", "Top_ARO", "Allocation_Proportions",
        "AMR_Gene_Family;ARO"
    ]
}

def get_report_type(filename):
    if "homscan_MAP.tsv" in filename: return "homscan_MAP"
    if "homscan.tsv" in filename: return "homscan"
    if "varscan.tsv" in filename: return "varscan"
    return None

def parse_amr_key(key_string):
    match = re.search(r';(ARO_\d+.*)', key_string)
    if match:
        split_point = match.start()
        family = key_string[:split_point]
        aro_part = match.group(1)
    else:
        try:
            family, aro_part = key_string.split(';', 1)
        except ValueError:
            return key_string, "N/A", "N/A"

    if aro_part.startswith("multiple"):
        return family, "multiple", "N/A"
    
    if '__' in aro_part:
        aro_id_block, misc_block = aro_part.split('__', 1)
        return family, aro_id_block, misc_block
    else:
        return family, aro_part, "N/A"

def reformat_and_write_summary(source_path, destination_path, report_type):
    header = FINAL_HEADERS.get(report_type)
    
    if not source_path.exists() or os.path.getsize(source_path) == 0:
        try:
            with open(destination_path, 'w') as f:
                f.write('#' + '\t'.join(header) + '\n')
        except IOError: pass
        return

    try:
        df = pd.read_csv(source_path, sep='\t')
        
        if df.empty:
            with open(destination_path, 'w') as f:
                f.write('#' + '\t'.join(header) + '\n')
            return

        if '#AMR_Gene_Family;ARO' in df.columns:
            df.rename(columns={'#AMR_Gene_Family;ARO': 'AMR_Gene_Family;ARO'}, inplace=True)

        key_col = 'AMR_Gene_Family;ARO'
        if key_col in df.columns:
            df[['AMR_Gene_Family', 'ARO', 'CARD_Short_Name_Misc']] = df[key_col].apply(
                lambda x: pd.Series(parse_amr_key(x))
            )
        else:
             print(BColors.red(f"Error: Key column '{key_col}' not found in {source_path.name}."))
             return

        final_df = pd.DataFrame()
        for col in header:
            if col in df.columns:
                final_df[col] = df[col]
            else:
                final_df[col] = "N/A"
        
        with open(destination_path, 'w', newline='') as f_out:
            f_out.write('#' + '\t'.join(final_df.columns) + '\n')
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

    files_to_process = [
        ('homscan', f"{filename_prefix}_homscan.tsv"),
        ('homscan', f"{filename_prefix}_homscan_MAP.tsv"),
        ('varscan', f"{filename_prefix}_varscan.tsv")
    ]
    
    for sub_dir, filename in files_to_process:
        source_path = main_output_dir / 'tmp' / sub_dir / filename
        destination_path = main_output_dir / filename
        report_type = get_report_type(filename)
        
        if report_type:
            reformat_and_write_summary(source_path, destination_path, report_type)

    print(BColors.green("\n\n--- Consolidation Complete ---"))

if __name__ == "__main__":
    main()