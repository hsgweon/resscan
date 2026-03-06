#!/usr/bin/env python3

# resscan/resscan_curate_metadata.py
import argparse
import os
import sys
import pandas as pd
import logging
import csv

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

def setup_logging(output_path):
    """Sets up logging."""
    # Create a log file in the same directory as the output file
    output_dir = os.path.dirname(os.path.abspath(output_path))
    log_file = os.path.join(output_dir, "metadata_curation.log")
    
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    # Also log to console if it's an error
    console = logging.StreamHandler()
    console.setLevel(logging.ERROR)
    logging.getLogger('').addHandler(console)
    return log_file

def load_rules(rules_path):
    """Loads curation rules from a CSV file."""
    print(BColors.cyan(f"--- Loading curation rules from: {rules_path} ---"))
    rules = []
    try:
        with open(rules_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            for row in reader:
                if not row or not row[0].strip() or row[0].startswith('#'):
                    continue
                
                if len(row) >= 2:
                    rules.append((row[0].strip(), row[1].strip()))
                else:
                    print(BColors.yellow(f"Warning: Skipping malformed rule line: {row}"))
    except IOError as e:
        print(BColors.red(f"Error reading rules file: {e}"), file=sys.stderr)
        return None
    
    print(BColors.green(f"--- Loaded {len(rules)} curation rules."))
    logging.info(f"Loaded {len(rules)} rules from {rules_path}")
    return rules

def apply_curation(df, rules):
    """
    Applies rules to the DataFrame.
    Matches 'Card_Short_Name' against the rules.
    Updates 'AMR_Gene_Family'.
    """
    changes_made = 0
        
    def curation_logic(row):
        target_name = str(row['Card_Short_Name'])
        current_family = row['AMR_Gene_Family']
        
        for pattern, new_family in rules:
            if pattern in target_name:
                if new_family != current_family:
                    return new_family, True # Changed
                return current_family, False # Matched but already same
        
        return current_family, False # No match

    # Apply logic
    print(BColors.cyan("--- Applying rules to metadata... ---"))
    
    # Apply returns a DataFrame/Series of tuples, we split them
    results = df.apply(curation_logic, axis=1, result_type='expand')
    
    # Update the dataframe
    df['AMR_Gene_Family'] = results[0]
    
    # Count True values in the second column of results
    changes_made = results[1].sum()

    return df, changes_made

def main():
    parser = argparse.ArgumentParser(
        description="Curates a ResScan metadata file by applying custom renaming rules to the AMR_Gene_Family column."
    )
    parser.add_argument(
        "-m", "--metadata",
        required=True,
        help="Path to the ResScan metadata file (e.g., resscan_DB_CARD_NR_metadata.txt)."
    )
    parser.add_argument(
        "-r", "--rules",
        required=True,
        help="Path to a CSV file with curation rules. Format: 'Name_Substring,New_Family'."
    )
    parser.add_argument(
        "-o", "--output",
        help="Path for the curated metadata file. If omitted, overwrites the input file."
    )
    args = parser.parse_args()

    # Determine output path
    output_path = args.output if args.output else args.metadata
    
    # Setup Log
    log_file = setup_logging(output_path)
    logging.info(f"Started curation. Input: {args.metadata}")

    # 1. Load Rules
    rules = load_rules(args.rules)
    if rules is None:
        sys.exit(1)

    # 2. Load Metadata
    print(BColors.cyan(f"--- Loading metadata from: {args.metadata} ---"))
    try:
        # Ensure we read as TSV
        df = pd.read_csv(args.metadata, sep='\t')
        
        # Validation: Check if columns exist
        if 'Card_Short_Name' not in df.columns or 'AMR_Gene_Family' not in df.columns:
            print(BColors.red("Error: Metadata file missing required columns 'Card_Short_Name' or 'AMR_Gene_Family'."), file=sys.stderr)
            sys.exit(1)
            
    except Exception as e:
        print(BColors.red(f"Error parsing metadata file: {e}"), file=sys.stderr)
        logging.error(f"Error parsing metadata: {e}")
        sys.exit(1)

    # 3. Apply Rules
    df, changes_count = apply_curation(df, rules)

    print(BColors.green(f"--- Curation complete. Updated {changes_count} entries. ---"))
    logging.info(f"Curation complete. Updated {changes_count} entries.")

    # 4. Save
    print(BColors.cyan(f"--- Saving curated metadata to: {output_path} ---"))
    try:
        df.to_csv(output_path, sep='\t', index=False)
        print(BColors.green("--- Successfully saved file."))
        print(BColors.green(f"--- Log written to: {log_file}"))
    except IOError as e:
        print(BColors.red(f"Error writing output file: {e}"), file=sys.stderr)
        logging.error(f"Error writing output: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()