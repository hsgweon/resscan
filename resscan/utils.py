#!/usr/bin/env python3

# resscan/utils.py
import csv
import os
import sys


class BColors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    OKCYAN = '\033[96m'
    GREEN = '\033[92m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    GREY = '\033[90m'
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
    @classmethod
    def blue(cls, text): return cls._colorize(cls.BLUE, text)


def load_metadata(metadata_path):
    """Loads metadata to map Sequence_ID to family, length, and annotation fields."""
    print(BColors.cyan(f"--- Loading metadata from: {metadata_path} ---"))
    if not os.path.exists(metadata_path):
        print(BColors.red(f"Error: Metadata file not found at '{metadata_path}'"), file=sys.stderr)
        sys.exit(1)
    aro_to_family, aro_to_length, aro_to_annotation = {}, {}, {}
    try:
        with open(metadata_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            if 'SeqNucLength' not in reader.fieldnames:
                print(BColors.red("Error: Metadata file must contain a 'SeqNucLength' column for normalisation."), file=sys.stderr)
                sys.exit(1)
            for row in reader:
                sid = row['Sequence_ID']
                aro_to_family[sid] = row['AMR_Gene_Family']
                try:
                    aro_to_length[sid] = int(row['SeqNucLength'])
                except (ValueError, TypeError):
                    aro_to_length[sid] = 0
                aro_to_annotation[sid] = {
                    'ARO_Name': row.get('ARO_Name', 'N/A'),
                    'Drug_Class': row.get('Drug_Class', 'N/A'),
                    'Resistance_Mechanism': row.get('Resistance_Mechanism', 'N/A'),
                }
    except (Exception, KeyError) as e:
        print(BColors.red(f"Error reading or parsing metadata file: {e}"), file=sys.stderr)
        sys.exit(1)
    print(BColors.green(f"--- Loaded metadata for {len(aro_to_family)} AROs."))
    return aro_to_family, aro_to_length, aro_to_annotation


def load_uscg_metrics(uscg_report_path):
    """Loads the Overall RPK and FPK values from the USCG report."""
    print(BColors.cyan(f"--- Loading USCG metrics from: {uscg_report_path} ---"))
    if not os.path.exists(uscg_report_path):
        print(BColors.red(f"Error: USCG report file not found at '{uscg_report_path}'"), file=sys.stderr)
        return None, None
    uscg_rpk, uscg_fpk = None, None
    try:
        with open(uscg_report_path, 'r') as f:
            for line in f:
                if line.startswith("Overall_RPK_Across_All_USCGs"):
                    uscg_rpk = float(line.strip().split('\t')[1])
                elif line.startswith("Overall_FPK_Across_All_USCGs"):
                    uscg_fpk = float(line.strip().split('\t')[1])
    except (IOError, ValueError, IndexError) as e:
        print(BColors.red(f"Error reading or parsing USCG report file: {e}"), file=sys.stderr)
        return None, None
    if uscg_rpk is None:
        print(BColors.yellow("Warning: 'Overall_RPK_Across_All_USCGs' not found in the USCG report."), file=sys.stderr)
    if uscg_fpk is None:
        print(BColors.yellow("Warning: 'Overall_FPK_Across_All_USCGs' not found in the USCG report. FPK-based normalization will be skipped."), file=sys.stderr)
    return uscg_rpk, uscg_fpk


def load_sequencing_metrics(filepath):
    """Loads total bases, total reads, and total fragments from the metrics file."""
    print(BColors.cyan(f"--- Loading sequencing metrics from: {filepath} ---"))
    if not os.path.exists(filepath):
        print(BColors.red(f"Error: Sequencing metrics file not found at '{filepath}'"), file=sys.stderr)
        sys.exit(1)
    try:
        with open(filepath, 'r') as f:
            lines = [l.strip() for l in f.readlines() if l.strip()]
        total_bases = int(lines[0])
        total_reads = int(lines[1]) if len(lines) > 1 else None
        total_fragments = int(lines[2]) if len(lines) > 2 else total_reads
        msg = f"--- Found total bases: {total_bases:,}"
        if total_reads is not None:
            msg += f" | total reads: {total_reads:,}"
        if total_fragments is not None:
            msg += f" | total fragments: {total_fragments:,}"
        print(BColors.green(msg + " ---"))
        return total_bases, total_reads, total_fragments
    except (IOError, ValueError, IndexError) as e:
        print(BColors.red(f"Error reading or parsing sequencing metrics file: {e}"), file=sys.stderr)
        sys.exit(1)
