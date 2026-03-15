#!/usr/bin/env python3

# resscan/scgscan_calculate_total_bases.py
import argparse
import os
import sys
import gzip
import bz2
import itertools
from concurrent.futures import ProcessPoolExecutor

class BColors:
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
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

def open_fastq_file(filepath):
    """Opens a FASTQ file, transparently handling .gz, .bz2, or uncompressed."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    elif filepath.endswith('.bz2'):
        return bz2.open(filepath, 'rt')
    else:
        return open(filepath, 'r')

def count_bases_in_file(filepath):
    """Worker function to count bases in a single file."""
    if not os.path.exists(filepath):
        return (filepath, 0, f"File not found")
    
    local_total = 0
    try:
        with open_fastq_file(filepath) as f:
            for line in itertools.islice(f, 1, None, 4):
                local_total += len(line.rstrip())
        return (filepath, local_total, None)
    except Exception as e:
        return (filepath, 0, str(e))

def main():
    parser = argparse.ArgumentParser(
        description="Calculates total bases from FASTQ files using multiple CPU cores."
    )
    parser.add_argument("-i", "--input-files", required=True, help="Comma-delimited list of FASTQ files.")
    parser.add_argument("-o", "--output-file", required=True, help="Path to output file.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of parallel processes (default: 4).")
    
    args = parser.parse_args()
    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]

    if not input_files_list:
        print(BColors.red("Error: No input files provided."), file=sys.stderr)
        sys.exit(1)

    print(BColors.cyan(f"--- Starting parallel base count with {args.threads} threads ---"))

    total_bases = 0
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(count_bases_in_file, input_files_list)
        
        for filepath, count, error in results:
            if error:
                print(BColors.red(f"Error processing {os.path.basename(filepath)}: {error}"), file=sys.stderr)
                sys.exit(1)
            else:
                print(BColors.green(f"  + {os.path.basename(filepath)}: {count:,} bases"))
                total_bases += count

    print(BColors.cyan(f"--- Total bases calculated: {total_bases:,} ---"))

    # Write output
    try:
        os.makedirs(os.path.dirname(os.path.abspath(args.output_file)), exist_ok=True)
        with open(args.output_file, 'w') as f_out:
            f_out.write(str(total_bases))
    except IOError as e:
        print(BColors.red(f"Error writing output: {e}"), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()