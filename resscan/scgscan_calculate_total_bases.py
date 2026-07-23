#!/usr/bin/env python3

# resscan/scgscan_calculate_total_bases.py
import argparse
import os
import sys
import gzip
import bz2
import itertools
from concurrent.futures import ProcessPoolExecutor
try:
    from resscan.utils import BColors
except ImportError:
    from utils import BColors

def open_fastq_file(filepath):
    """Opens a FASTQ file, transparently handling .gz, .bz2, or uncompressed."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    elif filepath.endswith('.bz2'):
        return bz2.open(filepath, 'rt')
    else:
        return open(filepath, 'r')

def count_bases_in_file(filepath):
    """Worker function to count bases and reads in a single file."""
    if not os.path.exists(filepath):
        return (filepath, 0, 0, f"File not found")

    local_bases = 0
    local_reads = 0
    try:
        with open_fastq_file(filepath) as f:
            for line in itertools.islice(f, 1, None, 4):
                stripped = line.rstrip()
                local_bases += len(stripped)
                local_reads += 1
        return (filepath, local_bases, local_reads, None)
    except Exception as e:
        return (filepath, 0, 0, str(e))

def main():
    parser = argparse.ArgumentParser(
        description="Calculates total bases from FASTQ files using multiple CPU cores."
    )
    parser.add_argument("-i", "--input-files", required=True, help="Comma-delimited list of FASTQ files.")
    parser.add_argument("-o", "--output-file", required=True, help="Path to output file.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of parallel processes (default: 4).")
    parser.add_argument("--mate-count", type=int, default=1, help="Number of mates per fragment (e.g. 2 for paired-end). Default: 1")
    
    args = parser.parse_args()
    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]

    if not input_files_list:
        print(BColors.red("Error: No input files provided."), file=sys.stderr)
        sys.exit(1)

    print(BColors.cyan(f"--- Starting parallel base count with {args.threads} threads ---"))

    total_bases = 0
    total_reads = 0
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(count_bases_in_file, input_files_list)

        for filepath, bases, reads, error in results:
            if error:
                print(BColors.red(f"Error processing {os.path.basename(filepath)}: {error}"), file=sys.stderr)
                sys.exit(1)
            else:
                print(BColors.green(f"  + {os.path.basename(filepath)}: {bases:,} bases, {reads:,} reads"))
                total_bases += bases
                total_reads += reads

    total_fragments = total_reads // args.mate_count
    print(BColors.cyan(f"--- Total bases: {total_bases:,} | Total reads: {total_reads:,} | Total fragments: {total_fragments:,} (mate-count: {args.mate_count}) ---"))

    try:
        os.makedirs(os.path.dirname(os.path.abspath(args.output_file)), exist_ok=True)
        with open(args.output_file, 'w') as f_out:
            f_out.write(f"{total_bases}\n{total_reads}\n{total_fragments}\n")
    except IOError as e:
        print(BColors.red(f"Error writing output: {e}"), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()