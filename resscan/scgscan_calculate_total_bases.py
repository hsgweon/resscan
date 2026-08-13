#!/usr/bin/env python3

# resscan/scgscan_calculate_total_bases.py
import argparse
import os
import sys
import gzip
import bz2
import math
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

def count_bases_in_file(task):
    """
    Worker function to count bases and reads in a single file.

    Also accumulates the effective-length offset. A read of length r that must
    present m aligned bases gathers evidence for a gene over a window of
    (L - 2m + r), so it contributes (r - 2m) to the difference between the
    window and the gene length. Summed over eligible reads this gives the
    per-sample offset used to correct per-kilobase abundance. At a fraction of
    0.5, m = r/2 and every term is zero.
    """
    filepath, min_aln_len, min_aln_frac = task
    if not os.path.exists(filepath):
        return (filepath, 0, 0, 0, 0, "File not found")

    local_bases = 0
    local_reads = 0
    local_offset = 0
    local_eligible = 0
    try:
        with open_fastq_file(filepath) as f:
            for line in itertools.islice(f, 1, None, 4):
                read_length = len(line.rstrip())
                local_bases += read_length
                local_reads += 1

                required = 0
                if min_aln_frac > 0 and read_length > 0:
                    required = math.ceil(read_length * min_aln_frac)
                required = max(int(min_aln_len), required)

                if read_length >= required:
                    local_offset += read_length - 2 * required
                    local_eligible += 1
        return (filepath, local_bases, local_reads, local_offset, local_eligible, None)
    except Exception as e:
        return (filepath, 0, 0, 0, 0, str(e))

def main():
    parser = argparse.ArgumentParser(
        description="Calculates total bases from FASTQ files using multiple CPU cores."
    )
    parser.add_argument("-i", "--input-files", required=True, help="Comma-delimited list of FASTQ files.")
    parser.add_argument("-o", "--output-file", required=True, help="Path to output file.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of parallel processes (default: 4).")
    parser.add_argument("--mate-count", type=int, default=1, help="Number of mates per fragment (e.g. 2 for paired-end). Default: 1")
    parser.add_argument("--min-aln-len", type=int, default=40, help="Alignment length floor used by HomScan, needed for the effective-length offset. Default: 40")
    parser.add_argument("--min-aln-frac", type=float, default=0.5, help="Alignment length fraction used by HomScan, needed for the effective-length offset. Default: 0.5")

    args = parser.parse_args()
    input_files_list = [f.strip() for f in args.input_files.split(',') if f.strip()]

    if not input_files_list:
        print(BColors.red("Error: No input files provided."), file=sys.stderr)
        sys.exit(1)

    print(BColors.cyan(f"--- Starting parallel base count with {args.threads} threads ---"))

    total_bases = 0
    total_reads = 0
    total_offset = 0
    total_eligible = 0
    tasks = [(f, args.min_aln_len, args.min_aln_frac) for f in input_files_list]
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        results = executor.map(count_bases_in_file, tasks)

        for filepath, bases, reads, offset, eligible, error in results:
            if error:
                print(BColors.red(f"Error processing {os.path.basename(filepath)}: {error}"), file=sys.stderr)
                sys.exit(1)
            else:
                print(BColors.green(f"  + {os.path.basename(filepath)}: {bases:,} bases, {reads:,} reads"))
                total_bases += bases
                total_reads += reads
                total_offset += offset
                total_eligible += eligible

    total_fragments = total_reads // args.mate_count
    print(BColors.cyan(f"--- Total bases: {total_bases:,} | Total reads: {total_reads:,} | Total fragments: {total_fragments:,} (mate-count: {args.mate_count}) ---"))

    mean_offset = (total_offset / total_eligible) if total_eligible > 0 else 0.0
    if total_eligible < total_reads:
        print(BColors.yellow(f"--- Warning: {total_reads - total_eligible:,} of {total_reads:,} reads are shorter than the "
                             f"alignment length floor ({args.min_aln_len} bp) and can never be assigned to a gene. ---"))
    print(BColors.cyan(f"--- Effective-length offset: {mean_offset:+.2f} bp "
                       f"(0.00 means per-kilobase abundance is unbiased) ---"))

    try:
        os.makedirs(os.path.dirname(os.path.abspath(args.output_file)), exist_ok=True)
        with open(args.output_file, 'w') as f_out:
            f_out.write(f"{total_bases}\n{total_reads}\n{total_fragments}\n{mean_offset}\n{total_eligible}\n")
    except IOError as e:
        print(BColors.red(f"Error writing output: {e}"), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()