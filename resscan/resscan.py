#!/usr/bin/env python3

# resscan/resscan.py
import argparse
import logging
import subprocess
import sys
import time
from pathlib import Path
import os
import platform
import re
import shutil
import glob

import importlib.metadata
from importlib import resources
try:
    from resscan.utils import BColors as Colors, describe_database
except ImportError:
    from utils import BColors as Colors, describe_database

__version__ = "1.2.1"

class ColoredFormatter(logging.Formatter):
    """Custom formatter to add colours to log messages for console output."""
    LOG_FORMAT = "%(message)s"
    FORMATS = {
        logging.DEBUG: Colors.GREY + LOG_FORMAT + Colors.ENDC,
        logging.INFO: Colors.GREEN + LOG_FORMAT + Colors.ENDC,
        logging.WARNING: Colors.WARNING + LOG_FORMAT + Colors.ENDC,
        logging.ERROR: Colors.FAIL + Colors.BOLD + LOG_FORMAT + Colors.ENDC,
        logging.CRITICAL: Colors.FAIL + Colors.BOLD + Colors.UNDERLINE + LOG_FORMAT + Colors.ENDC,
    }
    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)

def get_py_dependency_versions():
    """Retrieves the versions of specified Python packages."""
    dependencies = ['pip', 'setuptools', 'pandas', 'numpy']
    versions = {}
    for package in dependencies:
        try:
            versions[package] = importlib.metadata.version(package)
        except Exception:
            versions[package] = "Not Found"
    return versions

def get_tool_versions():
    """Retrieves the versions of external command-line tools."""
    versions = {}
    try:
        result = subprocess.run(['bwa'], capture_output=True, text=True, check=False)
        version_line = next((line for line in result.stderr.splitlines() if line.startswith('Version:')), None)
        versions['BWA'] = version_line.split(':', 1)[1].strip() if version_line else 'Unknown'
    except Exception:
        versions['BWA'] = 'Not Found'
    try:
        result = subprocess.run(['diamond', '--version'], capture_output=True, text=True, check=True)
        versions['Diamond'] = result.stdout.strip().split()[-1]
    except Exception:
        versions['Diamond'] = 'Not Found'
    return versions

def setup_logging(log_file):
    """Sets up logging to a file (plain) and console (coloured)."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers = []
    file_handler = logging.FileHandler(log_file, mode='w')
    file_formatter = logging.Formatter("%(asctime)s [%(levelname)s] - %(message)s")
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(ColoredFormatter())
    logger.addHandler(console_handler)

def log_boxed_header(title, color_class):
    """Logs a formatted, coloured, boxed section header."""
    padding = 2
    width = len(title) + padding * 2
    border_line = f"+{'-' * width}+"
    title_line = f"|{' ' * padding}{title.upper()}{' ' * padding}|"
    
    header_text = f"{border_line}\n{title_line}\n{border_line}"
    
    print(f"\n{color_class}{Colors.BOLD}{header_text}{Colors.ENDC}\n")
    
    logging.getLogger().handlers[0].handle(logging.makeLogRecord({
        'msg': f"\n--- {title.upper()} ---\n", 'levelname': 'INFO', 'levelno': logging.INFO,
        'pathname': '', 'lineno': 0, 'exc_info': None, 'func': None
    }))

def run_command(command):
    """Runs a command, logs its output in real-time, and returns success status."""
    command_str = ' '.join(map(str, command))
    print(f"{Colors.CYAN}Running Command:{Colors.ENDC} {command_str}")
    logging.getLogger().handlers[0].handle(logging.makeLogRecord({
        'msg': f"Running Command: {command_str}", 'levelname': 'INFO', 'levelno': logging.INFO,
        'pathname': '', 'lineno': 0, 'exc_info': None, 'func': None
    }))
    try:
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1, universal_newlines=True
        )
        with process.stdout:
            for line in iter(process.stdout.readline, ''):
                print(f"{Colors.GREY}\t{line.strip()}{Colors.ENDC}")
                logging.getLogger().handlers[0].handle(logging.makeLogRecord({
                    'msg': line.strip(), 'levelname': 'INFO', 'levelno': logging.INFO,
                    'pathname': '', 'lineno': 0, 'exc_info': None, 'func': None
                }))
        return_code = process.wait()
        if return_code != 0:
            logging.error(f"Command failed with exit code {return_code}: {command_str}")
            return False
    except FileNotFoundError:
        logging.error(f"Error: The command '{command[0]}' was not found.")
        return False
    except Exception as e:
        logging.error(f"An unexpected error occurred while running command: {e}")
        return False
    logging.info("Command finished successfully.")
    return True

def file_exists_and_is_not_empty(file_path):
    """Checks if a file exists and is not empty."""
    return file_path.is_file() and file_path.stat().st_size > 0

def format_duration(seconds):
    """Formats a duration in seconds into a human-readable string."""
    if seconds < 60:
        return f"{seconds:.2f} seconds"
    
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    
    parts = []
    if hours > 0:
        parts.append(f"{int(hours)} hours")
    if minutes > 0:
        parts.append(f"{int(minutes)} minutes")
    if seconds > 0:
        parts.append(f"{int(seconds)} seconds")
        
    return ", ".join(parts)

def get_data_path(filename):
    """Finds the absolute path to a data file packaged with the installation."""
    with resources.as_file(resources.files('resscan').joinpath(filename)) as path:
        return path

MATE_SUFFIX = re.compile(r'^(?P<stem>.+)[._-]R?(?P<mate>[12])(?:_001)?$')
FASTQ_EXTENSIONS = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')


def run_stem_of(path):
    """
    Run identity of a FASTQ file: its name with the mate marker removed.

    Files sharing a stem are the mates of one sequencing run. The stem is
    matched greedily so the marker binds to the last candidate in the name,
    otherwise a sample called something like CH_SR1_DSS matches on its own
    embedded '_R1'.
    """
    base = os.path.basename(str(path))
    for ext in FASTQ_EXTENSIONS:
        if base.endswith(ext):
            base = base[:-len(ext)]
            break
    match = MATE_SUFFIX.match(base)
    return match.group('stem') if match else base


def validate_run_structure(parsed_runs):
    """
    Rejects input that cannot describe real sequencing runs.

    A DNA fragment is read from one end or from both, so a run has one mate or
    two; there is no third mate. A larger count means several runs were passed
    as one, which ResScan cannot detect from the reads themselves: it would
    take the count at face value and derive the wrong number of fragments,
    silently distorting every fragment-normalised metric.
    """
    for index, run in enumerate(parsed_runs, start=1):
        if len(run) <= 2:
            continue

        stems = []
        for f in run:
            stem = run_stem_of(f)
            if stem not in stems:
                stems.append(stem)

        print(f"{Colors.FAIL}Error: run {index} lists {len(run)} files. A sequencing run has "
              f"one mate (single-end) or two (paired-end); there is no third mate.{Colors.ENDC}")
        for f in run:
            print(f"{Colors.FAIL}         {f}{Colors.ENDC}")

        if len(stems) > 1:
            print(f"{Colors.WARNING}\nThese look like {len(stems)} separate sequencing runs of the same "
                  f"sample, which must be separated by ';' rather than ','.{Colors.ENDC}")
            groups = []
            for stem in stems:
                groups.append(",".join(str(f) for f in run if run_stem_of(f) == stem))
            print(f"{Colors.WARNING}Try instead:{Colors.ENDC}\n  -i '" + ";".join(groups) + "'")
        else:
            print(f"{Colors.WARNING}\nAll {len(run)} files appear to belong to one run. If some are "
                  f"unpaired or singleton reads, they cannot be counted as mates of a pair; "
                  f"quantify them separately.{Colors.ENDC}")
        sys.exit(1)


def validate_inputs_exist(parsed_runs):
    """
    Checks every input file up front, reporting all problems at once.

    Without this the first sign of a mistyped path is BWA failing, after the
    reference index has been built; in a large batch that is a slow way to
    learn about a typo, and every remaining bad path has to be found one run
    at a time.
    """
    missing, empty = [], []
    for run in parsed_runs:
        for f in run:
            path = Path(f)
            if not path.is_file():
                missing.append(f)
            elif path.stat().st_size == 0:
                empty.append(f)

    if missing or empty:
        for f in missing:
            print(f"{Colors.FAIL}Error: input file not found: {f}{Colors.ENDC}")
        for f in empty:
            print(f"{Colors.FAIL}Error: input file is empty: {f}{Colors.ENDC}")
        print(f"{Colors.WARNING}\nCheck the paths above; they are resolved relative to the "
              f"directory you launched from.{Colors.ENDC}")
        sys.exit(1)


def main():
    start_time = time.monotonic()

    parser = argparse.ArgumentParser(
        description=f"resscan pipeline v{__version__}. A wrapper script to run the resscan analysis pipeline.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    # I/O Arguments
    parser.add_argument("-i", "--input-fastqs", required=True, type=str, help="Input FASTQ files. Use commas to separate mates within a run and semicolons to separate multiple runs. E.g. run1_R1.fq,run1_R2.fq;run2_R1.fq,run2_R2.fq")
    parser.add_argument("-o", "--output-prefix", required=True, help="Prefix for output files and the main output directory.")
    
    # Database Arguments
    parser.add_argument(
        "--card",
        dest="card_db_dir",
        required=True, 
        help="Directory containing the CARD database files: 'resscan_DB_CARD_NR_all_nuc.fasta' and 'resscan_DB_CARD_NR_metadata.txt'."
    )
    parser.add_argument("--db-scg", default=None, help="Path to a custom Single Copy Genes FASTA file (overrides packaged DB).")
    parser.add_argument("--db-scg-lengths", default=None, help="Path to a custom SCG gene lengths TSV file (overrides packaged DB).")
    
    # Performance & Filtering Arguments
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads to use (default: 8).")
    parser.add_argument("--homscan-pid-cutoff", type=float, default=0.95, help="Minimum percent identity for HOMSCAN hits (0.0-1.0 scale). Default: 0.95")
    parser.add_argument("--varscan-pid-cutoff", type=float, default=0.95, help="Minimum nucleotide percent identity for VARSCAN hits (0.0-1.0 scale). Default: 0.95")
    parser.add_argument("--consensus-cutoff", type=float, default=0.8, help="Consensus cutoff for homscan family assignment. Default: 0.80")
    parser.add_argument("--homscan-min-aln-len", type=int, default=40, help="Absolute floor (bp) on the aligned length of a HOMSCAN hit. Applied together with --homscan-min-aln-frac; the stricter governs. Reads shorter than this can never be assigned. Default: 40")
    parser.add_argument("--homscan-min-aln-frac", type=float, default=0.5, help="Minimum fraction of a read's own length that must align for a HOMSCAN hit. At 0.5 per-kilobase abundance is unbiased at any read length. Set to 0 to use --homscan-min-aln-len alone. Default: 0.5")
    
    # Gene Type & PID Type Arguments
    parser.add_argument("--homscan-gene-types", default='H', help="Comma-delimited list of gene types for homscan (e.g., 'H,K'). Default: 'H'")
    parser.add_argument("--varscan-gene-types", default='V,R', help="Comma-delimited list of variant types for varscan (e.g., 'V,R,O'). Default: 'V,R'")
    parser.add_argument("--homscan-pid-type", choices=['protein', 'nucleotide'], default='protein', help="PID type to use for homscan filtering. Default: protein")
    
    # MAP Resolver Arguments
    map_group = parser.add_argument_group('MAP Resolver Arguments', 'Options for the Maximum A Posteriori abundance resolver.')
    map_group.add_argument("--map-priors-file", default=None, help="Path to a tab-separated file of priors for the MAP resolver.")
    map_group.add_argument("--map-metric-column", default='RPKPMC', help="The numeric column from the homscan report to use for MAP abundance resolution. Default: RPKPMC. (MAP allocates ambiguous reads in proportion, so per-sample normalisations — the …G/…M/…PC/…PMC scalings — cancel; this choice effectively only selects read- vs fragment-based and per-kb vs raw counting.)")
    map_group.add_argument("--map-base-prior", type=float, default=1.0, help="Baseline prior 'pseudo-count' for genes NOT in the priors file. Default: 1.0")
    map_group.add_argument("--map-prior-strength", type=float, default=1.0, help="Multiplier for the influence of the priors file. Default: 1.0")

    # Control Flow Arguments
    parser.add_argument("--skip-mapping", action="store_true", help="Skip mapping and start from homscan/varscan.")
    parser.add_argument("--skip-to-tabulate", action="store_true", help="Skip mapping and SAM processing; re-run tabulation (including MAP) and consolidation only. Useful for re-running with different MAP priors.")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode, creating detailed logs for sub-scripts.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output directory.")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args()

    # Setup Paths
    runs = [r.strip() for r in args.input_fastqs.split(';') if r.strip()]
    parsed_runs = [[m.strip() for m in run.split(',') if m.strip()] for run in runs]
    validate_run_structure(parsed_runs)
    validate_inputs_exist(parsed_runs)
    mate_counts = [len(r) for r in parsed_runs]
    if len(set(mate_counts)) > 1:
        print(f"{Colors.FAIL}Error: All runs must have the same number of mates. Found: {mate_counts}{Colors.ENDC}")
        sys.exit(1)
    mate_count = mate_counts[0] if mate_counts else 1
    fastq_files = [Path(f) for run in parsed_runs for f in run]
    input_fastqs_str = ",".join(map(str, fastq_files))
    output_prefix_path = Path(args.output_prefix).resolve()
    output_prefix_name = output_prefix_path.name
    main_output_dir = output_prefix_path

    # Resolve Database Paths from the user-provided directory
    card_db_dir = Path(args.card_db_dir).resolve()
    if not card_db_dir.is_dir():
        print(f"{Colors.FAIL}Error: The provided CARD database directory does not exist: {card_db_dir}{Colors.ENDC}")
        sys.exit(1)

    db_card = card_db_dir / "resscan_DB_CARD_NR_all_nuc.fasta"
    db_card_metadata = card_db_dir / "resscan_DB_CARD_NR_metadata.txt"

    if not db_card.is_file():
        print(f"{Colors.FAIL}Error: Required CARD FASTA file not found at: {db_card}{Colors.ENDC}")
        print(f"{Colors.FAIL}Please ensure the directory contains 'resscan_DB_CARD_NR_all_nuc.fasta'{Colors.ENDC}")
        sys.exit(1)
    
    if not db_card_metadata.is_file():
        print(f"{Colors.FAIL}Error: Required CARD metadata file not found at: {db_card_metadata}{Colors.ENDC}")
        print(f"{Colors.FAIL}Please ensure the directory contains 'resscan_DB_CARD_NR_metadata.txt'{Colors.ENDC}")
        sys.exit(1)

    # Resolve SCG Database (Packaged or Custom)
    db_scg = args.db_scg or get_data_path('databases/resscan_DB_SCG/SCG.fasta')
    db_scg_lengths = args.db_scg_lengths or get_data_path('databases/resscan_DB_SCG/SCG_gene_lengths.tsv')

    # Handle Existing Directory and --overwrite Flag
    if main_output_dir.exists():
        if not args.overwrite:
            print(f"{Colors.FAIL}Error: Output directory '{main_output_dir}' already exists.{Colors.ENDC}")
            print(f"{Colors.FAIL}Please remove it or use the --overwrite flag to proceed.{Colors.ENDC}")
            sys.exit(1)
        else:
            if args.skip_to_tabulate:
                print(f"{Colors.WARNING}Warning: --overwrite and --skip-to-tabulate are set.{Colors.ENDC}")
                print(f"{Colors.WARNING}Cleaning up previous tabulation results...{Colors.ENDC}")
                files_to_remove = [
                    main_output_dir / "tmp" / "homscan" / f"{output_prefix_name}_homscan.tsv",
                    main_output_dir / "tmp" / "homscan" / f"{output_prefix_name}_homscan_detailed.tsv",
                    main_output_dir / "tmp" / "varscan" / f"{output_prefix_name}_varscan.tsv",
                    main_output_dir / f"{output_prefix_name}_homscan.tsv",
                    main_output_dir / f"{output_prefix_name}_homscan_detailed.tsv",
                    main_output_dir / f"{output_prefix_name}_varscan.tsv",
                ]
                for f in files_to_remove:
                    if f.exists(): os.remove(f)
            elif args.skip_mapping and not args.skip_to_tabulate:
                print(f"{Colors.WARNING}Warning: --overwrite and --skip-mapping are set.{Colors.ENDC}")
                print(f"{Colors.WARNING}Cleaning up previous analysis results but preserving mapping data...{Colors.ENDC}")
                dirs_to_remove = [
                    main_output_dir / "tmp" / "homscan", main_output_dir / "tmp" / "varscan",
                    main_output_dir / "logs", main_output_dir / f"{output_prefix_name}_homscan_html",
                    main_output_dir / f"{output_prefix_name}_varscan_html"
                ]
                for d in dirs_to_remove:
                    if d.exists(): shutil.rmtree(d)
                for f in glob.glob(str(main_output_dir / f"{output_prefix_name}_*.tsv")):
                    os.remove(f)
            else:
                print(f"{Colors.WARNING}Warning: --overwrite flag is set. The entire directory '{main_output_dir}' will be removed.{Colors.ENDC}")
                try:
                    shutil.rmtree(main_output_dir)
                except OSError as e:
                    print(f"{Colors.FAIL}Error: Could not remove directory {main_output_dir}. Reason: {e}{Colors.ENDC}")
                    sys.exit(1)

    # Setup Logging and Directories
    log_dir = main_output_dir / "logs"
    setup_logging(log_dir / f"{output_prefix_name}_run.log")
    
    # Write Run Details
    run_details_file = log_dir / f"{output_prefix_name}_run_details.txt"
    print(f"{Colors.BLUE}{Colors.BOLD}--- resscan pipeline starting (v{__version__}) ---{Colors.ENDC}")
    db_info = describe_database(card_db_dir, fasta_path=db_card, metadata_path=db_card_metadata)

    with open(run_details_file, "w") as f:
        f.write(f"ResScan Pipeline Run Details\n" + "-"*30 + "\n")
        f.write(f"Pipeline Version: {__version__}\n")
        f.write(f"Date and Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Command-line: {' '.join(sys.argv)}\n\n")
        f.write("--- System & Tool Information ---\n")
        f.write(f"Python Version: {platform.python_version()}\n")
        for tool, version in get_tool_versions().items(): f.write(f"{tool} Version: {version}\n")
        for dep, version in get_py_dependency_versions().items(): f.write(f"{dep} Version: {version}\n")
        f.write(f"\n--- Databases Used ---\n")
        f.write(f"CARD FASTA: {db_card}\n")
        f.write(f"CARD Metadata: {db_card_metadata}\n")
        f.write(f"SCG FASTA: {db_scg}\n")
        for key, value in db_info.items():
            f.write(f"{key}: {value}\n")

    # A compact provenance record, written beside the result tables rather than
    # inside logs/, so that it travels with the results when they are shared.
    provenance_file = main_output_dir / f"{output_prefix_name}_provenance.txt"
    with open(provenance_file, "w") as f:
        f.write("# ResScan run provenance. Records the software and reference database\n")
        f.write("# that produced the result tables in this directory.\n")
        f.write(f"ResScan_Version\t{__version__}\n")
        f.write(f"Run_Date\t{time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Command_Line\t{' '.join(sys.argv)}\n")
        f.write(f"Python_Version\t{platform.python_version()}\n")
        for tool, version in get_tool_versions().items():
            f.write(f"{tool}_Version\t{version}\n")
        for key, value in db_info.items():
            f.write(f"Database_{key}\t{value}\n")
        f.write(f"SCG_Database\t{db_scg}\n")
        f.write(f"Param_homscan_min_aln_len\t{args.homscan_min_aln_len}\n")
        f.write(f"Param_homscan_min_aln_frac\t{args.homscan_min_aln_frac}\n")
        f.write(f"Param_homscan_pid_cutoff\t{args.homscan_pid_cutoff}\n")
        f.write(f"Param_homscan_pid_type\t{args.homscan_pid_type}\n")
        f.write(f"Param_varscan_pid_cutoff\t{args.varscan_pid_cutoff}\n")
        f.write(f"Param_consensus_cutoff\t{args.consensus_cutoff}\n")
        f.write(f"Param_map_metric_column\t{args.map_metric_column}\n")

    # Define all sub-directories and script paths
    src_dir = Path(__file__).parent.resolve()
    tmp_dir = main_output_dir / "tmp"
    resscan_tmp = tmp_dir / "resscan"
    scgscan_tmp = tmp_dir / "scgscan"
    homscan_tmp = tmp_dir / "homscan"
    varscan_tmp = tmp_dir / "varscan"
    homscan_html_dir = main_output_dir / f"{output_prefix_name}_homscan_html"
    varscan_html_dir = main_output_dir / f"{output_prefix_name}_varscan_html"

    for d in [resscan_tmp, scgscan_tmp, homscan_tmp, varscan_tmp, homscan_html_dir, varscan_html_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # ==========================================================================
    # MAPPING
    # ==========================================================================
    if not args.skip_mapping and not args.skip_to_tabulate:
        log_boxed_header("Mapping", Colors.HEADER)
        # BWA map for CARD
        amr_map_cmd = [sys.executable, str(src_dir / "resscan_map_reads_bwa.py"), "-i", input_fastqs_str, "--db", str(db_card), "--tmp-dir", str(resscan_tmp), "-t", str(args.threads)]
        if not run_command(amr_map_cmd): sys.exit(1)

        # DIAMOND map for SCG
        scg_map_cmd = [sys.executable, str(src_dir / "scgscan_map_reads_diamond.py"), "-i", input_fastqs_str, "--db", str(db_scg), "--tmp-dir", str(scgscan_tmp), "-t", str(args.threads)]
        if not run_command(scg_map_cmd): sys.exit(1)

        diamond_outputs = [scgscan_tmp / f"{f.name}.diamond.tsv" for f in fastq_files]
        uscg_report = scgscan_tmp / f"{output_prefix_name}_uscg_report.tsv"
        uscg_table = scgscan_tmp / f"{output_prefix_name}_uscg.tsv"
        
        valid_diamond_outputs = [str(p) for p in diamond_outputs if file_exists_and_is_not_empty(p)]

        if valid_diamond_outputs:
            scg_quant_cmd = [sys.executable, str(src_dir / "scgscan_quantify_from_diamond.py"), "-i", ",".join(valid_diamond_outputs), "-l", str(db_scg_lengths), "-o", str(uscg_report), "--output-table", str(uscg_table), "-e", "1e-5"]
            if not run_command(scg_quant_cmd): sys.exit(1)
        else:
            logging.warning("No valid DIAMOND hits found in any files. Skipping SCG quantification.")
            uscg_report.touch()

        total_bases_file = scgscan_tmp / f"{output_prefix_name}_total_bases.txt"
        calc_bases_cmd = [sys.executable, str(src_dir / "scgscan_calculate_total_bases.py"), "-i", input_fastqs_str, "-o", str(total_bases_file), "-t", str(args.threads), "--mate-count", str(mate_count), "--min-aln-len", str(args.homscan_min_aln_len), "--min-aln-frac", str(args.homscan_min_aln_frac)]
        if not run_command(calc_bases_cmd): sys.exit(1)
    else:
        if args.skip_to_tabulate:
            print(f"{Colors.WARNING}Skipping mapping and SAM processing as requested by --skip-to-tabulate.{Colors.ENDC}")
            total_bases_file = scgscan_tmp / f"{output_prefix_name}_total_bases.txt"
            uscg_report = scgscan_tmp / f"{output_prefix_name}_uscg_report.tsv"
            if not total_bases_file.is_file() or not uscg_report.is_file():
                logging.error("--skip-to-tabulate requires total_bases_file and uscg_report from a previous mapping run.")
                sys.exit(1)
        else:
            print(f"{Colors.WARNING}Skipping mapping steps as requested by --skip-mapping.{Colors.ENDC}")
            total_bases_file = scgscan_tmp / f"{output_prefix_name}_total_bases.txt"
            uscg_report = scgscan_tmp / f"{output_prefix_name}_uscg_report.tsv"
            if not total_bases_file.is_file() or not uscg_report.is_file():
                logging.error("--skip-mapping was used, but required files from the mapping step are missing.")
                sys.exit(1)

    # Define dynamic file paths for downstream steps
    sam_files = [resscan_tmp / f"{f.name}.sam" for f in fastq_files]
    sam_files_str = ",".join(map(str, sam_files))
    homscan_hits_files = [homscan_tmp / f"{s.name}_hits.tsv" for s in sam_files]
    homscan_hits_files_str = ",".join(map(str, homscan_hits_files))

    # ==========================================================================
    # HOMSCAN
    # ==========================================================================
    log_boxed_header("HomScan", Colors.HEADER)

    if not args.skip_to_tabulate:
        valid_sam_files = [str(p) for p in sam_files if file_exists_and_is_not_empty(p)]
        if valid_sam_files:
            homscan_process_cmd = [
                sys.executable, str(src_dir / "homscan_process_sam.py"), "-i", ",".join(valid_sam_files), "--db", str(db_card),
                "--tmp-dir", str(homscan_tmp), "--output-prefix", output_prefix_name,
                "--gene-types", args.homscan_gene_types,
                "--pid-cutoff", str(args.homscan_pid_cutoff),
                "--pid-type", args.homscan_pid_type,
                "--min-aln-len", str(args.homscan_min_aln_len),
                "--min-aln-frac", str(args.homscan_min_aln_frac)
            ]
            if args.debug: homscan_process_cmd.append("--debug")
            if not run_command(homscan_process_cmd): sys.exit(1)
        else:
            logging.warning("No non-empty SAM files found. Skipping homscan process step.")

    valid_hits_files = [str(p) for p in homscan_hits_files if file_exists_and_is_not_empty(p)]
    if valid_hits_files:
        homscan_tabulate_cmd = [
            sys.executable, str(src_dir / "homscan_tabulate_and_normalise.py"), "-i", ",".join(valid_hits_files),
            "--metadata", str(db_card_metadata), "--uscg-report", str(uscg_report),
            "--total-bases-file", str(total_bases_file), "--pid-cutoff", str(args.homscan_pid_cutoff),
            "--consensus", str(args.consensus_cutoff), "--tmp-dir", str(homscan_tmp),
            "--output-prefix", output_prefix_name,
            "--pid-type", args.homscan_pid_type,
            "--map-metric-column", args.map_metric_column,
            "--map-base-prior", str(args.map_base_prior),
            "--map-prior-strength", str(args.map_prior_strength),
        ]
        if args.map_priors_file:
            homscan_tabulate_cmd.extend(["--map-priors-file", args.map_priors_file])
        if not run_command(homscan_tabulate_cmd): sys.exit(1)

        if not args.skip_to_tabulate:
            homscan_visualise_cmd = [
                sys.executable, str(src_dir / "homscan_visualise.py"), "-i", ",".join(valid_hits_files),
                "--metadata", str(db_card_metadata), "--db", str(db_card),
                "-o", str(homscan_html_dir), "--pid-cutoff", str(args.homscan_pid_cutoff),
                "--pid-type", args.homscan_pid_type
            ]
            if not run_command(homscan_visualise_cmd): sys.exit(1)
    else:
        logging.warning("Homscan hits not found or empty. Skipping homscan tabulation.")


    # ==========================================================================
    # VARSCAN
    # ==========================================================================
    log_boxed_header("VarScan", Colors.HEADER)

    if not args.skip_to_tabulate:
        valid_sam_files = [str(p) for p in sam_files if file_exists_and_is_not_empty(p)]
        if valid_sam_files:
            varscan_process_cmd = [
                sys.executable, str(src_dir / "varscan_process_sam.py"), "-i", ",".join(valid_sam_files), "--db", str(db_card),
                "--tmp-dir", str(varscan_tmp), "--output-prefix", output_prefix_name,
                "--gene-types", args.varscan_gene_types
            ]
            if args.debug: varscan_process_cmd.append("--debug")
            if not run_command(varscan_process_cmd): sys.exit(1)
        else:
            logging.warning("No non-empty SAM files found. Skipping varscan process step.")

    variant_hits = varscan_tmp / f"{output_prefix_name}_variant_hits.tsv"
    if file_exists_and_is_not_empty(variant_hits):
        varscan_tabulate_cmd = [
            sys.executable, str(src_dir / "varscan_tabulate_and_normalise.py"), "-i", str(variant_hits),
            "--metadata", str(db_card_metadata), "--uscg-report", str(uscg_report),
            "--total-bases-file", str(total_bases_file), "--tmp-dir", str(varscan_tmp),
            "--output-prefix", output_prefix_name, "--pid-cutoff", str(args.varscan_pid_cutoff)
        ]
        if not run_command(varscan_tabulate_cmd): sys.exit(1)

        if not args.skip_to_tabulate:
            variant_alignments = varscan_tmp / f"{output_prefix_name}_variant_alignments.txt"
            if file_exists_and_is_not_empty(variant_alignments):
                varscan_visualise_cmd = [
                    sys.executable, str(src_dir / "varscan_visualise.py"), "--variant-hits", str(variant_hits),
                    "--variant-alignments", str(variant_alignments), "--metadata", str(db_card_metadata),
                    "-o", str(varscan_html_dir),
                    "--pid-cutoff", str(args.varscan_pid_cutoff)
                ]
                if not run_command(varscan_visualise_cmd): sys.exit(1)
            else:
                logging.warning("Variant alignments not found. Skipping varscan visualisation.")
    else:
        logging.warning("Variant hits not found. Skipping varscan tabulation.")


    # ==========================================================================
    # FINAL
    # ==========================================================================
    log_boxed_header("Final Report Consolidation", Colors.BLUE)
    final_cmd = [sys.executable, str(src_dir / "resscan_consolidate_all.py"), "--output-prefix", str(output_prefix_path)]
    if not run_command(final_cmd): sys.exit(1)

    end_time = time.monotonic()
    duration_seconds = end_time - start_time
    
    print(f"\n{Colors.BLUE}{Colors.BOLD}--- resscan pipeline finished successfully! ---{Colors.ENDC}")
    print(f"{Colors.BLUE}Total execution time: {format_duration(duration_seconds)}{Colors.ENDC}")

    print(f"\n{Colors.BOLD}--- Key Output Files ---{Colors.ENDC}")

    def summarize_files(title, directory, pattern):
        files = sorted(glob.glob(str(directory / pattern)))
        if files:
            print(f"\n{Colors.GREEN}{title}:{Colors.ENDC}")
            for f in files:
                print(f"  - {f}")
    
    summarize_files("Final Summary Tables", main_output_dir, f"{output_prefix_name}_*.tsv")
    if not args.skip_to_tabulate:
        summarize_files("HomScan Visualisations", homscan_html_dir, "*.html")
        summarize_files("VarScan Visualisations", varscan_html_dir, "*.html")
    summarize_files("Log Files", log_dir, "*")


if __name__ == "__main__":
    main()
