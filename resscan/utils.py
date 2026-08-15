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

    # A zero denominator is legitimate for samples that contain little or no
    # cellular material (extracellular DNA, filtered water, a blank). The
    # per-cell metrics are undefined rather than zero in that case, so they are
    # reported as NA. Say so plainly: otherwise a whole column of NA looks like
    # a pipeline failure.
    if (uscg_rpk is None or uscg_rpk <= 1e-9) and (uscg_fpk is None or uscg_fpk <= 1e-9):
        print(BColors.yellow(
            "Warning: no single-copy marker genes were detected in this sample, so the "
            "per-cell denominator is zero. Gene-length- and depth-normalised metrics "
            "(RPK, FPK, RPKG, FPKG, RPKM, FPKM) are unaffected, but all per-cell metrics "
            "(RPKPC/FPKPC, RPKPMC/FPKPMC, RPKPGC/FPKPGC) will be reported as NA. This is "
            "expected for samples with little or no cellular biomass."), file=sys.stderr)

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


def load_effective_length_offset(filepath):
    """
    Loads the mean effective-length offset (in bp) from the metrics file.

    This is the fourth line, written from v1.1.0 onwards. Metrics files produced
    by earlier versions have no such line; zero is returned in that case, which
    leaves per-kilobase abundance unchanged.
    """
    try:
        with open(filepath, 'r') as f:
            lines = [l.strip() for l in f.readlines() if l.strip()]
        if len(lines) > 3:
            return float(lines[3])
    except (IOError, ValueError, IndexError):
        pass
    return 0.0


DB_INFO_FILENAME = "resscan_DB_CARD_INFO.txt"


def md5_of_file(path, chunk_size=1 << 20):
    """MD5 checksum of a file, or None if it cannot be read."""
    import hashlib
    try:
        h = hashlib.md5()
        with open(path, 'rb') as f:
            for chunk in iter(lambda: f.read(chunk_size), b''):
                h.update(chunk)
        return h.hexdigest()
    except (IOError, OSError):
        return None


def read_card_version(input_dir):
    """
    Best-effort CARD release version, read from card.json in the raw download.

    CARD ships a '_version' key in that file. It is absent from older releases
    and from downloads that omit card.json, so a missing value is not an error.
    """
    import json
    path = os.path.join(input_dir, "card.json")
    if not os.path.exists(path):
        return None
    try:
        with open(path, 'r', encoding='utf-8', errors='replace') as f:
            data = json.load(f)
        if isinstance(data, dict):
            for key in ("_version", "version", "_CARD_version"):
                if data.get(key):
                    return str(data[key])
    except (IOError, ValueError):
        pass
    return None


def describe_database(db_dir, fasta_path=None, metadata_path=None):
    """
    Identifies the reference database in use.

    Prefers the manifest written at build time, which carries the CARD release
    and sequence counts. Databases built before manifests existed have none, so
    the files are checksummed directly instead: a database is then still
    identified unambiguously, just without a human-readable version.
    """
    info = {}
    manifest = os.path.join(str(db_dir), DB_INFO_FILENAME)
    if os.path.exists(manifest):
        try:
            with open(manifest, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#') or '\t' not in line:
                        continue
                    key, _, value = line.partition('\t')
                    info[key.strip()] = value.strip()
            info["Manifest"] = "present"
        except IOError:
            info["Manifest"] = "unreadable"
    else:
        info["Manifest"] = "absent (database predates build manifests; identified by checksum)"

    info.setdefault("Directory", str(db_dir))
    for label, path in (("Reference_FASTA", fasta_path), ("Reference_Metadata", metadata_path)):
        if not path:
            continue
        key = f"{label}_MD5"
        if key not in info:
            digest = md5_of_file(str(path))
            if digest:
                info[key] = digest
    return info
