#!/usr/bin/env python3

# resscan/resscan_build_db.py
import argparse
import os
import sys
import shutil
import re
import csv
import logging
import datetime
from collections import defaultdict

# --- BColors Class ---
class BColors:
    """A helper class to add color to terminal output."""
    HEADER = '\033[95m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
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

# --- Global Constants ---
class ResScanFiles:
    COMMON_DB_PREFIX = "resscan_DB_CARD"
    
    EXPECTED_INPUTS = {
        "snps": "snps.txt",
        "homolog": "nucleotide_fasta_protein_homolog_model.fasta",
        "knockout": "nucleotide_fasta_protein_knockout_model.fasta",
        "variant": "nucleotide_fasta_protein_variant_model.fasta",
        "overexpression": "nucleotide_fasta_protein_overexpression_model.fasta",
        "rrna": "nucleotide_fasta_rRNA_gene_variant_model.fasta",
        "aro_index": "aro_index.tsv"
    }

    class Stage1:
        README_SUFFIX = "README"
        FASTA_HOMOLOG = "homolog_model"
        FASTA_KNOCKOUT = "knockout_model"
        FASTA_VARIANT = "protein_variant_model"
        FASTA_OVEREXPRESSION = "protein_overexpression_model"
        FASTA_RRNA = "rRNA_gene_variant_model"
        FASTA_COMBINED = "all"
        METADATA = "metadata"
        METADATA_NON_VARIANT_IGNORED = "metadata_ignored_non_variant"
        METADATA_VARIANT_IGNORED = "metadata_ignored_variant"
    class Stage2:
        NR_METADATA = "NR_metadata"
        NR_FASTA = "NR_all"
        NR_FASTA_NUC = "NR_all_nuc" # Nucleotide - Required by ResScan

# --- Data Structures ---
class Config:
    def __init__(self):
        self.path_snps_file = None
        self.path_protein_homolog_model = None
        self.path_protein_knockout_model = None
        self.path_protein_variant_model = None
        self.path_protein_overexpression_model = None
        self.path_rRNA_gene_variant_model = None
        self.path_aro_index = None
        self.output_dir_name = None
        self.overwrite_outputs = False
        self.keep_intermediate = False 

class FastaEntry:
    def __init__(self, header="", sequence=""):
        self.header = header
        self.sequence = sequence

class ProteinFastaHeaderInfo:
    def __init__(self):
        self.aro_number = "N/A"
        self.name = "N/A"
        self.card_short_name = "N/A"
        self.parsed_successfully = False

class SNPSFileLine:
    def __init__(self):
        self.accession_aro = ""
        self.name = ""
        self.model_type = ""
        self.parameter_type = ""
        self.mutations_str = ""
        self.card_short_name = ""
        self.source = ""
        self.citation = ""
        self.original_line = ""
        self.line_number = 0
        self.is_protein_pathway_candidate = False
        self.is_rrna_pathway_candidate = False

class ParsedProteinMutation:
    def __init__(self):
        self.original_aa = ''
        self.new_aa = ''
        self.position = 0
        self.is_frameshift = False
        self.is_nonsense = False
        self.original_token = ''

class ParsedNucleotideMutation:
    def __init__(self):
        self.original_base = ''
        self.new_base = ''
        self.position = 0
        self.original_token = ''

class MetadataEntry:
    def __init__(self, line=None, is_header=False):
        self.sequence_id = ""
        self.aro_number = ""
        self.card_short_name = ""
        self.model_type = ""
        self.parameter_type = ""
        self.mutation_string_in_log = ""
        self.nucleotide_mutation_position_s = ""
        self.seq_nuc_length = 0
        self.original_line = line if line is not None else ""

        if line and not is_header:
            cols = split_string(line, '\t')
            if len(cols) >= 8:
                self.sequence_id = cols[0]
                self.aro_number = cols[1]
                self.card_short_name = cols[2]
                self.model_type = cols[3]
                self.parameter_type = cols[4]
                self.mutation_string_in_log = cols[5]
                self.nucleotide_mutation_position_s = cols[6]
                try: self.seq_nuc_length = int(cols[7])
                except (ValueError, IndexError): self.seq_nuc_length = 0
            else:
                self.sequence_id = "ERROR_PARSING"

class AroIndexEntry:
    def __init__(self, line=None, line_num=0):
        self.aro_accession = "N/A"
        self.cvterm_id = ""
        self.model_sequence_id = ""
        self.model_id = ""
        self.model_name = ""
        self.aro_name = ""
        self.protein_accession = ""
        self.dna_accession = ""
        self.amr_gene_family = ""
        self.drug_class = ""
        self.resistance_mechanism = ""
        self.card_short_name = ""

        if line:
            cols = split_string(line, '\t')
            if len(cols) >= 12:
                colon_pos = cols[0].rfind(':')
                self.aro_accession = cols[0][colon_pos + 1:] if colon_pos != -1 and colon_pos + 1 < len(cols[0]) else "N/A"
                self.cvterm_id = cols[1]
                self.model_sequence_id = cols[2]
                self.model_id = cols[3]
                self.model_name = cols[4]
                self.aro_name = cols[5]
                self.protein_accession = cols[6]
                self.dna_accession = cols[7]
                self.amr_gene_family = cols[8]
                self.drug_class = cols[9]
                self.resistance_mechanism = cols[10]
                self.card_short_name = cols[11]

# --- Helper Functions ---
def trim_string(s): return s.strip()

def split_string(s, delimiter): return [trim_string(token) for token in s.split(delimiter)]

def to_uppercase_str(s): return s.upper()

def file_exists(path): return os.path.exists(path) and os.path.isfile(path)

def format_output_filename(suffix_part, ext):
    if not suffix_part:
        return f"{ResScanFiles.COMMON_DB_PREFIX}{ext}"
    else:
        return f"{ResScanFiles.COMMON_DB_PREFIX}_{suffix_part}{ext}"

def sanitize_for_id(s):
    if not s or s == "NA" or s == "Wildtype": 
        return ""
    s = re.sub(r'[ ()\[\]/|*,.;:]', '_', s)
    s = re.sub(r'__+', '_', s) 
    s = s.strip('_') 
    return s

def get_model_code(model_type):
    if "homolog" in model_type: return 'H'
    if "knockout" in model_type: return 'K'
    if "protein variant" in model_type: return 'V'
    if "overexpression" in model_type: return 'O'
    if "rRNA" in model_type: return 'R'
    return 'U'

def generate_unique_id(aro, name, model_type, mutation, id_tracker):
    mut_sanitized = sanitize_for_id(mutation)
    mutation_part = f"__{mut_sanitized}" if mut_sanitized else ""
    base_id = f"ARO_{sanitize_for_id(aro)}__{sanitize_for_id(name)}__{get_model_code(model_type)}{mutation_part}"
    count = id_tracker[base_id]
    id_tracker[base_id] += 1
    return f"{base_id}_{count}" if count > 0 else base_id

# --- Codon Translation ---
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

REVERSE_CODON_TABLE = {
    'A': "GCT", 'C': "TGT", 'D': "GAT", 'E': "GAA", 'F': "TTT",
    'G': "GGT", 'H': "CAT", 'I': "ATT", 'K': "AAA", 'L': "TTA",
    'M': "ATG", 'N': "AAT", 'P': "CCT", 'Q': "CAA", 'R': "CGT",
    'S': "TCT", 'T': "ACT", 'V': "GTT", 'W': "TGG", 'Y': "TAT",
    '*': "TAA"
}

def translate_codon(codon):
    if len(codon) != 3 or '-' in codon or 'N' in codon.upper(): return 'X'
    return CODON_TABLE.get(codon.upper(), 'X')

def translate_nucleotide_sequence(nuc_seq):
    prot_seq = []
    err_msgs = []
    valid_nuc_seq = ""

    for char in to_uppercase_str(nuc_seq):
        if char in "ACGT": valid_nuc_seq += char
        elif char == 'U': valid_nuc_seq += 'T'
        elif not char.isspace(): err_msgs.append(f"Invalid char '{char}'.")

    if len(valid_nuc_seq) < 3:
        err_msgs.append("Sequence too short for translation.")
        return "", "; ".join(err_msgs)

    for i in range(0, len(valid_nuc_seq) - 2, 3):
        codon = valid_nuc_seq[i:i+3]
        aa = translate_codon(codon)
        prot_seq.append(aa)
        if aa == '*': break 

    return "".join(prot_seq), "; ".join(err_msgs)

def join_set_to_string(s_set, delimiter):
    if not s_set: return "NA"
    return delimiter.join(sorted(list(s_set)))

def join_vector_to_delimited_string(s_vec, delimiter):
    if not s_vec: return "NA"
    return delimiter.join(s_vec)

def scan_input_directory(input_dir):
    if not os.path.isdir(input_dir):
        raise RuntimeError(f"Error: Input directory does not exist: {input_dir}")

    config = Config()
    
    def get_path(key):
        filename = ResScanFiles.EXPECTED_INPUTS[key]
        path = os.path.join(input_dir, filename)
        if not os.path.exists(path):
            raise RuntimeError(f"Mandatory file missing in input directory: {filename}")
        return path

    config.path_snps_file = get_path("snps")
    config.path_protein_homolog_model = get_path("homolog")
    config.path_protein_knockout_model = get_path("knockout")
    config.path_protein_variant_model = get_path("variant")
    config.path_protein_overexpression_model = get_path("overexpression")
    config.path_rRNA_gene_variant_model = get_path("rrna")
    config.path_aro_index = get_path("aro_index")

    return config

def setup_logging(output_dir):
    """Sets up logging to a file for auditing purposes."""
    log_file = os.path.join(output_dir, "build_db.log")
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logging.info("--- ResScan Database Build Started ---")
    logging.info(f"Command line: {' '.join(sys.argv)}")
    print(BColors.green(f"Info: Auditing log initialized at: {log_file}"))

def setup_output_directory(dir_path, overwrite):
    if os.path.exists(dir_path):
        if not os.path.isdir(dir_path):
            raise RuntimeError(f"Error: Output path exists but is not a directory: {dir_path}")
        if overwrite:
            print(BColors.yellow("Info: Removing existing output directory."))
            shutil.rmtree(dir_path)
        else:
            raise RuntimeError(f"Error: Output directory '{dir_path}' exists. Use --overwrite to replace it.")
    os.makedirs(dir_path)
    print(BColors.green(f"Info: Created output directory: {dir_path}"))
    return True

def read_next_fasta_entry(fasta_file_handle):
    entry = FastaEntry()
    line = fasta_file_handle.readline()
    while line:
        line = trim_string(line)
        if line and line.startswith('>'):
            entry.header = line
            break
        line = fasta_file_handle.readline()

    if not entry.header: return None

    while True:
        pos = fasta_file_handle.tell()
        line = fasta_file_handle.readline()
        if not line or line.startswith('>'):
            fasta_file_handle.seek(pos)
            break
        entry.sequence += trim_string(line)

    entry.sequence = to_uppercase_str(entry.sequence)
    return entry

def parse_protein_fasta_header(header_line):
    info = ProteinFastaHeaderInfo()
    content = header_line[1:]
    
    match_aro = re.search(r'\|ARO:(\d+)', content)
    if match_aro: info.aro_number = match_aro.group(1)

    last_pipe_pos = content.rfind('|')
    name_part = trim_string(content[last_pipe_pos + 1:]) if last_pipe_pos != -1 else content
    
    if name_part:
        info.name = name_part
        end_pos = len(name_part)
        space_pos = name_part.find(' ')
        bracket_pos = name_part.find('[')
        if space_pos != -1: end_pos = min(end_pos, space_pos)
        if bracket_pos != -1: end_pos = min(end_pos, bracket_pos)
        info.card_short_name = trim_string(name_part[:end_pos])

    info.parsed_successfully = (info.aro_number != "N/A" or info.card_short_name != "N/A")
    return info

def process_protein_model_fasta(in_path, suffix, model_type, metadata_file, ignored_file, db_dir, id_tracker, sequence_to_nucleotide_map):
    try:
        in_fasta = open(in_path, 'r')
    except IOError:
        msg = f"SYSTEM_ERROR\tCould not open: {in_path}"
        ignored_file.write(msg + "\n")
        logging.error(msg)
        return

    out_prot_path = os.path.join(db_dir, format_output_filename(suffix, ".fasta"))
    out_nuc_path = os.path.join(db_dir, format_output_filename(suffix + "_nuc", ".fasta"))

    try:
        out_prot_fasta = open(out_prot_path, 'w')
        out_nuc_fasta = open(out_nuc_path, 'w')
    except IOError as e:
        msg = f"SYSTEM_ERROR\tCould not create output files for {suffix}: {e}"
        ignored_file.write(msg + "\n")
        logging.error(msg)
        in_fasta.close()
        return

    count = 0
    while True:
        entry = read_next_fasta_entry(in_fasta)
        if not entry: break

        prot_seq, err = translate_nucleotide_sequence(entry.sequence)
        h_info = parse_protein_fasta_header(entry.header)

        if h_info.parsed_successfully and prot_seq and not err:
            seq_id = generate_unique_id(h_info.aro_number, h_info.card_short_name, model_type, "", id_tracker)
            out_prot_fasta.write(f">{seq_id}\n{prot_seq}\n")
            out_nuc_fasta.write(f">{seq_id}\n{entry.sequence}\n")
            metadata_file.write(f"{seq_id}\t{h_info.aro_number}\t{h_info.card_short_name}\t"
                                f"{model_type}\tNA\tNA\tNA\t{len(entry.sequence)}\n")
            sequence_to_nucleotide_map[prot_seq] = entry.sequence
            count += 1
        else:
            ignored_file.write(f"{entry.header}\tParse/translation failed: {err}\n")
    
    in_fasta.close()
    out_prot_fasta.close()
    out_nuc_fasta.close()
    logging.info(f"Processed {count} sequences for {model_type}")

def load_fasta_to_map(path):
    fasta_map = {}
    try:
        with open(path, 'r') as f:
            while True:
                entry = read_next_fasta_entry(f)
                if not entry: break
                match_aro = re.search(r'\|ARO:(\d+)', entry.header)
                if match_aro: fasta_map[match_aro.group(1)] = entry
    except IOError:
        return False, fasta_map
    return True, fasta_map

def parse_snps_file_line(line, line_num):
    entry = SNPSFileLine()
    entry.original_line = line
    entry.line_number = line_num
    cols = split_string(line, '\t')
    if len(cols) >= 6:
        entry.accession_aro = cols[0]
        entry.name = cols[1]
        entry.model_type = cols[2]
        entry.parameter_type = cols[3]
        entry.mutations_str = cols[4]
        entry.card_short_name = cols[5]
        
        is_prot = "protein" in entry.model_type
        is_rrna = "rRNA" in entry.model_type
        is_var = "variant" in entry.parameter_type or "mutation" in entry.parameter_type
        
        entry.is_protein_pathway_candidate = is_prot and is_var
        entry.is_rrna_pathway_candidate = is_rrna and is_var
    return entry

def parse_protein_mutations_string(mut_str):
    muts = []
    if not mut_str or trim_string(mut_str) == "NA":
        return True, muts 

    for token in split_string(mut_str, ','):
        if not token: continue
        p_mut = ParsedProteinMutation()
        p_mut.original_token = token 

        # Frameshift: (AA)(Pos)fs
        match_fs = re.match(r'([A-Z])(\d+)fs', token)
        if match_fs:
            p_mut.original_aa = match_fs.group(1)[0]
            p_mut.position = int(match_fs.group(2))
            p_mut.is_frameshift = True
            muts.append(p_mut)
            continue

        # Nonsense: (AA)(Pos)Ter
        match_ter = re.match(r'([A-Z])(\d+)Ter', token)
        if match_ter:
            p_mut.original_aa = match_ter.group(1)[0]
            p_mut.position = int(match_ter.group(2))
            p_mut.new_aa = '*'
            p_mut.is_nonsense = True
            muts.append(p_mut)
            continue

        # Substitution: (OriginalAA)(Pos)(NewAA)
        match_sub = re.match(r'([A-Z])(\d+)([A-Z*])', token)
        if match_sub:
            p_mut.original_aa = match_sub.group(1)[0]
            p_mut.position = int(match_sub.group(2))
            p_mut.new_aa = match_sub.group(3)[0]
            muts.append(p_mut)
            continue
        
        return False, []

    return True, muts

def parse_nucleotide_mutations_string(mut_str):
    muts = []
    if not mut_str or trim_string(mut_str) == "NA":
        return True, muts

    for token in split_string(mut_str, ','):
        if not token: continue
        n_mut = ParsedNucleotideMutation()
        n_mut.original_token = token
        
        match_nuc = re.match(r'([acgtuACGTU])(\d+)([acgtuACGTU])', token)
        if match_nuc:
            n_mut.original_base = to_uppercase_str(match_nuc.group(1))[0]
            n_mut.position = int(match_nuc.group(2))
            n_mut.new_base = to_uppercase_str(match_nuc.group(3))[0]
            if n_mut.original_base == 'U': n_mut.original_base = 'T'
            if n_mut.new_base == 'U': n_mut.new_base = 'T'
            muts.append(n_mut)
            continue
        
        return False, []

    return True, muts

def process_protein_variant_entry(snps, source_map, out_prot_fasta, out_nuc_fasta, metadata_file, ignored_file, id_tracker, sequence_to_nucleotide_map):
    source_entry = source_map.get(snps.accession_aro)
    if not source_entry or not source_entry.sequence:
        ignored_file.write(f"{snps.original_line}\tError: ARO not found or empty sequence in source FASTA\n")
        return

    def write_entry(prot_seq, nuc_seq, mut_tag, nuc_pos):
        if not prot_seq: return
        seq_id = generate_unique_id(snps.accession_aro, snps.card_short_name, snps.model_type, mut_tag, id_tracker)
        out_prot_fasta.write(f">{seq_id}\n{prot_seq}\n")
        out_nuc_fasta.write(f">{seq_id}\n{nuc_seq}\n")
        metadata_file.write(f"{seq_id}\t{snps.accession_aro}\t{snps.card_short_name}\t"
                            f"{snps.model_type}\t{snps.parameter_type}\t{mut_tag}\t"
                            f"{nuc_pos}\t{len(nuc_seq)}\n")
        sequence_to_nucleotide_map[prot_seq] = nuc_seq

    parse_success, mutations = parse_protein_mutations_string(snps.mutations_str)
    if not parse_success:
        ignored_file.write(f"{snps.original_line}\tError: Could not parse protein mutations string\n")
        return

    if not mutations: # Wildtype
        prot_seq, err = translate_nucleotide_sequence(source_entry.sequence)
        if not prot_seq or err:
            ignored_file.write(f"{snps.original_line}\tError: Wildtype translation failed: {err}\n")
            return
        write_entry(prot_seq, source_entry.sequence, "", "NA")
        return

    # Handle frameshift mutations
    if mutations[0].is_frameshift:
        fs = mutations[0]
        nuc_pos_0 = (fs.position - 1) * 3
        for del_len in [1, 2]:
            if nuc_pos_0 >= 0 and nuc_pos_0 + del_len <= len(source_entry.sequence):
                nt_del = list(source_entry.sequence)
                del nt_del[nuc_pos_0 : nuc_pos_0 + del_len]
                nt_del_str = "".join(nt_del)
                
                prot_del, fs_err = translate_nucleotide_sequence(nt_del_str)
                if prot_del and not fs_err:
                    write_entry(prot_del, nt_del_str, f"{fs.original_token}_{del_len}del", str(nuc_pos_0 + 1))
                else:
                    ignored_file.write(f"{snps.original_line}\tError: Frameshift translation failed for {fs.original_token}_{del_len}del: {fs_err}\n")
            else:
                ignored_file.write(f"{snps.original_line}\tError: Frameshift position out of bounds for {fs.original_token}\n")
        return

    # Handle substitution/nonsense mutations
    current_nuc_list = list(source_entry.sequence)
    nuc_pos_vec = []
    
    orig_prot, err_trans = translate_nucleotide_sequence("".join(current_nuc_list))
    if not orig_prot:
        ignored_file.write(f"{snps.original_line}\tError: Initial translation failed for substitution: {err_trans}\n")
        return

    for mut in mutations:
        codon_start_pos = (mut.position - 1) * 3
        if mut.position <= 0 or codon_start_pos + 3 > len(current_nuc_list):
            ignored_file.write(f"{snps.original_line}\tError: Mutation position out of bounds for {mut.original_token}\n")
            return

        current_prot_check, _ = translate_nucleotide_sequence("".join(current_nuc_list))
        if mut.position > len(current_prot_check) or current_prot_check[mut.position-1] != mut.original_aa:
            if mut.position > len(current_prot_check) or current_prot_check[mut.position-1] != mut.new_aa:
                found_aa = current_prot_check[mut.position-1] if mut.position <= len(current_prot_check) else 'End-of-Sequence'
                ignored_file.write(f"{snps.original_line}\tError: AA mismatch for {mut.original_token}. Expected '{mut.original_aa}', found '{found_aa}'\n")
                return

        if mut.new_aa in REVERSE_CODON_TABLE:
            new_codon = REVERSE_CODON_TABLE[mut.new_aa]
            for i in range(3):
                current_nuc_list[codon_start_pos + i] = new_codon[i]
            nuc_pos_vec.append(str(codon_start_pos + 1))
        else:
            ignored_file.write(f"{snps.original_line}\tError: No reverse codon found for new AA '{mut.new_aa}' in {mut.original_token}\n")
            return

    final_nuc_seq = "".join(current_nuc_list)
    final_prot, final_err = translate_nucleotide_sequence(final_nuc_seq)
    if not final_prot or final_err:
        ignored_file.write(f"{snps.original_line}\tError: Final translation failed for substitution: {final_err}\n")
        return
    
    write_entry(final_prot, final_nuc_seq, snps.mutations_str, join_vector_to_delimited_string(nuc_pos_vec, ","))

def process_rrna_variant_entry(snps, rrna_map, out_fasta, metadata_file, ignored_file, id_tracker, sequence_to_nucleotide_map):
    source_entry = rrna_map.get(snps.accession_aro)
    if not source_entry or not source_entry.sequence:
        ignored_file.write(f"{snps.original_line}\tError: ARO not found or empty sequence in source FASTA\n")
        return

    def write_entry(nuc_seq, mut_tag, nuc_pos):
        seq_id = generate_unique_id(snps.accession_aro, snps.card_short_name, snps.model_type, mut_tag, id_tracker)
        out_fasta.write(f">{seq_id}\n{nuc_seq}\n")
        metadata_file.write(f"{seq_id}\t{snps.accession_aro}\t{snps.card_short_name}\t"
                            f"{snps.model_type}\t{snps.parameter_type}\t{mut_tag}\t"
                            f"{nuc_pos}\t{len(nuc_seq)}\n")
        sequence_to_nucleotide_map[nuc_seq] = nuc_seq

    parse_success, mutations = parse_nucleotide_mutations_string(snps.mutations_str)
    if not parse_success:
        ignored_file.write(f"{snps.original_line}\tError: Could not parse rRNA mutations string\n")
        return

    if not mutations: # Wildtype
        write_entry(source_entry.sequence, "", "NA")
        return

    current_nuc_list = list(source_entry.sequence)
    nuc_pos_vec = []
    for mut in mutations:
        if mut.position <= 0 or mut.position > len(current_nuc_list):
            ignored_file.write(f"{snps.original_line}\tError: rRNA Mutation position out of bounds for {mut.original_token}\n")
            return
        
        if current_nuc_list[mut.position - 1] != mut.original_base and current_nuc_list[mut.position - 1] != mut.new_base:
            ignored_file.write(f"{snps.original_line}\tError: rRNA Mutation mismatch for {mut.original_token}. Expected '{mut.original_base}', found '{current_nuc_list[mut.position-1]}'\n")
            return
        
        current_nuc_list[mut.position - 1] = mut.new_base
        nuc_pos_vec.append(str(mut.position))
    
    final_nuc_seq = "".join(current_nuc_list)
    write_entry(final_nuc_seq, snps.mutations_str, join_vector_to_delimited_string(nuc_pos_vec, ","))

def load_all_metadata(path):
    entries = []
    try:
        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader) 
            for row in reader:
                entries.append(MetadataEntry('\t'.join(row)))
    except IOError:
        raise RuntimeError(f"Failed to load Stage 1 metadata from {path}")
    return entries

def load_all_fasta_sequences(path):
    fasta_map = {}
    try:
        with open(path, 'r') as f:
            while True:
                entry = read_next_fasta_entry(f)
                if not entry: break
                if entry.header and len(entry.header) > 1:
                    fasta_map[trim_string(entry.header[1:])] = entry
    except IOError:
        raise RuntimeError(f"Failed to load combined Stage 1 FASTA from {path}")
    return fasta_map

def load_aro_index_file(path):
    aro_index_map = {}
    try:
        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader) 
            for i, row in enumerate(reader):
                entry = AroIndexEntry('\t'.join(row), i + 1)
                if entry.aro_accession != "N/A":
                    aro_index_map[entry.aro_accession] = entry
    except IOError:
        raise RuntimeError(f"Failed to load ARO index file from {path}")
    return aro_index_map


def cleanup_intermediate_files(db_dir):
    print(BColors.cyan("\n--- Cleaning up intermediate files ---"))
    suffixes_to_remove = [
        ResScanFiles.Stage1.FASTA_HOMOLOG,
        ResScanFiles.Stage1.FASTA_KNOCKOUT,
        ResScanFiles.Stage1.FASTA_VARIANT,
        ResScanFiles.Stage1.FASTA_OVEREXPRESSION,
        ResScanFiles.Stage1.FASTA_RRNA,
        ResScanFiles.Stage1.FASTA_COMBINED,
        ResScanFiles.Stage1.METADATA,
        ResScanFiles.Stage2.NR_FASTA 
    ]

    files_removed = 0
    for suffix in suffixes_to_remove:
        for ext in [".fasta", ".txt"]:
            filepath = os.path.join(db_dir, format_output_filename(suffix, ext))
            if os.path.exists(filepath):
                os.remove(filepath)
                files_removed += 1
        
        if "model" in suffix: 
            filepath_nuc = os.path.join(db_dir, format_output_filename(suffix + "_nuc", ".fasta"))
            if os.path.exists(filepath_nuc):
                os.remove(filepath_nuc)
                files_removed += 1

    print(BColors.green(f"Removed {files_removed} intermediate files."))


def run_first_stage(config, db_dir, sequence_to_nucleotide_map):
    print(BColors.cyan("\n--- Stage 1: Processing CARD Models and SNPs ---"))
    logging.info("Starting Stage 1")
    
    stage1_id_tracker = defaultdict(int)

    def open_output_stream(suffix, ext, header):
        path = os.path.join(db_dir, format_output_filename(suffix, ext))
        stream = open(path, 'w')
        if header: stream.write(header + '\n')
        return stream

    metadata_file = open_output_stream(
        ResScanFiles.Stage1.METADATA, ".txt", 
        "Sequence_ID\tARO_Number\tCard_Short_Name\tModel_Type\tParameter_Type\tMutation_String\tNucleotide_Mutation_Position\tSeqNucLength"
    )
    
    # Non-variant stays simple
    non_variant_ignored = open_output_stream(
        ResScanFiles.Stage1.METADATA_NON_VARIANT_IGNORED, ".txt", 
        "Original_Header\tIssue_Description"
    )
    
    # Variant now adopts the snps.txt column structure plus the reason
    variant_header = "Accession_ARO\tName\tModel_Type\tParameter_Type\tMutation\tCard_Short_Name\tSource\tCitation\tIssue_Description"
    variant_ignored = open_output_stream(
        ResScanFiles.Stage1.METADATA_VARIANT_IGNORED, ".txt", 
        variant_header
    )

    print(BColors.cyan("Processing homolog models..."))
    process_protein_model_fasta(config.path_protein_homolog_model, ResScanFiles.Stage1.FASTA_HOMOLOG, "homolog", metadata_file, non_variant_ignored, db_dir, stage1_id_tracker, sequence_to_nucleotide_map)
    
    print(BColors.cyan("Processing knockout models..."))
    process_protein_model_fasta(config.path_protein_knockout_model, ResScanFiles.Stage1.FASTA_KNOCKOUT, "knockout", metadata_file, non_variant_ignored, db_dir, stage1_id_tracker, sequence_to_nucleotide_map)
    
    print(BColors.cyan("Loading variant FASTA files..."))
    success, protein_variant_map = load_fasta_to_map(config.path_protein_variant_model)
    if not success: raise RuntimeError("Failed to load Protein Variant FASTA")
    success, protein_overexpression_map = load_fasta_to_map(config.path_protein_overexpression_model)
    if not success: raise RuntimeError("Failed to load Protein Overexpression FASTA")
    success, rrna_variant_map = load_fasta_to_map(config.path_rRNA_gene_variant_model)
    if not success: raise RuntimeError("Failed to load rRNA Variant FASTA")
    
    prot_var_fasta = open_output_stream(ResScanFiles.Stage1.FASTA_VARIANT, ".fasta", "")
    prot_var_nuc_fasta = open_output_stream(ResScanFiles.Stage1.FASTA_VARIANT + "_nuc", ".fasta", "")
    prot_over_fasta = open_output_stream(ResScanFiles.Stage1.FASTA_OVEREXPRESSION, ".fasta", "")
    prot_over_nuc_fasta = open_output_stream(ResScanFiles.Stage1.FASTA_OVEREXPRESSION + "_nuc", ".fasta", "")
    rrna_fasta = open_output_stream(ResScanFiles.Stage1.FASTA_RRNA, ".fasta", "")

    print(BColors.cyan("Processing snps.txt for variants..."))
    try:
        with open(config.path_snps_file, 'r') as snps_file:
            csv_reader = csv.reader(snps_file, delimiter='\t')
            next(csv_reader) 
            for i, row in enumerate(csv_reader):
                line = '\t'.join(row)
                entry = parse_snps_file_line(line, i + 1)
                if entry.is_protein_pathway_candidate:
                    if "variant" in entry.model_type:
                        process_protein_variant_entry(entry, protein_variant_map, prot_var_fasta, prot_var_nuc_fasta, metadata_file, variant_ignored, stage1_id_tracker, sequence_to_nucleotide_map)
                    elif "overexpression" in entry.model_type:
                        process_protein_variant_entry(entry, protein_overexpression_map, prot_over_fasta, prot_over_nuc_fasta, metadata_file, variant_ignored, stage1_id_tracker, sequence_to_nucleotide_map)
                elif entry.is_rrna_pathway_candidate:
                    process_rrna_variant_entry(entry, rrna_variant_map, rrna_fasta, metadata_file, variant_ignored, stage1_id_tracker, sequence_to_nucleotide_map)
                elif entry.accession_aro != "": 
                    variant_ignored.write(f"{entry.original_line}\tSkipped: Model/Parameter Type not applicable.\n")
    except IOError:
        raise RuntimeError(f"Failed to open snps file: {config.path_snps_file}")
    
    metadata_file.close()
    non_variant_ignored.close()
    variant_ignored.close()
    prot_var_fasta.close()
    prot_var_nuc_fasta.close()
    prot_over_fasta.close()
    prot_over_nuc_fasta.close()
    rrna_fasta.close()
    print(BColors.green("--- Stage 1 Complete ---"))
    logging.info("Stage 1 Complete")



def run_second_stage(config, db_dir, sequence_to_nucleotide_map):
    print(BColors.cyan("\n--- Stage 2: Removing Redundancy and Finalizing Outputs ---"))
    logging.info("Starting Stage 2: Redundancy Removal")

    combined_fasta_path = os.path.join(db_dir, format_output_filename(ResScanFiles.Stage1.FASTA_COMBINED, ".fasta"))
    
    with open(combined_fasta_path, 'w') as combined_stream:
        for suffix in [ResScanFiles.Stage1.FASTA_HOMOLOG, ResScanFiles.Stage1.FASTA_KNOCKOUT, ResScanFiles.Stage1.FASTA_VARIANT, ResScanFiles.Stage1.FASTA_OVEREXPRESSION, ResScanFiles.Stage1.FASTA_RRNA]:
            individual_fasta_path = os.path.join(db_dir, format_output_filename(suffix, ".fasta"))
            if os.path.exists(individual_fasta_path):
                with open(individual_fasta_path, 'r') as individual_stream:
                    shutil.copyfileobj(individual_stream, combined_stream)
    
    all_metadata = load_all_metadata(os.path.join(db_dir, format_output_filename(ResScanFiles.Stage1.METADATA, ".txt")))
    id_to_fasta = load_all_fasta_sequences(combined_fasta_path)
    aro_index_map = load_aro_index_file(config.path_aro_index)

    id_to_metadata = {}
    sequence_to_ids = defaultdict(list)
    for meta in all_metadata:
        if meta.sequence_id == "ERROR_PARSING": continue
        id_to_metadata[meta.sequence_id] = meta
        if meta.sequence_id in id_to_fasta:
            sequence_to_ids[id_to_fasta[meta.sequence_id].sequence].append(meta.sequence_id)
    
    logging.info(f"Consolidating {len(id_to_metadata)} metadata entries and {len(id_to_fasta)} sequences.")

    nr_meta_path = os.path.join(db_dir, format_output_filename(ResScanFiles.Stage2.NR_METADATA, ".txt"))
    nr_fasta_path = os.path.join(db_dir, format_output_filename(ResScanFiles.Stage2.NR_FASTA, ".fasta"))
    nr_nuc_fasta_path = os.path.join(db_dir, format_output_filename(ResScanFiles.Stage2.NR_FASTA_NUC, ".fasta"))

    with open(nr_meta_path, 'w', newline='') as nr_meta_stream, \
         open(nr_fasta_path, 'w') as nr_fasta_stream, \
         open(nr_nuc_fasta_path, 'w') as nr_nuc_fasta_stream:
        
        nr_meta_writer = csv.writer(nr_meta_stream, delimiter='\t')
        nr_meta_writer.writerow(["Sequence_ID", "ARO_Number", "Card_Short_Name", "Model_Type", "Parameter_Type",
                                 "Mutation_String", "Nucleotide_Mutation_Position", "SeqNucLength",
                                 "Consolidated_Stage1_IDs", "Consolidated_ARO_Numbers", "Consolidated_Mutation_Strings",
                                 "Consolidated_Parameter_Types", "Redundancy_Count", "AMR_Gene_Family", "ARO_Name",
                                 "Drug_Class", "Resistance_Mechanism"])
        
        final_id_tracker = defaultdict(int)

        for unique_sequence, ids in sequence_to_ids.items():
            if not ids: continue

            rep_meta = id_to_metadata[ids[0]] 

            consolidated_aros = set()
            consolidated_mutations = set()
            consolidated_params = set()
            consolidated_families = set()
            consolidated_aro_names = set()
            consolidated_drugs = set()
            consolidated_mechs = set()

            for id_val in ids:
                meta = id_to_metadata[id_val]
                if meta.aro_number != "NA":
                    consolidated_aros.add(meta.aro_number)
                    if meta.aro_number in aro_index_map:
                        aro_data = aro_index_map[meta.aro_number]
                        if aro_data.amr_gene_family and aro_data.amr_gene_family != "NA": consolidated_families.add(aro_data.amr_gene_family)
                        if aro_data.aro_name and aro_data.aro_name != "NA": consolidated_aro_names.add(aro_data.aro_name)
                        if aro_data.drug_class and aro_data.drug_class != "NA": consolidated_drugs.add(aro_data.drug_class)
                        if aro_data.resistance_mechanism and aro_data.resistance_mechanism != "NA": consolidated_mechs.add(aro_data.resistance_mechanism)
                if meta.mutation_string_in_log not in ["NA", "Wildtype"]:
                    consolidated_mutations.add(meta.mutation_string_in_log)
                if meta.parameter_type != "NA":
                    consolidated_params.add(meta.parameter_type)

            consolidated_aros_str = join_set_to_string(consolidated_aros, ";")
            consolidated_muts_str = join_set_to_string(consolidated_mutations, ";")
            consolidated_muts_sanitized_for_header = sanitize_for_id(consolidated_muts_str)
            mutation_part_for_header = f"__{consolidated_muts_sanitized_for_header}" if consolidated_muts_sanitized_for_header else ""

            base_id = (f"ARO_{sanitize_for_id(consolidated_aros_str)}__{sanitize_for_id(rep_meta.card_short_name)}__"
                       f"{get_model_code(rep_meta.model_type)}{mutation_part_for_header}")
            
            count = final_id_tracker[base_id]
            final_id_tracker[base_id] += 1
            final_id = f"{base_id}_{count}" if count > 0 else base_id

            if count > 0:
                warning_msg = f"Warning: Generated identical base ID '{base_id}' for a new unique sequence. Appending suffix: {final_id}"
                sys.stderr.write(BColors.yellow(warning_msg + "\n"))
                logging.warning(warning_msg)

            nr_fasta_stream.write(f">{final_id}\n{unique_sequence}\n")
            
            nuc_seq_for_protein = sequence_to_nucleotide_map.get(unique_sequence)
            if nuc_seq_for_protein:
                nr_nuc_fasta_stream.write(f">{final_id}\n{nuc_seq_for_protein}\n")
            else:
                logging.warning(f"Could not find corresponding nucleotide sequence for final ID: {final_id}")
            
            nr_meta_writer.writerow([
                final_id, rep_meta.aro_number, rep_meta.card_short_name, rep_meta.model_type, rep_meta.parameter_type,
                rep_meta.mutation_string_in_log, rep_meta.nucleotide_mutation_position_s, rep_meta.seq_nuc_length,
                join_vector_to_delimited_string(ids, ";"),
                consolidated_aros_str,
                (consolidated_muts_str if consolidated_muts_str != "NA" else ""), 
                join_set_to_string(consolidated_params, ";"),
                len(ids),
                join_set_to_string(consolidated_families, ";"),
                join_set_to_string(consolidated_aro_names, ";"),
                join_set_to_string(consolidated_drugs, ";"),
                join_set_to_string(consolidated_mechs, ";")
            ])
    print(BColors.green("--- Stage 2 Complete ---"))
    logging.info("Stage 2 Complete")

def write_readme_file(config, db_dir):
    readme_path = os.path.join(db_dir, format_output_filename(ResScanFiles.Stage1.README_SUFFIX, ".txt"))
    try:
        with open(readme_path, 'w') as readme_file:
            readme_file.write("ResScan Database Build Documentation\n")
            readme_file.write("==================================\n\n")
            readme_file.write(f"Generated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            readme_file.write("This document describes the final output files.\n\n")
            
            readme_file.write("--- Final FASTA File ---\n")
            readme_file.write(f"**{format_output_filename(ResScanFiles.Stage2.NR_FASTA_NUC, '.fasta')}**\n")
            readme_file.write("   - Contains the non-redundant nucleotide sequences for all CARD models.\n")
            readme_file.write("   - This is the primary reference database used by ResScan for alignment.\n\n")
            
            readme_file.write("**Header Format**\n")
            readme_file.write("Format: `>ARO_[AROs]__[GeneName]__[ModelCode][__Mutations]`\n")
            readme_file.write("Example: `>ARO_3000597__AAC_6_-Iad__V__S181G` (with mutation)\n")
            readme_file.write("Example: `>ARO_3000001__blaTEM-1__H` (wildtype, no mutation part)\n")
            readme_file.write(" - `[AROs]`: Consolidated ARO numbers for this sequence.\n")
            readme_file.write(" - `[GeneName]`: The CARD short name of a representative entry.\n")
            readme_file.write(" - `[ModelCode]`: Single character for model type (H:homolog, K:knockout, V:variant, O:overexpression, R:rRNA).\n")
            readme_file.write(" - `[__Mutations]`: Optional. Consolidated mutations, if present. Omitted for wildtype sequences.\n\n")
            
            readme_file.write(f"--- Final Metadata File: {format_output_filename(ResScanFiles.Stage2.NR_METADATA, '.txt')} ---\n")
            readme_file.write("This file contains consolidated metadata for each unique sequence.\n\n")
            readme_file.write("Columns:\n")
            readme_file.write("1.  **Sequence_ID**: The unique ID that matches the FASTA headers. This is the primary key.\n")
            readme_file.write("2-8. **ARO_Number, etc.**: Metadata from a single representative pre-consolidation entry.\n")
            readme_file.write("9.  **Consolidated_Stage1_IDs**: A semicolon-separated list of all original IDs from Stage 1 that were merged into this single unique sequence. Provides a full audit trail.\n")
            readme_file.write("10. **Consolidated_ARO_Numbers**: All unique ARO numbers for this sequence.\n")
            readme_file.write("11. **Consolidated_Mutation_Strings**: All unique mutations that result in this sequence (empty if wildtype).\n")
            readme_file.write("12-17. Additional consolidated metadata (Parameter Types, Redundancy Count, Gene Family, etc.).\n")
            
        print(BColors.green(f"README file written to: {readme_path}"))
    except IOError as e:
        print(BColors.red(f"Error writing README file: {e}"), file=sys.stderr)



def write_html_ignored_report(db_dir):
    """Generates a detailed HTML report with distinct structures for different ignored categories."""
    html_path = os.path.join(db_dir, "ignored_sequences_report.html")
    non_var_path = os.path.join(db_dir, format_output_filename(ResScanFiles.Stage1.METADATA_NON_VARIANT_IGNORED, ".txt"))
    var_path = os.path.join(db_dir, format_output_filename(ResScanFiles.Stage1.METADATA_VARIANT_IGNORED, ".txt"))

    css = """
    <style>
        body { font-family: 'Segoe UI', sans-serif; line-height: 1.4; color: #333; max-width: 1600px; margin: 0 auto; padding: 20px; background-color: #f8f9fa; }
        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
        h2 { color: #d35400; margin-top: 40px; padding: 5px 10px; background: #fff; border-left: 5px solid #d35400; }
        .table-container { overflow-x: auto; background: #fff; box-shadow: 0 4px 6px rgba(0,0,0,0.1); border-radius: 8px; margin: 20px 0; }
        table { width: 100%; border-collapse: collapse; font-size: 0.85em; }
        th { background-color: #2c3e50; color: white; text-align: left; padding: 12px; position: sticky; top: 0; }
        td { padding: 10px; border-bottom: 1px solid #eee; word-break: break-word; vertical-align: top; }
        tr:nth-child(even) { background-color: #fcfcfc; }
        tr:hover { background-color: #f1f4f6; }
        .reason-cell { font-weight: bold; color: #c0392b; }
        .no-data { font-style: italic; color: #95a5a6; text-align: center; padding: 40px; }
        .summary { margin-bottom: 20px; padding: 15px; background: #eef2f7; border-radius: 6px; border: 1px solid #d1d9e1; display: inline-block; }
    </style>
    """

    def get_data(file_path):
        data = []
        if os.path.exists(file_path):
            with open(file_path, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                headers = next(reader, None)
                for row in reader:
                    if row: data.append(row)
        return data

    non_var_data = get_data(non_var_path)
    var_data = get_data(var_path)

    with open(html_path, 'w') as h:
        h.write(f"<!DOCTYPE html><html><head><meta charset='utf-8'><title>ResScan Ignored Report</title>{css}</head><body>")
        h.write("<h1>ResScan Database Build: Ignored Sequences Report</h1>")
        h.write(f"<div class='summary'><strong>Non-Variant Ignored:</strong> {len(non_var_data)} | <strong>Variant Ignored:</strong> {len(var_data)}</div>")

        # Table 1: Non-Variant
        h.write("<h2>1. Non-Variant Models (Translation/Header Issues)</h2>")
        if non_var_data:
            h.write("<div class='table-container'><table><thead><tr><th>Original Header</th><th>Issue Description</th></tr></thead><tbody>")
            for row in non_var_data:
                h.write(f"<tr><td>{row[0]}</td><td class='reason-cell'>{row[1]}</td></tr>")
            h.write("</tbody></table></div>")
        else:
            h.write("<p class='no-data'>No non-variant sequences were ignored.</p>")

        # Table 2: Variant (Detailed)
        h.write("<h2>2. Variant & rRNA Models (Detailed SNP Data)</h2>")
        if var_data:
            v_headers = ["ARO", "Name", "Model", "Param", "Mutation", "Short Name", "Source", "Citation", "Issue Description"]
            h.write("<div class='table-container'><table><thead><tr>")
            for head in v_headers: h.write(f"<th>{head}</th>")
            h.write("</tr></thead><tbody>")
            for row in var_data:
                h.write("<tr>")
                for i, cell in enumerate(row):
                    cls = " class='reason-cell'" if i == len(row)-1 else ""
                    h.write(f"<td{cls}>{cell}</td>")
                h.write("</tr>")
            h.write("</tbody></table></div>")
        else:
            h.write("<p class='no-data'>No variant entries were ignored.</p>")

        h.write(f"<p style='color: #7f8c8d; font-size: 0.8em;'>Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>")
        h.write("</body></html>")


def main():
    parser = argparse.ArgumentParser(description="resscan_build_db - Processes CARD data from a directory into ResScan format.")
    parser.add_argument("-i", "--input-dir", required=True, help="Path to the folder containing downloaded CARD files.")
    parser.add_argument("-d", "--output-dir", required=True, help="Name of the output directory.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output directory.")
    parser.add_argument("--keep-intermediate", action="store_true", help="Keep intermediate build files (Stage 1 outputs and Protein FASTA).")
    
    args = parser.parse_args()

    try:
        config = scan_input_directory(args.input_dir)
        config.output_dir_name = args.output_dir
        config.overwrite_outputs = args.overwrite
        config.keep_intermediate = args.keep_intermediate

        setup_output_directory(config.output_dir_name, config.overwrite_outputs)
        setup_logging(config.output_dir_name)
        
        db_dir = config.output_dir_name

        print(BColors.cyan(f"Starting ResScan Database Build..."))
        
        sequence_to_nucleotide_map = {}
        
        run_first_stage(config, db_dir, sequence_to_nucleotide_map)
        run_second_stage(config, db_dir, sequence_to_nucleotide_map)
        write_readme_file(config, db_dir)
        write_html_ignored_report(db_dir)
        
        if not config.keep_intermediate:
            cleanup_intermediate_files(db_dir)
        else:
            print(BColors.yellow("\n--- Skipping cleanup (intermediate files preserved) ---"))
        
        print(BColors.green("\n--- All Stages Complete ---"))
        print(BColors.green(f"Final processed files are in '{config.output_dir_name}' directory."))
        logging.info("ResScan Database Build Completed Successfully")

    except Exception as e:
        sys.stderr.write(BColors.red(f"Critical Error: {e}\n"))
        logging.error(f"Critical Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()