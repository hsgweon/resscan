"""
Tests for how alignment records are turned into evidence.

Each test builds a tiny reference database and hand-written SAM records for one
specific situation, so a failure points at one rule rather than at the pipeline.
"""
import csv
import random
import subprocess
import sys
from pathlib import Path

import pytest

PKG = Path(__file__).resolve().parent.parent / "resscan"
sys.path.insert(0, str(PKG.parent))

from resscan import varscan_process_sam as vps          # noqa: E402
from resscan import scgscan_quantify_from_diamond as sq  # noqa: E402

CODONS = {'S': 'TCG', 'L': 'TTA', 'D': 'GAC', 'N': 'AAT'}


def revcomp(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def make_gene(length_codons=120, seed=1):
    """A random open reading frame with no stop codons, S at codon 10 and D at codon 60."""
    rng = random.Random(seed)
    sense = [c for c, aa in vps.CODON_TABLE.items() if aa != '_']
    codons = [rng.choice(sense) for _ in range(length_codons)]
    codons[9], codons[59] = CODONS['S'], CODONS['D']
    return codons


def with_mutations(codons, **changes):
    out = list(codons)
    for position, residue in changes.items():
        out[int(position[1:]) - 1] = CODONS[residue]
    return "".join(out)


@pytest.fixture
def variant_db(tmp_path):
    wild = make_gene()
    refs = {
        "ARO_1__Test_gene__V__S10L": with_mutations(wild, p10='L'),
        "ARO_1__Test_gene__V__S10L_D60N": with_mutations(wild, p10='L', p60='N'),
    }
    fasta = tmp_path / "db.fasta"
    fasta.write_text("".join(f">{name}\n{seq}\n" for name, seq in refs.items()))
    return fasta, refs


def sam_line(qname, flag, rname, pos, cigar, seq, nm, score):
    return "\t".join([qname, str(flag), rname, str(pos), "0", cigar, "*", "0", "0",
                      seq, "*" if seq == "*" else "I" * len(seq), f"NM:i:{nm}", f"AS:i:{score}"])


def run_varscan(tmp_path, fasta, lines):
    sam = tmp_path / "reads.sam"
    sam.write_text("\n".join(lines) + "\n")
    sequences = vps.parse_fasta_db(str(fasta))
    mutations = vps.parse_mutation_info_from_db(sequences)
    out = tmp_path / "hits.tsv"
    vps.process_sam_files_for_variants(
        [str(sam)], sequences, mutations, str(out), False, str(tmp_path / "debug.tsv"),
        str(tmp_path / "aln.txt"), 100, False, {'V', 'R'})
    with open(out) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def test_single_substitution_read_survives_arbitrary_primary(tmp_path, variant_db):
    """
    A read with S10L that stops before codon 60 matches both references
    perfectly. If the aligner makes the combination its primary, the read must
    still be counted for the single-substitution model.
    """
    fasta, refs = variant_db
    read = refs["ARO_1__Test_gene__V__S10L"][0:150]          # codons 1-50
    hits = run_varscan(tmp_path, fasta, [
        sam_line("r1", 0, "ARO_1__Test_gene__V__S10L_D60N", 1, "150M", read, 0, 150),
        sam_line("r1", 256, "ARO_1__Test_gene__V__S10L", 1, "150M", "*", 0, 150),
    ])
    assert [h["ARO_matched"] for h in hits] == ["ARO_1__Test_gene__V__S10L"]


def test_double_mutant_read_counted_once_for_combination(tmp_path, variant_db):
    """A read carrying both substitutions confirms both models but counts once, for the pair."""
    fasta, refs = variant_db
    read = refs["ARO_1__Test_gene__V__S10L_D60N"][0:210]     # codons 1-70
    hits = run_varscan(tmp_path, fasta, [
        sam_line("r1", 0, "ARO_1__Test_gene__V__S10L", 1, "210M", read, 3, 195),
        sam_line("r1", 256, "ARO_1__Test_gene__V__S10L_D60N", 1, "210M", "*", 0, 210),
    ])
    assert [h["ARO_matched"] for h in hits] == ["ARO_1__Test_gene__V__S10L_D60N"]


def test_wild_type_read_confirms_nothing(tmp_path, variant_db):
    fasta, refs = variant_db
    read = with_mutations(make_gene())[0:150]
    nm = sum(a != b for a, b in zip(read, refs["ARO_1__Test_gene__V__S10L"]))
    hits = run_varscan(tmp_path, fasta, [
        sam_line("r1", 0, "ARO_1__Test_gene__V__S10L", 1, "150M", read, nm, 140),
        sam_line("r1", 256, "ARO_1__Test_gene__V__S10L_D60N", 1, "150M", "*", nm, 140),
    ])
    assert hits == []


def test_secondary_on_opposite_strand_is_reoriented(tmp_path, variant_db):
    """A secondary on the other strand must be read from the reverse complement of the primary's SEQ."""
    fasta, refs = variant_db
    read = refs["ARO_1__Test_gene__V__S10L"][0:150]
    hits = run_varscan(tmp_path, fasta, [
        # Primary on a reference VarScan does not track, stored reverse-complemented.
        sam_line("r1", 16, "ARO_9__Other__H", 1, "150M", revcomp(read), 0, 150),
        sam_line("r1", 256, "ARO_1__Test_gene__V__S10L", 1, "150M", "*", 0, 150),
    ])
    assert [h["ARO_matched"] for h in hits] == ["ARO_1__Test_gene__V__S10L"]


def test_reads_are_grouped_per_read_not_per_file(tmp_path, variant_db):
    fasta, refs = variant_db
    read = refs["ARO_1__Test_gene__V__S10L"][0:150]
    hits = run_varscan(tmp_path, fasta, [
        sam_line("r1", 0, "ARO_1__Test_gene__V__S10L", 1, "150M", read, 0, 150),
        sam_line("r2", 0, "ARO_1__Test_gene__V__S10L", 1, "150M", read, 0, 150),
    ])
    assert sorted(h["query_id"] for h in hits) == ["r1_1", "r2_1"]


def test_clipping_rule_matches_reference_ends():
    """A clip is legitimate only where the read runs off the end of the reference."""
    assert vps.is_valid_alignment("20S130M", 1, 130, 500)        # overhangs the start
    assert vps.is_valid_alignment("130M20S", 371, 130, 500)      # overhangs the end
    assert not vps.is_valid_alignment("20S130M", 200, 130, 500)  # clipped inside the gene
    assert not vps.is_valid_alignment("130M20S", 1, 130, 500)


def test_homscan_ignores_hard_clipped_supplementary_sequence(tmp_path):
    """A supplementary record written after the primary must not replace the full read."""
    gene = with_mutations(make_gene(seed=7))
    fasta = tmp_path / "db.fasta"
    fasta.write_text(f">ARO_5__GeneA__H\n{gene}\n>ARO_6__GeneB__H\n{with_mutations(make_gene(seed=8))}\n")
    read = gene[len(gene) - 100:] + "ACGT" * 12 + "AC"           # runs 50 bp past the gene end
    sam = tmp_path / "reads.fq.sam"
    sam.write_text("\n".join([
        sam_line("r1", 0, "ARO_5__GeneA__H", len(gene) - 99, "100M50S", read, 0, 100),
        sam_line("r1", 2048, "ARO_6__GeneB__H", 1, "100H50M", read[100:], 20, 30),
    ]) + "\n")
    subprocess.run([sys.executable, str(PKG / "homscan_process_sam.py"), "-i", str(sam), "--db", str(fasta),
                    "--tmp-dir", str(tmp_path), "--output-prefix", "t", "--pid-cutoff", "0.95"],
                   check=True, capture_output=True)
    with open(tmp_path / "reads.fq.sam_hits.tsv") as handle:
        hits = list(csv.DictReader(handle, delimiter='\t'))
    assert [(h["ARO_matched"], h["protein_pid"]) for h in hits] == [("ARO_5__GeneA__H", "100.00")]


def test_marker_density_is_per_kilobase_of_nucleotides(tmp_path):
    """
    A marker and an AMR gene of equal physical length that attract equal reads
    are both present once per cell, so their densities must be equal.
    """
    report = tmp_path / "report.tsv"
    rpk, fpk = sq.write_report(str(report), {"geneA": 90}, {"geneA": set(range(90))}, {"geneA": 300}, 90)
    amr_rpk = 90 / (900 / 1000.0)                                # same reads over the same 900 nt
    assert rpk == pytest.approx(amr_rpk)
    assert fpk == pytest.approx(amr_rpk)


def _diamond_line(read, gene, qstart, qend, qlen=None, evalue="1e-20"):
    fields = [read, gene, "100", str((abs(qend - qstart) + 1) // 3), "0", "0", str(qstart), str(qend), "1", "50", evalue, "100"]
    if qlen is not None:
        fields.append(str(qlen))
    return "\t".join(fields)


def test_marker_reads_need_half_the_read_aligned(tmp_path):
    """Markers are counted by the rule used for AMR genes: at least half of the read must align."""
    hits = tmp_path / "hits.tsv"
    hits.write_text("\n".join([
        _diamond_line("r1", "rplB___001", 1, 150, 150),        # whole read
        _diamond_line("r2", "rplB___001", 150, 76, 150),       # exactly half, reverse frame
        _diamond_line("r3", "rplB___001", 1, 72, 150),         # 72 of 150 bases: an end overlap
        _diamond_line("r4", "rplB___001", 1, 45, 60),          # short read, 45 of 60 aligned
        _diamond_line("r5", "rplB___001", 1, 36, 60),          # above half the read but below the 40-base floor
        _diamond_line("r6", "rplB___001", 1, 150, 150, evalue="1e-3"),
    ]) + "\n")
    counts, _fragments, total, passing, too_short = sq.process_diamond_file(str(hits), 1e-5, 40, 0.5)
    assert counts["rplB"] == 3 and passing == 3 and too_short == 2 and total == 6
    counts, *_ = sq.process_diamond_file(str(hits), 1e-5)       # rule off: e-value only
    assert counts["rplB"] == 5


def test_marker_files_without_read_length_are_counted_on_evalue_alone(tmp_path):
    hits = tmp_path / "old.tsv"
    hits.write_text(_diamond_line("r1", "rplB___001", 1, 30) + "\n")
    counts, _fragments, _total, passing, too_short = sq.process_diamond_file(str(hits), 1e-5, 40, 0.5)
    assert counts["rplB"] == 1 and passing == 1 and too_short == 0


def test_marker_density_uses_the_effective_length_offset(tmp_path):
    report = tmp_path / "report.tsv"
    rpk, fpk = sq.write_report(str(report), {"geneA": 90}, {"geneA": set(range(90))}, {"geneA": 300}, 90, eff_len_offset=-100.0)
    assert rpk == pytest.approx(90 / (800 / 1000.0))             # reads: 900 nt - 100
    assert fpk == pytest.approx(90 / (900 / 1000.0))             # fragments: plain length, as for AMR genes


def test_packaged_marker_set():
    """35 markers; the five that attract reads from related genes are not among them."""
    db = PKG / "databases" / "resscan_DB_SCG"
    lengths = [line.split("\t")[0] for line in (db / "SCG_gene_lengths.tsv").read_text().splitlines()[1:] if line.strip()]
    in_fasta = {line[1:].split("___")[0] for line in open(db / "SCG.fasta") if line.startswith(">")}
    assert len(lengths) == 35 and set(lengths) == in_fasta
    assert not {"rpoA", "rpoB", "cysS", "hisS", "valS"} & in_fasta
