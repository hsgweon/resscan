"""Tests for the HTML reports: where mutation markers go, and what each page embeds."""
import json
import re
import subprocess
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent / "resscan"
sys.path.insert(0, str(PKG.parent))

from resscan import varscan_visualise as vv  # noqa: E402


def test_marker_positions_ignore_substitution_like_entry_names():
    """'P2x' in PBP2x, 'P1a' in PBP1a and 'E2A' in GE2A are part of the name, not mutations."""
    for name, mutation, expected in [
        ("ARO_3003043__Spne_PBP2x_AMX__V__M339F", "M339F", [339]),
        ("ARO_3003041__Spne_PBP1a_AMX__V__S575T", "S575T", [575]),
        ("ARO_3003438__Efac_EFTu_GE2A__V__G275A", "G275A", [275]),
        ("ARO_3004174__Ctra_23S_MAC__R__a2663c_t2663c", "a2663c_t2663c", [2663, 2663]),
    ]:
        header = f"> SRR1.2_1 | maps to | {name} | mutation: {mutation}"
        assert vv.parse_header_info(header) == (name.split("__")[2], expected)


def test_marker_line_marks_the_codon_and_skips_gap_columns():
    ref = "ATG" + "---" + "AAACCCGGG"            # an insertion in the read after codon 1
    marker = vv.generate_marker_line(ref, 'V', [3], 1)   # codon 3 = reference bases 7-9
    assert marker == "." * 9 + "^^^" + "..."


def test_mismatch_highlighting_marks_only_aligned_differences():
    out = vv.highlight_mismatches("ACGTACGT", "..GAAC-.")
    assert out == '..G<span class="mismatch">A</span>AC<span class="mismatch">-</span>.'
    assert re.sub(r"<.*?>", "", out) == "..GAAC-."


def test_homscan_page_embeds_only_its_own_reads(tmp_path):
    db = tmp_path / "db.fasta"
    db.write_text(">ARO_1__geneA__H\n" + "ACGT" * 100 + "\n>ARO_2__geneB__H\n" + "TTGA" * 100 + "\n")
    meta = tmp_path / "meta.tsv"
    meta.write_text("Sequence_ID\tAMR_Gene_Family\nARO_1__geneA__H\tfamily A\nARO_2__geneB__H\tfamily B\n")
    cols = ["ARO_matched", "query_id", "nucleotide_pid", "nucleotide_mismatches", "nucleotide_denominator",
            "protein_pid", "protein_mismatches", "protein_denominator", "position_on_ref", "read_length",
            "target_length", "uniqueness"]
    rows = [["ARO_1__geneA__H", "readA_1", "100.00", "0", "150", "100.00", "0", "50", "1", "150", "400", "UNIQUE"],
            ["ARO_2__geneB__H", "readB_1", "100.00", "0", "150", "100.00", "0", "50", "1", "150", "400", "UNIQUE"],
            ["ARO_1__geneA__H", "shared_1", "100.00", "0", "150", "100.00", "0", "50", "10", "150", "400", "NOT_UNIQUE"],
            ["ARO_2__geneB__H", "shared_1", "100.00", "0", "150", "100.00", "0", "50", "10", "150", "400", "NOT_UNIQUE"]]
    hits = tmp_path / "hits.tsv"
    hits.write_text("\n".join("\t".join(r) for r in [cols] + rows) + "\n")
    out = tmp_path / "html"
    subprocess.run([sys.executable, str(PKG / "homscan_visualise.py"), "-i", str(hits), "--db", str(db),
                    "--metadata", str(meta), "-o", str(out)], check=True, capture_output=True)
    page = (out / "family_a.html").read_text()
    embedded = json.loads(re.search(r"const allReadHits = (\{.*?\});\n", page).group(1))
    assert set(embedded["reads"]) == {"readA_1", "shared_1"}
    # a read shown here still lists its match in the other family
    names = [reference for reference, _family in embedded["references"]]
    shared = embedded["reads"]["shared_1"]
    assert shared[0] == 2
    assert {names[row[0]] for row in shared[1:]} == {"ARO_1__geneA__H", "ARO_2__geneB__H"}


def test_homscan_depth_and_coverage_use_every_read_even_when_the_pile_is_capped():
    from resscan import homscan_visualise as hv
    hits = [{"query_id": f"r{i}_1", "ARO_matched": "ARO_1__geneA__H", "position_on_ref": str(1 + i % 200),
             "nucleotide_denominator": "100", "nucleotide_pid": "100.00", "protein_pid": "100.00",
             "uniqueness": "UNIQUE"} for i in range(1000)]
    depth = hv.compute_depth(400, hits)
    assert sum(depth) == 1000 * 100                      # every aligned base of every read is counted
    section, stats, drawn = hv.generate_single_aro_plot_html("ARO_1__geneA__H", 400, hits, len(hits), "protein", 50)
    assert len(drawn) == 50 and section.count('<rect class="r') == 50
    assert stats["reads"] == 1000 and stats["lateral_coverage_pct"] == round(299 / 400 * 100, 2)
    assert "50 of 1,000 shown" in section


def test_page_read_budget_lowers_the_limit_only_for_pages_with_many_references():
    from resscan import homscan_visualise as hv
    assert hv.reads_drawn_limit(500, 2) == 500
    assert hv.reads_drawn_limit(500, 230) == hv.PAGE_READ_BUDGET // 230
    assert hv.reads_drawn_limit(500, 5000) == hv.MIN_READS_DRAWN
    assert hv.reads_drawn_limit(0, 230) == 0             # 0 means draw everything


def test_final_table_rows_are_split_into_own_and_group_rows(tmp_path):
    from resscan import homscan_visualise as hv
    table = tmp_path / "homscan.tsv"
    table.write_text(
        "#AMR_Gene_Family\tARO\tRead_Count\tRPKPC\tTop_ARO\tAllocation_Proportions\n"
        "fam\tARO_1__geneA__H\t7\t0.0100\tARO_1__geneA__H\tARO_1__geneA__H:100.0%\n"
        "fam\tmultiple\t90\t0.2000\tARO_1\tARO_1:60.0%,ARO_2:40.0%\n")
    summary = hv.load_summary_table(str(table))
    own, groups = hv.final_rows_for(summary, "ARO_1__geneA__H")
    assert own["reads"] == "7" and [g["share"] for g in groups] == ["60.0%"]
    own, groups = hv.final_rows_for(summary, "ARO_2__geneB__H")
    assert own is None and groups[0]["group_size"] == 2


def test_varscan_block_view_window_and_site_evidence():
    ref = "ATG" * 20 + "CAA" + "GGC" * 20                 # codon 21 = CAA (Q) at reference bases 61-63
    que = "." * 30 + ref[30:90] + "." * (len(ref) - 90)
    block = {"header": "> read_1 | maps to | ARO_9__Test_PBP2x__V__E21Q | mutation: E21Q",
             "refseq": ref, "queseq": que, "raw_ref_line": "REFSEQ:    " + ref, "aro": "ARO_9__Test_PBP2x__V__E21Q"}
    view = vv.build_block_view(block, flank=10)
    assert (view["aligned_from"], view["aligned_to"]) == (31, 90)
    assert (view["window_ref_from"], view["window_ref_to"]) == (20, 100)
    assert view["marker_window"].index("^^^") == 60 - 20   # base 61 is column 60 of the reference
    assert view["sites"] == [{"label": "E21Q", "reference": "CAA", "read": "CAA", "reference_aa": "Q", "read_aa": "Q"}]


def test_index_page_links_only_to_pages_that_exist(tmp_path):
    from resscan import resscan_report_index as ri
    hom_dir = tmp_path / "s_homscan_html"
    hom_dir.mkdir()
    (hom_dir / "family_a.html").write_text("x")
    rows = [{"AMR_Gene_Family": "family A", "ARO": "ARO_1__geneA__H", "ARO_Name": "geneA", "Drug_Class": "x;y | x",
             "Read_Count": "10", "Fragment_Count": "6", "Lateral_Coverage_%": "50", "RPKPC": "0.5", "RPKG": "3",
             "Top_ARO": "ARO_1__geneA__H", "Allocation_Proportions": "ARO_1__geneA__H:100.0%"},
            {"AMR_Gene_Family": "family B", "ARO": "multiple", "ARO_Name": "geneB | geneC", "Drug_Class": "z",
             "Read_Count": "4", "Fragment_Count": "2", "Lateral_Coverage_%": "20", "RPKPC": "0.1", "RPKG": "1",
             "Top_ARO": "ARO_2", "Allocation_Proportions": "ARO_2:70.0%,ARO_3:30.0%"}]
    page = ri.build_index("s", rows, [], [], [], {}, str(hom_dir), str(tmp_path / "none"), "s_homscan_html",
                          "s_varscan_html", {"ARO_2": ("geneB", "family B")})
    assert '<a href="s_homscan_html/family_a.html">family A</a>' in page
    assert "family_b.html" not in page                    # no page on disk, so no link
    assert "geneB" in page and "group of 2" in page and "x; y" in page


def test_long_family_names_give_short_distinct_file_names():
    from resscan.report_common import sanitize_filename
    assert sanitize_filename("TEM beta-lactamase") == "tem_beta-lactamase.html"
    long_a = ";".join(["antibiotic efflux regulatory protein"] * 12) + ";variant A"
    long_b = ";".join(["antibiotic efflux regulatory protein"] * 12) + ";variant B"
    assert len(sanitize_filename(long_a).encode()) < 255
    assert sanitize_filename(long_a) != sanitize_filename(long_b)
