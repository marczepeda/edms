"""
Tests for edms.bio.pe — prime editing design/annotation logic.

Covers the pure sequence-math helpers (get_codons, get_codon_frames,
found_list_in_order) with hand-verified expectations, and the
pandas-DataFrame-oriented functions (find_enzyme_sites, pegRNA_outcome,
sensor_designer, merge, group_pe, shared_sequences, print_shared_sequences*)
using small synthetic DataFrames built from hand-verifiable DNA sequences.
"""
import pandas as pd
import pytest
from Bio.Seq import Seq

from edms.bio import pe


# --------------------------------------------------------------------------- #
# get_codons / get_codon_frames
# --------------------------------------------------------------------------- #

def test_get_codons_frame0():
    assert pe.get_codons("ATGCTGTGA") == ["ATG", "CTG", "TGA"]


def test_get_codons_frame1():
    # frame=1 drops the first base; trailing partial codon (1 nt) is dropped too
    assert pe.get_codons("ATGCTGTGA", frame=1) == ["TGC", "TGT"]


def test_get_codons_frame2():
    assert pe.get_codons("ATGCTGTGA", frame=2) == ["GCT", "GTG"]


def test_get_codons_empty_on_short_sequence():
    # len < 3 -> range(frame, len-2, 3) is empty
    assert pe.get_codons("AT") == []


def test_get_codon_frames():
    assert pe.get_codon_frames("ATGCTGTGA") == [
        ["ATG", "CTG", "TGA"],
        ["TGC", "TGT"],
        ["GCT", "GTG"],
    ]


# --------------------------------------------------------------------------- #
# found_list_in_order
# --------------------------------------------------------------------------- #

def test_found_list_in_order_found_at_start():
    assert pe.found_list_in_order(["A", "B", "C", "D"], ["A", "B"]) == 0


def test_found_list_in_order_found_in_middle():
    assert pe.found_list_in_order(["A", "B", "C", "D"], ["B", "C"]) == 1


def test_found_list_in_order_runs_past_end_returns_negative_one():
    assert pe.found_list_in_order(["A", "B", "C", "D"], ["C", "D", "E"]) == -1


def test_found_list_in_order_middle_mismatch_returns_negative_one():
    # 3-item sub_ls where the *middle* item is wrong is correctly rejected
    assert pe.found_list_in_order(["A", "B", "C", "D", "E"], ["C", "X", "E"]) == -1


def test_found_list_in_order_out_of_order_returns_negative_one():
    # 'C' then 'B' never occurs consecutively in that order in main_ls
    assert pe.found_list_in_order(["A", "B", "C", "D"], ["C", "B"]) == -1


def test_found_list_in_order_trailing_mismatch_returns_negative_one():
    # 'C','D' match but the trailing 'X' does not match 'E' -> not found
    assert pe.found_list_in_order(["A", "B", "C", "D", "E"], ["C", "D", "X"]) == -1


def test_found_list_in_order_single_item():
    assert pe.found_list_in_order(["A", "B", "C"], ["C"]) == 2


# --------------------------------------------------------------------------- #
# find_enzyme_sites
# --------------------------------------------------------------------------- #

def test_find_enzyme_sites_counts_fwd_and_rc_hits():
    # Esp3I recognition site: CGTCTC (fwd), GAGACG (rc)
    df = pd.DataFrame({
        "Oligonucleotide": [
            "AAACGTCTCAAA",  # one fwd hit
            "AAAAAAAAAA",    # no hits
            "AAAGAGACGAAA",  # one rc hit
        ]
    })
    out = pe.find_enzyme_sites(df, enzyme="Esp3I")
    assert list(out["Esp3I"]) == [1, 0, 1]
    assert out["Esp3I_fwd_i"].iloc[0] == [3]
    assert out["Esp3I_fwd_i"].iloc[1] == []
    assert out["Esp3I_rc_i"].iloc[2] == [3]


def test_find_enzyme_sites_counts_multiple_hits_in_one_row():
    df = pd.DataFrame({"Oligonucleotide": ["CGTCTCAAACGTCTC"]})
    out = pe.find_enzyme_sites(df, enzyme="Esp3I")
    assert out["Esp3I"].iloc[0] == 2
    assert out["Esp3I_fwd_i"].iloc[0] == [0, 9]


# --------------------------------------------------------------------------- #
# pegRNA_outcome
# --------------------------------------------------------------------------- #

def test_pegRNA_outcome_plus_strand_correct_geometry():
    ref = "ATCGGCTAACGTTAGCATGCTAGGGACCTGATCGTACGGA"
    spacer = ref[5:25]
    rc = str(Seq(ref).reverse_complement())
    pbs = rc[10:20]
    rtt = rc[20:30]
    edit_seq = ref[:20] + str(Seq(rtt + pbs).reverse_complement()) + ref[30:]

    df = pd.DataFrame({
        "Strand": ["+"],
        "Spacer_sequence": [spacer],
        "PBS_sequence": [pbs],
        "RTT_sequence": [rtt],
        "Reference_sequence": [ref],
        "Edit_sequence": [edit_seq],
    })
    out = pe.pegRNA_outcome(df, detailed=True)
    assert bool(out["Programmed_sequence_match"].iloc[0]) is True
    assert bool(out["Spacer_reference_match"].iloc[0]) is True
    assert bool(out["PBS_reference_match"].iloc[0]) is True
    assert bool(out["Edit_junction_match"].iloc[0]) is True
    assert out["Spacer_reference_matches"].iloc[0] == 1
    assert out["PBS_reference_matches"].iloc[0] == 1


def test_pegRNA_outcome_minus_strand_correct_geometry():
    ref = "ATCGGCTAACGTTAGCATGCTAGGGACCTGATCGTACGGA"
    rc = str(Seq(ref).reverse_complement())
    spacer = rc[5:25]
    rtt = ref[10:20]
    pbs = ref[20:30]
    edit_seq = ref[:10] + rtt + pbs + ref[30:]

    df = pd.DataFrame({
        "Strand": ["-"],
        "Spacer_sequence": [spacer],
        "PBS_sequence": [pbs],
        "RTT_sequence": [rtt],
        "Reference_sequence": [ref],
        "Edit_sequence": [edit_seq],
    })
    out = pe.pegRNA_outcome(df, detailed=True)
    assert bool(out["Programmed_sequence_match"].iloc[0]) is True


def test_pegRNA_outcome_detects_broken_pbs():
    ref = "ATCGGCTAACGTTAGCATGCTAGGGACCTGATCGTACGGA"
    rc = str(Seq(ref).reverse_complement())
    spacer = rc[5:25]
    rtt = ref[10:20]
    pbs = ref[20:30]
    edit_seq = ref[:10] + rtt + pbs + ref[30:]

    df = pd.DataFrame({
        "Strand": ["-"],
        "Spacer_sequence": [spacer],
        "PBS_sequence": ["AAAAAAAAAA"],  # corrupted PBS: won't match reference
        "RTT_sequence": [rtt],
        "Reference_sequence": [ref],
        "Edit_sequence": [edit_seq],
    })
    out = pe.pegRNA_outcome(df, detailed=True)
    assert bool(out["Programmed_sequence_match"].iloc[0]) is False
    assert bool(out["PBS_reference_match"].iloc[0]) is False


def test_pegRNA_outcome_missing_required_column_raises():
    df = pd.DataFrame({"Strand": ["+"]})
    with pytest.raises(ValueError):
        pe.pegRNA_outcome(df)


def test_pegRNA_outcome_not_detailed_only_adds_summary_col():
    ref = "ATCGGCTAACGTTAGCATGCTAGGGACCTGATCGTACGGA"
    rc = str(Seq(ref).reverse_complement())
    spacer = rc[5:25]
    rtt = ref[10:20]
    pbs = ref[20:30]
    edit_seq = ref[:10] + rtt + pbs + ref[30:]
    df = pd.DataFrame({
        "Strand": ["-"], "Spacer_sequence": [spacer], "PBS_sequence": [pbs],
        "RTT_sequence": [rtt], "Reference_sequence": [ref], "Edit_sequence": [edit_seq],
    })
    out = pe.pegRNA_outcome(df, detailed=False)
    assert "Programmed_sequence_match" in out.columns
    assert "Spacer_reference_match" not in out.columns


# --------------------------------------------------------------------------- #
# sensor_designer
# --------------------------------------------------------------------------- #

def test_sensor_designer_forward_plus_strand():
    ref = "ATCGGCTAACGTTAGCATGCTAGGGACCTGATCGTACGGA"
    spacer = ref[5:25]
    df = pd.DataFrame({"Spacer_sequence": [spacer], "Strand": ["+"], "Reference_sequence": [ref]})
    out = pe.sensor_designer(df, sensor_length=20, before_spacer=5, sensor_orientation="forward")
    assert out["Sensor_sequence"].iloc[0] == ref[0:20]


def test_sensor_designer_revcom_plus_strand():
    ref = "ATCGGCTAACGTTAGCATGCTAGGGACCTGATCGTACGGA"
    spacer = ref[5:25]
    df = pd.DataFrame({"Spacer_sequence": [spacer], "Strand": ["+"], "Reference_sequence": [ref]})
    out = pe.sensor_designer(df, sensor_length=20, before_spacer=5, sensor_orientation="revcom")
    assert out["Sensor_sequence"].iloc[0] == str(Seq(ref[0:20]).reverse_complement())


def test_sensor_designer_minus_strand():
    ref = "ATCGGCTAACGTTAGCATGCTAGGGACCTGATCGTACGGA"
    rc = str(Seq(ref).reverse_complement())
    spacer = rc[5:25]
    df = pd.DataFrame({"Spacer_sequence": [spacer], "Strand": ["-"], "Reference_sequence": [ref]})
    out = pe.sensor_designer(df, sensor_length=20, before_spacer=5, sensor_orientation="forward")
    assert out["Sensor_sequence"].iloc[0] == rc[0:20]


def test_sensor_designer_invalid_sensor_length_raises():
    df = pd.DataFrame({"Spacer_sequence": ["AAA"], "Strand": ["+"], "Reference_sequence": ["AAA"]})
    with pytest.raises(ValueError):
        pe.sensor_designer(df, sensor_length=0)


def test_sensor_designer_invalid_before_spacer_raises():
    df = pd.DataFrame({"Spacer_sequence": ["AAA"], "Strand": ["+"], "Reference_sequence": ["AAA"]})
    with pytest.raises(ValueError):
        pe.sensor_designer(df, before_spacer=0)


def test_sensor_designer_odd_sensor_length_gets_incremented():
    ref = "ATCGGCTAACGTTAGCATGCTAGGGACCTGATCGTACGGA"
    spacer = ref[5:25]
    df = pd.DataFrame({"Spacer_sequence": [spacer], "Strand": ["+"], "Reference_sequence": [ref]})
    out = pe.sensor_designer(df, sensor_length=19, before_spacer=5, sensor_orientation="forward")
    # 19 -> bumped up to 20 internally
    assert len(out["Sensor_sequence"].iloc[0]) == 20


# --------------------------------------------------------------------------- #
# merge
# --------------------------------------------------------------------------- #

def test_merge_joins_and_assigns_ngRNA_group_cycling():
    epegRNAs = pd.DataFrame({"pegRNA_number": [1, 2, 3], "Spacer_sequence": ["AAA", "BBB", "CCC"]})
    ngRNAs = pd.DataFrame({"pegRNA_number": [1, 1, 1, 2, 2, 3], "ngRNA_number": [1, 2, 3, 1, 2, 1]})
    out = pe.merge(epegRNAs, ngRNAs, ngRNAs_groups_max=2)
    assert list(out["pegRNA_number"]) == [1, 1, 1, 2, 2, 3]
    assert list(out["ngRNA_number"]) == [1, 2, 3, 1, 2, 1]
    # group = count % max + 1, starting count at 1 and incrementing per row for that pegRNA
    assert list(out["ngRNA_group"]) == [2, 1, 2, 2, 1, 2]


def test_merge_drops_ngRNAs_without_matching_pegRNA():
    epegRNAs = pd.DataFrame({"pegRNA_number": [1], "Spacer_sequence": ["AAA"]})
    ngRNAs = pd.DataFrame({"pegRNA_number": [1, 99], "ngRNA_number": [1, 1]})
    out = pe.merge(epegRNAs, ngRNAs)
    assert list(out["pegRNA_number"]) == [1]


# --------------------------------------------------------------------------- #
# group_pe
# --------------------------------------------------------------------------- #

def test_group_pe_pairs_shared_spacer_and_ngRNA():
    df = pd.DataFrame({
        "epegRNA": ["e1", "e2"],
        "ngRNA": ["n1", "n1"],
        "Spacer_sequence_epegRNA": ["SPACER1", "SPACER1"],
        "RTT_sequence": ["AAAAAAAAAA", "AAAAAAAAAT"],
        "PBS_sequence": ["GGGGG", "GGGGG"],
        "extra": ["x1", "x2"],
    })
    out = pe.group_pe(df, other_cols=["extra"])
    assert set(out["epegRNA"]) == {"e1", "e2"}
    assert set(out["epegRNA_compare"]) == {"e1", "e2"}
    # every pairing has itself excluded (epegRNA != epegRNA_compare)
    assert all(out["epegRNA"] != out["epegRNA_compare"])
    # match/mismatch score defaults: match=2, mismatch=-1; 9 matches + 1 mismatch => score=17
    # mismatches reported = len(RTT) - alignment.score = 10 - 17 = -7
    assert list(out["RTT_alignments_mismatches"]) == [-7, -7]
    assert list(out["PBS_lengths"]) == ["(5,5)", "(5,5)"]


def test_group_pe_no_shared_spacer_yields_empty():
    df = pd.DataFrame({
        "epegRNA": ["e1", "e2"],
        "ngRNA": ["n1", "n2"],
        "Spacer_sequence_epegRNA": ["SPACER1", "SPACER2"],
        "RTT_sequence": ["AAAAAAAAAA", "CCCCCCCCCC"],
        "PBS_sequence": ["GGGGG", "TTTTT"],
        "extra": ["x1", "x2"],
    })
    out = pe.group_pe(df, other_cols=["extra"])
    assert out.empty


# --------------------------------------------------------------------------- #
# shared_sequences
# --------------------------------------------------------------------------- #

def test_shared_sequences_reduces_to_unique_spacer_pbs_groups():
    df = pd.DataFrame({
        "Target_name": ["gene"] * 4,
        "pegRNA_number": [1, 2, 3, 4],
        "Strand": ["+"] * 4,
        "Edit": ["A10G", "A11G", "C20T", "C21T"],
        "Spacer_sequence": ["SP1", "SP1", "SP2", "SP2"],
        "PBS_sequence": ["PBS1", "PBS1", "PBS2", "PBS2"],
        "RTT_length": [10, 12, 10, 12],
    })
    out = pe.shared_sequences(df, hist_plot=False)
    assert len(out) == 2
    row0 = out[out["Spacer_sequence"] == "SP1"].iloc[0]
    assert row0["pegRNA_numbers"] == [1, 2]
    assert row0["Edits"] == ["A10G", "A11G"]
    assert row0["AA_numbers"] == [10, 11]
    assert row0["AA_numbers_min"] == 10
    assert row0["AA_numbers_max"] == 11
    assert bool(row0["AA_numbers_continuous"]) is True

    row1 = out[out["Spacer_sequence"] == "SP2"].iloc[0]
    assert row1["AA_numbers_min"] == 20
    assert row1["AA_numbers_max"] == 21


def test_shared_sequences_detects_noncontinuous_aa_numbers():
    df = pd.DataFrame({
        "Target_name": ["gene"] * 2,
        "pegRNA_number": [1, 2],
        "Strand": ["+"] * 2,
        "Edit": ["A10G", "A15G"],  # gap between 10 and 15
        "Spacer_sequence": ["SP1", "SP1"],
        "PBS_sequence": ["PBS1", "PBS1"],
        "RTT_length": [10, 10],
    })
    out = pe.shared_sequences(df, hist_plot=False)
    assert bool(out["AA_numbers_continuous"].iloc[0]) is False


# --------------------------------------------------------------------------- #
# print_shared_sequences / print_shared_sequences_mutant
# --------------------------------------------------------------------------- #

def test_print_shared_sequences_outputs_expected_text(capsys):
    dc = {
        "A": pd.DataFrame({"Spacer_sequence": ["SP1", "SP2"], "PBS_sequence": ["PBS1", "PBS2"]}),
        "B": pd.DataFrame({"Spacer_sequence": ["SP3", "SP4"], "PBS_sequence": ["PBS3", "PBS4"]}),
    }
    pe.print_shared_sequences(dc)
    out = capsys.readouterr().out
    assert "A_spacer" in out
    assert "B_PBS" in out
    assert "SP1" in out and "PBS3" in out


def test_print_shared_sequences_mutant_outputs_expected_text(capsys):
    dc = {"A": pd.DataFrame({"Spacer_sequence": ["SP1"], "PBS_sequence": ["PBS1"], "Priority_mut": ["M1"]})}
    pe.print_shared_sequences_mutant(dc)
    out = capsys.readouterr().out
    assert "A_mutant" in out
    assert "M1" in out
