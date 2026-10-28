"""
Tests for edms.bio.pe — prime editing design/annotation logic.

Covers the pure sequence-math helpers (get_codons, get_codon_frames,
found_list_in_order) with hand-verified expectations, and the
pandas-DataFrame-oriented functions (find_enzyme_sites, pegRNA_outcome,
sensor_designer, merge, group_pe, shared_sequences, print_shared_sequences*)
using small synthetic DataFrames built from hand-verifiable DNA sequences.
"""
import numpy as np
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


# --------------------------------------------------------------------------- #
# pegRNA_molecule_key / drop_duplicate_pegRNAs
# --------------------------------------------------------------------------- #

def _peg(spacer="ACGTACGTACGTACGTACGT", rtt="AAACCC", pbs="GGGTTT", edit="R1C", **kw):
    row = {"Spacer_sequence": spacer, "Scaffold_sequence": "GTTTAAGAGC",
           "RTT_sequence": rtt, "PBS_sequence": pbs, "Edit": edit}
    row.update(kw)
    return row


class TestPegRNAMoleculeKey:
    def test_identical_sequences_give_identical_keys(self):
        df = pd.DataFrame([_peg(edit="R1C"), _peg(edit="R1H")])
        key = pe.pegRNA_molecule_key(df)
        assert key.iloc[0] == key.iloc[1]  # same molecule, different edit label

    def test_a_different_rtt_gives_a_different_key(self):
        df = pd.DataFrame([_peg(rtt="AAACCC"), _peg(rtt="AAACCG")])
        key = pe.pegRNA_molecule_key(df)
        assert key.iloc[0] != key.iloc[1]

    def test_case_is_ignored(self):
        df = pd.DataFrame([_peg(rtt="aaaccc"), _peg(rtt="AAACCC")])
        key = pe.pegRNA_molecule_key(df)
        assert key.iloc[0] == key.iloc[1]

    def test_linker_is_included_when_present(self):
        df = pd.DataFrame([_peg(Linker_sequence="AAAA"), _peg(Linker_sequence="TTTT")])
        key = pe.pegRNA_molecule_key(df)
        assert key.iloc[0] != key.iloc[1]

    def test_missing_sequence_columns_raises(self):
        with pytest.raises(ValueError, match="none of the sequence columns"):
            pe.pegRNA_molecule_key(pd.DataFrame({"Edit": ["R1C"]}))


class TestDropDuplicatePegRNAs:
    def test_removes_duplicates_within_the_frame(self):
        df = pd.DataFrame([_peg(edit="R1C"), _peg(edit="R1H"), _peg(rtt="AAACCG", edit="R1K")])
        out = pe.drop_duplicate_pegRNAs(df)
        assert len(out) == 2
        assert out.iloc[0]["Edit"] == "R1C"  # first occurrence is kept

    def test_removes_rows_already_present_in_existing(self):
        existing = pd.DataFrame([_peg(rtt="AAACCC")])
        new = pd.DataFrame([_peg(rtt="AAACCC", edit="R1H"), _peg(rtt="AAACCG", edit="R1K")])
        out = pe.drop_duplicate_pegRNAs(new, existing=existing)
        assert list(out["Edit"]) == ["R1K"]

    def test_keeps_everything_when_nothing_collides(self):
        existing = pd.DataFrame([_peg(rtt="TTTTTT")])
        new = pd.DataFrame([_peg(rtt="AAACCC"), _peg(rtt="AAACCG")])
        assert len(pe.drop_duplicate_pegRNAs(new, existing=existing)) == 2

    def test_empty_input_is_returned_unchanged(self):
        empty = pd.DataFrame(columns=["Spacer_sequence", "Scaffold_sequence",
                                      "RTT_sequence", "PBS_sequence", "Edit"])
        assert len(pe.drop_duplicate_pegRNAs(empty)) == 0

    def test_empty_existing_is_ignored(self):
        new = pd.DataFrame([_peg(rtt="AAACCC")])
        empty = pd.DataFrame(columns=new.columns)
        assert len(pe.drop_duplicate_pegRNAs(new, existing=empty)) == 1

    def test_index_is_reset(self):
        df = pd.DataFrame([_peg(edit="R1C"), _peg(edit="R1H"), _peg(rtt="AAACCG", edit="R1K")])
        out = pe.drop_duplicate_pegRNAs(df)
        assert list(out.index) == list(range(len(out)))


# --------------------------------------------------------------------------- #
# codon_nt_difference / sort_pegRNAs
# --------------------------------------------------------------------------- #

class TestCodonNtDifference:
    def test_identical_codons(self):
        assert pe.codon_nt_difference("GTG", "GTG") == 0

    def test_counts_each_differing_position(self):
        assert pe.codon_nt_difference("CAC", "GGC") == 2
        assert pe.codon_nt_difference("CAC", "GGA") == 3

    def test_case_insensitive(self):
        assert pe.codon_nt_difference("cac", "GGC") == 2

    def test_missing_or_short_codon_is_nan(self):
        assert np.isnan(pe.codon_nt_difference("", "GGC"))
        assert np.isnan(pe.codon_nt_difference("GT", "GTG"))
        assert np.isnan(pe.codon_nt_difference(np.nan, "GTG"))


def _sortable(edit, spacer, pbs, ref_codon, edit_codon, usage, number,
              silent_ref="AAA", silent_edit="AAG", silent_usage=0.5, aa=1):
    return {"Edit": edit, "AA_number": aa, "Spacer_sequence": spacer, "PBS_sequence": pbs,
            "Target_reference_codon": ref_codon, "Target_edit_codon": edit_codon,
            "Target_edit_codon_usage": usage, "Silent_reference_codon": silent_ref,
            "Silent_edit_codon": silent_edit, "Silent_edit_codon_usage": silent_usage,
            "pegRNA_number": number}


class TestSortPegRNAs:
    def test_adds_the_nucleotide_difference_columns(self):
        df = pd.DataFrame([_sortable("Y6V", "S1", "P1", "TAC", "GTG", 0.47, 1)])
        out = pe.sort_pegRNAs(df)
        assert out.loc[0, "Target_codon_nt_difference"] == 3   # TAC -> GTG
        assert out.loc[0, "Silent_codon_nt_difference"] == 1   # AAA -> AAG

    def test_orders_by_nucleotide_difference_then_usage(self):
        # shuffled on input; TAC->GTC differs at 2, the rest at 3.
        # pegRNA_number deliberately contradicts the usage order, so falling back to it
        # would give GTT, GTA, GTG and fail -- only a usage-based tiebreak passes.
        df = pd.DataFrame([
            _sortable("Y6V", "S1", "P1", "TAC", "GTC", 0.24, 1),
            _sortable("Y6V", "S1", "P1", "TAC", "GTA", 0.11, 8),
            _sortable("Y6V", "S1", "P1", "TAC", "GTG", 0.47, 9),
            _sortable("Y6V", "S1", "P1", "TAC", "GTT", 0.18, 7),
        ])
        out = pe.sort_pegRNAs(df)
        assert list(out["Target_edit_codon"]) == ["GTG", "GTT", "GTA", "GTC"]

    def test_keeps_each_edit_spacer_pbs_group_contiguous(self):
        df = pd.DataFrame([
            _sortable("Y6V", "S1", "P1", "TAC", "GTG", 0.47, 1),
            _sortable("Y6V", "S2", "P2", "TAC", "GTT", 0.18, 2),
            _sortable("Y6V", "S1", "P1", "TAC", "GTC", 0.24, 3),
            _sortable("Y6V", "S2", "P2", "TAC", "GTA", 0.11, 4),
        ])
        out = pe.sort_pegRNAs(df)
        codes = out[["Edit", "Spacer_sequence", "PBS_sequence"]].astype(str).agg("|".join, axis=1)
        runs = (codes != codes.shift()).cumsum().nunique()
        assert runs == codes.nunique()   # no group is split across the file

    def test_a_revised_pegrna_lands_beside_its_group(self):
        # the revised row (different RTT, fewest differences) must sort into its own group
        df = pd.DataFrame([
            _sortable("A1B", "Sx", "Px", "AAA", "CCC", 0.5, 10),
            _sortable("Y6V", "S1", "P1", "TAC", "GTG", 0.47, 1),
            _sortable("Y6V", "S1", "P1", "TAC", "GTC", 0.24, 99),   # the revised one
            _sortable("A1B", "Sx", "Px", "AAA", "CCG", 0.4, 11),
        ])
        out = pe.sort_pegRNAs(df).reset_index(drop=True)
        y6v = out.index[out["Edit"] == "Y6V"].tolist()
        assert y6v == list(range(min(y6v), min(y6v) + len(y6v)))   # contiguous
        assert out.loc[y6v[-1], "pegRNA_number"] == 99             # ranked last in its group

    def test_amino_acid_order_takes_precedence(self):
        df = pd.DataFrame([
            _sortable("Q9X", "S1", "P1", "AAA", "CCC", 0.5, 1, aa=9),
            _sortable("A1X", "S1", "P1", "AAA", "CCC", 0.5, 2, aa=1),
        ])
        out = pe.sort_pegRNAs(df)
        assert list(out["AA_number"]) == [1, 9]

    def test_empty_frame_is_returned_unchanged(self):
        empty = pd.DataFrame(columns=["Edit", "Spacer_sequence", "PBS_sequence"])
        assert len(pe.sort_pegRNAs(empty)) == 0

    def test_tolerates_missing_optional_columns(self):
        df = pd.DataFrame([{"Edit": "B", "pegRNA_number": 2}, {"Edit": "A", "pegRNA_number": 1}])
        out = pe.sort_pegRNAs(df)
        assert list(out["Edit"]) == ["A", "B"]

    def test_index_is_reset(self):
        df = pd.DataFrame([
            _sortable("Y6V", "S1", "P1", "TAC", "GTC", 0.24, 2),
            _sortable("Y6V", "S1", "P1", "TAC", "GTG", 0.47, 1),
        ])
        assert list(pe.sort_pegRNAs(df).index) == [0, 1]


# --------------------------------------------------------------------------- #
# refresh_codon_annotations
# --------------------------------------------------------------------------- #

# in-frame reference; codon 3 ([9,12)) stands in for the intended edit and
# codon 5 ([15,18)) for the silent mutation
REFRESH_REF = "ATGGCTCTGACCAAGCATGAATTCGTGCGT"


def _annotated(edit_sequence):
    """A row whose codon annotations describe REFRESH_REF, not necessarily edit_sequence."""
    return {
        "Reference_sequence": REFRESH_REF,
        "Edit_sequence": edit_sequence,
        "Target_reference_codon": "ACC", "Target_reference_AA": "T",
        "Target_edit_codon": "ACC", "Target_edit_AA": "T",
        "Target_edit_codon_usage": 0.36, "Target_edit_codon_index": 9,
        "Silent_reference_codon": "CAT", "Silent_reference_AA": "H",
        "Silent_edit_codon": "CAT", "Silent_edit_AA": "H",
        "Silent_edit_codon_usage": 0.41, "Silent_edit_codon_index": 15,
    }


class TestRefreshCodonAnnotations:
    def test_rereads_the_edit_codon_from_the_rewritten_sequence(self):
        # codon 3 rewritten ACC -> ACG (still Thr)
        edit = REFRESH_REF[:9] + "ACG" + REFRESH_REF[12:]
        out = pe.refresh_codon_annotations(pd.DataFrame([_annotated(edit)]))
        assert out.loc[0, "Target_edit_codon"] == "ACG"
        assert out.loc[0, "Target_edit_AA"] == "T"
        assert out.loc[0, "Target_edit_codon_usage"] == 0.12   # ACG human usage

    def test_rereads_the_silent_codon(self):
        edit = REFRESH_REF[:15] + "CAC" + REFRESH_REF[18:]     # CAT -> CAC, still His
        out = pe.refresh_codon_annotations(pd.DataFrame([_annotated(edit)]))
        assert out.loc[0, "Silent_edit_codon"] == "CAC"
        assert out.loc[0, "Silent_edit_AA"] == "H"
        assert out.loc[0, "Silent_edit_codon_usage"] == 0.59

    def test_leaves_the_reference_columns_alone(self):
        edit = REFRESH_REF[:9] + "ACG" + REFRESH_REF[12:]
        out = pe.refresh_codon_annotations(pd.DataFrame([_annotated(edit)]))
        assert out.loc[0, "Target_reference_codon"] == "ACC"   # Reference_sequence never changes
        assert out.loc[0, "Silent_reference_codon"] == "CAT"

    def test_recomputes_the_change_indexes(self):
        edit = REFRESH_REF[:9] + "ACG" + REFRESH_REF[12:]      # one base differs, at 11
        out = pe.refresh_codon_annotations(pd.DataFrame([_annotated(edit)]))
        assert out.loc[0, "Edit_sequence_change_indexes"] == "11"
        assert out.loc[0, "Reference_sequence_change_indexes"] == "11"

    def test_unchanged_sequence_leaves_everything_as_it_was(self):
        out = pe.refresh_codon_annotations(pd.DataFrame([_annotated(REFRESH_REF)]))
        assert out.loc[0, "Target_edit_codon"] == "ACC"
        assert out.loc[0, "Silent_edit_codon"] == "CAT"
        assert out.loc[0, "Edit_sequence_change_indexes"] == ""

    def test_row_without_a_recorded_index_is_left_alone(self):
        row = _annotated(REFRESH_REF)
        row["Silent_edit_codon_index"] = np.nan
        row["Silent_edit_codon"] = "GGG"          # deliberately wrong, but unindexed
        out = pe.refresh_codon_annotations(pd.DataFrame([row]))
        assert out.loc[0, "Silent_edit_codon"] == "GGG"

    def test_index_running_off_the_end_blanks_the_codon(self):
        row = _annotated(REFRESH_REF)
        row["Silent_edit_codon_index"] = len(REFRESH_REF) - 1   # only 1 base available
        out = pe.refresh_codon_annotations(pd.DataFrame([row]))
        assert out.loc[0, "Silent_edit_codon"] == ""
        assert out.loc[0, "Silent_edit_AA"] == ""

    def test_empty_frame_is_returned_unchanged(self):
        assert len(pe.refresh_codon_annotations(pd.DataFrame(columns=["Edit_sequence"]))) == 0

    def test_frame_without_codon_columns_is_tolerated(self):
        df = pd.DataFrame([{"Edit_sequence": REFRESH_REF, "Reference_sequence": REFRESH_REF}])
        out = pe.refresh_codon_annotations(df)
        assert out.loc[0, "Edit_sequence_change_indexes"] == ""

    def test_refresh_is_idempotent(self):
        edit = REFRESH_REF[:9] + "ACG" + REFRESH_REF[12:]
        once = pe.refresh_codon_annotations(pd.DataFrame([_annotated(edit)]))
        twice = pe.refresh_codon_annotations(once.copy())
        pd.testing.assert_frame_equal(once, twice)


# --------------------------------------------------------------------------- #
# merge_pegRNAs_ngRNAs
# --------------------------------------------------------------------------- #

MERGE_ALLELE = "ATGGCTCTGACCAAGCATGAATTCGTGCGTAGCATCCCAGATAACCAGTAC"


def _peg_row(number, allele=MERGE_ALLELE, **kw):
    row = {"pegRNA_number": number, "Spacer_sequence": "AAACCCGGGTTTAAACCCGG",
           "Strand": "+", "Edit_sequence": allele, "Reference_sequence": MERGE_ALLELE,
           "Target_name": "demo", "Annotation": "PAM_disrupted", "Extension_sequence": "ACGT"}
    row.update(kw)
    return row


def _ng_row(number, spacer, allele=MERGE_ALLELE, **kw):
    row = {"pegRNA_number": number, "ngRNA_number": number * 10,
           "Spacer_sequence": spacer, "Strand": "-", "Edit_sequence": allele,
           "Reference_sequence": MERGE_ALLELE, "Target_name": "demo",
           "Annotation": "PE3b-seed", "Extension_sequence": ""}
    row.update(kw)
    return row


class TestMergePegRNAsNgRNAs:
    def test_one_row_per_pegrna_ngrna_pair(self):
        peg = pd.DataFrame([_peg_row(1)])
        ng = pd.DataFrame([_ng_row(1, MERGE_ALLELE[3:23]), _ng_row(1, MERGE_ALLELE[6:26])])
        out = pe.merge_pegRNAs_ngRNAs(peg, ng)
        assert len(out) == 2                       # the pegRNA is repeated per ngRNA
        assert out['pegRNA_number'].tolist() == [1, 1]

    def test_no_ngrna_group_is_assigned(self):
        peg = pd.DataFrame([_peg_row(1)])
        ng = pd.DataFrame([_ng_row(1, MERGE_ALLELE[3:23])])
        assert 'ngRNA_group' not in pe.merge_pegRNAs_ngRNAs(peg, ng).columns

    def test_columns_that_differ_get_the_prefix(self):
        peg = pd.DataFrame([_peg_row(1)])
        ng = pd.DataFrame([_ng_row(1, MERGE_ALLELE[3:23])])
        out = pe.merge_pegRNAs_ngRNAs(peg, ng)
        assert 'ngRNA_Spacer_sequence' in out.columns and 'ngRNA_Strand' in out.columns
        assert out.loc[0, 'Spacer_sequence'] == peg.loc[0, 'Spacer_sequence']   # pegRNA's kept

    def test_identical_columns_are_dropped(self):
        peg = pd.DataFrame([_peg_row(1)])
        ng = pd.DataFrame([_ng_row(1, MERGE_ALLELE[3:23])])
        out = pe.merge_pegRNAs_ngRNAs(peg, ng)
        for col in ('Reference_sequence', 'Target_name'):
            assert col in out.columns and f'ngRNA_{col}' not in out.columns

    def test_columns_empty_on_the_ngrna_side_are_dropped(self):
        peg = pd.DataFrame([_peg_row(1)])
        ng = pd.DataFrame([_ng_row(1, MERGE_ALLELE[3:23])])          # Extension_sequence is ''
        out = pe.merge_pegRNAs_ngRNAs(peg, ng)
        assert 'ngRNA_Extension_sequence' not in out.columns

    def test_no_pandas_suffixes_survive(self):
        peg = pd.DataFrame([_peg_row(1)])
        ng = pd.DataFrame([_ng_row(1, MERGE_ALLELE[3:23])])
        out = pe.merge_pegRNAs_ngRNAs(peg, ng)
        assert not [c for c in out.columns if c.endswith(('_x', '_y'))]

    def test_pairs_the_ngrna_cannot_bind_are_dropped(self):
        peg = pd.DataFrame([_peg_row(1)])
        ng = pd.DataFrame([_ng_row(1, MERGE_ALLELE[3:23]),
                           _ng_row(1, 'TTTTTTTTTTTTTTTTTTTT')])      # absent from the allele
        out = pe.merge_pegRNAs_ngRNAs(peg, ng)
        assert len(out) == 1
        assert out.loc[0, 'ngRNA_Spacer_sequence'] == MERGE_ALLELE[3:23]

    def test_a_reverse_complement_spacer_still_counts_as_bindable(self):
        rc = str(Seq(MERGE_ALLELE[3:23]).reverse_complement())
        peg = pd.DataFrame([_peg_row(1)])
        ng = pd.DataFrame([_ng_row(1, rc)])
        assert len(pe.merge_pegRNAs_ngRNAs(peg, ng)) == 1

    def test_grna_type_is_dropped_from_both_sides(self):
        peg = pd.DataFrame([_peg_row(1, gRNA_type='pegRNA')])
        ng = pd.DataFrame([_ng_row(1, MERGE_ALLELE[3:23], gRNA_type='ngRNA')])
        out = pe.merge_pegRNAs_ngRNAs(peg, ng)
        assert 'gRNA_type' not in out.columns and 'ngRNA_gRNA_type' not in out.columns
        assert len(out) == 1                       # the merge itself is unaffected

    def test_missing_join_key_raises(self):
        peg = pd.DataFrame([{"Spacer_sequence": "A"}])
        ng = pd.DataFrame([_ng_row(1, MERGE_ALLELE[3:23])])
        with pytest.raises(ValueError, match="pegRNA_number"):
            pe.merge_pegRNAs_ngRNAs(peg, ng)


# --------------------------------------------------------------------------- #
# refresh_ngRNA_spacers
# --------------------------------------------------------------------------- #

class TestRefreshNgRNASpacers:
    def _frame(self, spacer, strand='+'):
        return pd.DataFrame([{"pegRNA_number": 1, "Spacer_sequence": spacer, "Strand": strand,
                              "Annotation": "PE3", "Edit_sequence": MERGE_ALLELE}])

    def test_recuts_a_spacer_whose_allele_was_rewritten(self):
        after = MERGE_ALLELE[:10] + 'ttt' + MERGE_ALLELE[13:]
        out = pe.refresh_ngRNA_spacers(self._frame(MERGE_ALLELE[5:25]), {1: MERGE_ALLELE}, {1: after})
        assert out.loc[0, 'Spacer_sequence'] == after[5:25]

    def test_updates_the_stored_allele_too(self):
        after = MERGE_ALLELE[:10] + 'ttt' + MERGE_ALLELE[13:]
        out = pe.refresh_ngRNA_spacers(self._frame(MERGE_ALLELE[5:25]), {1: MERGE_ALLELE}, {1: after})
        assert out.loc[0, 'Edit_sequence'] == after

    def test_a_reverse_complement_spacer_is_recut_in_its_own_orientation(self):
        rc = str(Seq(MERGE_ALLELE[5:25]).reverse_complement())
        after = MERGE_ALLELE[:10] + 'ttt' + MERGE_ALLELE[13:]
        out = pe.refresh_ngRNA_spacers(self._frame(rc, strand='-'), {1: MERGE_ALLELE}, {1: after})
        assert out.loc[0, 'Spacer_sequence'] == str(Seq(after[5:25]).reverse_complement())

    def test_annotation_becomes_pe3b_when_an_edit_lands_in_the_spacer(self):
        after = MERGE_ALLELE[:10] + 'ttt' + MERGE_ALLELE[13:]
        out = pe.refresh_ngRNA_spacers(self._frame(MERGE_ALLELE[5:25]), {1: MERGE_ALLELE}, {1: after})
        assert out.loc[0, 'Annotation'].startswith('PE3b')

    def test_a_spacer_outside_the_change_keeps_its_sequence_but_gets_the_new_allele(self):
        after = MERGE_ALLELE[:45] + 'ttt' + MERGE_ALLELE[48:]   # far from the spacer
        out = pe.refresh_ngRNA_spacers(self._frame(MERGE_ALLELE[0:20]), {1: MERGE_ALLELE}, {1: after})
        assert out.loc[0, 'Spacer_sequence'] == MERGE_ALLELE[0:20]   # spacer untouched
        assert out.loc[0, 'Edit_sequence'] == after                  # but the allele moved

    def test_unchanged_allele_is_a_no_op(self):
        before = self._frame(MERGE_ALLELE[5:25])
        out = pe.refresh_ngRNA_spacers(before, {1: MERGE_ALLELE}, {1: MERGE_ALLELE})
        assert out.loc[0, 'Spacer_sequence'] == MERGE_ALLELE[5:25]

    def test_pegrna_absent_from_the_maps_is_skipped(self):
        out = pe.refresh_ngRNA_spacers(self._frame(MERGE_ALLELE[5:25]), {}, {})
        assert out.loc[0, 'Spacer_sequence'] == MERGE_ALLELE[5:25]

    def test_alleles_of_different_length_are_skipped(self):
        # the swap only substitutes, so a length change means these cannot be aligned
        out = pe.refresh_ngRNA_spacers(self._frame(MERGE_ALLELE[5:25]),
                                       {1: MERGE_ALLELE}, {1: MERGE_ALLELE[:-3]})
        assert out.loc[0, 'Spacer_sequence'] == MERGE_ALLELE[5:25]

    def test_empty_frame_is_returned_unchanged(self):
        assert len(pe.refresh_ngRNA_spacers(pd.DataFrame(columns=['pegRNA_number']), {}, {})) == 0
