"""
Tests for edms.bio.clone — molecular cloning / GG-cloning oligo design, UMI generation,
PCR master-mix calculations, PCR simulation, and off-target alignment.
"""
import numpy as np
import pandas as pd
import pytest
from Bio.Seq import Seq

from edms.bio import clone


# ---------------------------------------------------------------------------
# ord_form()
# ---------------------------------------------------------------------------

class TestOrdForm:
    def test_renames_and_computes_scale_bp(self):
        df = pd.DataFrame({
            "oid_top": ["a1", "a2"],
            "oseq_top": ["A" * 59, "A" * 60],
        })
        out = clone.ord_form(df=df, id="id", seq="seq", suf="_top", pre="o")
        assert list(out.columns) == ["Oligo Name", "Sequence", "Scale (µmol)", "bp"]
        assert out["Oligo Name"].tolist() == ["a1", "a2"]
        # boundary: len < 60 -> 0.025 scale; len >= 60 -> 0.05 scale
        assert out["Scale (µmol)"].tolist() == [0.025, 0.05]
        assert out["bp"].tolist() == [59, 60]


# ---------------------------------------------------------------------------
# tb()
# ---------------------------------------------------------------------------

class TestTb:
    def test_tG_adds_5prime_G_when_missing(self):
        df = pd.DataFrame({"ID": [1], "Spacer": ["ATCGA"]})
        out = clone.tb(df=df, id="ID", seq="Spacer", t5="CACC", t3="GTTT",
                        b5="AAAC", b3="GGGG", tG=True, pre="o")
        assert out["oSpacer_top"].iloc[0] == "CACC" + "G" + "ATCGA" + "GTTT"
        expected_bot = "AAAC" + str(Seq("GATCGA").reverse_complement()) + "GGGG"
        assert out["oSpacer_bot"].iloc[0] == expected_bot
        assert out["oID_top"].iloc[0] == "o1_top"
        assert out["oID_bot"].iloc[0] == "o1_bot"

    def test_tG_does_not_add_G_when_already_present(self):
        df = pd.DataFrame({"ID": [1], "Spacer": ["GATCGA"]})
        out = clone.tb(df=df, id="ID", seq="Spacer", t5="CACC", t3="", b5="AAAC", b3="", tG=True, pre="o")
        assert out["oSpacer_top"].iloc[0] == "CACC" + "GATCGA"
        expected_bot = "AAAC" + str(Seq("GATCGA").reverse_complement())
        assert out["oSpacer_bot"].iloc[0] == expected_bot

    def test_tG_false_never_adds_G(self):
        df = pd.DataFrame({"ID": [1], "Spacer": ["ATCGA"]})
        out = clone.tb(df=df, id="ID", seq="Spacer", t5="CACC", t3="", b5="AAAC", b3="", tG=False, pre="o")
        assert out["oSpacer_top"].iloc[0] == "CACC" + "ATCGA"
        expected_bot = "AAAC" + str(Seq("ATCGA").reverse_complement())
        assert out["oSpacer_bot"].iloc[0] == expected_bot

    def test_multiple_rows(self):
        df = pd.DataFrame({"ID": [1, 2], "Spacer": ["ATCGA", "GATCC"]})
        out = clone.tb(df=df, id="ID", seq="Spacer", t5="", t3="", b5="", b3="", tG=False, pre="o")
        assert out["oSpacer_top"].tolist() == ["ATCGA", "GATCC"]
        assert out["oSpacer_bot"].tolist() == [str(Seq("ATCGA").reverse_complement()),
                                                str(Seq("GATCC").reverse_complement())]


# ---------------------------------------------------------------------------
# Scaffold options
# ---------------------------------------------------------------------------

FLIP_AND_EXTEND_TOP = "TATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC"
FLIP_AND_EXTEND_PEG_BOT = "GGTGCCACTTGGGTGGGATGCTCGCAGGCATTCAGCCAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCATAGCTCTTAAAC"
FLIP_AND_EXTEND_NG_BOT = "AAAAGCACCGACTCGGTGCCACTTGGGTGGGATGCTCGCAGGCATTCAGCCAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCATAGCTCTTAAAC"


class TestScaffoldOverhangs:
    def test_flip_and_extend_overhangs_are_unchanged(self):
        # These were hardcoded before the overhangs were taken from the scaffold sequence.
        assert clone.SGRNA_SCAFFOLD_OVERHANGS["SpCas9 Flip and Extend"] == {"t3": "", "b5": "AAAC"}
        assert clone.PEGRNA_SCAFFOLD_OVERHANGS["SpCas9 Flip and Extend"] == {
            "spacer_t3": "GTTTAAGAGC", "spacer_b5": "", "extension_t5": "", "extension_b3": "GCACCGACTC"}

    def test_mod_7_shares_the_sgRNA_overhangs(self):
        assert clone.SGRNA_SCAFFOLD_OVERHANGS["SpCas9 mod_7"] == clone.SGRNA_SCAFFOLD_OVERHANGS["SpCas9 Flip and Extend"]

    def test_mod_7_differs_only_in_the_spacer_t3(self):
        mod_7 = clone.PEGRNA_SCAFFOLD_OVERHANGS["SpCas9 mod_7"]
        assert mod_7 == {"spacer_t3": "GTTTCAGAGC", "spacer_b5": "", "extension_t5": "", "extension_b3": "GCACCGACTC"}

    def test_overhangs_are_the_ends_of_the_scaffold(self):
        for scaffold, sequence in clone.SCAFFOLD_SEQUENCES.items():
            assert clone.PEGRNA_SCAFFOLD_OVERHANGS[scaffold]["spacer_t3"] == sequence[:10]
            assert clone.PEGRNA_SCAFFOLD_OVERHANGS[scaffold]["extension_b3"] == str(Seq(sequence[-10:]).reverse_complement())
            assert clone.SGRNA_SCAFFOLD_OVERHANGS[scaffold]["b5"] == str(Seq(sequence[:4]).reverse_complement())


class TestGetScaffoldOligos:
    def test_flip_and_extend_pegRNA_oligos_are_unchanged(self):
        out = clone.get_scaffold_oligos(scaffold="SpCas9 Flip and Extend", pre="pegRNA")
        assert out["Oligo Name"].tolist() == ["pegRNA_scaffold_top", "pegRNA_scaffold_bot"]
        assert out["Sequence"].tolist() == [FLIP_AND_EXTEND_TOP, FLIP_AND_EXTEND_PEG_BOT]
        assert out["bp"].tolist() == [99, 99]

    def test_flip_and_extend_ngRNA_oligos_are_unchanged(self):
        out = clone.get_scaffold_oligos(scaffold="SpCas9 Flip and Extend", pre="ngRNA")
        assert out["Sequence"].tolist() == [FLIP_AND_EXTEND_TOP, FLIP_AND_EXTEND_NG_BOT]
        assert out["bp"].tolist() == [99, 113]

    def test_mod_7_oligos_follow_its_own_sequence(self):
        sequence = clone.SCAFFOLD_SEQUENCES["SpCas9 mod_7"]
        peg = clone.get_scaffold_oligos(scaffold="SpCas9 mod_7", pre="pegRNA")
        ng = clone.get_scaffold_oligos(scaffold="SpCas9 mod_7", pre="ngRNA")
        assert peg["Sequence"].tolist() == [sequence[10:], str(Seq(sequence[:-10]).reverse_complement())]
        assert ng["Sequence"].tolist() == [sequence[10:], "AAAA" + str(Seq(sequence).reverse_complement())]
        assert peg["bp"].tolist() == [66, 66]  # shorter scaffold than Flip and Extend

    def test_scaffold_name_is_case_insensitive(self):
        assert clone.get_scaffold_oligos(scaffold="spcas9 mod_7", pre="pegRNA")["Sequence"].tolist() == \
               clone.get_scaffold_oligos(scaffold="SpCas9 mod_7", pre="pegRNA")["Sequence"].tolist()

    def test_unrecognized_scaffold_raises(self):
        with pytest.raises(ValueError, match="not an option"):
            clone.get_scaffold_oligos(scaffold="SpCas9", pre="pegRNA")


class TestScaffoldChoiceReachesTheOrder:
    def test_epegRNAs_orders_the_chosen_scaffold(self):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["GTCGATCGAT"], "RTT_sequence": ["AAAA"],
                           "PBS_sequence": ["CCCC"], "Linker_sequence": ["GGGG"]})
        out = clone.epegRNAs(df=df, id="ID", order=True, order_scaffold=True, scaffold="SpCas9 mod_7")
        ordered = out.loc[out["Oligo Name"] == "pegRNA_scaffold_top", "Sequence"].iloc[0]
        assert ordered == clone.SCAFFOLD_SEQUENCES["SpCas9 mod_7"][10:]

    def test_ngRNAs_orders_the_chosen_scaffold(self):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["GTCGATCGAT"]})
        out = clone.ngRNAs(df=df, id="ID", order=True, order_scaffold=True, vector_scaffold=False, scaffold="SpCas9 mod_7")
        ordered = out.loc[out["Oligo Name"] == "ngRNA_scaffold_top", "Sequence"].iloc[0]
        assert ordered == clone.SCAFFOLD_SEQUENCES["SpCas9 mod_7"][10:]

    def test_mod_7_spacer_oligo_carries_its_own_first_10_nt(self):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["GTCGATCGAT"]})
        out = clone.ngRNAs(df=df, id="ID", order=False, vector_scaffold=False, scaffold="SpCas9 mod_7")
        assert out["ns_Spacer_sequence_top"].iloc[0].endswith("GTTTCAGAGC")


# ---------------------------------------------------------------------------
# sgRNAs()
# ---------------------------------------------------------------------------

class TestSgRNAs:
    def test_order_false_returns_top_bot_columns(self):
        df = pd.DataFrame({"ID": [1, 2], "Spacer_sequence": ["ATCGATCGAT", "GATCGATCGA"]})
        out = clone.sgRNAs(df=df, id="ID", order=False)
        assert "oSpacer_sequence_top" in out.columns
        assert "oSpacer_sequence_bot" in out.columns
        assert len(out) == 2

    def test_order_true_produces_ordering_format(self):
        df = pd.DataFrame({"ID": [1, 2], "Spacer_sequence": ["ATCGATCGAT", "GATCGATCGA"]})
        out = clone.sgRNAs(df=df, id="ID", order=True)
        # order format concatenates top and bot ord_form outputs -> 2 rows per input row
        assert list(out.columns) == ["Oligo Name", "Sequence", "Scale (µmol)", "bp"]
        assert len(out) == 4
        assert set(out["Oligo Name"]) == {"o1_top", "o2_top", "o1_bot", "o2_bot"}


# ---------------------------------------------------------------------------
# get_overhangs()
# ---------------------------------------------------------------------------

PEGRNA_PARTS = lambda promoter=None, scaffold=None, motif=None: [
    (promoter, clone.PEGRNA_PROMOTER_OVERHANGS),
    (scaffold, clone.PEGRNA_SCAFFOLD_OVERHANGS),
    (motif, clone.EPEGRNA_MOTIF_OVERHANGS),
]


class TestGetOverhangs:
    def test_fills_unspecified_overhangs_from_chosen_parts(self):
        out = clone.get_overhangs(parts=PEGRNA_PARTS("7SK", "SpCas9 Flip and Extend", "tevo2"),
                                  spacer_t5=None, spacer_t3=None, spacer_b5=None, spacer_b3=None,
                                  extension_t5=None, extension_t3=None, extension_b5=None, extension_b3=None)
        assert out == {"spacer_t5": "TCTC", "spacer_t3": "GTTTAAGAGC", "spacer_b5": "", "spacer_b3": "",
                       "extension_t5": "", "extension_t3": "", "extension_b5": "GACG", "extension_b3": "GCACCGACTC"}

    def test_hU6_and_tevopreQ1(self):
        out = clone.get_overhangs(parts=PEGRNA_PARTS(promoter="hU6", motif="tevopreQ1"),
                                  spacer_t5=None, extension_b5=None)
        assert out == {"spacer_t5": "CACC", "extension_b5": "CGCG"}

    def test_sgRNA_parts_use_their_own_scaffold_overhangs(self):
        out = clone.get_overhangs(parts=[("hU6", clone.SGRNA_PROMOTER_OVERHANGS),
                                         ("SpCas9 Flip and Extend", clone.SGRNA_SCAFFOLD_OVERHANGS)],
                                  t5=None, t3=None, b5=None, b3=None)
        assert out == {"t5": "CACC", "t3": "", "b5": "AAAC", "b3": ""}

    def test_explicit_overhang_overrides_the_chosen_part(self):
        out = clone.get_overhangs(parts=PEGRNA_PARTS(promoter="7SK"), spacer_t5="AAAA", spacer_b3=None)
        assert out == {"spacer_t5": "AAAA", "spacer_b3": ""}

    def test_empty_string_overrides_rather_than_falls_back(self):
        out = clone.get_overhangs(parts=PEGRNA_PARTS(promoter="7SK"), spacer_t5="")
        assert out == {"spacer_t5": ""}

    def test_part_name_is_case_insensitive(self):
        assert clone.get_overhangs(parts=PEGRNA_PARTS(promoter="7sk"), spacer_t5=None) == {"spacer_t5": "TCTC"}
        assert clone.get_overhangs(parts=PEGRNA_PARTS(motif="TEVO2"), extension_b5=None) == {"extension_b5": "GACG"}

    def test_part_left_as_none_contributes_nothing(self):
        assert clone.get_overhangs(parts=PEGRNA_PARTS(), spacer_t5=None, spacer_t3=None) == {"spacer_t5": "", "spacer_t3": ""}

    def test_unrecognized_name_raises(self):
        with pytest.raises(ValueError, match="not an option"):
            clone.get_overhangs(parts=PEGRNA_PARTS(promoter="U6"), spacer_t5=None)

    def test_prints_the_defaults_it_supplied(self, capsys):
        clone.get_overhangs(parts=PEGRNA_PARTS("7SK", "SpCas9 Flip and Extend", "tevo2"),
                            spacer_t5=None, spacer_t3=None, spacer_b5=None, spacer_b3=None,
                            extension_t5=None, extension_t3=None, extension_b5=None, extension_b3=None)
        out = capsys.readouterr().out
        assert "7SK overhangs: spacer_t5='TCTC', spacer_b3=''" in out
        assert "SpCas9 Flip and Extend overhangs:" in out
        assert "tevo2 overhangs:" in out

    def test_says_nothing_when_every_overhang_was_specified(self, capsys):
        clone.get_overhangs(parts=PEGRNA_PARTS(promoter="7SK"), spacer_t5="AAAA", spacer_b3="TTTT")
        assert capsys.readouterr().out == ""

    def test_reports_only_the_overhangs_it_actually_supplied(self, capsys):
        clone.get_overhangs(parts=PEGRNA_PARTS(promoter="7SK"), spacer_t5="AAAA", spacer_b3=None)
        out = capsys.readouterr().out
        assert "spacer_b3=''" in out
        assert "spacer_t5" not in out


# ---------------------------------------------------------------------------
# sgRNAs()
# ---------------------------------------------------------------------------

class TestSgRNAsOverhangs:
    def test_default_promoter_and_scaffold_overhangs(self):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["GTCGATCGAT"]})
        out = clone.sgRNAs(df=df, id="ID", order=False)
        assert out["oSpacer_sequence_top"].iloc[0] == "CACC" + "GTCGATCGAT"
        assert out["oSpacer_sequence_bot"].iloc[0] == "AAAC" + str(Seq("GTCGATCGAT").reverse_complement())

    def test_7SK_promoter(self):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["GTCGATCGAT"]})
        out = clone.sgRNAs(df=df, id="ID", order=False, promoter="7SK")
        assert out["oSpacer_sequence_top"].iloc[0].startswith("TCTC")

    def test_overhangs_override_named_parts(self):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["GTCGATCGAT"]})
        out = clone.sgRNAs(df=df, id="ID", order=False, t5="AAAA", b5="TTTT")
        assert out["oSpacer_sequence_top"].iloc[0].startswith("AAAA")
        assert out["oSpacer_sequence_bot"].iloc[0].startswith("TTTT")


# ---------------------------------------------------------------------------
# epegRNAs()
# ---------------------------------------------------------------------------

class TestEpegRNAs:
    def test_make_extension_concatenates_rtt_pbs_linker(self):
        df = pd.DataFrame({
            "ID": [1],
            "Spacer_sequence": ["ATCGATCGAT"],
            "RTT_sequence": ["AAAA"],
            "PBS_sequence": ["CCCC"],
            "Linker_sequence": ["GGGG"],
        })
        out = clone.epegRNAs(df=df, id="ID", make_extension=True, order=False, order_scaffold=False)
        assert out["Extension_sequence"].iloc[0] == "AAAACCCCGGGG"

    def test_make_extension_warns_on_null_values(self, capsys):
        df = pd.DataFrame({
            "ID": [1],
            "Spacer_sequence": ["ATCGATCGAT"],
            "RTT_sequence": [None],
            "PBS_sequence": ["CCCC"],
            "Linker_sequence": ["GGGG"],
        })
        out = clone.epegRNAs(df=df, id="ID", make_extension=True, order=False, order_scaffold=False)
        captured = capsys.readouterr()
        assert "Null values found" in captured.out
        assert out["Extension_sequence"].iloc[0] == "CCCCGGGG"  # None filled as empty string

    def test_default_promoter_and_motif_overhangs(self):
        df = pd.DataFrame({
            "ID": [1],
            "Spacer_sequence": ["GTCGATCGAT"],
            "RTT_sequence": ["AAAA"],
            "PBS_sequence": ["CCCC"],
            "Linker_sequence": ["GGGG"],
        })
        out = clone.epegRNAs(df=df, id="ID", order=False)
        assert out["es_Spacer_sequence_top"].iloc[0].startswith("TCTC")  # 7SK
        assert out["es_Spacer_sequence_top"].iloc[0].endswith("GTTTAAGAGC")  # SpCas9 Flip and Extend
        assert out["ee_Extension_sequence_bot"].iloc[0].startswith("GACG")  # tevo2
        assert out["ee_Extension_sequence_bot"].iloc[0].endswith("GCACCGACTC")  # SpCas9 Flip and Extend

    def test_named_parts_select_their_overhangs(self):
        df = pd.DataFrame({
            "ID": [1],
            "Spacer_sequence": ["GTCGATCGAT"],
            "RTT_sequence": ["AAAA"],
            "PBS_sequence": ["CCCC"],
            "Linker_sequence": ["GGGG"],
        })
        out = clone.epegRNAs(df=df, id="ID", order=False, promoter="hU6", epegRNA_motif="tevopreQ1")
        assert out["es_Spacer_sequence_top"].iloc[0].startswith("CACC")
        assert out["ee_Extension_sequence_bot"].iloc[0].startswith("CGCG")

    def test_overhangs_override_named_parts(self):
        df = pd.DataFrame({
            "ID": [1],
            "Spacer_sequence": ["GTCGATCGAT"],
            "RTT_sequence": ["AAAA"],
            "PBS_sequence": ["CCCC"],
            "Linker_sequence": ["GGGG"],
        })
        out = clone.epegRNAs(df=df, id="ID", order=False, promoter="7SK", epegRNA_motif="tevo2",
                             spacer_t5="AAAA", extension_b5="TTTT")
        assert out["es_Spacer_sequence_top"].iloc[0].startswith("AAAA")
        assert out["ee_Extension_sequence_bot"].iloc[0].startswith("TTTT")

    def test_order_true_without_scaffold(self):
        df = pd.DataFrame({
            "ID": [1],
            "Spacer_sequence": ["ATCGATCGAT"],
            "RTT_sequence": ["AAAA"],
            "PBS_sequence": ["CCCC"],
            "Linker_sequence": ["GGGG"],
        })
        out = clone.epegRNAs(df=df, id="ID", make_extension=True, order=True, order_scaffold=False)
        # spacer top/bot + extension top/bot = 4 rows
        assert len(out) == 4
        assert set(out["Oligo Name"]) == {"es_1_top", "es_1_bot", "ee_1_top", "ee_1_bot"}

    def test_order_true_with_scaffold_adds_scaffold_rows(self):
        # Fixed bug: an unconditional (not gated on `order`) leftover block used to call
        # ord_form() looking for a column literally named after the 100-nt scaffold DNA
        # sequence, which never existed, crashing with a KeyError whenever order_scaffold=True.
        # The correct scaffold-handling logic (a hardcoded pd.DataFrame of scaffold rows,
        # added only when order==True) is what actually runs now.
        df = pd.DataFrame({
            "ID": [1],
            "Spacer_sequence": ["ATCGATCGAT"],
            "RTT_sequence": ["AAAA"],
            "PBS_sequence": ["CCCC"],
            "Linker_sequence": ["GGGG"],
        })
        out = clone.epegRNAs(df=df, id="ID", make_extension=True, order=True, order_scaffold=True)
        assert len(out) == 6
        assert "pegRNA_scaffold_top" in out["Oligo Name"].tolist()
        assert "pegRNA_scaffold_bot" in out["Oligo Name"].tolist()


# ---------------------------------------------------------------------------
# ngRNAs()
# ---------------------------------------------------------------------------

class TestNgRNAs:
    def test_order_false(self):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["ATCGATCGAT"]})
        out = clone.ngRNAs(df=df, id="ID", order=False)
        assert "ns_Spacer_sequence_top" in out.columns

    def test_default_is_the_vector_scaffold_sgRNA_style_oligos(self):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["GTCGATCGAT"]})
        out = clone.ngRNAs(df=df, id="ID", order=False)
        assert out["ns_Spacer_sequence_top"].iloc[0] == "TCTC" + "GTCGATCGAT"  # 7SK, no scaffold on the top oligo
        assert out["ns_Spacer_sequence_bot"].iloc[0] == "AAAC" + str(Seq("GTCGATCGAT").reverse_complement())

    def test_vector_scaffold_matches_sgRNAs_apart_from_the_promoter(self):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["GTCGATCGAT"]})
        ng = clone.ngRNAs(df=df.copy(), id="ID", order=False, promoter="hU6")
        sg = clone.sgRNAs(df=df.copy(), id="ID", order=False)
        assert ng["ns_Spacer_sequence_top"].iloc[0] == sg["oSpacer_sequence_top"].iloc[0]
        assert ng["ns_Spacer_sequence_bot"].iloc[0] == sg["oSpacer_sequence_bot"].iloc[0]

    def test_both_scaffolds_give_the_same_vector_scaffold_overhangs(self):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["GTCGATCGAT"]})
        mod_7 = clone.ngRNAs(df=df.copy(), id="ID", order=False, scaffold="SpCas9 mod_7")
        flip = clone.ngRNAs(df=df.copy(), id="ID", order=False, scaffold="SpCas9 Flip and Extend")
        assert mod_7["ns_Spacer_sequence_top"].iloc[0] == flip["ns_Spacer_sequence_top"].iloc[0]
        assert mod_7["ns_Spacer_sequence_bot"].iloc[0] == flip["ns_Spacer_sequence_bot"].iloc[0]

    def test_no_vector_scaffold_restores_the_GG_fragment_overhangs(self):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["GTCGATCGAT"]})
        out = clone.ngRNAs(df=df, id="ID", order=False, vector_scaffold=False, scaffold="SpCas9 Flip and Extend")
        assert out["ns_Spacer_sequence_top"].iloc[0] == "TCTC" + "GTCGATCGAT" + "GTTTAAGAGC"

    def test_hU6_promoter_and_overhang_override(self):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["GTCGATCGAT"]})
        assert clone.ngRNAs(df=df.copy(), id="ID", order=False, promoter="hU6")["ns_Spacer_sequence_top"].iloc[0].startswith("CACC")
        assert clone.ngRNAs(df=df.copy(), id="ID", order=False, spacer_t5="AAAA")["ns_Spacer_sequence_top"].iloc[0].startswith("AAAA")

    def test_order_true_with_scaffold(self):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["ATCGATCGAT"]})
        out = clone.ngRNAs(df=df, id="ID", order=True, order_scaffold=True, vector_scaffold=False)
        assert len(out) == 4  # spacer top/bot + scaffold top/bot
        assert "ngRNA_scaffold_top" in out["Oligo Name"].tolist()
        assert "ngRNA_scaffold_bot" in out["Oligo Name"].tolist()

    def test_order_scaffold_is_ignored_under_vector_scaffold(self, capsys):
        df = pd.DataFrame({"ID": [1], "Spacer_sequence": ["ATCGATCGAT"]})
        out = clone.ngRNAs(df=df, id="ID", order=True, order_scaffold=True)
        assert "Did not order the scaffold" in capsys.readouterr().out
        assert len(out) == 2  # spacer top/bot only
        assert "ngRNA_scaffold_top" not in out["Oligo Name"].tolist()


# ---------------------------------------------------------------------------
# Library pooled GG cloning: epegRNA_pool() / dms_pool() — real default resource CSVs
# ---------------------------------------------------------------------------

SCAFFOLD = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"


class TestEpegRNAPool:
    def test_end_to_end_with_default_resources(self):
        df = pd.DataFrame({
            "Subpool": [1, 2],
            "Spacer_sequence": ["GGCCCAGACTGAGCACGTGA", "AGCTGACGTAGCATGCATGA"],
            "Scaffold_sequence": [SCAFFOLD, SCAFFOLD],
            "RTT_sequence": ["ATCGATCG", "GATCGATC"],
            "PBS_sequence": ["CGTGCTCAGTCTG", "CGTGCTCAGTCTG"],
            "Linker_sequence": ["AAAA", "CCCC"],
        })
        out = clone.epegRNA_pool(
            df=df, barcode="Subpool",
            fwd_homology_t5_val="7SK", rev_homology_t3_val="tevopreQ1",
            file=None, return_df=True,
        )
        assert len(out) == 2
        assert out["Extension_sequence"].tolist() == ["ATCGATCGCGTGCTCAGTCTGAAAA", "GATCGATCCGTGCTCAGTCTGCCCC"]
        # Every oligo must contain its own spacer, scaffold, and extension as substrings
        for _, row in out.iterrows():
            assert row["Spacer_sequence"] in row["Oligonucleotide"]
            assert row["Scaffold_sequence"] in row["Oligonucleotide"]
            assert row["Extension_sequence"] in row["Oligonucleotide"]
            assert row["Oligonucleotide_length"] == len(row["Oligonucleotide"])
        # GG-cloning oligo should carry exactly 2 Esp3I recognition sites (fwd + rev), the
        # standard signature of a correctly-formed Type-IIS cloning insert.
        assert (out["Esp3I"] == 2).all()


class TestDmsPool:
    def test_end_to_end_with_default_resources(self):
        from edms.utils import load_resource_csv
        dms_pcr = load_resource_csv(filename="dms_pcr.csv")
        row = dms_pcr[dms_pcr["Forward Homology Name"] == "FOXA1-S165"].iloc[0]
        flank5 = row["Forward Homology Sequence Trim"]
        flank3 = row["Reverse Homology Sequence Trim"]
        middle = "AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTT"
        template = flank5 + middle + flank3

        df = pd.DataFrame({
            "Subpool": [1, 2],
            "Edit_sequence_with_silent_mutations": [template, template],
        })
        out = clone.dms_pool(
            df=df, barcode="Subpool",
            fwd_homology_t5_val="FOXA1-S165", rev_homology_t3_val="FOXA1-P205",
            file=None, return_df=True,
        )
        assert len(out) == 2
        # the trimmed template strips the flanks used to locate the insert boundaries
        assert (out["Edit_sequence_with_silent_mutations_trimmed"] == middle).all()
        for _, r in out.iterrows():
            assert middle in r["Oligonucleotide"]
        assert (out["Esp3I"] == 2).all()


# ---------------------------------------------------------------------------
# UMI helpers
# ---------------------------------------------------------------------------

class TestGenerateSequences:
    def test_length_1(self):
        assert sorted(clone.generate_sequences(length=1)) == ["A", "C", "G", "T"]

    def test_length_2_exhaustive(self):
        seqs = clone.generate_sequences(length=2)
        assert len(seqs) == 16
        assert len(set(seqs)) == 16  # all unique
        for s in seqs:
            assert len(s) == 2
            assert set(s) <= set("ACGT")


class TestFilterGC:
    def test_hand_computed_boundaries(self):
        # AATT: GC=0/4=0.0 (excluded, below lower bound)
        # GGCC: GC=4/4=1.0 (excluded, above upper bound)
        # ATGC: GC=2/4=0.5 (included, strictly within (0.4, 0.6))
        result = clone.filter_GC(["AATT", "GGCC", "ATGC"], GC_fract=(0.4, 0.6))
        assert result == ["ATGC"]

    def test_bounds_are_exclusive(self):
        # exactly at bound 0.5 should be excluded from a (0.5, 1.0) range and from (0.0, 0.5) range
        # since both comparisons use strict inequality
        assert clone.filter_GC(["ATGC"], GC_fract=(0.5, 1.0)) == []
        assert clone.filter_GC(["ATGC"], GC_fract=(0.0, 0.5)) == []


class TestShuffle:
    def test_shuffle_preserves_multiset(self):
        ls = ["A", "B", "C", "D", "E"]
        shuffled = clone.shuffle(ls)
        assert sorted(shuffled) == sorted(ls)

    def test_shuffle_does_not_mutate_input(self):
        ls = ["A", "B", "C", "D", "E"]
        ls_copy = ls.copy()
        clone.shuffle(ls)
        assert ls == ls_copy


class TestEncodeSequences:
    def test_base_mapping(self):
        enc = clone.encode_sequences(["ACGT", "TGCA"])
        expected = np.array([[0, 1, 2, 3], [3, 2, 1, 0]], dtype=np.uint8)
        np.testing.assert_array_equal(enc, expected)


class TestFastFilterByHamming:
    def test_hand_computed_filtering(self):
        # AAAA kept first. AAAT: distance to AAAA = 1, not > 1 -> excluded.
        # TTTT: distance to AAAA = 4 > 1 -> kept. AATT: distance to AAAA=2>1 and to TTTT=2>1 -> kept.
        result = clone.fast_filter_by_hamming(["AAAA", "AAAT", "TTTT", "AATT"], min_distance=1)
        assert result == ["AAAA", "TTTT", "AATT"]

    def test_empty_input(self):
        assert clone.fast_filter_by_hamming([], min_distance=1) == []

    def test_first_sequence_always_kept(self):
        result = clone.fast_filter_by_hamming(["GGGG"], min_distance=3)
        assert result == ["GGGG"]


class TestCountCsvRows:
    def test_counts_rows_excluding_header(self, tmp_path):
        p = tmp_path / "test.csv"
        p.write_text("col1,col2\n1,2\n3,4\n5,6\n")
        assert clone.count_csv_rows(str(p)) == 3

    def test_header_only(self, tmp_path):
        p = tmp_path / "empty.csv"
        p.write_text("col1,col2\n")
        assert clone.count_csv_rows(str(p)) == 0


class TestUmi:
    def test_umi_end_to_end_from_scratch(self, tmp_path):
        clone.umi(length=4, GC_fract=(0.0, 1.0), hamming=1, nrows=100, pt=None, dir=str(tmp_path))
        csvs = sorted(tmp_path.glob("*.csv"))
        assert len(csvs) >= 2  # the raw-generated file + at least one hamming-filtered file
        final = [c for c in csvs if "hamming" in c.name]
        assert len(final) >= 1
        df = pd.read_csv(final[-1])
        assert "UMI_sequence" in df.columns
        assert len(df) > 0
        # all sequences should be length 4 and drawn from ACGT
        for seq in df["UMI_sequence"]:
            assert len(seq) == 4
            assert set(seq) <= set("ACGT")
        # verify actual pairwise hamming distances satisfy the (strict) filter criterion
        seqs = df["UMI_sequence"].tolist()
        for i in range(len(seqs)):
            for j in range(i + 1, len(seqs)):
                dist = sum(a != b for a, b in zip(seqs[i], seqs[j]))
                assert dist > 1


# ---------------------------------------------------------------------------
# pcr_mm()
# ---------------------------------------------------------------------------

class TestPcrMm:
    def test_hand_computed_master_mix(self):
        primers = pd.Series({("FwdA", "RevA"): 2})
        out = clone.pcr_mm(primers=primers, template_uL=1)
        table = out[("FwdA", "RevA")]
        uL = table["uL"].tolist()
        uL_mm = table["uL MM"].tolist()
        # [water, Q5 buffer, dNTPs, fwd primer, rev primer, template, Q5 polymerase, total]
        expected_uL = [15.75, 5.00, 0.50, 1.25, 1.25, 1.00, 0.25, 25.00]
        expected_uL_mm = [round(v * 2 * 1.1, 2) for v in expected_uL]
        assert uL == pytest.approx(expected_uL)
        assert uL_mm == pytest.approx(expected_uL_mm)

    def test_total_row_equals_total_uL(self):
        primers = pd.Series({("P1", "P2"): 1})
        out = clone.pcr_mm(primers=primers, template_uL=2, total_uL=30)
        table = out[("P1", "P2")]
        assert table["uL"].iloc[-1] == pytest.approx(30.0)


# ---------------------------------------------------------------------------
# pcr_sim()
# ---------------------------------------------------------------------------

class TestPcrSim:
    def test_basic_product_without_extensions(self):
        template = "GGGG" + "ATCGATCG" + "TTTT"
        rev_bind = str(Seq("TTTT").reverse_complement())  # matches 3' end of template
        df = pd.DataFrame({
            "template": [template],
            "fwd": ["GGGG"],
            "rev": [rev_bind],
        })
        out = clone.pcr_sim(df=df, template_col="template", fwd_bind_col="fwd", rev_bind_col="rev")
        # product = fwd primer + middle region + revcomp(rev primer) = GGGG + ATCGATCG + TTTT
        assert out["PCR Product"].iloc[0] == "GGGGATCGATCGTTTT"

    def test_primer_in_product_false_excludes_primers(self):
        template = "GGGG" + "ATCGATCG" + "TTTT"
        rev_bind = str(Seq("TTTT").reverse_complement())
        df = pd.DataFrame({"template": [template], "fwd": ["GGGG"], "rev": [rev_bind]})
        out = clone.pcr_sim(df=df, template_col="template", fwd_bind_col="fwd", rev_bind_col="rev",
                             primer_in_product=False)
        assert out["PCR Product"].iloc[0] == "ATCGATCG"

    def test_fwd_and_rev_extensions(self):
        template = "GGGG" + "ATCGATCG" + "TTTT"
        rev_bind = str(Seq("TTTT").reverse_complement())
        df = pd.DataFrame({
            "template": [template], "fwd": ["GGGG"], "rev": [rev_bind],
            "fwd_ext": ["XXXX"], "rev_ext": ["YYYY"],
        })
        out = clone.pcr_sim(df=df, template_col="template", fwd_bind_col="fwd", rev_bind_col="rev",
                             fwd_ext_col="fwd_ext", rev_ext_col="rev_ext")
        fwd_full = "XXXX" + "GGGG"
        rev_full = "YYYY" + rev_bind
        expected = fwd_full + "ATCGATCG" + str(Seq(rev_full).reverse_complement())
        assert out["PCR Product"].iloc[0] == expected


# ---------------------------------------------------------------------------
# off_targets()
# ---------------------------------------------------------------------------

class TestOffTargets:
    def test_best_alignment_partner_and_score(self):
        df = pd.DataFrame({"seq": ["AAAAAAAAAA", "AAAAAAAAAT", "GGGGGGGGGG"]})
        out = clone.off_targets(df=df, col="seq", dir=None, return_df=True)
        # AAAAAAAAAA best matches AAAAAAAAAT (9 matches * 2 - 1 mismatch = 17), far better than GGGGGGGGGG (all mismatch -> -10)
        row0 = out.iloc[0]
        assert row0["Off Target Sequence"] == "AAAAAAAAAT"
        assert row0["Best Alignment Score"] == pytest.approx(17.0)
        row2 = out.iloc[2]
        assert row2["Best Alignment Score"] == pytest.approx(-10.0)

    def test_missing_column_raises(self):
        df = pd.DataFrame({"other": ["AAAA"]})
        with pytest.raises(ValueError):
            clone.off_targets(df=df, col="seq", dir=None, return_df=True)
