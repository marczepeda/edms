"""
Tests for edms.bio.signature

Covers:
- SNV / Indel / Signature dataclasses
- parse_signature_literal() grammar (valid literals, malformed input, edge cases)
- left_align_indels() leftmost-normalization behavior (including two previously-fixed bugs)
- concat_gapped_from_aligned() using both a minimal fake alignment object and real
  Bio.Align.PairwiseAligner alignments
- signature_from_alignment() end-to-end signature construction from an alignment
- expand_signature_units() per-nucleotide unit expansion
- compare_signature_with_n_extra_nt_or_less() / is_reference_match_with_n_extra_nt_or_less()
  boundary behavior at n and n+1 extra nucleotides
"""
from collections import Counter

import pytest
from Bio import Align

from edms.bio.signature import (
    SNV,
    Indel,
    Signature,
    SignatureParseError,
    parse_signature_literal,
    concat_gapped_from_aligned,
    left_align_indels,
    signature_from_alignment,
    expand_signature_units,
    compare_signature_with_n_extra_nt_or_less,
    is_reference_match_with_n_extra_nt_or_less,
)


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
class _FakeAlignment:
    """Minimal stand-in for a Bio.Align Alignment: concat_gapped_from_aligned()
    and signature_from_alignment() only ever touch the `.aligned` attribute,
    which they expect in the same shape Bio.Align produces:
    (ref_blocks, query_blocks), each a sequence of (start, end) pairs.
    """

    def __init__(self, ref_blocks, query_blocks):
        self.aligned = (ref_blocks, query_blocks)


def _global_aligner():
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    aligner.mismatch_score = -1
    aligner.match_score = 1
    return aligner


def _apply_indel(ref: str, ind: Indel) -> str:
    """Reconstruct the mutated sequence implied by a single Indel against ref,
    used to check whether a left-aligned Indel is truly equivalent to the
    original one (same resulting sequence)."""
    if ind.dellen > 0:
        return ref[: ind.pos] + ref[ind.pos + ind.dellen :]
    return ref[: ind.pos] + ind.ins + ref[ind.pos :]


# --------------------------------------------------------------------------- #
# Dataclasses
# --------------------------------------------------------------------------- #
class TestDataclasses:
    def test_snv_fields_and_equality(self):
        a = SNV(pos=5, ref="A", alt="T")
        b = SNV(pos=5, ref="A", alt="T")
        assert a == b
        assert a.pos == 5 and a.ref == "A" and a.alt == "T"

    def test_indel_fields_and_equality(self):
        a = Indel(pos=3, ins="", dellen=2)
        b = Indel(pos=3, ins="", dellen=2)
        assert a == b

    def test_snv_frozen_immutable(self):
        s = SNV(pos=1, ref="A", alt="G")
        with pytest.raises(Exception):
            s.pos = 2  # frozen dataclass -> raises dataclasses.FrozenInstanceError

    def test_indel_frozen_immutable(self):
        i = Indel(pos=1, ins="A", dellen=0)
        with pytest.raises(Exception):
            i.pos = 2

    def test_signature_holds_tuples(self):
        sig = Signature(snvs=(SNV(pos=1, ref="A", alt="G"),), indels=(Indel(pos=2, ins="T", dellen=0),))
        assert sig.snvs == (SNV(pos=1, ref="A", alt="G"),)
        assert sig.indels == (Indel(pos=2, ins="T", dellen=0),)

    def test_signature_equality(self):
        sig1 = Signature(snvs=(SNV(pos=1, ref="A", alt="G"),), indels=())
        sig2 = Signature(snvs=(SNV(pos=1, ref="A", alt="G"),), indels=())
        assert sig1 == sig2


# --------------------------------------------------------------------------- #
# parse_signature_literal
# --------------------------------------------------------------------------- #
class TestParseSignatureLiteral:
    def test_valid_signature_with_snv_and_indel(self):
        lit = (
            "Signature(snvs=(SNV(pos=3, ref='T', alt='A'),), "
            "indels=(Indel(pos=4, ins='GG', dellen=0),))"
        )
        result = parse_signature_literal(lit)
        assert isinstance(result, Signature)
        assert result.snvs == (SNV(pos=3, ref="T", alt="A"),)
        assert result.indels == (Indel(pos=4, ins="GG", dellen=0),)

    def test_valid_signature_multiple_snvs_and_indels(self):
        lit = (
            "Signature(snvs=(SNV(pos=1, ref='A', alt='C'), SNV(pos=9, ref='G', alt='T')), "
            "indels=(Indel(pos=5, ins='', dellen=3), Indel(pos=12, ins='AT', dellen=0)))"
        )
        result = parse_signature_literal(lit)
        assert result == Signature(
            snvs=(SNV(pos=1, ref="A", alt="C"), SNV(pos=9, ref="G", alt="T")),
            indels=(Indel(pos=5, ins="", dellen=3), Indel(pos=12, ins="AT", dellen=0)),
        )

    def test_valid_empty_signature(self):
        result = parse_signature_literal("Signature(snvs=(), indels=())")
        assert result == Signature(snvs=(), indels=())

    def test_signature_with_none_constant_allowed(self):
        # ast.Constant allows None per the code (isinstance check includes `n.value is None`)
        # SNV.ref accepts None fine since dataclass doesn't validate types.
        lit = "Signature(snvs=(SNV(pos=1, ref=None, alt='A'),), indels=())"
        result = parse_signature_literal(lit)
        assert result.snvs == (SNV(pos=1, ref=None, alt="A"),)

    def test_not_valid_python_returns_text_unchanged(self):
        text = "not valid python ((("
        result = parse_signature_literal(text)
        assert result == text

    def test_valid_expr_but_not_signature_returns_text_unchanged_tuple(self):
        text = "(1, 2, 3)"
        result = parse_signature_literal(text)
        assert result == text

    def test_valid_expr_but_not_signature_returns_text_unchanged_string(self):
        text = "'hello'"
        result = parse_signature_literal(text)
        assert result == text

    def test_bare_snv_ctor_returns_text_unchanged(self):
        # SNV(...) alone is a valid, allowed constructor call, but the top-level
        # result is not a Signature instance, so the original text is returned.
        text = "SNV(pos=1, ref='A', alt='T')"
        result = parse_signature_literal(text)
        assert result == text

    def test_disallowed_constructor_raises(self):
        with pytest.raises(SignatureParseError, match="Disallowed constructor"):
            parse_signature_literal("Foo(x=1)")

    def test_positional_args_raise(self):
        with pytest.raises(SignatureParseError, match="Positional arguments"):
            parse_signature_literal("SNV(1, 'A', 'T')")

    def test_positional_args_in_signature_raise(self):
        with pytest.raises(SignatureParseError, match="Positional arguments"):
            parse_signature_literal("Signature((), ())")

    def test_disallowed_syntax_list_raises(self):
        with pytest.raises(SignatureParseError, match="Disallowed syntax"):
            parse_signature_literal("Signature(snvs=[], indels=())")

    def test_disallowed_syntax_binop_raises(self):
        # A syntactically valid python expression whose node type (BinOp) is
        # not handled by conv() and not a SyntaxError, so it propagates as
        # SignatureParseError rather than returning the original text.
        with pytest.raises(SignatureParseError, match="Disallowed syntax"):
            parse_signature_literal("1 + 1")

    def test_disallowed_constant_float_raises(self):
        with pytest.raises(SignatureParseError, match="Disallowed constant"):
            parse_signature_literal("SNV(pos=1.5, ref='A', alt='T')")

    def test_missing_required_kwarg_raises_typeerror(self):
        # Not a SignatureParseError: the dataclass constructor itself rejects
        # the call before any Signature-specific validation happens.
        with pytest.raises(TypeError):
            parse_signature_literal("Signature(snvs=())")


# --------------------------------------------------------------------------- #
# left_align_indels
# --------------------------------------------------------------------------- #
class TestLeftAlignIndels:
    # Homopolymer run of A's at indices 0-3 in "AAAATCG"
    REF_HOMOPOLYMER = "AAAATCG"

    def test_deletion_in_homopolymer_normalizes_to_leftmost(self):
        # Deleting any single "A" out of the AAAA run yields the same result;
        # convention says normalize to the leftmost position (0).
        ind = Indel(pos=3, ins="", dellen=1)
        result = left_align_indels([ind], self.REF_HOMOPOLYMER)
        assert result == [Indel(pos=0, ins="", dellen=1)]
        # sanity: confirm it is truly equivalent to the original
        assert _apply_indel(self.REF_HOMOPOLYMER, ind) == _apply_indel(self.REF_HOMOPOLYMER, result[0])

    def test_deletion_middle_of_homopolymer_normalizes_to_leftmost(self):
        ind = Indel(pos=2, ins="", dellen=1)
        result = left_align_indels([ind], self.REF_HOMOPOLYMER)
        assert result == [Indel(pos=0, ins="", dellen=1)]

    def test_deletion_already_at_leftmost_unchanged(self):
        ind = Indel(pos=0, ins="", dellen=1)
        result = left_align_indels([ind], self.REF_HOMOPOLYMER)
        assert result == [Indel(pos=0, ins="", dellen=1)]

    def test_deletion_outside_repeat_context_unchanged(self):
        # Deleting "T" at pos 4: ref[3] == 'A' != 'T', so no shift should occur.
        ind = Indel(pos=4, ins="", dellen=1)
        result = left_align_indels([ind], self.REF_HOMOPOLYMER)
        assert result == [Indel(pos=4, ins="", dellen=1)]

    def test_single_base_insertion_in_homopolymer_normalizes_to_leftmost(self):
        # Inserting a single "A" right after the AAAA run (pos=4, before "T")
        # is equivalent to inserting it before the run (pos=0).
        ind = Indel(pos=4, ins="A", dellen=0)
        result = left_align_indels([ind], self.REF_HOMOPOLYMER)
        assert result == [Indel(pos=0, ins="A", dellen=0)]
        assert _apply_indel(self.REF_HOMOPOLYMER, ind) == _apply_indel(self.REF_HOMOPOLYMER, result[0])

    def test_single_base_insertion_already_leftmost_unchanged(self):
        ind = Indel(pos=0, ins="A", dellen=0)
        result = left_align_indels([ind], self.REF_HOMOPOLYMER)
        assert result == [Indel(pos=0, ins="A", dellen=0)]

    def test_insertion_not_matching_repeat_unchanged(self):
        # Inserting "C" right after the run: ref[3] == 'A' != 'C', no shift.
        ind = Indel(pos=4, ins="C", dellen=0)
        result = left_align_indels([ind], self.REF_HOMOPOLYMER)
        assert result == [Indel(pos=4, ins="C", dellen=0)]

    def test_noop_indel_passed_through_unchanged(self):
        # ins="" and dellen=0 hits neither branch; the function's else-clause
        # appends it unmodified.
        ind = Indel(pos=5, ins="", dellen=0)
        result = left_align_indels([ind], self.REF_HOMOPOLYMER)
        assert result == [ind]

    def test_two_deletions_from_same_run_merge_after_left_align(self):
        # Deleting one A from different offsets within the same homopolymer
        # run both normalize to pos=0 and should be de-duplicated.
        ind1 = Indel(pos=1, ins="", dellen=1)
        ind2 = Indel(pos=3, ins="", dellen=1)
        result = left_align_indels([ind1, ind2], self.REF_HOMOPOLYMER)
        assert result == [Indel(pos=0, ins="", dellen=1)]

    def test_result_sorted_by_pos_dellen_ins(self):
        ind_a = Indel(pos=6, ins="", dellen=1)  # 'G' at pos6, no shift context
        ind_b = Indel(pos=0, ins="A", dellen=0)  # already leftmost
        result = left_align_indels([ind_a, ind_b], self.REF_HOMOPOLYMER)
        assert [ (r.pos) for r in result ] == sorted(r.pos for r in result)

    # ---- Fixed bugs -------------------------------------------------------- #
    def test_multi_nt_deletion_in_period2_repeat_fully_left_aligned(self):
        """
        Fixed bug: signature.py left_align_indels(), deletion branch.

        `del_seq = ref[ind.pos:ind.pos+ind.dellen]` used to be computed once from the
        ORIGINAL position and never recomputed as `start` decreased, so the while-loop
        compared ref[start-1] against a *fixed* trailing base (`del_seq[-1]`) rather than
        the trailing base of the *current* sliding window (`ref[start-1+ind.dellen]`).
        This only happened to be correct for homopolymer runs (period 1), where every base
        in the run is identical anyway. The trailing base is now recomputed each iteration
        from the current `start`.

        For a period-2 tandem repeat ref="GATATC" (the "AT" repeat occupies indices 1-4),
        deleting 2 nt starting at pos=3 (removing ref[3:5]=="AT") is equivalent to deleting
        at pos=2 ("TA") or pos=1 ("AT") -- all three yield "GATC". The true leftmost
        position is 1.
        """
        ref = "GATATC"
        ind = Indel(pos=3, ins="", dellen=2)
        result = left_align_indels([ind], ref)

        assert len(result) == 1
        expected = Indel(pos=1, ins="", dellen=2)
        assert _apply_indel(ref, expected) == _apply_indel(ref, ind)  # sanity: truly equivalent

        assert result[0] == expected

    def test_multi_nt_insertion_left_align_preserves_sequence(self):
        """
        Fixed bug: signature.py left_align_indels(), insertion branch.

        When shifting an insertion left across a repeat, the code used to decrement `start`
        but leave `ind.ins` untouched. For a single-base insertion this was harmless
        (rotating a 1-character string is a no-op), but for a multi-character insertion the
        shifted Indel decoded to a *different* sequence than the original -- i.e. silent
        data corruption, not just an alignment-convention nicety. The insertion is now
        rotated (its last base moved to the front) each time it shifts left by one, which
        keeps the decoded sequence identical.

        ref = "GATATC" (an "AT" tandem repeat at indices 1-4).
        Inserting "AT" at pos=5 (right before the final "C") gives:
            ref[:5] + "AT" + ref[5:] == "GATAT" + "AT" + "C" == "GATATATC"
        The function now left-aligns this to pos=1 with the insertion rotated to "AT" (since
        rotating "AT" by one is still "AT" here), which decodes to:
            ref[:1] + "AT" + ref[1:] == "G" + "AT" + "ATATC" == "GATATATC"
        -- identical to the original sequence.
        """
        ref = "GATATC"
        ind = Indel(pos=5, ins="AT", dellen=0)
        result = left_align_indels([ind], ref)

        assert len(result) == 1
        original_seq = _apply_indel(ref, ind)
        shifted_seq = _apply_indel(ref, result[0])

        assert original_seq == "GATATATC"
        assert shifted_seq == original_seq
        assert result[0] == Indel(pos=1, ins="AT", dellen=0)


# --------------------------------------------------------------------------- #
# concat_gapped_from_aligned
# --------------------------------------------------------------------------- #
class TestConcatGappedFromAligned:
    def test_snv_only_real_global_alignment(self):
        ref = "ACGTACGT"
        query = "ACGAACGT"  # SNV at pos 3 T->A
        aln = _global_aligner().align(ref, query)[0]
        ref_g, query_g = concat_gapped_from_aligned(ref, query, aln)
        assert ref_g == "ACGTACGT"
        assert query_g == "ACGAACGT"

    def test_deletion_real_global_alignment(self):
        ref = "ACGTGGACGT"
        query = "ACGTACGT"
        aln = _global_aligner().align(ref, query)[0]
        ref_g, query_g = concat_gapped_from_aligned(ref, query, aln)
        assert ref_g == "ACGTGGACGT"
        assert query_g == "ACGT--ACGT"

    def test_insertion_real_global_alignment(self):
        ref = "ACGTACGT"
        query = "ACGTGGACGT"
        aln = _global_aligner().align(ref, query)[0]
        ref_g, query_g = concat_gapped_from_aligned(ref, query, aln)
        assert ref_g == "ACGT--ACGT"
        assert query_g == "ACGTGGACGT"

    def test_diagonal_mismatch_defensive_branch(self):
        # Two aligned blocks with an equal-length gap on both ref and query
        # between them exercises the "both advanced" defensive branch, which
        # treats the gap region as an unaligned diagonal (mismatches).
        ref = "ACGTACGT"
        query = "ACGTTTGT"
        aln = _FakeAlignment([(0, 4), (6, 8)], [(0, 4), (6, 8)])
        ref_g, query_g = concat_gapped_from_aligned(ref, query, aln)
        assert ref_g == ref
        assert query_g == query

    def test_local_alignment_excludes_terminal_gaps_when_disabled(self):
        ref = "TTACGTAA"
        query = "ACGT"
        aln = _FakeAlignment([(2, 6)], [(0, 4)])
        ref_g, query_g = concat_gapped_from_aligned(ref, query, aln, include_terminal_gaps=False)
        assert ref_g == "ACGT"
        assert query_g == "ACGT"

    def test_local_alignment_includes_terminal_gaps_when_enabled(self):
        ref = "TTACGTAA"
        query = "ACGT"
        aln = _FakeAlignment([(2, 6)], [(0, 4)])
        ref_g, query_g = concat_gapped_from_aligned(ref, query, aln, include_terminal_gaps=True)
        assert ref_g == "TTACGTAA"
        assert query_g == "--ACGT--"

    def test_mismatched_block_counts_raises_assertion(self):
        aln = _FakeAlignment([(0, 4), (6, 8)], [(0, 4)])
        with pytest.raises(AssertionError):
            concat_gapped_from_aligned("ACGTACGT", "ACGT", aln)


# --------------------------------------------------------------------------- #
# signature_from_alignment
# --------------------------------------------------------------------------- #
class TestSignatureFromAlignment:
    def test_single_snv(self):
        ref = "ACGTACGT"
        query = "ACGAACGT"
        aln = _global_aligner().align(ref, query)[0]
        sig = signature_from_alignment(aln, ref, query)
        assert sig == Signature(snvs=(SNV(pos=3, ref="T", alt="A"),), indels=())

    def test_multi_nt_deletion_compressed_into_single_indel(self):
        ref = "ACGTGGACGT"
        query = "ACGTACGT"
        aln = _global_aligner().align(ref, query)[0]
        sig = signature_from_alignment(aln, ref, query)
        assert sig == Signature(snvs=(), indels=(Indel(pos=4, ins="", dellen=2),))

    def test_multi_nt_insertion_compressed_into_single_indel(self):
        ref = "ACGTACGT"
        query = "ACGTGGACGT"
        aln = _global_aligner().align(ref, query)[0]
        sig = signature_from_alignment(aln, ref, query)
        assert sig == Signature(snvs=(), indels=(Indel(pos=4, ins="GG", dellen=0),))

    def test_perfect_match_gives_empty_signature(self):
        ref = "ACGTACGT"
        query = "ACGTACGT"
        aln = _global_aligner().align(ref, query)[0]
        sig = signature_from_alignment(aln, ref, query)
        assert sig == Signature(snvs=(), indels=())

    def test_left_alignment_is_applied_for_homopolymer_deletion(self):
        # The aligner naturally reports the gap at the start of the run here,
        # but this test targets the explicit left_align_indels() call inside
        # signature_from_alignment() directly: we hand-construct a fake
        # alignment whose gap sits at the RIGHTMOST-equivalent position within
        # a homopolymer run, and confirm signature_from_alignment() still
        # normalizes it to the leftmost position (pos=1), not just wherever
        # the aligner happened to place the gap.
        ref = "GAAAATCG"    # 'A' run at indices 1-4
        query = "GAAATCG"   # one 'A' deleted
        # ref blocks: [0,4) and [5,8); query blocks: [0,4) and [4,7)
        # -> gap deletes ref[4], the LAST 'A' of the run (rightmost-equivalent choice)
        aln = _FakeAlignment([(0, 4), (5, 8)], [(0, 4), (4, 7)])
        sig = signature_from_alignment(aln, ref, query)
        assert sig == Signature(snvs=(), indels=(Indel(pos=1, ins="", dellen=1),))

    def test_snv_and_deletion_combined(self):
        ref = "ACGTGGACGTT"
        query = "ACGTACGTA"  # deletion of "GG" (pos4-6) plus SNV T->A at the end
        aln = _global_aligner().align(ref, query)[0]
        sig = signature_from_alignment(aln, ref, query)
        assert Indel(pos=4, ins="", dellen=2) in sig.indels
        assert any(s.ref == "T" and s.alt == "A" for s in sig.snvs)


# --------------------------------------------------------------------------- #
# expand_signature_units
# --------------------------------------------------------------------------- #
class TestExpandSignatureUnits:
    def test_snv_expands_to_one_unit(self):
        sig = Signature(snvs=(SNV(pos=5, ref="a", alt="t"),), indels=())
        units = expand_signature_units(sig)
        # ref/alt should be upper-cased in the unit
        assert units == Counter({("SNV", 5, "A", "T"): 1})

    def test_deletion_expands_per_nucleotide(self):
        sig = Signature(snvs=(), indels=(Indel(pos=4, ins="", dellen=3),))
        units = expand_signature_units(sig)
        assert units == Counter({("DEL", 4): 1, ("DEL", 5): 1, ("DEL", 6): 1})

    def test_insertion_expands_per_nucleotide_with_index_and_base(self):
        sig = Signature(snvs=(), indels=(Indel(pos=10, ins="gg", dellen=0),))
        units = expand_signature_units(sig)
        assert units == Counter({("INS", 10, 0, "G"): 1, ("INS", 10, 1, "G"): 1})

    def test_noop_indel_contributes_no_units(self):
        sig = Signature(snvs=(), indels=(Indel(pos=1, ins="", dellen=0),))
        units = expand_signature_units(sig)
        assert units == Counter()

    def test_combined_signature(self):
        sig = Signature(
            snvs=(SNV(pos=1, ref="A", alt="C"),),
            indels=(Indel(pos=4, ins="", dellen=2), Indel(pos=10, ins="GG", dellen=0)),
        )
        units = expand_signature_units(sig)
        expected = Counter(
            {
                ("SNV", 1, "A", "C"): 1,
                ("DEL", 4): 1,
                ("DEL", 5): 1,
                ("INS", 10, 0, "G"): 1,
                ("INS", 10, 1, "G"): 1,
            }
        )
        assert units == expected


# --------------------------------------------------------------------------- #
# compare_signature_with_n_extra_nt_or_less / is_reference_match_with_n_extra_nt_or_less
# --------------------------------------------------------------------------- #
class TestCompareSignatureWithNExtraNtOrLess:
    REF_SIG = Signature(snvs=(SNV(pos=5, ref="A", alt="T"),), indels=())

    def _query_with_n_extra_snvs(self, n):
        extras = tuple(SNV(pos=100 + i, ref="C", alt="G") for i in range(n))
        return Signature(snvs=(SNV(pos=5, ref="A", alt="T"),) + extras, indels=())

    def test_exact_match_zero_extra(self):
        query = self._query_with_n_extra_snvs(0)
        result = compare_signature_with_n_extra_nt_or_less(query, self.REF_SIG, n_extra_nt=0)
        assert result["is_match_with_0_extra_nt_or_less"] is True
        assert result["n_extra_nt"] == 0
        assert not result["missing_from_query"]
        assert not result["extra_in_query"]

    def test_boundary_extra_equal_to_n_passes(self):
        # exactly n=2 extra units, n_extra_nt=2 -> pass (uses <=)
        query = self._query_with_n_extra_snvs(2)
        result = compare_signature_with_n_extra_nt_or_less(query, self.REF_SIG, n_extra_nt=2)
        assert result["n_extra_nt"] == 2
        assert result["is_match_with_2_extra_nt_or_less"] is True

    def test_boundary_extra_one_more_than_n_fails(self):
        # same query (2 extra units) but n_extra_nt=1 -> fail
        query = self._query_with_n_extra_snvs(2)
        result = compare_signature_with_n_extra_nt_or_less(query, self.REF_SIG, n_extra_nt=1)
        assert result["n_extra_nt"] == 2
        assert result["is_match_with_1_extra_nt_or_less"] is False

    def test_boundary_extra_equals_n_plus_one_relative_to_threshold(self):
        # n_extra_nt=3 (threshold) with exactly 3 extra units -> pass;
        # with 4 extra units -> fail. Confirms hard cutoff at the boundary.
        query_at = self._query_with_n_extra_snvs(3)
        query_over = self._query_with_n_extra_snvs(4)

        result_at = compare_signature_with_n_extra_nt_or_less(query_at, self.REF_SIG, n_extra_nt=3)
        result_over = compare_signature_with_n_extra_nt_or_less(query_over, self.REF_SIG, n_extra_nt=3)

        assert result_at["is_match_with_3_extra_nt_or_less"] is True
        assert result_over["is_match_with_3_extra_nt_or_less"] is False

    def test_missing_reference_unit_fails_regardless_of_n(self):
        # query lacks the required reference SNV -> should fail even with a huge n_extra_nt
        query_missing = Signature(snvs=(SNV(pos=100, ref="C", alt="G"),), indels=())
        result = compare_signature_with_n_extra_nt_or_less(query_missing, self.REF_SIG, n_extra_nt=1000)
        assert result["is_match_with_1000_extra_nt_or_less"] is False
        assert result["missing_from_query"] == Counter({("SNV", 5, "A", "T"): 1})

    def test_default_n_extra_nt_is_one(self):
        query = self._query_with_n_extra_snvs(1)
        result = compare_signature_with_n_extra_nt_or_less(query, self.REF_SIG)
        assert "is_match_with_1_extra_nt_or_less" in result
        assert result["is_match_with_1_extra_nt_or_less"] is True

    def test_accepts_counter_inputs_directly(self):
        q_units = expand_signature_units(self._query_with_n_extra_snvs(1))
        r_units = expand_signature_units(self.REF_SIG)
        result = compare_signature_with_n_extra_nt_or_less(q_units, r_units, n_extra_nt=1)
        assert result["is_match_with_1_extra_nt_or_less"] is True

    def test_query_wrong_type_raises_valueerror(self):
        with pytest.raises(ValueError, match="Query must be"):
            compare_signature_with_n_extra_nt_or_less("not a signature", self.REF_SIG)

    def test_reference_wrong_type_raises_valueerror(self):
        with pytest.raises(ValueError, match="Reference must be"):
            compare_signature_with_n_extra_nt_or_less(self.REF_SIG, 12345)


class TestIsReferenceMatchWithNExtraNtOrLess:
    REF_SIG = Signature(snvs=(SNV(pos=5, ref="A", alt="T"),), indels=())

    def test_matches_compare_function_true_case(self):
        query = Signature(snvs=(SNV(pos=5, ref="A", alt="T"), SNV(pos=6, ref="C", alt="G")), indels=())
        assert is_reference_match_with_n_extra_nt_or_less(query, self.REF_SIG, n_extra_nt=1) is True

    def test_matches_compare_function_false_case(self):
        query = Signature(
            snvs=(
                SNV(pos=5, ref="A", alt="T"),
                SNV(pos=6, ref="C", alt="G"),
                SNV(pos=7, ref="C", alt="G"),
            ),
            indels=(),
        )
        assert is_reference_match_with_n_extra_nt_or_less(query, self.REF_SIG, n_extra_nt=1) is False

    def test_default_n_extra_nt(self):
        query = Signature(snvs=(SNV(pos=5, ref="A", alt="T"),), indels=())
        assert is_reference_match_with_n_extra_nt_or_less(query, self.REF_SIG) is True
