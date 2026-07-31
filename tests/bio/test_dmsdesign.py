"""Tests for edms.bio.dmsdesign.

dmsdesign.py used to run argparse parsing and logging/file setup at *module
import time*. That was fixed by moving all of that into ``_configure_cli()``,
which is only invoked from ``main()``. These tests exercise the pure/DataFrame
-light helper functions directly (with hand-computed expected values) and run
a handful of ``main()`` integration tests through a temporary CSV input and
``-out`` directory.
"""
import csv
import sys

import pytest

import edms.bio.dmsdesign as dmsdesign


# ---------------------------------------------------------------------------
# codon_nt_difference
# ---------------------------------------------------------------------------

class TestCodonNtDifference:
    def test_identical_codons(self):
        assert dmsdesign.codon_nt_difference('ATG', 'ATG') == 0

    def test_one_base_difference(self):
        assert dmsdesign.codon_nt_difference('ATG', 'ATC') == 1

    def test_all_bases_different(self):
        assert dmsdesign.codon_nt_difference('ATG', 'CCC') == 3

    def test_case_insensitive(self):
        assert dmsdesign.codon_nt_difference('atg', 'ATG') == 0
        assert dmsdesign.codon_nt_difference('gga', 'GGT') == 1


# ---------------------------------------------------------------------------
# strip_saturation_parentheses / validate_saturation_input
# ---------------------------------------------------------------------------

class TestStripSaturationParentheses:
    def test_splits_left_middle_right(self):
        left, middle, right = dmsdesign.strip_saturation_parentheses('AAA(BBB)CCC')
        assert (left, middle, right) == ('AAA', 'BBB', 'CCC')

    def test_no_left_flank(self):
        left, middle, right = dmsdesign.strip_saturation_parentheses('(GAT)CCG')
        assert (left, middle, right) == ('', 'GAT', 'CCG')

    def test_no_right_flank(self):
        left, middle, right = dmsdesign.strip_saturation_parentheses('ATG(GAT)')
        assert (left, middle, right) == ('ATG', 'GAT', '')


class TestValidateSaturationInput:
    def test_valid_input_does_not_raise(self):
        # Should return normally (no SystemExit).
        dmsdesign.validate_saturation_input('ATG(GATCCG)TAG')

    def test_missing_parentheses_exits(self):
        with pytest.raises(SystemExit) as exc_info:
            dmsdesign.validate_saturation_input('ATGGATCCGTAG')
        assert exc_info.value.code == 1

    def test_more_than_one_parenthesized_region_exits(self):
        with pytest.raises(SystemExit) as exc_info:
            dmsdesign.validate_saturation_input('AT(G(GATCCG)TAG)')
        assert exc_info.value.code == 1

    def test_operator_characters_inside_region_are_invalid(self):
        # For saturation input the parenthesized region must be plain ACGT;
        # a substitution operator inside it is treated as an invalid character.
        with pytest.raises(SystemExit) as exc_info:
            dmsdesign.validate_saturation_input('ATG(GA/T)CCG')
        assert exc_info.value.code == 1


# ---------------------------------------------------------------------------
# process_sequence
# ---------------------------------------------------------------------------

class TestProcessSequence:
    def test_substitution(self):
        result = dmsdesign.process_sequence('ATG(G/A)CCG')
        assert result['reference_sequence'] == 'ATGGCCG'
        assert result['edit_sequence'] == 'ATGaCCG'
        assert result['editnumber_sequence'] == 'ATG1CCG'
        # Keyed on edit text + start position (so identical-text edits at different
        # positions don't collide); see test_duplicate_edit_annotation_text_collision_bug.
        assert result['editformat2sequence'] == {'(G/A)__at_3': ['G', 'a', 1]}
        assert result['editnumber2sequence'] == {1: ['G', 'a']}
        assert result['edit_span_length_w_ref'] == 1
        assert result['edit_span_length_w_edit'] == 1
        assert result['edit_start_in_ref'] == 3
        assert result['edit_stop_in_ref_rev'] == 3

    def test_insertion(self):
        result = dmsdesign.process_sequence('ATG(+ATC)CCG')
        assert result['reference_sequence'] == 'ATGCCG'
        assert result['edit_sequence'] == 'ATGatcCCG'
        assert result['editformat2sequence'] == {'(+ATC)__at_3': ['', 'atc', 1]}
        assert result['edit_span_length_w_ref'] == 0
        assert result['edit_span_length_w_edit'] == 3

    def test_deletion(self):
        result = dmsdesign.process_sequence('ATGGCC(-GCC)TAG')
        assert result['reference_sequence'] == 'ATGGCCGCCTAG'
        assert result['edit_sequence'] == 'ATGGCCTAG'
        assert result['editformat2sequence'] == {'(-GCC)__at_6': ['GCC', '', 1]}
        assert result['edit_span_length_w_ref'] == 3
        assert result['edit_span_length_w_edit'] == 0

    def test_combination_edit(self):
        # Example straight from the module's own CLI help text.
        result = dmsdesign.process_sequence('A(-TG)CTGTGAT(G/C)TCGTGATG(+atcg)A')
        assert result['reference_sequence'] == 'ATGCTGTGATGTCGTGATGA'
        assert result['edit_sequence'] == 'ACTGTGATcTCGTGATGatcgA'
        assert result['editnumber2sequence'] == {
            1: ['TG', ''],
            2: ['G', 'c'],
            3: ['', 'atcg'],
        }
        assert result['editnumber_sequence'] == 'A1CTGTGAT2TCGTGATG3A'

    def test_invalid_character_exits(self):
        with pytest.raises(SystemExit) as exc_info:
            dmsdesign.process_sequence('ATGN(G/A)CCG')
        assert exc_info.value.code == 1

    def test_unmatched_parenthesis_exits(self):
        with pytest.raises(SystemExit) as exc_info:
            dmsdesign.process_sequence('ATG(G/ACCG')
        assert exc_info.value.code == 1

    def test_nested_parentheses_exits(self):
        with pytest.raises(SystemExit) as exc_info:
            dmsdesign.process_sequence('ATG((G/A))CCG')
        assert exc_info.value.code == 1

    def test_empty_parentheses_exits(self):
        with pytest.raises(SystemExit) as exc_info:
            dmsdesign.process_sequence('ATG()CCG')
        assert exc_info.value.code == 1

    def test_multiple_operators_in_one_region_exits(self):
        with pytest.raises(SystemExit) as exc_info:
            dmsdesign.process_sequence('ATG(G/A/T)CCG')
        assert exc_info.value.code == 1

    def test_duplicate_edit_annotation_text_collision_bug(self):
        # Fixed bug: editformat2sequence used to be keyed by the literal '(ref/edit)' text of
        # each edit, so two distinct edit regions sharing identical annotation text collided and
        # editnumber_sequence mis-numbered every occurrence with the last edit's counter.
        # editformat2sequence is now keyed on text + position, and the numbered-sequence
        # rendering substitutes each edit occurrence by its own absolute position.
        result = dmsdesign.process_sequence('AA(A/T)CC(A/T)GG')
        # Two distinct parenthesized edit regions were supplied; both the
        # per-edit lookup table and the numbered-sequence rendering should
        # distinguish them.
        assert len(result['editformat2sequence']) == 2
        assert result['editnumber_sequence'] == 'AA1CC2GG'


# ---------------------------------------------------------------------------
# edit_type_from_target
# ---------------------------------------------------------------------------

class TestEditTypeFromTarget:
    def test_substitution(self):
        assert dmsdesign.edit_type_from_target('ATG(G/A)CCG') == 'substitution'

    def test_insertion(self):
        assert dmsdesign.edit_type_from_target('ATG(+ATC)CCG') == 'insertion'

    def test_deletion(self):
        assert dmsdesign.edit_type_from_target('ATG(-ATC)CCG') == 'deletion'

    def test_combination(self):
        result = dmsdesign.edit_type_from_target('A(-TG)CTGTGAT(G/C)TCGTGATG(+atcg)A')
        assert result == 'substitution& insertion& deletion'


# ---------------------------------------------------------------------------
# changed_positions
# ---------------------------------------------------------------------------

class TestChangedPositions:
    def test_substitution_lowercase_position(self):
        # 'a' at index 3 is the lowercase edited base emitted by process_sequence.
        positions = dmsdesign.changed_positions('ATGGCCG', 'ATGaCCG')
        assert positions == [3]

    def test_insertion_lowercase_positions(self):
        positions = dmsdesign.changed_positions('ATGCCG', 'ATGatcCCG')
        assert positions == [3, 4, 5]

    def test_deletion_falls_back_to_difflib(self):
        # A pure deletion leaves no lowercase bases in edit_sequence, so the
        # difflib fallback path is exercised.
        positions = dmsdesign.changed_positions('ATGGCCGCCTAG', 'ATGGCCTAG')
        assert positions == [6]


# ---------------------------------------------------------------------------
# candidate_silent_codons
# ---------------------------------------------------------------------------

class TestCandidateSilentCodons:
    # reference/edit pair: 13 G's, with the middle base substituted (lowercase
    # 'a' at index 6) mimicking the output of process_sequence for a G->A
    # substitution flanked by Gly codons on both sides.
    REF = 'GGGGGGGGGGGGG'
    EDIT = 'GGGGGGaGGGGGG'

    def test_close_mode_orders_by_distance_then_position(self):
        candidates = dmsdesign.candidate_silent_codons(self.EDIT, self.REF, [6], mode='close')
        starts = [c[1] for c in candidates]
        # Codon starting at 3 (distance 2) is closest, then 9 (distance 4),
        # then 0 (distance 5). The codon at 6 overlaps the intended edit and
        # must never be offered as a silent-mutation candidate.
        assert starts == [3, 9, 0]
        assert 6 not in starts

    def test_synonymous_choices_sorted_by_nt_diff_then_usage(self):
        candidates = dmsdesign.candidate_silent_codons(self.EDIT, self.REF, [6], mode='close')
        _, _, _, codon_ref, synonymous = candidates[0]
        assert codon_ref == 'GGG'
        # All GGN synonyms of GGG differ by exactly 1 nt, so ties break on
        # descending codon usage: GGC(0.34) > GGA(0.25) > GGT(0.16).
        assert [c[0] for c in synonymous] == ['GGC', 'GGA', 'GGT']

    def test_upstream_mode_excludes_codons_after_mutation(self):
        candidates = dmsdesign.candidate_silent_codons(self.EDIT, self.REF, [6], mode='upstream')
        starts = [c[1] for c in candidates]
        assert starts == [3, 0]
        assert 9 not in starts

    def test_downstream_mode_excludes_codons_before_mutation(self):
        candidates = dmsdesign.candidate_silent_codons(self.EDIT, self.REF, [6], mode='downstream')
        starts = [c[1] for c in candidates]
        assert starts == [9]

    def test_allowed_range_restricts_candidates(self):
        candidates = dmsdesign.candidate_silent_codons(
            self.EDIT, self.REF, [6], allowed_range=(0, 3), mode='close'
        )
        starts = [c[1] for c in candidates]
        assert starts == [0]

    def test_no_mutation_positions_uses_sequence_midpoint(self):
        # With an empty mutation_positions list, the function should not
        # crash and should center on len(seq)/2 without raising.
        candidates = dmsdesign.candidate_silent_codons(self.EDIT, self.REF, [], mode='close')
        assert len(candidates) == 4  # all 4 codons are eligible now


# ---------------------------------------------------------------------------
# choose_distributed_candidates
# ---------------------------------------------------------------------------

class TestChooseDistributedCandidates:
    EDIT = 'GGGGGGaGGGGGG'
    REF = 'GGGGGGGGGGGGG'

    def _candidates(self):
        return dmsdesign.candidate_silent_codons(self.EDIT, self.REF, [6], mode='close')

    def test_zero_requested_returns_empty(self):
        assert dmsdesign.choose_distributed_candidates(self._candidates(), 0, None, self.EDIT) == []

    def test_empty_candidates_returns_empty(self):
        assert dmsdesign.choose_distributed_candidates([], 3, None, self.EDIT) == []

    def test_two_requested_spread_across_full_sequence(self):
        selected = dmsdesign.choose_distributed_candidates(self._candidates(), 2, None, self.EDIT)
        starts = sorted(c[1] for c in selected)
        # Targets are the two ends (0 and 12); nearest available codon starts
        # are 0 and 9 (6 is excluded as the intended-edit codon).
        assert starts == [0, 9]

    def test_one_requested_targets_midpoint(self):
        selected = dmsdesign.choose_distributed_candidates(self._candidates(), 1, None, self.EDIT)
        assert len(selected) == 1
        assert selected[0][1] == 3


# ---------------------------------------------------------------------------
# add_silent_mutations
# ---------------------------------------------------------------------------

class TestAddSilentMutations:
    EDIT = 'GGGGGGaGGGGGG'
    REF = 'GGGGGGGGGGGGG'

    def test_zero_silent_mutations_is_a_no_op(self):
        edited, annotation, positions, added = dmsdesign.add_silent_mutations(
            self.EDIT, self.REF, 0, mode='close'
        )
        assert edited == self.EDIT
        assert annotation == ''
        assert positions == ''
        assert added == 0

    def test_close_mode_two_mutations(self):
        edited, annotation, positions, added = dmsdesign.add_silent_mutations(
            self.EDIT, self.REF, 2, mode='close'
        )
        assert added == 2
        assert edited == 'GGGggcaGGggcG'
        assert annotation == 'G2GGG>ggc;G4GGG>ggc'
        assert positions == '3-5;9-11'

    def test_upstream_mode_one_mutation(self):
        edited, annotation, positions, added = dmsdesign.add_silent_mutations(
            self.EDIT, self.REF, 1, mode='upstream'
        )
        assert added == 1
        assert edited == 'GGGggcaGGGGGG'
        assert annotation == 'G2GGG>ggc'

    def test_downstream_mode_one_mutation(self):
        edited, annotation, positions, added = dmsdesign.add_silent_mutations(
            self.EDIT, self.REF, 1, mode='downstream'
        )
        assert added == 1
        assert edited == 'GGGGGGaGGggcG'
        assert annotation == 'G4GGG>ggc'

    def test_distribute_mode_two_mutations(self):
        edited, annotation, positions, added = dmsdesign.add_silent_mutations(
            self.EDIT, self.REF, 2, mode='distribute'
        )
        assert added == 2
        assert edited == 'ggcGGGaGGggcG'
        assert positions == '0-2;9-11'

    def test_invalid_mode_raises_value_error(self):
        with pytest.raises(ValueError):
            dmsdesign.add_silent_mutations(self.EDIT, self.REF, 1, mode='barcode')

    def test_stops_early_when_candidates_exhausted(self):
        # Only 3 eligible codons exist (0, 3, 9); requesting 10 should just
        # return as many as are available rather than erroring.
        edited, annotation, positions, added = dmsdesign.add_silent_mutations(
            self.EDIT, self.REF, 10, mode='close'
        )
        assert added == 3


# ---------------------------------------------------------------------------
# barcode_candidate_pool / barcode_capacity / cyclic_barcode_indices /
# apply_silent_barcode
# ---------------------------------------------------------------------------

class TestBarcodePipeline:
    EDIT = 'GGGGGGaGGGGGG'
    REF = 'GGGGGGGGGGGGG'

    def test_pool_picks_closest_n_candidates(self):
        pool = dmsdesign.barcode_candidate_pool(self.EDIT, self.REF, 2)
        assert [c['codon_start'] for c in pool] == [3, 9]
        assert pool[0]['choices'] == ['GGC', 'GGA', 'GGT']
        assert pool[0]['aa'] == 'G'

    def test_pool_empty_when_n_silent_not_positive(self):
        assert dmsdesign.barcode_candidate_pool(self.EDIT, self.REF, 0) == []
        assert dmsdesign.barcode_candidate_pool(self.EDIT, self.REF, -1) == []

    def test_capacity_multiplies_choice_counts(self):
        pool = dmsdesign.barcode_candidate_pool(self.EDIT, self.REF, 2)
        assert dmsdesign.barcode_capacity(pool) == 9  # 3 choices * 3 choices

    def test_capacity_zero_for_empty_pool(self):
        assert dmsdesign.barcode_capacity([]) == 0

    def test_cyclic_indices_mixed_radix_counting(self):
        pool = dmsdesign.barcode_candidate_pool(self.EDIT, self.REF, 2)
        assert dmsdesign.cyclic_barcode_indices(0, pool) == [0, 0]
        assert dmsdesign.cyclic_barcode_indices(1, pool) == [1, 0]
        assert dmsdesign.cyclic_barcode_indices(3, pool) == [0, 1]
        assert dmsdesign.cyclic_barcode_indices(4, pool) == [1, 1]

    def test_apply_silent_barcode_index_zero_uses_top_choice_everywhere(self):
        edited, annotation, positions, added, pattern, barcode_index, reuse = (
            dmsdesign.apply_silent_barcode(self.EDIT, self.REF, 2, barcode_id=0)
        )
        assert edited == 'GGGggcaGGggcG'
        assert added == 2
        assert barcode_index == 0
        assert reuse == 0
        assert pattern == '2:GGG>ggc;4:GGG>ggc'

    def test_apply_silent_barcode_index_one_varies_first_codon(self):
        edited, annotation, positions, added, pattern, barcode_index, reuse = (
            dmsdesign.apply_silent_barcode(self.EDIT, self.REF, 2, barcode_id=1)
        )
        assert edited == 'GGGggaaGGggcG'
        assert barcode_index == 1
        assert reuse == 0

    def test_apply_silent_barcode_wraps_and_tracks_reuse_count(self):
        # capacity is 9, so barcode_id=10 should map to the same variant as
        # barcode_id=1 but with reuse_count incremented.
        edited1, _, _, _, _, idx1, reuse1 = dmsdesign.apply_silent_barcode(
            self.EDIT, self.REF, 2, barcode_id=1
        )
        edited2, _, _, _, _, idx2, reuse2 = dmsdesign.apply_silent_barcode(
            self.EDIT, self.REF, 2, barcode_id=1 + 9
        )
        assert edited1 == edited2
        assert idx1 == idx2
        assert reuse1 == 0
        assert reuse2 == 1

    def test_apply_silent_barcode_zero_capacity_returns_defaults(self):
        result = dmsdesign.apply_silent_barcode(self.EDIT, self.REF, 0, barcode_id=0)
        assert result == (self.EDIT, '', '', 0, '', '', '')


# ---------------------------------------------------------------------------
# sorted_codons_by_difference / sorted_codons_by_usage
# ---------------------------------------------------------------------------

class TestSortedCodons:
    def test_sorted_by_difference_prefers_most_nt_changes(self):
        # GGG (ref) vs GGC (1 diff), GAT (2 diff), TTT (3 diff... not same aa,
        # but the function only cares about nt distance from codon_ref).
        entries = [['GGC', 0.34], ['GAT', 0.46], ['TTT', 0.45]]
        result = dmsdesign.sorted_codons_by_difference('GGG', entries)
        assert [c[0] for c in result] == ['TTT', 'GAT', 'GGC']

    def test_sorted_by_difference_ties_break_on_usage_then_alpha(self):
        entries = [['GGC', 0.34], ['GGA', 0.25], ['GGT', 0.16]]
        result = dmsdesign.sorted_codons_by_difference('GGG', entries)
        assert [c[0] for c in result] == ['GGC', 'GGA', 'GGT']

    def test_sorted_by_usage_descending(self):
        entries = [['GGC', 0.34], ['GGA', 0.25], ['GGT', 0.16], ['GGG', 0.25]]
        result = dmsdesign.sorted_codons_by_usage(entries)
        assert [c[0] for c in result] == ['GGC', 'GGA', 'GGG', 'GGT']


# ---------------------------------------------------------------------------
# saturating_mutagenesis_input_sequences
# ---------------------------------------------------------------------------

class TestSaturatingMutagenesisInputSequences:
    def test_base_mode_generates_3_variants_per_position(self):
        names, seqs, out_names, ranges = dmsdesign.saturating_mutagenesis_input_sequences(
            't2', 'ATG(GA)CCG', 'base'
        )
        assert len(names) == len(seqs) == len(out_names) == len(ranges) == 6  # 2 positions * 3 bases
        assert out_names[0] == 't2_1_GtoA'
        assert seqs[0] == 'ATG(G/A)ACCG'
        assert ranges[0] == (3, 5)
        # internal name == output name for 'base' mode (no __internal_ suffix).
        assert names[0] == out_names[0]

    def test_aa_subs_mode_generates_substitutions_for_every_other_aa(self):
        names, seqs, out_names, ranges = dmsdesign.saturating_mutagenesis_input_sequences(
            't1', 'ATG(GAT)CCG', 'aa_subs'
        )
        # 1 codon * 19 other amino acids, each with 1-4 synonymous codon
        # choices depending on the target amino acid's codon family size.
        assert len(seqs) == len(names) > 19
        assert out_names[0] == 't1_1_DtoG'
        assert seqs[0] == 'ATG(GAT/GGC)CCG'
        assert names[0] == 't1_1_DtoG__internal_1_G_GGC'
        assert ranges[0] == (3, 6)
        # No wild-type (Asp) substitution should be generated.
        assert not any(name.endswith('_DtoD') for name in out_names)

    def test_aa_ins_mode_inserts_after_each_codon(self):
        names, seqs, out_names, ranges = dmsdesign.saturating_mutagenesis_input_sequences(
            't1', 'ATG(GAT)CCG', 'aa_ins'
        )
        assert out_names[0].startswith('t1_1_DtoD')
        # The reference codon stays outside the parens; only the inserted
        # codon is wrapped in (+...).
        assert seqs[0] == 'ATGGAT(+GGC)CCG'
        # Stop codons ('X') must never be inserted.
        assert not any(name.endswith('toX') for name in out_names)

    def test_aa_dels_mode_deletes_whole_codons(self):
        names, seqs, out_names, ranges = dmsdesign.saturating_mutagenesis_input_sequences(
            't3', 'ATG(GATGAT)CCG', 'aa_dels'
        )
        assert len(seqs) == 2  # one deletion per full codon in the region
        assert out_names[0] == 't3_1_DDtoD'
        assert seqs[0] == 'ATG(-GAT)GATCCG'
        assert out_names[1] == 't3_2_DPtoP'
        assert seqs[1] == 'ATGGAT(-GAT)CCG'

    def test_aa_silent_mode_only_produces_synonymous_codons(self):
        names, seqs, out_names, ranges = dmsdesign.saturating_mutagenesis_input_sequences(
            't4', 'ATG(GGG)CCG', 'aa_silent'
        )
        # GGN family has 4 codons, so 3 synonymous alternatives to GGG.
        assert len(seqs) == 3
        assert all(name.endswith('_GtoG') for name in out_names)
        assert seqs[0] == 'ATG(GGG/GGC)CCG'

    def test_unknown_sm_type_yields_no_designs(self):
        names, seqs, out_names, ranges = dmsdesign.saturating_mutagenesis_input_sequences(
            't5', 'ATG(GAT)CCG', 'not_a_real_mode'
        )
        assert names == seqs == out_names == ranges == []

    def test_invalid_saturation_syntax_exits(self):
        with pytest.raises(SystemExit):
            dmsdesign.saturating_mutagenesis_input_sequences('t6', 'ATGGATCCG', 'base')


# ---------------------------------------------------------------------------
# read_targets
# ---------------------------------------------------------------------------

class TestReadTargets:
    def test_reads_csv_and_skips_header_and_uppercases_sequence(self, tmp_path):
        p = tmp_path / 'targets.csv'
        p.write_text('target_name,target_sequence\ntarget_01,atg(g/a)ccg\n')
        targets = dmsdesign.read_targets(str(p))
        assert targets == [('target_01', 'ATG(G/A)CCG')]

    def test_reads_txt_and_skips_header(self, tmp_path):
        p = tmp_path / 'targets.txt'
        p.write_text('header line ignored\ntarget_01\tATG(G/A)CCG\n')
        targets = dmsdesign.read_targets(str(p))
        assert targets == [('target_01', 'ATG(G/A)CCG')]

    def test_csv_row_with_too_few_columns_exits(self, tmp_path):
        p = tmp_path / 'bad.csv'
        p.write_text('target_name,target_sequence\nonlyname\n')
        with pytest.raises(SystemExit) as exc_info:
            dmsdesign.read_targets(str(p))
        assert exc_info.value.code == 1

    def test_unsupported_extension_exits(self, tmp_path):
        p = tmp_path / 'targets.dat'
        p.write_text('header\nname,ATG\n')
        with pytest.raises(SystemExit) as exc_info:
            dmsdesign.read_targets(str(p))
        assert exc_info.value.code == 1

    def test_blank_lines_in_txt_are_skipped(self, tmp_path):
        p = tmp_path / 'targets.txt'
        p.write_text('header\ntarget_01\tATG(G/A)CCG\n\ntarget_02\tATG(G/C)CCG\n')
        targets = dmsdesign.read_targets(str(p))
        assert targets == [
            ('target_01', 'ATG(G/A)CCG'),
            ('target_02', 'ATG(G/C)CCG'),
        ]


# ---------------------------------------------------------------------------
# build_designs (depends on module-level globals set by _configure_cli)
# ---------------------------------------------------------------------------

class TestBuildDesigns:
    def test_non_saturation_target(self, tmp_path, monkeypatch):
        p = tmp_path / 'targets.csv'
        p.write_text('target_name,target_sequence\ntarget_01,ATG(G/A)CCG\n')
        monkeypatch.setattr(dmsdesign, 'file_in', str(p))
        monkeypatch.setattr(dmsdesign, 'saturation_mutagenesis', False)

        designs = dmsdesign.build_designs()
        assert len(designs) == 1
        design = designs[0]
        assert design['target_name'] == 'target_01'
        assert design['internal_target_name'] == 'target_01'
        assert design['saturation_mutagenesis'] == ''
        assert design['silent_allowed_range'] is None
        assert design['reference_sequence'] == 'ATGGCCG'
        assert design['edit_sequence'] == 'ATGaCCG'

    def test_saturation_target_expands_to_multiple_designs(self, tmp_path, monkeypatch):
        p = tmp_path / 'targets.csv'
        p.write_text('target_name,target_sequence\ntarget_sat,ATG(GA)CCG\n')
        monkeypatch.setattr(dmsdesign, 'file_in', str(p))
        monkeypatch.setattr(dmsdesign, 'saturation_mutagenesis', 'base')

        designs = dmsdesign.build_designs()
        assert len(designs) == 6  # 2 positions * 3 alternate bases
        first = designs[0]
        assert first['target_name'] == 'target_sat_1_GtoA'
        assert first['internal_target_name'] == 'target_sat_1_GtoA'
        assert first['saturation_mutagenesis'] == 'base'
        assert first['silent_allowed_range'] == (3, 5)


# ---------------------------------------------------------------------------
# main() integration tests
# ---------------------------------------------------------------------------

def _read_output_csv(out_dir):
    import glob
    files = glob.glob(str(out_dir / '*_DMSDesign.csv'))
    assert len(files) == 1, 'expected exactly one DMSDesign.csv output file'
    with open(files[0], newline='') as handle:
        return list(csv.DictReader(handle))


EXPECTED_COLUMNS = {
    'Target_name', 'Internal_target_name', 'Target_sequence', 'Design_number', 'Edit_type',
    'Saturation_mutagenesis', 'Reference_sequence', 'Edit_sequence', 'Edit_sequence_with_silent_mutations',
    'Edit_start_in_ref', 'Edit_span_length_ref', 'Edit_span_length_edit',
    'Silent_mutations_requested', 'Silent_mutation_mode', 'Silent_mutations_added',
    'Silent_mutation_annotation', 'Silent_mutation_positions',
    'Silent_mutation_allowed_region', 'Silent_barcode_id', 'Silent_barcode_pattern', 'Silent_barcode_reuse_count',
}


class TestMainIntegration:
    def test_substitution_input(self, tmp_path, monkeypatch):
        input_csv = tmp_path / 'input.csv'
        input_csv.write_text(
            'target_name,target_sequence\n'
            'target_01,ATGTGCTGTGATGGTAT(G/A)CCGGCGTAGTAATCGTAGC\n'
        )
        out_dir = tmp_path / 'out_sub'
        monkeypatch.setattr(sys, 'argv', ['dmsdesign.py', '-f', str(input_csv), '-out', str(out_dir)])
        dmsdesign.main()

        rows = _read_output_csv(out_dir)
        assert len(rows) == 1
        row = rows[0]
        assert EXPECTED_COLUMNS.issubset(row.keys())
        assert row['Target_name'] == 'target_01'
        assert row['Edit_type'] == 'substitution'
        assert row['Reference_sequence'] == 'ATGTGCTGTGATGGTATGCCGGCGTAGTAATCGTAGC'
        assert row['Edit_sequence'] == 'ATGTGCTGTGATGGTATaCCGGCGTAGTAATCGTAGC'
        # No silent mutations were requested (default 0).
        assert row['Silent_mutations_added'] == '0'

    def test_insertion_input_with_silent_barcode_mutations(self, tmp_path, monkeypatch):
        input_csv = tmp_path / 'input.csv'
        input_csv.write_text(
            'target_name,target_sequence\n'
            'target_02,ATGTGCTGTGATGGTATG(+ATCTCGATGA)CCGGCGTAGTAATCGTAGC\n'
        )
        out_dir = tmp_path / 'out_ins'
        monkeypatch.setattr(
            sys, 'argv',
            ['dmsdesign.py', '-f', str(input_csv), '-out', str(out_dir), '-silent_mut', '2'],
        )
        dmsdesign.main()

        rows = _read_output_csv(out_dir)
        assert len(rows) == 1
        row = rows[0]
        assert row['Edit_type'] == 'insertion'
        assert row['Silent_mutation_mode'] == 'barcode'  # default mode
        assert row['Silent_mutations_added'] == '2'
        assert row['Silent_barcode_id'] == '0'
        # Silent-mutated sequence must differ from the plain edit sequence
        # but still be the same length.
        assert row['Edit_sequence_with_silent_mutations'] != row['Edit_sequence']
        assert len(row['Edit_sequence_with_silent_mutations']) == len(row['Edit_sequence'])

    def test_deletion_input_with_distribute_mode(self, tmp_path, monkeypatch):
        input_csv = tmp_path / 'input.csv'
        input_csv.write_text(
            'target_name,target_sequence\n'
            'target_03,ATGTGCTGTGATGG(-TATGCCG)GCGTAGTAATCGTAGC\n'
        )
        out_dir = tmp_path / 'out_del'
        monkeypatch.setattr(
            sys, 'argv',
            ['dmsdesign.py', '-f', str(input_csv), '-out', str(out_dir),
             '-silent_mut', '1', '-silent_mut_mode', 'distribute'],
        )
        dmsdesign.main()

        rows = _read_output_csv(out_dir)
        assert len(rows) == 1
        row = rows[0]
        assert row['Edit_type'] == 'deletion'
        assert row['Silent_mutation_mode'] == 'distribute'
        assert row['Silent_mutations_added'] == '1'
        # Non-barcode modes leave the barcode-specific columns blank.
        assert row['Silent_barcode_id'] == ''
        assert row['Silent_barcode_pattern'] == ''

    def test_saturation_base_mode_expands_designs(self, tmp_path, monkeypatch):
        input_csv = tmp_path / 'input.csv'
        input_csv.write_text(
            'target_name,target_sequence\n'
            'target_sat,ATGTGCTGT(GAT)GGTATGCCGGCGTAGTAATCGTAGC\n'
        )
        out_dir = tmp_path / 'out_sat'
        monkeypatch.setattr(
            sys, 'argv',
            ['dmsdesign.py', '-f', str(input_csv), '-out', str(out_dir), '-sat_mut', 'base'],
        )
        dmsdesign.main()

        rows = _read_output_csv(out_dir)
        assert len(rows) == 9  # 3-base region * 3 alternate bases each
        assert all(row['Saturation_mutagenesis'] == 'base' for row in rows)
        names = {row['Target_name'] for row in rows}
        assert 'target_sat_1_GtoA' in names
        assert rows[0]['Silent_mutation_allowed_region'] == '9-11'

    def test_empty_input_file_exits(self, tmp_path, monkeypatch):
        input_csv = tmp_path / 'input.csv'
        input_csv.write_text('target_name,target_sequence\n')  # header only, no data rows
        out_dir = tmp_path / 'out_empty'
        monkeypatch.setattr(sys, 'argv', ['dmsdesign.py', '-f', str(input_csv), '-out', str(out_dir)])

        with pytest.raises(SystemExit) as exc_info:
            dmsdesign.main()
        assert exc_info.value.code == 1
