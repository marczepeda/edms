"""
Tests for edms.bio.primedesign — pegRNA/ngRNA design (the PrimeDesign tool ported into edms).

primedesign.py used to run argparse.parse_args() and logging/output-dir setup at MODULE IMPORT
TIME. This was fixed by moving all CLI-parsing into `_configure_cli()`, called at the top of
`main()`. The module now imports cleanly with no side effects; `main()` re-parses `sys.argv` via
`_configure_cli()` on every call, so tests drive it by monkeypatching `sys.argv` and calling
`edms.bio.primedesign.main()` directly (verified-working pattern from the task brief).

The smaller functions defined before `main()` (iupac2bases, gc_content, reverse_complement,
codon_nt_difference, silent_codon_in_allowed_range, silent_mutation_edit_sequence_in_allowed_range,
process_sequence, saturating_mutagenesis_input_sequences) are normal importable top-level
functions with no import-time side effects, and are tested directly with hand-computed
expectations.
"""
import sys

import pandas as pd
import pytest

from edms.bio.primedesign import (
    iupac2bases,
    gc_content,
    reverse_complement,
    codon_nt_difference,
    silent_codon_in_allowed_range,
    silent_mutation_edit_sequence_in_allowed_range,
    process_sequence,
    saturating_mutagenesis_input_sequences,
    ranked_silent_mutations,
    pam_spacer_annotation,
    select_silent_mutations,
    extend_rtt_for_silent_mutations,
    with_no_silent_mutation_option,
    target_aa_number,
    codon_annotation,
    codon_columns,
    SILENT_CODON_COLUMNS,
    NUCLEOTIDE_CHANGE_COLUMNS,
    nucleotide_change_columns,
    codon_dict,
)


# ---------------------------------------------------------------------------
# iupac2bases()
# ---------------------------------------------------------------------------

class TestIupac2Bases:
    def test_plain_bases_map_to_themselves(self):
        assert iupac2bases('A') == 'A'
        assert iupac2bases('T') == 'T'
        assert iupac2bases('C') == 'C'
        assert iupac2bases('G') == 'G'

    def test_degenerate_codes(self):
        assert iupac2bases('N') == '[ACTG]'
        assert iupac2bases('R') == '[AG]'
        assert iupac2bases('Y') == '[CT]'
        assert iupac2bases('W') == '[AT]'

    def test_structural_characters_pass_through(self):
        assert iupac2bases('(') == '('
        assert iupac2bases(')') == ')'
        assert iupac2bases('+') == '+'
        assert iupac2bases('-') == '-'
        assert iupac2bases('/') == '/'

    def test_invalid_symbol_exits(self):
        with pytest.raises(SystemExit) as excinfo:
            iupac2bases('Q')
        assert excinfo.value.code == 1


# ---------------------------------------------------------------------------
# gc_content()
# ---------------------------------------------------------------------------

class TestGcContent:
    def test_hand_computed(self):
        assert gc_content('ATCG') == "0.50"
        assert gc_content('AATT') == "0.00"
        assert gc_content('GGCC') == "1.00"

    def test_case_insensitive(self):
        assert gc_content('atcg') == gc_content('ATCG')


# ---------------------------------------------------------------------------
# reverse_complement()
# ---------------------------------------------------------------------------

class TestReverseComplement:
    def test_simple_dna(self):
        assert reverse_complement('ATCG') == 'CGAT'

    def test_lowercase(self):
        assert reverse_complement('atcg') == 'cgat'

    def test_bracket_and_structural_symbols_flip(self):
        # brackets swap direction (open<->close) and the whole string reverses;
        # '(' <-> ')' also swap under this scheme.
        assert reverse_complement('ATG[CGG]') == '[CCG]CAT'

    def test_edit_annotation_symbols_preserved(self):
        assert reverse_complement('(A/T)') == '(A/T)'  # palindromic under paren-swap + revcomp
        # '(' <-> ')' swap AND the string reverses, so '(+AT)' -> complement-in-place ')+TA(' -> reversed '(AT+)'
        assert reverse_complement('(+AT)') == '(AT+)'


# ---------------------------------------------------------------------------
# codon_nt_difference()
# ---------------------------------------------------------------------------

class TestCodonNtDifference:
    def test_identical_codons(self):
        assert codon_nt_difference('ATG', 'ATG') == 0

    def test_one_nt_difference(self):
        assert codon_nt_difference('ATG', 'ATC') == 1

    def test_all_different(self):
        assert codon_nt_difference('ATG', 'CCC') == 3

    def test_case_insensitive(self):
        assert codon_nt_difference('atg', 'ATC') == 1


# ---------------------------------------------------------------------------
# silent_codon_in_allowed_range() / silent_mutation_edit_sequence_in_allowed_range()
# ---------------------------------------------------------------------------

class TestSilentCodonInAllowedRange:
    def test_no_allowed_range_always_true(self):
        assert silent_codon_in_allowed_range({}, 100, 200) is True

    def test_within_range(self):
        entry = {'silent_mutation_allowed_range': (5, 11)}
        assert silent_codon_in_allowed_range(entry, 5, 8) is True

    def test_starts_before_range(self):
        entry = {'silent_mutation_allowed_range': (5, 11)}
        assert silent_codon_in_allowed_range(entry, 4, 8) is False

    def test_ends_after_range(self):
        entry = {'silent_mutation_allowed_range': (5, 11)}
        assert silent_codon_in_allowed_range(entry, 8, 12) is False


class TestSilentMutationEditSequenceInAllowedRange:
    def test_no_allowed_range_always_true(self):
        assert silent_mutation_edit_sequence_in_allowed_range({}, 'AAA', 'CCC') is True

    def test_empty_silent_sequence_always_true(self):
        entry = {'silent_mutation_allowed_range': (3, 6)}
        assert silent_mutation_edit_sequence_in_allowed_range(entry, 'AAATTTAAA', '') is True

    def test_mismatched_lengths_false(self):
        entry = {'silent_mutation_allowed_range': (3, 6)}
        assert silent_mutation_edit_sequence_in_allowed_range(entry, 'AAATTTAAA', 'AAATTT') is False

    def test_change_within_allowed_codon_true(self):
        entry = {'silent_mutation_allowed_range': (3, 6)}
        # differs only within codon [3:6]
        assert silent_mutation_edit_sequence_in_allowed_range(entry, 'AAATTTAAA', 'AAAGGGAAA') is True

    def test_change_outside_allowed_codon_false(self):
        entry = {'silent_mutation_allowed_range': (3, 6)}
        # differs within codon [6:9], outside the allowed [3,6) range
        assert silent_mutation_edit_sequence_in_allowed_range(entry, 'AAATTTAAA', 'AAATTTGAA') is False


# ---------------------------------------------------------------------------
# process_sequence()
# ---------------------------------------------------------------------------

class TestProcessSequenceSubstitution:
    def test_basic_substitution(self):
        (editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence,
         editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit,
         edit_start_in_ref, edit_stop_in_ref_rev) = process_sequence('AAA(G/A)CCC')

        # Keyed on edit text + start position (so identical-text edits at different
        # positions don't collide); see test_duplicate_edit_annotation_text_collision_bug.
        assert editformat2sequence == {'(G/A)__at_3': ['G', 'a', 1]}
        assert editnumber2sequence == {1: ['G', 'a']}
        assert reference_sequence == 'AAAGCCC'
        assert edit_sequence == 'AAAaCCC'
        assert editnumber_sequence == 'AAA1CCC'
        assert edit_span_length_w_ref == 1
        assert edit_span_length_w_edit == 1
        assert edit_start_in_ref == 3   # index of '(' in the raw input string
        assert edit_stop_in_ref_rev == 3  # index of ')' in the reversed raw input string


class TestProcessSequenceInsertion:
    def test_basic_insertion(self):
        (editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence,
         editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit,
         *_rest) = process_sequence('AAAG(+ATCTCGATGA)CCC')

        assert editformat2sequence == {'(+ATCTCGATGA)__at_4': ['', 'atctcgatga', 1]}
        assert reference_sequence == 'AAAGCCC'                     # insertion contributes nothing to reference
        assert edit_sequence == 'AAAGatctcgatgaCCC'
        assert edit_span_length_w_ref == 0
        assert edit_span_length_w_edit == 10


class TestProcessSequenceDeletion:
    def test_basic_deletion(self):
        (editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence,
         editnumber_sequence, edit_span_length_w_ref, edit_span_length_w_edit,
         *_rest) = process_sequence('AAA(-TATGCCG)GCC')

        assert editformat2sequence == {'(-TATGCCG)__at_3': ['TATGCCG', '', 1]}
        assert reference_sequence == 'AAATATGCCGGCC'
        assert edit_sequence == 'AAAGCC'                            # deletion contributes nothing to edit
        assert edit_span_length_w_ref == 7
        assert edit_span_length_w_edit == 0


class TestProcessSequenceCombination:
    def test_combination_edit_matches_docstring_example(self):
        # From the module's own docstring:
        #   Reference: ATGCTGTGAT G TCGTGATG    A
        #   Edit:      A--CTGTGAT C TCGTGATGatcgA
        #   Format:    A(-TG)CTGTGAT(G/C)TCGTGATG(+atcg)A
        (editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence,
         editnumber_sequence, *_rest) = process_sequence('A(-TG)CTGTGAT(G/C)TCGTGATG(+atcg)A')

        assert reference_sequence == 'ATGCTGTGATGTCGTGATGA'
        assert edit_sequence == 'ACTGTGATcTCGTGATGatcgA'
        assert editnumber2sequence == {1: ['TG', ''], 2: ['G', 'c'], 3: ['', 'atcg']}
        assert editnumber_sequence == 'A(1)CTGTGAT(2)TCGTGATG(3)A'.replace('(1)', '1').replace('(2)', '2').replace('(3)', '3')


class TestProcessSequenceValidation:
    def test_invalid_character_exits(self):
        with pytest.raises(SystemExit) as excinfo:
            process_sequence('ATGXYZ(G/A)CCGG')
        assert excinfo.value.code == 1

    def test_unbalanced_parentheses_exits(self):
        with pytest.raises(SystemExit):
            process_sequence('ATG(G/A')

    def test_nested_parentheses_exits(self):
        with pytest.raises(SystemExit):
            process_sequence('ATG((G/A))CCGG')

    def test_empty_parentheses_exits(self):
        with pytest.raises(SystemExit):
            process_sequence('ATG()CCGG')

    def test_multiple_annotations_in_one_set_exits(self):
        with pytest.raises(SystemExit):
            process_sequence('ATG(G/A/C)CCGG')

    def test_no_parentheses_exits(self):
        # zero parenthesis sets: count('(') == count(')') == 0, but the code requires > 0
        with pytest.raises(SystemExit):
            process_sequence('ATGATGATG')

    def test_duplicate_edit_annotation_text_collision_bug(self):
        # Fixed bug: editformat2sequence used to be keyed by the literal edit annotation TEXT
        # (e.g. '(A/T)'), not by position. Two distinct edit regions sharing identical annotation
        # text (e.g. two separate (A/T) substitutions at different positions) collided in the
        # dict, and the naive str.replace(edit_text, ...) calls building editnumber_sequence
        # stamped every occurrence of that text with the same (wrong, last-seen) edit number.
        # editformat2sequence is now keyed on text + position, and each edit occurrence is
        # substituted back in by its own absolute position -- so editnumber_sequence correctly
        # distinguishes the two occurrences ('AA1CC2GG', not 'AA2CC2GG'), which also fixes the
        # downstream PAM-relocation `.replace()` calls in main() (lines ~701-703, 794-800, 835-841)
        # that rely on editnumber_sequence/editnumber2sequence to relocate PAM-adjacent edits.
        (editformat2sequence, editnumber2sequence, reference_sequence, edit_sequence,
         editnumber_sequence, *_rest) = process_sequence('AA(A/T)CC(A/T)GG')

        assert len(editformat2sequence) == 2
        assert editnumber_sequence == 'AA1CC2GG'
        assert editnumber2sequence == {1: ['A', 't'], 2: ['A', 't']}


# ---------------------------------------------------------------------------
# saturating_mutagenesis_input_sequences()
# ---------------------------------------------------------------------------

class TestSaturatingMutagenesisInputSequences:
    def test_base_mode_generates_3_substitutions_per_position(self):
        names, seqs, out_names, ranges = saturating_mutagenesis_input_sequences('t1', 'AAA(TGC)CCC', 'base')
        # 3 nt in the target region * 3 alternate bases each = 9 designs
        assert len(names) == 9
        assert len(seqs) == 9
        assert len(out_names) == 9
        # every allowed range should be (3, 6) -- the span of "TGC" within "AAA(TGC)CCC"
        assert all(r == (3, 6) for r in ranges)
        # first design substitutes position 1 (0-indexed 0 of "TGC" -> 'T') to 'A'
        assert names[0] == 't1_1_TtoA'
        assert seqs[0] == 'AAA(T/A)GCCCC'

    def test_base_mode_names_are_unique(self):
        names, *_ = saturating_mutagenesis_input_sequences('t1', 'AAA(TGC)CCC', 'base')
        assert len(names) == len(set(names))

    def test_aa_subs_mode_produces_all_other_amino_acids(self):
        names, seqs, out_names, ranges = saturating_mutagenesis_input_sequences('t1', 'AAA(TGC)CCC', 'aa_subs')
        # TGC = Cys; substitutions should target every other amino acid (not Cys itself)
        assert len(names) > 0
        assert all('_1_CtoC' not in n for n in out_names)  # never substitutes Cys for Cys
        assert all(s.count('(') == 1 and s.count(')') == 1 for s in seqs)  # exactly one edit region each
        # internal names must be unique even though several rows share the same output_target_name
        assert len(names) == len(set(names))
        assert all(r == (3, 6) for r in ranges)

    def test_aa_dels_mode_produces_one_design_per_codon(self):
        # sequence_to_edit = "TGCGGA" -> 2 codons -> 1 deletion design per codon position
        names, seqs, out_names, ranges = saturating_mutagenesis_input_sequences('t1', 'AAA(TGCGGA)CCC', 'aa_dels')
        assert len(names) == 2
        for s in seqs:
            assert '(-' in s

    def test_invalid_format_exits(self):
        with pytest.raises(SystemExit):
            saturating_mutagenesis_input_sequences('t1', 'AAA(TGC)(CCC)', 'base')  # 2 parenthesis sets


# ---------------------------------------------------------------------------
# main() end-to-end integration tests
# ---------------------------------------------------------------------------
# Fixed, pre-validated sequences (found by brute-force random search against the real PrimeDesign
# algorithm with the default pe_format 'NNNNNNNNNNNNNNNNN/NNN[NGG]') that are long enough and
# happen to contain a usable PAM near the edit, so main() actually finds >=1 pegRNA design.

SUBSTITUTION_SEQ = (
    "AACTTTAAGAAATTATGTGCATGCCTTCAAGACCCAGAGACCTAATCATAGCGCTCCTCA"
    "(T/G)"
    "TTTGGCTCATACGCATCTGGGTCTTCGGCTTGAAATTGAGGGCAACCACGTGACTACTTC"
)
INSERTION_SEQ = (
    "AACTTTAAGAAATTATGTGCATGCCTTCAAGACCCAGAGACCTAATCATAGCGCTCCTCA"
    "(+TACGA)"
    "TTTGGCTCATACGCATCTGGGTCTTCGGCTTGAAATTGAGGGCAACCACGTGACTACTTC"
)
DELETION_SEQ = (
    "ACCTATAAGATTGTCGTTCGCGGATTACATTAAATAACATCGTTGTGGTAAGCGG"
    "(-GAAAG)"
    "CATTTGTGTCGTAGAAAATTGGGTGATGAGCGCGGTTCTAACAAGTAATAATGATAAGCC"
)
# Same as SUBSTITUTION_SEQ, but with the edit moved 3nt left so it sits on an in-frame codon
# boundary suitable for saturation mutagenesis (target region = codon "TCA").
SAT_MUT_SEQ = (
    "AACTTTAAGAAATTATGTGCATGCCTTCAAGACCCAGAGACCTAATCATAGCGCTCC"
    "(TCA)"
    "TTTGGCTCATACGCATCTGGGTCTTCGGCTTGAAATTGAGGGCAACCACGTGACTACTTC"
)

EXPECTED_COLUMNS = [
    'Target_name', 'Target_sequence', 'pegRNA_number', 'gRNA_type', 'Spacer_sequence',
    'Spacer_GC_content', 'PAM_sequence', 'Extension_sequence', 'Strand', 'Annotation',
    'pegRNA-to-edit_distance', 'Nick_index', 'ngRNA-to-pegRNA_distance', 'PBS_length',
    'PBS_GC_content', 'RTT_length', 'RTT_GC_content', 'First_extension_nucleotide',
    'Spacer_sequence_order_TOP', 'Spacer_sequence_order_BOTTOM',
    'pegRNA_extension_sequence_order_TOP', 'pegRNA_extension_sequence_order_BOTTOM',
    'Edit_type', 'Reference_sequence', 'Edit_sequence', 'Silent_mutation_relative_to_edit',
]


def _expected_columns():
    # pulled from the module rather than duplicated, so the two cannot drift
    return EXPECTED_COLUMNS + SILENT_CODON_COLUMNS + NUCLEOTIDE_CHANGE_COLUMNS


def _write_input_csv(tmp_path, rows):
    p = tmp_path / "input.csv"
    lines = ["target_name,target_sequence"]
    for name, seq in rows:
        lines.append(f"{name},{seq}")
    p.write_text("\n".join(lines) + "\n")
    return p


def _run_main(monkeypatch, tmp_path, input_csv, extra_args=None):
    argv = ['primedesign.py', '-f', str(input_csv), '-out', str(tmp_path / 'out')]
    if extra_args:
        argv += extra_args
    monkeypatch.setattr(sys, 'argv', argv)
    import edms.bio.primedesign as primedesign
    primedesign.main()
    out_csvs = list((tmp_path / 'out').glob('*.csv'))
    assert len(out_csvs) == 1, f"expected exactly one output csv, found {out_csvs}"
    return pd.read_csv(out_csvs[0])


class TestMainSubstitution:
    def test_produces_nonempty_output_with_expected_columns(self, tmp_path, monkeypatch):
        input_csv = _write_input_csv(tmp_path, [("target_01", SUBSTITUTION_SEQ)])
        df = _run_main(monkeypatch, tmp_path, input_csv)
        assert list(df.columns) == _expected_columns()
        assert len(df) > 0
        assert set(df['Target_name'].unique()) == {"target_01"}
        # Edit_type is only populated on pegRNA rows (ngRNA rows carry NaN there)
        assert set(df['Edit_type'].dropna().unique()) <= {"substitution"}
        assert set(df['gRNA_type'].unique()) <= {"pegRNA", "ngRNA"}
        assert (df['gRNA_type'] == 'pegRNA').any()


class TestMainInsertion:
    def test_produces_nonempty_output(self, tmp_path, monkeypatch):
        input_csv = _write_input_csv(tmp_path, [("target_02", INSERTION_SEQ)])
        df = _run_main(monkeypatch, tmp_path, input_csv)
        assert len(df) > 0
        assert set(df['Edit_type'].dropna().unique()) <= {"insertion"}


class TestMainDeletion:
    def test_produces_nonempty_output(self, tmp_path, monkeypatch):
        input_csv = _write_input_csv(tmp_path, [("target_03", DELETION_SEQ)])
        df = _run_main(monkeypatch, tmp_path, input_csv)
        assert len(df) > 0
        assert set(df['Edit_type'].dropna().unique()) <= {"deletion"}


class TestMainSaturationMutagenesis:
    def test_sat_mut_aa_subs_produces_many_targets(self, tmp_path, monkeypatch):
        input_csv = _write_input_csv(tmp_path, [("target_04", SAT_MUT_SEQ)])
        df = _run_main(monkeypatch, tmp_path, input_csv, extra_args=['-sat_mut', 'aa_subs'])
        assert len(df) > 0
        # saturation mutagenesis fans a single input out into many named sub-targets
        assert df['Target_name'].nunique() > 1
        assert all(name.startswith('target_04_1_') for name in df['Target_name'].unique())


class TestMainSilentMutation:
    def test_silent_mut_flag_annotates_relative_position(self, tmp_path, monkeypatch):
        input_csv = _write_input_csv(tmp_path, [("target_05", SUBSTITUTION_SEQ)])
        df = _run_main(monkeypatch, tmp_path, input_csv, extra_args=['-silent_mut'])
        assert len(df) > 0
        # with -silent_mut, at least some designs should carry a non-null relative-position tag
        assert df['Silent_mutation_relative_to_edit'].notna().any()
        assert set(df['Silent_mutation_relative_to_edit'].dropna().unique()) <= {"upstream", "downstream", "overlap"}


class TestMainErrorHandling:
    def test_missing_required_file_flag_exits(self, tmp_path, monkeypatch):
        # argparse itself should reject a call with no -f/--file argument
        monkeypatch.setattr(sys, 'argv', ['primedesign.py', '-out', str(tmp_path / 'out')])
        import edms.bio.primedesign as primedesign
        with pytest.raises(SystemExit):
            primedesign.main()

    def test_unsupported_extension_exits(self, tmp_path, monkeypatch):
        bad_input = tmp_path / "input.tsv"
        bad_input.write_text("target_name\ttarget_sequence\ntarget_01\tAAA(T/G)CCC\n")
        monkeypatch.setattr(sys, 'argv', ['primedesign.py', '-f', str(bad_input), '-out', str(tmp_path / 'out')])
        import edms.bio.primedesign as primedesign
        with pytest.raises(SystemExit):
            primedesign.main()


# ---------------------------------------------------------------------------
# Silent mutation selection: ranked_silent_mutations()
# ---------------------------------------------------------------------------

# 20-codon in-frame reference. With the nick at 30 and the SpCas9 pe_format, the
# spacer bases the RTT can reach are [30,33) (codon 10, AGC) and the PAM is
# [33,36) (codon 11, CGG -- matches NGG and has 5 synonyms, so it can be broken
# silently). Codon 15 (CAG, [45,48)) stands in for the intended edit.
RANK_REF = ("ATG" "GCT" "CTG" "ACC" "AAG" "CAT" "GAA" "TTC" "GTG" "CGT"
            "AGC" "CGG" "CCA" "GAT" "AAC" "CAG" "TAC" "ATT" "GGA" "TCT")
NO_RANGE = {'silent_mutation_allowed_range': None}

# SpCas9 'NNNNNNNNNNNNNNNNN/NNN[NGG]' offsets, as main() computes them
CUT_IDX, SPACER_END_IDX, PAM_START_IDX, PAM_END_IDX = 17, 20, 20, 23


def _rank_plus(entry=NO_RANGE, nick=30, rtt=30, edit=(45, 48)):
    return ranked_silent_mutations(entry, RANK_REF, nick, rtt, '+', '[ACGT]GG', edit[0], edit[1])


class TestRankedSilentMutations:
    def test_pam_codon_is_CGG_as_the_fixture_assumes(self):
        # guards the fixture itself: the other tests are meaningless if this drifts
        assert RANK_REF[33:36] == 'CGG'
        assert RANK_REF[30:33] == 'AGC'

    def test_every_candidate_is_synonymous(self):
        for c in _rank_plus():
            assert codon_dict[c['codon_ref']][1] == codon_dict[c['codon_edit']][1]

    def test_every_candidate_lies_wholly_inside_the_rtt(self):
        for c in _rank_plus():
            assert c['codon_start_abs'] >= 30
            assert c['codon_end_abs'] <= 60

    def test_codon_carrying_the_intended_edit_is_excluded(self):
        for c in _rank_plus():
            assert not (c['codon_start_abs'] < 48 and c['codon_end_abs'] > 45)

    def test_pam_disruptors_rank_first(self):
        candidates = _rank_plus()
        disruptors = [c for c in candidates if c['disrupts_pam']]
        assert disruptors, 'fixture should offer PAM-disrupting silent codons'
        assert all(c['disrupts_pam'] for c in candidates[:len(disruptors)])

    def test_pam_disruptors_are_exactly_the_arg_codons_that_break_NGG(self):
        disruptors = {c['codon_edit'] for c in _rank_plus() if c['disrupts_pam']}
        # CGG -> AGG still matches NGG, so it must NOT count as a disruptor
        assert disruptors == {'AGA', 'CGC', 'CGA', 'CGT'}
        assert 'AGG' not in disruptors

    def test_ranking_is_monotonic_on_the_full_key(self):
        key = lambda c: (c['disrupts_pam'], c['disrupts_spacer'], c['nt_difference'], c['usage_edit'])
        candidates = _rank_plus()
        assert all(key(a) >= key(b) for a, b in zip(candidates, candidates[1:]))

    def test_within_the_pam_tier_most_nt_differences_come_first(self):
        disruptors = [c for c in _rank_plus() if c['disrupts_pam']]
        # CGG->AGA differs at 2 positions; the other three differ at 1
        assert disruptors[0]['codon_edit'] == 'AGA'
        assert disruptors[0]['nt_difference'] == 2

    def test_codon_usage_breaks_ties_within_a_tier(self):
        one_nt = [c for c in _rank_plus() if c['disrupts_pam'] and c['nt_difference'] == 1]
        usages = [c['usage_edit'] for c in one_nt]
        assert usages == sorted(usages, reverse=True)
        assert [c['codon_edit'] for c in one_nt] == ['CGC', 'CGA', 'CGT']

    def test_spacer_disruptors_rank_above_plain_candidates(self):
        candidates = _rank_plus()
        rest = [c for c in candidates if not c['disrupts_pam']]
        spacer = [i for i, c in enumerate(rest) if c['disrupts_spacer']]
        assert spacer == list(range(len(spacer)))

    def test_only_the_reachable_spacer_codon_counts_as_spacer_disruption(self):
        for c in _rank_plus():
            if c['disrupts_spacer']:
                assert c['codon_start_abs'] == 30  # the one codon overlapping [30,33)

    def test_allowed_range_restricts_the_search(self):
        entry = {'silent_mutation_allowed_range': (36, 45)}
        candidates = _rank_plus(entry=entry)
        assert candidates
        for c in candidates:
            assert c['codon_start_abs'] >= 36 and c['codon_end_abs'] <= 45
        # the PAM codon is outside that range, so nothing can disrupt the PAM
        assert not any(c['disrupts_pam'] for c in candidates)

    def test_single_codon_amino_acids_yield_nothing(self):
        # TGG (Trp) is the only codon for W, so a PAM made of it cannot be broken silently
        ref = RANK_REF[:33] + 'TGG' + RANK_REF[36:]
        candidates = ranked_silent_mutations(NO_RANGE, ref, 30, 30, '+', '[ACGT]GG', 45, 48)
        assert not any(c['codon_start_abs'] == 33 for c in candidates)

    def test_edit_that_already_broke_the_pam_yields_no_pam_credit(self):
        # allele whose PAM is already destroyed: no candidate should claim to disrupt it
        ref = RANK_REF[:33] + 'CAA' + RANK_REF[36:]
        candidates = ranked_silent_mutations(NO_RANGE, ref, 30, 30, '+', '[ACGT]GG', 45, 48)
        assert candidates
        assert not any(c['disrupts_pam'] for c in candidates)

    def test_minus_strand_searches_the_mirrored_window(self):
        candidates = ranked_silent_mutations(NO_RANGE, RANK_REF, 40, 30, '-', '[ACGT]GG', 0, 0)
        assert candidates
        for c in candidates:
            assert c['codon_start_abs'] >= 10 and c['codon_end_abs'] <= 40
            assert codon_dict[c['codon_ref']][1] == codon_dict[c['codon_edit']][1]

    def test_applying_a_candidate_changes_only_its_own_codon(self):
        for c in _rank_plus():
            new = c['silent_mutation_edit_sequence']
            assert len(new) == len(RANK_REF)
            differing = [i for i, (a, b) in enumerate(zip(RANK_REF.upper(), new.upper())) if a != b]
            assert all(c['codon_start_abs'] <= i < c['codon_end_abs'] for i in differing)


# ---------------------------------------------------------------------------
# pam_spacer_annotation()
# ---------------------------------------------------------------------------

def _annotate(allele, strand='+', nick=30, reference=None):
    return pam_spacer_annotation(reference if reference is not None else RANK_REF, allele,
                                 nick, nick, strand, '[ACGT]GG',
                                 CUT_IDX, SPACER_END_IDX, PAM_START_IDX, PAM_END_IDX)


class TestPamSpacerAnnotation:
    def test_unchanged_allele_is_intact(self):
        assert _annotate(RANK_REF) == 'PAM_&_Spacer_Intact'

    def test_broken_pam_reports_pam_disrupted(self):
        allele = RANK_REF[:33] + 'CAA' + RANK_REF[36:]
        assert _annotate(allele) == 'PAM_disrupted'

    def test_pam_that_still_matches_NGG_is_not_disrupted(self):
        # CGG -> AGG is a changed PAM that still matches, so it is not a disruption
        allele = RANK_REF[:33] + 'AGG' + RANK_REF[36:]
        assert _annotate(allele) == 'PAM_&_Spacer_Intact'

    def test_changed_spacer_with_intact_pam_reports_spacer_disrupted(self):
        allele = RANK_REF[:30] + 'TCT' + RANK_REF[33:]
        assert _annotate(allele) == 'Spacer_disrupted'

    def test_pam_disruption_outranks_spacer_disruption(self):
        allele = RANK_REF[:30] + 'TCT' + 'CAA' + RANK_REF[36:]
        assert _annotate(allele) == 'PAM_disrupted'

    def test_changes_outside_pam_and_spacer_are_ignored(self):
        allele = RANK_REF[:48] + 'TAT' + RANK_REF[51:]
        assert _annotate(allele) == 'PAM_&_Spacer_Intact'

    def test_minus_strand_reads_the_mirrored_offsets(self):
        # nick at 40: PAM is [34,37) read as the reverse complement, spacer tail [37,40)
        ref = RANK_REF[:34] + 'CCT' + RANK_REF[37:]      # revcomp('CCT') == 'AGG', an NGG
        assert _annotate(ref, strand='-', nick=40, reference=ref) == 'PAM_&_Spacer_Intact'
        broken = ref[:34] + 'TTA' + ref[37:]
        assert _annotate(broken, strand='-', nick=40, reference=ref) == 'PAM_disrupted'
        spacer = ref[:37] + ('AAA' if ref[37:40] != 'AAA' else 'TTT') + ref[40:]
        assert _annotate(spacer, strand='-', nick=40, reference=ref) == 'Spacer_disrupted'


# ---------------------------------------------------------------------------
# target_aa_number() / codon_annotation() / codon_columns()
# ---------------------------------------------------------------------------

class TestTargetAaNumber:
    def test_defaults_to_one_with_no_index_and_no_range(self):
        assert target_aa_number({}, 0) == 1
        assert target_aa_number({}, 3) == 2

    def test_non_saturation_counts_from_the_start_of_the_target_sequence(self):
        # cBioPortal convention: aa_index names codon 1 of the whole sequence.
        # FOXA1_S194Afs*127 has its edit at nt 190 -> codon 63 -> 63 + 131 = 194
        entry = {'silent_mutation_allowed_range': None, 'aa_index': 131}
        assert target_aa_number(entry, 0) == 131
        assert target_aa_number(entry, 189) == 194

    def test_saturation_counts_from_the_start_of_the_parenthesized_region(self):
        # FOXA1 convention: index names the first residue of the window (flank5 = 111 nt)
        entry = {'silent_mutation_allowed_range': (111, 417), 'aa_index': 168}
        assert target_aa_number(entry, 111) == 168
        assert target_aa_number(entry, 114) == 169

    def test_codon_before_the_window_numbers_below_the_index(self):
        entry = {'silent_mutation_allowed_range': (111, 417), 'aa_index': 168}
        assert target_aa_number(entry, 108) == 167


class TestCodonAnnotation:
    def test_returns_codon_amino_acid_usage_number_and_index(self):
        assert codon_annotation('ATGGCTCTG', 3, 42) == ['GCT', 'A', codon_dict['GCT'][2], 42, 3]

    def test_blank_for_a_truncated_codon(self):
        assert codon_annotation('ATGGC', 3, 2) == ['', '', '', '', '']

    def test_blank_for_an_ambiguous_codon(self):
        assert codon_annotation('ATGGNT', 3, 2) == ['', '', '', '', '']


class TestCodonColumns:
    def _entry(self):
        return {'silent_mutation_allowed_range': None, 'aa_index': 1,
                'reference_sequence': RANK_REF,
                'edit_sequence': RANK_REF[:45] + 'GAG' + RANK_REF[48:],  # CAG(Q) -> GAG(E)
                'edit_start_in_ref': 45}

    def test_column_count_matches_the_declared_names(self):
        assert len(codon_columns(self._entry(), None)) == len(SILENT_CODON_COLUMNS) == 20

    def test_silent_half_is_blank_when_no_silent_mutation_was_made(self):
        assert codon_columns(self._entry(), None)[10:] == [''] * 10

    def test_target_half_reports_the_edited_codon(self):
        cols = codon_columns(self._entry(), None)
        assert cols[0] == 'CAG' and cols[1] == 'Q'      # reference codon
        assert cols[5] == 'GAG' and cols[6] == 'E'      # edited codon
        assert cols[3] == cols[8] == 16                  # same residue, aa_index 1

    def test_silent_half_reports_a_synonymous_pair(self):
        entry = self._entry()
        selected = ranked_silent_mutations(entry, entry['edit_sequence'], 30, 30, '+',
                                           '[ACGT]GG', 45, 48)[0]
        cols = codon_columns(entry, selected)
        assert cols[10] == selected['codon_ref'] and cols[15] == selected['codon_edit']
        assert cols[11] == cols[16]                      # amino acid unchanged
        assert cols[13] == cols[18]                      # same residue number


# ---------------------------------------------------------------------------
# nucleotide_change_columns()
# ---------------------------------------------------------------------------

NT_REF = 'ATGCATTACGGCTAAGCCATGGTTCAGACCTGA'


def _changes(reference, edit):
    """Return the two index lists as ints."""
    cols = nucleotide_change_columns(reference, edit)
    return ([int(x) for x in cols[0].split(';') if x],
            [int(x) for x in cols[1].split(';') if x])


class TestNucleotideChangeColumns:
    def test_identical_sequences_report_nothing(self):
        assert nucleotide_change_columns(NT_REF, NT_REF) == ['', '']

    def test_returns_one_value_per_declared_column(self):
        assert len(nucleotide_change_columns(NT_REF, NT_REF)) == len(NUCLEOTIDE_CHANGE_COLUMNS) == 2

    def test_single_substitution_is_reported_in_both_coordinate_systems(self):
        edit = NT_REF[:10] + ('T' if NT_REF[10] != 'T' else 'A') + NT_REF[11:]
        assert _changes(NT_REF, edit) == ([10], [10])

    def test_separate_substitutions_are_all_reported(self):
        # stands in for an intended edit plus a silent mutation elsewhere
        edit = list(NT_REF)
        for i in (4, 20, 21):
            edit[i] = 'G' if NT_REF[i] != 'G' else 'C'
        edit_positions, reference_positions = _changes(NT_REF, ''.join(edit))
        assert edit_positions == reference_positions == [4, 20, 21]

    def test_case_is_ignored(self):
        # silent mutations are written in lowercase by PrimeDesign
        edit = NT_REF[:10].lower() + NT_REF[10:]
        assert nucleotide_change_columns(NT_REF, edit) == ['', '']

    def test_deletion_has_reference_positions_but_no_edit_positions(self):
        edit = NT_REF[:12] + NT_REF[17:]          # drops 5 reference bases
        edit_positions, reference_positions = _changes(NT_REF, edit)
        assert edit_positions == []
        assert len(reference_positions) == 5
        # alignment-agnostic: removing the reported reference bases must give the edit
        kept = ''.join(c for i, c in enumerate(NT_REF) if i not in set(reference_positions))
        assert kept == edit

    def test_insertion_has_edit_positions_but_no_reference_positions(self):
        edit = NT_REF[:12] + 'GGCTT' + NT_REF[12:]
        edit_positions, reference_positions = _changes(NT_REF, edit)
        assert reference_positions == []
        assert len(edit_positions) == 5
        kept = ''.join(c for i, c in enumerate(edit) if i not in set(edit_positions))
        assert kept == NT_REF

    def test_deletion_plus_a_distant_substitution(self):
        # the realistic case: an indel edit and a silent mutation far from it
        edit = NT_REF[:4] + ('G' if NT_REF[4] != 'G' else 'C') + NT_REF[5:12] + NT_REF[17:]
        edit_positions, reference_positions = _changes(NT_REF, edit)
        assert 4 in edit_positions and 4 in reference_positions
        assert len(reference_positions) == len(edit_positions) + 5   # 5 deleted bases

    def test_indexes_are_sorted_and_zero_based(self):
        edit = list(NT_REF)
        for i in (0, 7, len(NT_REF) - 1):
            edit[i] = 'G' if NT_REF[i] != 'G' else 'C'
        edit_positions, _ = _changes(NT_REF, ''.join(edit))
        assert edit_positions == sorted(edit_positions)
        assert edit_positions[0] == 0 and edit_positions[-1] == len(NT_REF) - 1

    def test_values_contain_no_commas(self):
        # the output is written as naive comma-joined CSV, so a comma would corrupt it
        edit = list(NT_REF)
        for i in (2, 9, 15):
            edit[i] = 'G' if NT_REF[i] != 'G' else 'C'
        for value in nucleotide_change_columns(NT_REF, ''.join(edit)):
            assert ',' not in value


# ---------------------------------------------------------------------------
# number_of_pegrnas = 0 means "all"
# ---------------------------------------------------------------------------

class TestMainNumberOfPegrnas:
    def _count(self, tmp_path, monkeypatch, extra):
        tmp_path.mkdir(parents=True, exist_ok=True)
        input_csv = _write_input_csv(tmp_path, [("target_01", SAT_MUT_SEQ)])
        df = _run_main(monkeypatch, tmp_path, input_csv, extra_args=extra)
        return df[df['gRNA_type'] == 'pegRNA']

    def test_zero_means_all_not_none(self, tmp_path, monkeypatch):
        # The selection slices the ranked pegRNAs, so a bare [:0] would silently return
        # nothing. This fixture yields at most one pegRNA per design, so 0 and 1 give the
        # same count here; what matters is that 0 is not read as "keep zero".
        capped = self._count(tmp_path / 'capped', monkeypatch,
                             ['-sat_mut', 'aa_subs', '-n_pegrnas', '1'])
        every = self._count(tmp_path / 'every', monkeypatch,
                            ['-sat_mut', 'aa_subs', '-n_pegrnas', '0'])
        assert len(every) > 0
        assert len(every) >= len(capped)

    def test_an_explicit_cap_limits_pegrnas_per_design(self, tmp_path, monkeypatch):
        capped = self._count(tmp_path / 'capped', monkeypatch,
                             ['-sat_mut', 'aa_subs', '-n_pegrnas', '1'])
        per_design = capped.groupby(['Target_name', 'Target_sequence'])['pegRNA_number'].nunique()
        assert per_design.max() == 1

    def test_zero_is_the_default(self, tmp_path, monkeypatch):
        explicit = self._count(tmp_path / 'explicit', monkeypatch,
                               ['-sat_mut', 'aa_subs', '-n_pegrnas', '0'])
        default = self._count(tmp_path / 'default', monkeypatch, ['-sat_mut', 'aa_subs'])
        assert len(default) == len(explicit)


# ---------------------------------------------------------------------------
# non-saturation runs emit every PAM- or spacer-disrupting silent mutation
# ---------------------------------------------------------------------------

class TestMainSilentMutationOptions:
    def _pegrnas(self, tmp_path, monkeypatch, seq, extra=None):
        tmp_path.mkdir(parents=True, exist_ok=True)
        input_csv = _write_input_csv(tmp_path, [("target_01", seq)])
        df = _run_main(monkeypatch, tmp_path, input_csv, extra_args=(extra or []) + ['-silent_mut'])
        return df[df['gRNA_type'] == 'pegRNA']

    def _plain(self, tmp_path, monkeypatch, seq):
        """Same run without -silent_mut, as the comparison point."""
        tmp_path.mkdir(parents=True, exist_ok=True)
        input_csv = _write_input_csv(tmp_path, [("target_01", seq)])
        df = _run_main(monkeypatch, tmp_path, input_csv)
        return df[df['gRNA_type'] == 'pegRNA']

    def test_emits_more_pegrnas_than_without_silent_mutations(self, tmp_path, monkeypatch):
        # the edit is fixed, so any extra rows can only come from silent mutation options
        plain = self._plain(tmp_path / 'plain', monkeypatch, SUBSTITUTION_SEQ)
        options = self._pegrnas(tmp_path / 'options', monkeypatch, SUBSTITUTION_SEQ)
        assert len(options) > len(plain)

    def test_every_design_also_gets_a_no_silent_mutation_pegrna(self, tmp_path, monkeypatch):
        # pandas reads the literal 'None' written for these rows as NaN
        peg = self._pegrnas(tmp_path / 'cover', monkeypatch, SUBSTITUTION_SEQ)
        assert len(peg) > 0
        design = ['Spacer_sequence', 'RTT_length', 'PBS_length']
        none_rows = peg[peg['Silent_mutation_relative_to_edit'].isna()].groupby(design).size()
        assert len(none_rows) == peg.groupby(design).ngroups   # one for every design
        assert (none_rows == 1).all()                          # and exactly one each

    def test_most_designs_also_offer_a_silent_mutation(self, tmp_path, monkeypatch):
        peg = self._pegrnas(tmp_path / 'silent', monkeypatch, SUBSTITUTION_SEQ)
        assert peg['Silent_mutation_relative_to_edit'].notna().mean() > 0.5

    def test_pegrnas_without_a_silent_mutation_have_the_shortest_rtts(self):
        # the ones that miss out should be the short RTTs, where the arm leaves no room
        short = ranked_silent_mutations(NO_RANGE, RANK_REF, 30, 18, '+', '[ACGT]GG', 45, 48,
                                        homology_downstream=15)
        long_enough = ranked_silent_mutations(NO_RANGE, RANK_REF, 30, 30, '+', '[ACGT]GG', 45, 48,
                                              homology_downstream=15)
        assert len(short) < len(long_enough)

    def test_options_are_distinct_molecules(self, tmp_path, monkeypatch):
        # one row per option, never the same extension twice for a spacer
        peg = self._pegrnas(tmp_path / 'distinct', monkeypatch, SUBSTITUTION_SEQ)
        assert not peg.duplicated(subset=['Spacer_sequence', 'Extension_sequence']).any()

    def test_maximum_limits_options_per_design(self, tmp_path, monkeypatch):
        # at most the maximum, plus the one no-silent-mutation design
        capped = self._pegrnas(tmp_path / 'cap2', monkeypatch, SUBSTITUTION_SEQ,
                               extra=['-n_silent_mut', '2'])
        per_design = capped.groupby(['Spacer_sequence', 'RTT_length', 'PBS_length']).size()
        assert per_design.max() <= 3

    def test_uncapped_emits_more_than_capped(self, tmp_path, monkeypatch):
        capped = self._pegrnas(tmp_path / 'c1', monkeypatch, SUBSTITUTION_SEQ,
                               extra=['-n_silent_mut', '1'])
        every = self._pegrnas(tmp_path / 'c0', monkeypatch, SUBSTITUTION_SEQ)
        assert len(every) > len(capped)

    def test_options_for_one_design_differ_only_in_the_silent_codon(self, tmp_path, monkeypatch):
        peg = self._pegrnas(tmp_path / 'group', monkeypatch, SUBSTITUTION_SEQ)
        grouped = peg.groupby(['Spacer_sequence', 'RTT_length', 'PBS_length'])
        multi = [sub for _, sub in grouped if len(sub) > 1]
        assert multi, 'expected at least one design with several silent options'
        for sub in multi:
            assert sub['Silent_edit_codon_index'].nunique() + sub['Silent_edit_codon'].nunique() > 1


    def test_saturation_mutagenesis_still_emits_one_option(self, tmp_path, monkeypatch):
        # enumeration is deliberately limited to the non-saturation path
        peg = self._pegrnas(tmp_path / 'sat', monkeypatch, SAT_MUT_SEQ, extra=['-sat_mut', 'aa_subs'])
        per_design = peg.groupby(['Target_name', 'Target_sequence', 'Spacer_sequence',
                                  'RTT_length', 'PBS_length']).size()
        assert per_design.max() == 1


# ---------------------------------------------------------------------------
# select_silent_mutations()
# ---------------------------------------------------------------------------

def _cand(name, pam=False, spacer=False):
    return {'codon_edit': name, 'disrupts_pam': pam, 'disrupts_spacer': spacer}


def _cand2(name, pam=False, spacer=False, distance=0):
    return {'codon_edit': name, 'disrupts_pam': pam, 'disrupts_spacer': spacer,
            'distance_to_pam': distance}


class TestSelectSilentMutations:
    def test_default_keeps_every_disruptor(self):
        candidates = [_cand2('A', pam=True), _cand2('B', spacer=True), _cand2('C', pam=True)]
        assert [c['codon_edit'] for c in select_silent_mutations(candidates)] == ['A', 'B', 'C']

    def test_default_adds_the_non_disruptor_closest_to_the_pam(self):
        candidates = [_cand2('pam', pam=True),
                      _cand2('far', distance=30), _cand2('near', distance=3), _cand2('mid', distance=12)]
        assert [c['codon_edit'] for c in select_silent_mutations(candidates)] == ['pam', 'near']

    def test_default_keeps_only_one_non_disruptor(self):
        candidates = [_cand2('a', distance=9), _cand2('b', distance=3), _cand2('c', distance=21)]
        assert len(select_silent_mutations(candidates)) == 1

    def test_ties_on_distance_keep_the_better_ranked(self):
        # the list arrives ranked, so the first of an equal-distance pair wins
        candidates = [_cand2('first', distance=6), _cand2('second', distance=6)]
        assert select_silent_mutations(candidates)[0]['codon_edit'] == 'first'

    def test_default_with_only_disruptors_adds_nothing_extra(self):
        candidates = [_cand2('A', pam=True), _cand2('B', spacer=True)]
        assert len(select_silent_mutations(candidates)) == 2

    def test_a_maximum_takes_the_top_of_the_full_ranked_list(self):
        candidates = [_cand2(str(i), distance=i) for i in range(6)]
        assert [c['codon_edit'] for c in select_silent_mutations(candidates, 3)] == ['0', '1', '2']

    def test_a_maximum_can_exceed_what_the_curated_set_gives(self):
        candidates = [_cand2('pam', pam=True)] + [_cand2(f'n{i}', distance=i + 1) for i in range(5)]
        assert len(select_silent_mutations(candidates)) == 2          # curated: disruptor + nearest
        assert len(select_silent_mutations(candidates, 4)) == 4       # explicit maximum overrides

    def test_no_candidates_means_one_pegrna_with_no_silent_mutation(self):
        assert select_silent_mutations([]) == [None]
        assert select_silent_mutations([], 3) == [None]

    def test_does_not_mutate_the_input(self):
        candidates = [_cand2('A', pam=True), _cand2('B')]
        before = list(candidates)
        select_silent_mutations(candidates)
        assert candidates == before


class TestWithNoSilentMutationOption:
    def test_appends_the_no_silent_design(self):
        options = with_no_silent_mutation_option([_cand2('A', pam=True)])
        assert len(options) == 2 and options[-1] is None

    def test_does_not_double_up_when_there_is_nothing_to_add(self):
        assert with_no_silent_mutation_option([None]) == [None]

    def test_keeps_the_silent_options_first_and_in_order(self):
        options = with_no_silent_mutation_option([_cand2('A'), _cand2('B')])
        assert [c['codon_edit'] for c in options[:-1]] == ['A', 'B']


class TestDistanceToPam:
    def test_candidates_carry_their_distance_to_the_pam(self):
        candidates = ranked_silent_mutations(NO_RANGE, RANK_REF, 30, 30, '+', '[ACGT]GG', 45, 48)
        assert all('distance_to_pam' in c for c in candidates)

    def test_a_codon_overlapping_the_pam_is_distance_zero(self):
        # the PAM sits at [33,36); the codon there must report 0
        candidates = ranked_silent_mutations(NO_RANGE, RANK_REF, 30, 30, '+', '[ACGT]GG', 45, 48)
        at_pam = [c for c in candidates if c['codon_start_abs'] == 33]
        assert at_pam and all(c['distance_to_pam'] == 0 for c in at_pam)

    def test_distance_grows_with_separation_from_the_pam(self):
        candidates = ranked_silent_mutations(NO_RANGE, RANK_REF, 30, 30, '+', '[ACGT]GG', 45, 48)
        by_start = {c['codon_start_abs']: c['distance_to_pam'] for c in candidates}
        far = max(by_start)
        assert by_start[far] > 0 and by_start[far] >= by_start[min(by_start)]



# ---------------------------------------------------------------------------
# homology arm exclusion / extend_rtt_for_silent_mutations()
# ---------------------------------------------------------------------------

# RANK_REF with the codon at [30,33) replaced by TGG (Trp), the one amino acid with
# no synonym, so the first codon in the window can never supply a candidate
NO_SYNONYM_REF = RANK_REF[:30] + 'TGG' + RANK_REF[33:]


class TestHomologyArmExclusion:
    def test_arm_is_excluded_from_the_search(self):
        # RTT [30,60); a 15 nt arm leaves [30,45)
        kept = ranked_silent_mutations(NO_RANGE, RANK_REF, 30, 30, '+', '[ACGT]GG', 45, 48,
                                       homology_downstream=15)
        assert kept
        assert all(c['codon_end_abs'] <= 45 for c in kept)

    def test_zero_arm_searches_the_whole_rtt(self):
        whole = ranked_silent_mutations(NO_RANGE, RANK_REF, 30, 30, '+', '[ACGT]GG', 45, 48)
        kept = ranked_silent_mutations(NO_RANGE, RANK_REF, 30, 30, '+', '[ACGT]GG', 45, 48,
                                       homology_downstream=15)
        assert len(whole) > len(kept)

    def test_pam_disruptors_survive_the_exclusion(self):
        # the PAM sits next to the nick, at the far end from the arm, so it must be unaffected
        kept = ranked_silent_mutations(NO_RANGE, RANK_REF, 30, 30, '+', '[ACGT]GG', 45, 48,
                                       homology_downstream=15)
        assert any(c['disrupts_pam'] for c in kept)

    def test_minus_strand_excludes_the_mirrored_arm(self):
        # RTT [10,40); a 15 nt arm leaves [25,40)
        kept = ranked_silent_mutations(NO_RANGE, RANK_REF, 40, 30, '-', '[ACGT]GG', 0, 0,
                                       homology_downstream=15)
        assert kept
        assert all(c['codon_start_abs'] >= 25 for c in kept)


class TestExtendRttForSilentMutations:
    def _extend(self, reference, rtt_length, maximum, strand='+', hd=15):
        return extend_rtt_for_silent_mutations(NO_RANGE, reference, 30, rtt_length, strand,
                                               '[ACGT]GG', 45, 48, hd, maximum)

    def test_leaves_the_rtt_alone_when_a_candidate_already_fits(self):
        length, candidates = self._extend(RANK_REF, 30, 50)
        assert length == 30 and candidates

    def test_extends_until_a_candidate_appears(self):
        # rtt 18 with a 15 nt arm leaves only [30,33), which is the no-synonym codon
        assert ranked_silent_mutations(NO_RANGE, NO_SYNONYM_REF, 30, 18, '+', '[ACGT]GG', 45, 48, 15) == []
        length, candidates = self._extend(NO_SYNONYM_REF, 18, 50)
        assert length > 18 and candidates
        assert all(c['codon_end_abs'] <= 30 + length - 15 for c in candidates)

    def test_extension_is_the_shortest_that_works(self):
        # from 19 the window ends mid-codon, so the shortest workable extension is 2, not a
        # whole codon: stepping a codon at a time would overshoot to 22 and fail this
        length, _ = self._extend(NO_SYNONYM_REF, 19, 50)
        assert length == 21
        for shorter in range(19, length):
            assert ranked_silent_mutations(NO_RANGE, NO_SYNONYM_REF, 30, shorter, '+',
                                           '[ACGT]GG', 45, 48, 15) == []

    def test_never_reaches_the_maximum_length(self):
        length, _ = self._extend(NO_SYNONYM_REF, 18, 24)
        assert length < 24

    def test_returns_the_original_length_when_no_extension_helps(self):
        # a maximum that leaves no room to grow must not lengthen the pegRNA for nothing
        length, candidates = self._extend(NO_SYNONYM_REF, 18, 19)
        assert length == 18 and candidates == []
