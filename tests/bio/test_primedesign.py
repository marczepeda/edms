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
        assert list(df.columns) == EXPECTED_COLUMNS
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
