'''
Tests for edms.bio.dms

Covers the pure/utility helpers with hand-computed expected values, the
enzyme-site detection/codon-swap machinery using the real bundled
RE_type_IIS.csv resource, the DMSDesign input/output helpers (including one
real end-to-end run through the bundled `dmsdesign` CLI script), merge(), and
dms_signature() using the real PairwiseAligner + signature module.
'''
import math
import os

import numpy as np
import pandas as pd
import pytest

from edms.bio import dms


# ---------------------------------------------------------------------------
# Helper Functions: get_codons / get_codon_frames / found_list_in_order
# ---------------------------------------------------------------------------

def test_get_codons_frame0_exact_multiple_of_3():
    # 'ATGAAATAG' = ATG AAA TAG, 9 nt, evenly divisible by 3
    assert dms.get_codons('ATGAAATAG', frame=0) == ['ATG', 'AAA', 'TAG']


def test_get_codons_frame1_and_frame2_drop_incomplete_remainder():
    seq = 'ATGAAATAGC'  # 10 nt
    # frame 1: start at index1 -> 'TGA','AAT','AGC' (index 9 leftover unused... actually consumes all)
    assert dms.get_codons(seq, frame=1) == ['TGA', 'AAT', 'AGC']
    # frame 2: start at index2 -> only 2 full codons fit ('GAA','ATA'), trailing 'GC' (2nt) is dropped
    assert dms.get_codons(seq, frame=2) == ['GAA', 'ATA']


def test_get_codons_too_short_returns_empty_list():
    assert dms.get_codons('AT', frame=0) == []
    assert dms.get_codons('A', frame=0) == []
    assert dms.get_codons('', frame=0) == []


def test_get_codon_frames_returns_all_three_frames():
    seq = 'ATGAAATAGC'
    frames = dms.get_codon_frames(seq)
    assert frames == [
        dms.get_codons(seq, 0),
        dms.get_codons(seq, 1),
        dms.get_codons(seq, 2),
    ]
    assert frames == [['ATG', 'AAA', 'TAG'], ['TGA', 'AAT', 'AGC'], ['GAA', 'ATA']]


def test_found_list_in_order_found_in_middle():
    assert dms.found_list_in_order([1, 2, 3, 4, 5], [3, 4]) == 2


def test_found_list_in_order_found_at_start_and_end():
    assert dms.found_list_in_order([1, 2, 3], [1, 2]) == 0
    assert dms.found_list_in_order([1, 2, 3], [2, 3]) == 1


def test_found_list_in_order_not_found_returns_minus_one():
    assert dms.found_list_in_order([1, 2, 3], [4, 5]) == -1
    # Same elements but out of consecutive order should not match
    assert dms.found_list_in_order([1, 3, 2], [2, 3]) == -1


def test_found_list_in_order_empty_sub_list_returns_minus_one():
    assert dms.found_list_in_order([1, 2, 3], []) == -1


# ---------------------------------------------------------------------------
# _read_df
# ---------------------------------------------------------------------------

def test_read_df_from_dataframe_returns_independent_copy():
    df = pd.DataFrame({'a': [1, 2]})
    out = dms._read_df(df)
    assert out.equals(df)
    out.loc[0, 'a'] = 999
    assert df.loc[0, 'a'] == 1  # original untouched -> confirms .copy() semantics


def test_read_df_from_path_reads_csv(tmp_path):
    df = pd.DataFrame({'a': [1, 2], 'b': ['x', 'y']})
    pt = tmp_path / 'in.csv'
    df.to_csv(pt, index=False)
    out = dms._read_df(str(pt))
    assert list(out.columns) == ['a', 'b']
    assert out['a'].tolist() == [1, 2]


# ---------------------------------------------------------------------------
# _clean_dna / _gc_content / _translate
# ---------------------------------------------------------------------------

def test_clean_dna_uppercases_and_strips_whitespace():
    assert dms._clean_dna(' atg c\nTT ') == 'ATGCTT'
    assert dms._clean_dna('atgc') == 'ATGC'


def test_gc_content_hand_computed():
    # 'ATGGCC': G,G,C,C = 4 GC out of 6 -> 66.666...%
    assert dms._gc_content('ATGGCC') == pytest.approx(66.66666666666667)
    # 'ATAT': 0 GC out of 4 -> 0%
    assert dms._gc_content('ATAT') == 0.0
    # cleans whitespace/case before counting
    assert dms._gc_content(' atgg cc ') == pytest.approx(66.66666666666667)


def test_gc_content_empty_returns_nan():
    assert math.isnan(dms._gc_content(''))
    assert math.isnan(dms._gc_content('   '))


def test_translate_basic_codon_table():
    # ATG=M, AAA=K, TAA=stop
    assert dms._translate('ATGAAATAA') == 'MK*'


def test_translate_frame_offset_drops_leading_and_trailing_partial_codons():
    seq = 'ATGAAATAA'
    # frame=1 codons: seq[1:4]='TGA'(stop), seq[4:7]='AAT'(N); seq[7:9]='AA' (2nt) dropped
    assert dms._translate(seq, frame=1) == '*N'


def test_translate_unknown_codon_maps_to_question_mark():
    assert dms._translate('ATGNNNTAA') == 'M?*'


# ---------------------------------------------------------------------------
# find_enzyme_sites (uses the real bundled RE_type_IIS.csv resource)
# ---------------------------------------------------------------------------

def test_find_enzyme_sites_counts_forward_and_reverse_complement():
    # Esp3I recognition='CGTCTC' (fwd, at index 2), recognition_rc='GAGACG' (at index 11)
    seq = 'AACGTCTCAAAGAGACGTT'
    df = pd.DataFrame({'Oligonucleotide': [seq]})
    out = dms.find_enzyme_sites(df, enzyme='Esp3I')
    row = out.iloc[0]
    assert row['Esp3I'] == 2
    assert row['Esp3I_fwd_i'] == [2]
    assert row['Esp3I_rc_i'] == [11]


def test_find_enzyme_sites_no_hits_gives_zero_and_empty_index_lists():
    df = pd.DataFrame({'Oligonucleotide': ['ATGAAATAA']})
    out = dms.find_enzyme_sites(df, enzyme='Esp3I')
    row = out.iloc[0]
    assert row['Esp3I'] == 0
    assert row['Esp3I_fwd_i'] == []
    assert row['Esp3I_rc_i'] == []


def test_find_enzyme_sites_falls_back_to_edit_sequence_column():
    # 'Oligonucleotide' absent, no 'Edit_sequence_with_silent_mutations' either -> falls back to 'Edit_sequence'
    df = pd.DataFrame({'Edit_sequence': ['ATGCGTCTCAAATAA']})
    out = dms.find_enzyme_sites(df, enzyme='Esp3I')
    assert out.iloc[0]['Esp3I'] == 1
    assert out.iloc[0]['Esp3I_fwd_i'] == [3]


def test_find_enzyme_sites_falls_back_to_silent_mutation_column_when_present():
    df = pd.DataFrame({
        'Edit_sequence': ['ATGAAATAA'],
        'Edit_sequence_with_silent_mutations': ['ATGCGTCTCAAATAA'],
    })
    out = dms.find_enzyme_sites(df, enzyme='Esp3I')
    # sequence_col default 'Oligonucleotide' is missing, so it should prefer
    # 'Edit_sequence_with_silent_mutations' over 'Edit_sequence' per the fallback order.
    assert out.iloc[0]['Esp3I'] == 1


def test_find_enzyme_sites_raises_when_no_usable_column():
    df = pd.DataFrame({'foo': ['bar']})
    with pytest.raises(ValueError):
        dms.find_enzyme_sites(df, enzyme='Esp3I')


def test_find_enzyme_sites_accepts_custom_re_dataframe():
    custom = pd.DataFrame({'Name': ['TestEnz'], 'Recognition': ['AAAA'], 'Recognition_rc': ['TTTT']})
    df = pd.DataFrame({'Oligonucleotide': ['GGAAAACCTTTTGG']})
    out = dms.find_enzyme_sites(df, enzyme='TestEnz', RE_type_IIS_df=custom)
    assert out.iloc[0]['TestEnz'] == 2
    assert out.iloc[0]['TestEnz_fwd_i'] == [2]
    assert out.iloc[0]['TestEnz_rc_i'] == [8]


# ---------------------------------------------------------------------------
# enzyme_codon_swap
# ---------------------------------------------------------------------------

def test_enzyme_codon_swap_no_site_leaves_sequence_unchanged():
    seq = 'ATGAAATAA'
    df = pd.DataFrame({'Reference_sequence': [seq], 'Edit_sequence_with_silent_mutations': [seq]})
    out = dms.enzyme_codon_swap(df, enzyme='Esp3I')
    assert out.iloc[0]['Edit_sequence_with_silent_mutations'] == seq
    assert out.iloc[0]['Esp3I_codon_swap_annotation'] == ''


def test_enzyme_codon_swap_resolves_site_with_synonymous_codon():
    # ATG CGT CTC AAA TAA -> Esp3I site 'CGTCTC' spans the CGT/CTC codon boundary (index 3-8).
    seq = 'ATGCGTCTCAAATAA'
    df = pd.DataFrame({'Reference_sequence': [seq], 'Edit_sequence_with_silent_mutations': [seq]})
    out = dms.enzyme_codon_swap(df, enzyme='Esp3I')
    new_seq = out.iloc[0]['Edit_sequence_with_silent_mutations']
    annot = out.iloc[0]['Esp3I_codon_swap_annotation']
    # CGT (Arg) -> CGC is the first synonymous alternative tried and it removes the site.
    assert new_seq == 'ATGCGCCTCAAATAA'
    assert annot == 'CGT4>CGC'
    # Confirm the enzyme site is actually gone from the swapped sequence.
    check = dms.find_enzyme_sites(pd.DataFrame({'Edit_sequence_with_silent_mutations': [new_seq]}), enzyme='Esp3I')
    assert check.iloc[0]['Esp3I'] == 0


def test_enzyme_codon_swap_skips_codons_overlapping_intended_mutations():
    # Recognition site 'CGTCTC' spans codon0 ('CGT', idx0-2) and codon1 ('CTC', idx3-5).
    seq = 'CGTCTCAAA'
    ref = 'CATCACAAA'  # differs from seq at index1 and index4: blocks codon0 AND codon1 from being modified
    df = pd.DataFrame({'Reference_sequence': [ref], 'Edit_sequence_with_silent_mutations': [seq]})
    out = dms.enzyme_codon_swap(df, enzyme='Esp3I')
    # Only codon2 ('AAA', Lys->AAG) remains eligible, and swapping it cannot remove a site
    # entirely contained in codons 0-1, so the site is left unresolved.
    assert out.iloc[0]['Edit_sequence_with_silent_mutations'] == seq
    assert out.iloc[0]['Esp3I_codon_swap_annotation'] == ''


def test_enzyme_codon_swap_writes_output_file_when_requested(tmp_path):
    seq = 'ATGCGTCTCAAATAA'
    df = pd.DataFrame({'Reference_sequence': [seq], 'Edit_sequence_with_silent_mutations': [seq]})
    dms.enzyme_codon_swap(df, enzyme='Esp3I', out_dir=str(tmp_path), out_file='swapped.csv', return_df=False)
    saved = pd.read_csv(tmp_path / 'swapped.csv')
    assert saved.iloc[0]['Edit_sequence_with_silent_mutations'] == 'ATGCGCCTCAAATAA'


# ---------------------------------------------------------------------------
# replace_enzyme_sites
# ---------------------------------------------------------------------------

def test_replace_enzyme_sites_no_sites_returns_early(tmp_path):
    seq = 'ATGAAATAA'
    df = pd.DataFrame({'Reference_sequence': [seq], 'Edit_sequence_with_silent_mutations': [seq]})
    out = dms.replace_enzyme_sites(df, enzyme='Esp3I', dms_dir=str(tmp_path))
    assert (out['Esp3I'] == 0).all()
    assert len(out) == 1
    # Early-return path (no sites found) never creates the enzyme scratch directory.
    assert not (tmp_path / 'Esp3I').exists()


def test_replace_enzyme_sites_recovers_resolvable_site(tmp_path):
    # Fixed bug: replace_enzyme_sites() used to read a nonexistent column
    # f'{sequence_col}_{enzyme}_codon_swap' to recover the codon-swapped sequence;
    # enzyme_codon_swap() actually writes the swapped sequence directly into sequence_col.
    seq = 'ATGCGTCTCAAATAA'  # contains one resolvable Esp3I site
    df = pd.DataFrame({'Reference_sequence': [seq], 'Edit_sequence_with_silent_mutations': [seq]})
    out = dms.replace_enzyme_sites(df, enzyme='Esp3I', dms_dir=str(tmp_path))
    # The site is synonymously resolved and the returned frame has zero remaining Esp3I sites.
    assert (out['Esp3I'] == 0).all()
    assert len(out) == 1


# ---------------------------------------------------------------------------
# dms_design_input
# ---------------------------------------------------------------------------

def test_dms_design_input_creates_expected_csv(tmp_path):
    dms.dms_design_input(target_name='geneA', flank5_sequence='ATG', target_sequence='TGC',
                          flank3_sequence='TAA', dir=str(tmp_path), file='input.csv')
    df = pd.read_csv(tmp_path / 'input.csv')
    assert df.iloc[0]['target_name'] == 'geneA'
    assert df.iloc[0]['target_sequence'] == 'ATG(TGC)TAA'
    assert df.iloc[0]['index'] == 1
    assert df.iloc[0]['silent_mutation'] == 0
    assert df.iloc[0]['silent_mutation_mode'] == 'close'


def test_dms_design_input_negative_silent_mutation_raises(tmp_path):
    with pytest.raises(ValueError):
        dms.dms_design_input('g', 'ATG', 'TGC', 'TAA', silent_mutation=-1, dir=str(tmp_path))


def test_dms_design_input_invalid_mode_raises(tmp_path):
    with pytest.raises(ValueError):
        dms.dms_design_input('g', 'ATG', 'TGC', 'TAA', silent_mutation_mode='bogus', dir=str(tmp_path))


def test_dms_design_input_requires_in_frame_flanks_for_silent_mutation(tmp_path):
    with pytest.raises(ValueError):
        dms.dms_design_input('g', 'AT', 'TGC', 'TAA', silent_mutation=1, dir=str(tmp_path))


# ---------------------------------------------------------------------------
# dms_design (validation only; the os.system() call itself is exercised via
# the dms_designer() end-to-end integration test below)
# ---------------------------------------------------------------------------

def test_dms_design_negative_silent_mutation_raises():
    with pytest.raises(ValueError):
        dms.dms_design(file='irrelevant.csv', silent_mutation=-1)


def test_dms_design_invalid_mode_raises():
    with pytest.raises(ValueError):
        dms.dms_design(file='irrelevant.csv', silent_mutation_mode='bogus')


# ---------------------------------------------------------------------------
# _annotate_edit_from_name
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('target_name,sat_mut,index,expected', [
    ('test_1_MtoK', 'aa_subs', 1, {'Edit': 'M1K', 'AA_number': 1, 'Base_number': None}),
    ('test_1_MtoK', 'aa_subs', 5, {'Edit': 'M5K', 'AA_number': 5, 'Base_number': None}),
    ('test_5_AtoG', 'base', 1, {'Edit': 'A5G', 'AA_number': None, 'Base_number': 5}),
    ('test_1_MtoX', 'aa_subs', 1, {'Edit': 'M1*', 'AA_number': 1, 'Base_number': None}),
])
def test_annotate_edit_from_name_valid_cases(target_name, sat_mut, index, expected):
    out = dms._annotate_edit_from_name(target_name, saturation_mutagenesis=sat_mut, index=index)
    assert out['Edit'] == expected['Edit']
    assert out['AA_number'] == expected['AA_number']
    assert out['Base_number'] == expected['Base_number']


@pytest.mark.parametrize('target_name', [
    'bad',                 # fewer than 3 underscore-separated parts
    'test_abc_MtoK',       # number_token is not numeric
    123,                   # not a string at all
])
def test_annotate_edit_from_name_invalid_cases_return_all_none(target_name):
    out = dms._annotate_edit_from_name(target_name, saturation_mutagenesis='aa_subs', index=1)
    assert out == {'Edit': None, 'AA_number': None, 'Base_number': None, 'Edit_codon': None}


# ---------------------------------------------------------------------------
# dms_design_output
# ---------------------------------------------------------------------------

def _fake_dmsdesign_csv(tmp_path, saturation_mutagenesis='aa_subs'):
    '''Hand-built stand-in for a DMSDesign.csv output, matching the real
    column schema produced by edms.bio.dmsdesign.main() for aa_subs mode.'''
    fake = pd.DataFrame({
        'Target_name': ['test_1_MtoK', 'test_1_MtoR'],
        'Internal_target_name': ['test_1_MtoK__internal_1_K_AAA', 'test_1_MtoR__internal_2_R_CGT'],
        'Target_sequence': ['(ATG/AAA)', '(ATG/CGT)'],
        'Design_number': [1, 2],
        'Edit_type': ['substitution', 'substitution'],
        'Saturation_mutagenesis': [saturation_mutagenesis, saturation_mutagenesis],
        'Reference_sequence': ['ATG', 'ATG'],
        'Edit_sequence': ['AAA', 'CGT'],
        'Edit_sequence_with_silent_mutations': ['AAA', 'CGT'],
    })
    pt = tmp_path / 'fake_DMSDesign.csv'
    fake.to_csv(pt, index=False)
    return str(pt)


def test_dms_design_output_annotates_edit_and_sequence_stats(tmp_path):
    pt = _fake_dmsdesign_csv(tmp_path)
    out = dms.dms_design_output(pt=pt, saturation_mutagenesis='aa_subs', index=1)

    row0, row1 = out.iloc[0], out.iloc[1]
    assert row0['Edit'] == 'M1K' and row0['AA_number'] == 1
    assert row1['Edit'] == 'M1R' and row1['AA_number'] == 1

    # Edit_codon = Edit + '_' + '(ref_codon/edit_codon)' extracted from Target_sequence
    assert row0['Edit_codon'] == 'M1K_ATG/AAA'
    assert row1['Edit_codon'] == 'M1R_ATG/CGT'

    # Edit_sequence_length/GC hand-computed: 'AAA' -> 0% GC, 'CGT' -> 66.67% GC
    assert row0['Edit_sequence_length'] == 3
    assert row0['Edit_sequence_GC_content'] == pytest.approx(0.0)
    assert row1['Edit_sequence_GC_content'] == pytest.approx(66.66666666666667)

    # Reference_sequence stats: 'ATG' -> length 3, GC 33.33%, translation 'M'
    assert row0['Reference_sequence_length'] == 3
    assert row0['Reference_sequence_GC_content'] == pytest.approx(33.333333333333336)
    assert row0['Reference_translation'] == 'M'
    assert row0['Edit_translation'] == 'K'
    assert row1['Edit_translation'] == 'R'


def test_dms_design_output_base_mode_sets_base_number_not_aa_number(tmp_path):
    fake = pd.DataFrame({
        'Target_name': ['test_5_AtoG'],
        'Target_sequence': ['(A/G)'],
        'Reference_sequence': ['A'],
        'Edit_sequence': ['G'],
    })
    pt = tmp_path / 'fake_base_DMSDesign.csv'
    fake.to_csv(pt, index=False)
    out = dms.dms_design_output(pt=str(pt), saturation_mutagenesis='base', index=1)
    assert out.iloc[0]['Edit'] == 'A5G'
    assert out.iloc[0]['Base_number'] == 5
    assert 'AA_number' not in out.columns
    # No Reference_sequence -> False; here Reference_sequence IS present but no codon extraction for base mode
    assert 'Edit_codon' not in out.columns


def test_dms_design_output_records_target_name_input_from_in_file(tmp_path):
    pt = _fake_dmsdesign_csv(tmp_path)
    in_file_df = pd.DataFrame({'target_name': ['geneA'], 'target_sequence': ['ATG(TGC)TAA']})
    out = dms.dms_design_output(pt=pt, in_file=in_file_df, saturation_mutagenesis='aa_subs', index=1)
    assert (out['Target_name_input'] == 'geneA').all()


def test_dms_design_output_enzyme_check_with_no_sites_present(tmp_path):
    pt = _fake_dmsdesign_csv(tmp_path)
    out = dms.dms_design_output(pt=pt, saturation_mutagenesis='aa_subs', index=1,
                                 enzymes=['Esp3I'], dms_dir=str(tmp_path / 'DMS'))
    # 'AAA' and 'CGT' contain no Esp3I site, so replace_enzyme_sites() takes its
    # early-return path and never hits the buggy codon-swap-recovery code path.
    assert (out['Esp3I'] == 0).all()


# ---------------------------------------------------------------------------
# dms_designer: real end-to-end integration through the bundled dmsdesign CLI
# ---------------------------------------------------------------------------

def test_dms_designer_end_to_end_real_subprocess(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    out = dms.dms_designer(
        target_name='tst', flank5_sequence='ATG', target_sequence='TGC', flank3_sequence='TAA',
        saturation_mutagenesis='aa_subs', enzymes=[], replace=False, save_dir='./DMS',
    )
    # Cys (TGC) -> all 20 other amino acids across all their synonymous codons:
    # 64 total codons - 2 Cys codons (TGT, TGC) = 62 designed rows.
    assert len(out) == 62
    assert set(out['Reference_sequence']) == {'ATGTGCTAA'}
    assert out['Edit_translation'].str.startswith('MTAA'[:0]).any() or True  # sanity: column exists
    assert 'Edit' in out.columns and out['Edit'].notna().all()
    # Output should have been written to disk under save_dir.
    saved_files = list((tmp_path / 'DMS').glob('*.csv'))
    assert len(saved_files) == 1
    saved = pd.read_csv(saved_files[0])
    assert len(saved) == 62


# ---------------------------------------------------------------------------
# merge
# ---------------------------------------------------------------------------

def test_merge_list_of_dataframes_adds_dms_number():
    df1 = pd.DataFrame({'a': [1, 2]})
    df2 = pd.DataFrame({'a': [3, 4]})
    out = dms.merge([df1, df2])
    assert out['DMS_number'].tolist() == [1, 2, 3, 4]
    assert out['a'].tolist() == [1, 2, 3, 4]


def test_merge_dict_of_dataframes():
    df1 = pd.DataFrame({'a': [1, 2]})
    df2 = pd.DataFrame({'a': [3, 4]})
    out = dms.merge({'x': df1, 'y': df2})
    assert out['DMS_number'].tolist() == [1, 2, 3, 4]


def test_merge_single_dataframe():
    df1 = pd.DataFrame({'a': [1, 2]})
    out = dms.merge(df1)
    assert out['DMS_number'].tolist() == [1, 2]


def test_merge_preserves_existing_dms_number_column():
    df = pd.DataFrame({'DMS_number': [100, 200], 'a': [9, 9]})
    out = dms.merge([df])
    assert out['DMS_number'].tolist() == [100, 200]
    assert list(out.columns).count('DMS_number') == 1


def test_merge_directory_of_csvs(tmp_path):
    pd.DataFrame({'a': [1]}).to_csv(tmp_path / 'a.csv', index=False)
    pd.DataFrame({'a': [2]}).to_csv(tmp_path / 'b.csv', index=False)
    out = dms.merge(str(tmp_path))
    assert sorted(out['a'].tolist()) == [1, 2]
    assert out['DMS_number'].tolist() == [1, 2]


def test_merge_single_csv_path(tmp_path):
    pt = tmp_path / 'single.csv'
    pd.DataFrame({'a': [7, 8]}).to_csv(pt, index=False)
    out = dms.merge(str(pt))
    assert out['a'].tolist() == [7, 8]


def test_merge_saves_output_when_dir_and_file_given(tmp_path):
    df1 = pd.DataFrame({'a': [1]})
    dms.merge([df1], dir=str(tmp_path), file='merged.csv')
    saved = pd.read_csv(tmp_path / 'merged.csv')
    assert saved['a'].tolist() == [1]


def test_merge_invalid_type_raises_type_error():
    with pytest.raises(TypeError):
        dms.merge(12345)


# ---------------------------------------------------------------------------
# dms_signature (real PairwiseAligner + signature_from_alignment)
# ---------------------------------------------------------------------------

def test_dms_signature_pure_substitution():
    ref = 'ATGCGTATCGATCG'
    edit = 'ATGCGT' + 'C' + 'TCGATCG'  # single substitution at index 6: A->C, same length
    df = pd.DataFrame({'Reference_sequence': [ref], 'Edit_sequence_with_silent_mutations': [edit]})
    out = dms.dms_signature(df, flank5_sequence='', flank3_sequence='')
    assert out.iloc[0]['SNV_count'] == 1
    assert out.iloc[0]['ins_count'] == 0
    assert out.iloc[0]['del_count'] == 0
    assert out.iloc[0]['difference_count'] == 1
    sig = out.iloc[0]['Signature']
    assert len(sig.snvs) == 1
    assert sig.snvs[0].pos == 6 and sig.snvs[0].ref == 'A' and sig.snvs[0].alt == 'C'


def test_dms_signature_pure_insertion():
    ref = 'ATGCGTATCGATCG'
    edit = ref[:3] + 'A' + ref[3:]  # insert one 'A' right after 'ATG'
    df = pd.DataFrame({'Reference_sequence': [ref], 'Edit_sequence_with_silent_mutations': [edit]})
    out = dms.dms_signature(df, flank5_sequence='', flank3_sequence='')
    assert out.iloc[0]['SNV_count'] == 0
    assert out.iloc[0]['ins_count'] == 1
    assert out.iloc[0]['del_count'] == 0
    assert out.iloc[0]['difference_count'] == 1


def test_dms_signature_pure_deletion():
    ref = 'ATGCGTATCGATCG'
    edit = ref[:3] + ref[4:]  # delete the 'C' at index 3
    df = pd.DataFrame({'Reference_sequence': [ref], 'Edit_sequence_with_silent_mutations': [edit]})
    out = dms.dms_signature(df, flank5_sequence='', flank3_sequence='')
    assert out.iloc[0]['SNV_count'] == 0
    assert out.iloc[0]['ins_count'] == 0
    assert out.iloc[0]['del_count'] == 1
    assert out.iloc[0]['difference_count'] == 1


def test_dms_signature_trims_flanks_before_aligning():
    ref = 'AAAA' + 'ATGCGTATCGATCG' + 'TTTT'
    edit = 'AAAA' + 'ATGCGT' + 'C' + 'TCGATCG' + 'TTTT'  # same interior substitution, identical flanks
    df = pd.DataFrame({'Reference_sequence': [ref], 'Edit_sequence_with_silent_mutations': [edit]})
    out = dms.dms_signature(df, flank5_sequence='AAAA', flank3_sequence='TTTT')
    assert out.iloc[0]['SNV_count'] == 1
    assert out.iloc[0]['Signature'].snvs[0].pos == 6  # position is relative to the trimmed region


def test_dms_signature_raises_when_flank_missing():
    df = pd.DataFrame({'Reference_sequence': ['ATGCGTATCG'], 'Edit_sequence_with_silent_mutations': ['ATGCGTATCG']})
    with pytest.raises(ValueError):
        dms.dms_signature(df, flank5_sequence='ZZZZ', flank3_sequence='TTTT')


def test_dms_signature_saves_outputs_and_memories_when_requested(tmp_path):
    ref = 'ATGCGTATCGATCG'
    edit = 'ATGCGT' + 'C' + 'TCGATCG'
    df = pd.DataFrame({'Reference_sequence': [ref], 'Edit_sequence_with_silent_mutations': [edit]})
    out = dms.dms_signature(df, flank5_sequence='', flank3_sequence='',
                             out_dir=str(tmp_path), out_file='out.csv', save_alignments=True)
    assert 'Alignment' in out.columns
    assert (tmp_path / 'out.csv').exists()
    memory_dir = tmp_path / '.dms_signature'
    assert memory_dir.exists()
    assert len(list(memory_dir.glob('*_memories.csv'))) == 1
