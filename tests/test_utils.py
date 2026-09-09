'''
Tests for edms.utils

Covers apply_filters(), the shared filtering step behind `edms pe pipeline` and
`edms dms pipeline`.
'''
import pandas as pd
import pytest

from edms import utils


def _df():
    return pd.DataFrame({
        'Edit': ['A1G', 'A1G', 'A2T', 'A3C'],
        'difference_count': [1, 1, 2, 3],
        'Reference_sequence': ['ACGTTCAAG', 'ACGTTCAAG', 'ACGTTCAAG', 'GGGGGGGGG'],
    })


def test_apply_filters_without_filters_returns_the_dataframe_untouched():
    df = _df()
    assert utils.apply_filters(df) is df


def test_apply_filters_applies_queries_in_order():
    out = utils.apply_filters(_df(), query=["Reference_sequence.str.contains('ACGTTCAAG')",
                                            'difference_count >= 2'])
    assert list(out['Edit']) == ['A2T']


def test_apply_filters_drop_duplicates_keeps_first_occurrence():
    out = utils.apply_filters(_df(), drop_duplicates=['Edit'])
    assert list(out['Edit']) == ['A1G', 'A2T', 'A3C']


def test_apply_filters_accepts_bare_strings_for_both_arguments():
    out = utils.apply_filters(_df(), query='difference_count == 1', drop_duplicates='Edit')
    assert list(out['Edit']) == ['A1G']


def test_apply_filters_resets_the_index():
    out = utils.apply_filters(_df(), query='difference_count >= 2')
    assert list(out.index) == [0, 1]


def test_apply_filters_reports_row_counts_per_step(capsys):
    utils.apply_filters(_df(), query=['difference_count >= 2'], drop_duplicates=['Edit'],
                        label='pegRNAs')
    out = capsys.readouterr().out
    assert 'Filtering: 4 pegRNAs' in out
    assert '4 -> 2 pegRNAs' in out


def test_apply_filters_raises_on_an_unknown_column_and_lists_what_is_available():
    with pytest.raises(ValueError) as excinfo:
        utils.apply_filters(_df(), query=['Esp3I_codon_swap_recovered == False'])
    message = str(excinfo.value)
    assert 'Filter query failed' in message
    assert 'Available columns' in message
    assert 'difference_count' in message


def test_apply_filters_raises_when_a_query_empties_the_dataframe():
    with pytest.raises(ValueError, match='removed every row'):
        utils.apply_filters(_df(), query=['difference_count > 99'])


def test_apply_filters_allow_empty_returns_an_empty_dataframe_instead():
    out = utils.apply_filters(_df(), query=['difference_count > 99'], allow_empty=True)
    assert len(out) == 0


def test_apply_filters_raises_on_a_missing_drop_duplicates_column():
    with pytest.raises(ValueError, match='not found'):
        utils.apply_filters(_df(), drop_duplicates=['Nope'])
