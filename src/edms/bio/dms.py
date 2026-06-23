'''
Module: dms.py
Author: Marc Zepeda
Created: 2026-06-19
Description: Deep Mutational Scan

Usage:
[Biological Dictionaries]
- dna_aa_codon_table: DNA to AA codon table
- aa_dna_codon_table: AA to DNA codon table

[Helper Functions]
- get_codons(): returns all codons within a specified frame for a nucleotide sequence
- get_codon_frames(): returns all codon frames for a nucleotide sequence
- found_list_in_order(): returns index of sub_ls found in consecutive order in main_ls or -1
- find_enzyme_sites(): find enzyme sites in DMS edit sequences or oligonucleotides
- enzyme_codon_swap(): synonymously modify DMS edit sequences to disrupt an RE recognition site
- replace_enzyme_sites(): save DMS sequences with enzyme sites before codon swap, save recovered silent-codon-swapped sequences after codon swap, and save unresolved sequences to lost

[DMSDesign]
- dms_design_input(): creates and checks DMSDesign saturation mutagenesis input file with optional silent mutation mode metadata
- dms_design(): run DMSDesign
- dms_design_output(): finishes annotations from DMSDesign output
- dms_designer(): execute DMSDesign for EDMS workflow
- merge(): combine one or more DMSDesign outputs into one edit-sequence library

[DMS]
- dms_signature(): create signatures for DMS outcomes using alignments
'''

# Import packages
import os
import re
import datetime
from edms import config
import pandas as pd
import numpy as np
from typing import Literal
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .signature import signature_from_alignment
from ..gen import io as io
from ..gen import tidy as t
from ..utils import memory_timer, load_resource_csv, mkdir

# Biological Dictionaries
''' dna_aa_codon_table: DNA to AA codon table '''
dna_aa_codon_table = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

''' aa_dna_codon_table: AA to DNA codon table '''
aa_dna_codon_table = {}
for _codon, _aa in dna_aa_codon_table.items():
    aa_dna_codon_table.setdefault(_aa, []).append(_codon)

# Helper Functions
def get_codons(sequence: str, frame: int = 0) -> list[str]:
    ''' get_codons(): returns all codons within a specified frame for a nucleotide sequence. '''
    return [sequence[i:i + 3] for i in range(frame, len(sequence) - 2, 3)]


def get_codon_frames(sequence: str) -> list[list[str]]:
    ''' get_codon_frames(): returns all codon frames for a nucleotide sequence. '''
    return [get_codons(sequence, frame) for frame in range(3)]


def found_list_in_order(main_ls: list, sub_ls: list) -> int:
    ''' found_list_in_order(): returns index of sub_ls found in consecutive order in main_ls or -1. '''
    if len(sub_ls) == 0:
        return -1
    for i in range(len(main_ls) - len(sub_ls) + 1):
        if main_ls[i:i + len(sub_ls)] == sub_ls:
            return i
    return -1


def _read_df(df: pd.DataFrame | str, literal_eval: bool = True) -> pd.DataFrame:
    '''Read dataframe from path if needed.'''
    if isinstance(df, str):
        return io.get(pt=df, literal_eval=literal_eval)
    return df.copy()


def _clean_dna(sequence: str) -> str:
    '''Uppercase DNA sequence and remove whitespace.'''
    return re.sub(r'\s+', '', str(sequence)).upper()


def _gc_content(sequence: str) -> float:
    '''Return percent GC content.'''
    sequence = _clean_dna(sequence)
    if len(sequence) == 0:
        return np.nan
    return 100 * (sequence.count('G') + sequence.count('C')) / len(sequence)


def _translate(sequence: str, frame: int = 0) -> str:
    '''Translate an in-frame DNA sequence using dna_aa_codon_table.'''
    codons = get_codons(_clean_dna(sequence), frame=frame)
    return ''.join([dna_aa_codon_table.get(codon, '?') for codon in codons])


def find_enzyme_sites(df: pd.DataFrame | str, enzyme: str, sequence_col: str = 'Oligonucleotide',
                      RE_type_IIS_df: pd.DataFrame | str = None, literal_eval: bool = False) -> pd.DataFrame:
    '''
    find_enzyme_sites(): find enzyme sites in DMS edit sequences or oligonucleotides.

    Parameters:
    df (pd.DataFrame | str): DataFrame or file path
    enzyme (str): Enzyme name (e.g. Esp3I, BsaI, BspMI)
    sequence_col (str, optional): column containing sequences to search (Default: Oligonucleotide)
    RE_type_IIS_df (pd.DataFrame | str, optional): Type IIS RE information
    literal_eval (bool, optional): convert string representations when loading files
    '''
    df = _read_df(df=df, literal_eval=literal_eval)
    if sequence_col not in df.columns:
        if 'Edit_sequence_with_silent_mutations' in df.columns:
            sequence_col = 'Edit_sequence_with_silent_mutations'
        elif 'Edit_sequence' in df.columns:
            sequence_col = 'Edit_sequence'
        else:
            raise ValueError(f'{sequence_col} not found in df columns.')

    if isinstance(RE_type_IIS_df, str):
        RE_type_IIS_df = io.get(pt=RE_type_IIS_df, literal_eval=literal_eval)
    elif RE_type_IIS_df is None:
        RE_type_IIS_df = load_resource_csv(filename='RE_type_IIS.csv')

    rec = RE_type_IIS_df[RE_type_IIS_df['Name'] == enzyme]['Recognition'].values[0]
    rec_rc = RE_type_IIS_df[RE_type_IIS_df['Name'] == enzyme]['Recognition_rc'].values[0]
    enzyme_sites_fwd = [t.find_all(_clean_dna(seq), rec) for seq in df[sequence_col]]
    enzyme_sites_rc = [t.find_all(_clean_dna(seq), rec_rc) for seq in df[sequence_col]]
    df[enzyme] = [len(fwd) + len(rc) for fwd, rc in zip(enzyme_sites_fwd, enzyme_sites_rc)]
    df[f'{enzyme}_fwd_i'] = enzyme_sites_fwd
    df[f'{enzyme}_rc_i'] = enzyme_sites_rc
    return df


def enzyme_codon_swap(df: pd.DataFrame | str, enzyme: str, sequence_col: str = 'Edit_sequence_with_silent_mutations',
                      reference_col: str = 'Reference_sequence', RE_type_IIS_df: pd.DataFrame | str = None,
                      out_dir: str = None, out_file: str = None, return_df: bool = True,
                      literal_eval: bool = False, comments: bool = False) -> pd.DataFrame:
    '''
    enzyme_codon_swap(): synonymously modify DMS edit sequences to disrupt an RE recognition site.

    Notes:
    - This assumes `sequence_col` is coding and starts in-frame.
    - Codons overlapping non-synonymous intended mutation positions are not modified.
    - If no synonymous codon removes the site, the sequence is left unchanged.
    '''
    df = _read_df(df=df, literal_eval=literal_eval)
    if sequence_col not in df.columns:
        if 'Edit_sequence_with_silent_mutations' in df.columns:
            sequence_col = 'Edit_sequence_with_silent_mutations'
        elif 'Edit_sequence' in df.columns:
            sequence_col = 'Edit_sequence'
        else:
            raise ValueError(f'{sequence_col} not found in df columns.')

    df = find_enzyme_sites(df=df, enzyme=enzyme, sequence_col=sequence_col,
                           RE_type_IIS_df=RE_type_IIS_df, literal_eval=literal_eval)
    changed_sequences = []
    annotations = []
    for _, row in df.iterrows():
        seq = _clean_dna(row[sequence_col])
        ref = _clean_dna(row[reference_col]) if reference_col in df.columns else seq
        if row[enzyme] == 0:
            changed_sequences.append(seq)
            annotations.append('')
            continue

        mutation_positions = {i for i, (a, b) in enumerate(zip(ref, seq)) if a != b}
        best_seq = seq
        best_annot = ''
        for codon_start in range(0, len(seq) - 2, 3):
            codon = seq[codon_start:codon_start + 3]
            if codon not in dna_aa_codon_table:
                continue
            if any(codon_start <= pos < codon_start + 3 for pos in mutation_positions):
                continue
            aa = dna_aa_codon_table[codon]
            for codon_alt in aa_dna_codon_table.get(aa, []):
                if codon_alt == codon:
                    continue
                trial = seq[:codon_start] + codon_alt + seq[codon_start + 3:]
                trial_df = pd.DataFrame({sequence_col: [trial]})
                trial_df = find_enzyme_sites(trial_df, enzyme=enzyme, sequence_col=sequence_col,
                                             RE_type_IIS_df=RE_type_IIS_df, literal_eval=literal_eval)
                if int(trial_df.iloc[0][enzyme]) == 0:
                    best_seq = trial
                    best_annot = f'{codon}{codon_start + 1}>{codon_alt}'
                    break
            if best_annot:
                break
        if comments and best_annot:
            print(f'{enzyme} disrupted by {best_annot}')
        changed_sequences.append(best_seq)
        annotations.append(best_annot)

    df[sequence_col] = changed_sequences
    df[f'{enzyme}_codon_swap_annotation'] = annotations
    if out_dir is not None and out_file is not None:
        io.save(obj=df, dir=out_dir, file=out_file)
    if return_df:
        return df

def replace_enzyme_sites(df: pd.DataFrame | str, enzyme: str,
                         sequence_col: str = 'Edit_sequence_with_silent_mutations',
                         reference_col: str = 'Reference_sequence',
                         dms_dir: str = './DMS',
                         RE_type_IIS_df: pd.DataFrame | str = None,
                         literal_eval: bool = False,
                         comments: bool = False) -> pd.DataFrame:
    '''
    replace_enzyme_sites(): save DMS sequences with enzyme sites before codon swap,
    save recovered silent-codon-swapped sequences after codon swap, and save unresolved
    sequences to lost.
    '''
    df = _read_df(df=df, literal_eval=literal_eval)

    df = find_enzyme_sites(
        df=df, enzyme=enzyme, sequence_col=sequence_col,
        RE_type_IIS_df=RE_type_IIS_df, literal_eval=literal_eval
    )

    df_clean = df[df[enzyme] == 0].copy()
    df_enzyme_before = df[df[enzyme] != 0].copy()

    if len(df_enzyme_before) == 0:
        return df

    before_dir = f'{dms_dir}/{enzyme}/codon_swap_before'
    after_dir = f'{dms_dir}/{enzyme}/codon_swap_after'
    lost_dir = f'{dms_dir}/{enzyme}/lost'

    io.save(
        obj=df_enzyme_before,
        dir=before_dir,
        file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_{enzyme}_before.csv'
    )

    df_swap = enzyme_codon_swap(
        df=df_enzyme_before,
        enzyme=enzyme,
        sequence_col=sequence_col,
        reference_col=reference_col,
        RE_type_IIS_df=RE_type_IIS_df,
        literal_eval=literal_eval,
        comments=comments,
    )

    df_swap[f'{enzyme}_codon_swap_recovered'] = df_swap[f'{enzyme}_codon_swap_annotation'].astype(str) != ''
    df_swap.loc[df_swap[f'{enzyme}_codon_swap_recovered'], sequence_col] = df_swap.loc[
        df_swap[f'{enzyme}_codon_swap_recovered'], f'{sequence_col}_{enzyme}_codon_swap'
    ]

    df_swap = find_enzyme_sites(
        df=df_swap,
        enzyme=enzyme,
        sequence_col=sequence_col,
        RE_type_IIS_df=RE_type_IIS_df,
        literal_eval=literal_eval,
    )

    df_after = df_swap[
        (df_swap[f'{enzyme}_codon_swap_recovered']) &
        (df_swap[enzyme] == 0)
    ].copy()

    df_lost = df_swap[df_swap[enzyme] != 0].copy()

    if len(df_after) > 0:
        io.save(
            obj=df_after,
            dir=after_dir,
            file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_{enzyme}_after.csv'
        )

    if len(df_lost) > 0:
        io.save(
            obj=df_lost,
            dir=lost_dir,
            file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_{enzyme}_lost.csv'
        )

    df_clean[f'{enzyme}_codon_swap_recovered'] = np.nan
    df_clean[f'{enzyme}_codon_swap_annotation'] = np.nan
    out = pd.concat([df_clean, df_after], ignore_index=True)

    out = find_enzyme_sites(
        df=out,
        enzyme=enzyme,
        sequence_col=sequence_col,
        RE_type_IIS_df=RE_type_IIS_df,
        literal_eval=literal_eval,
    )
    out = out[out[enzyme] == 0].reset_index(drop=True)

    if comments:
        print(f'{enzyme}: {len(df_enzyme_before)} DMS sequences had enzyme sites before codon swap.')
        print(f'{enzyme}: {len(df_after)} DMS sequences recovered by silent codon swap.')
        print(f'{enzyme}: {len(df_lost)} DMS sequences lost due to unresolved enzyme sites.')

    return out

# DMSDesign
def dms_design_input(target_name: str, flank5_sequence: str, target_sequence: str, flank3_sequence: str,
                     index: int = 1, silent_mutation: int = 0, silent_mutation_mode: str = 'close',
                     dir: str = '.', file: str = 'dms_design_input.csv'):
    '''
    dms_design_input(): creates and checks DMSDesign saturation mutagenesis input file.

    Parameters:
    target_name (str): name of target
    flank5_sequence (str): in-frame sequence 5' of saturation mutagenesis region
    target_sequence (str): in-frame saturation mutagenesis region
    flank3_sequence (str): in-frame sequence 3' of saturation mutagenesis region
    index (int, optional): 1st amino acid or base in target_sequence index (Default: 1)
    silent_mutation (int, optional): number of silent codon changes to add (Default: 0)
    silent_mutation_mode (str, optional): silent mutation placement mode; close, upstream, downstream, distribute, or barcode (Default: close)
    dir (str, optional): output directory
    file (str, optional): output filename
    '''
    if silent_mutation < 0:
        raise ValueError('silent_mutation must be an integer greater than or equal to 0.')
    if silent_mutation_mode not in ['close', 'upstream', 'downstream', 'distribute', 'barcode']:
        raise ValueError("silent_mutation_mode must be one of: 'close', 'upstream', 'downstream', 'distribute', 'barcode'.")
    if silent_mutation > 0:
        for name, seq in [('flank5_sequence', flank5_sequence), ('target_sequence', target_sequence), ('flank3_sequence', flank3_sequence)]:
            if len(seq) % 3 != 0:
                raise ValueError(f'Length of {name} ({len(seq)}) must be divisible by 3 when silent_mutation > 0.')

    io.save(obj=pd.DataFrame({'target_name': [target_name],
                              'target_sequence': [f'{flank5_sequence}({target_sequence}){flank3_sequence}'],
                              'index': [index],
                              'silent_mutation': [silent_mutation],
                              'silent_mutation_mode': [silent_mutation_mode]}),
            dir=dir, file=file)


def dms_design(file: str, saturation_mutagenesis: Literal['aa', 'aa_subs', 'aa_ins', 'aa_dels', 'aa_silent', 'base'] = None,
               silent_mutation: int = 0, silent_mutation_mode: str = 'close',
               out_dir: str = './DATETIMESTAMP_DMSDesign'):
    '''
    dms_design(): run DMSDesign.

    Parameters:
    file (str): input file with target_name,target_sequence columns
    saturation_mutagenesis (str, optional): aa, aa_subs, aa_ins, aa_dels, aa_silent, or base
    silent_mutation (int, optional): number of silent codon changes to add (Default: 0)
    silent_mutation_mode (str, optional): silent mutation placement mode; close, upstream, downstream, distribute, or barcode (Default: close)
    out_dir (str, optional): DMSDesign output directory
    '''
    if silent_mutation < 0:
        raise ValueError('silent_mutation must be an integer greater than or equal to 0.')
    if silent_mutation_mode not in ['close', 'upstream', 'downstream', 'distribute', 'barcode']:
        raise ValueError("silent_mutation_mode must be one of: 'close', 'upstream', 'downstream', 'distribute', 'barcode'.")
    cmd = f'python -m edms.bio.dmsdesign -f {file}'
    if saturation_mutagenesis is not None:
        cmd += f' -sat_mut {saturation_mutagenesis}'
    if silent_mutation != 0:
        cmd += f' -silent_mut {silent_mutation}'
        cmd += f' -silent_mut_mode {silent_mutation_mode}'
    if out_dir != './DATETIMESTAMP_DMSDesign':
        cmd += f' -out {out_dir}'
    print(cmd)
    os.system(cmd)


def _annotate_edit_from_name(target_name: str, saturation_mutagenesis: str, index: int = 1) -> dict:
    '''Infer Edit, AA_number/Base_number, and Edit_codon from DMSDesign target names.'''
    out = {'Edit': None, 'AA_number': None, 'Base_number': None, 'Edit_codon': None}
    if not isinstance(target_name, str):
        return out
    parts = target_name.split('_')
    if len(parts) < 3:
        return out
    number_token = parts[-2]
    change_token = parts[-1]
    if not re.fullmatch(r'-?\d+', number_token):
        return out
    number = int(number_token) + index - 1
    if 'to' in change_token:
        before, after = change_token.split('to', 1)
        before = before.replace('X', '*')
        after = after.replace('X', '*')
        out['Edit'] = f'{before}{number}{after}'
    if saturation_mutagenesis and 'aa' in saturation_mutagenesis:
        out['AA_number'] = number
    elif saturation_mutagenesis == 'base':
        out['Base_number'] = number
    return out


def dms_design_output(pt: str, in_file: pd.DataFrame | str = None, saturation_mutagenesis: str = None,
                      index: int = 1, enzymes: list[str] = None, sequence_col: str = 'Edit_sequence_with_silent_mutations',
                      replace: bool = True, dms_dir: str = './DMS', literal_eval: bool = False, comments: bool = False) -> pd.DataFrame:
    '''
    dms_design_output(): finishes annotations from DMSDesign output.

    Parameters:
    pt (str): path to DMSDesign output csv
    in_file (pd.DataFrame | str, optional): original DMSDesign input file
    saturation_mutagenesis (str, optional): aa, aa_subs, aa_ins, aa_dels, aa_silent, or base
    index (int, optional): 1st amino acid or base in target sequence index (Default: 1)
    enzymes (list[str], optional): Type IIS enzymes to check in edit sequences
    sequence_col (str, optional): sequence column used for enzyme checks
    replace (bool, optional): synonymously replace DMS edit sequences with enzyme sites and remove unresolved sequences (Default: True)
    dms_dir (str, optional): directory to save DMS files (Default: './DMS')
    literal_eval (bool, optional): convert string representations when loading files
    comments (bool, optional): print replacement summary comments (Default: False)
    '''
    df = io.get(pt=pt, literal_eval=literal_eval)
    if in_file is not None:
        in_file_df = _read_df(in_file, literal_eval=literal_eval)
        if saturation_mutagenesis is not None and 'target_name' in in_file_df.columns and len(in_file_df) == 1:
            df['Target_name_input'] = in_file_df.iloc[0]['target_name']

    if saturation_mutagenesis is not None:
        annots = [_annotate_edit_from_name(name, saturation_mutagenesis=saturation_mutagenesis, index=index)
                  for name in df['Target_name']]
        df['Edit'] = [a['Edit'] for a in annots]
        if saturation_mutagenesis and 'aa' in saturation_mutagenesis:
            df['AA_number'] = [a['AA_number'] for a in annots]
            if 'Target_sequence' in df.columns:
                df['Edit_codon'] = df['Edit'].astype(str) + '_' + df['Target_sequence'].str.extract(r'\(([ACGTacgt]{3}/[ACGTacgt]{3})\)')[0].fillna('')
        elif saturation_mutagenesis == 'base':
            df['Base_number'] = [a['Base_number'] for a in annots]

    seq_col = sequence_col if sequence_col in df.columns else 'Edit_sequence'
    df['Edit_sequence_length'] = [_clean_dna(seq).__len__() for seq in df[seq_col]]
    df['Edit_sequence_GC_content'] = [_gc_content(seq) for seq in df[seq_col]]
    if 'Reference_sequence' in df.columns:
        df['Reference_sequence_length'] = [_clean_dna(seq).__len__() for seq in df['Reference_sequence']]
        df['Reference_sequence_GC_content'] = [_gc_content(seq) for seq in df['Reference_sequence']]
        df['Reference_translation'] = [_translate(seq) for seq in df['Reference_sequence']]
        df['Edit_translation'] = [_translate(seq) for seq in df[seq_col]]

    if enzymes is not None:
        for enzyme in enzymes:
            df = find_enzyme_sites(
                df=df,
                enzyme=enzyme,
                sequence_col=seq_col,
                literal_eval=literal_eval,
            )

            if replace:
                df = replace_enzyme_sites(
                    df=df,
                    enzyme=enzyme,
                    sequence_col=seq_col,
                    reference_col='Reference_sequence',
                    dms_dir=dms_dir,
                    literal_eval=literal_eval,
                    comments=comments,
                )

    cols = ['Target_name', 'Internal_target_name', 'Design_number', 'Edit', 'AA_number', 'Base_number', 'Edit_type',
            'Saturation_mutagenesis', 'Reference_sequence', 'Edit_sequence', 'Edit_sequence_with_silent_mutations',
            'Silent_mutations_requested', 'Silent_mutation_mode', 'Silent_mutations_added', 'Silent_mutation_annotation',
            'Silent_mutation_positions', 'Silent_mutation_allowed_region',
            'Edit_sequence_length', 'Edit_sequence_GC_content', 'Reference_sequence_length', 'Reference_sequence_GC_content',
            'Reference_translation', 'Edit_translation', 'Edit_codon', 'Target_sequence']
    cols = [c for c in cols if c in df.columns]
    for c in df.columns:
        if c not in cols:
            cols.append(c)
    return df[cols]


def dms_designer(in_file: str = None, target_name: str = None, flank5_sequence: str = None,
                 target_sequence: str = None, flank3_sequence: str = None, index: int = 1,
                 silent_mutation: int = 0, silent_mutation_mode: str = 'barcode', saturation_mutagenesis: str = None,
                 enzymes: list[str] = ['Esp3I'], replace: bool = True,
                 out_dir: str = './DMSDesign/DATETIMESTAMP_DMSDesign',
                 save_dir: str = './DMS', save_file: str = None) -> pd.DataFrame:
    '''
    dms_designer(): execute DMSDesign for EDMS workflow.

    Parameters:
    in_file (str, optional): DMSDesign input file path
    target_name (str, optional): name of target (required if in_file not provided)
    flank5_sequence (str, optional): in-frame sequence 5' of saturation mutagenesis region (required if in_file not provided)
    target_sequence (str, optional): in-frame saturation mutagenesis region (required if in_file not provided)
    flank3_sequence (str, optional): in-frame sequence 3' of saturation mutagenesis region (required if in_file not provided)
    index (int, optional): 1st amino acid or base in target_sequence index (Default: 1)
    silent_mutation (int, optional): number of silent codon changes to add (Default: 0)
    silent_mutation_mode (str, optional): silent mutation placement mode; close, upstream, downstream, distribute, or barcode (Default: barcode)
    saturation_mutagenesis (str, optional): aa, aa_subs, aa_ins, aa_dels, aa_silent, or base
    enzymes (list[str], optional): Type IIS enzymes to check in edit sequences (Default: ['Esp3I'])
    replace (bool, optional): synonymously replace DMS edit sequences with enzyme sites and remove unresolved sequences (Default: True)
    out_dir (str, optional): DMSDesign output directory
    save_dir (str, optional): directory to save final annotated DMSDesign output
    save_file (str, optional): filename for final annotated DMSDesign output (if None, generated from target_name or in_file)

    *** Example saturation_mutagenesis.TXT file *** ---------------------------------------
    |											|
    |	target	ATGTGC(TGTGATGGTATGCCGGCGTAGTAA)TCGTAG   1                              |
    |											|
    ---------------------------------------------------------------------------------------

    *** Example saturation_mutagenesis.CSV file *** ---------------------------------------
    |											|
    |	target,ATGTGC(TGTGATGGTATGCCGGCGTAGTAA)TCGTAG,1		                        |
    |											|
    ---------------------------------------------------------------------------------------

    *** Example not_saturation_mutagenesis.TXT file *** -----------------------------------
    |											|
    |	target_01_substitution	ATGTGCTGTGATGGTAT(G/A)CCGGCGTAGTAATCGTAGC   1           |
    |	target_01_insertion	ATGTGCTGTGATGGTATG(+ATCTCGATGA)CCGGCGTAGTAATCGTAGC  1   |
    |	target_01_deletion	ATGTGCTGTGATGG(-TATGCCG)GCGTAGTAATCGTAGC    1           |
    |											|
    ---------------------------------------------------------------------------------------

    *** Example not_saturation_mutagenesis.CSV file *** -----------------------------------
    |											|
    |	target_01_substitution,ATGTGCTGTGATGGTAT(G/A)CCGGCGTAGTAATCGTAGC,1		|
    |	target_01_insertion,ATGTGCTGTGATGGTATG(+ATCTCGATGA)CCGGCGTAGTAATCGTAGC,1	|
    |	target_01_deletion,ATGTGCTGTGATGG(-TATGCCG)GCGTAGTAATCGTAGC,1			|
    |											|
    ---------------------------------------------------------------------------------------

    *** Formatting different DNA edits *** ------------------------------------------------
    |											|
    |	Substitution edit:	Format: (reference/edit)	Example:(G/A)		|
    |	Insertion edit:		Format: (+insertion)		Example:(+ATCG)		|
    |	Deletion edit:		Format: (-deletion)		Example:(-ATCG)		|
    |											|
    ---------------------------------------------------------------------------------------

    *** Combination edit example *** ------------------------------------------------------
    |											|
    |	Reference:			ATGCTGTGAT G TCGTGATG    A			|
    |	Edit:				A--CTGTGAT C TCGTGATGatcgA			|
    |	Sequence format:	A(-TG)CTGTGAT(G/C)TCGTGATG(+atcg)A			|
    |											|
    ---------------------------------------------------------------------------------------
    '''
    if silent_mutation < 0:
        raise ValueError('silent_mutation must be an integer greater than or equal to 0.')
    if silent_mutation_mode not in ['close', 'upstream', 'downstream', 'distribute', 'barcode']:
        raise ValueError("silent_mutation_mode must be one of: 'close', 'upstream', 'downstream', 'distribute', 'barcode'.")

    generated_input = False
    if in_file is None:
        if None in [target_name, flank5_sequence, target_sequence, flank3_sequence]:
            raise ValueError('Provide either in_file or target_name, flank5_sequence, target_sequence, and flank3_sequence.')
        in_file = f'./{"_".join(target_name.split(" "))}.csv'
        dms_design_input(target_name=target_name, flank5_sequence=flank5_sequence, target_sequence=target_sequence,
                         flank3_sequence=flank3_sequence, index=index, silent_mutation=silent_mutation,
                         silent_mutation_mode=silent_mutation_mode, dir='.', file=os.path.basename(in_file))
        generated_input = True
    elif silent_mutation > 0:
        in_file_df = io.get(pt=in_file)
        if saturation_mutagenesis is not None and len(in_file_df) > 0:
            seq = in_file_df.iloc[0]['target_sequence']
            flank5_sequence = seq.split('(')[0]
            target_sequence = seq.split('(')[1].split(')')[0]
            flank3_sequence = seq.split(')')[1]
            for name, value in [('flank5_sequence', flank5_sequence), ('target_sequence', target_sequence), ('flank3_sequence', flank3_sequence)]:
                if len(value) % 3 != 0:
                    raise ValueError(f'Length of {name} ({len(value)}) must be divisible by 3 when silent_mutation > 0.')
        else:
            print("Warning: Manually verify that target sequences in 'in_file' are in-frame when 'silent_mutation' > 0.")

    dms_design(file=in_file, saturation_mutagenesis=saturation_mutagenesis,
               silent_mutation=silent_mutation, silent_mutation_mode=silent_mutation_mode, out_dir=out_dir)

    search_dir = './DMSDesign' if out_dir == './DMSDesign/DATETIMESTAMP_DMSDesign' else out_dir
    output_candidates = sorted([os.path.join(search_dir, file) for file in io.relative_paths(search_dir) if 'DMSDesign.csv' in file], reverse=True)
    if len(output_candidates) == 0:
        raise FileNotFoundError(f'No DMSDesign.csv output found in {search_dir}.')
    pt = output_candidates[0]
    df = dms_design_output(pt=pt, in_file=in_file, saturation_mutagenesis=saturation_mutagenesis,
                           index=index, enzymes=enzymes, replace=replace, dms_dir=save_dir)

    if save_file is None:
        if target_name is None and not generated_input:
            target_name = os.path.splitext(os.path.basename(in_file))[0]
        save_file = f'{"_".join(str(target_name).split(" "))}_DMSDesign.csv'
    io.save(obj=df, dir=save_dir, file=save_file)
    return df


def merge(dms_designs: str | dict | list | pd.DataFrame, dir: str = None, file: str = None,
          literal_eval: bool = False) -> pd.DataFrame:
    '''
    merge(): combine one or more DMSDesign outputs into one edit-sequence library.

    Parameters:
    dms_designs (str | dict | list | pd.DataFrame): dataframe, path, directory, list of paths/dataframes, or dict of dataframes
    dir (str, optional): output directory
    file (str, optional): output filename
    literal_eval (bool, optional): convert string representations when loading files
    '''
    frames = []
    if isinstance(dms_designs, pd.DataFrame):
        frames = [dms_designs]
    elif isinstance(dms_designs, dict):
        frames = [value for value in dms_designs.values()]
    elif isinstance(dms_designs, list):
        frames = [io.get(pt=x, literal_eval=literal_eval) if isinstance(x, str) else x for x in dms_designs]
    elif isinstance(dms_designs, str):
        if os.path.isdir(dms_designs):
            paths = [os.path.join(dms_designs, x) for x in os.listdir(dms_designs) if x.endswith('.csv')]
            frames = [io.get(pt=x, literal_eval=literal_eval) for x in sorted(paths)]
        else:
            frames = [io.get(pt=dms_designs, literal_eval=literal_eval)]
    else:
        raise TypeError('dms_designs must be a dataframe, path, directory, list, or dict.')

    out = pd.concat(frames, ignore_index=True)
    if 'DMS_number' not in out.columns:
        out.insert(0, 'DMS_number', np.arange(1, len(out) + 1))
    if dir is not None and file is not None:
        io.save(obj=out, dir=dir, file=file)
    return out


# DMS
def dms_signature(df: pd.DataFrame | str, config_key: str=None,
                flank5_sequence: str = None, flank3_sequence: str = None, flank5_length: int=0, flank3_length: int=0,
                reference_sequence: str = 'Reference_sequence', edit_sequence: str = 'Edit_sequence_with_silent_mutations',
                match_score: float = 2, mismatch_score: float = -1, open_gap_score: float = -10, extend_gap_score: float = -0.1,
                out_dir: str=None, out_file: str=None, save_alignments: bool=False, return_df: bool=True, literal_eval: bool = False) -> pd.DataFrame:
    '''
    dms_signature(): create signatures for DMS outcomes using alignments.

    Parameters:
    df (pd.DataFrame | str): DMS design dataframe or path
    config_key (str, optional): config file key (FWD primer_REV primer) with 'motif5' (flank5_sequence) & 'motif3' (flank3_sequence)
    flank5_sequence (str, optional): flank5 sequence
    flank3_sequence (str, optional): flank3 sequence
    flank5_length (int, optional): length of flank5 sequence to include in alignment if provided (Default: 0)
    flank3_length (int, optional): length of flank3 sequence to include in alignment if provided (Default: 0)
    reference_sequence (str, optional): column name for reference sequences (Default: 'Reference_sequence')
    edit_sequence (str, optional): column name for post RTT sequences (Default: 'Edit_sequence_with_silent_mutations')
    match_score (float, optional): alignment match score (Default: 2)
    mismatch_score (float, optional): alignment mismatch score (Default: -1)
    open_gap_score (float, optional): alignment gap open score (Default: -10)
    extend_gap_score (float, optional): alignment gap extension score (Default: -0.1)
    out_dir (str, optional): output directory
    out_file (str, optional): output filename
    save_alignments (bool, optional): save alignments (Default: False, save memory)
    return_df (bool, optional): return dataframe (Default: True)
    literal_eval (bool, optional): convert string representations when loading files (Default: False)
    '''
    # Initialize timer; memory reporting
    memory_timer(reset=True)
    memories = []

    df = _read_df(df, literal_eval=literal_eval)

    # High sequence homology; punish gaps
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = match_score  # Score for a match
    aligner.mismatch_score = mismatch_score  # Penalty for a mismatch; applied to both strands
    aligner.open_gap_score = open_gap_score  # Penalty for opening a gap; applied to both strands
    aligner.extend_gap_score = extend_gap_score  # Penalty for extending a gap; applied to both strands

    # Use config_key to get flank5_sequence and flank3_sequence if provided
    if flank5_sequence is None and flank3_sequence is None:
        config_key = config.get_info(id=config_key)
        flank5_sequence = config_key['motif5']
        flank3_sequence = config_key['motif3']

    if save_alignments:
        alignments_list = []
    signatures_list = []
    for i, (reference_seq, edit_seq) in enumerate(t.zip_cols(df=df, cols=[reference_sequence, edit_sequence])):
        reference_seq = _clean_dna(reference_seq)
        edit_seq = _clean_dna(edit_seq)

        # Trim flanks if provided
        if flank5_sequence is not None and flank3_sequence is not None:
            if reference_seq.find(flank5_sequence)==-1 or reference_seq.rfind(flank3_sequence)==-1:
                raise(ValueError(f"Flank5 or Flank3 sequences were not found in reference sequence for pegRNAs row {i}.\nPlease check the flank5 ({flank5_sequence}) and flank3 ({flank3_sequence}) sequences.\nReference sequence: {reference_seq}"))
            if edit_seq.find(flank5_sequence)==-1 or edit_seq.rfind(flank3_sequence)==-1:
                raise(ValueError(f"Flank5 or Flank3 sequences were not found in edit sequence for pegRNAs row {i}.\nPlease check the flank5 ({flank5_sequence}) and flank3 ({flank3_sequence}) sequences.\nEdit sequence: {edit_seq}"))
            
            # Trim flanks from reference and edit sequences
            reference_seq = reference_seq[reference_seq.find(flank5_sequence)+len(flank5_sequence)-flank5_length : reference_seq.rfind(flank3_sequence)+flank3_length]
            edit_seq = edit_seq[edit_seq.find(flank5_sequence)+len(flank5_sequence)-flank5_length : edit_seq.rfind(flank3_sequence)+flank3_length]

        # Create and append alignment
        alignment = aligner.align(reference_seq, edit_seq)[0]
        if save_alignments:
            alignments_list.append(alignment)

        # Create and append signature
        signatures_list.append(signature_from_alignment(ref_seq=reference_seq, query_seq=edit_seq, alignment=alignment))
    
    # Create Alignment and Signature columns
    if save_alignments:
        df['Alignment'] = alignments_list
    df['Signature'] = signatures_list

    # Count # of SNVs, insertions, deletions in Signature
    snvs_ls = []
    ins_ls = []
    dels_ls = []
    for sign in df['Signature']: # Iterate through signatures

        # Count SNVs, insertions, deletions
        snvs_ls.append(len([snv for snv in sign.snvs]))
        ins_ls.append(sum([len(ind.ins) for ind in sign.indels]))
        dels_ls.append(sum([ind.dellen for ind in sign.indels]))

    # Create count columns
    df['SNV_count'] = snvs_ls
    df['ins_count'] = ins_ls
    df['del_count'] = dels_ls
    df['difference_count'] = df['SNV_count'] + df['ins_count'] + df['del_count']

    # Save & Return
    memories.append(memory_timer(task=f"dms_signature(): {len(df)} out of {len(df)}"))
    if out_dir is not None and out_file is not None:
        io.save(obj=pd.DataFrame(memories, columns=['Task','Memory, MB','Time, s']),
                dir=os.path.join(out_dir,f'.dms_signature'),
                file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_memories.csv')
        io.save(obj=df, dir=out_dir, file=out_file)
    if return_df==True: 
        return df