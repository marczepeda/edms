'''
Module: clone.py
Author: Marc Zepeda
Created: 2024-05-29 
Description: Molecular cloning 

Usage:
[Individual GG Cloning]
- ord_form(): Sigma Alrich ordering formatter
- tb(): designs top & bottom oligonucleotides
- sgRNAs(): design GG cloning oligonucleotides for cutting and base editing sgRNAs
- epegRNAs(): design GG cloning oligonucleotides for prime editing epegRNAs
- ngRNAs(): design GG cloning oligonucleotides for prime editing ngRNAs

[Library GG Cloning]
- pe_pcr1: Dataframe of PE library PCR1 primers
- RE_typeIIS: Dataframe of TypeIIS restriction enzymes
- generate_sequences(): recursively generates all possible sequences of A, T, C, G of the specified length
- filter_GC(): filters sequences based on GC content
- shuffle(): randomly reorganizes a list
- pe_twist_oligos(): makes twist oligonucleotides for prime editing

[Master Mix]
- pcr_mm(): NEB Q5 PCR master mix calculations

[Simulation]
- pcr_sim(): returns dataframe with simulated pcr product 
'''

# Import packages
import pandas as pd
import numpy as np
import random
from Bio.Seq import Seq
from ..gen import tidy as t

# Individual GG cloning
def ord_form(df:pd.DataFrame,id:str,seq:str,suf:str,pre:str):
    ''' 
    ord_form(): Sigma Alrich ordering formatter
    
    Parameters:
    df (dataframe): pandas dataframe
    id (str): id column name
    seq (str): oligonucleotide sequence
    suf (str): suffix for oligonucleotide category
    pre (str): prefix for oligonucleotide category
    
    Dependencies: pandas
    '''
    ord = df[[(pre+id+suf),(pre+seq+suf)]]
    ord = ord.rename(columns={(pre+id+suf):'Oligo Name',(pre+seq+suf):'Sequence'})
    scale = []
    bp = []
    for s in ord['Sequence']:
        if len(s)<60: scale.append(0.025)
        else: scale.append(0.05)
        bp.append(len(s))
    ord['Scale (µmol)']=scale
    ord['bp']=bp
    return ord

def tb(df:pd.DataFrame,id:str,seq:str,t5:str,t3:str,b5:str,b3:str,tG:bool,pre:str):
    ''' 
    tb(): designs top & bottom oligonucleotides
    
    Parameters:
    df (datframe): Dataframe with sequences
    id (str): id column name
    seq (str): sequence column name
    t5 (str): top oligonucleotide 5' overhang
    t3 (str): top oligonucleotide 3' overhang
    b5 (str): bottom oligonucleotide 5' overhang
    b3 (str): bottom oligonucleotide 3' overhang
    tG (bool): add 5' G to spacer if needed
    pre (str): prefix for ids and id column

    Dependencies: pandas & Bio.Seq
    '''
    top_ids=[]
    bot_ids=[]
    top_seqs=[]
    bot_seqs=[]
    for i,s in enumerate(df[seq]):
        top_ids.append(pre+str(df.iloc[i][id])+'_t')
        bot_ids.append(pre+str(df.iloc[i][id])+'_b')
        if (tG==True)&(s[0]!='G'): top_seqs.append(t5+'G'+s+t3)
        else: top_seqs.append(t5+s+t3)
        bot_seqs.append(b5+str(Seq(s).reverse_complement()+b3))
    df[pre+id+'_t']=top_ids
    df[pre+id+'_b']=bot_ids
    df[pre+seq+'_t']=top_seqs
    df[pre+seq+'_b']=bot_seqs
    
    return df

def sgRNAs(df:pd.DataFrame,id:str,spacer='Spacer_sequence',t5='CACC',t3='',b5='AAAC',b3='',tG=True,order=True):
    ''' 
    sgRNAs(): design GG cloning oligonucleotides for cutting and base editing sgRNAs
    
    Parameters:
    df (dataframe): Dataframe with sequence information for epegRNAs
    id (str): id column name
    spacer (str): spacer column name (Default: Spacer_sequence)
    t5 (str): top oligonucleotide 5' overhang
    t3 (str): top oligonucleotide 3' overhang
    b5 (str): bottom oligonucleotide 5' overhang (revcom)
    b3 (str): bottom oligonucleotide 3' overhang (revcom)
    tG (bool): add 5' G to spacer if needed (Default: True)
    order (bool): order format
    
    Dependencies: pandas, top_bot(), & ord_form()
    '''
    df=tb(df=df,id=id,seq=spacer,t5=t5,t3=t3,b5=b5,b3=b3,tG=tG,pre='o') # Make top and bottom oligos for spacer inserts
    if order==True: return pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_t',pre='o'), # Sigma order format
                                      ord_form(df=df,id=id,seq=spacer,suf='_b',pre='o')]).reset_index(drop=True)
    else: return df # Original dataframe with top and bottom oligos

def epegRNAs(df: pd.DataFrame,id: str,tG=True, order=True,make_extension=True,
             spacer='Spacer_sequence',spacer_t5='CACC',spacer_t3='GTTTAAGAGC',spacer_b5='',spacer_b3='',
             extension='Extension_sequence',RTT='RTT_sequence',PBS='PBS_sequence',linker='Linker_sequence',extension_t5='',extension_t3='',extension_b5='CGCG',extension_b3='GCACCGACTC'):
    ''' 
    epegRNAs(): design GG cloning oligonucleotides for prime editing epegRNAs
    
    Parameters:
    df (dataframe): Dataframe with sequence information for epegRNAs
    id (str): id column name
    tG (bool, optional): add 5' G to spacer if needed (Default: True)
    order (bool, optional): order format (Default: True)
    make_extension (bool, optional): concatenate RTT, PBS, and linker to make extension sequence (Default: True)
    spacer (str, optional): epegRNA spacer column name (Default: Spacer_sequence)
        _t5 (str, optional): top oligonucleotide 5' overhang
        _t3 (str, optional): top oligonucleotide 3' overhang
        _b5 (str, optional): bottom oligonucleotide 5' overhang
        _b3 (str, optional): bottom oligonucleotide 3' overhang
    extension (str, optional): epegRNA extension name (Default: Extension_sequence)
        _t5 (str, optional): top oligonucleotide 5' overhang
        _t3 (str, optional): top oligonucleotide 3' overhang
        _b5 (str, optional): bottom oligonucleotide 5' overhang
        _b3 (str, optional): bottom oligonucleotide 3' overhang
    RTT (str, optional): epegRNA reverse transcripase template column name (Default: RTT_sequence)
    PBS (str, optional): epegRNA primer binding site column name (Default: PBS_sequence)
    linker (str, optional): epegRNA linker column name(Default: Linker_sequence)
    
    Assumptions:
    1. epegRNA scaffold: GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
    2. epegRNA motif: tevoPreQ1 (CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA)
    
    Dependencies: pandas, top_bot(), & ord_form()
    ''' 
    if make_extension==True: df[extension] = df[RTT]+df[PBS]+df[linker] # Make extension by concatenating RTT, PBS, and linker
    else: print(f'Warning: Did not make extension sequence!\nMake sure "{extension}" column includes RTT+PBS+linker for epegRNAs.')
    df=tb(df=df,id=id,seq=spacer,t5=spacer_t5,t3=spacer_t3,b5=spacer_b5,b3=spacer_b3,tG=tG,pre='ps_') # Make top and bottom oligos for spacer inserts
    df=tb(df=df,id=id,seq=extension,t5=extension_t5,t3=extension_t3,b5=extension_b5,b3=extension_b3,tG=False,pre='pe_') # Make top and bottom oligos for extension inserts
    if order==True: return pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_t',pre='ps_'), # Sigma order format
                                      ord_form(df=df,id=id,seq=spacer,suf='_b',pre='ps_'),
                                      ord_form(df=df,id=id,seq=extension,suf='_t',pre='pe_'),
                                      ord_form(df=df,id=id,seq=extension,suf='_b',pre='pe_')]).reset_index(drop=True)
    else: return df # Original dataframe with top and bottom oligos

def ngRNAs(df: pd.DataFrame,id: str,tG=True, order=True,
           spacer='Spacer_sequence',ngRNA_sp_t5='CACC',ngRNA_sp_t3='GTTTAAGAGC',ngRNA_sp_b5='',ngRNA_sp_b3=''):
    ''' 
    ngRNAs(): design GG cloning oligonucleotides for prime editing ngRNAs
    
    Parameters:
    df (dataframe): Dataframe with spacers
    id (str): id column
    tG (bool, optional): add 5' G to spacer if needed (Default: True)
    order (bool, optional): order format (Default: True)
    spacer (str, optional): ngRNA spacer column name
    ngRNA_sp_t5 (str, optional): top oligonucleotide 5' overhang
    ngRNA_sp_t3 (str, optional): top oligonucleotide 3' overhang
    ngRNA_sp_b5 (str, optional): bottom oligonucleotide 5' overhang
    ngRNA_sp_b3 (str, optional): bottom oligonucleotide 3' overhang}
    
    Assumptions:
    1. ngRNA scaffold: GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
    
    Dependencies: pandas, top_bot(), & ord_form()
'''
    df=tb(df=df,id=id,seq=spacer,t5=ngRNA_sp_t5,t3=ngRNA_sp_t3,b5=ngRNA_sp_b5,b3=ngRNA_sp_b3,tG=tG,pre='ns_') # Make top and bottom oligos for spacer inserts
    if order==True: return pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_t',pre='ns_'), # Sigma order format
                                      ord_form(df=df,id=id,seq=spacer,suf='_b',pre='ns_')]).reset_index(drop=True)
    
    else: return df # Original dataframe with top and bottom oligos

# Library GG cloning
''' pe_pcr1: Dataframe of PE library PCR1 primers
'''
pe_pcr1 = pd.DataFrame({
  0: {
    'Subpool Number': 1,
    'Forward Barcode': 'CGAATTCCGT',
    'Reverse Barcode': 'GCTTAGAACA',
    'Homology Arm 5': 'GGTCTCTGTTT',
    'Homology Arm 3': 'CGCGGGAGACC',
    'Forward PCR1 Template': 'CGAATTCCGTGGTCTCTGTTT',
    'Reverse PCR1 Template': 'TGTTCTAAGCGGTCTCCCGCG',
    'Forward PCR1 Primer ID': 'fpMUZ187.1',
    'Forward PCR1 Primer': 'CGAATTCCGTGGTCTCTG',
    'Reverse PCR1 Primer ID': 'rpMUZ187.1',
    'Reverse PCR1 Primer': 'TGTTCTAAGCGGTCTCCC',
    'Forward PCR1 Primer NEB Tm': 63,
    'Reverse PCR1 Primer NEB Tm': 64,
    'PCR1 Primers NEB Ta': 63,
    'Forward PCR1 Primer IDT Tm': 59.5,
    'Reverse PCR1 Primer IDT Tm': 60.6,
    'Forward PCR1 Primer IDT Th': 31.7,
    'Reverse PCR1 Primer IDT Th': 17.5
  },
  1: {
    'Subpool Number': 2,
    'Forward Barcode': 'GTTTATCGGG',
    'Reverse Barcode': 'ACTGACTGTA',
    'Homology Arm 5': 'GGTCTCTGTTT',
    'Homology Arm 3': 'CGCGGGAGACC',
    'Forward PCR1 Template': 'GTTTATCGGGGGTCTCTGTTT',
    'Reverse PCR1 Template': 'TACAGTCAGTGGTCTCCCGCG',
    'Forward PCR1 Primer ID': 'fpMUZ187.2',
    'Forward PCR1 Primer': 'GTTTATCGGGGGTCTCTG',
    'Reverse PCR1 Primer ID': 'rpMUZ187.2',
    'Reverse PCR1 Primer': 'TACAGTCAGTGGTCTCCC',
    'Forward PCR1 Primer NEB Tm': 62,
    'Reverse PCR1 Primer NEB Tm': 64,
    'PCR1 Primers NEB Ta': 63,
    'Forward PCR1 Primer IDT Tm': 58.7,
    'Reverse PCR1 Primer IDT Tm': 59.6,
    'Forward PCR1 Primer IDT Th': 2.8,
    'Reverse PCR1 Primer IDT Th': 36.4
  },
  2: {
    'Subpool Number': 3,
    'Forward Barcode': 'ACCGATATGG',
    'Reverse Barcode': 'CTCGTCATAG',
    'Homology Arm 5': 'GGTCTCTGTTT',
    'Homology Arm 3': 'CGCGGGAGACC',
    'Forward PCR1 Template': 'ACCGATATGGGGTCTCTGTTT',
    'Reverse PCR1 Template': 'CTATGACGAGGGTCTCCCGCG',
    'Forward PCR1 Primer ID': 'fpMUZ187.3',
    'Forward PCR1 Primer': 'ACCGATATGGGGTCTCTG',
    'Reverse PCR1 Primer ID': 'rpMUZ187.3',
    'Reverse PCR1 Primer': 'CTATGACGAGGGTCTCCC',
    'Forward PCR1 Primer NEB Tm': 63,
    'Reverse PCR1 Primer NEB Tm': 63,
    'PCR1 Primers NEB Ta': 63,
    'Forward PCR1 Primer IDT Tm': 59.8,
    'Reverse PCR1 Primer IDT Tm': 59.9,
    'Forward PCR1 Primer IDT Th': 46.2,
    'Reverse PCR1 Primer IDT Th': 45.3
  },
  3: {
    'Subpool Number': 4,
    'Forward Barcode': 'GAGGTCTTTC',
    'Reverse Barcode': 'CACAACATAC',
    'Homology Arm 5': 'GGTCTCTGTTT',
    'Homology Arm 3': 'CGCGGGAGACC',
    'Forward PCR1 Template': 'GAGGTCTTTCGGTCTCTGTTT',
    'Reverse PCR1 Template': 'GTATGTTGTGGGTCTCCCGCG',
    'Forward PCR1 Primer ID': 'fpMUZ187.4',
    'Forward PCR1 Primer': 'GAGGTCTTTCGGTCTCTG',
    'Reverse PCR1 Primer ID': 'rpMUZ187.4',
    'Reverse PCR1 Primer': 'GTATGTTGTGGGTCTCCC',
    'Forward PCR1 Primer NEB Tm': 62,
    'Reverse PCR1 Primer NEB Tm': 63,
    'PCR1 Primers NEB Ta': 63,
    'Forward PCR1 Primer IDT Tm': 58.5,
    'Reverse PCR1 Primer IDT Tm': 58.9,
    'Forward PCR1 Primer IDT Th': 27.1,
    'Reverse PCR1 Primer IDT Th': 42.6
  },
  4: {
    'Subpool Number': 5,
    'Forward Barcode': 'GAGTAGCTCA',
    'Reverse Barcode': 'ATGTACCCAA',
    'Homology Arm 5': 'GGTCTCTGTTT',
    'Homology Arm 3': 'CGCGGGAGACC',
    'Forward PCR1 Template': 'GAGTAGCTCAGGTCTCTGTTT',
    'Reverse PCR1 Template': 'TTGGGTACATGGTCTCCCGCG',
    'Forward PCR1 Primer ID': 'fpMUZ187.5',
    'Forward PCR1 Primer': 'GAGTAGCTCAGGTCTCTG',
    'Reverse PCR1 Primer ID': 'rpMUZ187.5',
    'Reverse PCR1 Primer': 'TTGGGTACATGGTCTCCC',
    'Forward PCR1 Primer NEB Tm': 61,
    'Reverse PCR1 Primer NEB Tm': 65,
    'PCR1 Primers NEB Ta': 63,
    'Forward PCR1 Primer IDT Tm': 57.6,
    'Reverse PCR1 Primer IDT Tm': 60.4,
    'Forward PCR1 Primer IDT Th': 39.5,
    'Reverse PCR1 Primer IDT Th': 42.5
  },
  5: {
    'Subpool Number': 6,
    'Forward Barcode': 'GGATGCATGA',
    'Reverse Barcode': 'CATCAAGCTT',
    'Homology Arm 5': 'GGTCTCTGTTT',
    'Homology Arm 3': 'CGCGGGAGACC',
    'Forward PCR1 Template': 'GGATGCATGAGGTCTCTGTTT',
    'Reverse PCR1 Template': 'AAGCTTGATGGGTCTCCCGCG',
    'Forward PCR1 Primer ID': 'fpMUZ187.6',
    'Forward PCR1 Primer': 'GGATGCATGAGGTCTCTG',
    'Reverse PCR1 Primer ID': 'rpMUZ187.6',
    'Reverse PCR1 Primer': 'AAGCTTGATGGGTCTCCC',
    'Forward PCR1 Primer NEB Tm': 62,
    'Reverse PCR1 Primer NEB Tm': 65,
    'PCR1 Primers NEB Ta': 63,
    'Forward PCR1 Primer IDT Tm': 59.0,
    'Reverse PCR1 Primer IDT Tm': 61.1,
    'Forward PCR1 Primer IDT Th': 31.9,
    'Reverse PCR1 Primer IDT Th': 42.6
  },
  6: {
    'Subpool Number': 7,
    'Forward Barcode': 'ATGAGGACGA',
    'Reverse Barcode': 'CACATAAAGG',
    'Homology Arm 5': 'GGTCTCTGTTT',
    'Homology Arm 3': 'CGCGGGAGACC',
    'Forward PCR1 Template': 'ATGAGGACGAGGTCTCTGTTT',
    'Reverse PCR1 Template': 'CCTTTATGTGGGTCTCCCGCG',
    'Forward PCR1 Primer ID': 'fpMUZ187.7',
    'Forward PCR1 Primer': 'ATGAGGACGAGGTCTCTG',
    'Reverse PCR1 Primer ID': 'rpMUZ187.7',
    'Reverse PCR1 Primer': 'CCTTTATGTGGGTCTCCC',
    'Forward PCR1 Primer NEB Tm': 63,
    'Reverse PCR1 Primer NEB Tm': 63,
    'PCR1 Primers NEB Ta': 63,
    'Forward PCR1 Primer IDT Tm': 59.7,
    'Reverse PCR1 Primer IDT Tm': 58.7,
    'Forward PCR1 Primer IDT Th': 34.0,
    'Reverse PCR1 Primer IDT Th': 42.6
  },
  7: {
    'Subpool Number': 8,
    'Forward Barcode': 'GGTAGACACG',
    'Reverse Barcode': 'TCGACTTAGA',
    'Homology Arm 5': 'GGTCTCTGTTT',
    'Homology Arm 3': 'CGCGGGAGACC',
    'Forward PCR1 Template': 'GGTAGACACGGGTCTCTGTTT',
    'Reverse PCR1 Template': 'TCTAAGTCGAGGTCTCCCGCG',
    'Forward PCR1 Primer ID': 'fpMUZ187.8',
    'Forward PCR1 Primer': 'GGTAGACACGGGTCTCTG',
    'Reverse PCR1 Primer ID': 'rpMUZ187.8',
    'Reverse PCR1 Primer': 'TCTAAGTCGAGGTCTCCC',
    'Forward PCR1 Primer NEB Tm': 64,
    'Reverse PCR1 Primer NEB Tm': 63,
    'PCR1 Primers NEB Ta': 63,
    'Forward PCR1 Primer IDT Tm': 60.4,
    'Reverse PCR1 Primer IDT Tm': 59.0,
    'Forward PCR1 Primer IDT Th': 52.1,
    'Reverse PCR1 Primer IDT Th': 24.3
  },
  8: {
    'Subpool Number': 9,
    'Forward Barcode': 'AGGCATGACT',
    'Reverse Barcode': 'GTTACAAGTC',
    'Homology Arm 5': 'GGTCTCTGTTT',
    'Homology Arm 3': 'CGCGGGAGACC',
    'Forward PCR1 Template': 'AGGCATGACTGGTCTCTGTTT',
    'Reverse PCR1 Template': 'GACTTGTAACGGTCTCCCGCG',
    'Forward PCR1 Primer ID': 'fpMUZ187.9',
    'Forward PCR1 Primer': 'AGGCATGACTGGTCTCTG',
    'Reverse PCR1 Primer ID': 'rpMUZ187.9',
    'Reverse PCR1 Primer': 'GACTTGTAACGGTCTCCC',
    'Forward PCR1 Primer NEB Tm': 64,
    'Reverse PCR1 Primer NEB Tm': 62,
    'PCR1 Primers NEB Ta': 63,
    'Forward PCR1 Primer IDT Tm': 60.7,
    'Reverse PCR1 Primer IDT Tm': 58.8,
    'Forward PCR1 Primer IDT Th': 31.9,
    'Reverse PCR1 Primer IDT Th': 30.9
  }
}).T

''' RE_typeIIS: Dataframe of TypeIIS restriction enzymes
'''
RE_typeIIS = pd.DataFrame({
  0: {
      'Name': 'Esp3I',
      'Sequence_t5': 'N(CGTCTC)N|N',
      'Sequence_b3': 'N(GCAGAG)NNNNN|N',
      'Recognition': 'CGTCTC',
      'Recognition_rc': 'GAGACG'
  },
  1: {
      'Name': 'BsaI',
      'Sequence_t5': 'N(GGTCTC)N|N',
      'Sequence_b3': 'N(CCAGAG)NNNNN|N',
      'Recognition': 'GGTCTC',
      'Recognition_rc': 'GAGACC'
  },
  2: {
      'Name': 'BspMI',
      'Sequence_t5': 'N(ACCTGC)NNNN|N',
      'Sequence_b3': 'N(TGGACG)NNNNNNNN|N',
      'Recognition': 'ACCTGC',
      'Recognition_rc': 'GCAGGT'
  }
}).T

def generate_sequences(length:int, current_sequence:str=""):
    """
    generate_sequences(): recursively generates all possible sequences of A, T, C, G of the specified length

    Parameters:
    length (int): the length of each unique molecular identifier in the list
    current_sequence (str, recursion): the current sequence being built for the unique molecular identifier
    
    Dependencies: generate_sequences()
    """
    bases = ['A', 'T', 'C', 'G']

    # Base case: if the current sequence length matches the desired length
    if len(current_sequence) == length: return [current_sequence]

    # Recursive case: extend the current sequence with each base and recurse
    sequences = []
    for base in bases: sequences.extend(generate_sequences(length, current_sequence + base))
    
    # Return final list containing all unique molecular identifiers of specified length
    return sequences

def filter_GC(sequences: list, GC_fract: tuple):
    '''
    filter_GC(): filters sequences based on GC content
    
    Parameters:
    sequences (list): list of sequences (str)
    GC_fract (tuple): Pair of GC content boundaries written fractions (Ex: (0.4,0.6))
    '''
    return [sequence for sequence in sequences if ((len(t.find_all(sequence,'G'))+len(t.find_all(sequence,'C')))/len(sequence)>GC_fract[0])&((len(t.find_all(sequence,'G'))+len(t.find_all(sequence,'C')))/len(sequence)<GC_fract[1])]

def shuffle(ls: list):
    """
    shuffle(): randomly reorganizes a list

    Parameters:
    ls (list): the list to be shuffled.
    
    Dependencies: random
    """
    ls2 = ls[:]  # Create a copy of the list to avoid modifying the original
    random.shuffle(ls2)
    return ls2

def pe_twist_oligos(df: pd.DataFrame,id_pre:str,tG=True, make_extension=True,UMI_length:int=8,UMI_GC_fract:tuple=(0.4,0.6),
                    fwd_barcode_t5='Forward Barcode',rev_barcode_t3='Reverse Barcode',
                    homology_arm_t5='Homology Arm 5',homology_arm_t3='Homology Arm 3',
                    ngRNA_hU6_gg_insert='GTTTAGAGACGATCGACGTCTCACACC',epegRNA_gg_insert='GTTTAAGAGCAGGTGCTAGACCTGCGTCGGTGC',
                    ngRNA_spacer='Spacer_sequence_ngRNA',epegRNA_spacer='Spacer_sequence_epegRNA',
                    epegRNA_extension='Extension_sequence',epegRNA_RTT='RTT_sequence',
                    epegRNA_PBS='PBS_sequence',epegRNA_linker='Linker_sequence',
                    epegRNA_pbs_length='PBS_length',ngRNA_group='ngRNA_group'):
    ''' 
    pe_twist_oligos(): makes twist oligonucleotides for prime editing
    
    Parameters:
    df (dataframe): Dataframe with sequence information for epegRNAs & corresponding ngRNAs
    id_pre (str): Prefix for ID column
    tG (bool, optional): add 5' G to spacer if needed (Default: True)
    UMI_length (int, optional): unique molecular identifier length (Default: 8; 4^8=65,536 UMIs before GC content filtering; adds 16 nt total)
    UMI_GC_fract (tuple, optional): unique molecular identifier GC content boundaries written fractions (Default: (0.4,0.6))
    make_extension (bool, optional): concatenate RTT, PBS, and linker to make extension sequence (Default: True)
    fwd_barcode_t5 (bool, optional): forward barcode column name (Default: Forward Barcode)
    rev_barcode_t3 (bool, optional): reverse barcode column name (Default: Reverse Barcode)
    homology_arm_t5 (bool, optional): homology arm t5 column name (Default: Homology Arm 5)
    homology_arm_t3 (bool, optional): homology arm t5 column name (Default: Homology Arm 3)
    ngRNA_hU6_gg_insert (str, optional): ngRNA scaffold to hU6 Golden Gate insert sequence (Default:  ngRNA_scafold_5nt - Esp3I(R) - random_nt - Esp3I(F) - hU6; GTTTAGAGACGATCGACGTCTCACACC)
    epegRNA_gg_insert (str, optional): epegRNA scaffold Golden Gate insert sequence (Default: epegRNA_scaffold_8nt - BspMI (R) - random_5nt - BspMI - epegRNA_scaffold_8nt; GTTTAAGAGCAGGTGCTAGACCTGCGTCGGTGC)
    ngRNA_spacer (str, optional): ngRNA spacer column name (Default: ngRNA_Spacer_sequence)
    epegRNA_spacer (str, optional): epegRNA spacer column name (Default: epegRNA_Spacer_sequence)
    epegRNA_extension (str, optional): epegRNA extension name (Default: Extension_sequence)
    epegRNA_RTT (str, optional): epegRNA reverse transcripase template column name (Default: RTT_sequence)
    epegRNA_PBS (str, optional): epegRNA primer binding site column name (Default: PBS_sequence)
    epegRNA_linker (str, optional): epegRNA linker column name (Default: Linker_sequence)
    
    Assumptions:
    1. Oligo Template: FWD Barcode - BsaI - mU6 - ngRNA_spacer - ngRNA_scaffold_5nt - Esp3I(R) - random_5nt - Esp3I(F) - hU6 - epegRNA_spacer - epegRNA_scaffold_8nt - BspMI (R) - random_5nt - BspMI - epegRNA_scaffold_8nt - epegRNA_extension - tevopreQ1_motif_5nt - BsaI - REV Barcode
    2. epegRNA motif: tevoPreQ1 (CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA)
    
    Dependencies: pandas,random,generate_sequences(),filter_GC(),shuffle(),pe_pcr1
    '''
    # Make extension by concatenating RTT, PBS, and linker
    if make_extension==True: df[epegRNA_extension] = df[epegRNA_RTT]+df[epegRNA_PBS]+df[epegRNA_linker]
    else: print(f'Warning: Did not make extension sequence!\nMake sure "{epegRNA_extension}" column includes RTT+PBS+linker for epegRNAs.')

    # Assign subpool barcodes to different PBS lengths and ngRNA groups
    pbs_lengths = sorted(df[epegRNA_pbs_length].value_counts().keys())
    ngRNA_groups = sorted(df[ngRNA_group].value_counts().keys())
    pbs_lengths_ls = [val for val in pbs_lengths for _ in np.arange(len(ngRNA_groups))] 
    ngRNA_groups_ls = ngRNA_groups*len(pbs_lengths)
    pcr1_barcodes = pe_pcr1.iloc[:len(pbs_lengths_ls)]
    pcr1_barcodes[epegRNA_pbs_length]=pbs_lengths_ls
    pcr1_barcodes[ngRNA_group]=ngRNA_groups_ls
    df = pd.merge(left=df,right=pcr1_barcodes,on=[epegRNA_pbs_length,ngRNA_group])

    # Assign UMI to each twist oligo
    UMI_sequences = filter_GC(generate_sequences(length=UMI_length),UMI_GC_fract)
    if len(UMI_sequences)<df.shape[0]: KeyError(f'UMI_length={UMI_length} results in {len(UMI_sequences)} UMIs, which is less than {df.shape[0]} twist oligonucleotides!')
    df['Forward UMI']=shuffle(ls=UMI_sequences)[:df.shape[0]]
    df['Reverse UMI']=shuffle(ls=UMI_sequences)[:df.shape[0]]

    # Make twist oligo & determine length
    df = df.sort_values(by=[epegRNA_pbs_length,ngRNA_group]).reset_index(drop=True)
    df[f'{epegRNA_spacer}_nt1']=[s[0] for s in df[epegRNA_spacer]]
    df[f'{ngRNA_spacer}_nt1']=[s[0] for s in df[ngRNA_spacer]]
    twist_oligos = []
    twist_products = []
    for i,(epegRNA_spacer_nt1,ngRNA_spacer_nt1) in enumerate(t.zip_cols(df,[f'{epegRNA_spacer}_nt1',f'{ngRNA_spacer}_nt1'])):
        if tG==True: # Append 5'G to spacer if not already present
            if epegRNA_spacer_nt1=='G' and ngRNA_spacer_nt1=='G':
                twist_oligos.append(df.iloc[i]['Forward UMI']+
                                    df.iloc[i][fwd_barcode_t5]+df.iloc[i][homology_arm_t5]+
                                    df.iloc[i][ngRNA_spacer]+ngRNA_hU6_gg_insert+
                                    df.iloc[i][epegRNA_spacer]+epegRNA_gg_insert+
                                    df.iloc[i][epegRNA_extension]+
                                    df.iloc[i][homology_arm_t3]+df.iloc[i][rev_barcode_t3]+
                                    df.iloc[i]['Reverse UMI'])
                twist_products.append(df.iloc[i][homology_arm_t5]+
                                      df.iloc[i][ngRNA_spacer]+ngRNA_hU6_gg_insert+
                                      df.iloc[i][epegRNA_spacer]+epegRNA_gg_insert+
                                      df.iloc[i][epegRNA_extension]+
                                      df.iloc[i][homology_arm_t3])
            elif ngRNA_spacer_nt1=='G':
                twist_oligos.append(df.iloc[i]['Forward UMI']+
                                    df.iloc[i][fwd_barcode_t5]+df.iloc[i][homology_arm_t5]+
                                    df.iloc[i][ngRNA_spacer]+ngRNA_hU6_gg_insert+
                                    'G'+df.iloc[i][epegRNA_spacer]+epegRNA_gg_insert+
                                    df.iloc[i][epegRNA_extension]+
                                    df.iloc[i][homology_arm_t3]+df.iloc[i][rev_barcode_t3]+
                                    df.iloc[i]['Reverse UMI'])
                twist_products.append(df.iloc[i][homology_arm_t5]+
                                      df.iloc[i][ngRNA_spacer]+ngRNA_hU6_gg_insert+
                                      'G'+df.iloc[i][epegRNA_spacer]+epegRNA_gg_insert+
                                      df.iloc[i][epegRNA_extension]+
                                      df.iloc[i][homology_arm_t3])
            elif epegRNA_spacer_nt1=='G':
                twist_oligos.append(df.iloc[i]['Forward UMI']+
                                    df.iloc[i][fwd_barcode_t5]+df.iloc[i][homology_arm_t5]+
                                    'G'+df.iloc[i][ngRNA_spacer]+ngRNA_hU6_gg_insert+
                                    df.iloc[i][epegRNA_spacer]+epegRNA_gg_insert+
                                    df.iloc[i][epegRNA_extension]+
                                    df.iloc[i][homology_arm_t3]+df.iloc[i][rev_barcode_t3]+
                                    df.iloc[i]['Reverse UMI'])
                twist_products.append(df.iloc[i][homology_arm_t5]+
                                      'G'+df.iloc[i][ngRNA_spacer]+ngRNA_hU6_gg_insert+
                                      df.iloc[i][epegRNA_spacer]+epegRNA_gg_insert+
                                      df.iloc[i][epegRNA_extension]+
                                      df.iloc[i][homology_arm_t3])
            else:
                twist_oligos.append(df.iloc[i]['Forward UMI']+
                                    df.iloc[i][fwd_barcode_t5]+df.iloc[i][homology_arm_t5]+
                                    'G'+df.iloc[i][ngRNA_spacer]+ngRNA_hU6_gg_insert+
                                    'G'+df.iloc[i][epegRNA_spacer]+epegRNA_gg_insert+
                                    df.iloc[i][epegRNA_extension]+
                                    df.iloc[i][homology_arm_t3]+df.iloc[i][rev_barcode_t3]+
                                    df.iloc[i]['Reverse UMI'])
                twist_products.append(df.iloc[i][homology_arm_t5]+
                                      'G'+df.iloc[i][ngRNA_spacer]+ngRNA_hU6_gg_insert+
                                      'G'+df.iloc[i][epegRNA_spacer]+epegRNA_gg_insert+
                                      df.iloc[i][epegRNA_extension]+
                                      df.iloc[i][homology_arm_t3])
        else: # Do not append 5'G to spacer if not already present
            twist_oligos.append(df.iloc[i]['Forward UMI']+
                                df.iloc[i][fwd_barcode_t5]+df.iloc[i][homology_arm_t5]+
                                df.iloc[i][ngRNA_spacer]+ngRNA_hU6_gg_insert+
                                df.iloc[i][epegRNA_spacer]+epegRNA_gg_insert+
                                df.iloc[i][epegRNA_extension]+
                                df.iloc[i][homology_arm_t3]+df.iloc[i][rev_barcode_t3]+
                                df.iloc[i]['Reverse UMI'])
            twist_products.append(df.iloc[i][homology_arm_t5]+
                                  df.iloc[i][ngRNA_spacer]+ngRNA_hU6_gg_insert+
                                  df.iloc[i][epegRNA_spacer]+epegRNA_gg_insert+
                                  df.iloc[i][epegRNA_extension]+
                                  df.iloc[i][homology_arm_t3])
    df['ngRNA_hU6_GG_insert'] = ngRNA_hU6_gg_insert
    df['epegRNA_GG_insert'] = epegRNA_gg_insert
    df['Twist_oligo'] = twist_oligos
    df['Twist_oligo_length']=[len(twist) for twist in df['Twist_oligo']]
    df['ID'] = [f'{id_pre}.{i+1}' for i in range(len(df))]

    # Check for 2 recognition sites per enzyme
    for (enzyme,recognition,recognition_rc) in t.zip_cols(df=RE_typeIIS,cols=['Name','Recognition','Recognition_rc']):
        enzyme_sites = [len(t.find_all(twist_product,recognition)) + # Check forward direction
                        len(t.find_all(twist_product,recognition_rc)) # Check reverse direction
                        for twist_product in twist_products] # Iterate through Twist product
        df[enzyme] = enzyme_sites
        print(f"{enzyme} does not have 2 reconition sites for {[id for (id,enzyme_site) in t.zip_cols(df=df,cols=['ID',enzyme]) if enzyme_site!=2]}")
    
    return df

# Master Mix
def pcr_mm(primers: pd.Series, template_uL: int, template='1-2 ng/uL template',
           Q5_mm_x_stock=5,dNTP_mM_stock=10,fwd_uM_stock=10,rev_uM_stock=10,Q5_U_uL_stock=2,
           Q5_mm_x_desired=1,dNTP_mM_desired=0.2,fwd_uM_desired=0.5,rev_uM_desired=0.5,Q5_U_uL_desired=0.02,
           total_uL=25,mm_x=1.1):
    '''
    pcr_mm(): NEB Q5 PCR master mix calculations
    
    Parameters:
    primers (Series): value_counts() for primers
    template_uL (int): template uL per reaction
    template (str, optional): template name (Default: '1-2 ng/uL template')
    Q5_mm_x_stock (int, optional): Q5 reaction master mix stock (Default: 5)
    dNTP_mM_stock (int, optional): [dNTP] stock in mM (Default: 10)
    fwd_uM_stock (int, optional): [FWD Primer] stock in mM (Default: 10)
    rev_uM_stock (int, optional): [REV Primer] stock in mM (Default: 10)
    Q5_U_uL_stock (int, optional): [Q5 Polymerase] stock in U/uL (Default: 2)
    Q5_mm_x_desired (int, optional): Q5 reaction master mix desired (Default: 1)
    dNTP_mM_desired (int, optional): [dNTP] desired in mM (Default: 0.2)
    fwd_uM_desired (float, optional): [FWD Primer] desired in mM (Default: 0.5)
    rev_uM_desired (float, optional): [REV Primer] desired in mM (Default: 0.5)
    Q5_U_uL_desired (float, optional): [Q5 Polymerase] desired in U/uL (Default: 0.02)
    total_uL (int, optional): total uL per reaction (Default: 25)
    mm_x (float, optional): master mix multiplier (Default: 1.1)

    Dependencies: pandas
    '''
    pcr_mm_dc = dict()
    for i,(pcr1_fwd,pcr1_rev) in enumerate(primers.keys()):
        pcr_mm_dc[(pcr1_fwd,pcr1_rev)] = pd.DataFrame({'Component':['Nuclease-free H2O',f'{Q5_mm_x_stock}x Q5 Reaction Buffer','dNTPs',pcr1_fwd,pcr1_rev,template,'Q5 Polymerase','Total'],
                                                       'Stock':['',Q5_mm_x_stock,dNTP_mM_stock,fwd_uM_stock,rev_uM_stock,'',Q5_U_uL_stock,''],
                                                       'Desired':['',Q5_mm_x_desired,dNTP_mM_desired,fwd_uM_desired,rev_uM_desired,'',Q5_U_uL_desired,''],
                                                       'Unit':['','x','mM','uM','uM','','U/uL',''],
                                                       'uL': [round(total_uL-sum([Q5_mm_x_desired/Q5_mm_x_stock,dNTP_mM_desired/dNTP_mM_stock,fwd_uM_desired/fwd_uM_stock,rev_uM_desired/rev_uM_stock,template_uL/total_uL,Q5_U_uL_desired/Q5_U_uL_stock]*total_uL),2),
                                                              round(Q5_mm_x_desired/Q5_mm_x_stock*total_uL,2),
                                                              round(dNTP_mM_desired/dNTP_mM_stock*total_uL,2),
                                                              round(fwd_uM_desired/fwd_uM_stock*total_uL,2),
                                                              round(rev_uM_desired/rev_uM_stock*total_uL,2),
                                                              round(template_uL,2),
                                                              round(Q5_U_uL_desired/Q5_U_uL_stock*total_uL,2),
                                                              round(total_uL,2)],
                                                       'uL MM': [round((total_uL-sum([Q5_mm_x_desired/Q5_mm_x_stock,dNTP_mM_desired/dNTP_mM_stock,fwd_uM_desired/fwd_uM_stock,rev_uM_desired/rev_uM_stock,template_uL/total_uL,Q5_U_uL_desired/Q5_U_uL_stock]*total_uL))*primers.iloc[i]*mm_x,2),
                                                                 round(Q5_mm_x_desired/Q5_mm_x_stock*total_uL*primers.iloc[i]*mm_x,2),
                                                                 round(dNTP_mM_desired/dNTP_mM_stock*total_uL*primers.iloc[i]*mm_x,2),
                                                                 round(fwd_uM_desired/fwd_uM_stock*total_uL*primers.iloc[i]*mm_x,2),
                                                                 round(rev_uM_desired/rev_uM_stock*total_uL*primers.iloc[i]*mm_x,2),
                                                                 round(template_uL*primers.iloc[i]*mm_x,2),
                                                                 round(Q5_U_uL_desired/Q5_U_uL_stock*total_uL*primers.iloc[i]*mm_x,2),
                                                                 round(total_uL*primers.iloc[i]*mm_x,2)]
                                                     })
    return pcr_mm_dc

# Simulation
def pcr_sim(df: pd.DataFrame,template_col: str, fwd_bind_col: str, rev_bind_col: str,
            fwd_ext_col: str=None, rev_ext_col: str=None, product_col='PCR Product'):
    '''
    pcr_sim(): returns dataframe with simulated pcr product 
    
    Parameters:
    df (dataframe): dataframe with template & primers
    template_col (str): template column name
    fwd_bind_col (str): fwd primer binding region column name 
    rev_bind_col (str): rev primer binding region column name 
    fwd_ext_col (str, optional): fwd primer extension region column name (Default: None)
    rev_ext_col (str, optional): rev primer extension region column name (Default: None)
    product_col (str, optional): pcr product column name (Default: 'PCR Product')

    Dependencies: pandas,Bio.Seq,tidy
    '''
    pcr_product_ls = []

    if fwd_ext_col is not None and rev_ext_col is not None: # FWD & REV primers have extension regions
        for (template,fwd_bind,rev_bind,fwd_ext,rev_ext) in t.zip_cols(df=df,cols=[template_col,fwd_bind_col,rev_bind_col,fwd_ext_col,rev_ext_col]):
            fwd = fwd_ext + fwd_bind
            rev = rev_ext + rev_bind
            rc_rev_bind = ''.join(Seq(rev_bind).reverse_complement())
            rc_rev = ''.join(Seq(rev).reverse_complement())
            pcr_product_ls.append(fwd+ # fwd primer
                                  template[template.find(fwd_bind)+len(fwd_bind):template.find(rc_rev_bind)]+ # template between primers
                                  rc_rev) # reverse complement of reverse primer
    
    elif fwd_ext_col is not None: # FWD primers have extension regions
        for (template,fwd_bind,rev_bind,fwd_ext) in t.zip_cols(df=df,cols=[template_col,fwd_bind_col,rev_bind_col,fwd_ext_col]):
            fwd = fwd_ext + fwd_bind
            rev = rev_bind
            rc_rev_bind = ''.join(Seq(rev_bind).reverse_complement())
            rc_rev = ''.join(Seq(rev).reverse_complement())
            pcr_product_ls.append(fwd+ # fwd primer
                                  template[template.find(fwd_bind)+len(fwd_bind):template.find(rc_rev_bind)]+ # template between primers
                                  rc_rev) # reverse complement of reverse primer
    
    elif rev_ext_col is not None: # REV primers have extension regions
        for (template,fwd_bind,rev_bind,rev_ext) in t.zip_cols(df=df,cols=[template_col,fwd_bind_col,rev_bind_col,rev_ext_col]):
            fwd = fwd_bind
            rev = rev_ext + rev_bind
            rc_rev_bind = ''.join(Seq(rev_bind).reverse_complement())
            rc_rev = ''.join(Seq(rev).reverse_complement())
            pcr_product_ls.append(fwd+ # fwd primer
                                  template[template.find(fwd_bind)+len(fwd_bind):template.find(rc_rev_bind)]+ # template between primers
                                  rc_rev) # reverse complement of reverse primer
    
    else: # FWD and REV primers do not have extension regions
        for (template,fwd_bind,rev_bind) in t.zip_cols(df=df,cols=[template_col,fwd_bind_col,rev_bind_col]):
            fwd = fwd_bind
            rev = rev_bind
            rc_rev_bind = ''.join(Seq(rev_bind).reverse_complement())
            rc_rev = ''.join(Seq(rev).reverse_complement())
            pcr_product_ls.append(fwd+ # fwd primer
                                  template[template.find(fwd_bind)+len(fwd_bind):template.find(rc_rev_bind)]+ # template between primers
                                  rc_rev) # reverse complement of reverse primer
            
    df[product_col]=pcr_product_ls
    return df