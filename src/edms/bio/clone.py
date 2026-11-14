'''
Module: clone.py
Author: Marc Zepeda
Created: 2024-05-29 
Description: Molecular cloning 

Usage:
[Individual GG Cloning]
- ord_form(): Sigma Alrich ordering formatter
- tb(): designs top & bottom oligonucleotides
- get_overhangs(): fill in GG cloning overhangs from the chosen promoter, scaffold, & epegRNA motif
- get_scaffold_oligos(): top & bottom oligonucleotides for ordering the scaffold
- sgRNAs(): design GG cloning oligonucleotides for cutting and base editing sgRNAs
- epegRNAs(): design GG cloning oligonucleotides for prime editing epegRNAs
- ngRNAs(): design GG cloning oligonucleotides for prime editing ngRNAs

[Library GG Cloning]
- epegRNA_pool(): design GG cloning oligonucleotides for pooled prime editing epegRNAs
- dms_pool(): design GG cloning oligonucleotides for deep mutational scanning (DMS) libraries

[UMIs]
- generate_sequences(): recursively generates all possible sequences of A, T, C, G of the specified length
- filter_GC(): filters sequences based on GC content
- shuffle(): randomly reorganizes a list
- encode_sequences(): Convert sequences to integer arrays for fast comparison.
- hamming_matrix(): pairwise Hamming distances between two encoded sequence sets
- min_hamming_to(): minimum Hamming distance from each sequence to any sequence in a reference set
- read_sequences(): read a list of DNA sequences from a file path, dataframe, series, or iterable
- min_hamming_within(): minimum Hamming distance from each sequence to any other sequence in the list
- gc_fractions(): GC content of each sequence as a fraction of its length
- max_homopolymers(): longest run of a single repeated base in each sequence
- fast_filter_by_hamming(): fastFilter sequences such that all retained sequences have a Hamming distance 'min_distance' from each other, using NumPy for speed.
- count_csv_rows(): counts number of rows in a csv file (subtracting 1 for header)
- umi_report(): print summary statistics for a set of UMIs
- umi(): generates unique molecular identifiers (UMIs) of specified length, GC content, and Hamming distance

[Master Mix]
- pcr_mm(): NEB Q5 PCR master mix calculations
- rxn_mm(): master mix for any reaction given stock and desired concentrations
- anneal_mm(): phosphorylate and anneal oligo duplexes
- golden_gate(): Esp3I/T4 ligase assembly from a vector and fragments
- gibson(): fragment volumes for an equimolar NEBuilder assembly

[Simulation]
- pcr_sim(): returns dataframe with simulated pcr product
- off_targets(): Find off-target sequences for a list of sequences using pairwise alignment.
'''

# Import packages
import pandas as pd
import numpy as np
import os
import random
import time
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import datetime
from ..gen import io
from ..gen import tidy as t
from ..bio import fastq as fq
from ..utils import load_resource_csv

# Individual GG cloning
def ord_form(df:pd.DataFrame, id:str, seq:str, suf:str, pre:str) -> pd.DataFrame:
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

def tb(df:pd.DataFrame, id:str, seq:str, t5:str, t3:str, 
       b5:str, b3:str, tG:bool, pre:str) -> pd.DataFrame:
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
        top_ids.append(pre+str(df.iloc[i][id])+'_top')
        bot_ids.append(pre+str(df.iloc[i][id])+'_bot')
        if (tG==True)&(s[0]!='G'): 
            top_seqs.append(t5+'G'+s+t3)
            bot_seqs.append(b5+str(Seq('G'+s).reverse_complement()+b3))
        else: 
            top_seqs.append(t5+s+t3)
            bot_seqs.append(b5+str(Seq(s).reverse_complement()+b3))
    df[pre+id+'_top']=top_ids
    df[pre+id+'_bot']=bot_ids
    df[pre+seq+'_top']=top_seqs
    df[pre+seq+'_bot']=bot_seqs
    
    return df

SCAFFOLD_SEQUENCES = {'SpCas9 Flip and Extend': 'GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC',
                      'SpCas9 Flip and Extend + com-modified': 'GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC',
                      'SpCas9 mod_7': 'GTTTCAGAGCTAGAAATAGCAAGTTGAAATAAGACTAGTTCGTTATCACGCCGAAAGGCGGGCACCGAGTCGGTGC'}

# Overhangs contributed by the promoter & scaffold for sgRNAs, whose spacer is cloned into a vector that already carries the scaffold
SGRNA_PROMOTER_OVERHANGS = {'hU6': {'t5': 'CACC', 'b3': ''},
                            '7SK': {'t5': 'TCTC', 'b3': ''}}
SGRNA_SCAFFOLD_OVERHANGS = {scaffold: {'t3': '', # The vector carries the scaffold, so the spacer only needs the scaffold's first 4 nt on the bottom oligo
                                       'b5': str(Seq(sequence[:4]).reverse_complement())}
                            for scaffold,sequence in SCAFFOLD_SEQUENCES.items()}

# Overhangs contributed by the promoter, scaffold, & epegRNA motif for peg/ngRNAs, whose scaffold is a separate GG fragment
PEGRNA_PROMOTER_OVERHANGS = {'hU6': {'spacer_t5': 'CACC', 'spacer_b3': ''},
                             '7SK': {'spacer_t5': 'TCTC', 'spacer_b3': ''}}
PEGRNA_SCAFFOLD_OVERHANGS = {scaffold: {'spacer_t3': sequence[:10], 'spacer_b5': '', # The spacer insert carries the scaffold's first 10 nt
                                        'extension_t5': '', 'extension_b3': str(Seq(sequence[-10:]).reverse_complement())} # The extension insert carries the scaffold's last 10 nt
                             for scaffold,sequence in SCAFFOLD_SEQUENCES.items()}
EPEGRNA_MOTIF_OVERHANGS = {'tevopreQ1': {'extension_b5': 'CGCG', 'extension_t3': ''},
                           'tevo2': {'extension_b5': 'GACG', 'extension_t3': ''}}

# ngRNAs cloned into a vector that already carries the scaffold take the sgRNA overhangs under the spacer_ names
NGRNA_VECTOR_PROMOTER_OVERHANGS = {promoter: {f'spacer_{overhang}': sequence for overhang,sequence in overhangs.items()}
                                   for promoter,overhangs in SGRNA_PROMOTER_OVERHANGS.items()}
NGRNA_VECTOR_SCAFFOLD_OVERHANGS = {scaffold: {f'spacer_{overhang}': sequence for overhang,sequence in overhangs.items()}
                                   for scaffold,overhangs in SGRNA_SCAFFOLD_OVERHANGS.items()}

def get_scaffold_oligos(scaffold: str, pre: str) -> pd.DataFrame:
    ''' 
    get_scaffold_oligos(): top & bottom oligonucleotides for ordering the scaffold
    
    Parameters:
    scaffold (str): scaffold name (i.e., SpCas9 Flip and Extend, SpCas9 mod_7)
    pre (str): prefix for the oligonucleotide names (i.e., pegRNA, ngRNA)
    
    Dependencies: pandas & Bio.Seq
    '''
    scaffold_names = {name.lower(): name for name in SCAFFOLD_SEQUENCES} # Match the scaffold name regardless of case
    if scaffold is None or scaffold.lower() not in scaffold_names:
        raise ValueError(f"'{scaffold}' is not an option; choose from {list(SCAFFOLD_SEQUENCES.keys())} to order the scaffold.")
    sequence = SCAFFOLD_SEQUENCES[scaffold_names[scaffold.lower()]]

    top = sequence[10:] # The spacer insert carries the first 10 nt
    if pre=='ngRNA': bot = 'AAAA'+str(Seq(sequence).reverse_complement()) # ngRNAs have no extension insert, so the scaffold carries the vector's overhang
    else: bot = str(Seq(sequence[:-10]).reverse_complement()) # The extension insert carries the last 10 nt

    return pd.DataFrame({'Oligo Name': [f'{pre}_scaffold_top',f'{pre}_scaffold_bot'],
                         'Sequence': [top,bot],
                         'Scale (µmol)': [0.05, 0.05],
                         'bp': [len(top),len(bot)]})

def get_overhangs(parts: list, **overhangs) -> dict:
    ''' 
    get_overhangs(): fill in GG cloning overhangs from the chosen promoter, scaffold, & epegRNA motif
    
    Parameters:
    parts (list): (name, options) tuples pairing each chosen part with its overhang options
        (i.e., [('hU6', SGRNA_PROMOTER_OVERHANGS), ('SpCas9 Flip and Extend', SGRNA_SCAFFOLD_OVERHANGS)])
    overhangs: overhangs keyed by name; None takes the chosen part's overhang, anything else overrides it
    
    Dependencies: none
    '''
    for name,options in parts:
        if name is None: continue # Part was not chosen; its overhangs stay empty unless specified
        option_names = {option.lower(): option for option in options} # Match the part name regardless of case
        if name.lower() not in option_names:
            raise ValueError(f"'{name}' is not an option; choose from {list(options.keys())} or specify the overhangs directly.")
        part = option_names[name.lower()]
        used = {overhang: sequence for overhang,sequence in options[part].items() if overhangs.get(overhang,'') is None} # Overhangs this part supplies that were not specified
        if used: print(f"{part} overhangs: " + ', '.join(f"{overhang}='{sequence}'" for overhang,sequence in used.items()))
        overhangs.update(used)

    return {overhang: '' if sequence is None else sequence for overhang,sequence in overhangs.items()}

def sgRNAs(df:pd.DataFrame | str,id:str, spacer: str='Spacer_sequence', 
           promoter: str='hU6', scaffold: str='SpCas9 Flip and Extend',
           t5: str=None, t3: str=None,b5: str=None,b3: str=None,tG: bool=True, 
           order: bool=True, file:str=None) -> pd.DataFrame:
    ''' 
    sgRNAs(): design GG cloning oligonucleotides for cutting and base editing sgRNAs
    
    Parameters:
    df (dataframe | str): Dataframe with sequence information for sgRNAs (or file path)
    id (str): id column name
    spacer (str): spacer column name (Default: Spacer_sequence)
    promoter (str, optional): promoter upstream of the spacer (Default: hU6)
        hU6: CACC t5, '' b3
        7SK: TCTC t5, '' b3
    scaffold (str, optional): scaffold downstream of the spacer (Default: SpCas9 Flip and Extend)
        SpCas9 Flip and Extend: '' t3, AAAC b5
        SpCas9 mod_7: '' t3, AAAC b5
    t5 (str): top oligonucleotide 5' overhang (Default: from promoter)
    t3 (str): top oligonucleotide 3' overhang (Default: from scaffold)
    b5 (str): bottom oligonucleotide 5' overhang (revcom) (Default: from scaffold)
    b3 (str): bottom oligonucleotide 3' overhang (revcom) (Default: from promoter)
    tG (bool): add 5' G to spacer if needed (Default: True)
    order (bool): order format
    file (str, optional): output file path
    
    Assumptions:
    1. sgRNA scaffold: SpCas9 Flip and Extend (GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC),
       SpCas9 Flip and Extend + com-modified [required for VLPs] (GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC),
       or SpCas9 mod_7 (GTTTCAGAGCTAGAAATAGCAAGTTGAAATAAGACTAGTTCGTTATCACGCCGAAAGGCGGGCACCGAGTCGGTGC)
    
    Dependencies: pandas, io, top_bot(), ord_form(), & get_overhangs()
    '''
    if type(df)==str: # Get sgRNAs dataframe from file path if needed
        df = io.get(pt=df)

    # Fill in the overhangs that were not specified from the chosen promoter & scaffold
    overhangs = get_overhangs(parts=[(promoter,SGRNA_PROMOTER_OVERHANGS),(scaffold,SGRNA_SCAFFOLD_OVERHANGS)],
                              t5=t5, t3=t3, b5=b5, b3=b3)

    df=tb(df=df,id=id,seq=spacer,t5=overhangs['t5'],t3=overhangs['t3'],b5=overhangs['b5'],b3=overhangs['b3'],tG=tG,pre='o') # Make top and bottom oligos for spacer inserts
    if order==True: # Sigma order format (or original dataframe with top and bottom oligos)
        df = pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_top',pre='o'), 
                        ord_form(df=df,id=id,seq=spacer,suf='_bot',pre='o')]).reset_index(drop=True)
    
    # Save & return dataframe
    io.save(obj=df, file=file)  
    return df

def epegRNAs(df: pd.DataFrame | str, id: str, tG: str=True, order: bool=True, make_extension: bool=True,
             promoter: str='7SK', scaffold: str='SpCas9 Flip and Extend', epegRNA_motif: str='tevo2',
             spacer: str='Spacer_sequence', spacer_t5: str=None, spacer_t3: str=None, 
             spacer_b5: str=None, spacer_b3: str=None, extension: str='Extension_sequence', RTT: str='RTT_sequence',
             PBS: str='PBS_sequence', linker: str='Linker_sequence', extension_t5: str=None, extension_t3: str=None,
             extension_b5: str=None, extension_b3: str=None,
             order_scaffold: bool=False, file:str=None) -> pd.DataFrame:
    ''' 
    epegRNAs(): design GG cloning oligonucleotides for prime editing epegRNAs
    
    Parameters:
    df (dataframe | str): Dataframe with sequence information for epegRNAs (or file path)
    id (str): id column name
    tG (bool, optional): add 5' G to spacer if needed (Default: True)
    order (bool, optional): order format (Default: True)
    make_extension (bool, optional): concatenate RTT, PBS, and linker to make extension sequence (Default: True)
    promoter (str, optional): promoter upstream of the spacer (Default: 7SK)
        hU6: CACC spacer_t5, '' spacer_b3
        7SK: TCTC spacer_t5, '' spacer_b3
    scaffold (str, optional): scaffold between the spacer & extension (Default: SpCas9 Flip and Extend)
        SpCas9 Flip and Extend: GTTTAAGAGC spacer_t3, '' spacer_b5, '' extension_t5, GCACCGACTC extension_b3
        SpCas9 mod_7: GTTTCAGAGC spacer_t3, '' spacer_b5, '' extension_t5, GCACCGACTC extension_b3
    epegRNA_motif (str, optional): 3' hairpin motif downstream of the extension (Default: tevo2)
        tevopreQ1: CGCG extension_b5, '' extension_t3
        tevo2: GACG extension_b5, '' extension_t3
    spacer (str, optional): epegRNA spacer column name (Default: Spacer_sequence)
        _t5 (str, optional): top oligonucleotide 5' overhang (Default: from promoter)
        _t3 (str, optional): top oligonucleotide 3' overhang (Default: from scaffold)
        _b5 (str, optional): bottom oligonucleotide 5' overhang (Default: from scaffold)
        _b3 (str, optional): bottom oligonucleotide 3' overhang (Default: from promoter)
    extension (str, optional): epegRNA extension name (Default: Extension_sequence)
        _t5 (str, optional): top oligonucleotide 5' overhang (Default: from scaffold)
        _t3 (str, optional): top oligonucleotide 3' overhang (Default: from epegRNA_motif)
        _b5 (str, optional): bottom oligonucleotide 5' overhang (Default: from epegRNA_motif)
        _b3 (str, optional): bottom oligonucleotide 3' overhang (Default: from scaffold)
    RTT (str, optional): epegRNA reverse transcripase template column name (Default: RTT_sequence)
    PBS (str, optional): epegRNA primer binding site column name (Default: PBS_sequence)
    linker (str, optional): epegRNA linker column name(Default: Linker_sequence)
    order_scaffold (bool, optional): order top and bottom oligonucleotide for scaffold sequence (Default: False)
    file (str, optional): output file path
    
    Assumptions:
    1. epegRNA scaffold: SpCas9 Flip and Extend (GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC)
       SpCas9 Flip and Extend + com-modified [required for VLPs] (GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC),
       or SpCas9 mod_7 (GTTTCAGAGCTAGAAATAGCAAGTTGAAATAAGACTAGTTCGTTATCACGCCGAAAGGCGGGCACCGAGTCGGTGC)
    2. epegRNA motif: tevo2 (CGTCGGTTCTACCTAGTTACGACGTTAAACCAACTAGGA) or tevopreQ1 (CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA)
    
    Dependencies: pandas, top_bot(), ord_form(), & get_overhangs()
    '''
    if type(df)==str: # Get epegRNAs dataframe from file path if needed
        df = io.get(pt=df)
    
    # Fill in the overhangs that were not specified from the chosen promoter, scaffold, & epegRNA motif
    overhangs = get_overhangs(parts=[(promoter,PEGRNA_PROMOTER_OVERHANGS),(scaffold,PEGRNA_SCAFFOLD_OVERHANGS),(epegRNA_motif,EPEGRNA_MOTIF_OVERHANGS)],
                              spacer_t5=spacer_t5, spacer_t3=spacer_t3, spacer_b5=spacer_b5, spacer_b3=spacer_b3,
                              extension_t5=extension_t5, extension_t3=extension_t3, extension_b5=extension_b5, extension_b3=extension_b3)
    
    if make_extension==True: 
        if df[RTT].isnull().any():
            print(f'Warning: Null values found in "{RTT}" column for some epegRNAs.')
        if df[PBS].isnull().any():
            print(f'Warning: Null values found in "{PBS}" column for some epegRNAs.')
        if df[linker].isnull().any():
            print(f'Warning: Null values found in "{linker}" column for some epegRNAs.')
        df[extension] = [rtt+pbs+linker for rtt,pbs,linker in zip(df[RTT].fillna(""), df[PBS].fillna(""), df[linker].fillna(""))] # Make extension by concatenating RTT, PBS, and linker
    else: print(f'Warning: Did not make extension sequence!\nMake sure "{extension}" column includes RTT+PBS+linker for epegRNAs.')
    df=tb(df=df,id=id,seq=spacer,t5=overhangs['spacer_t5'],t3=overhangs['spacer_t3'],b5=overhangs['spacer_b5'],b3=overhangs['spacer_b3'],tG=tG,pre='es_') # Make top and bottom oligos for spacer inserts
    df=tb(df=df,id=id,seq=extension,t5=overhangs['extension_t5'],t3=overhangs['extension_t3'],b5=overhangs['extension_b5'],b3=overhangs['extension_b3'],tG=False,pre='ee_') # Make top and bottom oligos for extension inserts
    if order==True: # Sigma order format (or original dataframe with top and bottom oligos)
        if order_scaffold==True: # Include the scaffold sequence in the order
            df = pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_top',pre='es_'),
                            ord_form(df=df,id=id,seq=spacer,suf='_bot',pre='es_'),
                            ord_form(df=df,id=id,seq=extension,suf='_top',pre='ee_'),
                            ord_form(df=df,id=id,seq=extension,suf='_bot',pre='ee_'),
                            get_scaffold_oligos(scaffold=scaffold, pre='pegRNA')
                            ]).reset_index(drop=True)
        else: # Do not include the scaffold sequence in the order
            df = pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_top',pre='es_'),
                            ord_form(df=df,id=id,seq=spacer,suf='_bot',pre='es_'),
                            ord_form(df=df,id=id,seq=extension,suf='_top',pre='ee_'),
                            ord_form(df=df,id=id,seq=extension,suf='_bot',pre='ee_')
                            ]).reset_index(drop=True)

    # Save & return dataframe
    io.save(obj=df, file=file)  
    return df

def ngRNAs(df: pd.DataFrame | str, id: str, tG: bool=True, order: bool=True, vector_scaffold: bool=True,
           promoter: str='7SK', scaffold: str='SpCas9 mod_7',
           spacer: str='Spacer_sequence', spacer_t5: str=None, spacer_t3: str=None,
           spacer_b5: str=None, spacer_b3: str=None, order_scaffold: bool=False, 
           file:str=None) -> pd.DataFrame:
    ''' 
    ngRNAs(): design GG cloning oligonucleotides for prime editing ngRNAs
    
    Parameters:
    df (dataframe | str): Dataframe with spacers (or file path)
    id (str): id column
    tG (bool, optional): add 5' G to spacer if needed (Default: True)
    order (bool, optional): order format (Default: True)
    vector_scaffold (bool, optional): the vector already carries the scaffold, so the spacer is cloned like an sgRNA (Default: True)
        True: the scaffold contributes '' spacer_t3 & AAAC spacer_b5, and there is no scaffold to order
        False: the scaffold is a separate GG fragment, so it contributes its own ends & can be ordered
    promoter (str, optional): promoter upstream of the spacer (Default: 7SK)
        hU6: CACC spacer_t5, '' spacer_b3
        7SK: TCTC spacer_t5, '' spacer_b3
    scaffold (str, optional): scaffold downstream of the spacer (Default: SpCas9 mod_7)
        vector_scaffold=True, both options: '' spacer_t3, AAAC spacer_b5
        vector_scaffold=False, SpCas9 Flip and Extend: GTTTAAGAGC spacer_t3, '' spacer_b5
        vector_scaffold=False, SpCas9 mod_7: GTTTCAGAGC spacer_t3, '' spacer_b5
    spacer (str, optional): ngRNA spacer column name
    spacer_t5 (str, optional): top oligonucleotide 5' overhang (Default: from promoter)
    spacer_t3 (str, optional): top oligonucleotide 3' overhang (Default: from scaffold)
    spacer_b5 (str, optional): bottom oligonucleotide 5' overhang (Default: from scaffold)
    spacer_b3 (str, optional): bottom oligonucleotide 3' overhang (Default: from promoter)
    order_scaffold (bool, optional): order top and bottom oligonucleotide for scaffold sequence (Default: False)
        ignored when vector_scaffold=True, since the vector carries the scaffold
    file (str, optional): output file path
    
    Assumptions:
    1. ngRNA scaffold: SpCas9 Flip and Extend (GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC)
       SpCas9 Flip and Extend + com-modified [required for VLPs] (GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC),
       or SpCas9 mod_7 (GTTTCAGAGCTAGAAATAGCAAGTTGAAATAAGACTAGTTCGTTATCACGCCGAAAGGCGGGCACCGAGTCGGTGC)
    
    Dependencies: pandas, top_bot(), ord_form(), & get_overhangs()
    '''
    if type(df)==str: # Get ngRNAs dataframe from file path if needed
        df = io.get(pt=df)
    
    # Fill in the overhangs that were not specified from the chosen promoter & scaffold
    if vector_scaffold==True: parts = [(promoter,NGRNA_VECTOR_PROMOTER_OVERHANGS),(scaffold,NGRNA_VECTOR_SCAFFOLD_OVERHANGS)] # The vector carries the scaffold, so the spacer is cloned like an sgRNA
    else: parts = [(promoter,PEGRNA_PROMOTER_OVERHANGS),(scaffold,PEGRNA_SCAFFOLD_OVERHANGS)] # The scaffold is a separate GG fragment
    overhangs = get_overhangs(parts=parts,
                              spacer_t5=spacer_t5, spacer_t3=spacer_t3, spacer_b5=spacer_b5, spacer_b3=spacer_b3)
    
    if (order_scaffold==True)&(vector_scaffold==True): # The vector carries the scaffold, so there is nothing to order
        print('Warning: Did not order the scaffold oligonucleotides!\nThe vector carries the scaffold when vector_scaffold=True.')
        order_scaffold = False
    
    df=tb(df=df,id=id,seq=spacer,t5=overhangs['spacer_t5'],t3=overhangs['spacer_t3'],b5=overhangs['spacer_b5'],b3=overhangs['spacer_b3'],tG=tG,pre='ns_') # Make top and bottom oligos for spacer inserts
    if order==True: # Sigma order format (or original dataframe with top and bottom oligos)
        if order_scaffold==True: # Include the scaffold sequence in the order
            df = pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_top',pre='ns_'),
                            ord_form(df=df,id=id,seq=spacer,suf='_bot',pre='ns_'),
                            get_scaffold_oligos(scaffold=scaffold, pre='ngRNA')
                            ]).reset_index(drop=True)
        else: # Do not include the scaffold sequence in the order
            df = pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_top',pre='ns_'),
                            ord_form(df=df,id=id,seq=spacer,suf='_bot',pre='ns_')
                            ]).reset_index(drop=True)

    # Save & return dataframe
    io.save(obj=df, file=file)  
    return df

# Library GG cloning
def epegRNA_pool(
    df: pd.DataFrame | str, barcode: str=None, 
    tG:bool=True,
    make_extension:bool=True,
    UMI_df: pd.DataFrame | str=None,
    PCR_df: pd.DataFrame| str=None,
    RE_type_IIS_df: pd.DataFrame| str=None,
    UMI_i: int=0, enzymes: list=['Esp3I'], barcode_i:int=1, 
    fwd_barcode_t5:str='Forward Barcode', rev_barcode_t3:str='Reverse Barcode',
    fwd_RE_t5:str='Forward RE Name', rev_RE_t3:str='Reverse RE Name',
    fwd_homology_t5:str='Forward Homology Name', rev_homology_t3:str='Reverse Homology Name',
    fwd_RE_t5_val:str='Esp3I', rev_RE_t3_val:str='Esp3I',
    fwd_homology_t5_val:str='7SK', rev_homology_t3_val:str='tevo2',
    epegRNA_spacer:str='Spacer_sequence', epegRNA_scaffold:str='Scaffold_sequence',
    epegRNA_extension:str='Extension_sequence', epegRNA_RTT:str='RTT_sequence',
    epegRNA_PBS:str='PBS_sequence', epegRNA_linker:str='Linker_sequence',
    file:str=None, return_df:bool=True) -> pd.DataFrame:
    ''' 
    epegRNA_pool(): design GG cloning oligonucleotides for pooled prime editing epegRNAs
    
    Parameters:
    df (dataframe | str): Dataframe with sequence information for epegRNAs (or file path)
    barcode (str, optional): subpool barcode column name in df (Default: None); disregarded when barcode_i = 0
    tG (bool, optional): add 5' G(s) to spacer if needed (Default: True); leading G yields optimal expression and G(s) may be added to avoid creating an Esp3I site at promoter-spacer junction
    make_extension (bool, optional): concatenate RTT, PBS, and linker to make extension sequence (Default: True)
    UMI_df (dataframe | str, optional): Dataframe with UMI sequences (or file path)
    PCR_df (dataframe | str, optional): Dataframe with PCR primer and subpool barcode information (or file path)
    RE_type_IIS_df (dataframe | str, optional): Dataframe with Type IIS restriction enzyme information (or file path)
    UMI_i (int, optional): UMI start index (Default: 0)
    enzymes (list, optional): list of Type IIS restriction enzymes to check for (Default: ['Esp3I'])
    barcode_i (int, optional): subpool barcode start index (Default: 1); 0 omits the forward and reverse subpool barcodes (barcode is disregarded) and uses 20 nt UMIs instead of 15 nt UMIs
    fwd_barcode_t5 (str, optional): forward barcode column name (Default: 'Forward Barcode')
    rev_barcode_t3 (str, optional): reverse barcode column name (Default: 'Reverse Barcode')
    fwd_RE_t5 (str, optional): forward restriction enzyme column name (Default: 'Forward RE Name')
    rev_RE_t3 (str, optional): reverse restriction enzyme column name (Default: 'Reverse RE Name')
    fwd_homology_t5 (str, optional): forward homology column name (Default: 'Forward Homology Name')
    rev_homology_t3 (str, optional): reverse homology column name (Default: 'Reverse Homology Name')
    fwd_RE_t5_val (str, optional): forward restriction enzyme column value (Default: 'Esp3I')
    rev_RE_t3_val (str, optional): reverse restriction enzyme column value (Default: 'Esp3I')
    fwd_homology_t5_val (str, optional): forward homology column value (Default: '7SK')
    rev_homology_t3_val (str, optional): reverse homology column value (Default: 'tevopreQ1')
    epegRNA_spacer (str, optional): epegRNA spacer column name (Default: 'Spacer_sequence')
    epegRNA_scaffold (str, optional): epegRNA scaffold sequence column name (Default: 'Scaffold_sequence')
    epegRNA_extension (str, optional): epegRNA extension name (Default: 'Extension_sequence')
    epegRNA_RTT (str, optional): epegRNA reverse transcripase template column name (Default: 'RTT_sequence')
    epegRNA_PBS (str, optional): epegRNA primer binding site column name (Default: 'PBS_sequence')
    epegRNA_linker (str, optional): epegRNA linker column name (Default: 'Linker_sequence')
    file (str, optional): output file path (Default: None)
    return_df (bool, optional): return dataframe (Default: True)

    barcode_i > 0 Assumptions:
    1. Oligo Template: FWD UMI - FWD Barcode - FWD RE - FWD Homology - epegRNA_spacer - epegRNA_scaffold - epegRNA_extension - RE Homology - REV RE - REV Barcode - REV UMI
    2. epegRNA motif: tevoPreQ1 (CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA) or tevo2 (CGTCGGTTCTACCTAGTTACGACGTTAAACCAACTAGGA)
    3. epegRNA scaffold: SpCas9 Flip and Extend (GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC)
       SpCas9 Flip and Extend + com-modified [required for VLPs] (GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC),
       or SpCas9 mod_7 (GTTTCAGAGCTAGAAATAGCAAGTTGAAATAAGACTAGTTCGTTATCACGCCGAAAGGCGGGCACCGAGTCGGTGC)
    
    barcode_i = 0 Assumptions:
    1. Oligo Template: FWD UMI - FWD RE - FWD Homology - epegRNA_spacer - epegRNA_scaffold - epegRNA_extension - RE Homology - REV RE - REV UMI
    2. epegRNA motif: tevoPreQ1 (CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA) or tevo2 (CGTCGGTTCTACCTAGTTACGACGTTAAACCAACTAGGA)
    3. epegRNA scaffold: SpCas9 Flip and Extend (GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC)
       SpCas9 Flip and Extend + com-modified [required for VLPs] (GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC),
       or SpCas9 mod_7 (GTTTCAGAGCTAGAAATAGCAAGTTGAAATAAGACTAGTTCGTTATCACGCCGAAAGGCGGGCACCGAGTCGGTGC)
    
    Dependencies: pandas,
    '''
    # Get dataframes from file paths if needed
    if isinstance(df, str):
        df = io.get(pt=df)
    else:
        df = df.copy()

    if isinstance(UMI_df, str):
        UMI_df = io.get(pt=UMI_df)
    elif isinstance(UMI_df, pd.DataFrame):
        UMI_df = UMI_df.copy()

    if isinstance(PCR_df, str):
        PCR_df = io.get(pt=PCR_df)
    elif isinstance(PCR_df, pd.DataFrame):
        PCR_df = PCR_df.copy()

    if isinstance(RE_type_IIS_df, str):
        RE_type_IIS_df = io.get(pt=RE_type_IIS_df)
    elif isinstance(RE_type_IIS_df, pd.DataFrame):
        RE_type_IIS_df = RE_type_IIS_df.copy()
    
    # Load UMI and PCR dataframes from resources if not provided
    subpool = barcode_i != 0 # barcode_i = 0 => no subpool barcodes; 20 nt UMIs take their place
    if subpool == False: barcode = None # Disregard subpool barcode column
    save_time = False
    if UMI_df is None and PCR_df is None and RE_type_IIS_df is None and enzymes == ['Esp3I']: # Default UMI, PCR, Enzyme dataframes, and enzymes (save time)
        if subpool: UMI_df = load_resource_csv(filename='UMI_15_hamming_4_yield_18687_Esp3I_0.csv')
        else: UMI_df = load_resource_csv(filename='UMI_10_hamming_5_yield_73_2_combinations_interleaved_Esp3I_0.csv')
        PCR_df = load_resource_csv(filename='edms_pcr.csv')
        RE_type_IIS_df = load_resource_csv(filename='RE_type_IIS.csv')
        save_time = True
    
    else: # Not default (checking for RE sites)
        if UMI_df is None: # Get UMI dataframe from resources if not provided
            if subpool: UMI_df = load_resource_csv(filename='UMI_15_hamming_4_yield_19356.csv')
            else: UMI_df = load_resource_csv(filename='UMI_10_hamming_5_yield_73_2_combinations_interleaved_Esp3I_0.csv')
        
        if PCR_df is None: # Get PCR dataframe from resources if not provided
            PCR_df = load_resource_csv(filename='edms_pcr.csv')
        else:
            PCR_df = PCR_df.copy() # Make copy of provided PCR dataframe to avoid modifying original
            if not all(col in PCR_df.columns for col in [fwd_barcode_t5, rev_barcode_t3, fwd_RE_t5, rev_RE_t3, fwd_homology_t5, rev_homology_t3, fwd_RE_t5.replace("Name", "Sequence"), rev_RE_t3.replace("Name", "Sequence"), fwd_homology_t5.replace("Name", "Sequence"), rev_homology_t3.replace("Name", "Sequence")]):
                raise ValueError(f'PCR dataframe must include columns: {fwd_barcode_t5}, {rev_barcode_t3}, {fwd_RE_t5}, {rev_RE_t3}, {fwd_homology_t5}, and {rev_homology_t3} to filter for specific RE and homology sequences!\nAlso, need columns with RE and homology sequences: {fwd_RE_t5.replace("Name", "Sequence")}, {rev_RE_t3.replace("Name", "Sequence")}, {fwd_homology_t5.replace("Name", "Sequence")}, and {rev_homology_t3.replace("Name", "Sequence")}!')
        
        if RE_type_IIS_df is None: # Get from resources if not provided
            RE_type_IIS_df = load_resource_csv(filename='RE_type_IIS.csv')

    # Make extension by concatenating RTT, PBS, and linker
    if make_extension==True: df[epegRNA_extension] = df[epegRNA_RTT]+df[epegRNA_PBS]+df[epegRNA_linker]
    else: print(f'Warning: Did not make extension sequence!\nMake sure "{epegRNA_extension}" column includes RTT+PBS+linker for epegRNAs.')

    # Print PCR_df options
    print(f'PCR_df options:\n{PCR_df[[fwd_RE_t5, rev_RE_t3, fwd_homology_t5, rev_homology_t3]].drop_duplicates().reset_index(drop=True)}\n')

    # Filter PCR_df for specific RE and homology sequences
    PCR_df = PCR_df[
        (PCR_df[fwd_RE_t5] == fwd_RE_t5_val) &
        (PCR_df[rev_RE_t3] == rev_RE_t3_val) &
        (PCR_df[fwd_homology_t5] == fwd_homology_t5_val) &
        (PCR_df[rev_homology_t3] == rev_homology_t3_val)
    ]

    # Assign subpool barcodes
    if subpool:
        barcodes = sorted(df[barcode].unique())
        PCR_barcodes = PCR_df.iloc[barcode_i:barcode_i+len(barcodes)] # Get corresponding barcodes from PCR dataframe starting at barcode_i
        PCR_barcodes[barcode] = barcodes
        df = pd.merge(left=df,right=PCR_barcodes,on=barcode)
        print(f'{len(PCR_barcodes)} Barcodes: index = [{barcode_i}:{barcode_i+PCR_barcodes.shape[0]}]')    
    else: # barcode_i = 0: no subpool barcodes, keep the RE and homology sequences from the first PCR dataframe row
        PCR_barcodes = PCR_df.iloc[:1].copy()
        PCR_barcodes[fwd_barcode_t5] = ''
        PCR_barcodes[rev_barcode_t3] = ''
        for col,val in PCR_barcodes.iloc[0].items(): df[col] = val
        print('0 Barcodes: no subpool barcodes (20 nt UMIs instead)')

    # Check for RE sites in UMI sequences + subpool barcodes (discard UMIs with RE sites)
    if save_time == False: # Not default UMI, PCR, Enzyme dataframes, and enzymes (do not save time)
        fwd_UMI_ls = []
        fwd_barcode_ls = []
        fwd_UMI_barcode_ls = []

        rev_UMI_ls = []
        rev_barcode_ls = []
        rev_barcode_UMI_ls = []

        for umi in UMI_df['UMI_sequence']: # Iterate through UMIs
            for fwd_barcode in PCR_barcodes[fwd_barcode_t5]: # Iterate through FWD barcodes
                fwd_UMI_ls.append(umi)
                fwd_barcode_ls.append(fwd_barcode)
                fwd_UMI_barcode_ls.append(umi+fwd_barcode)
            
            for rev_barcode in PCR_barcodes[rev_barcode_t3]: # Iterate through REV barcodes
                rev_UMI_ls.append(umi)
                rev_barcode_ls.append(rev_barcode)
                rev_barcode_UMI_ls.append(rev_barcode+umi)
        
        UMI_barcode_df = pd.DataFrame({'Forward UMI':fwd_UMI_ls,
                                    'Forward Barcode':fwd_barcode_ls,
                                    'Forward UMI+Barcode':fwd_UMI_barcode_ls,
                                    'Reverse UMI':rev_UMI_ls,
                                    'Reverse Barcode':rev_barcode_ls,
                                    'Reverse Barcode+UMI':rev_barcode_UMI_ls})

        for (enzyme,recognition,recognition_rc) in t.zip_cols(df=RE_type_IIS_df[RE_type_IIS_df['Name'].isin(enzymes)], # Just Esp3I by default
                                                            cols=['Name','Recognition','Recognition_rc']):
            # Check Forward UMI+Barcode for RE sites
            UMI_barcode_df[f'{enzyme} Forward UMI+Barcode'] = [len(t.find_all(oligo,recognition))+len(t.find_all(oligo,recognition_rc)) for oligo in UMI_barcode_df['Forward UMI+Barcode']] # Iterate through oligonucleotides
            
            # Check Forward UMI+Barcode for RE sites
            UMI_barcode_df[f'{enzyme} Reverse Barcode+UMI'] = [len(t.find_all(oligo,recognition))+len(t.find_all(oligo,recognition_rc)) for oligo in UMI_barcode_df['Reverse Barcode+UMI']] # Iterate through oligonucleotides

            # Count number of RE sites (should be 0)
            RE_ls = []
            for umi in UMI_df['UMI_sequence']:
                RE_ls.append(sum(UMI_barcode_df[(UMI_barcode_df['Forward UMI']==umi)][f'{enzyme} Forward UMI+Barcode']) + sum(UMI_barcode_df[(UMI_barcode_df['Reverse UMI']==umi)][f'{enzyme} Reverse Barcode+UMI']))
            UMI_df[enzyme] = RE_ls
        
        # Keep UMIs with no RE sites
        UMI_df = UMI_df[UMI_df[enzymes].sum(axis=1)==0]

    # Assign UMI to each oligo
    UMI_sequences = UMI_df['UMI_sequence'].tolist()
    if len(UMI_sequences[UMI_i:])<df.shape[0]*2: # Check if there are enough UMIs for the number of oligos
        ValueError(f'{len(UMI_sequences)} UMIs is less than 2x{df.shape[0]} oligonucleotides!')
    UMI_sequences = UMI_sequences[UMI_i:UMI_i+df.shape[0]*2] # Get UMIs starting at UMI_i
    print(f'{len(UMI_sequences)} UMIs: index = [{UMI_i}:{UMI_i+df.shape[0]*2}]')
    df['Forward UMI'] = UMI_sequences[:df.shape[0]] # Assign forward UMIs
    df['Reverse UMI'] = UMI_sequences[df.shape[0]:] # Assign reverse UMIs

    # Make oligo & determine length
    df[f'{epegRNA_spacer}_nt1']=[s[0] for s in df[epegRNA_spacer]]
    oligos = []
    for i,(epegRNA_spacer_nt1) in enumerate(df[f'{epegRNA_spacer}_nt1']):
        if tG==True: # Add 5'G(s) to spacer if needed for optimal expression or to avoid creating an Esp3I site at promoter-spacer junction
            if (epegRNA_spacer_nt1!='G'): # Spacer does not start with G, add 5'G(s)
                promoter_spacer_junction = df.iloc[i][fwd_homology_t5.replace("Name", "Sequence")]+'G'+df.iloc[i][epegRNA_spacer] # Check if adding 5'G creates an RE site at promoter-spacer junction
                if any(recognition in promoter_spacer_junction for recognition in RE_type_IIS_df[RE_type_IIS_df['Name'].isin(enzymes)]['Recognition'].tolist()) or any(recognition_rc in promoter_spacer_junction for recognition_rc in RE_type_IIS_df[RE_type_IIS_df['Name'].isin(enzymes)]['Recognition_rc'].tolist()): # It does, add 2 Gs
                    add_G = 'GG'
                else: # It does not, add 1 G
                    add_G = 'G'
            else: # Spacer already starts with G; still need to check if it creates an RE site at promoter-spacer junction
                promoter_spacer_junction = df.iloc[i][fwd_homology_t5.replace("Name", "Sequence")]+df.iloc[i][epegRNA_spacer] # Check if promoter-spacer junction creates an RE site
                if any(recognition in promoter_spacer_junction for recognition in RE_type_IIS_df[RE_type_IIS_df['Name'].isin(enzymes)]['Recognition'].tolist()) or any(recognition_rc in promoter_spacer_junction for recognition_rc in RE_type_IIS_df[RE_type_IIS_df['Name'].isin(enzymes)]['Recognition_rc'].tolist()): # It does, add 1 G
                    add_G = 'G'
                else: # It does not, do not add G
                    add_G = ''
        else: # Do not append 5'G to spacer
            add_G = ''
        # Make final oligonucleotide sequence based on template: FWD UMI - FWD Barcode - FWD RE - FWD Homology - epegRNA_spacer - epegRNA_scaffold - epegRNA_extension - RE Homology - REV RE - REV Barcode - REV UMI
        oligo = df.iloc[i]['Forward UMI']+df.iloc[i][fwd_barcode_t5]+df.iloc[i][fwd_RE_t5.replace("Name", "Sequence")]+df.iloc[i][fwd_homology_t5.replace("Name", "Sequence")]+add_G+ \
                        df.iloc[i][epegRNA_spacer]+df.iloc[i][epegRNA_scaffold]+df.iloc[i][epegRNA_extension]+ \
                        df.iloc[i][rev_homology_t3.replace("Name", "Sequence")]+df.iloc[i][rev_RE_t3.replace("Name", "Sequence")]+df.iloc[i][rev_barcode_t3]+df.iloc[i]['Reverse UMI']
        oligos.append(oligo.upper()) # Make sure oligos are uppercase
            
    df['Oligonucleotide'] = oligos
    df['Oligonucleotide_length']=[len(oligo) for oligo in df['Oligonucleotide']]

    # Check for 2 recognition sites per enzyme
    for (enzyme,recognition,recognition_rc) in t.zip_cols(df=RE_type_IIS_df[RE_type_IIS_df['Name'].isin(enzymes)], # Just Esp3I by default
                                                          cols=['Name','Recognition','Recognition_rc']):
        # Check forward direction for recognition sites
        enzyme_sites_fwd = [t.find_all(oligo,recognition) for oligo in oligos] # Iterate through oligonucleotides

        # Check reverse direction for recognition sites
        enzyme_sites_rc = [t.find_all(oligo,recognition_rc) for oligo in oligos] # Iterate through oligonucleotides
        
        # Sum recognition sites in both directions
        enzyme_sites = [len(enzyme_site_fwd)+len(enzyme_site_rc) for (enzyme_site_fwd,enzyme_site_rc) in zip(enzyme_sites_fwd,enzyme_sites_rc)]
        
        df[enzyme] = enzyme_sites
        df[f'{enzyme}_fwd_i'] = enzyme_sites_fwd
        df[f'{enzyme}_rc_i'] = enzyme_sites_rc
    
    # Save & return dataframe
    io.save(obj=df, file=file)
    if return_df:
        return df

def dms_pool(
    df: pd.DataFrame | str, barcode: str,
    fwd_homology_t5_val: str | list, rev_homology_t3_val: str | list,
    barcode_val: list | None = None,
    UMI_df: pd.DataFrame | str=None,
    PCR_df: pd.DataFrame| str=None,
    RE_type_IIS_df: pd.DataFrame| str=None,
    UMI_i: int=0, enzymes: list=['Esp3I'], barcode_i:int=0, 
    fwd_barcode_t5:str='Forward Barcode', rev_barcode_t3:str='Reverse Barcode',
    fwd_RE_t5:str='Forward RE Name', rev_RE_t3:str='Reverse RE Name',
    fwd_homology_t5:str='Forward Homology Name', rev_homology_t3:str='Reverse Homology Name',
    fwd_RE_t5_val:str='Esp3I', rev_RE_t3_val:str='Esp3I',
    template:str='Edit_sequence_with_silent_mutations',
    file:str=None, return_df:bool=True) -> pd.DataFrame:
    ''' 
    dms_pool(): design GG cloning oligonucleotides for deep mutational scanning (DMS) libraries
    
    Parameters:
    df (dataframe | str): Dataframe with sequence information for epegRNAs (or file path)
    barcode (str): subpool barcode column name in df
    fwd_homology_t5_val (str | list): forward homology column value in PCR_df (Examples: 'FOXA1-S165', 'FOXA1-W199', 'FOXA1-R233', or ['FOXA1-S165', 'FOXA1-W199', 'FOXA1-R233'])
    rev_homology_t3_val (str | list): reverse homology column value in PCR_df (Examples: 'FOXA1-P205', 'FOXA1-G239', 'FOXA1-P272', or ['FOXA1-P205', 'FOXA1-G239', 'FOXA1-P272'])
    barcode_val (list | None): subpool barcode column values in df[barcode] (Examples: ['1', '2', '3'] or None)
    UMI_df (dataframe | str, optional): Dataframe with UMI sequences (or file path)
    PCR_df (dataframe | str, optional): Dataframe with PCR primer and subpool barcode information (or file path)
    RE_type_IIS_df (dataframe | str, optional): Dataframe with Type IIS restriction enzyme information (or file path)
    UMI_i (int, optional): UMI start index (Default: 0)
    enzymes (list, optional): list of Type IIS restriction enzymes to check for (Default: ['Esp3I'])
    barcode_i (int, optional): subpool barcode start index (Default: 0)
    fwd_barcode_t5 (str, optional): forward barcode column name (Default: 'Forward Barcode')
    rev_barcode_t3 (str, optional): reverse barcode column name (Default: 'Reverse Barcode')
    fwd_RE_t5 (str, optional): forward restriction enzyme column name (Default: 'Forward RE Name')
    rev_RE_t3 (str, optional): reverse restriction enzyme column name (Default: 'Reverse RE Name')
    fwd_homology_t5 (str, optional): forward homology column name (Default: 'Forward Homology Name')
    rev_homology_t3 (str, optional): reverse homology column name (Default: 'Reverse Homology Name')
    fwd_RE_t5_val (str, optional): forward restriction enzyme column value (Default: 'Esp3I')
    rev_RE_t3_val (str, optional): reverse restriction enzyme column value (Default: 'Esp3I')
    template (str, optional): oligonucleotide template column name in df (Default: 'Edit_sequence_with_silent_mutations')
    file (str, optional): output file path (Default: None)
    return_df (bool, optional): return dataframe (Default: True)

    Assumptions:
    1. Oligo Template: FWD UMI - FWD Barcode - FWD RE - FWD Homology - Trimmed Template - RE Homology - REV RE - REV Barcode - REV UMI
    
    Dependencies: pandas
    '''
    # Get dataframes from file paths if needed
    if isinstance(df, str):
        df = io.get(pt=df)
    else:
        df = df.copy().reset_index(drop=True)

    if isinstance(UMI_df, str):
        UMI_df = io.get(pt=UMI_df)
    elif isinstance(UMI_df, pd.DataFrame):
        UMI_df = UMI_df.copy()

    if isinstance(PCR_df, str):
        PCR_df = io.get(pt=PCR_df)
    elif isinstance(PCR_df, pd.DataFrame):
        PCR_df = PCR_df.copy()

    if isinstance(RE_type_IIS_df, str):
        RE_type_IIS_df = io.get(pt=RE_type_IIS_df)
    elif isinstance(RE_type_IIS_df, pd.DataFrame):
        RE_type_IIS_df = RE_type_IIS_df.copy()
    
    # Load UMI and PCR dataframes from resources if not provided
    save_time = False
    if UMI_df is None and PCR_df is None and RE_type_IIS_df is None and enzymes == ['Esp3I']: # Default UMI, PCR, Enzyme dataframes, and enzymes (save time)
        UMI_df = load_resource_csv(filename='UMI_15_hamming_4_yield_18687_Esp3I_0.csv')
        PCR_df = load_resource_csv(filename='dms_pcr.csv')
        RE_type_IIS_df = load_resource_csv(filename='RE_type_IIS.csv')
        save_time = True
    
    else: # Not default (checking for RE sites)
        if UMI_df is None: # Get UMI dataframe from resources if not provided
            UMI_df = load_resource_csv(filename='UMI_15_hamming_4_yield_19356.csv')
        
        if PCR_df is None: # Get PCR dataframe from resources if not provided
            PCR_df = load_resource_csv(filename='dms_pcr.csv')
        else:
            PCR_df = PCR_df.copy() # Make copy of provided PCR dataframe to avoid modifying original
            if not all(col in PCR_df.columns for col in [fwd_barcode_t5, rev_barcode_t3, fwd_RE_t5, rev_RE_t3, fwd_homology_t5, rev_homology_t3, fwd_RE_t5.replace("Name", "Sequence"), rev_RE_t3.replace("Name", "Sequence"), fwd_homology_t5.replace("Name", "Sequence"), rev_homology_t3.replace("Name", "Sequence")]):
                raise ValueError(f'PCR dataframe must include columns: {fwd_barcode_t5}, {rev_barcode_t3}, {fwd_RE_t5}, {rev_RE_t3}, {fwd_homology_t5}, and {rev_homology_t3} to filter for specific RE and homology sequences!\nAlso, need columns with RE and homology sequences: {fwd_RE_t5.replace("Name", "Sequence")}, {rev_RE_t3.replace("Name", "Sequence")}, {fwd_homology_t5.replace("Name", "Sequence")}, and {rev_homology_t3.replace("Name", "Sequence")}!')
        
        if RE_type_IIS_df is None: # Get from resources if not provided
            RE_type_IIS_df = load_resource_csv(filename='RE_type_IIS.csv')

    # Print PCR_df options
    print(f'PCR_df options:\n{PCR_df[[fwd_RE_t5, rev_RE_t3, fwd_homology_t5, rev_homology_t3]].drop_duplicates().reset_index(drop=True)}\n')

    # Convert length 1 lists to strings for fwd_homology_t5_val and rev_homology_t3_val
    if isinstance(fwd_homology_t5_val, list) and len(fwd_homology_t5_val) == 1:
        fwd_homology_t5_val = fwd_homology_t5_val[0]
    if isinstance(rev_homology_t3_val, list) and len(rev_homology_t3_val) == 1:
        rev_homology_t3_val = rev_homology_t3_val[0]

    # Check that fwd_homology_t5_val and rev_homology_t3_val are both strings or both lists of the same length; barcode_val must be None (string) or list of the same length.
    if isinstance(fwd_homology_t5_val, str):
        if isinstance(rev_homology_t3_val, list):
            raise ValueError(f'If fwd_homology_t5_val is a string, rev_homology_t3_val and barcode_val must also be a string (not a list)!')
        elif barcode_val is not None:
            raise ValueError(f'If fwd_homology_t5_val is a string, barcode_val must be None!')
    elif isinstance(fwd_homology_t5_val, list):
        if isinstance(rev_homology_t3_val, str) or barcode_val is None:
            raise ValueError(f'If fwd_homology_t5_val is a list, rev_homology_t3_val and barcode_val must also be a list (not None)!')
        elif len(fwd_homology_t5_val) != len(rev_homology_t3_val) or len(fwd_homology_t5_val) != len(barcode_val):
            raise ValueError(f'fwd_homology_t5_val, rev_homology_t3_val, and barcode_val must be the same length if all are lists!')

    if barcode_val is not None: # Assign barcode_val to df[barcode] if provided
        if any(barcode_val[i] != df[barcode].unique()[i] for i in range(len(barcode_val))): # Check that barcode_val matches df[barcode] unique values
            raise ValueError(f'Barcode val does not exactly match the unique values in df[barcode]:\n\tBarcode_val: {barcode_val}\n\tUnique barcodes in df: {df[barcode].unique().tolist()}')


    # Filter PCR_df for specific RE and homology sequences
    if isinstance(fwd_homology_t5_val, str) & isinstance(rev_homology_t3_val, str): # All are strings
        PCR_df = PCR_df[
            (PCR_df[fwd_RE_t5] == fwd_RE_t5_val) &
            (PCR_df[rev_RE_t3] == rev_RE_t3_val) &
            (PCR_df[fwd_homology_t5] == fwd_homology_t5_val) &
            (PCR_df[rev_homology_t3] == rev_homology_t3_val)
        ].reset_index(drop=True)
        barcodes = sorted(df[barcode].unique())
        PCR_barcodes = PCR_df.iloc[barcode_i:barcode_i+len(barcodes)] # Get corresponding barcodes from PCR dataframe starting at barcode_i
        PCR_barcodes[barcode] = barcodes
        
    elif isinstance(fwd_homology_t5_val, list) & isinstance(rev_homology_t3_val, list): # Both are lists
        
        PCR_barcodes = pd.DataFrame()
        for i, (barcode_v, fwd_homology, rev_homology) in enumerate(zip(barcode_val,fwd_homology_t5_val, rev_homology_t3_val)):
            PCR_df_temp = PCR_df[
                (PCR_df[fwd_RE_t5] == fwd_RE_t5_val) &
                (PCR_df[rev_RE_t3] == rev_RE_t3_val) &
                (PCR_df[fwd_homology_t5] == fwd_homology) &
                (PCR_df[rev_homology_t3] == rev_homology)
            ].reset_index(drop=True).iloc[barcode_i+i:barcode_i+i+1] # Get corresponding barcodes from PCR dataframe starting at barcode_i
            PCR_df_temp[barcode] = barcode_v
            PCR_barcodes = pd.concat([PCR_barcodes, PCR_df_temp], ignore_index=True)

    # Assign subpool barcodes
    df = pd.merge(left=df,right=PCR_barcodes,on=barcode)
    print(f'{len(PCR_barcodes)} Barcodes: index = [{barcode_i}:{barcode_i+PCR_barcodes.shape[0]}]')    

    # Check for RE sites in UMI sequences + subpool barcodes (discard UMIs with RE sites)
    if save_time == False: # Not default UMI, PCR, Enzyme dataframes, and enzymes (do not save time)
        fwd_UMI_ls = []
        fwd_barcode_ls = []
        fwd_UMI_barcode_ls = []

        rev_UMI_ls = []
        rev_barcode_ls = []
        rev_barcode_UMI_ls = []

        for umi in UMI_df['UMI_sequence']: # Iterate through UMIs
            for fwd_barcode in PCR_barcodes[fwd_barcode_t5]: # Iterate through FWD barcodes
                fwd_UMI_ls.append(umi)
                fwd_barcode_ls.append(fwd_barcode)
                fwd_UMI_barcode_ls.append(umi+fwd_barcode)
            
            for rev_barcode in PCR_barcodes[rev_barcode_t3]: # Iterate through REV barcodes
                rev_UMI_ls.append(umi)
                rev_barcode_ls.append(rev_barcode)
                rev_barcode_UMI_ls.append(rev_barcode+umi)
        
        UMI_barcode_df = pd.DataFrame({'Forward UMI':fwd_UMI_ls,
                                    'Forward Barcode':fwd_barcode_ls,
                                    'Forward UMI+Barcode':fwd_UMI_barcode_ls,
                                    'Reverse UMI':rev_UMI_ls,
                                    'Reverse Barcode':rev_barcode_ls,
                                    'Reverse Barcode+UMI':rev_barcode_UMI_ls})

        for (enzyme,recognition,recognition_rc) in t.zip_cols(df=RE_type_IIS_df[RE_type_IIS_df['Name'].isin(enzymes)], # Just Esp3I by default
                                                            cols=['Name','Recognition','Recognition_rc']):
            # Check Forward UMI+Barcode for RE sites
            UMI_barcode_df[f'{enzyme} Forward UMI+Barcode'] = [len(t.find_all(oligo,recognition))+len(t.find_all(oligo,recognition_rc)) for oligo in UMI_barcode_df['Forward UMI+Barcode']] # Iterate through oligonucleotides
            
            # Check Forward UMI+Barcode for RE sites
            UMI_barcode_df[f'{enzyme} Reverse Barcode+UMI'] = [len(t.find_all(oligo,recognition))+len(t.find_all(oligo,recognition_rc)) for oligo in UMI_barcode_df['Reverse Barcode+UMI']] # Iterate through oligonucleotides

            # Count number of RE sites (should be 0)
            RE_ls = []
            for umi in UMI_df['UMI_sequence']:
                RE_ls.append(sum(UMI_barcode_df[(UMI_barcode_df['Forward UMI']==umi)][f'{enzyme} Forward UMI+Barcode']) + sum(UMI_barcode_df[(UMI_barcode_df['Reverse UMI']==umi)][f'{enzyme} Reverse Barcode+UMI']))
            UMI_df[enzyme] = RE_ls
        
        # Keep UMIs with no RE sites
        UMI_df = UMI_df[UMI_df[enzymes].sum(axis=1)==0]

    # Assign UMI to each oligo
    UMI_sequences = UMI_df['UMI_sequence'].tolist()
    if len(UMI_sequences[UMI_i:])<df.shape[0]*2: # Check if there are enough UMIs for the number of oligos
        ValueError(f'{len(UMI_sequences)} UMIs is less than 2x{df.shape[0]} oligonucleotides!')
    UMI_sequences = UMI_sequences[UMI_i:UMI_i+df.shape[0]*2] # Get UMIs starting at UMI_i
    print(f'{len(UMI_sequences)} UMIs: index = [{UMI_i}:{UMI_i+df.shape[0]*2}]')
    df['Forward UMI'] = UMI_sequences[:df.shape[0]] # Assign forward UMIs
    df['Reverse UMI'] = UMI_sequences[df.shape[0]:] # Assign reverse UMIs

    # Make oligo & determine length
    oligos = []
    trimmed_templates = []
    for i, (temp_seq, fwd_homology_t5_trim, rev_homology_t3_trim) in enumerate(t.zip_cols(
        df=df, cols=[template, fwd_homology_t5.replace("Name", "Sequence Trim"), rev_homology_t3.replace("Name", "Sequence Trim")])):

        # Trim flanks if provided
        if temp_seq.find(fwd_homology_t5_trim)==-1 or temp_seq.rfind(rev_homology_t3_trim)==-1:
            raise(ValueError(f"Flank5 or Flank3 sequences were not found in edit sequence for DMS row {i}.\nPlease check for flank5 ({fwd_homology_t5_trim}) and flank3 ({rev_homology_t3_trim}) sequences in the template sequence:\n{temp_seq}"))

        # Trim flanks from template sequences
        trimmed_temp_seq = temp_seq[temp_seq.find(fwd_homology_t5_trim)+len(fwd_homology_t5_trim) : temp_seq.rfind(rev_homology_t3_trim)]
        trimmed_templates.append(trimmed_temp_seq)

        # Make final oligonucleotide sequence based on template: FWD UMI - FWD Barcode - FWD RE - FWD Homology - Trimmed Template - RE Homology - REV RE - REV Barcode - REV UMI
        oligo = df.iloc[i]['Forward UMI']+df.iloc[i][fwd_barcode_t5]+df.iloc[i][fwd_RE_t5.replace("Name", "Sequence")]+df.iloc[i][fwd_homology_t5.replace("Name", "Sequence")]+ \
                trimmed_templates[i]+ \
                df.iloc[i][rev_homology_t3.replace("Name", "Sequence")]+df.iloc[i][rev_RE_t3.replace("Name", "Sequence")]+df.iloc[i][rev_barcode_t3]+df.iloc[i]['Reverse UMI']
        oligos.append(oligo.upper()) # Make sure oligos are uppercase

    df[f'{template}_trimmed'] = trimmed_templates
    df['Oligonucleotide'] = oligos
    df['Oligonucleotide_length']=[len(oligo) for oligo in df['Oligonucleotide']]

    # Check for 2 recognition sites per enzyme
    for (enzyme,recognition,recognition_rc) in t.zip_cols(df=RE_type_IIS_df[RE_type_IIS_df['Name'].isin(enzymes)], # Just Esp3I by default
                                                          cols=['Name','Recognition','Recognition_rc']):
        # Check forward direction for recognition sites
        enzyme_sites_fwd = [t.find_all(oligo,recognition) for oligo in oligos] # Iterate through oligonucleotides

        # Check reverse direction for recognition sites
        enzyme_sites_rc = [t.find_all(oligo,recognition_rc) for oligo in oligos] # Iterate through oligonucleotides
        
        # Sum recognition sites in both directions
        enzyme_sites = [len(enzyme_site_fwd)+len(enzyme_site_rc) for (enzyme_site_fwd,enzyme_site_rc) in zip(enzyme_sites_fwd,enzyme_sites_rc)]
        
        df[enzyme] = enzyme_sites
        df[f'{enzyme}_fwd_i'] = enzyme_sites_fwd
        df[f'{enzyme}_rc_i'] = enzyme_sites_rc
    
    # Save & return dataframe
    io.save(obj=df, file=file)
    if return_df:
        return df

# UMIs
def generate_sequences(length:int, current_sequence:str="") -> list:
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

def filter_GC(sequences: list, GC_fract: tuple) -> list:
    '''
    filter_GC(): filters sequences based on GC content
    
    Parameters:
    sequences (list): list of sequences (str)
    GC_fract (tuple): Pair of GC content boundaries written fractions (Ex: (0.4,0.6))
    '''
    return [sequence for sequence in sequences if ((len(t.find_all(sequence,'G'))+len(t.find_all(sequence,'C')))/len(sequence)>GC_fract[0])&((len(t.find_all(sequence,'G'))+len(t.find_all(sequence,'C')))/len(sequence)<GC_fract[1])]

def shuffle(ls: list) -> list:
    """
    shuffle(): randomly reorganizes a list

    Parameters:
    ls (list): the list to be shuffled.
    
    Dependencies: random
    """
    ls2 = ls[:]  # Create a copy of the list to avoid modifying the original
    random.shuffle(ls2)
    return ls2

def encode_sequences(seq_list: list) -> np.ndarray:
    """
    encode_sequences(): Convert sequences to integer arrays for fast comparison.

    Parameters:
    seq_list (list of str): List of equal-length DNA sequences (A/C/G/T).

    Dependencies: numpy
    """
    if len(seq_list) == 0:
        return np.empty((0, 0), dtype=np.uint8)

    lengths = {len(seq) for seq in seq_list}
    if len(lengths) > 1:
        raise ValueError(f"encode_sequences() requires equal-length sequences; found lengths {sorted(lengths)}.")

    # Pack every sequence into one ascii buffer, then map bases through a lookup table (much faster than a Python loop)
    ascii_arr = np.frombuffer(''.join(seq_list).encode('ascii'), dtype=np.uint8).reshape(len(seq_list), -1)
    lut = np.full(256, 255, dtype=np.uint8)
    for i, base in enumerate('ACGT'):
        lut[ord(base)] = i
    encoded = lut[ascii_arr]
    if np.any(encoded == 255):
        raise ValueError("encode_sequences() found a character outside A/C/G/T.")
    return encoded


def hamming_matrix(a_encoded: np.ndarray, b_encoded: np.ndarray) -> np.ndarray:
    """
    hamming_matrix(): pairwise Hamming distances between two encoded sequence sets.

    Parameters:
    a_encoded (ndarray): (m, length) encoded sequences from encode_sequences()
    b_encoded (ndarray): (n, length) encoded sequences from encode_sequences()

    Returns (m, n) ndarray of Hamming distances.

    Dependencies: numpy
    """
    distances = np.zeros((a_encoded.shape[0], b_encoded.shape[0]), dtype=np.int16)
    for position in range(a_encoded.shape[1]):  # accumulate per position to avoid an (m, n, length) temporary
        distances += (a_encoded[:, position][:, None] != b_encoded[None, :, position])
    return distances


def min_hamming_to(sequences: list, reference: list) -> np.ndarray:
    """
    min_hamming_to(): minimum Hamming distance from each sequence to any sequence in a reference set.

    Parameters:
    sequences (list of str): List of equal-length DNA sequences (A/C/G/T).
    reference (list of str): List of sequences of the same length to measure against.

    Returns one distance per sequence, or np.nan throughout if the reference set is empty.

    Dependencies: numpy, encode_sequences() & hamming_matrix()
    """
    if len(sequences) == 0:
        return np.empty(0)
    if not reference:
        return np.full(len(sequences), np.nan)

    encoded = encode_sequences(sequences)
    reference_encoded = encode_sequences(reference)
    if encoded.shape[1] != reference_encoded.shape[1]:
        raise ValueError(f"min_hamming_to() needs equal lengths; sequences are {encoded.shape[1]} nt and reference sequences are {reference_encoded.shape[1]} nt.")

    out = np.empty(len(sequences), dtype=np.int32)
    block = max(1, int(2e7 // len(reference)))  # cap each distance matrix at ~40 MB
    for start in range(0, len(sequences), block):
        stop = min(start + block, len(sequences))
        out[start:stop] = hamming_matrix(encoded[start:stop], reference_encoded).min(axis=1)
    return out


def read_sequences(obj, col: str = 'UMI_sequence', length: int = None) -> list:
    """
    read_sequences(): read a list of DNA sequences from a file path, dataframe, series, or iterable.

    Sequences are uppercased and stripped, blanks and duplicates are dropped, and order is preserved.

    Parameters:
    obj (str | DataFrame | Series | list | None): file path, dataframe, series or iterable of sequences;
        None returns an empty list
    col (str, optional): column name holding the sequences when obj is a file path or dataframe
        (Default: 'UMI_sequence')
    length (int, optional): raise if any sequence is not this length (Default: None, no check)

    Dependencies: pandas & io
    """
    if obj is None:
        return []
    if isinstance(obj, str):
        obj = io.get(pt=obj)

    if isinstance(obj, pd.DataFrame):
        if col not in obj.columns:
            raise ValueError(f"read_sequences() could not find column '{col}'; the file or dataframe has {list(obj.columns)}.")
        values = obj[col]
    elif isinstance(obj, pd.Series):
        values = obj
    else:
        values = pd.Series(list(obj), dtype=object)

    sequences = [str(value).strip().upper() for value in values.dropna()]
    sequences = list(dict.fromkeys(sequence for sequence in sequences if sequence))  # drop blanks & duplicates

    if length is not None:
        wrong = sorted({len(sequence) for sequence in sequences} - {length})
        if wrong:
            raise ValueError(f"read_sequences() expected sequences of length {length}; found lengths {wrong}. Hamming distances are only defined between equal-length sequences.")
    return sequences


def min_hamming_within(sequences: list, block: int = 256) -> np.ndarray:
    """
    min_hamming_within(): minimum Hamming distance from each sequence to any OTHER sequence in the list.

    Parameters:
    sequences (list of str): List of equal-length DNA sequences (A/C/G/T).
    block (int, optional): # of sequences compared per block, trades memory for speed (Default: 256)

    Dependencies: numpy, encode_sequences(), hamming_matrix()
    """
    n = len(sequences)
    if n < 2:
        return np.full(n, np.nan)

    encoded = encode_sequences(sequences)
    length = encoded.shape[1]
    out = np.empty(n, dtype=np.int32)
    for start in range(0, n, block):
        stop = min(start + block, n)
        distances = hamming_matrix(encoded[start:stop], encoded).astype(np.int32)
        np.fill_diagonal(distances[:, start:stop], length + 1)  # ignore each sequence's distance to itself
        out[start:stop] = distances.min(axis=1)
    return out


def gc_fractions(sequences: list) -> np.ndarray:
    """
    gc_fractions(): GC content of each sequence as a fraction of its length.

    Parameters:
    sequences (list of str): List of equal-length DNA sequences (A/C/G/T).

    Dependencies: numpy & encode_sequences()
    """
    if len(sequences) == 0:
        return np.empty(0)
    encoded = encode_sequences(sequences)
    return ((encoded == 1) | (encoded == 2)).sum(axis=1) / encoded.shape[1]  # C=1, G=2


def max_homopolymers(sequences: list) -> np.ndarray:
    """
    max_homopolymers(): longest run of a single repeated base in each sequence.

    Parameters:
    sequences (list of str): List of equal-length DNA sequences (A/C/G/T).

    Dependencies: numpy & encode_sequences()
    """
    if len(sequences) == 0:
        return np.empty(0, dtype=np.int16)
    encoded = encode_sequences(sequences)
    longest = np.ones(encoded.shape[0], dtype=np.int16)
    current = np.ones(encoded.shape[0], dtype=np.int16)
    for position in range(1, encoded.shape[1]):
        current = np.where(encoded[:, position] == encoded[:, position - 1], current + 1, 1)
        longest = np.maximum(longest, current)
    return longest


def fast_filter_by_hamming(sequences: list, min_distance: int, kept: list = None,
                           return_min_distances: bool = False):
    """
    fast_filter_by_hamming(): greedily filter sequences so that every retained sequence is more than
    'min_distance' away from every other retained sequence, using NumPy for speed.

    NOTE: the comparison is strictly greater than ('> min_distance'), so the smallest Hamming distance
    in the returned set is min_distance + 1. This is the long-standing behavior of this function.

    Parameters:
    sequences (list of str): List of equal-length DNA sequences (A/C/G/T) to consider, in priority order.
    min_distance (int): Hamming distance that retained sequences must exceed.
    kept (list, optional): sequences already retained from a previous call; they are kept as-is and
        included in the returned list (Default: None, start a new set)
    return_min_distances (bool, optional): also return, for each input sequence, its minimum Hamming
        distance to the retained set at the moment it was evaluated (np.nan if the set was empty)
        (Default: False)

    Dependencies: numpy, encode_sequences() & hamming_matrix()
    """
    kept_ls = list(kept) if kept is not None else []

    if not sequences:
        return (kept_ls, np.empty(0)) if return_min_distances else kept_ls

    encoded = encode_sequences(sequences)
    length = encoded.shape[1]
    min_distances = np.full(len(sequences), np.nan)

    # Buffer holding the encoded retained set; grows in place as sequences are added
    kept_encoded = encode_sequences(kept_ls) if kept_ls else np.empty((0, length), dtype=np.uint8)
    if kept_encoded.shape[0] and kept_encoded.shape[1] != length:
        raise ValueError("fast_filter_by_hamming() requires 'kept' sequences of the same length as 'sequences'.")
    buffer = np.empty((len(kept_ls) + len(sequences), length), dtype=np.uint8)
    buffer[:len(kept_ls)] = kept_encoded
    n_kept = len(kept_ls)

    # Vectorized prefilter: anything already too close to the pre-existing retained set can never be added
    if n_kept:
        prefilter = np.empty(len(sequences), dtype=np.int32)
        block = max(1, int(2e7 // n_kept))  # cap each distance matrix at ~40 MB
        for start in range(0, len(sequences), block):
            stop = min(start + block, len(sequences))
            prefilter[start:stop] = hamming_matrix(encoded[start:stop], buffer[:n_kept]).min(axis=1)
        min_distances[:] = prefilter
        candidates = np.flatnonzero(prefilter > min_distance)
    else:
        candidates = np.arange(len(sequences))

    # Sequential pass over the survivors, which only have to be checked against sequences added here
    n_start = n_kept
    for i in candidates:
        if n_kept > n_start:
            distances = np.count_nonzero(buffer[n_start:n_kept] != encoded[i], axis=1)
            closest = distances.min()
            if not np.isnan(min_distances[i]):
                closest = min(closest, min_distances[i])
            min_distances[i] = closest
            if closest <= min_distance:
                continue
        buffer[n_kept] = encoded[i]
        n_kept += 1
        kept_ls.append(sequences[i])

    return (kept_ls, min_distances) if return_min_distances else kept_ls

def count_csv_rows(pt:str) -> int:
    '''
    count_csv_rows(): counts number of rows in a csv file (subtracting 1 for header)

    Parameters:
    pt (str): file path to csv file
    '''
    with open(pt, 'r') as f:
        return sum(1 for line in f) - 1  # subtract 1 for header


def umi_report(kept: list, GC_fract: tuple, hamming: int, min_distances: np.ndarray = None,
               exclude: list = None, exclude_hamming: int = None):
    """
    umi_report(): print stage-independent summary statistics for a set of UMIs.

    Parameters:
    kept (list of str): retained UMI sequences
    GC_fract (tuple): pair of GC content boundaries written as fractions
    hamming (int): Hamming distance that retained UMIs were required to exceed
    min_distances (ndarray, optional): minimum Hamming distance within the retained set, one per UMI;
        recomputed here if not provided (Default: None)
    exclude (list of str, optional): reserved sequences the retained UMIs had to stay away from (Default: None)
    exclude_hamming (int, optional): Hamming distance the retained UMIs had to exceed from 'exclude' (Default: None)

    Dependencies: numpy, pandas, min_hamming_within(), min_hamming_to(), gc_fractions() & max_homopolymers()
    """
    if not kept:
        print('  No UMIs retained; nothing to report.')
        return

    print(f'\n--- UMI set summary (n = {len(kept)}) ---')

    # Hamming distance distribution (minimum distance from each UMI to any other retained UMI)
    if min_distances is None:
        min_distances = min_hamming_within(sequences=kept)
    valid = min_distances[~np.isnan(min_distances)]
    if valid.size:
        print('  Minimum Hamming distance to the nearest other UMI:')
        print(f'    min {int(valid.min())} | median {np.median(valid):.1f} | mean {valid.mean():.2f} | max {int(valid.max())}')
        counts = pd.Series(valid.astype(int)).value_counts().sort_index()
        for distance, count in counts.items():
            print(f'    d = {distance}: {count} UMIs ({count/len(kept)*100:.1f}%)')
        if valid.min() <= hamming:
            print(f'    WARNING: {int((valid <= hamming).sum())} UMIs are within the requested distance of {hamming}.')

    # Distance to the reserved set, if one was supplied
    if exclude:
        to_exclude = min_hamming_to(sequences=kept, reference=exclude)
        print(f'  Minimum Hamming distance to the {len(exclude)} excluded sequences:')
        print(f'    min {int(to_exclude.min())} | median {np.median(to_exclude):.1f} | mean {to_exclude.mean():.2f} | max {int(to_exclude.max())}')
        if exclude_hamming is not None and to_exclude.min() <= exclude_hamming:
            print(f'    WARNING: {int((to_exclude <= exclude_hamming).sum())} UMIs are within the requested distance of {exclude_hamming} of an excluded sequence.')

    # GC content distribution
    gc = gc_fractions(sequences=kept)
    print(f'  GC fraction (requested strictly between {GC_fract[0]} and {GC_fract[1]}):')
    print(f'    min {gc.min():.2f} | median {np.median(gc):.2f} | mean {gc.mean():.3f} | max {gc.max():.2f}')
    gc_counts = pd.Series(gc).round(3).value_counts().sort_index()
    for fraction, count in gc_counts.items():
        print(f'    GC = {fraction:.3f}: {count} UMIs ({count/len(kept)*100:.1f}%)')

    # Per-position base composition; a biased position shows up as a value far from 0.25
    encoded = encode_sequences(kept)
    print('  Per-position base frequency:')
    print(f'    {"pos":>4}  {"A":>6}  {"C":>6}  {"G":>6}  {"T":>6}')
    for position in range(encoded.shape[1]):
        frequencies = [np.count_nonzero(encoded[:, position] == code)/len(kept) for code in range(4)]
        print(f'    {position+1:>4}  ' + '  '.join(f'{frequency:>6.3f}' for frequency in frequencies))

    # Homopolymer runs
    runs = max_homopolymers(sequences=kept)
    print('  Longest homopolymer run:')
    run_counts = pd.Series(runs).value_counts().sort_index()
    for run, count in run_counts.items():
        print(f'    {run}x: {count} UMIs ({count/len(kept)*100:.1f}%)')
    print('-' * 40 + '\n')


def umi(length: int = 10, GC_fract: tuple = (0.4, 0.6), hamming: int = 4,
        nrows: int = 1000, mode: str = 'filter', pt: str = None,
        exclude=None, exclude_col: str = 'UMI_sequence', exclude_hamming: int = None,
        dir: str = '../out', file: str = None, verbose: bool = False) -> pd.DataFrame:
    '''
    umi(): generates unique molecular identifiers (UMIs) of specified length, GC content, and Hamming distance

    Parameters:
    length (int, optional): length of the unique molecular identifiers (Default: 10)
    GC_fract (tuple, optional): pair of GC content boundaries written as fractions; the bounds are
        exclusive, so (0.4, 0.6) at length 10 keeps only 5 G/C bases (Default: (0.4, 0.6))
    hamming (int, optional): Hamming distance that retained UMIs must exceed; the smallest distance in
        the output is therefore hamming + 1 (Default: 4)
    nrows (int, optional): # of UMIs to compare iteratively for hamming filtering (Default: 1000)
    mode (str, optional): 'filter' drops every UMI that fails the GC, exclusion or Hamming filter and
        writes only the survivors; 'keep' writes every candidate UMI annotated with GC_fract, GC_pass,
        min_hamming, max_homopolymer, keep and, when 'exclude' is given, min_hamming_exclude, so the
        thresholds can be revisited without regenerating (Default: 'filter')
    pt (str, optional): shuffled UMI candidate file path if already made, reused so a run can be
        resumed in the same order (Default: None)
    exclude (str | DataFrame | Series | list, optional): reserved UMI sequences the output must stay
        away from, as a file path, dataframe, series or list; useful for generating a new UMI set that
        will not be confused with one already in use (Default: None)
    exclude_col (str, optional): column holding the sequences when 'exclude' is a file path or
        dataframe (Default: 'UMI_sequence')
    exclude_hamming (int, optional): Hamming distance that retained UMIs must exceed from every
        excluded sequence; 0 removes exact matches only (Default: None, use 'hamming')
    dir (str, optional): save directory (Default: '../out')
    file (str, optional): output file name; built from the parameters if not provided (Default: None)
    verbose (bool, optional): print stage counts, timing, Hamming distance distribution, GC and
        per-position composition, and homopolymer runs (Default: False)

    Dependencies: pandas, numpy, os, time, datetime, io, generate_sequences(), shuffle(), gc_fractions(),
        max_homopolymers(), read_sequences(), min_hamming_to(), fast_filter_by_hamming(),
        min_hamming_within() & umi_report()
    '''
    modes = ('filter', 'keep')
    if mode not in modes:
        raise ValueError(f"umi() mode must be one of {modes}; got '{mode}'.")
    exclude_distance = hamming if exclude_hamming is None else exclude_hamming

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    start_time = time.time()

    def elapsed(message: str):
        if verbose:
            print(f'[{time.time()-start_time:7.1f}s] {message}')

    # Build the candidate pool, either from scratch or from a previously shuffled file
    if pt is None:
        sequences = generate_sequences(length=length)  # Generate all possible sequences of specified length
        elapsed(f'Generated {len(sequences)} sequences of length {length}.')
        candidates = shuffle(ls=sequences)
        io.save(obj=pd.DataFrame({'UMI_sequence': candidates}),
                file=os.path.join(dir, f'{timestamp}_UMI_{length}_candidates.csv'))
        elapsed(f'Shuffled and saved {len(candidates)} candidate sequences.')
    else:
        candidates = io.get(pt=pt)['UMI_sequence'].tolist()
        candidates_set = set(candidates)
        elapsed(f'Loaded {len(candidates)} candidate sequences from {pt}.')
        missing = [sequence for sequence in generate_sequences(length=length) if sequence not in candidates_set]
        if missing:  # Top up the pool if the loaded file predates a wider search space
            candidates.extend(shuffle(ls=missing))
            elapsed(f'Added {len(missing)} sequences of length {length} missing from {pt}.')

    # GC filter; in 'keep' mode nothing is dropped, the result is recorded instead
    gc = gc_fractions(sequences=candidates)
    gc_pass = (gc > GC_fract[0]) & (gc < GC_fract[1])
    print(f'{int(gc_pass.sum())} of {len(candidates)} sequences of length {length} have GC content strictly between {GC_fract[0]} and {GC_fract[1]}.')
    hamming_positions = np.flatnonzero(gc_pass)
    hamming_input = [candidates[position] for position in hamming_positions]
    elapsed(f'GC filter kept {len(hamming_input)} sequences.')

    # Exclusion filter; candidates too close to a reserved set can never be retained, so they are
    # removed before the greedy pass rather than competing for space in it
    exclude_sequences = read_sequences(exclude, col=exclude_col, length=length)
    exclude_distances = np.full(len(candidates), np.nan)
    if exclude_sequences:
        distances = min_hamming_to(sequences=hamming_input, reference=exclude_sequences)
        exclude_distances[hamming_positions] = distances
        too_close = distances <= exclude_distance
        hamming_positions = hamming_positions[~too_close]
        hamming_input = [candidate for candidate, close in zip(hamming_input, too_close) if not close]
        print(f'{int(too_close.sum())} of {len(too_close)} GC-passing sequences are within {exclude_distance} of the {len(exclude_sequences)} excluded sequences.')
        if not hamming_input:
            # A set produced by this function at distance d covers the whole pool at radius d: every
            # candidate it rejected was rejected for being within d of something it kept. Excluding a
            # complete previous run at the same distance therefore always leaves nothing behind.
            print(f'No candidates survive the exclusion filter. If the excluded sequences are a complete '
                  f'umi() run at hamming={exclude_distance}, they already cover the whole search space at '
                  f'that distance; exclude only the sequences actually in use, lower exclude_hamming, or '
                  f'raise length.')
        elapsed(f'Exclusion filter kept {len(hamming_input)} sequences.')

    # Hamming filter, processed in chunks of nrows so a long run checkpoints as it goes
    out_file = file if file is not None else f'{timestamp}_UMI_{length}_hamming_{hamming}_{mode}.csv'
    checkpoint = os.path.join(dir, out_file if mode == 'filter' else f'{timestamp}_UMI_{length}_hamming_{hamming}_checkpoint.csv')
    kept = []
    eval_distances = np.full(len(hamming_input), np.nan)
    for i in np.arange(start=0, stop=len(hamming_input), step=nrows):
        chunk = hamming_input[i:i+nrows]
        kept, chunk_distances = fast_filter_by_hamming(sequences=chunk, min_distance=hamming,
                                                       kept=kept, return_min_distances=True)
        eval_distances[i:i+len(chunk)] = chunk_distances

        # Checkpoint the retained set after each chunk so a long run can be interrupted
        print(f'Kept {len(kept)} sequences so far...')
        elapsed(f'Compared {min(i+nrows, len(hamming_input))} of {len(hamming_input)} sequences.')
        io.save(obj=pd.DataFrame({'UMI_sequence': kept}), file=checkpoint)

    elapsed(f'Hamming filter kept {len(kept)} of {len(hamming_input)} sequences.')

    # Recompute each retained UMI's distance to the final set; the value recorded during the greedy
    # pass only reflects the UMIs retained before it, which can overstate the true minimum
    kept_distances = min_hamming_within(sequences=kept) if len(kept) <= 20000 else None
    if kept_distances is None and verbose:
        print(f'  Skipped the final distance recomputation: {len(kept)} retained UMIs is above the 20000 limit.')
    if kept_distances is not None:
        elapsed('Recomputed minimum Hamming distances within the retained set.')

    # Assemble and save the output
    kept_set = set(kept)
    if mode == 'filter':
        df = pd.DataFrame({'UMI_sequence': kept})
    else:
        kept_index = {sequence: i for i, sequence in enumerate(kept)}
        min_hamming = np.full(len(candidates), np.nan)
        keep_flags = np.zeros(len(candidates), dtype=bool)
        min_hamming[hamming_positions] = eval_distances  # distance at the moment each candidate was evaluated
        for position, candidate in zip(hamming_positions, hamming_input):
            if candidate in kept_index:
                keep_flags[position] = True
                if kept_distances is not None:
                    min_hamming[position] = kept_distances[kept_index[candidate]]
        df = pd.DataFrame({'UMI_sequence': candidates,
                           'GC_fract': gc,
                           'GC_pass': gc_pass,
                           'min_hamming': min_hamming,
                           'max_homopolymer': max_homopolymers(sequences=candidates),
                           'keep': keep_flags})
        if exclude_sequences:  # distance to the reserved set; NaN where the GC filter already dropped the candidate
            df.insert(3, 'min_hamming_exclude', exclude_distances)
    io.save(obj=df, file=os.path.join(dir, out_file))
    if mode == 'keep' and os.path.exists(checkpoint):
        os.remove(checkpoint)  # the annotated table supersedes the checkpoint
    print(f'Kept {len(kept_set)} UMIs of length {length} (Hamming distance > {hamming}); wrote {len(df)} rows to {os.path.join(dir, out_file)}.')

    if verbose:
        umi_report(kept=kept, GC_fract=GC_fract, hamming=hamming, min_distances=kept_distances,
                   exclude=exclude_sequences, exclude_hamming=exclude_distance)
        print(f'Total time: {time.time()-start_time:.1f}s')

    return df


# Master Mix
def pcr_mm(primers: pd.Series, template_uL: int, template: str='1-2 ng/uL template',
           Q5_mm_x_stock: int=5, dNTP_mM_stock: int=10, fwd_uM_stock: int=10, rev_uM_stock: int=10, Q5_U_uL_stock: int=2,
           Q5_mm_x_desired: int=1, dNTP_mM_desired: float=0.2,fwd_uM_desired: float=0.5, rev_uM_desired: float=0.5, Q5_U_uL_desired: float=0.02,
           total_uL: int=25, mm_x: float=1.1) -> dict[pd.DataFrame]:
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
                                                     },index=pd.Index(list(np.arange(1,9)), name=f"{pcr1_fwd}_{pcr1_rev}"))
    return pcr_mm_dc

def rxn_mm(df: pd.DataFrame, rxn_uL: float, samples: int=1, mm_x: float=1.1,
           water: str='Nuclease-free water', component_col: str='Component',
           stock_col: str='Stock', desired_col: str='Desired',
           unit_col: str='Unit', uL_col: str='uL', mm_col: str='In MM',
           note_col: str='Note', file: str=None) -> pd.DataFrame:
    '''
    rxn_mm(): master mix for any reaction given stock and desired concentrations.

    Every reaction set up at the bench is the same table: components with a
    stock and a desired concentration, whose volume is rxn_uL * desired/stock,
    a few components added at a fixed volume, and water making up the rest. A
    PCR, a Golden Gate assembly and an oligo anneal differ only in the rows.

    Water is the balance rather than an input, so a component that changes
    cannot quietly push the reaction over volume; it comes back as an error
    instead.

    Parameters:
    df (dataframe): one row per component. A row with both a stock and a
        desired concentration is calculated; a row with a volume in uL_col is
        taken as given. A row with neither is an error.
    rxn_uL (float): total volume of one reaction
    samples (int, optional): reactions the master mix has to cover (Default: 1)
    mm_x (float, optional): master mix multiplier (Default: 1.1)
    water (str, optional): name of the balance row (Default: 'Nuclease-free water')
    component_col (str, optional): component name column (Default: 'Component')
    stock_col (str, optional): stock concentration column (Default: 'Stock')
    desired_col (str, optional): desired concentration column (Default: 'Desired')
    unit_col (str, optional): unit column (Default: 'Unit')
    uL_col (str, optional): fixed per-reaction volume column (Default: 'uL')
    mm_col (str, optional): optional boolean column; False keeps a component
        out of the master mix, for anything added per sample (Default: 'In MM')
    note_col (str, optional): optional note column carried through (Default: 'Note')
    file (str, optional): output file path

    Dependencies: pandas,io
    '''
    if component_col not in df.columns:
        raise KeyError(f"rxn_mm() needs a '{component_col}' column; df has "
                       f"{', '.join(map(str, df.columns))}")
    if rxn_uL <= 0:
        raise ValueError(f'rxn_uL must be greater than 0, got {rxn_uL}')
    if samples < 1:
        raise ValueError(f'samples must be at least 1, got {samples}')

    def _num(row, col):
        if col not in df.columns:
            return None
        value = row.get(col)
        if value is None or value == '' or pd.isna(value):
            return None
        return float(value)

    names, stocks, desireds, units, volumes, in_mm, notes = [], [], [], [], [], [], []
    for _, row in df.iterrows():
        name = row[component_col]
        stock, desired, fixed = (_num(row, stock_col), _num(row, desired_col),
                                 _num(row, uL_col))
        if fixed is not None:
            uL = fixed
        elif stock is not None and desired is not None:
            if stock <= 0:
                raise ValueError(f"'{name}' has a stock of {stock}; it must be "
                                 f'greater than 0')
            if desired > stock:
                raise ValueError(f"'{name}' wants {desired} from a stock of "
                                 f'{stock}; it cannot be reached by dilution')
            uL = rxn_uL * desired / stock
        else:
            raise ValueError(f"'{name}' needs either a volume in '{uL_col}' or "
                             f"both '{stock_col}' and '{desired_col}'")
        names.append(name)
        stocks.append(row.get(stock_col, ''))
        desireds.append(row.get(desired_col, ''))
        units.append(row.get(unit_col, '') if unit_col in df.columns else '')
        volumes.append(round(uL, 2))
        keep = row.get(mm_col) if mm_col in df.columns else True
        in_mm.append(True if keep is None or pd.isna(keep) else bool(keep))
        notes.append(row.get(note_col, '') if note_col in df.columns else '')

    water_uL = round(rxn_uL - sum(volumes), 2)
    if water_uL < 0:
        raise ValueError(f'The components come to {sum(volumes):.2f} uL, which '
                         f'is more than the {rxn_uL} uL reaction. Raise rxn_uL '
                         f'or lower a desired concentration.')

    # Water first, the way it goes into the tube; total last
    names = [water] + names + ['Total']
    stocks = [''] + stocks + ['']
    desireds = [''] + desireds + ['']
    units = [''] + units + ['']
    volumes = [water_uL] + volumes + [round(rxn_uL, 2)]
    in_mm = [True] + in_mm + [True]
    notes = [''] + notes + ['']

    scale = samples * mm_x
    mm_uL = [round(v * scale, 2) if keep else '-'
             for v, keep in zip(volumes, in_mm)]
    # The total is only what goes in the master mix tube, and is scaled from
    # the unrounded sum so it is the volume to expect rather than the sum of
    # nine rounded cells
    mm_uL[-1] = round(sum(v for v, keep in zip(volumes[:-1], in_mm[:-1]) if keep)
                      * scale, 2)
    mm_uL = pd.Series(mm_uL, dtype=object)

    def _blank(values):
        return pd.Series(['' if v is None or (not isinstance(v, str) and pd.isna(v))
                          else v for v in values], dtype=object)

    df_mm = pd.DataFrame({'Component': names,
                          'Stock': _blank(stocks),
                          'Desired': _blank(desireds),
                          'Unit': _blank(units),
                          'uL': volumes,
                          f'uL, MM ({samples} x {mm_x})': mm_uL,
                          'Note': _blank(notes)})
    if not any(n for n in notes if isinstance(n, str) and n):
        df_mm = df_mm.drop(columns='Note')
    if file is not None:
        io.save(obj=df_mm, file=file, id=True)
    return df_mm

# Simulation
def pcr_sim(df: pd.DataFrame | str,template_col: str, fwd_bind_col: str, rev_bind_col: str,
            fwd_ext_col: str=None, rev_ext_col: str=None, product_col: str='PCR Product',
            primer_in_product: bool=True, file:str=None) -> pd.DataFrame:
    '''
    pcr_sim(): returns dataframe with simulated pcr product 
    
    Parameters:
    df (dataframe | str): dataframe with template & primers (or file path)
    template_col (str): template column name
    fwd_bind_col (str): fwd primer binding region column name 
    rev_bind_col (str): rev primer binding region column name 
    fwd_ext_col (str, optional): fwd primer extension region column name (Default: None)
    rev_ext_col (str, optional): rev primer extension region column name (Default: None)
    product_col (str, optional): pcr product column name (Default: 'PCR Product')
    primer_in_product (bool, optional): include primer binding and extension regions in PCR product (Default: True)
    file (str, optional): output file path

    Dependencies: pandas,Bio.Seq,tidy
    '''
    if type(df)==str: # Get template & primers dataframe from file path if needed
        df = io.get(pt=df)

    pcr_product_ls = []
    
    if primer_in_product == False: # Don't include PCR primers in products; useful for alignments
       for (template,fwd_bind,rev_bind) in t.zip_cols(df=df,cols=[template_col,fwd_bind_col,rev_bind_col]):
           rc_rev_bind = ''.join(Seq(rev_bind).reverse_complement())
           pcr_product_ls.append(template[template.find(fwd_bind)+len(fwd_bind):template.find(rc_rev_bind)]) # template between primers
    elif fwd_ext_col is not None and rev_ext_col is not None: # FWD & REV primers have extension regions
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

    # Save & return dataframe
    io.save(obj=df, file=file) 
    return df

def off_targets(df: pd.DataFrame | str, col: str, match_score: float = 2, mismatch_score: float = -1, 
                open_gap_score: float = -10, extend_gap_score: float = -0.1, ckpt: int = 100,
                dir: str = None, return_df: bool = True) -> pd.DataFrame:
    '''
    off_targets(): Find off-target sequences for a list of sequences of the same length using pairwise alignment.
    
    Parameters:
    df (dataframe | str): DataFrame with sequences (or file path)
    col (str): Column name containing sequences to align.
    match_score (float, optional): Score for a match in the alignment (Default: 2).
    mismatch_score (float, optional): Penalty for a mismatch in the alignment (Default: -1).
    open_gap_score (float, optional): Penalty for opening a gap in the alignment (Default: -10).
    extend_gap_score (float, optional): Penalty for extending a gap in the alignment (Default: -0.1).
    ckpt (int, optional): Checkpoint interval for saving progress (Default: 100).
    dir (str, optional): Directory to save alignment checkpoints (Default: None).
    return_df (bool, optional): Whether to return a DataFrame with the results (Default: True).
    '''
    # High sequence homology; punish gaps
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = match_score  # Score for a match
    aligner.mismatch_score = mismatch_score  # Penalty for a mismatch; applied to both strands
    aligner.open_gap_score = open_gap_score  # Penalty for opening a gap; applied to both strands
    aligner.extend_gap_score = extend_gap_score  # Penalty for extending a gap; applied to both strands

    if type(df) == str:  # Get sequences dataframe from file path if needed
        df = io.get(pt=df)
    if col not in df.columns:  # Check if the specified column exists
        raise ValueError(f"Column '{col}' not found in the DataFrame. Available columns: {df.columns.tolist()}")
    seqs = df[col].tolist()  # Extract sequences from the specified column

    # Store target sequences, off-target sequences, their alignments, and their alignment scores
    off_target_seqs = []
    off_target_seqs_scores = []
    off_target_seqs_alignments = []
    for s,seq in enumerate(seqs): # Iterate though sequences
        if s==0: # Initial alignment status
            print(f'{s+1} out of {len(seqs)}') 
        elif s%(ckpt)==0: # Alignment status; save checkpoint
            print(f'{s+1} out of {len(seqs)}')
            if dir is not None:
                io.save(obj=pd.DataFrame({'Target Sequence': seqs[:s],
                                            'Off Target Sequence': off_target_seqs,
                                            'Best Alignment Score': off_target_seqs_scores,
                                            'Best Alignment Formatted': off_target_seqs_alignments}),
                        file=os.path.join(dir,'ckpt', f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_seqs_{s+1}.csv'))
            
        if seq is None: # Missing region (not applicable to count_alignments())
            continue

        seq_alignments_scores = []
        seq_alignments_aligned = []
        for ref in seqs: # Iterate though reference sequences
            if ref != seq:
                seq_alignment = aligner.align(seq.upper(),ref.upper()) # trim ngs sequence to reference sequence & align
                seq_alignments_scores.append(seq_alignment[0].score) # Save highest alignment score
                seq_alignments_aligned.append(seq_alignment.sequences[1]) # Save alignment matches
            else: # Skip self-alignment
                continue

        # Isolate maximum score alignment
        i = seq_alignments_scores.index(max(seq_alignments_scores))
        off_target_seqs.append(seq_alignments_aligned[i])
        off_target_seqs_scores.append(seq_alignments_scores[i])
        off_target_seqs_alignments.append(fq.format_alignment(a=seq, b=off_target_seqs[s], show=False, return_alignment=True))

    if dir is not None:  # Save final results
        io.save(obj=pd.DataFrame({'Target Sequence': seqs,
                                  'Off Target Sequence': off_target_seqs,
                                  'Best Alignment Score': off_target_seqs_scores,
                                  'Best Alignment Formatted': off_target_seqs_alignments}),
                file=os.path.join(dir, f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_seqs_{s+1}.csv'))
    
    if return_df:  # Return results as a DataFrame
        return pd.DataFrame({'Target Sequence': seqs,
                             'Off Target Sequence': off_target_seqs,
                             'Best Alignment Score': off_target_seqs_scores,
                             'Best Alignment Formatted': off_target_seqs_alignments})

def anneal_mm(oligos: int=1, rxn_uL: float=10, mm_x: float=1.05,
              top_uL: float=1.0, bottom_uL: float=1.0, buffer_stock: float=10,
              buffer_desired: float=1, pnk_uL: float=0.5,
              file: str=None) -> dict[str, pd.DataFrame]:
    '''
    anneal_mm(): phosphorylate and anneal oligo duplexes.

    Top and bottom oligos go in per tube, so they are listed at their
    per-reaction volume and left out of the master mix; everything else is
    shared. Returns the mix and the thermocycler program.

    Parameters:
    oligos (int, optional): duplexes to anneal (Default: 1)
    rxn_uL (float, optional): volume of one reaction (Default: 10)
    mm_x (float, optional): master mix multiplier (Default: 1.05)
    top_uL (float, optional): top oligo per reaction, uL (Default: 1.0)
    bottom_uL (float, optional): bottom oligo per reaction, uL (Default: 1.0)
    buffer_stock (float, optional): T4 ligase buffer stock, x (Default: 10)
    buffer_desired (float, optional): T4 ligase buffer in the reaction, x (Default: 1)
    pnk_uL (float, optional): T4 PNK per reaction, uL (Default: 0.5)
    file (str, optional): output file path

    Dependencies: pandas,io,rxn_mm()
    '''
    components = pd.DataFrame([
        {'Component': 'Top oligo (100 uM)', 'uL': top_uL, 'In MM': False},
        {'Component': 'Bottom oligo (100 uM)', 'uL': bottom_uL, 'In MM': False},
        {'Component': 'T4 DNA ligase buffer', 'Stock': buffer_stock,
         'Desired': buffer_desired, 'Unit': 'x'},
        {'Component': 'T4 PNK', 'uL': pnk_uL, 'Note': 'add last'},
    ])
    program = pd.DataFrame([
        {'Step': 1, 'Temperature': '37 C', 'Time': '30 min', 'Note': 'phosphorylate'},
        {'Step': 2, 'Temperature': '95 C', 'Time': '5 min', 'Note': 'denature'},
        {'Step': 3, 'Temperature': '4 C', 'Time': 'hold',
         'Note': 'ramp at 0.2 C/sec so the duplexes anneal rather than snap back'},
    ])
    tables = {'Anneal': rxn_mm(components, rxn_uL=rxn_uL, samples=oligos, mm_x=mm_x),
              'Anneal program': program}
    if file is not None:
        io.save(obj=tables, file=file, id=True)
    return tables

def golden_gate(inserts: int=3, samples: int=1, rxn_uL: float=10,
                mm_x: float=1.1, insert_uL: float=None, insert_stock: float=10,
                insert_desired: float=1, insert_unit: str='x',
                vector: str='Vector', vector_stock: float=100,
                vector_ng: float=100, buffer_stock: float=10,
                buffer_desired: float=1, atp_stock: float=10,
                atp_desired: float=1, enzyme: str='FastDigest Esp3I',
                enzyme_stock: float=20, enzyme_desired: float=1,
                ligase_stock: float=20, ligase_desired: float=1,
                shared_inserts: bool=False, cycles: int=15,
                file: str=None) -> dict[str, pd.DataFrame]:
    '''
    golden_gate(): Esp3I/T4 ligase assembly from a vector and fragments.

    Covers both routes into the same assembly: inserts that are annealed oligo
    duplexes, which go in at a dilution (insert_stock/insert_desired), and an
    insert that is a gel-purified PCR product, which goes in at a fixed volume
    (insert_uL). Returns the mix and the thermocycler program.

    The vector is the one component specified as an amount rather than a
    concentration, because that is how it is decided at the bench: so many
    nanograms of backbone, whatever the miniprep came back at.

    Parameters:
    inserts (int, optional): fragments going into the vector (Default: 3)
    samples (int, optional): assemblies the master mix has to cover (Default: 1)
    rxn_uL (float, optional): volume of one reaction (Default: 10)
    mm_x (float, optional): master mix multiplier (Default: 1.1)
    insert_uL (float, optional): fixed insert volume per reaction; set this for
        a gel-purified PCR product and the stock/desired pair is ignored. Zero
        or None means the insert is a dilution instead, so a template can carry
        the choice as a number (Default: None)
    insert_stock (float, optional): insert stock, in insert_unit (Default: 10)
    insert_desired (float, optional): insert in the reaction (Default: 1)
    insert_unit (str, optional): insert concentration unit (Default: 'x')
    vector (str, optional): backbone name (Default: 'Vector')
    vector_stock (float, optional): backbone concentration, ng/uL (Default: 100)
    vector_ng (float, optional): backbone per reaction, ng (Default: 100)
    buffer_stock (float, optional): FastDigest buffer stock, x (Default: 10)
    buffer_desired (float, optional): FastDigest buffer in the reaction (Default: 1)
    atp_stock (float, optional): ATP stock, mM (Default: 10)
    atp_desired (float, optional): ATP in the reaction, mM (Default: 1)
    enzyme (str, optional): type IIS enzyme (Default: 'FastDigest Esp3I')
    enzyme_stock (float, optional): enzyme stock, x (Default: 20)
    enzyme_desired (float, optional): enzyme in the reaction, x (Default: 1)
    ligase_stock (float, optional): T4 DNA ligase stock, x (Default: 20)
    ligase_desired (float, optional): T4 DNA ligase in the reaction, x (Default: 1)
    shared_inserts (bool, optional): True when every assembly takes the same
        fragments, which puts them in the master mix; False when each assembly
        has its own, which keeps them per tube (Default: False)
    cycles (int, optional): digest/ligate cycles (Default: 15)
    file (str, optional): output file path

    Dependencies: pandas,io,rxn_mm()
    '''
    if inserts < 1:
        raise ValueError(f'inserts must be at least 1, got {inserts}')
    if vector_stock <= 0:
        raise ValueError(f'vector_stock must be greater than 0, got {vector_stock}')

    rows = [{'Component': 'FastDigest buffer', 'Stock': buffer_stock,
             'Desired': buffer_desired, 'Unit': 'x'},
            {'Component': 'ATP', 'Stock': atp_stock, 'Desired': atp_desired,
             'Unit': 'mM'}]
    for i in range(1, inserts + 1):
        row = {'Component': f'Fragment {i}', 'In MM': shared_inserts}
        if insert_uL:
            row['uL'] = insert_uL
        else:
            row.update({'Stock': insert_stock, 'Desired': insert_desired,
                        'Unit': insert_unit})
        rows.append(row)
    rows.append({'Component': vector, 'uL': round(vector_ng / vector_stock, 2),
                 'Note': f'{vector_ng:g} ng at {vector_stock:g} ng/uL'})
    rows.append({'Component': enzyme, 'Stock': enzyme_stock,
                 'Desired': enzyme_desired, 'Unit': 'x'})
    rows.append({'Component': 'T4 DNA ligase', 'Stock': ligase_stock,
                 'Desired': ligase_desired, 'Unit': 'x'})

    program = pd.DataFrame([
        {'Step': 1, 'Temperature': '37 C', 'Time': '5 min', 'Note': 'digest'},
        {'Step': 2, 'Temperature': '16 C', 'Time': '5 min',
         'Note': f'ligate; repeat steps 1-2 for {cycles} cycles'},
        {'Step': 3, 'Temperature': '4 C', 'Time': 'hold', 'Note': ''},
    ])
    tables = {'Golden Gate': rxn_mm(pd.DataFrame(rows), rxn_uL=rxn_uL,
                                    samples=samples, mm_x=mm_x),
              'Golden Gate program': program}
    if file is not None:
        io.save(obj=tables, file=file, id=True)
    return tables

# Grams per mole per base pair of double-stranded DNA, the usual round number
NG_PER_PMOL_KB = 650

def gibson(df: pd.DataFrame, rxn_uL: float=15, mm_uL: float=None,
           pmol: float=0.02, pmol_small: float=0.05, small_kb: float=1.0,
           fragment_col: str='Fragment', kb_col: str='kb',
           conc_col: str='ng/uL', pmol_col: str='pmol',
           dilution_col: str='Dilution', mm: str='2x NEBuilder HiFi MM',
           file: str=None) -> pd.DataFrame:
    '''
    gibson(): fragment volumes for an equimolar NEBuilder assembly.

    An assembly is equimolar in fragments, not in mass, so a 5 kb backbone and
    a 0.3 kb insert at the same number of picomoles are fifteen-fold apart in
    nanograms. This converts each fragment's target picomoles to nanograms at
    650 ng per pmol per kb, then to a volume at its measured concentration,
    and gives the water that fills the reaction.

    Short fragments default to a higher pmol target than long ones, which is
    the usual way of giving a small insert its excess over the backbone.

    Parameters:
    df (dataframe): one row per fragment, with a length in kb and a measured
        concentration. A per-fragment pmol or dilution column overrides the
        defaults for that row.
    rxn_uL (float, optional): volume of one assembly (Default: 15)
    mm_uL (float, optional): 2x master mix volume; half the reaction when left
        unset, which is what a 2x mix means (Default: None)
    pmol (float, optional): target picomoles for a fragment at or above
        small_kb (Default: 0.02)
    pmol_small (float, optional): target picomoles below small_kb (Default: 0.05)
    small_kb (float, optional): length below which pmol_small applies (Default: 1.0)
    fragment_col (str, optional): fragment name column (Default: 'Fragment')
    kb_col (str, optional): fragment length column, kb (Default: 'kb')
    conc_col (str, optional): measured concentration column, ng/uL (Default: 'ng/uL')
    pmol_col (str, optional): optional per-fragment pmol column (Default: 'pmol')
    dilution_col (str, optional): optional per-fragment dilution factor, where 10
        means a 1-in-10 dilution of the stock (Default: 'Dilution')
    mm (str, optional): master mix name (Default: '2x NEBuilder HiFi MM')
    file (str, optional): output file path

    Dependencies: pandas,io,NG_PER_PMOL_KB
    '''
    for col in (fragment_col, kb_col, conc_col):
        if col not in df.columns:
            raise KeyError(f"gibson() needs a '{col}' column; df has "
                           f"{', '.join(map(str, df.columns))}")
    if len(df) == 0:
        raise ValueError('gibson() needs at least one fragment')
    if mm_uL is None:
        mm_uL = rxn_uL / 2
    if mm_uL >= rxn_uL:
        raise ValueError(f'mm_uL ({mm_uL}) must be less than rxn_uL ({rxn_uL})')

    names, kbs, pmols, ngs, concs, uLs = [], [], [], [], [], []
    for _, row in df.iterrows():
        name, kb = row[fragment_col], float(row[kb_col])
        if kb <= 0:
            raise ValueError(f"'{name}' has a length of {kb} kb; it must be "
                             f'greater than 0')
        conc = row[conc_col]
        if conc is None or pd.isna(conc) or float(conc) <= 0:
            raise ValueError(f"'{name}' has no measured concentration. Nanodrop "
                             f'the fragments before planning the assembly; the '
                             f'volumes cannot be guessed from the lengths alone.')
        conc = float(conc)
        want = row.get(pmol_col) if pmol_col in df.columns else None
        if want is None or pd.isna(want):
            want = pmol_small if kb < small_kb else pmol
        want = float(want)
        dilution = row.get(dilution_col) if dilution_col in df.columns else None
        dilution = 1.0 if dilution is None or pd.isna(dilution) else float(dilution)
        if dilution <= 0:
            raise ValueError(f"'{name}' has a dilution of {dilution}; it must "
                             f'be greater than 0')
        ng = want * kb * NG_PER_PMOL_KB
        names.append(name)
        kbs.append(round(kb, 3))
        pmols.append(round(want, 4))
        ngs.append(round(ng, 2))
        concs.append(round(conc / dilution, 2))
        uLs.append(round(ng / (conc / dilution), 2))

    water_uL = round(rxn_uL - mm_uL - sum(uLs), 2)
    if water_uL < 0:
        raise ValueError(
            f'The fragments come to {sum(uLs):.2f} uL, which leaves no room '
            f'beside {mm_uL} uL of master mix in a {rxn_uL} uL reaction. '
            f'Concentrate the fragments, lower the pmol target, or scale the '
            f'reaction up.')

    df_gibson = pd.DataFrame({
        'Fragment': names + [mm, 'Nuclease-free water', 'Total'],
        'kb': pd.Series(kbs + ['', '', ''], dtype=object),
        'pmol': pd.Series(pmols + ['', '', ''], dtype=object),
        'ng': pd.Series(ngs + ['', '', ''], dtype=object),
        'ng/uL': pd.Series(concs + ['', '', ''], dtype=object),
        'uL': uLs + [round(mm_uL, 2), water_uL, round(rxn_uL, 2)]})
    if file is not None:
        io.save(obj=df_gibson, file=file, id=True)
    return df_gibson
