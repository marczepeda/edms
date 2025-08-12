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
- epegRNA_pool(): design GG cloning oligonucleotides for pooled prime editing epegRNAs

[UMIs]
- generate_sequences(): recursively generates all possible sequences of A, T, C, G of the specified length
- filter_GC(): filters sequences based on GC content
- shuffle(): randomly reorganizes a list
- encode_sequences(): Convert sequences to integer arrays for fast comparison.
- fast_filter_by_hamming(): fastFilter sequences such that all retained sequences have a Hamming distance 'min_distance' from each other, using NumPy for speed.
- count_csv_rows(): counts number of rows in a csv file (subtracting 1 for header)
- UMI(): generates unique molecular identifiers (UMIs) of specified length, GC content, and Hamming distance

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
import datetime
from ..gen import io
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
        top_ids.append(pre+str(df.iloc[i][id])+'_top')
        bot_ids.append(pre+str(df.iloc[i][id])+'_bot')
        if (tG==True)&(s[0]!='G'): top_seqs.append(t5+'G'+s+t3)
        else: top_seqs.append(t5+s+t3)
        bot_seqs.append(b5+str(Seq(s).reverse_complement()+b3))
    df[pre+id+'_top']=top_ids
    df[pre+id+'_bot']=bot_ids
    df[pre+seq+'_top']=top_seqs
    df[pre+seq+'_bot']=bot_seqs
    
    return df

def sgRNAs(df:pd.DataFrame | str,id:str,spacer='Spacer_sequence',t5='CACC',t3='',b5='AAAC',b3='',tG=True,order=True,
           dir:str=None, file:str=None):
    ''' 
    sgRNAs(): design GG cloning oligonucleotides for cutting and base editing sgRNAs
    
    Parameters:
    df (dataframe | str): Dataframe with sequence information for sgRNAs (or file path)
    id (str): id column name
    spacer (str): spacer column name (Default: Spacer_sequence)
    t5 (str): top oligonucleotide 5' overhang
    t3 (str): top oligonucleotide 3' overhang
    b5 (str): bottom oligonucleotide 5' overhang (revcom)
    b3 (str): bottom oligonucleotide 3' overhang (revcom)
    tG (bool): add 5' G to spacer if needed (Default: True)
    order (bool): order format
    dir (str, optional): save directory
    file (str, optional): save file
    
    Dependencies: pandas, io, top_bot(), & ord_form()
    '''
    if type(df)==str: # Get sgRNAs dataframe from file path if needed
        df = io.get(pt=df)

    df=tb(df=df,id=id,seq=spacer,t5=t5,t3=t3,b5=b5,b3=b3,tG=tG,pre='o') # Make top and bottom oligos for spacer inserts
    if order==True: # Sigma order format (or original dataframe with top and bottom oligos)
        df = pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_top',pre='o'), 
                        ord_form(df=df,id=id,seq=spacer,suf='_bot',pre='o')]).reset_index(drop=True)
    
    # Save & return dataframe
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=df)  
    return df

def epegRNAs(df: pd.DataFrame | str,id: str,tG=True, order=True,make_extension=True,
             spacer='Spacer_sequence',spacer_t5='CACC',spacer_t3='GTTTAAGAGC',spacer_b5='',spacer_b3='',
             extension='Extension_sequence',RTT='RTT_sequence',PBS='PBS_sequence',linker='Linker_sequence',extension_t5='',extension_t3='',extension_b5='CGCG',extension_b3='GCACCGACTC',
             order_scaffold=False, dir:str=None, file:str=None):
    ''' 
    epegRNAs(): design GG cloning oligonucleotides for prime editing epegRNAs
    
    Parameters:
    df (dataframe | str): Dataframe with sequence information for epegRNAs (or file path)
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
    order_scaffold (bool, optional): order top and bottom oligonucleotide for scaffold sequence (Default: False)
    dir (str, optional): save directory
    file (str, optional): save file
    
    Assumptions:
    1. epegRNA scaffold: GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
    2. epegRNA motif: tevoPreQ1 (CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA)
    
    Dependencies: pandas, top_bot(), & ord_form()
    '''
    if type(df)==str: # Get epegRNAs dataframe from file path if needed
        df = io.get(pt=df)
    
    if make_extension==True: df[extension] = df[RTT]+df[PBS]+df[linker] # Make extension by concatenating RTT, PBS, and linker
    else: print(f'Warning: Did not make extension sequence!\nMake sure "{extension}" column includes RTT+PBS+linker for epegRNAs.')
    df=tb(df=df,id=id,seq=spacer,t5=spacer_t5,t3=spacer_t3,b5=spacer_b5,b3=spacer_b3,tG=tG,pre='es_') # Make top and bottom oligos for spacer inserts
    df=tb(df=df,id=id,seq=extension,t5=extension_t5,t3=extension_t3,b5=extension_b5,b3=extension_b3,tG=False,pre='ee_') # Make top and bottom oligos for extension inserts
    if order_scaffold==True: # Order top and bottom oligonucleotide for scaffold sequence
        df = pd.concat([ord_form(df=df,id=id,seq='TATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC',suf='_top',pre='scaffold_'),
                        ord_form(df=df,id=id,seq='GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC',suf='_bot',pre='scaffold_')]).reset_index(drop=True)
    if order==True: # Sigma order format (or original dataframe with top and bottom oligos)
        if order_scaffold==True: # Include the scaffold sequence in the order
            df = pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_top',pre='es_'),
                            ord_form(df=df,id=id,seq=spacer,suf='_bot',pre='es_'),
                            ord_form(df=df,id=id,seq=extension,suf='_top',pre='ee_'),
                            ord_form(df=df,id=id,seq=extension,suf='_bot',pre='ee_'),
                            pd.DataFrame({'Oligo Name': ['pegRNA_scaffold_top','pegRNA_scaffold_bot'],
                                          'Sequence': ['TATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC','GGTGCCACTTGGGTGGGATGCTCGCAGGCATTCAGCCAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCATAGCTCTTAAAC'],
                                          'Scale (µmol)': [0.05, 0.05],
                                          'bp': [99, 99]}
                            )]).reset_index(drop=True)
        else: # Do not include the scaffold sequence in the order
            df = pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_top',pre='es_'),
                            ord_form(df=df,id=id,seq=spacer,suf='_bot',pre='es_'),
                            ord_form(df=df,id=id,seq=extension,suf='_top',pre='ee_'),
                            ord_form(df=df,id=id,seq=extension,suf='_bot',pre='ee_')
                            ]).reset_index(drop=True)

    # Save & return dataframe
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=df)  
    return df

def ngRNAs(df: pd.DataFrame | str,id: str, tG=True, order=True,
           spacer='Spacer_sequence', spacer_t5='CACC',spacer_t3='GTTTAAGAGC',spacer_b5='',spacer_b3='',
           order_scaffold=False, dir:str=None, file:str=None):
    ''' 
    ngRNAs(): design GG cloning oligonucleotides for prime editing ngRNAs
    
    Parameters:
    df (dataframe | str): Dataframe with spacers (or file path)
    id (str): id column
    tG (bool, optional): add 5' G to spacer if needed (Default: True)
    order (bool, optional): order format (Default: True)
    spacer (str, optional): ngRNA spacer column name
    spacer_t5 (str, optional): top oligonucleotide 5' overhang
    spacer_t3 (str, optional): top oligonucleotide 3' overhang
    spacer_b5 (str, optional): bottom oligonucleotide 5' overhang
    spacer_b3 (str, optional): bottom oligonucleotide 3' overhang}
    order_scaffold (bool, optional): order top and bottom oligonucleotide for scaffold sequence (Default: False)
    dir (str, optional): save directory
    file (str, optional): save file
    
    Assumptions:
    1. ngRNA scaffold: GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
    
    Dependencies: pandas, top_bot(), & ord_form()
    '''
    if type(df)==str: # Get ngRNAs dataframe from file path if needed
        df = io.get(pt=df)
    
    df=tb(df=df,id=id,seq=spacer,t5=spacer_t5,t3=spacer_t3,b5=spacer_b5,b3=spacer_b3,tG=tG,pre='ns_') # Make top and bottom oligos for spacer inserts
    if order==True: # Sigma order format (or original dataframe with top and bottom oligos)
        if order_scaffold==True: # Include the scaffold sequence in the order
            df = pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_top',pre='ns_'),
                            ord_form(df=df,id=id,seq=spacer,suf='_bot',pre='ns_'),
                            pd.DataFrame({'Oligo Name': ['ngRNA_scaffold_top','ngRNA_scaffold_bot'],
                                          'Sequence': ['TATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC','AAAAGCACCGACTCGGTGCCACTTGGGTGGGATGCTCGCAGGCATTCAGCCAAGTTGATAACGGACTAGCCTTATTTAAACTTGCTATGCTGTTTCCAGCATAGCTCTTAAAC'],
                                          'Scale (µmol)': [0.05, 0.05],
                                          'bp': [99, 113]})
                            ]).reset_index(drop=True)
        else: # Do not include the scaffold sequence in the order
            df = pd.concat([ord_form(df=df,id=id,seq=spacer,suf='_top',pre='ns_'),
                            ord_form(df=df,id=id,seq=spacer,suf='_bot',pre='ns_')
                            ]).reset_index(drop=True)

    # Save & return dataframe
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=df)  
    return df

# Library GG cloning
def epegRNA_pool(df: pd.DataFrame | str, tG:bool=True, make_extension:bool=True,
                 UMI_df: pd.DataFrame | str='../resources/UMI_15_shuffle_hamming_5_100000.csv',
                 PCR_df: pd.DataFrame| str='../resources/edms_pcr.csv',
                 RE_type_IIS_df: pd.DataFrame| str='../resources/RE_type_IIS.csv',enzymes: list=['Esp3I'],
                 barcode:str='Barcode', barcode_i:int=0, fwd_barcode_t5:str='Forward Barcode',rev_barcode_t3:str='Reverse Barcode',
                 Esp3I_hU6:str='Esp3I_hU6',tevopreQ1_Esp3I:str='tevopreQ1_Esp3I',
                 epegRNA_spacer:str='Spacer_sequence',
                 epegRNA_extension:str='Extension_sequence',epegRNA_RTT:str='RTT_sequence',
                 epegRNA_PBS:str='PBS_sequence',epegRNA_linker:str='Linker_sequence',
                 dir:str=None, file:str=None, return_df:bool=False):
    ''' 
    epegRNA_pool(): design GG cloning oligonucleotides for pooled prime editing epegRNAs
    
    Parameters:
    df (dataframe | str): Dataframe with sequence information for epegRNAs & corresponding ngRNAs (or file path)
    tG (bool, optional): add 5' G to spacer if needed (Default: True)
    make_extension (bool, optional): concatenate RTT, PBS, and linker to make extension sequence (Default: True)
    UMI_df (dataframe | str, optional): Dataframe with UMI sequences (or file path) (Default: '../resources/UMI_15_shuffle_hamming_5_100000.csv')
    PCR_df (dataframe | str, optional): Dataframe with PCR primer and subpool barcode information (or file path) (Default: '../resources/edms_pcr.csv')
    RE_type_IIS_df (dataframe | str, optional): Dataframe with Type IIS restriction enzyme information (or file path) (Default: '../resources/RE_typeIIS.csv')
    enzymes (list, optional): list of Type IIS restriction enzymes to check for (Default: ['Esp3I'])
    barcode (str, optional): subpool barcode column name (Default: Barcode)
    barcode_i (int, optional): subpool barcode index (Default: 0)
    fwd_barcode_t5 (bool, optional): forward barcode column name (Default: Forward Barcode)
    rev_barcode_t3 (bool, optional): reverse barcode column name (Default: Reverse Barcode)
    Esp3I_hU6 (bool, optional): Esp3I_hU6 column name (Default: Esp3I_hU6)
    tevopreQ1_Esp3I (bool, optional): tevopreQ1_Esp3I column name (Default: tevopreQ1_Esp3I)
    epegRNA_spacer (str, optional): epegRNA spacer column name (Default: Spacer_sequence)
    epegRNA_extension (str, optional): epegRNA extension name (Default: Extension_sequence)
    epegRNA_RTT (str, optional): epegRNA reverse transcripase template column name (Default: RTT_sequence)
    epegRNA_PBS (str, optional): epegRNA primer binding site column name (Default: PBS_sequence)
    epegRNA_linker (str, optional): epegRNA linker column name (Default: Linker_sequence)
    dir (str, optional): save directory
    file (str, optional): save file

    Assumptions:
    1. Oligo Template: FWD UMI - FWD Barcode - Esp3I(F) - hU6 - epegRNA_spacer - epegRNA_scaffold - epegRNA_extension - tevopreQ1 motif - Esp3I(R) - REV Barcode - REV UMI
    2. epegRNA motif: tevoPreQ1 (CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA)
    
    Dependencies: pandas,
    '''
    # Get dataframes from file path if needed
    if type(df)==str:
        df = io.get(pt=df)
    if type(UMI_df)==str:
        UMI_df = io.get(pt=UMI_df)
    if type(PCR_df)==str:
        PCR_df = io.get(pt=PCR_df)
    if type(RE_type_IIS_df)==str:
        RE_type_IIS_df = io.get(pt=RE_type_IIS_df)

    # Make extension by concatenating RTT, PBS, and linker
    if make_extension==True: df[epegRNA_extension] = df[epegRNA_RTT]+df[epegRNA_PBS]+df[epegRNA_linker]
    else: print(f'Warning: Did not make extension sequence!\nMake sure "{epegRNA_extension}" column includes RTT+PBS+linker for epegRNAs.')

    # Assign subpool barcodes
    barcodes = sorted(df[barcode].unique())
    PCR_barcodes = PCR_df.iloc[barcode_i:barcode_i+len(barcodes)]
    PCR_barcodes[barcode] = barcodes
    df = pd.merge(left=df,right=PCR_barcodes,on=barcode)

    # Assign UMI to each twist oligo
    UMI_sequences = UMI_df['UMI_sequence'].tolist()
    if len(UMI_sequences)<df.shape[0]*2: 
        ValueError(f'{len(UMI_sequences)} UMIs is less than 2x{df.shape[0]} twist oligonucleotides!')

    # Make twist oligo & determine length
    df.sort_values(by=barcode,inplace=True,ignore_index=True)
    df[f'{epegRNA_spacer}_nt1']=[s[0] for s in df[epegRNA_spacer]]
    twist_oligos = []
    for i,(epegRNA_spacer_nt1) in enumerate(df[f'{epegRNA_spacer}_nt1']):
        if tG==True: # Append 5'G to spacer if not already present
            if epegRNA_spacer_nt1=='G':
                twist_oligos.append(df.iloc[i]['Forward UMI']+
                                    df.iloc[i][fwd_barcode_t5]+df.iloc[i][Esp3I_hU6]+
                                    df.iloc[i][epegRNA_spacer]+df.iloc[i]['Scaffold_sequence']+
                                    df.iloc[i][epegRNA_extension]+df.iloc[i]['Sensor_sequence']+
                                    df.iloc[i][tevopreQ1_Esp3I]+df.iloc[i][rev_barcode_t3]+
                                    df.iloc[i]['Reverse UMI'])
            else:
                twist_oligos.append(df.iloc[i]['Forward UMI']+
                                    df.iloc[i][fwd_barcode_t5]+df.iloc[i][Esp3I_hU6]+'G'+
                                    df.iloc[i][epegRNA_spacer]+df.iloc[i]['Scaffold_sequence']+
                                    df.iloc[i][epegRNA_extension]+df.iloc[i]['Sensor_sequence']+
                                    df.iloc[i][tevopreQ1_Esp3I]+df.iloc[i][rev_barcode_t3]+
                                    df.iloc[i]['Reverse UMI'])
        else: # Do not append 5'G to spacer if not already present
            twist_oligos.append(df.iloc[i]['Forward UMI']+
                                df.iloc[i][fwd_barcode_t5]+df.iloc[i][Esp3I_hU6]+
                                df.iloc[i][epegRNA_spacer]+df.iloc[i]['Scaffold_sequence']+
                                df.iloc[i][epegRNA_extension]+df.iloc[i]['Sensor_sequence']+
                                df.iloc[i][tevopreQ1_Esp3I]+df.iloc[i][rev_barcode_t3]+
                                df.iloc[i]['Reverse UMI'])
            
    df['Twist_oligo'] = twist_oligos
    df['Twist_oligo_length']=[len(twist) for twist in df['Twist_oligo']]

    # Check for 2 recognition sites per enzyme
    for (enzyme,recognition,recognition_rc) in t.zip_cols(df=RE_type_IIS_df[RE_type_IIS_df['Name'] in enzymes],cols=['Name','Recognition','Recognition_rc']):
        enzyme_sites = [len(t.find_all(twist_oligo,recognition)) + # Check forward direction
                        len(t.find_all(twist_oligo,recognition_rc)) # Check reverse direction
                        for twist_oligo in twist_oligos] # Iterate through Twist product
        df[enzyme] = enzyme_sites
    
    # Save & return dataframe
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=df)
    if return_df:
        return df

# UMIs
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

def encode_sequences(seq_list: list):
    """
    encode_sequences(): Convert sequences to integer arrays for fast comparison.
    
    Parameters:
    seq_list (list of str): List of DNA sequences (A/C/G/T).
    """
    base_to_int = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return np.array([[base_to_int[base] for base in seq] for seq in seq_list], dtype=np.uint8)

def fast_filter_by_hamming(sequences: list, min_distance: int, dir: str = '../out'):
    """
    fast_filter_by_hamming(): fastFilter sequences such that all retained sequences have a Hamming distance 'min_distance' from each other, using NumPy for speed.
    
    Parameters:
    sequences (list of str): List of equal-length DNA sequences (A/C/G/T).
    min_distance (int): Minimum required Hamming distance.
    dir (str, optional): save directory for intermediate results (Default: '../out')
    """
    if not sequences:
        return []

    encoded = encode_sequences(sequences)
    keep_indices = [0]
    
    for i in range(1, len(encoded)):
        current = encoded[i]
        kept = encoded[keep_indices]
        
        # Compute Hamming distances in batch
        distances = np.count_nonzero(kept != current, axis=1)
        
        if np.all(distances > min_distance):
            keep_indices.append(i)
        
        if len(keep_indices) % 1000 == 0: # Save iteratively every 1000 sequences
            print(f'Kept {len(keep_indices)} sequences so far...')
            io.save(dir=dir, 
                    file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_UMI_{len(sequences[0])}_hamming_{min_distance}_compare_{len(sequences)}_yield_{len(keep_indices)}.csv',
                    obj=pd.DataFrame({'UMI_sequence': [sequences[i] for i in keep_indices]}))
    
    return [sequences[i] for i in keep_indices]

def count_csv_rows(pt:str):
    '''
    count_csv_rows(): counts number of rows in a csv file (subtracting 1 for header)

    Parameters:
    pt (str): file path to csv file
    '''
    with open(pt, 'r') as f:
        return sum(1 for line in f) - 1  # subtract 1 for header

def UMI(length: int = 15, GC_fract: tuple = (0.4, 0.6), shuffle: bool = True,
        hamming: int = 5, nrows: int=1000, pt: str=None, dir: str = '../out'):
    '''
    UMI(): generates unique molecular identifiers (UMIs) of specified length, GC content, and Hamming distance
    
    Parameters:
    length (int, optional): length of the unique molecular identifiers (Default: 15)
    GC_fract (tuple, optional): pair of GC content boundaries written as fractions (Default: (0.4, 0.6))
    shuffle (bool, optional): randomly reorganize the list of UMIs (Default: True)
    hamming (int, optional): Minimum Hamming distance between UMIs (Default: 5)
    nrows (int, optional): # of UMIs to compare iteratively for hamming filtering (Default: 1000)
    pt (str, optional): Shuffled UMI file path if already made (Default: None)
    dir (str, optional): save directory
    '''
    if pt is None: # Start from scratch if no file path is provided
        sequences = generate_sequences(length=length) # Generate all possible sequences of specified length
        filtered_sequences = filter_GC(sequences=sequences, GC_fract=GC_fract) # Filter sequences based on GC content
        print(f'Generated {len(filtered_sequences)} sequences of length {length} with GC content between {GC_fract[0]} and {GC_fract[1]}.')

        if shuffle is True: # Shuffle the sequences if specified & save
            filtered_sequences = shuffle(ls=filtered_sequences)
            io.save(dir=dir, 
                    file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_UMI_{length}_shuffle.csv', 
                    obj=pd.DataFrame({'UMI_sequence': filtered_sequences}))
        else:
            io.save(dir=dir, 
                    file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_UMI_{length}.csv', 
                    obj=pd.DataFrame({'UMI_sequence': filtered_sequences}))
    
        #filtered_sequences = filtered_sequences[:nrows] # Limit to nrows for hamming filtering

    else: # Load from provided file path
        #filtered_sequences = io.get(pt=pt,nrows=nrows)['UMI_sequence'].tolist()
        filtered_sequences = io.get(pt=pt)['UMI_sequence'].tolist()
    
    filtered_sequences_save = [] # List to save filtered sequences iteratively
    for i in np.arange(start=0, stop=count_csv_rows(pt), step=nrows): # Process in chunks of nrows
        filtered_sequences_save.extend( # Filter sequences based on Hamming distance & compare to previously saved sequences
            fast_filter_by_hamming(sequences = filtered_sequences_save + filtered_sequences[i:i+nrows],
                                   min_distance = hamming,
                                   dir = dir)
        )
    #filtered_sequences = fast_filter_by_hamming(sequences=filtered_sequences, min_distance=hamming,dir=dir) # Filter sequences based on Hamming distance

    io.save(dir=dir, 
            file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_UMI_{length}_hamming_{hamming}_compare_{nrows}_yield_{len(filtered_sequences_save)}.csv', 
            obj=pd.DataFrame({'UMI_sequence': filtered_sequences_save}))

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
                                                     },index=pd.Index(list(np.arange(1,9)), name=f"{pcr1_fwd}_{pcr1_rev}"))
    return pcr_mm_dc

# Simulation
def pcr_sim(df: pd.DataFrame | str,template_col: str, fwd_bind_col: str, rev_bind_col: str,
            fwd_ext_col: str=None, rev_ext_col: str=None, product_col='PCR Product',
            dir:str=None, file:str=None):
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
    dir (str, optional): save directory
    file (str, optional): save file

    Dependencies: pandas,Bio.Seq,tidy
    '''
    if type(df)==str: # Get template & primers dataframe from file path if needed
        df = io.get(pt=df)

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

    # Save & return dataframe
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=df) 
    return df