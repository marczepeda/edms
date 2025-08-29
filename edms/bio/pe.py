''' 
Module: pe.py
Author: Marc Zepeda
Created: 2024-08-31
Description: Prime Editing

Usage:
[Biological Dictionaries]
- dna_aa_codon_table: DNA to AA codon table
- aa_dna_codon_table: AA to DNA codon table

[Helper Functions]
- get_codons(): returns all codons within a specified frame for a nucleotide sequence
- get_codon_frames(): returns all codon frames for a nucleotide sequence
- found_list_in_order(): returns index of sub_ls found consecutive order in main_ls or -1
- find_enzyme_sites(): find enzyme sites in pegRNAs or ngRNAs
- enzyme_codon_swap(): modify pegRNA RTT sequences to disrupt a RE recognition site
 
[PrimeDesign]
- prime_design_input(): creates and checks PrimeDesign saturation mutagenesis input file
- prime_design(): run PrimeDesign using Docker (NEED TO BE RUNNING DESKTOP APP)
- prime_design_output(): splits peg/ngRNAs from PrimeDesign output & finishes annotations
- prime_designer(): execute PrimeDesign for EDMS using Docker (NEED TO BE RUNNING DESKTOP APP)
- merge(): rejoins epeg/ngRNAs & creates ngRNA_groups

[pegRNA]
- epegRNA_linkers(): generate epegRNA linkers between PBS and 3' hairpin motif & finish annotations
- shared_sequences(): Reduce PE library into shared spacers and PBS sequences
- pilot_screen(): Create pilot screen for EDMS
- sensor_designer(): Design pegRNA sensors
- rtt_designer(): design all possible RTT for given spacer & PBS (WT, single insertions, & single deletions)
- pegRNA_outcome(): confirm that pegRNAs should create the predicted edit
- pegRNA_signature(): create signatures for pegRNA outcomes using alignments

[Comparing pegRNA libraries]
- print_shared_sequences(): prints spacer and PBS sequences from dictionary of shared_sequences libraries
- print_shared_sequences_mutant(): prints spacer and PBS sequences as well as priority mutant from dictionary of shared_sequences libraries

[Comparing pegRNAs]
- group_pe(): returns a dataframe containing groups of (epegRNA,ngRNA) pairs that share spacers and have similar PBS and performs pairwise alignment for RTT  
'''

# Import packages
import pandas as pd
import numpy as np
import os
import re
import datetime
from typing import Literal
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import math
from typing import Literal

from ..bio.signature import signature_from_alignment
from ..bio import pegLIT as pegLIT
from ..gen import io as io
from ..gen import tidy as t
from ..gen import plot as p
from ..dat import cosmic as co 
from ..dat import cvar
from ..bio import fastq as fq
from ..utils import memory_timer,load_resource_csv

# Biological Dictionaries
''' dna_aa_codon_table: DNA to AA codon table'''
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

''' aa_dna_codon_table: AA to DNA codon table'''
aa_dna_codon_table = {
    "F": ["TTT", "TTC"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "Y": ["TAT", "TAC"],
    "*": ["TAA", "TAG", "TGA"],  # Stop codons
    "C": ["TGT", "TGC"],
    "W": ["TGG"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "H": ["CAT", "CAC"],
    "Q": ["CAA", "CAG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "I": ["ATT", "ATC", "ATA"],
    "M": ["ATG"],  # Start codon
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "N": ["AAT", "AAC"],
    "K": ["AAA", "AAG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "G": ["GGT", "GGC", "GGA", "GGG"]
}

# Helper Functions 
def get_codons(sequence: str, frame: int=0) -> list[str]:
    ''' 
    get_codons(): returns all codons within a specified frame for a nucleotide sequence
    
    Parameters:
    sequence (str): nucletide sequence
    frame (int, optional): codon frame (0, 1, or 2)

    Dependencies:
    '''
    return [sequence[i:i+3] for i in range(frame, len(sequence) - 2, 3)]

def get_codon_frames(sequence: str) -> list[list[str]]:
    ''' 
    get_codon_frames(): returns all codon frames for a nucleotide sequence
    
    Parameters:
    seqeuence (str): nucleotide sequence

    Dependencies:
    ''' 
    return [get_codons(sequence,frame) for frame in range(3)]

def found_list_in_order(main_ls: list, sub_ls: list) -> int:
    ''' 
    found_list_in_order(): returns index of sub_ls found consecutive order in main_ls or -1
    
    Parameters:
    main_ls (list): search for it here
    sub_ls (list): find this list

    Dependencies:
    '''
    found=False # Initialize found variable
    for m,item in enumerate(main_ls): # Iterate through main_ls
        
        s=0 # Start index for sub_ls
        if item == sub_ls[0]: # If item matches sub_ls
            
            for sub_item in sub_ls: # Iterate through sub_ls
                try:
                    if sub_item == main_ls[m+s]: # Check sub_ls and main_ls item match
                        if s == 0: # Return index
                            index = m
                        
                        s+=1 # Increment s to check next item in main_ls 
                        
                    else: # If item does not match sub_ls, break
                        break
                    
                    if s+1 == len(sub_ls): # If last item in sub_ls; found True
                        found=True

                except IndexError: # End of main_ls reached
                    return -1
            
        if found==True: # Found all items in order
            return index 
    
    return -1 # Not found

def find_enzyme_sites(df: pd.DataFrame | str, enzyme: str, RE_type_IIS_df: pd.DataFrame | str = None, literal_eval: bool=True) -> pd.DataFrame:
    ''' 
    find_enzyme_sites(): find enzyme sites in pegRNAs or ngRNAs
    
    Parameters:
    df (pd.DataFrame | str): DataFrame with pegRNAs or ngRNAs or file path to DataFrame
    enzyme (str): Enzyme name (e.g. Esp3I, BsaI, BspMI, etc.)
    RE_type_IIS_df (pd.DataFrame | str, optional): DataFrame with Type IIS RE information (or file path)
    literal_eval (bool, optional): convert string representations (Default: True)
    '''
    # Get dataframes from file path if needed
    if type(df)==str:
        df = io.get(pt=df, literal_eval=literal_eval)

    if type(RE_type_IIS_df)==str:
        RE_type_IIS_df = io.get(pt=RE_type_IIS_df, literal_eval=literal_eval)
    elif RE_type_IIS_df is None: # Get from resources if not provided
        RE_type_IIS_df = load_resource_csv(filename='RE_type_IIS.csv')

    # Check forward & reverse direction for recognition sites on pegRNAs
    df_enzyme_sites_fwd = [t.find_all(oligo,RE_type_IIS_df[RE_type_IIS_df['Name']==enzyme]['Recognition'].values[0]) for oligo in df['Oligonucleotide']] # Iterate through oligonucleotides
    df_enzyme_sites_rc = [t.find_all(oligo,RE_type_IIS_df[RE_type_IIS_df['Name']==enzyme]['Recognition_rc'].values[0]) for oligo in df['Oligonucleotide']] # Iterate through oligonucleotides
    df_enzyme_sites = [len(enzyme_site_fwd)+len(enzyme_site_rc) for (enzyme_site_fwd,enzyme_site_rc) in zip(df_enzyme_sites_fwd,df_enzyme_sites_rc)] # Sum forward & reverse direction
    
    # Add enzyme sites to DataFrame & return
    df[enzyme] = df_enzyme_sites
    df[f'{enzyme}_fwd_i'] = df_enzyme_sites_fwd
    df[f'{enzyme}_rc_i'] = df_enzyme_sites_rc
    return df

def enzyme_codon_swap(pegRNAs: pd.DataFrame | str, in_file: pd.DataFrame | str, enzyme: str, 
                      RE_type_IIS_df: pd.DataFrame | str = None, out_dir: str = None, 
                      out_file: str = None, return_df: bool = True, literal_eval: bool=True, comments: bool=False) -> pd.DataFrame:
    '''
    enzyme_codon_swap(): modify pegRNA RTT sequences to disrupt a RE recognition site

    Parameters:
    pegRNAs (pd.DataFrame | str): pegRNAs DataFrame or file path to pegRNAs DataFrame
    in_file (dataframe | str): Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (column names required)
    enzyme (str): Enzyme name (e.g. Esp3I, BsaI, BspMI, etc.)
    RE_type_IIS_df (dataframe | str, optional): Dataframe with Type IIS RE information (or file path)
    out_dir (str, optional): output directory
    out_file (str, optional): output filename
    return_df (bool, optional): Return pegRNAs DataFrame (Default: True)
    literal_eval (bool, optional): convert string representations (Default: True)
    comments (bool, optional): Print comments (Default: False)
    '''
    # Initialize timer; memory reporting
    memory_timer(reset=True)
    memories = []

    # Get pegRNAs & PrimeDesign input DataFrame from file path if needed
    if type(pegRNAs)==str:
        pegRNAs = io.get(pt=pegRNAs, literal_eval=literal_eval)
    if type(in_file)==str:
        in_file = io.get(pt=in_file, literal_eval=literal_eval)
    if type(RE_type_IIS_df)==str:
        RE_type_IIS_df = io.get(pt=RE_type_IIS_df, literal_eval=literal_eval)
    elif RE_type_IIS_df is None: # Get from resources if not provided
            RE_type_IIS_df = load_resource_csv(filename='RE_type_IIS.csv')

    # Filter pegRNAs based on enzyme count
    pegRNAs = pegRNAs[pegRNAs[enzyme] == 1]

    enzyme_rtt = []
    for (spacer,scaffold,rtt,enzyme_fwd_i_ls,enzyme_rc_i_ls) in t.zip_cols(df=pegRNAs,cols=['Spacer_sequence','Scaffold_sequence','RTT_sequence',f'{enzyme}_fwd_i',f'{enzyme}_rc_i']):
        if len(enzyme_fwd_i_ls) == 1: 
            if (enzyme_fwd_i_ls[0] >= len(spacer)+ len(scaffold) - 1 - len(RE_type_IIS_df[RE_type_IIS_df['Name']==enzyme]['Recognition'])) & \
               (enzyme_fwd_i_ls[0] <= len(spacer) + len(scaffold) + len(rtt) - 1): # enzyme site is completely or partially in the RTT
                enzyme_rtt.append(True)
            else:
                enzyme_rtt.append(False)
        elif len(enzyme_rc_i_ls) == 1:
            if (enzyme_rc_i_ls[0] >= len(spacer)+ len(scaffold) - 1 - len(RE_type_IIS_df[RE_type_IIS_df['Name']==enzyme]['Recognition_rc'])) & \
               (enzyme_rc_i_ls[0] <= len(spacer) + len(scaffold) + len(rtt) - 1):
                enzyme_rtt.append(True)
            else:
                enzyme_rtt.append(False)
        else:
            raise ValueError(f"Multiple enzyme indices found for {enzyme} in pegRNA: {enzyme_fwd_i_ls} (length = {len(enzyme_fwd_i_ls)}) & {enzyme_rc_i_ls} (length = {len(enzyme_rc_i_ls)})")

    pegRNAs = pegRNAs[enzyme_rtt]

    # Get reference sequence & codons (+ reverse complement)
    target_sequence = in_file.iloc[0]['target_sequence'] 
    seq = Seq(target_sequence.split('(')[1].split(')')[0]) # Break apart target sequences
    if len(seq)%3 != 0: raise(ValueError(f"Length of target sequence ({len(seq)}) must divisible by 3. Check input file."))
    flank5 = Seq(target_sequence.split('(')[0])
    if len(flank5)%3 != 0: raise(ValueError(f"Length of 5' flank ({len(flank5)}) must divisible by 3. Check input file."))
    flank3 = Seq(target_sequence.split(')')[1])
    if len(flank3)%3 != 0: raise(ValueError(f"Length of 3' flank ({len(flank3)}) must divisible by 3. Check input file."))

    f5_seq_f3_nuc = flank5 + seq + flank3  # Join full nucleotide reference sequence
    rc_f5_seq_f3_nuc = Seq.reverse_complement(f5_seq_f3_nuc) # Full nucleotide reference reverse complement sequence
    
    codons = get_codons(seq) # Codons
    codons_flank5 = get_codons(flank5) # Codons in-frame flank 5
    codons_flank3 = get_codons(flank3) # Codons in-frame flank 3
    extended_codons = codons_flank5 + codons + codons_flank3 # Codons including flank 5 and flank 3

    # Get new RTT sequences
    new_rtt_ls = []
    for (strand,spacer,scaffold,rtt,pbs,rtt_length,enzyme_fwd_i_ls,enzyme_rc_i_ls) in t.zip_cols(df=pegRNAs,
                                                                                      cols=['Strand','Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence',
                                                                                            'RTT_length',f'{enzyme}_fwd_i',f'{enzyme}_rc_i']):
        if strand=='+': # Spacer: + strand; PBS & RTT: - strand
            
            # Find spacer in sequence
            spacer_j = f5_seq_f3_nuc.find(spacer)
            if spacer_j == -1:
                raise ValueError(f"Spacer sequence '{spacer}' not found in target sequence. Please check the input file.")
            elif spacer_j != f5_seq_f3_nuc.rfind(spacer):
                raise ValueError(f"Multiple matches found for spacer sequence '{spacer}'. Please check the input file.")

            # Find PBS in reverse complement sequence
            pbs_j = rc_f5_seq_f3_nuc.find(pbs)
            if pbs_j == -1:
                raise ValueError(f"PBS sequence '{pbs}' not found in target sequence. Please check the input file.")
            elif pbs_j != rc_f5_seq_f3_nuc.rfind(pbs):
                print(pbs,pbs_j,rc_f5_seq_f3_nuc.rfind(pbs))
                raise ValueError(f"Multiple matches found for PBS sequence '{pbs}'. Please check the input file.")

            # Obtain reverse complement WT RTT & edit RTT in-frame from + strand
            rc_rtt = Seq.reverse_complement(Seq(rtt)) # reverse complement of rtt (+ strand)
            rc_rtt_codon_frames = get_codon_frames(rc_rtt) # codons

            rtt_wt = rc_f5_seq_f3_nuc[pbs_j-int(rtt_length):pbs_j]
            rc_rtt_wt = Seq.reverse_complement(rtt_wt) # reverse complement of rtt wt (+ strand)
            rc_rtt_wt_codon_frames = get_codon_frames(rc_rtt_wt) # codons
            if comments==True:
                print(f"Extended Codons (Here): {extended_codons[math.floor((len(rc_f5_seq_f3_nuc)-pbs_j)/3)-1:math.floor((len(rc_f5_seq_f3_nuc)-pbs_j+int(rtt_length))/3)]}")
            for i,(rc_rtt_wt_codon_frame,rc_rtt_codon_frame) in enumerate(zip(rc_rtt_wt_codon_frames,rc_rtt_codon_frames)): # Search for in-frame nucleotide sequence
                if comments==True:
                    print(f"rc_rtt_wt_codon_frame: {rc_rtt_wt_codon_frame}")
                
                index = found_list_in_order(extended_codons[math.floor((len(rc_f5_seq_f3_nuc)-pbs_j)/3)-1:math.floor((len(rc_f5_seq_f3_nuc)-pbs_j+rtt_length)/3)],rc_rtt_wt_codon_frame)
                if index != -1: # Codon frame from reverse complement of rtt matches extended codons of in-frame nucleotide sequence
                    rc_rtt_inframe_nuc_codons_flank5 = rc_rtt[:i] # Save codon frame flank 5'
                    rc_rtt_inframe_nuc_codons = rc_rtt_codon_frame # Save codon frame
                    rc_rtt_inframe_nuc_codons_flank3 = rc_rtt[i+3*len(rc_rtt_codon_frame):] # Save codon frame flank 3'
                    rc_rtt_inframe_nuc = Seq('').join(rc_rtt_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rc_rtt_inframe_prot = Seq.translate(rc_rtt_inframe_nuc) # Translate to in-frame protein sequence
                    
                    found=True
                    break
            
            if comments==True:
                print(f'Strand: {strand}')    
                print(f'Nucleotides (WT): {rc_rtt_wt}')
                print(f'Nucleotides (Edit): {rc_rtt}')
                print(f'Nucleotides 5\' of Codons (Edit): {rc_rtt_inframe_nuc_codons_flank5}')
                print(f'Nucleotides Codons (Edit): {rc_rtt_inframe_nuc_codons}')
                print(f'Nucleotides 3\' of Codons (Edit): {rc_rtt_inframe_nuc_codons_flank3}')
                print(f'Nucleotides In-Frame (Edit): {rc_rtt_inframe_nuc}')
                print(f'Amino Acids In-Frame (Edit): {rc_rtt_inframe_prot}')
            
            if found==False:
                raise(ValueError("RTT was not found."))
            
            # Find enzyme site in reverse complement sequence codons
            enzyme_i = str(rc_rtt_inframe_nuc).upper().find(RE_type_IIS_df[RE_type_IIS_df['Name']==enzyme]['Recognition'].values[0])
            if enzyme_i == -1: # Try reverse complement enzyme site
                enzyme_i = str(rc_rtt_inframe_nuc).upper().find(RE_type_IIS_df[RE_type_IIS_df['Name']==enzyme]['Recognition_rc'].values[0])
            
            if enzyme_i != -1: # Found enzyme site or reverse complement enzyme site
                enzyme_codon_i = math.floor(enzyme_i/3)
            else:
                new_rtt_ls.append(None)
                continue
            
            if comments==True:
                print(f'Enzyme site index: {enzyme_i}')
                print(f'Enzyme codon index: {enzyme_codon_i}')

            # Change codon swap enzyme site & save new RTT sequence
            codons = [str(codon).lower() for codon in aa_dna_codon_table[str(rc_rtt_inframe_prot[enzyme_codon_i])] if str(codon).upper() != str(rc_rtt_inframe_nuc_codons[enzyme_codon_i]).upper()]
            if comments==True:
                print(f'Codons: {codons}')
            if len(codons)!=0:
                rc_rtt_inframe_nuc_codons[enzyme_codon_i] = codons[0]
                if comments==True:
                    print(f'Nucleotides In-Frame (New): {Seq('').join(rc_rtt_inframe_nuc_codons)}')
                    print(f'Amino Acid In-Frame (New): {Seq.translate(Seq('').join(rc_rtt_inframe_nuc_codons))}')
                    print(f'RTT (New): {str(Seq.reverse_complement(Seq(rc_rtt_inframe_nuc_codons_flank5)+Seq('').join(rc_rtt_inframe_nuc_codons)+Seq(rc_rtt_inframe_nuc_codons_flank3)))}')
                    print(f'RTT (Old): {rtt}')
                    print(f'RTT (WT): {rtt_wt}')
                new_rtt_ls.append(str(Seq.reverse_complement(Seq(rc_rtt_inframe_nuc_codons_flank5)+
                                                             Seq('').join(rc_rtt_inframe_nuc_codons)+
                                                             Seq(rc_rtt_inframe_nuc_codons_flank3))))
            else:
                new_rtt_ls.append(None)

        elif strand=='-': # Spacer: - strand; PBS & RTT: + strand
            
            # Find spacer in sequence
            spacer_j = rc_f5_seq_f3_nuc.find(spacer)
            if spacer_j == -1:
                raise ValueError(f"Spacer sequence '{spacer}' not found in target sequence. Please check the input file.")
            if spacer_j != rc_f5_seq_f3_nuc.rfind(spacer):
                raise ValueError(f"Multiple matches found for spacer sequence '{spacer}' not found in target sequence. Please check the input file.")

            # Find PBS in sequence
            pbs_j = f5_seq_f3_nuc.find(pbs)
            if pbs_j == -1:
                raise ValueError(f"PBS sequence '{pbs}' not found in target sequence. Please check the input file.")
            if pbs_j != f5_seq_f3_nuc.rfind(pbs):
                print(pbs,pbs_j,f5_seq_f3_nuc.rfind(pbs))
                raise ValueError(f"Multiple matches found for PBS sequence '{pbs}' not found in target sequence. Please check the input file.")

            # Obtain WT RTT & edit RTT in-frame from + strand
            rtt_codon_frames = get_codon_frames(rtt) # codons

            rtt_wt = f5_seq_f3_nuc[pbs_j-int(rtt_length):pbs_j]
            rtt_wt_codon_frames = get_codon_frames(rtt_wt) # codons
            if comments==True:
                print(f"Extended Codons (Here): {extended_codons[math.ceil((pbs_j-rtt_length)/3)-1:math.ceil(pbs_j/3)]}")
            for i,(rtt_wt_codon_frame,rtt_codon_frame) in enumerate(zip(rtt_wt_codon_frames,rtt_codon_frames)): # Search for in-frame nucleotide sequence
                if comments==True:
                    print(f"rtt_wt_codon_frame: {rtt_wt_codon_frame}")

                index = found_list_in_order(extended_codons[math.ceil((pbs_j-rtt_length)/3)-1:math.ceil(pbs_j/3)],rtt_wt_codon_frame)
                if index != -1: # Codon frame from rtt matches extended codons of in-frame nucleotide sequence
                    rtt_inframe_nuc_codons_flank5 = rtt[:i] # Save codon frame flank 5'
                    rtt_inframe_nuc_codons = rtt_codon_frame # Save codon frame
                    rtt_inframe_nuc_codons_flank3 = rtt[i+3*len(rtt_codon_frame):] # Save codon frame flank 3'
                    rtt_inframe_nuc = Seq('').join(rtt_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rtt_inframe_prot = Seq.translate(rtt_inframe_nuc) # Translate to in-frame protein sequence
                    found=True
                    break
            
            if comments==True:
                print(f'Strand: {strand}')    
                print(f'Nucleotides: {rtt_wt}')
                print(f'Nucleotides (Edit): {rtt}')
                print(f'Nucleotides 5\' of Codons (Edit): {rtt_inframe_nuc_codons_flank5}')
                print(f'Nucleotides Codons (Edit): {rtt_inframe_nuc_codons}')
                print(f'Nucleotides 3\' of Codons (Edit): {rtt_inframe_nuc_codons_flank3}')
                print(f'Nucleotides In-Frame (Edit): {rtt_inframe_nuc}')
                print(f'Amino Acids In-Frame (Edit): {rtt_inframe_prot}')
            
            if found==False:
                raise(ValueError("RTT was not found."))
            
            # Find enzyme site in reverse complement sequence codons
            enzyme_i = str(rtt_inframe_nuc).upper().find(RE_type_IIS_df[RE_type_IIS_df['Name']==enzyme]['Recognition'].values[0])
            if enzyme_i == -1: # Try reverse complement enzyme site
                enzyme_i = str(rtt_inframe_nuc).upper().find(RE_type_IIS_df[RE_type_IIS_df['Name']==enzyme]['Recognition_rc'].values[0])
            
            if enzyme_i != -1: # Found enzyme site or reverse complement enzyme site
                enzyme_codon_i = math.floor(enzyme_i/3)
            else:
                new_rtt_ls.append(None)
                continue
            
            if comments==True:
                print(f'Enzyme site index: {enzyme_i}')
                print(f'Enzyme codon index: {enzyme_codon_i}')

            # Change codon swap enzyme site & save new RTT sequence
            codons = [str(codon).lower() for codon in aa_dna_codon_table[str(rtt_inframe_prot[enzyme_codon_i])] if str(codon).upper() != str(rtt_inframe_nuc_codons[enzyme_codon_i]).upper()]
            if comments==True:
                print(f'Codons: {codons}')
            if len(codons)!=0:
                rtt_inframe_nuc_codons[enzyme_codon_i] = codons[0]
                if comments==True:
                    print(f'Nucleotides In-Frame (New): {Seq('').join(rtt_inframe_nuc_codons)}')
                    print(f'Amino Acid In-Frame (New): {Seq.translate(Seq('').join(rtt_inframe_nuc_codons))}')
                    print(f'RTT (New): {str(Seq(rtt_inframe_nuc_codons_flank5)+Seq('').join(rtt_inframe_nuc_codons)+Seq(rtt_inframe_nuc_codons_flank3))}')
                    print(f'RTT (Old): {rtt}')
                    print(f'RTT (WT): {rtt_wt}')
                new_rtt_ls.append(str(Seq(rtt_inframe_nuc_codons_flank5)+
                                      Seq('').join(rtt_inframe_nuc_codons)+
                                      Seq(rtt_inframe_nuc_codons_flank3)))
            else:
                new_rtt_ls.append(None)

    # Update pegRNAs DataFrame with new RTT sequences
    pegRNAs['RTT_sequence'] = new_rtt_ls
    pegRNAs = pegRNAs[pegRNAs['RTT_sequence'].isna()==False].reset_index(drop=True)  # Filter out None RTT sequences

    # Update extension & oligonucleotide sequence
    if 'Linker_sequence' in pegRNAs.columns:  # epegRNA?
        pegRNAs['Extension_sequence'] = pegRNAs['RTT_sequence'] + pegRNAs['PBS_sequence'] + pegRNAs['Linker_sequence']
        pegRNAs['Oligonucleotide'] = pegRNAs['Spacer_sequence'] + pegRNAs['Scaffold_sequence'] + pegRNAs['RTT_sequence'] + pegRNAs['PBS_sequence'] + pegRNAs['Linker_sequence']  
        pegRNAs['Oligonucleotide'] = [str(oligo).upper() for oligo in pegRNAs['Oligonucleotide']] # Convert to uppercase
            
    else: # pegRNA
        pegRNAs['Extension_sequence'] = pegRNAs['RTT_sequence'] + pegRNAs['PBS_sequence']
        pegRNAs['Oligonucleotide'] = pegRNAs['Spacer_sequence'] + pegRNAs['Scaffold_sequence'] + pegRNAs['RTT_sequence'] + pegRNAs['PBS_sequence'] 
        pegRNAs['Oligonucleotide'] = [str(oligo).upper() for oligo in pegRNAs['Oligonucleotide']] # Convert to uppercase                                                                                                                  cols=['Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence'])]

    # Save & Return
    memories.append(memory_timer(task=f"enzyme_codon_swap()"))
    if out_dir is not None and out_file is not None:
        io.save(dir=os.path.join(out_dir,f'.{enzyme}_codon_swap'),
                file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_memories.csv',
                obj=pd.DataFrame(memories, columns=['Task','Memory, MB','Time, s']))
        io.save(dir=out_dir,file=out_file,obj=pegRNAs)
    if return_df==True: return pegRNAs

# PrimeDesign
def prime_design_input(target_name: str, flank5_sequence: str, 
                     target_sequence: str, flank3_sequence: str,
                     aa_index: int=1,
                     dir: str='.', file: str='prime_design_input.csv'):
    ''' 
    prime_design_input(): creates and checks PrimeDesign saturation mutagenesis input file
    
    Parameters:
    target_name (str): name of target
    flank5_sequence: in-frame nucleotide sequence with 5' of saturation mutagensis region (length must be divisible by 3)
    target_sequence (str): in-frame nucleotide sequence for the saturation mutagensis region (length must be divisible by 3)
    flank3_sequence: in-frame nucleotide sequence with 3' of saturation mutagensis region (length must be divisible by 3)
    aa_index (int, optional): 1st amino acid in target sequence index (Default: 1)
    dir (str, optional): name of the output directory 
    file (str, optional): name of the output file
    
    Dependencies: pandas & io
    
    Reference: https://github.com/pinellolab/PrimeDesign/tree/master/PrimeDesign
    '''
    # Check PrimeDesign saturation mutagenesis input file
    if len(flank5_sequence)%3 != 0: raise(ValueError(f"Length of flank5_sequence ({len(flank5_sequence)}) must divisible by 3."))
    if len(target_sequence)%3 != 0: raise(ValueError(f"Length of target_sequence ({len(target_sequence)}) must divisible by 3."))
    if len(flank5_sequence)%3 != 0: raise(ValueError(f"Length of flank3_sequence ({len(flank3_sequence)}) must divisible by 3."))

    # Create PrimeDesign saturation mutagenesis input file
    io.save(dir=dir,
            file=file,
            obj=pd.DataFrame({'target_name': [target_name],
                              'target_sequence': [f"{flank5_sequence}({target_sequence}){flank3_sequence}"],
                              'aa_index': [aa_index]}))

def prime_design(file: str, pbs_length_list: list = [],rtt_length_list: list = [], nicking_distance_minimum: int = 0,
                nicking_distance_maximum: int = 100, filter_c1_extension: bool = False, silent_mutation: bool = False,
                genome_wide_design: bool = False, saturation_mutagenesis: str = None, number_of_pegrnas: int = 3, number_of_ngrnas: int = 3,
                nicking_distance_pooled: int = 75, homology_downstream: int = 10, pbs_length_pooled: int = 14, rtt_max_length_pooled: int = 50,
                out_dir: str = './DATETIMESTAMP_PrimeDesign'):
    ''' 
    prime_design(): run PrimeDesign
    
    Parameters:
    file (str): input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (column names required)
    pbs_length_list (list, optional): list of primer binding site (PBS) lengths for the pegRNA extension. Example: 12 13 14 15
    rtt_length_list (list, optional): list of reverse transcription (RT) template lengths for the pegRNA extension. Example: 10 15 20
    nicking_distance_minimum (int, optional): minimum nicking distance for designing ngRNAs. (Default: 0 bp)
    nicking_distance_maximum (int, optional): maximum nicking distance for designing ngRNAs. (Default: 100 bp)
    filter_c1_extension (bool, optional): filter against pegRNA extensions that start with a C base. (Default: False)
    silent_mutation (bool, optional): introduce silent mutation into PAM assuming sequence is in-frame. Currently only available with SpCas9. (Default: False)
    genome_wide_design (bool, optional): whether this is a genome-wide pooled design. This option designs a set of pegRNAs per input without ranging PBS and RTT parameters.
    saturation_mutagenesis (str, optional): saturation mutagenesis design with prime editing (Options: 'aa', 'base').
    number_of_pegrnas (int, optional): maximum number of pegRNAs to design for each input sequence. The pegRNAs are ranked by 1) PAM disrupted > PAM intact then 2) distance to edit. (Default: 3)
    number_of_ngrnas (int, optional): maximum number of ngRNAs to design for each input sequence. The ngRNAs are ranked by 1) PE3b-seed > PE3b-nonseed > PE3 then 2) deviation from nicking_distance_pooled. (Default: 3)
    nicking_distance_pooled (int, optional): the nicking distance between pegRNAs and ngRNAs for pooled designs. PE3b annotation is priority (PE3b seed -> PE3b non-seed), followed by nicking distance closest to this parameter. (Default: 75 bp)
    homology_downstream (int, optional): for pooled designs (genome_wide or saturation_mutagenesis needs to be indicated), this parameter determines the RT extension length downstream of an edit for pegRNA designs. (Default: 10)
    pbs_length_pooled (int, optional): the PBS length to design pegRNAs for pooled design applications. (Default: 14 nt)
    rtt_max_length_pooled (int, optional): maximum RTT length to design pegRNAs for pooled design applications. (Default: 50 nt)
    out_dir (str, optional): name of output directory (Default: ./DATETIMESTAMP_PrimeDesign)
    
    Dependencies: os, numpy, & https://github.com/pinellolab/PrimeDesign
    '''
    # Write PrimeDesign Command Line
    cmd = 'python -m edms.bio.primedesign'
    cmd += f' -f {file}' # Append required parameters
    if pbs_length_list: cmd += f' -pbs {" ".join([str(val) for val in pbs_length_list])}' # Append optional parameters
    if rtt_length_list: cmd += f' -rtt {" ".join([str(val) for val in rtt_length_list])}'
    if nicking_distance_minimum!=0: cmd += f' -nick_dist_min {str(nicking_distance_minimum)}' 
    if nicking_distance_maximum!=100: cmd += f' -nick_dist_max {str(nicking_distance_maximum)}'
    if filter_c1_extension: cmd += f' -filter_c1 {str(filter_c1_extension)}'
    if silent_mutation: cmd += f' -silent_mut'
    if genome_wide_design: cmd += f' -genome_wide'
    if saturation_mutagenesis: cmd += f' -sat_mut {saturation_mutagenesis}'
    if number_of_pegrnas!=3: cmd += f' -n_pegrnas {number_of_pegrnas}'
    if number_of_ngrnas!=3: cmd += f' -n_ngrnas {number_of_ngrnas}'
    if nicking_distance_pooled!=75: cmd += f' -nick_dist_pooled {nicking_distance_pooled}'
    if homology_downstream!=10: cmd += f' -homology_downstream {homology_downstream}'
    if pbs_length_pooled!=14: cmd += f' -pbs_pooled {pbs_length_pooled}'
    if rtt_max_length_pooled!=50: cmd += f' -rtt_pooled {rtt_max_length_pooled}'
    if out_dir!='./DATETIMESTAMP_PrimeDesign': cmd+= f' -out {out_dir}'
    print(cmd)
    
    os.system(cmd) # Execute PrimeDesign Command Line

def prime_design_output(pt: str, scaffold_sequence: str, in_file: pd.DataFrame | str, saturation_mutagenesis:str=None, 
                        aa_index: int=1, enzymes: list[str]=['Esp3I'], replace: bool=True) -> dict[pd.DataFrame]:
    ''' 
    prime_design_output(): splits peg/ngRNAs from PrimeDesign output & finishes annotations
    
    Parameters:
    pt (str): path to primeDesign output
    scaffold_sequence (str): sgRNA scaffold sequence
        SpCas9 flip + extend (shorter): GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC
        SpCas9 flip + extend + com-modified (required for VLPs): GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
    in_file (Dataframe | str): input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (column names required)
    saturation_mutagenesis (str, optional): saturation mutagenesis design with prime editing (Options: 'aa', 'base').
    aa_index (int, optional): 1st amino acid in target sequence index (Default: 1)
    enzymes (list, optional): list of type IIS RE enzymes (i.e., Esp3I, BsaI, BspMI) to check for in pegRNAs and ngRNAs (Default: ['Esp3I'])
    replace (bool, optional): replace pegRNAs and remove ngRNAs with RE enzyme sites (Default: True)
    
    Dependencies: io, numpy, & pandas
    '''
    # Get target_name from input file
    if type(in_file) == str: # Get in_file from file path if needed
        in_file = io.get(pt=in_file)
    target_name_in_file = in_file.iloc[0]['target_name']

    if saturation_mutagenesis: # Saturation mutagenesis mode

        # Get PrimeDesign output & seperate pegRNAs and ngRNAs
        primeDesign_output = io.get(pt)
        pegRNAs = primeDesign_output[primeDesign_output['gRNA_type']=='pegRNA'].reset_index(drop=True)
        ngRNAs = primeDesign_output[primeDesign_output['gRNA_type']=='ngRNA'].reset_index(drop=True)

        # Generate pegRNAs
        pegRNAs['Edit']=[str(target_name.split('_')[-1].split('to')[0]) + # AA Before
                        str(int(target_name.split('_')[-2]) + aa_index-1) + # AA Index
                        str(target_name.split('_')[-1].split('to')[1]) # AA After
                        for target_name in pegRNAs['Target_name']]
        pegRNAs['Target_name']=[target_name_in_file]*len(pegRNAs)
        pegRNAs['Scaffold_sequence']=[scaffold_sequence]*len(pegRNAs)
        pegRNAs['RTT_sequence']=[pegRNAs.iloc[i]['Extension_sequence'][0:int(pegRNAs.iloc[i]['RTT_length'])] for i in range(len(pegRNAs))]
        pegRNAs['PBS_sequence']=[pegRNAs.iloc[i]['Extension_sequence'][int(pegRNAs.iloc[i]['RTT_length']):]  for i in range(len(pegRNAs))]
        pegRNAs['AA Number'] = [int(re.findall(r'-?\d+',edit)[0]) if edit is not None else aa_index for edit in pegRNAs['Edit']]
        pegRNAs = t.reorder_cols(df=pegRNAs,
                                 cols=['pegRNA_number','gRNA_type','Strand','Edit','AA Number', # Important metadata
                                    'Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence',  # Sequence information
                                    'Target_name','Target_sequence','Spacer_GC_content','PAM_sequence','Extension_sequence','Annotation','pegRNA-to-edit_distance','Nick_index','PBS_length','PBS_GC_content','RTT_length','RTT_GC_content','First_extension_nucleotide'], # Less important metadata
                                 keep=False) 
        
        # Generate ngRNAs
        ngRNAs['Edit']=[str(target_name.split('_')[-1].split('to')[0]) + # AA Before
                        str(int(target_name.split('_')[-2]) + aa_index-1) + # AA Index
                        str(target_name.split('_')[-1].split('to')[1]) # AA After
                        for target_name in ngRNAs['Target_name']]
        ngRNAs['Target_name']=[target_name_in_file]*len(ngRNAs)
        ngRNAs['Scaffold_sequence']=[scaffold_sequence]*len(ngRNAs)
        ngRNAs['ngRNA_number']=list(np.arange(1,len(ngRNAs)+1))
        ngRNAs['AA Number'] = [int(re.findall(r'-?\d+',edit)[0]) if edit is not None else aa_index for edit in ngRNAs['Edit']]
        ngRNAs = t.reorder_cols(df=ngRNAs,
                                cols=['pegRNA_number','ngRNA_number','gRNA_type','Strand','Edit','AA Number', # Important metadata
                                    'Spacer_sequence','Scaffold_sequence',  # Sequence information
                                    'Target_name','Target_sequence','Spacer_GC_content','PAM_sequence','Annotation','Nick_index','ngRNA-to-pegRNA_distance'], # Less important metadata
                                keep=False) 
    
    else: # Not saturation mutagenesis mode
        
        # Get PrimeDesign output & seperate pegRNAs and ngRNAs
        primeDesign_output = io.get(pt)
        pegRNAs = primeDesign_output[primeDesign_output['gRNA_type']=='pegRNA'].reset_index(drop=True)
        ngRNAs = primeDesign_output[primeDesign_output['gRNA_type']=='ngRNA'].reset_index(drop=True)

        # Generate pegRNAs
        pegRNAs['Scaffold_sequence']=[scaffold_sequence]*len(pegRNAs)
        pegRNAs['RTT_sequence']=[pegRNAs.iloc[i]['Extension_sequence'][0:int(pegRNAs.iloc[i]['RTT_length'])] for i in range(len(pegRNAs))]
        pegRNAs['PBS_sequence']=[pegRNAs.iloc[i]['Extension_sequence'][int(pegRNAs.iloc[i]['RTT_length']):]  for i in range(len(pegRNAs))]
        pegRNAs['Target_name']=[target_name_in_file]*len(pegRNAs)
        pegRNAs = t.reorder_cols(df=pegRNAs,
                                cols=['Target_name','pegRNA_number','gRNA_type','Strand', # Important metadata
                                    'Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence',  # Sequence information
                                    'Target_sequence','Spacer_GC_content','PAM_sequence','Extension_sequence','Annotation','pegRNA-to-edit_distance','Nick_index','PBS_length','PBS_GC_content','RTT_length','RTT_GC_content','First_extension_nucleotide'], # Less important metadata
                                keep=False) 
        
        # Generate ngRNAs
        ngRNAs['Scaffold_sequence']=[scaffold_sequence]*len(ngRNAs)
        ngRNAs['Target_name']=[target_name_in_file]*len(ngRNAs)
        ngRNAs = t.reorder_cols(df=ngRNAs,
                                cols=['Target_name','pegRNA_number','gRNA_type','Strand', # Important metadata
                                    'Spacer_sequence','Scaffold_sequence',  # Sequence information
                                    'Target_sequence','Spacer_GC_content','PAM_sequence','Extension_sequence','Annotation','Nick_index','ngRNA-to-pegRNA_distance'], # Less important metadata
                                keep=False)
    
    # Temporarily make pegRNAs and ngRNAs oligonucleotides
    pegRNAs['Oligonucleotide'] = [str(spacer+scaffold+rtt+pbs).upper() for (spacer, scaffold, rtt, pbs) in t.zip_cols(df=pegRNAs,cols=['Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence'])]
    ngRNAs['Oligonucleotide'] = [str(spacer+scaffold).upper() for (spacer, scaffold) in t.zip_cols(df=ngRNAs,cols=['Spacer_sequence','Scaffold_sequence'])]
    
    # Check for 0 recognition sites per enzyme
    for enzyme in enzymes:
        # pegRNAs: Find recognition sites for enzymes
        pegRNAs = find_enzyme_sites(df=pegRNAs, enzyme=enzyme)
        pegRNAs_edits = list(pegRNAs['Edit'].unique()) # Get pegRNA edits
    
        if replace: # Replace pegRNAs with RE enzyme sites
            
            # Store pegRNAs with recognition sites for enzymes
            pegRNAs_enzyme = pegRNAs[pegRNAs[enzyme]!=0]
            io.save(dir=f'../pegRNAs/{enzyme}/codon_swap_before',
                    file=f'{int(pegRNAs_enzyme.iloc[0]['PBS_length'])}.csv',
                    obj=pegRNAs_enzyme)
            
            # Codon swap pegRNAs with enzyme recognition site
            pegRNAs_enzyme = enzyme_codon_swap(pegRNAs=pegRNAs_enzyme,in_file=in_file,enzyme=enzyme)
            io.save(dir=f'../pegRNAs/{enzyme}/codon_swap_after',
                    file=f'{int(pegRNAs_enzyme.iloc[0]['PBS_length'])}.csv',
                    obj=pegRNAs_enzyme)
            pegRNAs = pd.concat([pegRNAs,pegRNAs_enzyme],ignore_index=True)
            print(f"pegRNAs edits recovered by modifying {enzyme} recognition site: {list(pegRNAs_enzyme['Edit'].unique())}")

            # Recheck pegRNAs for RE recognition sites and drop those with recognition sites
            pegRNAs = find_enzyme_sites(df=pegRNAs, enzyme=enzyme)
            pegRNAs = pegRNAs[pegRNAs[enzyme]==0].sort_values(by='pegRNA_number').reset_index(drop=True)

            # Store removed edits
            pegRNAs_enzyme = pegRNAs[pegRNAs[enzyme]!=0]
            remove_pegRNAs_edits = pegRNAs[pegRNAs['Edit'].isin(pegRNAs_enzyme['Edit'])]['Edit'].unique()
            
            # Save lost edits
            lost_pegRNAs_edits = [remove_edit for remove_edit in remove_pegRNAs_edits if remove_edit not in pegRNAs_edits]
            if len(lost_pegRNAs_edits) > 0:
                print(f"pegRNA edits lost due to {enzyme} recognition site: {lost_pegRNAs_edits}")
                io.save(dir=f'../pegRNAs/{enzyme}/lost',
                            file=f'{int(pegRNAs_enzyme.iloc[0]['PBS_length'])}.csv',
                            obj=pegRNAs_enzyme[pegRNAs_enzyme['Edit'].isin(lost_pegRNAs_edits)])

            # Drop enzyme column
            pegRNAs.drop(columns=[enzyme,f'{enzyme}_fwd_i',f'{enzyme}_rc_i'],inplace=True)

        # ngRNAs: Find recognition sites for enzymes
        ngRNAs = find_enzyme_sites(df=ngRNAs, enzyme=enzyme)

        if replace: # REMOVE ngRNAs with RE enzyme sites
            
            # Store ngRNAs with recognition sites for enzymes
            ngRNAs_enzyme = ngRNAs[ngRNAs[enzyme]!=0]
            io.save(dir=f'../ngRNAs/{enzyme}/codon_swap_before',
                    file=f'{int(pegRNAs.iloc[0]['PBS_length'])}.csv',
                    obj=ngRNAs_enzyme)

            # Drop ngRNAs with RE recognition sites
            ngRNAs = ngRNAs[ngRNAs[enzyme]==0].reset_index(drop=True)

            # Store removed edits
            ngRNAs_enzyme = ngRNAs[ngRNAs[enzyme]!=0]
            remove_ngRNAs_edits = ngRNAs[ngRNAs['Edit'].isin(ngRNAs_enzyme['Edit'])]['Edit'].unique()
            
            # Save lost edits
            lost_ngRNAs_edits = [remove_edit for remove_edit in remove_ngRNAs_edits if remove_edit not in ngRNAs['Edit'].unique()]
            if len(lost_ngRNAs_edits) > 0:
                print(f"ngRNA edits lost due to {enzyme} recognition site: {lost_ngRNAs_edits}")
                io.save(dir=f'../ngRNAs/{enzyme}/lost',
                        file=f'{int(pegRNAs.iloc[0]['PBS_length'])}.csv',
                        obj=ngRNAs_enzyme[ngRNAs_enzyme['Edit'].isin(lost_ngRNAs_edits)])

            # Drop enzyme column
            ngRNAs.drop(columns=[enzyme,f'{enzyme}_fwd_i',f'{enzyme}_rc_i'],inplace=True)
    
    # Remove oligonucleotide column
    pegRNAs.drop(columns=['Oligonucleotide'], inplace=True)
    ngRNAs.drop(columns=['Oligonucleotide'], inplace=True)

    return pegRNAs,ngRNAs

def prime_designer(target_name: str, flank5_sequence: str, target_sequence: str, flank3_sequence: str,
                  pbs_length_pooled_ls: list = [11,13,15], rtt_max_length_pooled: int = 50, silent_mutation: bool = True,
                  number_of_pegrnas: int = 1, number_of_ngrnas: int = 3,
                  scaffold_sequence: str='GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC', 
                  aa_index: int=1, enzymes: list[str]=['Esp3I'], replace: bool=True):
    '''
    prime_designer(): execute PrimeDesign saturation mutagenesis for EDMS
    
    Parameters:
    target_name (str): name of target
    flank5_sequence (str): in-frame nucleotide sequence with 5' of saturation mutagensis region (length must be divisible by 3)
    target_sequence (str): in-frame nucleotide sequence for the saturation mutagensis region (length must be divisible by 3)
    flank3_sequence (str): in-frame nucleotide sequence with 3' of saturation mutagensis region (length must be divisible by 3)
    pbs_length_pooled_ls (list, optional): list of primer binding site (PBS) lengths for the pegRNA extension (Default: [11,13,15])
    rtt_max_length_pooled (int, optional): maximum RTT length to design pegRNAs for pooled design applications. (Default: 50 nt)
    silent_mutation (bool, optional): introduce silent mutation into PAM assuming sequence is in-frame (Default: True)
    number_of_pegrnas (int, optional): maximum number of pegRNAs to design for each input sequence. The pegRNAs are ranked by 1) PAM disrupted > PAM intact then 2) distance to edit. (Default: 1)
    number_of_ngrnas (int, optional): maximum number of ngRNAs to design for each input sequence. The ngRNAs are ranked by 1) PE3b-seed > PE3b-nonseed > PE3 then 2) deviation from nicking_distance_pooled. (Default: 3)
    scaffold_sequence (str, optional): sgRNA scaffold sequence (Default: SpCas9 flip + extend = GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC)
        Alternative option for VLPs: SpCas9 flip + extend + com-modified = GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
    aa_index (int, optional): 1st amino acid in target sequence index (Default: 1)
    enzymes (list, optional): list of type IIS RE enzymes (i.e., Esp3I, BsaI, BspMI) to check for in pegRNAs and ngRNAs (Default: ['Esp3I'])
    replace (bool, optional): replace pegRNAs and remove ngRNAs with RE sites (Default: True)

    Dependencies: prime_design_input(), prime_design(), & prime_design_output()
    '''
    # Create PrimeDesign input file
    prime_design_input(target_name=target_name, flank5_sequence=flank5_sequence, target_sequence=target_sequence, 
                     flank3_sequence=flank3_sequence, aa_index=aa_index, dir='.', file=f'{"_".join(target_name.split(" "))}.csv')

    # Iterate through PBS lengths
    pegRNAs=dict()
    ngRNAs=dict()
    for pbs_length_pooled in pbs_length_pooled_ls:

        # Run PrimeDesign in saturation mutatgenesis mode
        prime_design(file=f'{"_".join(target_name.split(" "))}.csv', silent_mutation=silent_mutation, saturation_mutagenesis="aa",
                     number_of_pegrnas=number_of_pegrnas, number_of_ngrnas=number_of_ngrnas, pbs_length_pooled=pbs_length_pooled, rtt_max_length_pooled=rtt_max_length_pooled)
        
        # Obtain pegRNAs and ngRNAs from PrimeDesign output
        pegRNAs[pbs_length_pooled],ngRNAs[pbs_length_pooled] = prime_design_output(
            pt=sorted([file for file in io.relative_paths('.') if "PrimeDesign.csv" in file], reverse= True)[0], 
            scaffold_sequence=scaffold_sequence, in_file=f'./{"_".join(target_name.split(" "))}.csv', 
            saturation_mutagenesis='aa', aa_index=aa_index, enzymes=enzymes, replace=replace)
    
    # Save pegRNAs and ngRNAs
    io.save_dir(dir='../pegRNAs', suf='.csv', dc=pegRNAs)
    io.save_dir(dir='../ngRNAs', suf='.csv', dc=ngRNAs)

def merge(epegRNAs: str | dict | pd.DataFrame, ngRNAs: str | dict | pd.DataFrame, ngRNAs_groups_max: int=3,
          epegRNA_suffix: str='_epegRNA', ngRNA_suffix: str='_ngRNA', dir: str=None, file: str=None, literal_eval: bool=True) -> pd.DataFrame:
    '''
    merge(): rejoins epeg/ngRNAs & creates ngRNA_groups
    
    Parameters:
    epegRNAs (dict or dataframe): dictionary containing epegRNA dataframes or epegRNA dataframe
    ngRNAs (dict or dataframe): dictionary containing ngRNA dataframes or ngRNA dataframe
    ngRNAs_group_max (int, optional): maximum # of ngRNAs per epegRNA (Default: 3)
    epegRNA_suffix (str, optional): Suffix for epegRNAs columns (Default: epegRNA_)
    ngRNA_suffix (str, optional): Suffix for ngRNAs columns (Default: ngRNA_)
    literal_eval (bool, optional): convert string representations (Default: True)
    
    Dependencies: tidy & pandas
    '''
    # Get if epegRNAs and ngRNAs from path if needed
    if isinstance(epegRNAs, str): 
        if os.path.isdir(epegRNAs): # directory
            epegRNAs = io.get_dir(dir=epegRNAs, literal_eval=literal_eval)
        elif os.path.isfile(epegRNAs): # file
            epegRNAs = io.get(pt=epegRNAs, literal_eval=literal_eval)
        else:
            raise(ValueError(f"'epegRNAs' does not exist or is not a file/directory.\n{epegRNAs}"))
    if isinstance(ngRNAs, str): 
        if os.path.isdir(ngRNAs): # directory
            ngRNAs = io.get_dir(dir=ngRNAs, literal_eval=literal_eval)
        elif os.path.isfile(ngRNAs): # file
            ngRNAs = io.get(pt=ngRNAs, literal_eval=literal_eval)
        else:
            raise(ValueError(f"'ngRNAs' does not exist or is not a file/directory.\n{ngRNAs}"))

    # Join dictionary of dataframes if needed
    if isinstance(epegRNAs,dict): epegRNAs = t.join(epegRNAs).reset_index(drop=True)
    if isinstance(ngRNAs,dict): ngRNAs = t.join(ngRNAs).drop_duplicates(subset='ngRNA_number').reset_index(drop=True)

    # Limit to ngRNAs that correspond to epegRNAs
    ngRNAs = ngRNAs[[True if pegRNA_num in set(epegRNAs['pegRNA_number']) else False 
                     for pegRNA_num in ngRNAs['pegRNA_number']]].reset_index(drop=True)

    # Merge epegRNAs & ngRNAs
    epeg_ngRNAs = pd.merge(left=epegRNAs,
                           right=ngRNAs,
                           on='pegRNA_number',
                           suffixes=(epegRNA_suffix,ngRNA_suffix)).reset_index(drop=True)
    
    ngRNAs_dc = {(pegRNA_num):1 for (pegRNA_num) in list(epeg_ngRNAs['pegRNA_number'].value_counts().keys())}
    ngRNA_group_ls = []
    for pegRNA_num in epeg_ngRNAs['pegRNA_number']:
        ngRNA_group_ls.append(ngRNAs_dc[pegRNA_num]%ngRNAs_groups_max+1)
        ngRNAs_dc[pegRNA_num]+=1
    epeg_ngRNAs['ngRNA_group']=ngRNA_group_ls
    
    # Save epeg_ngRNAs if dir and file are provided
    if dir is not None and file is not None:
        io.save(dir=dir, file=file, obj=epeg_ngRNAs)

    return epeg_ngRNAs

# pegRNA
def epegRNA_linkers(pegRNAs: str | pd.DataFrame, epegRNA_motif_sequence: str='CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA',
                    linker_pattern: str='NNNNNNNN', excluded_motifs: list=['Esp3I'],
                    ckpt_dir: str=None, ckpt_file=None, ckpt_pt: str='',
                    out_dir: str=None, out_file: str=None, literal_eval: bool=True) -> pd.DataFrame:
    ''' 
    epegRNA_linkers(): generate epegRNA linkers between PBS and 3' hairpin motif & finish annotations
    
    Parameters:
    pegRNAs (str | dataframe): pegRNAs DataFrame or file path
    epegRNA_motif_sequence (str, optional): epegRNA motif sequence (Optional, Default: tevopreQ1)
    linker_pattern (str, optional): epegRNA linker pattern (Default: NNNNNNNN)
    excluded_motifs (list, optional): list of motifs or type IIS RE enzymes (i.e., Esp3I, BsaI, BspMI) to exclude from linker generation (Default: ['Esp3I'])
    ckpt_dir (str, optional): Checkpoint directory
    ckpt_file (str, optional): Checkpoint file name
    ckpt_pt (str, optional): Previous ckpt path
    literal_eval (bool, optional): convert string representations (Default: True)
    
    Dependencies: pandas, pegLIT, & io
    '''
    if type(pegRNAs)==str: # Get pegRNAs dataframe from file path if needed
        pegRNAs = io.get(pt=pegRNAs,literal_eval=literal_eval)

    # Parse excluded_motifs
    if excluded_motifs is not None: # Check if excluded_motifs is a list 
        RE_type_IIS_df = load_resource_csv(filename='RE_type_IIS.csv')
        for motif in excluded_motifs: # Find type IIS RE and replace with recognition sequence (+ reverse complement)
            if motif in list(RE_type_IIS_df['Name']):
                excluded_motifs.remove(motif)
                excluded_motifs.append(RE_type_IIS_df[RE_type_IIS_df['Name']==motif]['Recognition'].values[0])
                excluded_motifs.append(RE_type_IIS_df[RE_type_IIS_df['Name']==motif]['Recognition_rc'].values[0])
    
    # Get or make ckpt DataFrame & linkers
    linkers = []
    if ckpt_dir is not None and ckpt_file is not None: # Save ckpts
        if ckpt_pt=='': 
            ckpt = pd.DataFrame(columns=['pegRNA_number','Linker_sequence'])
        else: 
            ckpt = io.get(pt=ckpt_pt)
            linkers = list(ckpt['Linker_sequence']) # Get linkers from ckpt
    else: 
        ckpt = '' # Don't save ckpts, length needs to 0.

    # Generate epegRNA linkers between PBS and 3' hairpin motif
    for i in range(len(pegRNAs)):
        if i>=len(ckpt):
            linkers.extend(pegLIT.pegLIT(seq_spacer=pegRNAs.iloc[i]['Spacer_sequence'],seq_scaffold=pegRNAs.iloc[i]['Scaffold_sequence'],
                                         seq_template=pegRNAs.iloc[i]['RTT_sequence'],seq_pbs=pegRNAs.iloc[i]['PBS_sequence'],
                                         seq_motif=epegRNA_motif_sequence,linker_pattern=linker_pattern,excluded_motifs=excluded_motifs))
            if ckpt_dir is not None and ckpt_file is not None: # Save ckpts
                ckpt = pd.concat([ckpt,pd.DataFrame({'pegRNA_number': [i], 'Linker_sequence': [linkers[i]]})])
                io.save(dir=ckpt_dir,file=ckpt_file,obj=ckpt)
            print(f'Status: {i} out of {len(pegRNAs)}')
    
    # Generate epegRNAs
    pegRNAs['Linker_sequence'] = linkers
    pegRNAs['Motif_sequence'] = [epegRNA_motif_sequence]*len(pegRNAs)
    epegRNAs = t.reorder_cols(df=pegRNAs,
                              cols=['pegRNA_number','gRNA_type','Strand','Edit', # Important metadata
                                    'Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence','Linker_sequence','Motif_sequence']) # Sequence information
    
    # Save epeg_ngRNAs if dir and file are provided
    if out_dir is not None and out_file is not None:
        io.save(dir=out_dir, file=out_file, obj=epegRNAs)
    
    return epegRNAs

def shared_sequences(pegRNAs: pd.DataFrame | str, hist_plot:bool=True, hist_dir: str=None, hist_file: str=None, literal_eval: bool=True, **kwargs) -> pd.DataFrame:
    ''' 
    shared_sequences(): Reduce PE library into shared spacers and PBS sequences
    
    Parameters:
    pegRNAs (dataframe | str): pegRNAs DataFrame (or file path)
    hist_plot (bool, optional): display histogram of reduced PE library (Default: True)
    hist_dir (str, optional): directory to save histogram
    hist_file (str, optional): file name to save histogram
    literal_eval (bool, optional): convert string representations (Default: True)

    Dependencies: pandas & plot
    '''
    # Get pegRNAs DataFrame from file path if needed
    if type(pegRNAs)==str:
        pegRNAs = io.get(pt=pegRNAs, literal_eval=literal_eval)

    # Reduce PE library to the set shared of spacers and PBS motifs
    shared = sorted({(pegRNAs.iloc[i]['Spacer_sequence'],pegRNAs.iloc[i]['PBS_sequence']) for i in range(len(pegRNAs))})
    shared_pegRNAs_lib = pd.DataFrame(columns=['Target_name','pegRNA_numbers','Strand','Edits','Spacer_sequence','PBS_sequence'])
    for (spacer,pbs) in shared:
        shared_pegRNAs = pegRNAs[(pegRNAs['Spacer_sequence']==spacer)&(pegRNAs['PBS_sequence']==pbs)]
        shared_pegRNAs_lib = pd.concat([shared_pegRNAs_lib,
                                        pd.DataFrame({'Target_name': [shared_pegRNAs.iloc[0]['Target_name']],
                                                      'pegRNA_numbers': [shared_pegRNAs['pegRNA_number'].to_list()],
                                                      'Strand': [shared_pegRNAs.iloc[0]['Strand']],
                                                      'Edits': [shared_pegRNAs['Edit'].to_list()],
                                                      'Spacer_sequence': [spacer],
                                                      'PBS_sequence': [pbs],
                                                      'RTT_lengths': [sorted(int(rtt) for rtt in set(shared_pegRNAs['RTT_length'].to_list()))]})]).reset_index(drop=True)
    
    # Find shared AAs within the reduced PE library
    aa_numbers_ls=[]
    aa_numbers_min_ls=[]
    aa_numbers_max_ls=[]
    continous_ls=[]
    for edits in shared_pegRNAs_lib['Edits']:
        aa_numbers = {int(edit[1:-1]) for edit in edits}
        aa_numbers_min = min(aa_numbers)
        aa_numbers_max = max(aa_numbers)
        if aa_numbers == set(range(aa_numbers_min,aa_numbers_max+1)): continous=True
        else: continous=False
        aa_numbers_ls.append(sorted(aa_numbers))
        aa_numbers_min_ls.append(aa_numbers_min)
        aa_numbers_max_ls.append(aa_numbers_max)
        continous_ls.append(continous)
    shared_pegRNAs_lib['AA_numbers']=aa_numbers_ls
    shared_pegRNAs_lib['AA_numbers_min']=aa_numbers_min_ls
    shared_pegRNAs_lib['AA_numbers_max']=aa_numbers_max_ls
    shared_pegRNAs_lib['AA_numbers_continuous']=continous_ls
    shared_pegRNAs_lib = shared_pegRNAs_lib.sort_values(by=['AA_numbers_min','AA_numbers_max']).reset_index(drop=True)

    if hist_plot: # Generate histogram
        shared_hist = pd.DataFrame()
        for i,aa_numbers in enumerate(shared_pegRNAs_lib['AA_numbers']):
            shared_hist = pd.concat([shared_hist,pd.DataFrame({'Group_Spacer_PBS': [f'{str(i)}_{shared_pegRNAs_lib.iloc[i]["Spacer_sequence"]}_{shared_pegRNAs_lib.iloc[i]["PBS_sequence"]}']*len(aa_numbers),
                                                               'AA_number': aa_numbers})]).reset_index(drop=True)
        p.dist(typ='hist',df=shared_hist,x='AA_number',cols='Group_Spacer_PBS',x_axis='AA number',title=f'Shared Spacers & PBS Sequences in the {shared_pegRNAs_lib.iloc[0]['Target_name']} PE Library',
               x_axis_dims=(min(shared_hist['AA_number']),max(shared_hist['AA_number'])),figsize=(10,2),bins=max(shared_hist['AA_number'])-min(shared_hist['AA_number'])+1,
               legend_loc='upper center',legend_bbox_to_anchor=(0.5, -.3),dir=hist_dir,file=hist_file,legend_ncol=2,**kwargs)

    return shared_pegRNAs_lib

def pilot_screen(pegRNAs_dir: str, mutations_pt: str, database: Literal['COSMIC','ClinVar']='COSMIC', literal_eval: bool=True):
    ''' 
    pilot_screen(): Determine pilot screen for EDMS
    
    Parameters:
    pegRNAs_dir (str): directory with pegRNAs from prime_designer() output
    mutations_pt (str): path to mutations file (COSMIC or ClinVar)
    database (str, optional): database to use for priority mutations (Default: 'COSMIC')
    literal_eval (bool, optional): convert string representations (Default: True)
    
    Dependencies: io, cosmic, cvar, shared_sequences(), priority_muts(), & priority_edits()
    '''
    # Get pegRNAs from prime_designer() output
    pegRNAs = io.get_dir(pegRNAs_dir,literal_eval=literal_eval)

    # Get mutations from COSMIC or ClinVar file
    if database=='COSMIC':
        mutations = co.mutations(io.get(pt=mutations_pt,literal_eval=literal_eval))
    elif database=='ClinVar':
        mutations = cvar.mutations(io.get(pt=mutations_pt,literal_eval=literal_eval))
    else: raise(ValueError(f"Database {database} not supported. Use 'COSMIC' or 'ClinVar'."))

    # Isolate shared spacer & PBS sequences
    pegRNAs_shared = dict()
    for key,pegRNAs_pbs in pegRNAs.items():
        pegRNAs_shared[key] = shared_sequences(pegRNAs=pegRNAs_pbs,
                                               hist_dir='../shared_sequences',
                                               hist_file=f'{key}.png',
                                               show=False)
    io.save_dir(dir='../shared_sequences',
                suf='.csv',
                dc=pegRNAs_shared)

    # Determine priority mutations for each shared spacer & PBS sequence
    pegRNAs_shared_muts = dict()
    for key,pegRNAs_shared_pbs in pegRNAs_shared.items():
        if database=='COSMIC': pegRNAs_shared_muts[key]=co.priority_muts(pegRNAs_shared=pegRNAs_shared_pbs,
                                                                         df_cosmic=mutations)
        elif database=='ClinVar': pegRNAs_shared_muts[key]=cvar.priority_muts(pegRNAs_shared=pegRNAs_shared_pbs,
                                                                              df_clinvar=mutations)
        else: raise(ValueError(f"Database {database} not supported. Use 'COSMIC' or 'ClinVar'."))

    io.save_dir(dir='../shared_sequences_muts',
                suf='.csv',
                dc=pegRNAs_shared_muts)

    # Determine priority edits for each shared spacer & PBS sequence
    pegRNAs_priority = dict()
    for key,pegRNAs_shared_pbs in pegRNAs_shared_muts.items():
        if database=='COSMIC': pegRNAs_priority[key]=co.priority_edits(pegRNAs=pegRNAs[key],
                                                                       pegRNAs_shared=pegRNAs_shared_pbs,
                                                                       df_cosmic=mutations)
        elif database=='ClinVar':  pegRNAs_shared_muts[key]=cvar.priority_edits(pegRNAs=pegRNAs[key],
                                                                                pegRNAs_shared=pegRNAs_shared_pbs,
                                                                                df_clinvar=mutations)
        else: raise(ValueError(f"Database {database} not supported. Use 'COSMIC' or 'ClinVar'."))
    
    io.save_dir(dir='../pegRNAs_priority',
                suf='.csv',
                dc=pegRNAs_priority)

def sensor_designer(pegRNAs: pd.DataFrame | str, in_file: str, sensor_length: int=60, before_spacer: int=5, sensor_orientation: Literal['revcom','forward']='revcom',
                    out_dir: str=None, out_file: str=None, return_df: bool=True, literal_eval: bool=True) -> pd.DataFrame:
    ''' 
    sensor_designer(): design pegRNA sensors
    
    Parameters:
    pegRNAs (dataframe | str): pegRNAs DataFrame (or file path)
    in_file (dataframe | str): Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (column names required)
    sensor_length (int, optional): Total length of the sensor in bp (Default = 60)
    before_spacer (int, optional): Amount of nucleotide context to put before the protospacer in the sensor (Default = 5)
    sensor_orientation (Literal, optional): Orientation of the sensor relative to the protospacer (Options: 'revcom' [Default b/c minimize recombination] or ’forward’).
    out_dir (str, optional): output directory
    out_file (str, optional): output filename
    return_df (bool, optional): return dataframe (Default: True)
    literal_eval (bool, optional): convert string representations (Default: True)

    Dependencies: io, pandas, Bio.Seq.Seq
    '''
    # Initialize timer; memory reporting
    memory_timer(reset=True)
    memories = []

    # Get pegRNAs & PrimeDesign input DataFrame from file path if needed
    if type(pegRNAs)==str:
        pegRNAs = io.get(pt=pegRNAs, literal_eval=literal_eval)
    if type(in_file)==str:
        in_file = io.get(pt=in_file, literal_eval=literal_eval)

    # Get reference sequence & codons (+ reverse complement)
    target_sequence = in_file.iloc[0]['target_sequence']
    seq = Seq(target_sequence.split('(')[1].split(')')[0]) # Break apart target sequences
    flank5 = Seq(target_sequence.split('(')[0])
    flank3 = Seq(target_sequence.split(')')[1])
    seq_nuc = flank5 + seq + flank3  # Join full nucleotide reference sequence
    rc_seq_nuc = Seq.reverse_complement(seq_nuc) # Full nucleotide reference reverse complement sequence
    
    # Check sensor_length
    if sensor_length <= 0:
        raise ValueError(f"Sensor length <= {sensor_length}")
    elif sensor_length % 2 != 0:
        print("Warning sensor length was not an even integer. Added 1.")
        sensor_length += 1
    
    # Check before_spacer
    if before_spacer <= 0:
        raise ValueError(f"Before spacer length <= {before_spacer}")
    
    # Find sensors
    sensors = []
    for spacer,strand in t.zip_cols(df=pegRNAs,cols=['Spacer_sequence','Strand']): # Iterate through spacers

        if strand=='+': # Spacer: + strand; PBS & RTT: - strand
            
            # Find spacer in sequence; compute spacer5 index
            spacer5 = seq_nuc.find(spacer)
            if spacer5 == -1:
                raise ValueError(f"Spacer sequence '{spacer}' not found in target sequence. Please check the input file.")
            elif spacer5 != seq_nuc.rfind(spacer):
                raise ValueError(f"Multiple matches found for spacer sequence '{spacer}'. Please check the input file.")

            # Assign start & end index for sensor
            start = spacer5 - before_spacer
            end = start + sensor_length
            sensor = seq_nuc[start:end]

        elif strand=='-': # Spacer: - strand; PBS & RTT: + strand
            
            # Find spacer in sequence; compute center index
            spacer5 = rc_seq_nuc.find(spacer)
            if spacer5 == -1:
                raise ValueError(f"Spacer sequence '{spacer}' not found in target sequence. Please check the input file.")
            if spacer5 != rc_seq_nuc.rfind(spacer):
                raise ValueError(f"Multiple matches found for spacer sequence '{spacer}' not found in target sequence. Please check the input file.")
            
            # Assign start & end index for sensor
            start = spacer5 - before_spacer
            end = start + sensor_length
            sensor = rc_seq_nuc[start:end]
            
        # Append sensor to list (revcom if specified)
        if sensor_orientation=='revcom':
            sensors.append(str(Seq.reverse_complement(Seq(sensor))))
        elif sensor_orientation=='forward':
            sensors.append(sensor)
        else:
            raise(ValueError(f"sensor_orientation = {sensor_orientation} was not 'revcom' or 'forward'."))

    # Add to dataframe
    pegRNAs['Sensor_sequence'] = sensors
    
    # Save & Return
    memories.append(memory_timer(task=f"sensors()"))
    if out_dir is not None and out_file is not None:
        io.save(dir=os.path.join(out_dir,f'.sensors'),
                file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_memories.csv',
                obj=pd.DataFrame(memories, columns=['Task','Memory, MB','Time, s']))
        io.save(dir=out_dir,file=out_file,obj=pegRNAs)
    if return_df==True: return pegRNAs

def rtt_designer(pegRNAs: pd.DataFrame | str, in_file: pd.DataFrame | str, rtt_length: int=39, 
                 include_WT: bool=False, enzymes: list[str]=['Esp3I'], replace: bool=True,
                 out_dir: str=None, out_file: str=None, return_df: bool=True, literal_eval: bool=True, comments: bool=False) -> pd.DataFrame:
    ''' 
    rtt_designer(): design all possible RTT for given spacer & PBS (WT, single insertions, & single deletions)
    
    Parameters:
    pegRNAs (dataframe | str): pegRNAs DataFrame (or file path)
    in_file (dataframe | str): Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (column names required)
    RTT_length (int, optional): Reverse transcriptase template length (bp)
    include_WT (bool, optional): include wildtype RTT (Default: False)
    enzymes (list, optional): list of type IIS RE enzymes (i.e., Esp3I, BsaI, BspMI) to check for in pegRNAs (Default: ['Esp3I'])
    replace (bool, optional): replace pegRNAs with RE enzyme sites (Default: True)
    out_dir (str, optional): output directory
    out_file (str, optional): output filename
    return_df (bool, optional): return dataframe (Default: True)
    literal_eval (bool, optional): convert string representations (Default: True)
    comments (bool, optional): print comments (Default: False)

    Dependencies: io, pandas, Bio.Seq.Seq, shared_sequences(), get_codons(), get_codon_frames(), found_list_in_order(), & aa_dna_codon_table
    '''
    # Initialize timer; memory reporting
    memory_timer(reset=True)
    memories = []

    # Get pegRNAs & PrimeDesign input DataFrame from file path if needed
    if type(pegRNAs)==str:
        pegRNAs = io.get(pt=pegRNAs,literal_eval=literal_eval)
    if type(in_file)==str:
        in_file = io.get(pt=in_file,literal_eval=literal_eval)

    # Get reference sequence & codons (+ reverse complement)
    target_sequence = in_file.iloc[0]['target_sequence'] 
    seq = Seq(target_sequence.split('(')[1].split(')')[0]) # Break apart target sequences
    if len(seq)%3 != 0: raise(ValueError(f"Length of target sequence ({len(seq)}) must divisible by 3. Check input file."))
    flank5 = Seq(target_sequence.split('(')[0])
    if len(flank5)%3 != 0: raise(ValueError(f"Length of 5' flank ({len(flank5)}) must divisible by 3. Check input file."))
    flank3 = Seq(target_sequence.split(')')[1])
    if len(flank3)%3 != 0: raise(ValueError(f"Length of 3' flank ({len(flank3)}) must divisible by 3. Check input file."))
    target_name = in_file.iloc[0]['target_name'] # Get target name

    aa_index = in_file.iloc[0]['aa_index']

    f5_seq_f3_nuc = flank5 + seq + flank3  # Join full nucleotide reference sequence
    rc_f5_seq_f3_nuc = Seq.reverse_complement(f5_seq_f3_nuc) # Full nucleotide reference reverse complement sequence
    seq_prot = Seq.translate(seq) # In-frame amino acid sequence
    f5_seq_f3_prot = Seq.translate(f5_seq_f3_nuc) # Full in-frame protein sequence (including flanks)
    
    codons = get_codons(seq) # Codons
    codons_flank5 = get_codons(flank5) # Codons in-frame flank 5
    codons_flank3 = get_codons(flank3) # Codons in-frame flank 3
    extended_codons = codons_flank5 + codons + codons_flank3 # Codons including flank 5 and flank 3
    extended_codons_nuc = Seq('').join(extended_codons) # Join codons into full in-frame nucleotide sequence
    extended_codons_prot = Seq.translate(extended_codons_nuc) # Translate to full in-frame protein sequence
    extended_codons_aa_indexes = list(np.arange(aa_index-len(codons_flank5),aa_index-len(codons_flank5)+len(extended_codons_prot))) # Obtain full in-frame amino acid indexes

    print(f'FWD Ref: {f5_seq_f3_nuc}')
    print(f'REV Ref: {rc_f5_seq_f3_nuc}')
    print(f'Nucleotides: {seq}')
    print(f'Amino Acids: {seq_prot}\n')

    # Obtain shared spacer and PBS sequences 
    shared_pegRNAs_lib = shared_sequences(pegRNAs=pegRNAs,hist_plot=False)

    # Obtain WT RTT, single insertions, and single deletions
    if include_WT==True:
        wildtypes = pd.DataFrame()
    insertions = pd.DataFrame()
    deletions = pd.DataFrame()
    for j,(spacer,pbs,strand) in enumerate(t.zip_cols(df=shared_pegRNAs_lib,cols=['Spacer_sequence','PBS_sequence','Strand'])): # Iterate through primer binding sites

        found = False # Boolean for RTT wildtype found
        if strand=='+': # Spacer: + strand; PBS & RTT: - strand
            
            # Find spacer in sequence
            spacer_j = f5_seq_f3_nuc.find(spacer)
            if spacer_j == -1:
                raise ValueError(f"Spacer sequence '{spacer}' not found in target sequence. Please check the input file.")
            elif spacer_j != f5_seq_f3_nuc.rfind(spacer):
                raise ValueError(f"Multiple matches found for spacer sequence '{spacer}'. Please check the input file.")

            # Find PBS in reverse complement sequence
            pbs_j = rc_f5_seq_f3_nuc.find(pbs)
            if pbs_j == -1:
                raise ValueError(f"PBS sequence '{pbs}' not found in target sequence. Please check the input file.")
            elif pbs_j != rc_f5_seq_f3_nuc.rfind(pbs):
                print(pbs,pbs_j,rc_f5_seq_f3_nuc.rfind(pbs))
                raise ValueError(f"Multiple matches found for PBS sequence '{pbs}'. Please check the input file.")

            # Obtain WT RTT from - strand
            rtt_wt = rc_f5_seq_f3_nuc[pbs_j-rtt_length:pbs_j]
            if include_WT==True:
                wildtypes = pd.concat([wildtypes,
                                    pd.DataFrame({'pegRNA_number': [j],
                                                    'gRNA_type': ['pegRNA'],
                                                    'Strand': [strand],
                                                    'Edit': [None],
                                                    'Spacer_sequence': [spacer],
                                                    'Scaffold_sequence': [pegRNAs.iloc[0]['Scaffold_sequence']],
                                                    'RTT_sequence': [str(rtt_wt)],
                                                    'PBS_sequence': [pbs],
                                                    'Target_name': [target_name],
                                                    'Target_sequence': [None],
                                                    'Spacer_GC_content': [None], 
                                                    'PAM_sequence': [None],
                                                    'Extension_sequence': [''.join([str(rtt_wt),pbs])], 
                                                    'Annotation': ['wildtype'], 
                                                    'pegRNA-to-edit_distance': [None],
                                                    'Nick_index': [None],
                                                    'ngRNA-to-pegRNA_distance': [None], 
                                                    'PBS_length': [len(pbs)],
                                                    'PBS_GC_content': [None], 
                                                    'RTT_length': [rtt_length], 
                                                    'RTT_GC_content': [None],
                                                    'First_extension_nucleotide': [rtt_wt[0]]})]).reset_index(drop=True)
            
            # Obtain reverse complement WT RTT in-frame from + strand
            rc_rtt_wt = Seq.reverse_complement(rtt_wt) # reverse complement of rtt (+ strand)
            rc_rtt_codon_frames = get_codon_frames(rc_rtt_wt) # codons
            if comments==True:
                print(f"Extended Codons (Here): {extended_codons[math.floor((len(rc_f5_seq_f3_nuc)-pbs_j)/3)-1:math.floor((len(rc_f5_seq_f3_nuc)-pbs_j+rtt_length)/3)]}")
            for i,rc_rtt_codon_frame in enumerate(rc_rtt_codon_frames): # Search for in-frame nucleotide sequence
                if comments==True:
                    print(f"rc_rtt_codon_frame: {rc_rtt_codon_frame}")
                
                index = found_list_in_order(extended_codons[math.floor((len(rc_f5_seq_f3_nuc)-pbs_j)/3)-1:math.floor((len(rc_f5_seq_f3_nuc)-pbs_j+rtt_length)/3)],rc_rtt_codon_frame)
                if index != -1: # Codon frame from reverse complement of rtt matches extended codons of in-frame nucleotide sequence
                    rc_rtt_wt_inframe_nuc_codons_flank5 = rc_rtt_wt[:i] # Save codon frame flank 5'
                    rc_rtt_wt_inframe_nuc_codons = rc_rtt_codon_frame # Save codon frame
                    rc_rtt_wt_inframe_nuc_codons_flank3 = rc_rtt_wt[i+3*len(rc_rtt_codon_frame):] # Save codon frame flank 3'
                    rc_rtt_wt_inframe_nuc = Seq('').join(rc_rtt_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rc_rtt_wt_inframe_prot = Seq.translate(rc_rtt_wt_inframe_nuc) # Translate to in-frame protein sequence
                    rc_rtt_wt_inframe_prot_indexes = extended_codons_aa_indexes[extended_codons_prot.find(rc_rtt_wt_inframe_prot):extended_codons_prot.find(rc_rtt_wt_inframe_prot)+len(rc_rtt_wt_inframe_prot)] # Obtain correponding aa indexes
                    rc_rtt_wt_inframe_prot_deletions = f5_seq_f3_prot[f5_seq_f3_prot.find(rc_rtt_wt_inframe_prot):f5_seq_f3_prot.find(rc_rtt_wt_inframe_prot)+len(rc_rtt_wt_inframe_prot)+1] # Store AAs for deletion names
                    found=True
                    break
            
            if comments==True:
                print(f'Strand: {strand}')    
                print(f'Nucleotides: {rc_rtt_wt}')
                print(f'Nucleotides 5\' of Codons: {rc_rtt_wt_inframe_nuc_codons_flank5}')
                print(f'Nucleotides Codons: {rc_rtt_wt_inframe_nuc_codons}')
                print(f'Nucleotides 3\' of Codons: {rc_rtt_wt_inframe_nuc_codons_flank3}')
                print(f'Nucleotides In-Frame: {rc_rtt_wt_inframe_nuc}')
                print(f'Amino Acids In-Frame: {rc_rtt_wt_inframe_prot}')
                print(f'Amino Acid #s In-Frame: {rc_rtt_wt_inframe_prot_indexes}\n')

            if found==False:
                raise(ValueError("RTT was not found."))

            # Obtain single insertion RTTs from - strand
            edits_in = []
            rtts_in = []
            for i in range(len(rc_rtt_wt_inframe_nuc_codons)): # Iterate through all in-frame codon positions
                for codon_table_aa,codon_table_dna in aa_dna_codon_table.items(): # Obtain all possible codon insertions
                    if codon_table_aa!='*': # Remove stop codons
                        edits_in.append(f'{rc_rtt_wt_inframe_prot[i]}{rc_rtt_wt_inframe_prot_indexes[i]}{rc_rtt_wt_inframe_prot[i]}{codon_table_aa}')
                        rtts_in.append(Seq.reverse_complement(Seq('').join([rc_rtt_wt_inframe_nuc_codons_flank5, # Codon frame flank 5'
                                                              Seq('').join(rc_rtt_wt_inframe_nuc_codons[:i+1]), # Codons before insertion
                                                              Seq(codon_table_dna[0]).lower(), # Insertion codon
                                                              Seq('').join(rc_rtt_wt_inframe_nuc_codons[i+1:]), # Codons after insertion
                                                              rc_rtt_wt_inframe_nuc_codons_flank3]))) # Codon frame flank 3'
            
            if comments==True:
                print(f'Insertions: {edits_in}')
                print(f'Insertion RTTs: {rtts_in}\n')

            insertions = pd.concat([insertions,
                                    pd.DataFrame({'pegRNA_number': [j]*len(edits_in),
                                                  'gRNA_type': ['pegRNA']*len(edits_in),
                                                  'Strand': [strand]*len(edits_in),
                                                  'Edit': edits_in,
                                                  'Spacer_sequence': [spacer]*len(edits_in),
                                                  'Scaffold_sequence': [pegRNAs.iloc[0]['Scaffold_sequence']]*len(edits_in),
                                                  'RTT_sequence': [str(rtt_in) for rtt_in in rtts_in],
                                                  'PBS_sequence': [pbs]*len(edits_in),
                                                  'Target_name': [target_name]*len(edits_in),
                                                  'Target_sequence': [None]*len(edits_in),
                                                  'Spacer_GC_content': [None]*len(edits_in),
                                                  'PAM_sequence': [None]*len(edits_in),
                                                  'Extension_sequence': [''.join([str(rtt_in),pbs]) for rtt_in in rtts_in], 
                                                  'Annotation': ['insertion']*len(edits_in),
                                                  'pegRNA-to-edit_distance': [None]*len(edits_in),
                                                  'Nick_index': [None]*len(edits_in),
                                                  'ngRNA-to-pegRNA_distance': [None]*len(edits_in),
                                                  'PBS_length': [len(pbs)]*len(edits_in),
                                                  'PBS_GC_content': [None]*len(edits_in),
                                                  'RTT_length': [len(rtt_in) for rtt_in in rtts_in], 
                                                  'RTT_GC_content': [None]*len(edits_in),
                                                  'First_extension_nucleotide': [rtt_in[0] for rtt_in in rtts_in]})]).reset_index(drop=True)

            # Obtain single deletion RTTs from - strand
            edits_del = [f'{aa}{rc_rtt_wt_inframe_prot_deletions[i+1]}{rc_rtt_wt_inframe_prot_indexes[i]}{rc_rtt_wt_inframe_prot_deletions[i+1]}' for i,aa in enumerate(rc_rtt_wt_inframe_prot) if i<len(rc_rtt_wt_inframe_prot)-1] # Don't want last AA
            rtts_del = [Seq.reverse_complement(Seq('').join([rc_rtt_wt_inframe_nuc_codons_flank5, # Codon frame flank 5'
                                               Seq('').join(rc_rtt_wt_inframe_nuc_codons[:i]), # Codons before deletion
                                               Seq('').join(rc_rtt_wt_inframe_nuc_codons[i+1:]), # Codons after deletion
                                               rc_rtt_wt_inframe_nuc_codons_flank3])) # Codon frame flank 3'
                                               for i in range(len(rc_rtt_wt_inframe_nuc_codons)) if i<len(rc_rtt_wt_inframe_nuc_codons)-1] # Don't want last AA
            
            if comments==True:
                print(f'Deletions: {edits_del}')
                print(f'Deletion RTTs: {rtts_del}\n\n')

            deletions = pd.concat([deletions,
                                   pd.DataFrame({'pegRNA_number': [j]*len(edits_del),
                                                 'gRNA_type': ['pegRNA']*len(edits_del),
                                                 'Strand': [strand]*len(edits_del),
                                                 'Edit': edits_del,
                                                 'Spacer_sequence': [spacer]*len(edits_del),
                                                 'Scaffold_sequence': [pegRNAs.iloc[0]['Scaffold_sequence']]*len(edits_del),
                                                 'RTT_sequence': [str(rtt_del) for rtt_del in rtts_del],
                                                 'PBS_sequence': [pbs]*len(edits_del),
                                                 'Target_name': [target_name]*len(edits_del),
                                                 'Target_sequence': [None]*len(edits_del),
                                                 'Spacer_GC_content': [None]*len(edits_del), 
                                                 'PAM_sequence': [None]*len(edits_del),
                                                 'Extension_sequence': [''.join([str(rtt_del),pbs]) for rtt_del in rtts_del], 
                                                 'Annotation': ['deletion']*len(edits_del), 
                                                 'pegRNA-to-edit_distance': [None]*len(edits_del),
                                                 'Nick_index': [None]*len(edits_del),
                                                 'ngRNA-to-pegRNA_distance': [None]*len(edits_del), 
                                                 'PBS_length': [len(pbs)]*len(edits_del),
                                                 'PBS_GC_content': [None]*len(edits_del),
                                                 'RTT_length': [len(rtt_del) for rtt_del in rtts_del], 
                                                 'RTT_GC_content': [None]*len(edits_del),
                                                 'First_extension_nucleotide': [rtt_del[0] for rtt_del in rtts_del]})]).reset_index(drop=True)
            
        elif strand=='-': # Spacer: - strand; PBS & RTT: + strand
            
            # Find spacer in sequence
            spacer_j = rc_f5_seq_f3_nuc.find(spacer)
            if spacer_j == -1:
                raise ValueError(f"Spacer sequence '{spacer}' not found in target sequence. Please check the input file.")
            if spacer_j != rc_f5_seq_f3_nuc.rfind(spacer):
                raise ValueError(f"Multiple matches found for spacer sequence '{spacer}' not found in target sequence. Please check the input file.")

            # Find PBS in sequence
            pbs_j = f5_seq_f3_nuc.find(pbs)
            if pbs_j == -1:
                raise ValueError(f"PBS sequence '{pbs}' not found in target sequence. Please check the input file.")
            if pbs_j != f5_seq_f3_nuc.rfind(pbs):
                print(pbs,pbs_j,f5_seq_f3_nuc.rfind(pbs))
                raise ValueError(f"Multiple matches found for PBS sequence '{pbs}' not found in target sequence. Please check the input file.")
            
            # Obtain WT RTT from + strand
            rtt_wt = f5_seq_f3_nuc[pbs_j-rtt_length:pbs_j]
            if include_WT==True:
                wildtypes = pd.concat([wildtypes,
                                    pd.DataFrame({'pegRNA_number': [j],
                                                    'gRNA_type': ['pegRNA'],
                                                    'Strand': [strand],
                                                    'Edit': [None],
                                                    'Spacer_sequence': [spacer],
                                                    'Scaffold_sequence': [pegRNAs.iloc[0]['Scaffold_sequence']],
                                                    'RTT_sequence': [str(rtt_wt)],
                                                    'PBS_sequence': [pbs],
                                                    'Target_name': [target_name],
                                                    'Target_sequence': [None],
                                                    'Spacer_GC_content': [None], 
                                                    'PAM_sequence': [None],
                                                    'Extension_sequence': [''.join([str(rtt_wt),pbs])], 
                                                    'Annotation': ['wildtype'], 
                                                    'pegRNA-to-edit_distance': [None],
                                                    'Nick_index': [None],
                                                    'ngRNA-to-pegRNA_distance': [None], 
                                                    'PBS_length': [len(pbs)],
                                                    'PBS_GC_content': [None], 
                                                    'RTT_length': [rtt_length], 
                                                    'RTT_GC_content': [None],
                                                    'First_extension_nucleotide': [rtt_wt[0]]})]).reset_index(drop=True)
            
            # Obtain WT RTT in-frame from + strand
            rtt_codon_frames = get_codon_frames(rtt_wt) # codons
            if comments==True:
                print(f"Extended Codons (Here): {extended_codons[math.ceil((pbs_j-rtt_length)/3)-1:math.ceil(pbs_j/3)]}")
            for i,rtt_codon_frame in enumerate(rtt_codon_frames): # Search for in-frame nucleotide sequence
                if comments==True:
                    print(f"rtt_codon_frame: {rtt_codon_frame}")

                index = found_list_in_order(extended_codons[math.ceil((pbs_j-rtt_length)/3)-1:math.ceil(pbs_j/3)],rtt_codon_frame)
                if index != -1: # Codon frame from rtt matches extended codons of in-frame nucleotide sequence
                    rtt_wt_inframe_nuc_codons_flank5 = rtt_wt[:i] # Save codon frame flank 5'
                    rtt_wt_inframe_nuc_codons = rtt_codon_frame # Save codon frame
                    rtt_wt_inframe_nuc_codons_flank3 = rtt_wt[i+3*len(rtt_codon_frame):] # Save codon frame flank 3'
                    rtt_wt_inframe_nuc = Seq('').join(rtt_codon_frame) # Join codon frame to make in-frame nucleotide sequence
                    rtt_wt_inframe_prot = Seq.translate(rtt_wt_inframe_nuc) # Translate to in-frame protein sequence
                    rtt_wt_inframe_prot_indexes = extended_codons_aa_indexes[extended_codons_prot.find(rtt_wt_inframe_prot):extended_codons_prot.find(rtt_wt_inframe_prot)+len(rtt_wt_inframe_prot)] # Obtain correponding aa indexes
                    rtt_wt_inframe_prot_deletions = f5_seq_f3_prot[f5_seq_f3_prot.find(rtt_wt_inframe_prot):f5_seq_f3_prot.find(rtt_wt_inframe_prot)+len(rtt_wt_inframe_prot)+1] # Store AAs for deletion names
                    found=True
                    break
            
            if comments==True:
                print(f'Strand: {strand}')
                print(f'Nucleotides: {rtt_wt}')
                print(f'Nucleotides 5\' of Codons: {rtt_wt_inframe_nuc_codons_flank5}')
                print(f'Nucleotides Codons: {rtt_wt_inframe_nuc_codons}')
                print(f'Nucleotides 3\' of Codons: {rtt_wt_inframe_nuc_codons_flank3}')
                print(f'Nucleotides In-Frame: {rtt_wt_inframe_nuc}')
                print(f'Amino Acids In-Frame: {rtt_wt_inframe_prot}')
                print(f'Amino Acid #s In-Frame: {rtt_wt_inframe_prot_indexes}\n')
            
            if found==False:
                raise(ValueError("RTT was not found."))
            
            # Obtain single insertion RTTs from + strand
            edits_in = []
            rtts_in = []
            for i in range(len(rtt_wt_inframe_nuc_codons)): # Iterate through all in-frame codon positions
                for codon_table_aa,codon_table_dna in aa_dna_codon_table.items(): # Obtain all possible codon insertions
                    if codon_table_aa!='*': # Remove stop codons
                        edits_in.append(f'{rtt_wt_inframe_prot[i]}{rtt_wt_inframe_prot_indexes[i]}{rtt_wt_inframe_prot[i]}{codon_table_aa}')
                        rtts_in.append(Seq('').join([rtt_wt_inframe_nuc_codons_flank5, # Codon frame flank 5'
                                       Seq('').join(rtt_wt_inframe_nuc_codons[:i+1]), # Codons before insertion
                                       Seq(codon_table_dna[0]).lower(), # Insertion codon
                                       Seq('').join(rtt_wt_inframe_nuc_codons[i+1:]), # Codons after insertion
                                       rtt_wt_inframe_nuc_codons_flank3])) # Codon frame flank 3'

            if comments==True:
                print(f'Insertions: {edits_in}')
                print(f'Insertion RTTs: {rtts_in}\n')

            insertions = pd.concat([insertions,
                                    pd.DataFrame({'pegRNA_number': [j]*len(edits_in),
                                                  'gRNA_type': ['pegRNA']*len(edits_in),
                                                  'Strand': [strand]*len(edits_in),
                                                  'Edit': edits_in,
                                                  'Spacer_sequence': [spacer]*len(edits_in),
                                                  'Scaffold_sequence': [pegRNAs.iloc[0]['Scaffold_sequence']]*len(edits_in),
                                                  'RTT_sequence': [str(rtt_in) for rtt_in in rtts_in],
                                                  'PBS_sequence': [pbs]*len(edits_in),
                                                  'Target_name': [target_name]*len(edits_in),
                                                  'Target_sequence': [None]*len(edits_in),
                                                  'Spacer_GC_content': [None]*len(edits_in),
                                                  'PAM_sequence': [None]*len(edits_in),
                                                  'Extension_sequence': [''.join([str(rtt_in),pbs]) for rtt_in in rtts_in], 
                                                  'Annotation': ['insertion']*len(edits_in),
                                                  'pegRNA-to-edit_distance': [None]*len(edits_in),
                                                  'Nick_index': [None]*len(edits_in),
                                                  'ngRNA-to-pegRNA_distance': [None]*len(edits_in),
                                                  'PBS_length': [len(pbs)]*len(edits_in),
                                                  'PBS_GC_content': [None]*len(edits_in),
                                                  'RTT_length': [len(rtt_in) for rtt_in in rtts_in], 
                                                  'RTT_GC_content': [None]*len(edits_in),
                                                  'First_extension_nucleotide': [rtt_in[0] for rtt_in in rtts_in]})]).reset_index(drop=True)

            # Obtain single deletion RTTs from + strand
            edits_del = [f'{aa}{rtt_wt_inframe_prot_deletions[i+1]}{rtt_wt_inframe_prot_indexes[i]}{rtt_wt_inframe_prot_deletions[i+1]}' for i,aa in enumerate(rtt_wt_inframe_prot) if i!=0] # Don't want first AA
            rtts_del = [Seq('').join([rtt_wt_inframe_nuc_codons_flank5, # Codon frame flank 5'
                        Seq('').join(rtt_wt_inframe_nuc_codons[:i]), # Codons before deletion
                        Seq('').join(rtt_wt_inframe_nuc_codons[i+1:]), # Codons after deletion
                        rtt_wt_inframe_nuc_codons_flank3]) # Codon frame flank 3'
                        for i in range(len(rtt_wt_inframe_nuc_codons)) if i!=0] # Don't want first AA
            
            if comments==True:
                print(f'Deletions: {edits_del}')
                print(f'Deletion RTTs: {rtts_del}\n\n')

            deletions = pd.concat([deletions,
                                   pd.DataFrame({'pegRNA_number': [j]*len(edits_del),
                                                 'gRNA_type': ['pegRNA']*len(edits_del),
                                                 'Strand': [strand]*len(edits_del),
                                                 'Edit': edits_del,
                                                 'Spacer_sequence': [spacer]*len(edits_del),
                                                 'Scaffold_sequence': [pegRNAs.iloc[0]['Scaffold_sequence']]*len(edits_del),
                                                 'RTT_sequence': [str(rtt_del) for rtt_del in rtts_del],
                                                 'PBS_sequence': [pbs]*len(edits_del),
                                                 'Target_name': [target_name]*len(edits_del),
                                                 'Target_sequence': [None]*len(edits_del),
                                                 'Spacer_GC_content': [None]*len(edits_del), 
                                                 'PAM_sequence': [None]*len(edits_del),
                                                 'Extension_sequence': [''.join([str(rtt_del),pbs]) for rtt_del in rtts_del], 
                                                 'Annotation': ['deletion']*len(edits_del), 
                                                 'pegRNA-to-edit_distance': [None]*len(edits_del),
                                                 'Nick_index': [None]*len(edits_del),
                                                 'ngRNA-to-pegRNA_distance': [None]*len(edits_del), 
                                                 'PBS_length': [len(pbs)]*len(edits_del),
                                                 'PBS_GC_content': [None]*len(edits_del),
                                                 'RTT_length': [len(rtt_del) for rtt_del in rtts_del], 
                                                 'RTT_GC_content': [None]*len(edits_del),
                                                 'First_extension_nucleotide': [rtt_del[0] for rtt_del in rtts_del]})]).reset_index(drop=True)

        else: 
            raise ValueError('Strand column can only have "+" and "-".')

    # Combine wildtype, substitution, insertion, deletion libraries
    if include_WT==True:
        pegRNAs = pd.concat([pegRNAs,wildtypes,insertions,deletions]).reset_index(drop=True)
    else:
        pegRNAs = pd.concat([pegRNAs,insertions,deletions]).reset_index(drop=True)

    # Remove pegRNAs that make edits outside in_file target sequence
    pegRNAs['AA Number'] = [int(re.findall(r'-?\d+',edit)[0]) if edit is not None else aa_index for edit in pegRNAs['Edit']]
    pegRNAs = pegRNAs[(pegRNAs['AA Number'] >= aa_index) & (pegRNAs['AA Number'] <= aa_index+len(seq)/3-1)]

    # Temporarily make pegRNAs and ngRNAs oligonucleotides
    pegRNAs['Oligonucleotide'] = [str(spacer+scaffold+rtt+pbs).upper() for (spacer, scaffold, rtt, pbs) in t.zip_cols(df=pegRNAs,cols=['Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence'])]
    
    # Check for 0 recognition sites per enzyme
    for enzyme in enzymes:
        pegRNAs = find_enzyme_sites(df=pegRNAs, enzyme=enzyme) # Find recognition sites for enzymes
        pegRNAs_edits = list(pegRNAs['Edit'].unique()) # Get pegRNA edits

        if replace: # Replace pegRNAs with RE enzyme sites
            
            # Store pegRNAs with recognition sites for enzymes
            pegRNAs_enzyme = pegRNAs[pegRNAs[enzyme]!=0]
            io.save(dir=f'../rtt_designer/{enzyme}/codon_swap_before',
                    file=f'{int(pegRNAs_enzyme.iloc[0]['PBS_length'])}.csv',
                    obj=pegRNAs_enzyme)
            
            # Codon swap pegRNAs with enzyme recognition site
            pegRNAs_enzyme = enzyme_codon_swap(pegRNAs=pegRNAs_enzyme,in_file=in_file,enzyme=enzyme)
            io.save(dir=f'../rtt_designer/{enzyme}/codon_swap_after',
                    file=f'{int(pegRNAs_enzyme.iloc[0]['PBS_length'])}.csv',
                    obj=pegRNAs_enzyme)
            pegRNAs = pd.concat([pegRNAs,pegRNAs_enzyme],ignore_index=True)
            print(f"pegRNA edits recovered by modifying {enzyme} recognition site: {list(pegRNAs_enzyme['Edit'].unique())}")

            # Recheck pegRNAs for RE recognition sites and drop those with recognition sites
            pegRNAs = find_enzyme_sites(df=pegRNAs, enzyme=enzyme)
            pegRNAs = pegRNAs[pegRNAs[enzyme]==0].reset_index(drop=True)

            # Store removed edits
            pegRNAs_enzyme = pegRNAs[pegRNAs[enzyme]!=0]
            remove_pegRNAs_edits = pegRNAs[pegRNAs['Edit'].isin(pegRNAs_enzyme['Edit'])]['Edit'].unique()
            # Save lost edits
            lost_pegRNAs_edits = [remove_edit for remove_edit in remove_pegRNAs_edits if remove_edit not in pegRNAs_edits]
            if len(lost_pegRNAs_edits) > 0:
                print(f"pegRNAs edits lost due to {enzyme} recognition site: {lost_pegRNAs_edits}")
                io.save(dir=f'../rtt_designer/{enzyme}/lost',
                            file=f'{int(pegRNAs_enzyme.iloc[0]['PBS_length'])}.csv',
                            obj=pegRNAs_enzyme[pegRNAs_enzyme['Edit'].isin(lost_pegRNAs_edits)])

            # Drop enzyme column
            pegRNAs.drop(columns=[enzyme,f'{enzyme}_fwd_i',f'{enzyme}_rc_i'],inplace=True)

    # Remove oligonucleotide column
    pegRNAs = pegRNAs.drop(columns=['Oligonucleotide'])

    # Save & Return
    memories.append(memory_timer(task=f"rtt_designer()"))
    if out_dir is not None and out_file is not None:
        io.save(dir=os.path.join(out_dir,f'.rtt_designer'),
                file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_memories.csv',
                obj=pd.DataFrame(memories, columns=['Task','Memory, MB','Time, s']))
        io.save(dir=out_dir,file=out_file,obj=pegRNAs)
    if return_df==True: return pegRNAs

def pegRNA_outcome(pegRNAs: pd.DataFrame | str, in_file: pd.DataFrame | str,
                   match_score: float = 2, mismatch_score: float = -1, open_gap_score: float = -10, extend_gap_score: float = -0.1,
                   out_dir: str=None, out_file: str=None, return_df: bool=True, literal_eval: bool=True) -> pd.DataFrame:
    ''' 
    pegRNA_outcome(): confirm that pegRNAs should create the predicted edit
    
    Parameters:
    pegRNAs (dataframe | str): pegRNAs DataFrame (or file path)
    in_file (dataframe | str): Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (column names required)
    match_score (float, optional): match score for pairwise alignment (Default: 2)
    mismatch_score (float, optional): mismatch score for pairwise alignment (Default: -1)
    open_gap_score (float, optional): open gap score for pairwise alignment (Default: -10)
    extend_gap_score (float, optional): extend gap score for pairwise alignment (Default: -0.1)
    out_dir (str, optional): output directory
    out_file (str, optional): output filename
    return_df (bool, optional): return dataframe (Default: True)
    literal_eval (bool, optional): convert string representations (Default: True)
    
    Dependencies: Bio, numpy, pandas, fastq, datetime, re, os, memory_timer(), io
    '''
    # Initialize timer; memory reporting
    memory_timer(reset=True)
    memories = []

    # Get pegRNAs & PrimeDesign input DataFrame from file path if needed
    if type(pegRNAs)==str:
        pegRNAs = io.get(pt=pegRNAs,literal_eval=literal_eval)
    if type(in_file)==str:
        in_file = io.get(pt=in_file,literal_eval=literal_eval)

    # Catch all stop codons that are written as "X" instead of "*"
    pegRNAs['Edit'] = pegRNAs['Edit'].replace('X', '*', regex=True)

   # Get reference sequence & codons (+ reverse complement)
    target_sequence = in_file.iloc[0]['target_sequence'] 
    seq = Seq(target_sequence.split('(')[1].split(')')[0]) # Break apart target sequences
    if len(seq)%3 != 0: raise(ValueError(f"Length of target sequence ({len(seq)}) must divisible by 3. Check input file."))
    flank5 = Seq(target_sequence.split('(')[0])
    if len(flank5)%3 != 0: raise(ValueError(f"Length of flank5 ({len(flank5)}) must divisible by 3. Check input file."))
    flank3 = Seq(target_sequence.split(')')[1])
    if len(flank3)%3 != 0: raise(ValueError(f"Length of flank3 ({len(flank3)}) must divisible by 3. Check input file."))

    aa_index = in_file.iloc[0]['aa_index']

    f5_seq_f3_nuc = flank5 + seq + flank3  # Join full nucleotide reference sequence
    rc_f5_seq_f3_nuc = Seq.reverse_complement(f5_seq_f3_nuc) # Full nucleotide reference reverse complement sequence
    seq_prot = Seq.translate(seq) # In-frame amino acid sequence
    f5_seq_f3_prot = Seq.translate(f5_seq_f3_nuc) # Full in-frame protein sequence (including flanks)
    
    aa_indexes = list(np.arange(aa_index,aa_index+len(seq_prot))) # In-frame amino acid indexes
    seq_prot_deletions = Seq.translate(seq)+Seq.translate(flank3[:3]) # In-frame amino acid sequence + next AA for deletion names
    
    print(f'FWD Ref: {f5_seq_f3_nuc}')
    print(f'REV Ref: {rc_f5_seq_f3_nuc}')
    print(f'Nucleotides: {seq}')
    print(f'Amino Acids: {seq_prot}\n')

    # Get expected edits in the pegRNA library 
    edits_substitutions=[]
    edits_insertions=[]
    edits_deletions=[]
    for i,aa in enumerate(seq_prot):
        edits_substitutions.extend([f'{aa}{str(aa_indexes[i])}{aa2}' for aa2 in aa_dna_codon_table if (aa2!='*')&(aa2!=aa)])
        edits_insertions.extend([f'{aa}{str(aa_indexes[i])}{aa}{aa2}' for aa2 in aa_dna_codon_table if aa2!='*'])
        edits_deletions.append(f'{aa}{seq_prot_deletions[i+1]}{str(aa_indexes[i])}{seq_prot_deletions[i+1]}')
    edits_substitutions_set = set(edits_substitutions)
    edits_insertions_set = set(edits_insertions)
    edits_deletions_set = set(edits_deletions)
    
    print(f'Expected Edits in pegRNA library...\nSubstitutions: {edits_substitutions}\nInsertions: {edits_insertions}\nDeletions: {edits_deletions}')
    print(f'All substitutions present: {edits_substitutions_set.issubset(set(pegRNAs["Edit"]))}; missing: {edits_substitutions_set-set(pegRNAs["Edit"])}')
    print(f'All insertions present: {edits_insertions_set.issubset(set(pegRNAs["Edit"]))}; missing: {edits_insertions_set-set(pegRNAs["Edit"])}')
    print(f'All deletions present: {edits_deletions_set.issubset(set(pegRNAs["Edit"]))}; missing: {edits_deletions_set-set(pegRNAs["Edit"])}\n')

    # Determine post_RTT_sequences
    post_RTT_sequences = [] # Store post RTT sequences
    for (strand,pbs,rtt,edit,aa_number) in t.zip_cols(df=pegRNAs,cols=['Strand','PBS_sequence','RTT_sequence','Edit','AA Number']): # Iterate through primer binding sites

        if strand=='+': # Spacer: + strand; PBS & RTT: - strand

            # Find reverse complement PBS in sequence
            rc_pbs = Seq.reverse_complement(Seq(pbs)) # reverse complement of pbs (+ strand)
            
            rc_pbs_j = f5_seq_f3_nuc.find(str(rc_pbs))
            if rc_pbs_j == -1:
                raise ValueError(f"PBS sequence '{pbs}' not found in target sequence. Please check the input file.")
            elif rc_pbs_j != f5_seq_f3_nuc.rfind(str(rc_pbs)):
                print(rc_pbs,rc_pbs_j,rc_f5_seq_f3_nuc.rfind(str(rc_pbs)))
                raise ValueError(f"Multiple matches found for PBS sequence '{str(rc_pbs)}'. Please check the input file.")

            # Replace change sequence using reverse complement RTT
            if len(edit)==len(str(aa_number))+2: # Substitution
                post_RTT_sequences.append(str(f5_seq_f3_nuc[:rc_pbs_j+len(rc_pbs)]+
                                            Seq.reverse_complement(Seq(rtt.upper()))+
                                            f5_seq_f3_nuc[rc_pbs_j+len(rc_pbs)+len(rtt):])) # Save post RTT sequence
            
            elif edit.find(str(aa_number))>len(edit)-edit.find(str(aa_number))-len(str(aa_number)): # Deletion
                post_RTT_sequences.append(str(f5_seq_f3_nuc[:rc_pbs_j+len(rc_pbs)]+
                                            Seq.reverse_complement(Seq(rtt.upper()))+
                                            f5_seq_f3_nuc[rc_pbs_j+len(rc_pbs)+len(rtt)+3*(edit.find(str(aa_number))-1):]))
            
            else: # Insertion
                post_RTT_sequences.append(str(f5_seq_f3_nuc[:rc_pbs_j+len(rc_pbs)]+
                                            Seq.reverse_complement(Seq(rtt.upper()))+
                                            f5_seq_f3_nuc[rc_pbs_j+len(rc_pbs)+len(rtt)-3*(len(edit)-edit.find(str(aa_number))-len(str(aa_number))-1):]))
            
        elif strand=='-': # Spacer: - strand; PBS & RTT: + strand

            # Find PBS in sequence
            pbs = Seq(pbs) # pbs (+ strand)
            pbs_j = f5_seq_f3_nuc.find(str(pbs))
            if pbs_j == -1:
                raise ValueError(f"PBS sequence '{pbs}' not found in target sequence. Please check the input file.")
            elif pbs_j != f5_seq_f3_nuc.rfind(str(pbs)):
                print(pbs,pbs_j,f5_seq_f3_nuc.rfind(str(pbs)))
                raise ValueError(f"Multiple matches found for PBS sequence '{pbs}'. Please check the input file.")
            
            # Replace change sequence using RTT
            if len(edit)==len(str(aa_number))+2: # Substitution
                post_RTT_sequences.append(str(f5_seq_f3_nuc[:pbs_j-len(rtt)]+
                                              Seq(rtt.upper())+
                                              f5_seq_f3_nuc[pbs_j:]))
            elif edit.find(str(aa_number))>len(edit)-edit.find(str(aa_number))-len(str(aa_number)): # Deletion
                post_RTT_sequences.append(str(f5_seq_f3_nuc[:pbs_j-len(rtt)-3*(edit.find(str(aa_number))-1)]+
                                            Seq(rtt.upper())+
                                            f5_seq_f3_nuc[pbs_j:]))
                
            else: # Insertion
                post_RTT_sequences.append(str(f5_seq_f3_nuc[:pbs_j-len(rtt)+3*(len(edit)-edit.find(str(aa_number))-len(str(aa_number))-1)]+
                                              Seq(rtt.upper())+
                                              f5_seq_f3_nuc[pbs_j:]))

        else: 
            raise(ValueError('Error: Strand column can only have "+" and "-".'))
        
    # Determine edit from post RTT sequences
    pegRNAs['post_RTT_sequence']=post_RTT_sequences
    
    # Check edits & assign multiple edit annotations if needed
    edit_check = []
    for post_RTT_sequence in pegRNAs['post_RTT_sequence']:
        if len(f5_seq_f3_nuc)!=len(post_RTT_sequence): # Indel
            edit = fq.find_indel(wt=f5_seq_f3_nuc, mut=post_RTT_sequence, res=int(aa_index-len(flank5)/3), show=False, 
                                 match_score=match_score, mismatch_score=mismatch_score,
                                 open_gap_score=open_gap_score, extend_gap_score=extend_gap_score)[0]
            
            # Check for additional edits
            aa_number = int(re.findall(r'-?\d+',edit)[0])
            i = int(aa_number-aa_index+len(flank5)/3)
            if edit.find(str(aa_number))>len(edit)-edit.find(str(aa_number))-len(str(aa_number)): # Deletion
                
                # Look for next AA(s) that match the deleted AA in the edit
                i_ls = [i]
                for j,aa_j in enumerate(f5_seq_f3_prot[i+1:]):
                    if aa_j==edit[0]:
                        i_ls.append(i+j+1)
                    else:
                        break
                
                if len(i_ls)>1: # Multiple edit annotations
                    edits = [f'{f5_seq_f3_prot[i]}{f5_seq_f3_prot[i+1]}{int(i+aa_index-len(flank5)/3)}{f5_seq_f3_prot[i+1]}' for i in i_ls]
                    edit_check.append(edits)
                else: # Single edit annotation
                    edit_check.append(edit)

            else: # Insertion
                
                # Look for next AA(s) that match the inserted AA in the edit
                i_ls = [i]
                for j,aa_j in enumerate(f5_seq_f3_prot[i+1:]):
                    if aa_j==edit[-1]:
                        i_ls.append(i+j+1)
                    else:
                        break
                
                if len(i_ls)>1: # Multiple edit annotations
                    edits = [f'{f5_seq_f3_prot[i]}{int(i+aa_index-len(flank5)/3)}{f5_seq_f3_prot[i]}{edit[-1]}' for i in i_ls]
                    edit_check.append(edits)
                else: 
                    edit_check.append(edit)
            
        else: # Substitution
            edit_check.append(fq.find_AA_edits(wt=str(Seq.translate(f5_seq_f3_nuc)), 
                                               res=int(aa_index-len(flank5)/3), 
                                               seq=str(Seq.translate(Seq(post_RTT_sequence)))))
    pegRNAs['Edit_check'] = edit_check
    
    # Compare Edit_check with Edit
    pegRNAs['Edit_check_match'] = [edit_check==edit if isinstance(edit_check,str) else edit in edit_check for (edit_check,edit) in t.zip_cols(df=pegRNAs,cols=['Edit_check','Edit'])]
    print(f"All pegRNAs passed edit check: {all(pegRNAs['Edit_check_match'])}")

    # Save & Return
    memories.append(memory_timer(task=f"pegRNA_outcome(): {len(pegRNAs)} out of {len(pegRNAs)}"))
    if out_dir is not None and out_file is not None:
        io.save(dir=os.path.join(out_dir,f'.pegRNA_outcome'),
                file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_memories.csv',
                obj=pd.DataFrame(memories, columns=['Task','Memory, MB','Time, s']))
        io.save(dir=out_dir,file=out_file,obj=pegRNAs)
    if return_df==True: return pegRNAs

def pegRNA_signature(pegRNAs: pd.DataFrame | str, in_file: pd.DataFrame | str, 
                     match_score: float = 2, mismatch_score: float = -1, open_gap_score: float = -10, extend_gap_score: float = -0.1,
                     out_dir: str=None, out_file: str=None, save_alignments: bool=False, return_df: bool=True, literal_eval: bool=True) -> pd.DataFrame:
    ''' 
    pegRNA_signature(): create signatures for pegRNA outcomes using alignments
    
    pegRNAs (dataframe | str): pegRNAs DataFrame (or file path)
    in_file (dataframe | str): Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence,aa_index (column names required)
    match_score (float, optional): match score for pairwise alignment (Default: 2)
    mismatch_score (float, optional): mismatch score for pairwise alignment (Default: -1)
    open_gap_score (float, optional): open gap score for pairwise alignment (Default: -10)
    extend_gap_score (float, optional): extend gap score for pairwise alignment (Default: -0.1)
    out_dir (str, optional): output directory
    out_file (str, optional): output filename
    save_alignments (bool, optional): save alignments (Default: False, save memory)
    return_df (bool, optional): return dataframe (Default: True)
    literal_eval (bool, optional): convert string representations (Default: True)
    '''
    # Initialize timer; memory reporting
    memory_timer(reset=True)
    memories = []

    # Get pegRNAs & PrimeDesign input DataFrame from file path if needed
    if type(pegRNAs)==str:
        pegRNAs = io.get(pt=pegRNAs,literal_eval=literal_eval)
    if type(in_file)==str:
        in_file = io.get(pt=in_file,literal_eval=literal_eval)

    # High sequence homology; punish gaps
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = match_score  # Score for a match
    aligner.mismatch_score = mismatch_score  # Penalty for a mismatch; applied to both strands
    aligner.open_gap_score = open_gap_score  # Penalty for opening a gap; applied to both strands
    aligner.extend_gap_score = extend_gap_score  # Penalty for extending a gap; applied to both strands

    # Get wt sequence
    target_sequence = in_file.iloc[0]['target_sequence'] 
    seq = target_sequence.split('(')[1].split(')')[0] # Break apart target sequences
    if len(seq)%3 != 0: raise(ValueError(f"Length of target sequence ({len(seq)}) must divisible by 3. Check input file."))
    flank5 = target_sequence.split('(')[0]
    if len(flank5)%3 != 0: raise(ValueError(f"Length of flank5 ({len(flank5)}) must divisible by 3. Check input file."))
    flank3 = target_sequence.split(')')[1]
    if len(flank3)%3 != 0: raise(ValueError(f"Length of flank3 ({len(flank3)}) must divisible by 3. Check input file."))

    # Create alignments and Signatures
    start_i = len(flank5)
    end_i = start_i + len(seq)
    pegRNAs['Alignment'] = [aligner.align(Seq(seq),
                                          Seq((post_RTT_sequence[start_i:end_i])))[0] for post_RTT_sequence in pegRNAs['post_RTT_sequence']]
    pegRNAs['Signature'] = [signature_from_alignment(ref_seq=seq,
                                                     query_seq=post_RTT_sequence[start_i:end_i],
                                                     alignment=alignment) for alignment,post_RTT_sequence in t.zip_cols(df=pegRNAs, cols=['Alignment', 'post_RTT_sequence'])]
    if save_alignments==False:
        pegRNAs.drop(columns=['Alignment'],inplace=True)

    # Save & Return
    memories.append(memory_timer(task=f"pegRNA_signature(): {len(pegRNAs)} out of {len(pegRNAs)}"))
    if out_dir is not None and out_file is not None:
        io.save(dir=os.path.join(out_dir,f'.pegRNA_signature'),
                file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_memories.csv',
                obj=pd.DataFrame(memories, columns=['Task','Memory, MB','Time, s']))
        io.save(dir=out_dir,file=out_file,obj=pegRNAs)
    if return_df==True: return pegRNAs

# Comparing pegRNA libraries
def print_shared_sequences(dc: dict):
    ''' 
    print_shared_sequences(): prints spacer and PBS sequences from dictionary of shared_sequences libraries
    
    Parameters:
    dc (dict): dictionary of shared_sequences() libraries

    Dependencies: pandas
    '''
    keys_a = sorted(dc.keys())

    text = f""
    for key in keys_a: 
        text += f"\t{key}_spacer\t\t{key}_PBS\t\t"

    for v in range(len(dc[keys_a[0]])):
        text += f"\n{v}:\t"
        for key in keys_a:
            text += f"{dc[key].iloc[v]['Spacer_sequence']}\t{dc[key].iloc[v]['PBS_sequence']}\t\t"
    print(text)

def print_shared_sequences_mutant(dc: dict):
    ''' 
    print_shared_sequences_mutant(): prints spacer and PBS sequences as well as priority mutant from dictionary of shared_sequences libraries
    
    Parameters:
    dc (dict): dictionary of shared_sequences() libraries with priority mutant

    Depedencies: pandas
    '''
    keys_a = sorted(dc.keys())

    text = f""
    for key in keys_a: 
        text += f"\t{key}_spacer\t\t{key}_PBS\t\t{key}_mutant"

    for v in range(len(dc[keys_a[0]])):
        text += f"\n{v}:\t"
        for key in keys_a:
            text += f"{dc[key].iloc[v]['Spacer_sequence']}\t{dc[key].iloc[v]['PBS_sequence']}\t{dc[key].iloc[v]['Priority_mut']}\t\t"
    print(text)

# Comparing pegRNAs
def group_pe(df: pd.DataFrame, other_cols: list, epegRNA_id_col: str='epegRNA', ngRNA_id_col: str='ngRNA',
             epegRNA_spacer_col: str='Spacer_sequence_epegRNA', epegRNA_RTT_col: str='RTT_sequence',epegRNA_PBS_col: str='PBS_sequence',
             match_score: float = 2, mismatch_score: float = -1, open_gap_score: float = -10, extend_gap_score: float = -0.1):
    '''
    group_pe(): returns a dataframe containing groups of (epegRNA,ngRNA) pairs that share spacers and have similar PBS and performs pairwise alignment for RTT
    
    Parameters:
    df (dataframe): dataframe
    other_cols (list): names of other column that will be retained
    epegRNA_id_col (str, optional): epegRNA id column name (Default: epegRNA)
    ngRNA_id_col (str, optional): ngRNA id column name (Default: ngRNA)
    epegRNA_spacer_col (str, optional): epegRNA spacer column name (Default: Spacer_sequence_epegRNA)
    epegRNA_RTT_col (str, optional): epegRNA reverse transcripase template column name (Default: RTT_sequence_epegRNA)
    epegRNA_PBS_col (str, optional): epegRNA primer binding site column name (Default: PBS_sequence_epegRNA
    match_score (float, optional): match score for pairwise alignment (Default: 2)
    mismatch_score (float, optional): mismatch score for pairwise alignment (Default: -1)
    open_gap_score (float, optional): open gap score for pairwise alignment (Default: -10)
    extend_gap_score (float, optional): extend gap score for pairwise alignment (Default: -0.1)
    
    Dependencies: pandas,Bio
    '''
    # High sequence homology; punish gaps
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = match_score  # Score for a match
    aligner.mismatch_score = mismatch_score  # Penalty for a mismatch; applied to both strands
    aligner.open_gap_score = open_gap_score  # Penalty for opening a gap; applied to both strands
    aligner.extend_gap_score = extend_gap_score  # Penalty for extending a gap; applied to both strands

    # Isolate desired columns
    other_cols.extend([epegRNA_id_col,ngRNA_id_col,epegRNA_spacer_col,epegRNA_RTT_col,epegRNA_PBS_col])
    df = df[other_cols]

    df_pairs = pd.DataFrame() # (epegRNA,ngRNA) pairs dataframe
    for epegRNA_id in list(df[epegRNA_id_col].value_counts().keys()): # Iterate through epegRNA ids
        
        # Split dataframe to isolate 1 epegRNA from others with the same spacer
        df_epegRNA = df[df[epegRNA_id_col]==epegRNA_id].reset_index(drop=True)
        df_others = df[(df[epegRNA_id_col]!=epegRNA_id)&(df[epegRNA_spacer_col]==df_epegRNA.iloc[0][epegRNA_spacer_col])].reset_index(drop=True)

        # Iterate through (epegRNA,ngRNA) pairs and isolate...
        for i,(ngRNA,epegRNA_RTT,epegRNA_PBS) in enumerate(t.zip_cols(df=df_epegRNA,cols=[ngRNA_id_col,epegRNA_RTT_col,epegRNA_PBS_col])):
            df_others = df_others[df_others[ngRNA_id_col]==ngRNA].reset_index(drop=True) # shared ngRNAs
            df_others = df_others[(df_others[epegRNA_PBS_col].str.contains(epegRNA_PBS))|(epegRNA_PBS in df_others[epegRNA_PBS_col])].reset_index(drop=True) # similar PBS
            
            if df_others.empty==False: # Only retain successful pairs
                df_others['PBS_lengths'] = [f'({len(epegRNA_PBS)},{len(other_epegRNA_PBS)})' for other_epegRNA_PBS in df_others[epegRNA_PBS_col]] # Get PBS lengths
                
                # Quantify mismatches in RTT alignments
                RTT_alignments = []
                RTT_alignments_mismatches = []
                for other_epegRNA_RTT in df_others[epegRNA_RTT_col]:
                    RTT_alignment = aligner.align(epegRNA_RTT,other_epegRNA_RTT)[0]
                    RTT_alignments.append(RTT_alignment)
                    RTT_alignments_mismatches.append(int(len(epegRNA_RTT)-RTT_alignment.score))
                df_others['RTT_alignment'] = RTT_alignments
                df_others['RTT_alignments_mismatches'] = RTT_alignments_mismatches
                
                series_df_epegRNA = pd.concat([pd.DataFrame([df_epegRNA.iloc[i]])]*(len(df_others))).reset_index(drop=True)

                df_pair = pd.concat([df_others,series_df_epegRNA.rename(columns=lambda col: f"{col}_compare")],axis=1) # Append compared (epegRNA,ngRNA)
                df_pairs = pd.concat([df_pairs,df_pair]).reset_index(drop=True) # Save (epegRNA,ngRNA) pairs to output dataframe
    
    return df_pairs