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
- pegRNA_outcome(): confirm that pegRNAs should create the predicted edits
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

def enzyme_codon_swap(pegRNAs: pd.DataFrame | str, enzyme: str, 
                      RE_type_IIS_df: pd.DataFrame | str = None, out_dir: str = None, 
                      out_file: str = None, return_df: bool = True, literal_eval: bool=True, comments: bool=False) -> pd.DataFrame:
    '''
    enzyme_codon_swap(): modify pegRNA RTT sequences to disrupt a RE recognition site

    Parameters:
    pegRNAs (pd.DataFrame | str): pegRNAs DataFrame or file path to pegRNAs DataFrame
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

    # Get pegRNAs & RE_type_IIS DataFrames from file path if needed
    if type(pegRNAs)==str:
        pegRNAs = io.get(pt=pegRNAs, literal_eval=literal_eval)
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

    # Get new RTT and edit sequences
    new_rtt_ls = []
    new_edit_sequence_ls = []
    for (strand,spacer,scaffold,rtt,pbs,rtt_length,enzyme_fwd_i_ls,enzyme_rc_i_ls,reference_sequence,edit_sequence) in t.zip_cols(df=pegRNAs,
        cols=['Strand','Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence', 'RTT_length',f'{enzyme}_fwd_i',f'{enzyme}_rc_i','Reference_sequence','Edit_sequence']):
        
        # Get reverse complement and codons of reference sequence
        reference_sequence = Seq(reference_sequence)
        rc_reference_sequence = Seq.reverse_complement(reference_sequence)
        reference_sequence_codons = get_codons(reference_sequence)

        if strand=='+': # Spacer: + strand; PBS & RTT: - strand
            
            # Find spacer in sequence
            spacer_j = reference_sequence.find(spacer)
            if spacer_j == -1:
                raise ValueError(f"Spacer sequence '{spacer}' not found in reference sequence: {reference_sequence}")
            elif spacer_j != reference_sequence.rfind(spacer):
                raise ValueError(f"Multiple matches found for spacer sequence '{spacer}'.")

            # Find PBS in reverse complement sequence
            pbs_j = rc_reference_sequence.find(pbs)
            if pbs_j == -1:
                raise ValueError(f"PBS sequence '{pbs}' not found in reference sequence: {reference_sequence}")
            elif pbs_j != rc_reference_sequence.rfind(pbs):
                print(pbs,pbs_j,rc_reference_sequence.rfind(pbs))
                raise ValueError(f"Multiple matches found for PBS sequence '{pbs}'.")

            # Obtain reverse complement WT RTT & edit RTT in-frame from + strand
            rc_rtt = Seq.reverse_complement(Seq(rtt)) # reverse complement of rtt (+ strand)
            rc_rtt_codon_frames = get_codon_frames(rc_rtt) # codons

            rtt_wt = rc_reference_sequence[pbs_j-int(rtt_length):pbs_j]
            rc_rtt_wt = Seq.reverse_complement(rtt_wt) # reverse complement of rtt wt (+ strand)
            rc_rtt_wt_codon_frames = get_codon_frames(rc_rtt_wt) # codons
            if comments==True:
                print(f"Reference Sequence Codons (Here): {reference_sequence_codons[math.floor((len(rc_reference_sequence)-pbs_j)/3)-1:math.floor((len(rc_reference_sequence)-pbs_j+int(rtt_length))/3)]}")
            for i,(rc_rtt_wt_codon_frame,rc_rtt_codon_frame) in enumerate(zip(rc_rtt_wt_codon_frames,rc_rtt_codon_frames)): # Search for in-frame nucleotide sequence
                if comments==True:
                    print(f"rc_rtt_wt_codon_frame: {rc_rtt_wt_codon_frame}")
                
                index = found_list_in_order(reference_sequence_codons[math.floor((len(rc_reference_sequence)-pbs_j)/3)-1:math.floor((len(rc_reference_sequence)-pbs_j+rtt_length)/3)],rc_rtt_wt_codon_frame)
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
                new_edit_sequence_ls.append(None)
                continue
            
            if comments==True:
                print(f'Enzyme site index: {enzyme_i}')
                print(f'Enzyme codon index: {enzyme_codon_i}')

            # Change codon swap enzyme site; save new RTT sequence and edit sequence
            codons = [str(codon).lower() for codon in aa_dna_codon_table[str(rc_rtt_inframe_prot[enzyme_codon_i])] if str(codon).upper() != str(rc_rtt_inframe_nuc_codons[enzyme_codon_i]).upper()]
            if comments==True:
                print(f'Codons: {codons}')
            if len(codons)!=0:
                rc_rtt_inframe_nuc_codons[enzyme_codon_i] = codons[0]
                new_rtt = str(Seq.reverse_complement(Seq(rc_rtt_inframe_nuc_codons_flank5)+Seq('').join(rc_rtt_inframe_nuc_codons)+Seq(rc_rtt_inframe_nuc_codons_flank3)))
                new_rtt_ls.append(new_rtt)
                edit_sequence = str(edit_sequence)
                rc_rtt = str(rc_rtt)
                rc_rtt_edit_sequence_i = edit_sequence.upper().find(rc_rtt.upper())
                if rc_rtt_edit_sequence_i == -1:
                    raise ValueError(f"RTT sequence (RC) '{rc_rtt}' not found in edit sequence '{edit_sequence}'.")
                new_edit_sequence = edit_sequence[:rc_rtt_edit_sequence_i] + str(Seq.reverse_complement(Seq(new_rtt))) + edit_sequence[rc_rtt_edit_sequence_i+len(new_rtt):]
                new_edit_sequence_ls.append(new_edit_sequence)

                if comments==True:
                    print(f'Nucleotides In-Frame (New): {Seq('').join(rc_rtt_inframe_nuc_codons)}')
                    print(f'Amino Acid In-Frame (New): {Seq.translate(Seq('').join(rc_rtt_inframe_nuc_codons))}')
                    print(f'RTT (New): {new_rtt}')
                    print(f'RTT (Old): {rtt}')
                    print(f'RTT (WT): {rtt_wt}')
                    print(f'Edit Sequence (New): {new_edit_sequence}')
                    print(f'Edit Sequence (Old): {edit_sequence}')
                
            else:
                new_rtt_ls.append(None)
                new_edit_sequence_ls.append(None)

        elif strand=='-': # Spacer: - strand; PBS & RTT: + strand
            
            # Find spacer in sequence
            spacer_j = rc_reference_sequence.find(spacer)
            if spacer_j == -1:
                raise ValueError(f"Spacer sequence '{spacer}' not found in reference sequence: {reference_sequence}")
            if spacer_j != rc_reference_sequence.rfind(spacer):
                raise ValueError(f"Multiple matches found for spacer sequence '{spacer}' not found in reference sequence. Please check the input file.")

            # Find PBS in sequence
            pbs_j = reference_sequence.find(pbs)
            if pbs_j == -1:
                raise ValueError(f"PBS sequence '{pbs}' not found in reference sequence: {reference_sequence}")
            if pbs_j != reference_sequence.rfind(pbs):
                print(pbs,pbs_j,reference_sequence.rfind(pbs))
                raise ValueError(f"Multiple matches found for PBS sequence '{pbs}' not found in reference sequence")

            # Obtain WT RTT & edit RTT in-frame from + strand
            rtt_codon_frames = get_codon_frames(rtt) # codons

            rtt_wt = reference_sequence[pbs_j-int(rtt_length):pbs_j]
            rtt_wt_codon_frames = get_codon_frames(rtt_wt) # codons
            if comments==True:
                print(f"Reference Sequence Codons (Here): {reference_sequence_codons[math.ceil((pbs_j-rtt_length)/3)-1:math.ceil(pbs_j/3)]}")
            for i,(rtt_wt_codon_frame,rtt_codon_frame) in enumerate(zip(rtt_wt_codon_frames,rtt_codon_frames)): # Search for in-frame nucleotide sequence
                if comments==True:
                    print(f"rtt_wt_codon_frame: {rtt_wt_codon_frame}")

                index = found_list_in_order(reference_sequence_codons[math.ceil((pbs_j-rtt_length)/3)-1:math.ceil(pbs_j/3)],rtt_wt_codon_frame)
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
                new_edit_sequence_ls.append(None)
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
                new_rtt = str(Seq(rtt_inframe_nuc_codons_flank5)+Seq('').join(rtt_inframe_nuc_codons)+Seq(rtt_inframe_nuc_codons_flank3))
                new_rtt_ls.append(new_rtt)
                edit_sequence = str(edit_sequence)
                rtt = str(rtt)
                rtt_edit_sequence_i = edit_sequence.upper().find(rtt.upper())
                if rtt_edit_sequence_i == -1:
                    raise ValueError(f"RTT sequence '{rtt}' not found in edit sequence '{edit_sequence}'.")
                new_edit_sequence = edit_sequence[:rtt_edit_sequence_i] + new_rtt + edit_sequence[rtt_edit_sequence_i+len(new_rtt):]
                new_edit_sequence_ls.append(new_edit_sequence)

                if comments==True:
                    print(f'Nucleotides In-Frame (New): {Seq('').join(rtt_inframe_nuc_codons)}')
                    print(f'Amino Acid In-Frame (New): {Seq.translate(Seq('').join(rtt_inframe_nuc_codons))}')
                    print(f'RTT (New): {new_rtt}')
                    print(f'RTT (Old): {rtt}')
                    print(f'RTT (WT): {rtt_wt}')
                    print(f'Edit Sequence (New): {new_edit_sequence}')
                    print(f'Edit Sequence (Old): {edit_sequence}')

            else:
                new_rtt_ls.append(None)
                new_edit_sequence_ls.append(None)

    # Update pegRNAs DataFrame with new RTT and edit sequences
    pegRNAs['RTT_sequence'] = new_rtt_ls
    pegRNAs['Edit_sequence'] = new_edit_sequence_ls
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
                    index: int=1, silent_mutation: bool=True,
                    dir: str='.', file: str='prime_design_input.csv'):
    ''' 
    prime_design_input(): creates and checks PrimeDesign saturation mutagenesis input file
    
    Parameters:
    target_name (str): name of target
    flank5_sequence: in-frame nucleotide sequence with 5' of saturation mutagensis region (length must be divisible by 3)
    target_sequence (str): in-frame nucleotide sequence for the saturation mutagensis region (length must be divisible by 3)
    flank3_sequence: in-frame nucleotide sequence with 3' of saturation mutagensis region (length must be divisible by 3)
    index (int, optional): 1st amino acid or base in target sequence index (Default: 1)
    silent_mutation (bool, optional): check that sequences are in-frame for silent mutation option (Default: True)
    dir (str, optional): name of the output directory 
    file (str, optional): name of the output file
    
    Dependencies: pandas & io
    
    Reference: https://github.com/pinellolab/PrimeDesign/tree/master/PrimeDesign
    '''
    # Check PrimeDesign saturation mutagenesis input file 
    if silent_mutation == True:
        if len(flank5_sequence)%3 != 0: raise(ValueError(f"Length of flank5_sequence ({len(flank5_sequence)}) must divisible by 3."))
        if len(target_sequence)%3 != 0: raise(ValueError(f"Length of target_sequence ({len(target_sequence)}) must divisible by 3."))
        if len(flank3_sequence)%3 != 0: raise(ValueError(f"Length of flank3_sequence ({len(flank3_sequence)}) must divisible by 3."))

    # Create PrimeDesign saturation mutagenesis input file
    io.save(dir=dir,
            file=file,
            obj=pd.DataFrame({'target_name': [target_name],
                              'target_sequence': [f"{flank5_sequence}({target_sequence}){flank3_sequence}"],
                              'index': [index]}))

def prime_design(file: str, pe_format: str = 'NNNNNNNNNNNNNNNNN/NNN[NGG]', pbs_length_list: list = [], rtt_length_list: list = [], 
                nicking_distance_minimum: int = 0, nicking_distance_maximum: int = 100, filter_c1_extension: bool = False, silent_mutation: bool = False,
                genome_wide_design: bool = False, saturation_mutagenesis: str = None, number_of_pegrnas: int = 3, number_of_ngrnas: int = 3,
                nicking_distance_pooled: int = 75, homology_downstream: int = 15, pbs_length_pooled: int = 14, rtt_max_length_pooled: int = 50,
                out_dir: str = './DATETIMESTAMP_PrimeDesign'):
    ''' 
    prime_design(): execute PrimeDesign (EDMS version)
    
    Parameters:
    file (str): input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence,index (column names required)
    pe_format (str, optional): Prime editing formatting including the spacer, cut index -> /, and protospacer adjacent motif (PAM) -> [PAM] (Default: NNNNNNNNNNNNNNNNN/NNN[NGG])
    pbs_length_list (list, optional): list of primer binding site (PBS) lengths for the pegRNA extension. Example: 12 13 14 15
    rtt_length_list (list, optional): list of reverse transcription (RT) template lengths for the pegRNA extension. Example: 10 15 20
    nicking_distance_minimum (int, optional): minimum nicking distance for designing ngRNAs. (Default: 0 bp)
    nicking_distance_maximum (int, optional): maximum nicking distance for designing ngRNAs. (Default: 100 bp)
    filter_c1_extension (bool, optional): filter against pegRNA extensions that start with a C base. (Default: False)
    silent_mutation (bool, optional): introduce silent mutation into PAM assuming sequence is in-frame. Currently only available with SpCas9. (Default: False)
    genome_wide_design (bool, optional): whether this is a genome-wide pooled design. This option designs a set of pegRNAs per input without ranging PBS and RTT parameters.
    saturation_mutagenesis (str, optional): saturation mutagenesis design with prime editing (Options: 'aa', 'aa_subs', 'aa_ins', 'aa_dels', 'base'). The 'aa' option makes all amino acid substitutions ('aa_subs'),  +1 amino acid insertions ('aa_ins'), and -1 amino acid deletions ('aa_dels'). The 'base' option makes DNA base changes. (Default: None)
    number_of_pegrnas (int, optional): maximum number of pegRNAs to design for each input sequence. The pegRNAs are ranked by 1) PAM disrupted > PAM intact then 2) distance to edit. (Default: 3)
    number_of_ngrnas (int, optional): maximum number of ngRNAs to design for each input sequence. The ngRNAs are ranked by 1) PE3b-seed > PE3b-nonseed > PE3 then 2) deviation from nicking_distance_pooled. (Default: 3)
    nicking_distance_pooled (int, optional): the nicking distance between pegRNAs and ngRNAs for pooled designs. PE3b annotation is priority (PE3b seed -> PE3b non-seed), followed by nicking distance closest to this parameter. (Default: 75 bp)
    homology_downstream (int, optional): this parameter determines the minimum RT extension length downstream of an edit for pegRNA designs. (Default: 15)
    pbs_length_pooled (int, optional): the PBS length to design pegRNAs for pooled design applications. (Default: 14 nt)
    rtt_max_length_pooled (int, optional): maximum RTT length to design pegRNAs for pooled design applications. (Default: 50 nt)
    out_dir (str, optional): name of output directory (Default: ./DATETIMESTAMP_PrimeDesign)
    
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

    Dependencies: os, numpy, & https://github.com/pinellolab/PrimeDesign
    '''
    # Write PrimeDesign Command Line
    cmd = 'python -m edms.bio.primedesign'
    cmd += f' -f {file}' # Append required parameters
    if pe_format!='NNNNNNNNNNNNNNNNN/NNN[NGG]': cmd += f' -pe_format {pe_format}'
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
    if homology_downstream!=15: cmd += f' -homology_downstream {homology_downstream}'
    if pbs_length_pooled!=14: cmd += f' -pbs_pooled {pbs_length_pooled}'
    if rtt_max_length_pooled!=50: cmd += f' -rtt_pooled {rtt_max_length_pooled}'
    if out_dir!='./DATETIMESTAMP_PrimeDesign': cmd+= f' -out {out_dir}'
    print(cmd)
    
    os.system(cmd) # Execute PrimeDesign Command Line

def prime_design_output(pt: str, scaffold_sequence: str, in_file: pd.DataFrame | str, saturation_mutagenesis:str=None, 
                        index: int=1, enzymes: list[str]=['Esp3I'], replace: bool=True) -> dict[pd.DataFrame]:
    ''' 
    prime_design_output(): splits peg/ngRNAs from PrimeDesign output & finishes annotations
    
    Parameters:
    pt (str): path to primeDesign output
    scaffold_sequence (str): sgRNA scaffold sequence
        SpCas9 flip + extend (shorter): GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC
        SpCas9 flip + extend + com-modified (required for VLPs): GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
    in_file (Dataframe | str): input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence (column names required)
    saturation_mutagenesis (str, optional): saturation mutagenesis design with prime editing (Options: 'aa', 'aa_subs', 'aa_ins', 'aa_dels', 'base'). The 'aa' option makes all amino acid substitutions ('aa_subs'),  +1 amino acid insertions ('aa_ins'), and -1 amino acid deletions ('aa_dels'). The 'base' option makes DNA base changes. (Default: 'aa')
    index (int, optional): 1st amino acid or base in target sequence index (Default: 1)
    enzymes (list, optional): list of type IIS RE enzymes (i.e., Esp3I, BsaI, BspMI) to check for in pegRNAs and ngRNAs (Default: ['Esp3I'])
    replace (bool, optional): replace pegRNAs and remove ngRNAs with RE enzyme sites (Default: True if saturation_mutagenesis is not None)
    
    Dependencies: io, numpy, & pandas
    '''
    if type(in_file) == str: # Get in_file from file path if needed
        in_file = io.get(pt=in_file)

    if saturation_mutagenesis is not None: # Saturation mutagenesis mode
        
        # Get target_name from input file
        target_name_in_file = in_file.iloc[0]['target_name']

        # Get PrimeDesign output & seperate pegRNAs and ngRNAs
        primeDesign_output = io.get(pt)
        pegRNAs = primeDesign_output[primeDesign_output['gRNA_type']=='pegRNA'].reset_index(drop=True)
        ngRNAs = primeDesign_output[primeDesign_output['gRNA_type']=='ngRNA'].reset_index(drop=True)

        # Generate pegRNAs
        pegRNAs['Edit']=[str(target_name.split('_')[-1].split('to')[0]) + # AA Before
                        str(int(target_name.split('_')[-2]) + index-1) + # AA Index
                        str(target_name.split('_')[-1].split('to')[1]) # AA After
                        for target_name in pegRNAs['Target_name']]
        pegRNAs['Target_name']=[target_name_in_file]*len(pegRNAs)
        pegRNAs['Scaffold_sequence']=[scaffold_sequence]*len(pegRNAs)
        pegRNAs['RTT_sequence']=[pegRNAs.iloc[i]['Extension_sequence'][0:int(pegRNAs.iloc[i]['RTT_length'])] for i in range(len(pegRNAs))]
        pegRNAs['PBS_sequence']=[pegRNAs.iloc[i]['Extension_sequence'][int(pegRNAs.iloc[i]['RTT_length']):].upper()  for i in range(len(pegRNAs))]
        if 'aa' in saturation_mutagenesis: 
            pegRNAs['AA_number'] = [int(re.findall(r'-?\d+',edit)[0]) if edit is not None else index for edit in pegRNAs['Edit']]
            cols=['pegRNA_number','gRNA_type','Strand','Edit','AA_number','Edit_type', # Important metadata
                'Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence',  # Sequence information
                'Target_name','Target_sequence','Spacer_GC_content','PAM_sequence','Extension_sequence','Annotation','pegRNA-to-edit_distance','Nick_index','PBS_length','PBS_GC_content','RTT_length','RTT_GC_content','First_extension_nucleotide','Reference_sequence', 'Edit_sequence', 'Silent_mutation_relative_to_edit'] # Less important metadata
        else: 
            pegRNAs['Base_number'] = [int(re.findall(r'-?\d+',edit)[0]) if edit is not None else index for edit in pegRNAs['Edit']]
            cols=['pegRNA_number','gRNA_type','Strand','Edit','Base_number','Edit_type', # Important metadata
                'Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence',  # Sequence information
                'Target_name','Target_sequence','Spacer_GC_content','PAM_sequence','Extension_sequence','Annotation','pegRNA-to-edit_distance','Nick_index','PBS_length','PBS_GC_content','RTT_length','RTT_GC_content','First_extension_nucleotide','Reference_sequence', 'Edit_sequence', 'Silent_mutation_relative_to_edit'] # Less important metadata
        for pegRNA_col in pegRNAs.columns: # Keep any additional columns in pegRNAs DataFrame
            if pegRNA_col not in cols and pegRNA_col not in ['ngRNA-to-pegRNA_distance','Spacer_sequence_order_TOP','Spacer_sequence_order_BOTTOM','pegRNA_extension_sequence_order_TOP','pegRNA_extension_sequence_order_BOTTOM']:
                cols.append(pegRNA_col)
        pegRNAs = t.reorder_cols(df=pegRNAs, cols=cols, keep=False) 
        
        # Generate ngRNAs
        ngRNAs['Edit']=[str(target_name.split('_')[-1].split('to')[0]) + # AA Before
                        str(int(target_name.split('_')[-2]) + index-1) + # AA Index
                        str(target_name.split('_')[-1].split('to')[1]) # AA After
                        for target_name in ngRNAs['Target_name']]
        ngRNAs['Target_name']=[target_name_in_file]*len(ngRNAs)
        ngRNAs['Scaffold_sequence']=[scaffold_sequence]*len(ngRNAs)
        ngRNAs['ngRNA_number']=list(np.arange(1,len(ngRNAs)+1))
        if 'aa' in saturation_mutagenesis: 
            ngRNAs['AA_number'] = [int(re.findall(r'-?\d+',edit)[0]) if edit is not None else index for edit in ngRNAs['Edit']]
            cols = ['pegRNA_number','ngRNA_number','gRNA_type','Strand','Edit','AA_number', # Important metadata
                'Spacer_sequence','Scaffold_sequence',  # Sequence information
                'Target_name','Target_sequence','Spacer_GC_content','PAM_sequence','Annotation','Nick_index','ngRNA-to-pegRNA_distance','Reference_sequence', 'Edit_sequence', 'Silent_mutation_relative_to_edit'] # Less important metadata
        else: 
            ngRNAs['Base_number'] = [int(re.findall(r'-?\d+',edit)[0]) if edit is not None else index for edit in ngRNAs['Edit']]
            cols = ['pegRNA_number','ngRNA_number','gRNA_type','Strand','Edit','Base_number', # Important metadata
                'Spacer_sequence','Scaffold_sequence',  # Sequence information
                'Target_name','Target_sequence','Spacer_GC_content','PAM_sequence','Annotation','Nick_index','ngRNA-to-pegRNA_distance','Reference_sequence', 'Edit_sequence', 'Silent_mutation_relative_to_edit'] # Less important metadata
        
        for ngRNA_col in ngRNAs.columns: # Keep any additional columns in ngRNAs DataFrame
            if ngRNA_col not in cols and ngRNA_col not in ['Spacer_sequence_order_TOP','Spacer_sequence_order_BOTTOM','pegRNA_extension_sequence_order_TOP','pegRNA_extension_sequence_order_BOTTOM']:
                cols.append(ngRNA_col)
        ngRNAs = t.reorder_cols(df=ngRNAs, cols=cols, keep=False)
    
    else: # Not saturation mutagenesis mode
        
        # Get PrimeDesign output & seperate pegRNAs and ngRNAs
        primeDesign_output = io.get(pt)
        pegRNAs = primeDesign_output[primeDesign_output['gRNA_type']=='pegRNA'].reset_index(drop=True)
        ngRNAs = primeDesign_output[primeDesign_output['gRNA_type']=='ngRNA'].reset_index(drop=True)

        # Generate pegRNAs
        pegRNAs['Scaffold_sequence']=[scaffold_sequence]*len(pegRNAs)
        pegRNAs['RTT_sequence']=[pegRNAs.iloc[i]['Extension_sequence'][0:int(pegRNAs.iloc[i]['RTT_length'])] for i in range(len(pegRNAs))]
        pegRNAs['PBS_sequence']=[pegRNAs.iloc[i]['Extension_sequence'][int(pegRNAs.iloc[i]['RTT_length']):].upper()  for i in range(len(pegRNAs))]
        cols = ['Target_name','pegRNA_number','gRNA_type','Strand','Edit_type', # Important metadata
                'Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence',  # Sequence information
                'Target_sequence','Spacer_GC_content','PAM_sequence','Extension_sequence','Annotation','pegRNA-to-edit_distance','Nick_index','PBS_length','PBS_GC_content','RTT_length','RTT_GC_content','First_extension_nucleotide','Reference_sequence', 'Edit_sequence', 'Silent_mutation_relative_to_edit'] # Less important metadata
        for pegRNA_col in pegRNAs.columns: # Keep any additional columns in pegRNAs DataFrame
            if pegRNA_col not in cols and pegRNA_col not in ['ngRNA-to-pegRNA_distance','Spacer_sequence_order_TOP','Spacer_sequence_order_BOTTOM','pegRNA_extension_sequence_order_TOP','pegRNA_extension_sequence_order_BOTTOM']:
                cols.append(pegRNA_col)
        pegRNAs = t.reorder_cols(df=pegRNAs, cols=cols, keep=False) 
        
        # Generate ngRNAs
        ngRNAs['Scaffold_sequence']=[scaffold_sequence]*len(ngRNAs)
        ngRNAs['ngRNA_number']=list(np.arange(1,len(ngRNAs)+1))
        cols = ['Target_name','ngRNA_number','gRNA_type','Strand', # Important metadata
                'Spacer_sequence','Scaffold_sequence',  # Sequence information
                'Target_sequence','Spacer_GC_content','PAM_sequence','Annotation','Nick_index','ngRNA-to-pegRNA_distance','Reference_sequence', 'Edit_sequence', 'Silent_mutation_relative_to_edit'] # Less important metadata
        for ngRNA_col in ngRNAs.columns: # Keep any additional columns in ngRNAs DataFrame
            if ngRNA_col not in cols and ngRNA_col not in ['Spacer_sequence_order_TOP','Spacer_sequence_order_BOTTOM','pegRNA_extension_sequence_order_TOP','pegRNA_extension_sequence_order_BOTTOM']:
                cols.append(ngRNA_col)
        ngRNAs = t.reorder_cols(df=ngRNAs, cols=cols, keep=False)
        
        # Set replace to False if not saturation mutagenesis
        replace = False
    
    # Temporarily make pegRNAs and ngRNAs oligonucleotides
    pegRNAs['Oligonucleotide'] = [str(spacer+scaffold+rtt+pbs).upper() for (spacer, scaffold, rtt, pbs) in t.zip_cols(df=pegRNAs,cols=['Spacer_sequence','Scaffold_sequence','RTT_sequence','PBS_sequence'])]
    ngRNAs['Oligonucleotide'] = [str(spacer+scaffold).upper() for (spacer, scaffold) in t.zip_cols(df=ngRNAs,cols=['Spacer_sequence','Scaffold_sequence'])]
    
    # Check for 0 recognition sites per enzyme
    for enzyme in enzymes:
        # pegRNAs: Find recognition sites for enzymes
        pegRNAs = find_enzyme_sites(df=pegRNAs, enzyme=enzyme)
        
        if replace: # Replace pegRNAs with RE enzyme sites
            
            pegRNAs_edits = list(pegRNAs['Edit'].unique()) # Get pegRNA edits
    
            # Store pegRNAs with recognition sites for enzymes
            pegRNAs_enzyme = pegRNAs[pegRNAs[enzyme]!=0]
            io.save(dir=f'../pegRNAs/{enzyme}/codon_swap_before',
                    file=f'{int(pegRNAs_enzyme.iloc[0]['PBS_length'])}.csv',
                    obj=pegRNAs_enzyme)
            
            # Codon swap pegRNAs with enzyme recognition site
            pegRNAs_enzyme = enzyme_codon_swap(pegRNAs=pegRNAs_enzyme,enzyme=enzyme,comments=True)
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

def prime_designer(in_file: str = None, target_name: str = None, flank5_sequence: str = None, target_sequence: str = None, 
                flank3_sequence: str = None, index: int=1, pe_format: str = 'NNNNNNNNNNNNNNNNN/NNN[NGG]', pbs_length_list: list = [],
                rtt_length_list: list = [], nicking_distance_minimum: int = 0,nicking_distance_maximum: int = 100, filter_c1_extension: bool = False,
                silent_mutation: bool = True, genome_wide_design: bool = False, saturation_mutagenesis: str = None,
                number_of_pegrnas: int = 3, number_of_ngrnas: int = 3, nicking_distance_pooled: int = 75, homology_downstream: int = 15,
                pbs_length_pooled_list: list = [11,13,15], rtt_max_length_pooled: int = 50,
                scaffold_sequence: str='GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC', 
                enzymes: list[str]=['Esp3I'], replace: bool=True):
    '''
    prime_designer(): execute PrimeDesign saturation mutagenesis (EDMS version)
    
    Parameters:
    in_file (str, required option 1): path to input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence,index (column names required). See examples below. (Default: None)
    target_name (str, required option 2): name of target
    flank5_sequence (str, required option 2): in-frame nucleotide sequence with 5' of saturation mutagensis region (length must be divisible by 3)
    target_sequence (str, required option 2): in-frame nucleotide sequence for the saturation mutagensis region (length must be divisible by 3)
    flank3_sequence (str, required option 2): in-frame nucleotide sequence with 3' of saturation mutagensis region (length must be divisible by 3)
    index (int, required option 2): 1st amino acid or base in target sequence index (Default: 1)
    
    pe_format (str, optional): Prime editing formatting including the spacer, cut index -> /, and protospacer adjacent motif (PAM) -> [PAM] (Default: NNNNNNNNNNNNNNNNN/NNN[NGG])
    pbs_length_list (list, optional): list of primer binding site (PBS) lengths for the pegRNA extension. Example: 12 13 14 15
    rtt_length_list (list, optional): list of reverse transcription (RT) template lengths for the pegRNA extension. Example: 10 15 20
    nicking_distance_minimum (int, optional): minimum nicking distance for designing ngRNAs. (Default: 0 bp)
    nicking_distance_maximum (int, optional): maximum nicking distance for designing ngRNAs. (Default: 100 bp)
    filter_c1_extension (bool, optional): filter against pegRNA extensions that start with a C base. (Default: False)
    silent_mutation (bool, optional): introduce silent mutation into PAM assuming sequence is in-frame. Currently only available with SpCas9. (Default: False)
    genome_wide_design (bool, optional): whether this is a genome-wide pooled design. This option designs a set of pegRNAs per input without ranging PBS and RTT parameters (Default: False).
    saturation_mutagenesis (str, optional): saturation mutagenesis design with prime editing (Options: 'aa', 'aa_subs', 'aa_ins', 'aa_dels', 'base'). The 'aa' option makes all amino acid substitutions ('aa_subs'),  +1 amino acid insertions ('aa_ins'), and -1 amino acid deletions ('aa_dels'). The 'base' option makes DNA base changes. (Default: None)
    number_of_pegrnas (int, optional): maximum number of pegRNAs to design for each input sequence. The pegRNAs are ranked by 1) PAM disrupted > PAM intact then 2) distance to edit. (Default: 3)
    number_of_ngrnas (int, optional): maximum number of ngRNAs to design for each input sequence. The ngRNAs are ranked by 1) PE3b-seed > PE3b-nonseed > PE3 then 2) deviation from nicking_distance_pooled. (Default: 3)
    nicking_distance_pooled (int, optional): the nicking distance between pegRNAs and ngRNAs for pooled designs. PE3b annotation is priority (PE3b seed -> PE3b non-seed), followed by nicking distance closest to this parameter. (Default: 75 bp)
    homology_downstream (int, optional): this parameter determines the minimum RT extension length downstream of an edit for pegRNA designs. (Default: 15)
    pbs_length_pooled_list (list, optional): List of PBS lengths to design pegRNAs for pooled design applications. (Default: [11, 13, 15])
    rtt_max_length_pooled (int, optional): maximum RTT length to design pegRNAs for pooled design applications. (Default: 50 nt)
    scaffold_sequence (str, optional): sgRNA scaffold sequence (Default: SpCas9 flip + extend = GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC)
        Alternative option for VLPs: SpCas9 flip + extend + com-modified = GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGGCTGAATGCCTGCGAGCATCCCACCCAAGTGGCACCGAGTCGGTGC
    enzymes (list, optional): list of type IIS RE enzymes (i.e., Esp3I, BsaI, BspMI) to check for in pegRNAs and ngRNAs (Default: ['Esp3I'])
    replace (bool, optional): replace pegRNAs and remove ngRNAs with RE sites (Default: True)

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
    
    Dependencies: prime_design_input(), prime_design(), & prime_design_output()
    '''
    if in_file is None: # Create PrimeDesign input file if needed
        prime_design_input(target_name=target_name, flank5_sequence=flank5_sequence, target_sequence=target_sequence, 
                        flank3_sequence=flank3_sequence, index=index, silent_mutation=silent_mutation, dir='.', file=f'{"_".join(target_name.split(" "))}.csv')

    elif (in_file is not None) & (silent_mutation == True): # Check that in_file target sequences are in-frame if silent_mutation is True
        in_file_df = io.get(pt=in_file)

        if saturation_mutagenesis is not None: # Saturation mutagenesis mode
            in_file_df_target_sequence = in_file_df.iloc[0]['target_sequence']
            flank5_sequence = in_file_df_target_sequence.split('(')[0]
            target_sequence = in_file_df_target_sequence.split('(')[1].split(')')[0]
            flank3_sequence = in_file_df_target_sequence.split(')')[1]

            if len(flank5_sequence)%3 != 0: raise(ValueError(f"Length of flank5_sequence ({len(flank5_sequence)}) must divisible by 3."))
            if len(target_sequence)%3 != 0: raise(ValueError(f"Length of target_sequence ({len(target_sequence)}) must divisible by 3."))
            if len(flank3_sequence)%3 != 0: raise(ValueError(f"Length of flank3_sequence ({len(flank3_sequence)}) must divisible by 3."))

        else: # Not saturation mutagenesis mode
            print("Warning: Manually verify that target sequences in 'in_file' are in-frame when 'silent_mutation' is set to True.")

    # Iterate through PBS lengths
    pegRNAs=dict()
    ngRNAs=dict()
    for pbs_length_pooled in pbs_length_pooled_list:

        if in_file is None:
            # Run PrimeDesign
            prime_design(file=f'{"_".join(target_name.split(" "))}.csv', silent_mutation=silent_mutation, saturation_mutagenesis=saturation_mutagenesis,
                        number_of_pegrnas=number_of_pegrnas, number_of_ngrnas=number_of_ngrnas, pbs_length_pooled=pbs_length_pooled, 
                        rtt_max_length_pooled=rtt_max_length_pooled, homology_downstream=homology_downstream, pe_format=pe_format,
                        pbs_length_list=pbs_length_list, rtt_length_list=rtt_length_list, nicking_distance_minimum=nicking_distance_minimum,
                        nicking_distance_maximum=nicking_distance_maximum, filter_c1_extension=filter_c1_extension, genome_wide_design=genome_wide_design,
                        nicking_distance_pooled=nicking_distance_pooled)

            # Obtain pegRNAs and ngRNAs from PrimeDesign output
            pegRNAs[pbs_length_pooled],ngRNAs[pbs_length_pooled] = prime_design_output(
                pt=sorted([file for file in io.relative_paths('.') if "PrimeDesign.csv" in file], reverse= True)[0], 
                scaffold_sequence=scaffold_sequence, in_file=f'./{"_".join(target_name.split(" "))}.csv', 
                saturation_mutagenesis=saturation_mutagenesis, index=index, enzymes=enzymes, replace=replace)
        
        else:
            # Run PrimeDesign
            prime_design(file=in_file, silent_mutation=silent_mutation, saturation_mutagenesis=saturation_mutagenesis,
                        number_of_pegrnas=number_of_pegrnas, number_of_ngrnas=number_of_ngrnas, pbs_length_pooled=pbs_length_pooled, 
                        rtt_max_length_pooled=rtt_max_length_pooled, homology_downstream=homology_downstream, pe_format=pe_format,
                        pbs_length_list=pbs_length_list, rtt_length_list=rtt_length_list, nicking_distance_minimum=nicking_distance_minimum,
                        nicking_distance_maximum=nicking_distance_maximum, filter_c1_extension=filter_c1_extension, genome_wide_design=genome_wide_design,
                        nicking_distance_pooled=nicking_distance_pooled)
            
            # Obtain pegRNAs and ngRNAs from PrimeDesign output
            pegRNAs[pbs_length_pooled],ngRNAs[pbs_length_pooled] = prime_design_output(
                pt=sorted([file for file in io.relative_paths('.') if "PrimeDesign.csv" in file], reverse= True)[0], 
                scaffold_sequence=scaffold_sequence, in_file=in_file, 
                saturation_mutagenesis=saturation_mutagenesis, index=index, enzymes=enzymes, replace=replace)
        
        if (saturation_mutagenesis is None) & (genome_wide_design==False): # Only run once if not saturation mutagenesis
            pegRNAs['pegRNAs'] = pegRNAs.pop(pbs_length_pooled) # Rename dictionary keys
            ngRNAs['ngRNAs'] = ngRNAs.pop(pbs_length_pooled)
            break 
    
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

def sensor_designer(pegRNAs: pd.DataFrame | str, sensor_length: int=60, before_spacer: int=5, sensor_orientation: Literal['revcom','forward']='revcom',
                    out_dir: str=None, out_file: str=None, return_df: bool=True, literal_eval: bool=True) -> pd.DataFrame:
    ''' 
    sensor_designer(): design pegRNA sensors
    
    Parameters:
    pegRNAs (dataframe | str): pegRNAs DataFrame (or file path)
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

    # Get pegRNAs from file path if needed
    if type(pegRNAs)==str:
        pegRNAs = io.get(pt=pegRNAs, literal_eval=literal_eval)

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
    for spacer,strand,reference_sequence in t.zip_cols(df=pegRNAs,cols=['Spacer_sequence','Strand','Reference_sequence']): # Iterate through spacers
        
        # Get reverse complement and codons of reference sequence
        reference_sequence = Seq(reference_sequence)
        rc_reference_sequence = Seq.reverse_complement(reference_sequence)

        if strand=='+': # Spacer: + strand; PBS & RTT: - strand
            
            # Find spacer in sequence; compute spacer5 index
            spacer5 = reference_sequence.find(spacer)
            if spacer5 == -1:
                raise ValueError(f"Spacer sequence '{spacer}' not found in target sequence. Please check the input file.")
            elif spacer5 != reference_sequence.rfind(spacer):
                raise ValueError(f"Multiple matches found for spacer sequence '{spacer}'. Please check the input file.")

            # Assign start & end index for sensor
            start = spacer5 - before_spacer
            end = start + sensor_length
            sensor = reference_sequence[start:end]

        elif strand=='-': # Spacer: - strand; PBS & RTT: + strand
            
            # Find spacer in sequence; compute center index
            spacer5 = rc_reference_sequence.find(spacer)
            if spacer5 == -1:
                raise ValueError(f"Spacer sequence '{spacer}' not found in reference sequence.")
            if spacer5 != rc_reference_sequence.rfind(spacer):
                raise ValueError(f"Multiple matches found for spacer sequence '{spacer}' not found in reference sequence.")
            
            # Assign start & end index for sensor
            start = spacer5 - before_spacer
            end = start + sensor_length
            sensor = rc_reference_sequence[start:end]
            
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

def pegRNA_outcome(pegRNAs: pd.DataFrame | str, in_file: pd.DataFrame | str = None,
                   match_score: float = 2, mismatch_score: float = -1, open_gap_score: float = -10, extend_gap_score: float = -0.1,
                   out_dir: str=None, out_file: str=None, return_df: bool=True, literal_eval: bool=True) -> pd.DataFrame:
    ''' 
    pegRNA_outcome(): confirm that pegRNAs should create the predicted edits
    
    Parameters:
    pegRNAs (dataframe | str): pegRNAs DataFrame (or file path)
    in_file (dataframe | str): Input file (.txt or .csv) with sequences for PrimeDesign. Format: target_name,target_sequence,index (column names required)
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

    # Verify all expected edits are present in pegRNA library (for saturation mutagenesis only)
    if in_file is not None:

        # Get reference sequence & codons (+ reverse complement)
        target_sequence = in_file.iloc[0]['target_sequence'] 
        seq = Seq(target_sequence.split('(')[1].split(')')[0]) # Break apart target sequences
        if len(seq)%3 != 0: raise(ValueError(f"Length of target sequence ({len(seq)}) must divisible by 3. Check input file."))
        flank5 = Seq(target_sequence.split('(')[0])
        if len(flank5)%3 != 0: raise(ValueError(f"Length of flank5 ({len(flank5)}) must divisible by 3. Check input file."))
        flank3 = Seq(target_sequence.split(')')[1])
        if len(flank3)%3 != 0: raise(ValueError(f"Length of flank3 ({len(flank3)}) must divisible by 3. Check input file."))

        index = in_file.iloc[0]['index']

        f5_seq_f3_nuc = flank5 + seq + flank3  # Join full nucleotide reference sequence
        rc_f5_seq_f3_nuc = Seq.reverse_complement(f5_seq_f3_nuc) # Full nucleotide reference reverse complement sequence
        seq_prot = Seq.translate(seq) # In-frame amino acid sequence
        f5_seq_f3_prot = Seq.translate(f5_seq_f3_nuc) # Full in-frame protein sequence (including flanks)
        
        indexes = list(np.arange(index,index+len(seq_prot))) # In-frame amino acid indexes
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
            edits_substitutions.extend([f'{aa}{str(indexes[i])}{aa2}' for aa2 in aa_dna_codon_table if (aa2!='*')&(aa2!=aa)])
            edits_insertions.extend([f'{aa}{str(indexes[i])}{aa}{aa2}' for aa2 in aa_dna_codon_table if aa2!='*'])
            edits_deletions.append(f'{aa}{seq_prot_deletions[i+1]}{str(indexes[i])}{seq_prot_deletions[i+1]}')
        edits_substitutions_set = set(edits_substitutions)
        edits_insertions_set = set(edits_insertions)
        edits_deletions_set = set(edits_deletions)
        
        print(f'Expected Edits in pegRNA library...\nSubstitutions: {edits_substitutions}\nInsertions: {edits_insertions}\nDeletions: {edits_deletions}')
        print(f'All substitutions present: {edits_substitutions_set.issubset(set(pegRNAs["Edit"]))}; missing: {edits_substitutions_set-set(pegRNAs["Edit"])}')
        print(f'All insertions present: {edits_insertions_set.issubset(set(pegRNAs["Edit"]))}; missing: {edits_insertions_set-set(pegRNAs["Edit"])}')
        print(f'All deletions present: {edits_deletions_set.issubset(set(pegRNAs["Edit"]))}; missing: {edits_deletions_set-set(pegRNAs["Edit"])}\n')

    # Determine post_RTT_sequences
    post_RTT_sequences = [] # Store post RTT sequences
    for (strand,pbs,rtt,edit,aa_number,reference_sequence) in t.zip_cols(df=pegRNAs,cols=['Strand','PBS_sequence','RTT_sequence','Edit','AA_number','Reference_sequence']): # Iterate through primer binding sites

        # Get reverse complement and codons of reference sequence
        reference_sequence = Seq(reference_sequence)
        rc_reference_sequence = Seq.reverse_complement(reference_sequence)

        if strand=='+': # Spacer: + strand; PBS & RTT: - strand

            # Find reverse complement PBS in sequence
            rc_pbs = Seq.reverse_complement(Seq(pbs)) # reverse complement of pbs (+ strand)
            
            rc_pbs_j = reference_sequence.find(str(rc_pbs))
            if rc_pbs_j == -1:
                raise ValueError(f"PBS sequence '{pbs}' not found in reference sequence.")
            elif rc_pbs_j != reference_sequence.rfind(str(rc_pbs)):
                print(rc_pbs,rc_pbs_j,rc_reference_sequence.rfind(str(rc_pbs)))
                raise ValueError(f"Multiple matches found for PBS sequence '{str(rc_pbs)}'.")

            # Replace change sequence using reverse complement RTT
            if len(edit)==len(str(aa_number))+2: # Substitution
                post_RTT_sequences.append(str(reference_sequence[:rc_pbs_j+len(rc_pbs)]+
                                            Seq.reverse_complement(Seq(rtt.upper()))+
                                            reference_sequence[rc_pbs_j+len(rc_pbs)+len(rtt):])) # Save post RTT sequence
            
            elif edit.find(str(aa_number))>len(edit)-edit.find(str(aa_number))-len(str(aa_number)): # Deletion
                post_RTT_sequences.append(str(reference_sequence[:rc_pbs_j+len(rc_pbs)]+
                                            Seq.reverse_complement(Seq(rtt.upper()))+
                                            reference_sequence[rc_pbs_j+len(rc_pbs)+len(rtt)+3*(edit.find(str(aa_number))-1):]))
            
            else: # Insertion
                post_RTT_sequences.append(str(reference_sequence[:rc_pbs_j+len(rc_pbs)]+
                                            Seq.reverse_complement(Seq(rtt.upper()))+
                                            reference_sequence[rc_pbs_j+len(rc_pbs)+len(rtt)-3*(len(edit)-edit.find(str(aa_number))-len(str(aa_number))-1):]))
            
        elif strand=='-': # Spacer: - strand; PBS & RTT: + strand

            # Find PBS in sequence
            pbs = Seq(pbs) # pbs (+ strand)
            pbs_j = reference_sequence.find(str(pbs))
            if pbs_j == -1:
                raise ValueError(f"PBS sequence '{pbs}' not found in reference sequence.")
            elif pbs_j != reference_sequence.rfind(str(pbs)):
                print(pbs,pbs_j,reference_sequence.rfind(str(pbs)))
                raise ValueError(f"Multiple matches found for PBS sequence '{pbs}'.")
            
            # Replace change sequence using RTT
            if len(edit)==len(str(aa_number))+2: # Substitution
                post_RTT_sequences.append(str(reference_sequence[:pbs_j-len(rtt)]+
                                              Seq(rtt.upper())+
                                              reference_sequence[pbs_j:]))
            elif edit.find(str(aa_number))>len(edit)-edit.find(str(aa_number))-len(str(aa_number)): # Deletion
                post_RTT_sequences.append(str(reference_sequence[:pbs_j-len(rtt)-3*(edit.find(str(aa_number))-1)]+
                                            Seq(rtt.upper())+
                                            reference_sequence[pbs_j:]))
                
            else: # Insertion
                post_RTT_sequences.append(str(reference_sequence[:pbs_j-len(rtt)+3*(len(edit)-edit.find(str(aa_number))-len(str(aa_number))-1)]+
                                              Seq(rtt.upper())+
                                              reference_sequence[pbs_j:]))

        else: 
            raise(ValueError('Error: Strand column can only have "+" and "-".'))
        
    # Determine edit from post RTT sequences
    pegRNAs['Post_RTT_sequence']=post_RTT_sequences
    
    # Check edits & assign multiple edit annotations if needed
    edit_check = []
    for post_RTT_sequence,reference_sequence in t.zip_cols(df=pegRNAs,cols=['Post_RTT_sequence','Reference_sequence']):
        if len(reference_sequence)!=len(post_RTT_sequence): # Indel
            edit = fq.find_indel(wt=reference_sequence, mut=post_RTT_sequence, res=int(index-len(flank5)/3), show=False, 
                                 match_score=match_score, mismatch_score=mismatch_score,
                                 open_gap_score=open_gap_score, extend_gap_score=extend_gap_score)[0]
            
            # Check for additional edits
            aa_number = int(re.findall(r'-?\d+',edit)[0])
            i = int(aa_number-index+len(flank5)/3)
            if edit.find(str(aa_number))>len(edit)-edit.find(str(aa_number))-len(str(aa_number)): # Deletion
                
                # Look for next AA(s) that match the deleted AA in the edit
                i_ls = [i]
                for j,aa_j in enumerate(f5_seq_f3_prot[i+1:]):
                    if aa_j==edit[0]:
                        i_ls.append(i+j+1)
                    else:
                        break
                
                if len(i_ls)>1: # Multiple edit annotations
                    edits = [f'{f5_seq_f3_prot[i]}{f5_seq_f3_prot[i+1]}{int(i+index-len(flank5)/3)}{f5_seq_f3_prot[i+1]}' for i in i_ls]
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
                    edits = [f'{f5_seq_f3_prot[i]}{int(i+index-len(flank5)/3)}{f5_seq_f3_prot[i]}{edit[-1]}' for i in i_ls]
                    edit_check.append(edits)
                else: 
                    edit_check.append(edit)
            
        else: # Substitution
            edit_check.append(fq.find_AA_edits(wt=str(Seq.translate(reference_sequence)), 
                                               res=int(index-len(flank5)/3), 
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

def pegRNA_signature(pegRNAs: pd.DataFrame | str, flank5_sequence: str | Seq, flank3_sequence: str | Seq, flank_length: int=15,
                    reference_sequence: str='Reference_sequence', edit_sequence: str='Edit_sequence',
                    match_score: float = 2, mismatch_score: float = -1, open_gap_score: float = -10, extend_gap_score: float = -0.1,
                    out_dir: str=None, out_file: str=None, save_alignments: bool=False, return_df: bool=True, literal_eval: bool=True) -> pd.DataFrame:
    ''' 
    pegRNA_signature(): create signatures for pegRNA outcomes using alignments
    
    Parameters:
    pegRNAs (dataframe | str): pegRNAs DataFrame (or file path)
    flank5_sequence (str | Seq): flank5 sequence
    flank3_sequence (str | Seq): flank3 sequence
    flank_length (int, optional): length of flank sequences to include in alignment (Default: 15)
    reference_sequence (str, optional): column name for reference sequences (Default: 'Reference_sequence')
    edit_sequence (str, optional): column name for post RTT sequences (Default: 'Edit_sequence')
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

    # Get pegRNAs DataFrame from file path if needed
    if type(pegRNAs)==str:
        pegRNAs = io.get(pt=pegRNAs,literal_eval=literal_eval)

    # High sequence homology; punish gaps
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = match_score  # Score for a match
    aligner.mismatch_score = mismatch_score  # Penalty for a mismatch; applied to both strands
    aligner.open_gap_score = open_gap_score  # Penalty for opening a gap; applied to both strands
    aligner.extend_gap_score = extend_gap_score  # Penalty for extending a gap; applied to both strands

    # Get alignments and signatures
    aligments_list = []
    signatures_list = []
    for i,(reference_seq, edit_seq) in enumerate(t.zip_cols(df=pegRNAs, cols=[reference_sequence,edit_sequence])): # Iterate through reference and edit sequences
        
        # Convert Seq objects to strings
        flank5_sequence = str(flank5_sequence)
        flank3_sequence = str(flank3_sequence)
        reference_seq = str(reference_seq)
        edit_seq = str(edit_seq)

        # Trim flanks from reference and edit sequences
        if reference_seq.find(flank5_sequence)==-1 or reference_seq.rfind(flank3_sequence)==-1:
            raise(ValueError(f"Flank5 or Flank3 sequences were not found in reference sequence for pegRNAs row {i}.\nPlease check the flank5 ({flank5_sequence}) and flank3 ({flank3_sequence}) sequences.\nReference sequence: {reference_seq}"))
        if edit_seq.find(flank5_sequence)==-1 or edit_seq.rfind(flank3_sequence)==-1:
            raise(ValueError(f"Flank5 or Flank3 sequences were not found in edit sequence for pegRNAs row {i}.\nPlease check the flank5 ({flank5_sequence}) and flank3 ({flank3_sequence}) sequences.\nEdit sequence: {edit_seq}"))
        reference_seq = reference_seq[reference_seq.find(flank5_sequence)+len(flank5_sequence)-flank_length : reference_seq.rfind(flank3_sequence)+flank_length]
        edit_seq = edit_seq[edit_seq.find(flank5_sequence)+len(flank5_sequence)-flank_length : edit_seq.rfind(flank3_sequence)+flank_length]

        # Create and append alignment
        alignment = aligner.align(Seq(reference_seq), Seq(edit_seq))[0]
        aligments_list.append(alignment)

        # Create and append signature
        signatures_list.append(signature_from_alignment(ref_seq=reference_seq, query_seq=edit_seq, alignment=alignment))

    # Create Alignment and Signature columns
    pegRNAs['Alignment'] = aligments_list
    pegRNAs['Signature'] = signatures_list
    
    # Drop Alignment column if not saving
    if save_alignments==False:
        pegRNAs.drop(columns=['Alignment'],inplace=True)

    # Count # of SNVs, insertions, deletions in Signature
    snvs_ls = []
    ins_ls = []
    dels_ls = []
    for sign in pegRNAs['Signature']: # Iterate through signatures

        # Count SNVs, insertions, deletions
        snvs_ls.append(len([snv for snv in sign.snvs]))
        ins_ls.append(sum([len(ind.ins) for ind in sign.indels]))
        dels_ls.append(sum([ind.dellen for ind in sign.indels]))

    # Create count columns
    pegRNAs['SNV_count'] = snvs_ls
    pegRNAs['ins_count'] = ins_ls
    pegRNAs['del_count'] = dels_ls
    pegRNAs['difference_count'] = pegRNAs['SNV_count'] + pegRNAs['ins_count'] + pegRNAs['del_count']

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