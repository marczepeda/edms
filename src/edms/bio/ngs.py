'''
Module: ngs.py
Author: Marc Zepeda
Created: 2024-11-12
Description: Next generation sequencing

Usage:
[NGS Thermocycler]
- group_boundaries(): Group a list of integers into segments where each segment contains consecutive numbers.
- min_sec(): Convert decimal minutes to a tuple of (minutes, seconds).
- thermocycler(): Creates a dictionary of thermocycler objects from a DataFrame.

[NGS PCR calculation]
- pcr_mm(): NEB Q5 PCR master mix calculations
- pcr_mm_ultra(): NEBNext Ultra II Q5 PCR master mix calculations
- pcrs(): generates NGS PCR plan automatically
- umis(): Returns the ug of gDNA, number of molecules, and reads needed for a given number of genotypes and samples.

[Hamming distance calculation]
- hamming_distance(): returns the Hamming distance between two sequences
- compute_distance_matrix(): returns pairwise Hamming distance matrix for a list of sequences
''' 
# Import Packages
import math
from typing import Literal
import numpy as np
import pandas as pd
from Bio.Seq import Seq
import os
from ..gen import io
from ..gen import tidy as t

# NGS Thermocycler
def group_boundaries(nums: list[int]) -> list[tuple[int, int]]:
    '''
    group_boundaries(nums): Group a list of integers into segments where each segment contains consecutive numbers.
    
    Parameters:
    nums (list[int]): A list of integers to be grouped.
    '''
    if not nums: # Check if the list is empty
        return []

    nums = sorted(set(nums))  # Remove duplicates and sort the list
    groups = [] # Initialize an empty list to hold the groups
    start = nums[0]

    for i in range(1, len(nums)): # Iterate through the list starting from the second element
        if nums[i] - nums[i - 1] > 1:
            groups.append((start, nums[i - 1]))
            start = nums[i]

    groups.append((start, nums[-1]))  # Close the final segment
    return groups

def min_sec(decimal_minutes: float) -> tuple[int, int]:
    '''
    min_sec(): Convert decimal minutes to a tuple of (minutes, seconds).

    Parameters:
    decimal_minutes (float): Decimal representation of minutes.
    '''
    minutes = int(decimal_minutes)
    seconds = int(round((decimal_minutes - minutes) * 60))
    return minutes, seconds

def thermocycler(df: pd.DataFrame, n: Literal['0','1.5','2'] = '1', cycles: int | str = None, 
                 pcr_fwd_col: str=None, pcr_rev_col: str=None, umi: bool=False) -> dict[pd.DataFrame]:
    """
    thermocycler(): Creates a dictionary of thermocycler objects from a DataFrame.

    Parameters:
    df (pd.DataFrame): DataFrame containing PCR data.
    n (Literal[1, 2, 3]): The PCR number to process (Default: 1).
    cycles (int | str): Number of cycles for PCR1 (Default: None -> 30).
    pcr_fwd_col (str, optional): PCR FWD column name (Default: None -> f'PCR{n} FWD').
    pcr_rev_col (str, optional): PCR REV column name (Default: None -> f'PCR{n} REV').
    
    Dependencies: math, pandas, typing.Literal, tidy, min_sec(), group_boundaries()
    """
    dc = dict() # Initialize an empty dictionary to hold thermocycler dataframes
    
    # Set cycles and primer columns based on the PCR number
    if umi==False: # Non-UMI
        if str(n)=='1': # PCR1
            if cycles is None:
                cycles = 30
            if pcr_fwd_col is None:
                pcr_fwd_col = 'PCR1 FWD'
            if pcr_rev_col is None:
                pcr_rev_col = 'PCR1 REV'
        elif str(n)=='2': # PCR2
            if cycles is None:
                cycles = 8
            if pcr_fwd_col is None:
                pcr_fwd_col = 'PCR2 FWD'
            if pcr_rev_col is None:
                pcr_rev_col = 'PCR2 REV'
            df[pcr_fwd_col] = 'PCR2-FWD'
        else: 
            raise ValueError("n must be 1 or 2")
        
    else: # UMI
        if str(n)=='1': # PCR1 (UMI PCR linear amplification)
            if cycles is None:
                cycles = 3
            if pcr_fwd_col is None:
                pcr_fwd_col = 'PCR1 FWD'
            if pcr_rev_col is None:
                pcr_rev_col = 'PCR1 REV'
        elif str(n)=='1.5': # PCR1.5 (UMI PCR exponential apmlification)
            if cycles is None:
                cycles = 30
            if pcr_fwd_col is None:
                pcr_fwd_col = 'PCR1 FWD'
            if pcr_rev_col is None:
                pcr_rev_col = 'PCR1 REV'
            df[pcr_fwd_col] = 'P5-FWD'
            df[pcr_rev_col] = 'P7-REV'
        elif str(n)=='2': # PCR2
            if cycles is None:
                cycles = 8
            if pcr_fwd_col is None:
                pcr_fwd_col = 'PCR2 FWD'
            if pcr_rev_col is None:
                pcr_rev_col = 'PCR2 REV'
            df[pcr_fwd_col] = 'PCR2-FWD' 
        else: 
            raise ValueError("n must be 1, 1.5 or 2")
        
    
    # Iterate through unique combinations of FWD and REV primers for the specified PCR number
    for (fwd,rev) in t.unique_tuples(df=df, cols=[pcr_fwd_col,pcr_rev_col]):
        title = f"{fwd}_{rev}: " # Initialize the thermocycler title for the set of primers

        # Group the IDs based on the set of primers; determine the temperature and anneal_time
        ls = group_boundaries(list(df[(df[pcr_fwd_col]==fwd) & (df[pcr_rev_col]==rev)]['ID'].keys()))
        tm = df[(df[pcr_fwd_col]==fwd) & (df[pcr_rev_col]==rev)].iloc[0][f'PCR{n} Tm']
        bp = df[(df[pcr_fwd_col]==fwd) & (df[pcr_rev_col]==rev)].iloc[0]['PCR2 bp']
        (min,sec) = min_sec(math.floor(bp/500)/2+0.5)
        if min == 0:
            anneal_time = f"{sec}s"
        else:
            if sec == 0:
                anneal_time = f"{min}min"
            else:
                anneal_time = f"{min}min {sec}s"
        
        # Append grouped IDs to the thermocycler title
        for (start, end) in ls:
            if start == end:
                title += f"{df.iloc[start][f'PCR{n} ID']}, "
            else:
                title += f"{df.iloc[start][f'PCR{n} ID']} -> {df.iloc[end][f'PCR{n} ID']}, "
        
        # Create a DataFrame for the thermocycler steps
        dc[f'{fwd}_{rev}_{tm}°C'] = pd.DataFrame(
            {'Temperature':['98°C', '98°C', f'{tm}°C', '72°C', '72°C', '4°C', ''],
             'Time':['30s', '10s', '30s', anneal_time, '2min', '∞', ''],
             'Repeat':['', f'{cycles} cycles', f'{cycles} cycles', f'{cycles} cycles', '', '', '']},
             index=pd.Index(['1','2','2','2','3','4',f'{title[:-2]}'], name="Step"))

    return dc

# NGS PCR calculation  
def pcr_mm(primers: pd.Series,  template: str, template_uL: int,
           Q5_mm_x_stock: int=5, dNTP_mM_stock: int=10, fwd_uM_stock: int=10, rev_uM_stock: int=10, Q5_U_uL_stock: int=2,
           Q5_mm_x_desired: int=1, dNTP_mM_desired: float=0.2,fwd_uM_desired: float=0.5, rev_uM_desired: float=0.5, Q5_U_uL_desired: float=0.02,
           total_uL: int=25, mm_x: float=1.1) -> dict[pd.DataFrame]:
    '''
    pcr_mm(): NEB Q5 PCR master mix calculations
    
    Parameters:
    primers (Series): (FWD primer, REV primer) pairs
    template (str): template name
    template_uL (int): template uL per reaction
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
    total_uL (int, optional): total uL per reaction (Default: 20)
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
                                            
def pcr_mm_ultra(primers: pd.Series, template: str, template_uL: int,
                 Q5_mm_x_stock: int=2, fwd_uM_stock: int=10, rev_uM_stock: int=10,
                 Q5_mm_x_desired: int=1,fwd_uM_desired: float=0.5, rev_uM_desired: float=0.5,
                 total_uL: int=20, mm_x: float=1.1) -> dict[pd.DataFrame]:
    '''
    pcr_mm_ultra(): NEBNext Ultra II Q5 PCR master mix calculations
    
    Parameters:
    primers (Series): (FWD primer, REV primer) pairs
    template (str): template name
    template_uL (int): template uL per reaction
    Q5_mm_x_stock (int, optional): Q5 reaction master mix stock (Default: 2)
    fwd_uM_stock (int, optional): [FWD Primer] stock in mM (Default: 10)
    rev_uM_stock (int, optional): [REV Primer] stock in mM (Default: 10)
    Q5_mm_x_desired (int, optional): Q5 reaction master mix desired (Default: 1)
    fwd_uM_desired (float, optional): [FWD Primer] desired in mM (Default: 0.5)
    rev_uM_desired (float, optional): [REV Primer] desired in mM (Default: 0.5)
    total_uL (int, optional): total uL per reaction (Default: 20)
    mm_x (float, optional): master mix multiplier (Default: 1.1)

    Dependencies: pandas
    '''
    pcr_mm_dc = dict()
    for i,(pcr1_fwd,pcr1_rev) in enumerate(primers.keys()):
        pcr_mm_dc[(pcr1_fwd,pcr1_rev)] = pd.DataFrame({'Component':['Nuclease-free H2O',f'NEBNext Ultra II Q5 {Q5_mm_x_stock}x MM',pcr1_fwd,pcr1_rev,template,'Total'],
                                                       'Stock':['',Q5_mm_x_stock,fwd_uM_stock,rev_uM_stock,'',''],
                                                       'Desired':['',Q5_mm_x_desired,fwd_uM_desired,rev_uM_desired,'',''],
                                                       'Unit':['','x','uM','uM','',''],
                                                       'uL': [round(total_uL-sum([Q5_mm_x_desired/Q5_mm_x_stock,fwd_uM_desired/fwd_uM_stock,rev_uM_desired/rev_uM_stock,template_uL/total_uL]*total_uL),2),
                                                              round(Q5_mm_x_desired/Q5_mm_x_stock*total_uL,2),
                                                              round(fwd_uM_desired/fwd_uM_stock*total_uL,2),
                                                              round(rev_uM_desired/rev_uM_stock*total_uL,2),
                                                              round(template_uL,2),
                                                              round(total_uL,2)],
                                                       'uL MM': [round((total_uL-sum([Q5_mm_x_desired/Q5_mm_x_stock,fwd_uM_desired/fwd_uM_stock,rev_uM_desired/rev_uM_stock,template_uL/total_uL]*total_uL))*primers.iloc[i]*mm_x,2),
                                                                 round(Q5_mm_x_desired/Q5_mm_x_stock*total_uL*primers.iloc[i]*mm_x,2),
                                                                 round(fwd_uM_desired/fwd_uM_stock*total_uL*primers.iloc[i]*mm_x,2),
                                                                 round(rev_uM_desired/rev_uM_stock*total_uL*primers.iloc[i]*mm_x,2),
                                                                 round(template_uL*primers.iloc[i]*mm_x,2),
                                                                 round(total_uL*primers.iloc[i]*mm_x,2)]
                                                     },index=pd.Index(list(np.arange(1,7)), name=f"{pcr1_fwd}_{pcr1_rev}"))
    return pcr_mm_dc

def pcrs(df: pd.DataFrame | str, dir:str=None, file:str=None, gDNA_id_col: str='ID', 
         pcr1_id_col: str='PCR1 ID', pcr1_fwd_col: str='PCR1 FWD', pcr1_rev_col: str='PCR1 REV', 
         pcr2_id_col: str='PCR2 ID', pcr2_fwd_col: str='PCR2 FWD', pcr2_rev_col: str='PCR2 REV', umi_col: str='UMI',
         Q5_mm_x_stock: int=5, dNTP_mM_stock: int=10, fwd_uM_stock: int=10, rev_uM_stock: int=10, Q5_U_uL_stock: int=2,
         Q5_mm_x_desired: int=1,dNTP_mM_desired: float=0.2, fwd_uM_desired: float=0.5, rev_uM_desired: float=0.5, Q5_U_uL_desired: float=0.02,
         pcr1_total_uL: int=20, pcr2_total_uL: int=20, mm_x: float=1.1, ultra: bool=False, 
         pcr1_cycles: int | str = None, pcr2_cycles: int | str = None, umi_cycles: int | str = None, pcr1_5_Tm: int | str = None) -> tuple[dict[pd.DataFrame]]:
    '''
    pcrs(): generates NGS PCR plan automatically
    
    Parameters:
    df (DataFrame | str): NGS samples dataframe (or file path)
    dir (str, optional): save directory
    file (str, optional): save file
    gDNA_id_col (str, optional): gDNA ID column name (Default: 'ID')
    pcr1_id_col (str, optional): PCR1 ID column name (Default: 'PCR1 ID')
    pcr1_fwd_col (str, optional): PCR1 FWD column name (Default: 'PCR1 FWD')
    pcr1_rev_col (str, optional): PCR1 REV column name (Default: 'PCR1 REV')
    pcr1_id_col (str, optional): PCR2 ID column name (Default: 'PCR2 ID')
    pcr1_fwd_col (str, optional): PCR2 FWD column name (Default: 'PCR2 FWD')
    pcr1_rev_col (str, optional): PCR2 REV column name (Default: 'PCR2 REV')
    umi_col (str, optional): UMI column name (Default: 'UMI')
    template_uL (int): template uL per reaction
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
    pcr1_total_uL (int, optional): total uL per reaction (Default: 20)
    pcr2_total_uL (int, optional): total uL per reaction (Default: 20)
    mm_x (float, optional): master mix multiplier (Default: 1.1)
    ultra (bool, optional): use NEB Ultra II reagents (Default: False)
    pcr1_cycles (int | str): Number of cycles for the PCR1 process (Default: None -> 30).
    pcr2_cycles (int | str): Number of cycles for the PCR2 process (Default: None -> 8).
    umi_cycles (int | str): Number of cycles for the UMI PCR process (Default: None -> 3).
    pcr1_5_Tm (int | str): Annealing temperature for PCR1.5 process (Default: None -> 65; P5-FWD & P7-REV primers).

    Dependencies: pandas,numpy,os,io,tidy,thermocycler(),pcr_mm(),pcr_mm_ultra()
    '''
    # Get samples dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)

    # Define plate axis method
    def plate(df: pd.DataFrame, group: Literal['96-well plate','8-strip_plate'], pcr: Literal['1','2']) -> pd.DataFrame:
        ''' 
        plate(): Creates a DataFrame representing a 96-well plate or 8-strip plate layout.

        Parameters:
        df (pd.DataFrame): DataFrame containing PCR data.
        group (Literal['96-well plate','8-strip_plate']): Type of plate layout to create.
        pcr (Literal['1','2']): The PCR number to process.
        '''
        if group == '96-well plate': # Define 96-well plate axis
            rows = ['A','B','C','D','E','F','G','H']
            cols = np.arange(1,13,1)
        elif group == '8-strip plate': # Define PCR strip axis
            rows = ['A','B','C','D','E','F']
            cols = np.arange(1,9,1)
        else:
            raise ValueError("group must be '96-well plate' or '8-strip plate'")

        if str(pcr)=='1':
            
            # Store gDNA and PCR locations on plate
            ls_df_fwd_rev = []
            primers_ls = []
            plate_ls = []
            row_ls = []
            col_ls = []

            plate_i = 1
            for (fwd,rev) in t.unique_tuples(df=df,cols=[pcr1_fwd_col,pcr1_rev_col]):
                df_fwd_rev = df[(df[pcr1_fwd_col]==fwd) & (df[pcr1_rev_col]==rev)]
                ls_df_fwd_rev.append(df_fwd_rev)
                row_i = 0
                col_i = 0
                for i in range(df_fwd_rev.shape[0]):
                    if col_i >= len(cols):
                        if row_i >= len(rows)-1:
                            row_i = 0
                            col_i = 0
                            plate_i += 1
                        else:
                            row_i += 1
                            col_i = 0
                    primers_ls.append(f"{fwd}_{rev}")
                    plate_ls.append(plate_i)
                    row_ls.append(rows[row_i])
                    col_ls.append(cols[col_i])
                    col_i += 1
                plate_i += 1
            
            df_plate = pd.concat(ls_df_fwd_rev,ignore_index=True) # Concatenate all PCR1 FWD/REV dataframes
            df_plate[f'{group} (PCR1)'] = [f"{plate}_{primers}" for (primers,plate) in zip(primers_ls,plate_ls)]
            df_plate['row'] = row_ls
            df_plate['column'] = col_ls

            return df_plate
        
        elif str(pcr)=='2':
            # Copy dataframe to avoid modifying original
            df2=df.copy()

            # Store gDNA and PCR locations on plate
            plate_ls = []
            row_ls = []
            col_ls = []

            plate_i = 1
            row_i = 0
            col_i = 0
            for i in range(df.shape[0]):
                if col_i >= len(cols):
                    if row_i >= len(rows)-1:
                        row_i = 0
                        col_i = 0
                        plate_i += 1
                    else:
                        row_i += 1
                        col_i = 0
                plate_ls.append(plate_i)
                row_ls.append(rows[row_i])
                col_ls.append(cols[col_i])
                col_i += 1

            df2[f'{group} (PCR2)'] = plate_ls
            df2['row'] = row_ls
            df2['column'] = col_ls

            return df2
        
        else:
            raise ValueError("pcr must be '1' or '2'")

    # Create 96-well and 8-strip plate layouts
    df_96_well_pcr1 = plate(df=df,group='96-well plate', pcr='1')
    df_96_well_pcr2 = plate(df=df,group='96-well plate', pcr='2')
    df_8_strip_pcr1 = plate(df=df,group='8-strip plate', pcr='1')
    df_8_strip_pcr2 = plate(df=df,group='8-strip plate', pcr='2')

    # Create pivot tables for gDNA, PCR1, and PCR2s
    pivots = {f"96-well_{gDNA_id_col}": pd.pivot_table(data=df_96_well_pcr1,values=gDNA_id_col,index=['96-well plate (PCR1)','row'],columns='column',aggfunc='first'),
            f"96-well_{pcr1_id_col}": pd.pivot_table(data=df_96_well_pcr1,values=pcr1_id_col,index=['96-well plate (PCR1)','row'],columns='column',aggfunc='first'),
            f"96-well_{pcr2_id_col}": pd.pivot_table(data=df_96_well_pcr2,values=pcr2_id_col,index=['96-well plate (PCR2)','row'],columns='column',aggfunc='first'),
            f"96-well_{pcr2_fwd_col}": pd.pivot_table(data=df_96_well_pcr2,values=pcr2_fwd_col,index=['96-well plate (PCR2)','row'],columns='column',aggfunc='first'),
            f"96-well_{pcr2_rev_col}": pd.pivot_table(data=df_96_well_pcr2,values=pcr2_rev_col,index=['96-well plate (PCR2)','row'],columns='column',aggfunc='first'),
            f"8-strip_{gDNA_id_col}": pd.pivot_table(data=df_8_strip_pcr1,values=gDNA_id_col,index=['8-strip plate (PCR1)','row'],columns='column',aggfunc='first'),
            f"8-strip_{pcr1_id_col}": pd.pivot_table(data=df_8_strip_pcr1,values=pcr1_id_col,index=['8-strip plate (PCR1)','row'],columns='column',aggfunc='first'),
            f"8-strip_{pcr2_id_col}": pd.pivot_table(data=df_8_strip_pcr2,values=pcr2_id_col,index=['8-strip plate (PCR2)','row'],columns='column',aggfunc='first'),
            f"8-strip_{pcr2_fwd_col}": pd.pivot_table(data=df_8_strip_pcr2,values=pcr2_fwd_col,index=['8-strip plate (PCR2)','row'],columns='column',aggfunc='first'),
            f"8-strip_{pcr2_rev_col}": pd.pivot_table(data=df_8_strip_pcr2,values=pcr2_rev_col,index=['8-strip plate (PCR2)','row'],columns='column',aggfunc='first')
            }

    # Get unique primer pairs and their value counts for PCR1 and PCR2
    df['PCR2 FWD MM'] = 'PCR2-FWD'
    pcr1_primers_vcs = t.vcs_ordered(df=df,cols=[pcr1_fwd_col,pcr1_rev_col])
    pcr2_primers_vcs = t.vcs_ordered(df=df,cols=['PCR2 FWD MM',pcr2_rev_col])

    # Split dataframe into UMI and non-UMI if UMI column exists
    if umi_col not in df.columns: # UMI column does not exist

        # Make PCR master mixes for PCR1 and PCR2
        if ultra:
            Q5_mm_x_stock = 2 # NEBNext Ultra II Q5 master mix stock
            pcr1_mms = pcr_mm_ultra(primers=pcr1_primers_vcs,template='gDNA Extract',template_uL=pcr1_total_uL-sum([Q5_mm_x_desired/Q5_mm_x_stock,fwd_uM_desired/fwd_uM_stock,rev_uM_desired/rev_uM_stock]*pcr1_total_uL),
                                    Q5_mm_x_stock=Q5_mm_x_stock,fwd_uM_stock=fwd_uM_stock,rev_uM_stock=rev_uM_stock,
                                    Q5_mm_x_desired=Q5_mm_x_desired,fwd_uM_desired=fwd_uM_desired,rev_uM_desired=rev_uM_desired,
                                    total_uL=pcr1_total_uL,mm_x=mm_x)
            pcr2_mms = pcr_mm_ultra(primers=pcr2_primers_vcs,template='PCR1 Product',template_uL=1,
                                    Q5_mm_x_stock=Q5_mm_x_stock,fwd_uM_stock=fwd_uM_stock,rev_uM_stock=rev_uM_stock,
                                    Q5_mm_x_desired=Q5_mm_x_desired,fwd_uM_desired=fwd_uM_desired,rev_uM_desired=rev_uM_desired,
                                    total_uL=pcr2_total_uL,mm_x=mm_x)
        else:
            pcr1_mms = pcr_mm(primers=pcr1_primers_vcs,template='gDNA Extract',template_uL=5,
                            Q5_mm_x_stock=Q5_mm_x_stock,dNTP_mM_stock=dNTP_mM_stock,fwd_uM_stock=fwd_uM_stock,rev_uM_stock=rev_uM_stock,
                            Q5_U_uL_stock=Q5_U_uL_stock,Q5_mm_x_desired=Q5_mm_x_desired,dNTP_mM_desired=dNTP_mM_desired,fwd_uM_desired=fwd_uM_desired,
                            rev_uM_desired=rev_uM_desired,Q5_U_uL_desired=Q5_U_uL_desired,total_uL=pcr1_total_uL,mm_x=mm_x)
            pcr2_mms = pcr_mm(primers=pcr2_primers_vcs,template='PCR1 Product',template_uL=1,
                            Q5_mm_x_stock=Q5_mm_x_stock,dNTP_mM_stock=dNTP_mM_stock,fwd_uM_stock=fwd_uM_stock,rev_uM_stock=rev_uM_stock,
                            Q5_U_uL_stock=Q5_U_uL_stock,Q5_mm_x_desired=Q5_mm_x_desired,dNTP_mM_desired=dNTP_mM_desired,fwd_uM_desired=fwd_uM_desired,
                            rev_uM_desired=rev_uM_desired,Q5_U_uL_desired=Q5_U_uL_desired,total_uL=pcr2_total_uL,mm_x=mm_x)

        # Create thermocycler objects for PCR1 and PCR2
        pcr1_thermo = thermocycler(df=df, n='1', cycles=pcr1_cycles, pcr_fwd_col=pcr1_fwd_col, pcr_rev_col=pcr1_rev_col)
        pcr2_thermo = thermocycler(df=df, n='2', cycles=pcr2_cycles, pcr_fwd_col=pcr2_fwd_col, pcr_rev_col=pcr2_rev_col)

        if dir is not None and file is not None: # Save file if dir & file are specified
            io.mkdir(dir=dir)
            with pd.ExcelWriter(os.path.join(dir,file)) as writer:
                sr = 0 # starting row
                for key,pivot in pivots.items():
                    pivot.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all pivots
                    pivot.to_excel(writer,sheet_name=key) # Pivot per sheet
                    sr += len(pivot)+2 # Skip 2 lines after each pivot
                for key,pcr1_mm in pcr1_mms.items():
                    pcr1_mm.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all PCR MMs
                    pcr1_mm.to_excel(writer,sheet_name='_'.join(key)) # PCR MM per sheet
                    sr += pcr1_mm.shape[0]+2 # Skip 2 lines after each PCR MM
                for key,pcr2_mm in pcr2_mms.items():
                    pcr2_mm.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all PCR MMs
                    pcr2_mm.to_excel(writer,sheet_name='_'.join(key)) # PCR MM per sheet
                    sr += pcr2_mm.shape[0]+2 # Skip 2 lines after each PCR MM
                for key,thermo in pcr1_thermo.items():
                    thermo.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all thermocyler objects
                    thermo.to_excel(writer,sheet_name=key) # Thermocyler object per sheet
                    sr += thermo.shape[0]+2 # Skip 2 lines after each thermocyler object
                for key,thermo in pcr2_thermo.items():
                    thermo.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all thermocyler objects
                    thermo.to_excel(writer,sheet_name=key) # Thermocyler object per sheet
                    sr += thermo.shape[0]+2 # Skip 2 lines after each thermocyler object
                
        return pivots,pcr1_mms,pcr2_mms,pcr1_thermo,pcr2_thermo

    else: # UMI column exists

        # Create 96-well and 8-strip plate layouts for UMI and non-UMI samples
        df_umi = df[df[umi_col]==True].copy().reset_index(drop=True)
        df_non_umi = df[df[umi_col]==False].copy().reset_index(drop=True)

        # Make PCR1.5 primers and get their value counts
        df_umi['PCR1.5 ID'] = df_umi[pcr1_id_col] + '_UMI'
        df_umi['PCR1.5 FWD'] = 'P5-FWD'
        df_umi['PCR1.5 REV'] = 'P7-REV'
        if pcr1_5_Tm is None:
            pcr1_5_Tm = 65
        df_umi['PCR1.5 Tm'] = pcr1_5_Tm
        
        pcr1_5_primers_vcs = t.vcs_ordered(df=df_umi,cols=['PCR1.5 FWD','PCR1.5 REV'])

        # Make PCR master mixes for PCR1 and PCR2
        if ultra:
            Q5_mm_x_stock = 2 # NEBNext Ultra II Q5 master mix stock
            pcr1_mms = pcr_mm_ultra(primers=pcr1_primers_vcs,template='gDNA Extract',template_uL=pcr1_total_uL-sum([Q5_mm_x_desired/Q5_mm_x_stock,fwd_uM_desired/fwd_uM_stock,rev_uM_desired/rev_uM_stock]*pcr1_total_uL),
                                    Q5_mm_x_stock=Q5_mm_x_stock,fwd_uM_stock=fwd_uM_stock,rev_uM_stock=rev_uM_stock,
                                    Q5_mm_x_desired=Q5_mm_x_desired,fwd_uM_desired=fwd_uM_desired,rev_uM_desired=rev_uM_desired,
                                    total_uL=pcr1_total_uL,mm_x=mm_x)
            pcr1_5_mms = pcr_mm_ultra(primers=pcr1_5_primers_vcs,template='PCR1 Product',template_uL=1,
                                    Q5_mm_x_stock=Q5_mm_x_stock,fwd_uM_stock=fwd_uM_stock,rev_uM_stock=rev_uM_stock,
                                    Q5_mm_x_desired=Q5_mm_x_desired,fwd_uM_desired=fwd_uM_desired,rev_uM_desired=rev_uM_desired,
                                    total_uL=pcr2_total_uL,mm_x=mm_x)
            pcr2_mms = pcr_mm_ultra(primers=pcr2_primers_vcs,template='PCR1.5 Product',template_uL=1,
                                    Q5_mm_x_stock=Q5_mm_x_stock,fwd_uM_stock=fwd_uM_stock,rev_uM_stock=rev_uM_stock,
                                    Q5_mm_x_desired=Q5_mm_x_desired,fwd_uM_desired=fwd_uM_desired,rev_uM_desired=rev_uM_desired,
                                    total_uL=pcr2_total_uL,mm_x=mm_x)
        else:
            pcr1_mms = pcr_mm(primers=pcr1_primers_vcs,template='gDNA Extract',template_uL=5,
                            Q5_mm_x_stock=Q5_mm_x_stock,dNTP_mM_stock=dNTP_mM_stock,fwd_uM_stock=fwd_uM_stock,rev_uM_stock=rev_uM_stock,
                            Q5_U_uL_stock=Q5_U_uL_stock,Q5_mm_x_desired=Q5_mm_x_desired,dNTP_mM_desired=dNTP_mM_desired,fwd_uM_desired=fwd_uM_desired,
                            rev_uM_desired=rev_uM_desired,Q5_U_uL_desired=Q5_U_uL_desired,total_uL=pcr1_total_uL,mm_x=mm_x)
            pcr1_5_mms = pcr_mm(primers=pcr1_5_primers_vcs,template='PCR1 Product',template_uL=1,
                            Q5_mm_x_stock=Q5_mm_x_stock,dNTP_mM_stock=dNTP_mM_stock,fwd_uM_stock=fwd_uM_stock,rev_uM_stock=rev_uM_stock,
                            Q5_U_uL_stock=Q5_U_uL_stock,Q5_mm_x_desired=Q5_mm_x_desired,dNTP_mM_desired=dNTP_mM_desired,fwd_uM_desired=fwd_uM_desired,
                            rev_uM_desired=rev_uM_desired,Q5_U_uL_desired=Q5_U_uL_desired,total_uL=pcr2_total_uL,mm_x=mm_x)
            pcr2_mms = pcr_mm(primers=pcr2_primers_vcs,template='PCR1.5 Product',template_uL=1,
                            Q5_mm_x_stock=Q5_mm_x_stock,dNTP_mM_stock=dNTP_mM_stock,fwd_uM_stock=fwd_uM_stock,rev_uM_stock=rev_uM_stock,
                            Q5_U_uL_stock=Q5_U_uL_stock,Q5_mm_x_desired=Q5_mm_x_desired,dNTP_mM_desired=dNTP_mM_desired,fwd_uM_desired=fwd_uM_desired,
                            rev_uM_desired=rev_uM_desired,Q5_U_uL_desired=Q5_U_uL_desired,total_uL=pcr2_total_uL,mm_x=mm_x)
        
        # Create thermocycler objects for PCR1 and PCR2
        # UMI specific
        umi_pcr1_thermo = thermocycler(df=df_umi, n='1', cycles=umi_cycles, pcr_fwd_col=pcr1_fwd_col, pcr_rev_col=pcr1_rev_col, umi=True)
        umi_pcr1_5_thermo = thermocycler(df=df_umi, n='1.5', cycles=pcr1_cycles, pcr_fwd_col='PCR1.5 FWD', pcr_rev_col='PCR1.5 REV', umi=True)

        # Non-UMI specific
        pcr1_thermo = thermocycler(df=df_non_umi, n='1', cycles=pcr1_cycles, pcr_fwd_col=pcr1_fwd_col, pcr_rev_col=pcr1_rev_col)

        # Common for UMI and non-UMI
        pcr2_thermo = thermocycler(df=df, n='2', cycles=pcr2_cycles, pcr_fwd_col=pcr2_fwd_col, pcr_rev_col=pcr2_rev_col)

        if dir is not None and file is not None: # Save file if dir & file are specified
            io.mkdir(dir=dir)
            with pd.ExcelWriter(os.path.join(dir,file)) as writer:
                sr = 0 # starting row
                for key,pivot in pivots.items():
                    pivot.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all pivots
                    pivot.to_excel(writer,sheet_name=key) # Pivot per sheet
                    sr += len(pivot)+2 # Skip 2 lines after each pivot
                for key,pcr1_mm in pcr1_mms.items():
                    pcr1_mm.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all PCR MMs
                    pcr1_mm.to_excel(writer,sheet_name='_'.join(key)) # PCR MM per sheet
                    sr += pcr1_mm.shape[0]+2 # Skip 2 lines after each PCR MM
                for key,pcr1_5_mm in pcr1_5_mms.items():
                    pcr1_5_mm.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all PCR MMs
                    pcr1_5_mm.to_excel(writer,sheet_name='_'.join(key)) # PCR MM per sheet
                    sr += pcr1_5_mm.shape[0]+2 # Skip 2 lines after each PCR MM
                for key,pcr2_mm in pcr2_mms.items():
                    pcr2_mm.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all PCR MMs
                    pcr2_mm.to_excel(writer,sheet_name='_'.join(key)) # PCR MM per sheet
                    sr += pcr2_mm.shape[0]+2 # Skip 2 lines after each PCR MM
                for key,thermo in umi_pcr1_thermo.items():
                    thermo.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all thermocyler objects
                    thermo.to_excel(writer,sheet_name=key) # Thermocyler object per sheet
                    sr += thermo.shape[0]+2 # Skip 2 lines after each thermocyler object
                for key,thermo in umi_pcr1_5_thermo.items():
                    thermo.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all thermocyler objects
                    thermo.to_excel(writer,sheet_name=key) # Thermocyler object per sheet
                    sr += thermo.shape[0]+2 # Skip 2 lines after each thermocyler object
                for key,thermo in pcr1_thermo.items():
                    thermo.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all thermocyler objects
                    thermo.to_excel(writer,sheet_name=key) # Thermocyler object per sheet
                    sr += thermo.shape[0]+2 # Skip 2 lines after each thermocyler object
                for key,thermo in pcr2_thermo.items():
                    thermo.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all thermocyler objects
                    thermo.to_excel(writer,sheet_name=key) # Thermocyler object per sheet
                    sr += thermo.shape[0]+2 # Skip 2 lines after each thermocyler object
                
        return pivots,pcr1_mms,pcr2_mms,umi_pcr1_thermo,umi_pcr1_5_thermo,pcr1_thermo,pcr2_thermo

def umis(genotypes: int, samples: int=1, cell_coverage: int=1000, ug_gDNA_per_cell: float=6*10**-6, 
         ploidy_per_cell: int=2, umi_coverage: float=5) -> tuple[float,int,int,int]:
    """
    umis(): Determine the ug of gDNA and reads required for genotyping with UMIs.

    Parameters:
    genotypes (int): The number of different genotypes in the sample.
    samples (int, optional): The number of samples to be processed (Default: 1).
    cell_coverage (int, optional): The number of cells per genotype (Default: 1000).
    ug_gDNA_per_cell (float, optional): The amount of genomic DNA per cell in micrograms (Default: 6x10^(-6) ug/cell).
    ploidy_per_cell (int, optional): The ploidy level of the cells (Default: 2 for diploid).
    umi_coverage (float): The desired UMI coverage (Default: 5).
    """
    # Calculate the ug of gDNA needed per sample
    ug = genotypes * cell_coverage * ug_gDNA_per_cell

    # Calculate the total number of molecules per sample
    molecules = genotypes * cell_coverage * ploidy_per_cell

    # Calculate reads needed per sample for desired UMI coverage
    reads_needed = molecules * umi_coverage
    reads_needed_samples = reads_needed * samples

    # Print and & return the results
    print(f"Per: {ug} ug gDNA; {molecules/10**6:.2f} M molecules; {reads_needed/10**6:.2f} M reads")
    print(f"All: {samples} samples; {reads_needed_samples/10**6:.2f} M reads")
    return ug,molecules,reads_needed,reads_needed_samples

# Hamming distance calculation
def hamming_distance(seq1: str | Seq, seq2: str | Seq) -> int:
    """
    hamming_distance(): returns the Hamming distance between two sequences

    Parameters:
    seq1 (str | Seq): sequence 1 
    seq1 (str | Seq): sequence 2
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length.")
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


def hamming_distance_matrix(df: pd.DataFrame | str, id: str, seqs: str, dir:str=None, file:str=None) -> pd.DataFrame:
    """
    hamming_distance_matrix(): compute pairwise Hamming distance matrix for a list of sequences stored in a dataframe

    Parameters:
    df (dataframe | str): pandas dataframe (or file path)
    id (str): id column name
    seqs (str): sequences column name
    dir (str, optional): save directory
    file (str, optional): save file

    Dependencies: pandas
    """
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)

    # Create empty hamming distance matrix
    n = len(df[seqs])
    matrix = [[0] * n for _ in range(n)]

    # Fill in hamming distance matrix
    for i in range(n):
        for j in range(n):
            if i <= j:
                dist = hamming_distance(df.iloc[i][seqs], df.iloc[j][seqs])
                matrix[i][j] = dist
                matrix[j][i] = dist  # since it's symmetric
    
    # Save & return hamming distance matrix
    df_matrix = pd.DataFrame(matrix,columns=df[id],index=df[seqs])
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=df_matrix.reset_index(drop=False))  
    return df_matrix