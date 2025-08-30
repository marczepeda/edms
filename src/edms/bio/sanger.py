'''
Module: sanger.py
Author:LJ Knight
Created: 2025-07-24
Description:Sanger sequencing

Usage:
[Sanger Thermocycler]
- group_boundaries(): Group a list of integers into segments where each segment contains consecutive numbers.
- min_sec(): Convert decimal minutes to a tuple of (minutes, seconds).
- thermocycler(): Creates a dictionary of thermocycler objects from a DataFrame.

[Sanger PCR calculation]
- pcr_mm(): NEB Q5 PCR master mix calculations
- pcr_mm_ultra(): NEBNext Ultra II Q5 PCR master mix calculations
- pcrs(): generates Sanger PCR plan automatically
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
import warnings
warnings.filterwarnings("ignore")

# Sanger Thermocycler
def group_boundaries(nums: list[int]) -> list[tuple[int,int]]:
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

def min_sec(decimal_minutes: float) -> tuple[int,int]:
    '''
    min_sec(): Convert decimal minutes to a tuple of (minutes, seconds).

    Parameters:
    decimal_minutes (float): Decimal representation of minutes.
    '''
    minutes = int(decimal_minutes)
    seconds = int(round((decimal_minutes - minutes) * 60))
    return minutes, seconds

def thermocycler(df: pd.DataFrame, n: Literal[1] = 1, cycles: int | str = None, pcr_fwd_col: str=None, pcr_rev_col: str=None) -> dict[pd.DataFrame]:
    """
    thermocycler(): Creates a dictionary of thermocycler objects from a DataFrame.

    Parameters:
    df (pd.DataFrame): DataFrame containing PCR data.
    n (Literal[1): The PCR number to process (Default: 1).
    cycles (int | str): Number of cycles for PCR1 (Default: None -> 30).
    pcr_fwd_col (str, optional): PCR FWD column name (Default: None -> f'PCR{n} FWD').
    pcr_rev_col (str, optional): PCR REV column name (Default: None -> f'PCR{n} REV').
    
    Dependencies: math, pandas, typing.Literal, tidy, min_sec(), group_boundaries()
    """
    dc = dict() # Initialize an empty dictionary to hold thermocycler dataframes
     # Set cycles and primer columns based on the PCR number
    if n==1: # PCR1
        if cycles is None:
            cycles = 30
        if pcr_fwd_col is None:
            pcr_fwd_col = 'PCR1 FWD'
        if pcr_rev_col is None:
            pcr_rev_col = 'PCR1 REV'
    else: 
        raise ValueError("n must be 1")
    
    # Iterate through unique combinations of FWD and REV primers for the specified PCR number
    for (fwd,rev) in t.unique_tuples(df=df, cols=[pcr_fwd_col,pcr_rev_col]):

        title = f"{fwd}_{rev}: " # Initialize the thermocycler title for the set of primers

        # Group the IDs based on the set of primers; determine the temperature and anneal_time
        ls = group_boundaries(list(df[(df[pcr_fwd_col]==fwd) & (df[pcr_rev_col]==rev)]['ID'].keys()))
        tm = df[(df[pcr_fwd_col]==fwd) & (df[pcr_rev_col]==rev)].iloc[0][f'PCR{n} Tm']
        bp = df[(df[pcr_fwd_col]==fwd) & (df[pcr_rev_col]==rev)].iloc[0]['PCR1 bp']
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

# Sanger calculation
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
         Q5_mm_x_stock: int=5, dNTP_mM_stock: int=10, fwd_uM_stock: int=10, rev_uM_stock: int=10, Q5_U_uL_stock: int=2,
         Q5_mm_x_desired: int=1,dNTP_mM_desired: float=0.2, fwd_uM_desired: float=0.5, rev_uM_desired: float=0.5, Q5_U_uL_desired: float=0.02,
         total_uL:int=20, mm_x: float=1.1, cycles: int | str = None, ultra: bool=False) -> tuple[dict[pd.DataFrame]]:
    '''
    pcrs(): generates Sanger PCR plan automatically
    
    Parameters:
    df (DataFrame | str): Sanger samples dataframe (or file path)
    dir (str, optional): save directory
    file (str, optional): save file
    gDNA_id_col (str, optional): gDNA ID column name (Default: 'ID')
    pcr1_id_col (str, optional): PCR1 ID column name (Default: 'PCR1 ID')
    pcr1_fwd_col (str, optional): PCR1 FWD column name (Default: 'PCR1 FWD')
    pcr1_rev_col (str, optional): PCR1 REV column name (Default: 'PCR1 REV')
    pcr1_id_col (str, optional): PCR2 ID column name (Default: 'PCR2 ID')
    pcr1_fwd_col (str, optional): PCR2 FWD column name (Default: 'PCR2 FWD')
    pcr1_rev_col (str, optional): PCR2 REV column name (Default: 'PCR2 REV')
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
    cycles (int | str): Number of cycles for the PCR1 process (Default: None -> 30).
    ultra (bool, optional): use NEB Ultra II reagents (Default: False)

    Dependencies: pandas,numpy,os,io,tidy,thermocycler(),pcr_mm(),pcr_mm_ultra()
    '''
    # Get samples dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)

    # Define plate axis method
    def plate(df: pd.DataFrame, group: Literal['96-well plate','8-strip_plate']) -> pd.DataFrame:
        ''' 
        plate(): Creates a DataFrame representing a 96-well plate or 8-strip plate layout.

        Parameters:
        df (pd.DataFrame): DataFrame containing PCR data.
        group (Literal['96-well plate','8-strip_plate']): Type of plate layout to create.
        '''
        if group == '96-well plate': # Define 96-well plate axis
            rows = ['A','B','C','D','E','F','G','H']
            cols = np.arange(1,13,1)
        elif group == '8-strip plate': # Define PCR strip axis
            rows = ['A','B','C','D','E','F']
            cols = np.arange(1,9,1)
        else:
            raise ValueError("group must be '96-well plate' or '8-strip plate'")

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
    
    # Create 96-well and 8-strip plate layouts
    df_96_well_pcr1 = plate(df=df,group='96-well plate')
    df_8_strip_pcr1 = plate(df=df,group='8-strip plate')

    # Create pivot tables for gDNA and PCR1
    pivots = {f"96-well_{gDNA_id_col}": pd.pivot_table(data=df_96_well_pcr1,values=gDNA_id_col,index=['96-well plate (PCR1)','row'],columns='column',aggfunc='first'),
            f"96-well_{pcr1_id_col}": pd.pivot_table(data=df_96_well_pcr1,values=pcr1_id_col,index=['96-well plate (PCR1)','row'],columns='column',aggfunc='first'),
            f"8-strip_{gDNA_id_col}": pd.pivot_table(data=df_8_strip_pcr1,values=gDNA_id_col,index=['8-strip plate (PCR1)','row'],columns='column',aggfunc='first'),
            f"8-strip_{pcr1_id_col}": pd.pivot_table(data=df_8_strip_pcr1,values=pcr1_id_col,index=['8-strip plate (PCR1)','row'],columns='column',aggfunc='first'),
            }
    
    # Get unique primer pairs and their value counts for PCR1
    pcr1_primers_vcs = t.vcs_ordered(df=df,cols=[pcr1_fwd_col,pcr1_rev_col])
    if ultra:
        Q5_mm_x_stock = 2 # NEBNext Ultra II Q5 master mix stock
        pcr1_mms = pcr_mm_ultra(primers=pcr1_primers_vcs,template='gDNA Extract',template_uL=total_uL-sum([Q5_mm_x_desired/Q5_mm_x_stock,fwd_uM_desired/fwd_uM_stock,rev_uM_desired/rev_uM_stock]*total_uL),
                                Q5_mm_x_stock=Q5_mm_x_stock,fwd_uM_stock=fwd_uM_stock,rev_uM_stock=rev_uM_stock,
                                Q5_mm_x_desired=Q5_mm_x_desired,fwd_uM_desired=fwd_uM_desired,rev_uM_desired=rev_uM_desired,
                                total_uL=total_uL,mm_x=mm_x)
        
    else:
        pcr1_mms = pcr_mm(primers=pcr1_primers_vcs,template='gDNA Extract',template_uL=5,
                          Q5_mm_x_stock=Q5_mm_x_stock,dNTP_mM_stock=dNTP_mM_stock,fwd_uM_stock=fwd_uM_stock,rev_uM_stock=rev_uM_stock,
                          Q5_U_uL_stock=Q5_U_uL_stock,Q5_mm_x_desired=Q5_mm_x_desired,dNTP_mM_desired=dNTP_mM_desired,fwd_uM_desired=fwd_uM_desired,
                          rev_uM_desired=rev_uM_desired,Q5_U_uL_desired=Q5_U_uL_desired,total_uL=total_uL,mm_x=mm_x)
        
    # Create thermocycler objects for PCR1 and PCR2
    pcr1_thermo = thermocycler(df=df, n=1, cycles=cycles, pcr_fwd_col=pcr1_fwd_col, pcr_rev_col=pcr1_rev_col)
        
    if dir is not None and file is not None: # Save file if dir & file are specified
        io.mkdir(dir=dir)
        with pd.ExcelWriter(os.path.join(dir,file)) as writer:
            sr = 0 # starting row
            for key,pivot in pivots.items():
                pivot.to_excel(writer,sheet_name='Sanger Plan',startrow=sr) # Sheet with all pivots
                pivot.to_excel(writer,sheet_name=key) # Pivot per sheet
                sr += len(pivot)+2 # Skip 2 lines after each pivot
            for key,pcr1_mm in pcr1_mms.items():
                pcr1_mm.to_excel(writer,sheet_name='Sanger Plan',startrow=sr) # Sheet with all PCR MMs
                pcr1_mm.to_excel(writer,sheet_name='_'.join(key)) # PCR MM per sheet
                sr += pcr1_mm.shape[0]+2 # Skip 2 lines after each PCR MM
            for key,thermo in pcr1_thermo.items():
                thermo.to_excel(writer,sheet_name='Sanger Plan',startrow=sr) # Sheet with all thermocyler objects
                thermo.to_excel(writer,sheet_name=key) # Thermocyler object per sheet
                sr += thermo.shape[0]+2 # Skip 2 lines after each thermocyler object
            
    return pivots,pcr1_mms,pcr1_thermo

