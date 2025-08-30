'''
Module: transfect.py
Author: Marc Zepeda
Created: 2024-11-07
Description: Transfection

Usage:
[Transfection calculation]
- PE3(): generates PE3 transfection plan for HEK293T cells (Default: 96-well plate in triplicate using L2000)
- virus(): generates transfection plan for virus production from HEK293T cells (Default: 6-well plate using L3000)
'''

# Import packages
import pandas as pd
import numpy as np
import os
import re
from ..gen import io
import warnings
warnings.filterwarnings("ignore")

# Transfection calculation
def PE3(plasmids: pd.DataFrame | str, epegRNAs: pd.DataFrame | str, ngRNAs: pd.DataFrame | str,
        dir:str=None, file:str=None, 
        pegRNA_number_col: str='pegRNA_number', epegRNAs_name_col: str='Name', ngRNAs_name_col: str='Name',
        plasmid_col: str='Plasmid', description_col: str='Description', colony_col: str='Colony', ng_uL_col: str='ng/uL',
        PE_plasmid: str='pMUZ86.7', reps: int=3, mm_x: float=1.1, epegRNA_ng: int=66, ngRNA_ng: int=22, PE_ng: int=200, well_uL:int=10) -> dict[pd.DataFrame]:
    '''
    PE3(): generates PE3 transfection plan for HEK293T cells (Default: 96-well plate in triplicate using L2000)

    Parameters:
    plasmids (DataFrame | str): all nanodrop concentrations (Default: google sheets format) 
    epegRNAs (DataFrame | str): epegRNAs from PrimeDesign
    ngRNAs (DataFrame | str): ngRNAs from PrimeDesign
    dir (optional): save directory
    file (optional): save file
    pegRNA_number_col (str, optional): pegRNA_number column name from PrimeDesign (Default: 'pegRNA_number')
    epegRNAs_name_col (str, optional): epegRNA column name (Default: 'Name')
    ngRNAs_name_col (str, optional): ngRNA column name (Default: 'Name')
    plasmid_col (str, optional): plasmid column name (Default: 'Plasmid')
    description_col (str, optional): plasmid column name (Default: 'Description')
    colony_col (str, optional): colony column name (Default: 'Colony')
    ng_uL_col (str, optional): ng/uL column name (Default: 'ng/uL')
    PE_plasmid (str, optional): search for PE plasmid name in plasmids (Default: 'pMUZ86.7')
    reps (int, optional): replicates (Default: 3)
    mm_x (float, optional): master mix multiplier (Default: 1.1)
    epegRNA_ng (int, optional): epegRNA ngs per well (Default: 66)
    ngRNA_ng (int, optional): ngRNA ngs per well (Default: 22)
    PE_ng (int, optional): PE ngs per well (Default: 200)
    well_uL (int, optional): uL transfection mix per well (Default: 10)

    Dependencies: pandas,numpy,os,io
    '''
    # Get dataframes from file path if needed
    if type(plasmids)==str:
        plasmids = io.get(pt=plasmids)
    if type(epegRNAs)==str:
        epegRNAs = io.get(pt=epegRNAs)
    if type(ngRNAs)==str:
        ngRNAs = io.get(pt=ngRNAs)
    
    # Determine transfection conditions (Tube As)
    tube_A = pd.DataFrame()
    for (pegRNA_number,epegRNA) in zip(epegRNAs[pegRNA_number_col],epegRNAs[epegRNAs_name_col]): # Retrieve epegRNA name & #
        for ngRNA in ngRNAs[ngRNAs[pegRNA_number_col]==pegRNA_number][ngRNAs_name_col]: # Retrieve corresponding ngRNAs
            # Get plasmid concentrations...
            epegRNA_df = plasmids[plasmids[plasmid_col]==epegRNA] # for epegRNAs
            ngRNA_df = plasmids[plasmids[plasmid_col]==ngRNA] # for ngRNAs
            
            # Change plasmid concentration column names...
            epegRNA_df = epegRNA_df.rename(columns={plasmid_col:'epegRNA', # for epegRNAs
                                                    colony_col:f'epegRNA {colony_col}',
                                                    description_col:f'epegRNA {description_col}',
                                                    ng_uL_col:f'epegRNA {ng_uL_col}'})
            ngRNA_df = ngRNA_df.rename(columns={plasmid_col:'ngRNA', # for ngRNAs
                                                colony_col:f'ngRNA {colony_col}',
                                                description_col:f'ngRNA {description_col}',
                                                'ng/uL':f'ngRNA {ng_uL_col}'})
            
            # Merge epegRNA & ngRNAs for each transfection condition (Tube As)
            epegRNA_df['ngRNA']=ngRNA
            tube_A = pd.concat([tube_A,pd.merge(left=epegRNA_df,
                                                right=ngRNA_df,
                                                on='ngRNA')])
    
    # Retrieve PE plasmid
    tube_A['PE'] = PE_plasmid
    tube_A[f'PE {colony_col}'] = plasmids[plasmids[plasmid_col]==PE_plasmid].iloc[0][colony_col]
    tube_A[f'PE {ng_uL_col}'] = plasmids[plasmids[plasmid_col]==PE_plasmid].iloc[0][ng_uL_col]

    # Calculate uLs for each reagent
    tube_A['epegRNA uL'] = [round(mm_x*reps*epegRNA_ng/ng_uL,2) for ng_uL in tube_A[f'epegRNA {ng_uL_col}']]
    tube_A['ngRNA uL'] = [round(mm_x*reps*ngRNA_ng/ng_uL,2) for ng_uL in tube_A[f'ngRNA {ng_uL_col}']]
    tube_A['PE uL'] = [round(mm_x*reps*PE_ng/ng_uL,2) for ng_uL in tube_A[f'PE {ng_uL_col}']]
    tube_A['Optimem uL'] = round(mm_x*reps*well_uL/2 - tube_A['epegRNA uL'] - tube_A['ngRNA uL'] - tube_A['PE uL'],2)
    tube_A['L2000 uL'] = 0
    tube_A['Total uL'] = round(mm_x*reps*well_uL/2,2)
    tube_A['Tube'] = [f'A{i+1}' for i in range(tube_A.shape[0])]

    # Append tube B calculations to transfection conditions (Tube As)
    tube_B = pd.DataFrame({'Optimem uL': [round(9/10*mm_x*reps*mm_x*tube_A.shape[0],2)],
                           'L2000 uL': [round(1/10*well_uL/2*mm_x*reps*mm_x*tube_A.shape[0],2)],
                           'Total uL': [round(well_uL/2*mm_x*reps*mm_x*tube_A.shape[0],2)],
                           'Tube': ['B']})
    transfection = pd.concat([tube_A,tube_B]).reset_index(drop=True)

    # Define 96-well plate axis
    rows_96_well = ['A','B','C','D','E','F','G','H']
    cols_96_well = np.arange(1,13,1)

    # Store transfection condition (Tube As) locations on 96-well plate (excluding outer wells)
    plate_ls = []
    row_ls = []
    col_ls = []

    plate_i = 1
    row_i = 1
    col_i = 1
    for i in range(tube_A.shape[0]):
        if col_i >= len(cols_96_well)-1:
            if row_i >= len(rows_96_well)-1-reps:
                row_i = 1
                col_i = 1
                plate_i += 1
            else:
                row_i += reps
                col_i = 1
        plate_ls.append(plate_i)
        row_ls.append(rows_96_well[row_i])
        col_ls.append(cols_96_well[col_i])
        col_i += 1

    tube_A['plate'] = plate_ls
    tube_A['row'] = row_ls
    tube_A['column'] = col_ls

    # Create pivot tables for tranfection & 96-well plates
    pivots = {'Transfection': transfection,
              '96-well plates': pd.pivot_table(data=tube_A,values='Tube',index=['plate','row'],columns='column',aggfunc='first'),
             }
    
    if dir is not None and file is not None: # Save file if dir & file are specified
        io.mkdir(dir=dir)
        with pd.ExcelWriter(os.path.join(dir,file)) as writer:
            sr = 0 # starting row
            for key,pivot in pivots.items():
                pivot.to_excel(writer,sheet_name='Tranfection Plan',startrow=sr) # Sheet with all pivots
                pivot.to_excel(writer,sheet_name=key) # Pivot per sheet
                sr += len(pivot)+2 # Skip 2 lines after each pivot
    return pivots

def virus(plasmids: pd.DataFrame | str, plasmid_col: str='Plasmid', description_col: str='Description', colony_col: str='Colony',
          ng_uL_col: str='ng/uL', VSVG_plasmid: str='pMUZ26.6',GagPol_plasmid: str='pMUZ26.7',
          reps: int=1, mm_x: float=1.1, VSVG_ng: int=750, GagPol_ng: int=1500, transfer_ng: int=750, well_uL: int=500,
          dir:str=None, file:str=None) -> pd.DataFrame:
    '''
    virus(): generates transfection plan for virus production from HEK293T cells (Default: 6-well plate using L3000)

    Parameters:
    plasmids (DataFrame): all nanodrop concentrations (Default: google sheets format) 
    lasmid_col (str, optional): plasmid column name (Default: 'Plasmid')
    description_col (str, optional): plasmid column name (Default: 'Description')
    colony_col (str, optional): colony column name (Default: 'Colony')
    ng_uL_col (str, optional): ng/uL column name (Default: 'ng/uL')
    VSVG_plasmid (str, optional): search for VSVG plasmid name in plasmids (Default: 'pMUZ26.6')
    GagPol_plasmid (str, optional): search for GagPol plasmid name in plasmids (Default: 'pMUZ26.7')
    reps (int, optional): replicates (Default: 1)
    mm_x (float, optional): master mix multiplier (Default: 1.1)
    VSVG_ng (int, optional): VSVG ngs per well (Default: 750)
    GagPol_ng (int, optional): ngRNA ngs per well (Default: 1500)
    transfer_ng (int, optional): PE ngs per well (Default: 750)
    well_uL (int, optional): uL transfection mix per well (Default: 500)
    dir (optional): save directory
    file (optional): save file

    Dependencies: pandas,re
    '''
    transfers_num = plasmids[(plasmids[plasmid_col]!=VSVG_plasmid)&(plasmids[plasmid_col]!=GagPol_plasmid)].shape[0]
    optimem_uL = round(transfers_num*mm_x*reps*well_uL/2,2)
    P3000_uL = round(transfers_num*mm_x*reps*6,2)
    L3000_uL = round(transfers_num*mm_x*reps*7,2)

    tube_ls = []
    order_ls = []
    ng_ls = []
    uL_ls = []

    for i,plasmid in enumerate(plasmids[plasmid_col]):
        if (plasmid==VSVG_plasmid):
            tube_ls.append('A')
            order_ls.append(1)
            ng_ls.append(round(transfers_num*mm_x*reps*VSVG_ng,2))
            uL_ls.append(round(transfers_num*mm_x*reps*VSVG_ng/plasmids.iloc[i][ng_uL_col],2))
        elif (plasmid==GagPol_plasmid):
            tube_ls.append('A')
            order_ls.append(1)
            ng_ls.append(round(transfers_num*mm_x*reps*GagPol_ng,2))
            uL_ls.append(round(transfers_num*mm_x*reps*GagPol_ng/plasmids.iloc[i][ng_uL_col],2))
        else:
            tube_ls.append(', '.join(re.findall(r'\d+\.\d+|\d+',plasmid)))
            order_ls.append(2)
            ng_ls.append(round(transfer_ng*reps,2))
            uL_ls.append(round(transfer_ng*reps/plasmids.iloc[i][ng_uL_col],2))

    plasmids['ng'] = ng_ls
    plasmids['uL'] = uL_ls
    plasmids['Tube'] = tube_ls
    plasmids['Order'] = order_ls
    
    # Save & return virus dataframe
    df_virus = pd.concat([pd.DataFrame({plasmid_col: ['N/A','N/A'],
                                        description_col: ['Optimem','P3000 Enhancer'],
                                        colony_col: ['N/A','N/A'],
                                        ng_uL_col: ['N/A','N/A'],
                                        'ng': ['N/A','N/A'],
                                        'uL': [optimem_uL,P3000_uL],
                                        'Tube': ['A','A'],
                                        'Order': [1,1]}),
                          plasmids.sort_values(by=['Order',plasmid_col]),
                          pd.DataFrame({plasmid_col: ['N/A','N/A'],
                                        description_col: ['Optimem','L3000'],
                                        colony_col: ['N/A','N/A'],
                                        ng_uL_col: ['N/A','N/A'],
                                        'ng': ['N/A','N/A'],
                                        'uL': [optimem_uL,L3000_uL],
                                        'Tube': ['B','B'],
                                        'Order': [3,3]})]).reset_index(drop=True)
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=df_virus,id=True) 
    return df_virus