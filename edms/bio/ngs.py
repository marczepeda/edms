### ngs.py ###
# Author: Marc Zepeda
# Date: 2024-11-12

# Import Packages
import numpy as np
import pandas as pd
import os
from ..gen import io
import warnings
warnings.filterwarnings("ignore")

# NGS PCR calculation methods
def pcr_mm(primers: pd.Series, template: str, template_uL: int,
           Q5_mm_x_stock=5,dNTP_mM_stock=10,fwd_uM_stock=10,rev_uM_stock=10,Q5_U_uL_stock=2,
           Q5_mm_x_desired=1,dNTP_mM_desired=0.2,fwd_uM_desired=0.5,rev_uM_desired=0.5,Q5_U_uL_desired=0.02,
           total_uL=20,mm_x=1.1):
    '''
    pcr_mm:
    
    Parameters:
    primers (Series): value_counts() for primers
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
                                                     })
    return pcr_mm_dc
                                            

def pcrs(samples: pd.DataFrame, dir:str=None, file:str=None, gDNA_id_col='ID', 
         pcr1_id_col='PCR1 ID', pcr1_fwd_col='PCR1 FWD', pcr1_rev_col='PCR1 REV', 
         pcr2_id_col='PCR2 ID', pcr2_fwd_col='PCR2 FWD', pcr2_rev_col='PCR2 REV',
         Q5_mm_x_stock=5,dNTP_mM_stock=10,fwd_uM_stock=10,rev_uM_stock=10,Q5_U_uL_stock=2,
         Q5_mm_x_desired=1,dNTP_mM_desired=0.2,fwd_uM_desired=0.5,rev_uM_desired=0.5,Q5_U_uL_desired=0.02,
         total_uL=20,mm_x=1.1):
    '''
    pcrs(): generates NGS PCR plan automatically (Default: 96-well plates excludingn outer wells)
    
    Parameters:
    samples (DataFrame): NGS samples dataframe
    dir (optional): save directory
    file (optional): save file
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

    Dependencies: pandas,numpy,os,io
    '''
    # Define 96-well plate axis
    rows_96_well = ['A','B','C','D','E','F','G','H']
    cols_96_well = np.arange(1,13,1)

    # Store gDNA and PCR locations on 96-well plate (excluding outer wells)
    plate_ls = []
    row_ls = []
    col_ls = []

    plate_i = 1
    row_i = 1
    col_i = 1
    for i in range(samples.shape[0]):
        if col_i >= len(cols_96_well)-1:
            if row_i >= len(rows_96_well)-2:
                row_i = 1
                col_i = 1
                plate_i += 1
            else:
                row_i += 1
                col_i = 1
        plate_ls.append(plate_i)
        row_ls.append(rows_96_well[row_i])
        col_ls.append(cols_96_well[col_i])
        col_i += 1

    samples['plate'] = plate_ls
    samples['row'] = row_ls
    samples['column'] = col_ls

    # Create pivot tables for gDNA, PCR1, and PCR2s
    pivots = {gDNA_id_col: pd.pivot_table(data=samples,values=gDNA_id_col,index=['plate','row'],columns='column',aggfunc='first'),
              pcr1_id_col: pd.pivot_table(data=samples,values=pcr1_id_col,index=['plate','row'],columns='column',aggfunc='first'),
              pcr1_fwd_col: pd.pivot_table(data=samples,values=pcr1_fwd_col,index=['plate','row'],columns='column',aggfunc='first'),
              pcr1_rev_col: pd.pivot_table(data=samples,values=pcr1_rev_col,index=['plate','row'],columns='column',aggfunc='first'),
              pcr2_id_col: pd.pivot_table(data=samples,values=pcr2_id_col,index=['plate','row'],columns='column',aggfunc='first'),
              pcr2_fwd_col: pd.pivot_table(data=samples,values=pcr2_fwd_col,index=['plate','row'],columns='column',aggfunc='first'),
              pcr2_rev_col: pd.pivot_table(data=samples,values=pcr2_rev_col,index=['plate','row'],columns='column',aggfunc='first')
              }
    
    # Create PCR master mixes for PCR1 and PCR2 primer pairs
    samples['PCR2 FWD MM'] = 'PCR2 FWD'
    pcr1_mms = pcr_mm(primers=samples[[pcr1_fwd_col,pcr1_rev_col]].value_counts(),template='gDNA Extract',template_uL=5)
    pcr2_mms = pcr_mm(primers=samples[['PCR2 FWD MM',pcr2_rev_col]].value_counts(),template='PCR1 Product',template_uL=1)

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
                pcr1_mm.to_excel(writer,sheet_name='_'.join(key)) # Pivot per sheet
                sr += pcr1_mm.shape[0]+2 # Skip 2 lines after each pivot
            for key,pcr2_mm in pcr2_mms.items():
                pcr2_mm.to_excel(writer,sheet_name='NGS Plan',startrow=sr) # Sheet with all PCR MMs
                pcr2_mm.to_excel(writer,sheet_name='_'.join(key)) # Pivot per sheet
                sr += pcr2_mm.shape[0]+2 # Skip 2 lines after each pivot
            
    
    return pivots,pcr1_mms,pcr2_mms
        
        