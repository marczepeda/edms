''' 
Module: qPCR.py
Author: Marc Zepeda
Created: 2024-09-06
Description: quantative Polymerase Chain Reaction

Usage:
[qPCR data retrieval and analysis]
- cfx_Cq(): retrieve RT-qPCR data from CFX Cq csv
- ddCq(): computes ΔΔCq mean and error for all samples holding target pairs constant
- find_cfx(): find a CFX export csv in a directory by its report name
- cfx_amp(): retrieve qPCR amplification curves from CFX Quantification Amplification Results csv as a tidy dataframe

[qPCR visualization]
- amp(): plot qPCR amplification (sigmoidal) curves (RFU vs. Cycle) with one line per well
'''
# Import packages
import itertools
import glob
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ..gen import io as io
from ..gen import plot as p

# qPCR data retrieval and analysis
def cfx_Cq(pt: str, sample_col:str='Sample', cols: list=['Well','Fluor','Target','Sample','Cq']) -> pd.DataFrame:
    ''' 
    cfx_Cq(): retrieve RT-qPCR data from CFX Cq csv
    
    Parameters:
    pt (str): path to Cq csv file
    sample_col (str, optional): column name with cDNA sample identifier (Default: Sample)
    cols (list): list of column names to retain (Default: ['Well','Fluor','Target','Sample','Cq'])
    
    Dependencies: io
    '''
    data = io.get(pt).dropna(subset=[sample_col])[cols]
    data[sample_col] = [int(cDNA) if type(cDNA)==float else cDNA for cDNA in data[sample_col]]
    return data

def ddCq(data: pd.DataFrame | str = None, sample_col:str='Sample', target_col:str='Target', Cq_col:str='Cq',
         file:str=None) -> pd.DataFrame:
    ''' 
    ddCq(): computes ΔΔCq mean and error for all samples holding target pairs constant
    
    Parameters:
    data (dataframe | str, optional): Cq pandas dataframe (or file path) (Default: None; find '*Quantification Cq Results*.csv' in the current directory)
    sample_col (str, optional): column name with cDNA sample identifier (Default: Sample)
    target_col (str, optional): column name with target identifier (Default: Target)
    Cq_col (str, optional): column name with Cq value (Default: Cq)
    file (str, optional): output file path

    Dependencies: pandas, numpy, itertools, & find_cfx()
    '''
    # Find CFX Cq export if needed
    if data is None:
        data = find_cfx(report='Quantification Cq Results')

    # Get dataframe from file path if needed
    if type(data)==str:
        data = cfx_Cq(pt=data,sample_col=sample_col,cols=[sample_col,target_col,Cq_col])

    # Get sample and target lists
    sample_ls = list(data[sample_col].value_counts().keys())
    target_ls = list(data[target_col].value_counts().keys())

    # Compute Cq mean and error for each set of samples & targets
    samples = []
    targets = []
    Cq_means = []
    Cq_errs = []
    for sample in sample_ls: # Isolate samples
        temp = data[data[sample_col]==sample] 
        for target in target_ls: # Isolate targets & compute
            samples.append(sample)
            targets.append(target)
            Cq_means.append(np.mean(temp[temp[target_col]==target][Cq_col].to_list()))
            Cq_errs.append(np.std(temp[temp[target_col]==target][Cq_col].to_list()))
    data2 = pd.DataFrame({'Sample':samples,'Target':targets,'Cq_mean':Cq_means,'Cq_err':Cq_errs})
    
    # Compute ΔCq mean and error for all target pairs within each set of samples
    samples = []
    target_pairs = []
    dCq_means = []
    dCq_errs = []
    for sample in sample_ls: # Isolate samples
        temp = data2[data2['Sample']==sample]
        for target_pair in list(itertools.combinations(target_ls,2)): # Isolate target pairs
            samples.append(sample)
            target_pairs.append(f'{target_pair[0]} ~ {target_pair[1]}')
            dCq_means.append(temp.iloc[0]['Cq_mean'] - temp.iloc[1]['Cq_mean'])
            dCq_errs.append(np.sqrt(temp.iloc[0]['Cq_err']**2+temp.iloc[1]['Cq_err']**2))
    data3 = pd.DataFrame({'Sample':samples,'Targets':target_pairs,'dCq_mean':dCq_means,'dCq_err':dCq_errs})
    
    # Compute ΔΔCq mean and error for all sample holding target pairs constant
    samples_pairs = []
    sample1s = []
    sample2s = []
    target_pairs = []
    target1s = []
    target2s = []
    ddCq_means = []
    ddCq_errs = []
    RQ_means = []
    RQ_errs = []
    for target_pair in list(data3['Targets'].value_counts().keys()): # Isolate target pairs
        temp=data3[data3['Targets']==target_pair]
        for i in range(len(temp)): # Isolate 1 sample to compare to the rest of the samples
            temp2 = temp.drop(i)
            for j in range(len(temp2)): # Iterate through the rest of the samples
                samples_pairs.append(f'{temp2.iloc[j]["Sample"]} ~ {temp.iloc[i]["Sample"]}')
                sample1s.append(temp2.iloc[j]["Sample"])
                sample2s.append(temp.iloc[i]["Sample"])
                target_pairs.append(target_pair)
                target1s.append(target_pair.split(' ~ ')[0])
                target2s.append(target_pair.split(' ~ ')[1])
                ddCq_means.append(temp2.iloc[j]['dCq_mean'] - temp.iloc[i]['dCq_mean'])
                ddCq_errs.append(np.sqrt(temp2.iloc[j]['dCq_err']**2+temp.iloc[i]['dCq_err']**2))
                RQ_means.append(2**(-ddCq_means[-1]))
                RQ_errs.append(np.abs(ddCq_means[-1]*np.log(2)*ddCq_errs[-1]))
    data4 = pd.DataFrame({'Samples':samples_pairs,'Sample 1':sample1s,'Sample 2': sample2s,'Targets':target_pairs,'Target 1':target1s,'Target 2':target2s,'ddCq_mean':ddCq_means,'ddCq_err':ddCq_errs,'RQ_mean':RQ_means,'RQ_err':RQ_errs})

    # Save & return analyzed qPCR data
    io.save(obj=data4, file=file) 
    return data4

def _norm_well(well: str) -> str:
    '''
    _norm_well(): normalize well names to the unpadded format (e.g., A01 -> A1)
    
    Parameters:
    well (str): well name
    '''
    match = re.fullmatch(r'([A-Za-z]+)0*(\d+)', str(well).strip())
    return f'{match.group(1).upper()}{int(match.group(2))}' if match else well

def find_cfx(report: str, dir: str = '.', required: bool = True) -> str | None:
    '''
    find_cfx(): find a CFX export csv in a directory by its report name
    
    Parameters:
    report (str): CFX report name in the file name (e.g., 'Quantification Amplification Results', 'Quantification Cq Results')
    dir (str, optional): directory to search (Default: current directory)
    required (bool, optional): raise an error if no file is found (Default: True); otherwise return None
    
    Dependencies: glob, os
    '''
    pts = sorted(glob.glob(os.path.join(glob.escape(dir), f'*{report}*.csv')))
    if len(pts) == 1:
        print(f"Found {report}: {pts[0]}")
        return pts[0]
    elif len(pts) > 1:
        raise ValueError(f"Multiple '{report}' csv files in {os.path.abspath(dir)}; specify one:\n" + '\n'.join(pts))
    elif required:
        raise FileNotFoundError(f"No '*{report}*.csv' file in {os.path.abspath(dir)}")
    return None

def cfx_amp(pt: str, annot: pd.DataFrame | str = None, wells: list = None, wells_exclude: list = None,
            baseline: tuple = None, drop_no_Cq: bool = False) -> pd.DataFrame:
    '''
    cfx_amp(): retrieve qPCR amplification curves from CFX Quantification Amplification Results csv as a tidy dataframe
    
    Parameters:
    pt (str): path to Quantification Amplification Results csv file (columns: Cycle, A1, A2, ...)
    annot (dataframe | str, optional): well annotations (or file path) with a Well column, such as the CFX
                                       Quantification Cq Results csv or a custom plate layout; other columns
                                       (e.g., Target, Sample, Content, Cq) are merged onto each well
    wells (list, optional): wells to keep (Default: None; keep all wells)
    wells_exclude (list, optional): wells to exclude (Default: None)
    baseline (tuple, optional): (start, end) cycles (inclusive) whose mean RFU is subtracted from each well
                                (Default: None; CFX exports are already baseline subtracted)
    drop_no_Cq (bool, optional): drop wells without a Cq value; requires annot with a Cq column (Default: False)
    
    Returns:
    pd.DataFrame: tidy dataframe with Well, Row, Column, Cycle, RFU, and annotation columns
    
    Dependencies: re, pandas, io, _norm_well()
    '''
    # Wide (Cycle x Well) -> tidy (Well, Cycle, RFU)
    data = io.get(pt)
    data = data.loc[:, ~data.columns.astype(str).str.startswith('Unnamed')]
    data = data.melt(id_vars='Cycle', var_name='Well', value_name='RFU')
    data['Well'] = data['Well'].map(_norm_well)

    # Filter wells
    if wells is not None:
        data = data[data['Well'].isin([_norm_well(w) for w in wells])]
    if wells_exclude is not None:
        data = data[~data['Well'].isin([_norm_well(w) for w in wells_exclude])]
    
    # Plate row & column (e.g., to facet into a plate view)
    data['Row'] = data['Well'].str.extract(r'^([A-Za-z]+)', expand=False)
    data['Column'] = data['Well'].str.extract(r'(\d+)$', expand=False).astype(int)

    # Baseline subtraction
    if baseline is not None:
        start, end = baseline
        base = data[(data['Cycle'] >= start) & (data['Cycle'] <= end)].groupby('Well')['RFU'].mean()
        data['RFU'] = data['RFU'] - data['Well'].map(base)

    # Merge well annotations (dropping index & empty columns)
    if annot is not None:
        if type(annot) == str:
            annot = io.get(annot)
        annot = annot.loc[:, ~annot.columns.astype(str).str.startswith('Unnamed')].dropna(axis=1, how='all').copy()
        if 'Well' not in annot.columns:
            raise ValueError(f"annot must have a 'Well' column; columns: {list(annot.columns)}")
        annot['Well'] = annot['Well'].map(_norm_well)
        annot = annot.drop(columns=[col for col in annot.columns if col in data.columns and col != 'Well'])
        data = data.merge(annot, on='Well', how='left')
    
    if drop_no_Cq:
        if 'Cq' not in data.columns:
            raise ValueError("drop_no_Cq requires annot with a 'Cq' column (e.g., Quantification Cq Results csv)")
        data = data.dropna(subset=['Cq'])

    return data.reset_index(drop=True)

# qPCR visualization
def amp(df: pd.DataFrame | str = None, annot: pd.DataFrame | str = None, wells: list = None, wells_exclude: list = None,
        baseline: tuple = None, drop_no_Cq: bool = False, threshold: float = None, threshold_color: str = 'black',
        cols: str = None, x_axis: str = 'Cycle', y_axis: str = 'RFU', y_axis_scale: str = 'linear',
        file: str = None, dpi: int = 0, transparent: bool = True, show: bool = True, **kwargs):
    '''
    amp(): plot qPCR amplification (sigmoidal) curves (RFU vs. Cycle) with one line per well
    
    Parameters:
    df (dataframe | str, optional): tidy amplification dataframe from cfx_amp() or path to CFX Quantification Amplification Results csv
                                    (Default: None; find '*Quantification Amplification Results*.csv' in the current directory)
    annot (dataframe | str, optional): well annotations (or file path) with a Well column, such as the CFX
                                       Quantification Cq Results csv or a custom plate layout
                                       (Default: None; use '*Quantification Cq Results*.csv' next to df if present)
    wells (list, optional): wells to keep (Default: None; keep all wells)
    wells_exclude (list, optional): wells to exclude (Default: None)
    baseline (tuple, optional): (start, end) cycles (inclusive) whose mean RFU is subtracted from each well (Default: None)
    drop_no_Cq (bool, optional): drop wells without a Cq value; requires annot with a Cq column (Default: False)
    threshold (float, optional): draw a horizontal threshold line at this RFU (Default: None)
    threshold_color (str, optional): threshold line color (Default: black)
    cols (str, optional): color column name (e.g., Well, Sample, Target, Row, Column; Default: None)
    x_axis (str, optional): x-axis name (Default: Cycle)
    y_axis (str, optional): y-axis name (Default: RFU)
    y_axis_scale (str, optional): y-axis scale; 'log' drops RFU <= 0 (Default: linear)
    file (str, optional): output plot file path
    dpi (int, optional): figure dpi (Default: 1200 for non-HTML, 150 for HTML)
    transparent (bool, optional): save static images with transparent background (Default: True)
    show (bool, optional): show plot (Default: True)
    **kwargs: plot.scat() parameters (e.g., facetx='Column', facety='Row' for a plate view)

    Dependencies: os, pandas, matplotlib, plot, find_cfx(), cfx_amp()
    '''
    # Find CFX exports if needed
    if df is None:
        df = find_cfx(report='Quantification Amplification Results')
    if type(df) == str and annot is None:
        annot = find_cfx(report='Quantification Cq Results', dir=os.path.dirname(df) or '.', required=False)

    # Get tidy dataframe from file path if needed
    if type(df) == str:
        df = cfx_amp(pt=df, annot=annot, wells=wells, wells_exclude=wells_exclude,
                     baseline=baseline, drop_no_Cq=drop_no_Cq)
    
    # Log scale can't show RFU <= 0
    if y_axis_scale == 'log':
        df = df[df['RFU'] > 0]
    
    # One line per well (units) without aggregating wells that share a color (estimator=None)
    plot_kwargs = dict(graph='line', df=df, x='Cycle', y='RFU', cols=cols, units='Well', estimator=None,
                       x_axis=x_axis, y_axis=y_axis, y_axis_scale=y_axis_scale, dpi=dpi, transparent=transparent, **kwargs)
    if threshold is None:
        return p.scat(file=file, show=show, **plot_kwargs)
    
    # Draw threshold on every panel before saving/showing
    fig, axes = p.scat(file=None, show=False, **plot_kwargs)
    for ax in axes.flat:
        ax.axhline(threshold, color=threshold_color, linestyle='--', linewidth=1, zorder=0)
    plt.figure(fig) # Re-register figure closed by scat()
    p._final_save_show(fig, file=file, dpi=dpi, transparent=transparent, icon='scatter', show=show)
    return fig, axes
