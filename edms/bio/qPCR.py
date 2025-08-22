''' 
Module: qPCR.py
Author: Marc Zepeda
Created: 2024-09-06
Description: quantative Polymerase Chain Reaction

Usage:
[qPCR data retrieval and analysis]
- cfx_Cq(): retrieve RT-qPCR data from CFX Cq csv
- ddCq(): computes ΔΔCq mean and error for all samples holding target pairs constant
'''
# Import packages
import itertools
import pandas as pd
import numpy as np
from ..gen import io as io

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

def ddCq(data: pd.DataFrame | str, sample_col:str='Sample', target_col:str='Target', Cq_col:str='Cq',
         dir:str=None, file:str=None) -> pd.DataFrame:
    ''' 
    ddCq(): computes ΔΔCq mean and error for all samples holding target pairs constant
    
    Parameters:
    data (dataframe | str): Cq pandas dataframe (or file path)
    sample_col (str, optional): column name with cDNA sample identifier (Default: Sample)
    target_col (str, optional): column name with target identifier (Default: Target)
    Cq_col (str, optional): column name with Cq value (Default: Cq)
    dir (str, optional): save directory
    file (str, optional): save file

    Dependencies: pandas, numpy, & itertools
    '''
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
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=data4) 
    return data4