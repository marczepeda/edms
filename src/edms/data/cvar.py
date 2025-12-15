''' 
Module: cvar.py
Author: Marc Zepeda
Created: 2024-09-16
Description: ClinVar

Usage:
[ClinVar database]
- mutations(): returns ClinVar mutations dataframe for a given gene
- prevalence(): returns list of mutations sorted by prevalence on ClinVar

[Prime editing]
- priority_muts: returns the shared sequences library dataframe with priority mutations
- priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
'''

# Import packages
import pandas as pd
import re
import ast
from ..gen import io

# ClinVar database
def mutations(df: pd.DataFrame | str, gene_name:str, dir:str=None, file:str=None) -> pd.DataFrame:
    ''' 
    mutations: returns ClinVar mutations dataframe for a given gene.
    
    Parameters:
    df (dataframe | str): ClinVar dataframe (or file path)
    gene_name (str): gene name
    dir (str, optional): save directory
    file (str, optional): save file

    Dependencies: re, io
    '''
    # Isolate mutations corresponding to gene
    if type(df)==str: # Get ClinVar dataframe from file path if needed
        df = io.get(pt=df)
    df = df[df['Gene(s)']==gene_name].reset_index(drop=True)
    df = df.dropna(subset='Protein change').reset_index(drop=True)

    # Determine mutation positions, AAs before, and AAs after
    befores = []
    nums = []
    afters = []
    for aa_mut in df['Protein change']:
        num = int(re.findall(r'-?\d+',aa_mut)[0])
        nums.append(num)
        befores.append(aa_mut.split(str(num))[0])
        afters.append(aa_mut.split(str(num))[1])
    df['AA_position']=nums
    df['AA_before']=befores
    df['AA_after']=afters

    # Save & return dataframe
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=df) 
    return df

def prevalence(df: pd.DataFrame) -> list:
    ''' 
    prevalence(): returns list of mutations sorted by prevalence on ClinVar
    
    Parameters:
    df (dataframe): ClinVar mutations dataframe
    
    Dependencies: pandas
    '''
    return list(df['Protein change'].value_counts().keys())

# Prime editing
def priority_muts(pegRNAs_shared: pd.DataFrame, df_clinvar: pd.DataFrame | str,
                  dir:str=None, file:str=None) -> pd.DataFrame:
    ''' 
    priority_muts: returns the shared sequences library dataframe with priority mutations
    
    Parameters:
    pegRNAs_shared (dataframe | str): pegRNAs shared sequences library dataframe (or file path)
    df_clinvar (dataframe | str): ClinVar mutations() dataframe (or file path)
    dir (str, optional): save directory
    file (str, optional): save file
    
    Dependencies: pandas, ast, prevalence()
    '''
    # Get pegRNAs shared sequences library and ClinVar mutation() dataframes from file paths if needed
    if type(pegRNAs_shared)==str:
        pegRNAs_shared = io.get(pt=pegRNAs_shared)
    if type(df_clinvar)==str:
        df_clinvar = io.get(pt=df_clinvar)
    
    # Check if ClinVar mutations dataframe has 'Protein change' column without NaN values
    if df_clinvar["Protein change"].isna().any() == True:
        df_clinvar = mutations(df=df_clinvar, gene_name=pegRNAs_shared['Target_name'][0]) # Get ClinVar mutations dataframe for the gene in pegRNAs shared sequences library

    # Get list of priority mutants
    mut_priority_ls = prevalence(df=df_clinvar)

    # Determine priority mutations for pegRNAs shared sequences library
    priority_muts = []
    mutants_used = []
    for e,edits in enumerate(pegRNAs_shared['Edits']): # Search available edits for shared spacer & PBS sequence
        if type(edits)==str:
            for m,mutant in enumerate(mut_priority_ls): # Iterate through most clinically-relevant mutations
                if (mutant in set(ast.literal_eval(edits)))&(mutant not in mutants_used): # Select a clinically-relevant mutation that has not been used
                    priority_muts.append(mutant)
                    mutants_used.append(mutant)
                    break
            if len(priority_muts)!=e+1: # All clinically-relevant mutations have been used
                for edit in ast.literal_eval(edits): # Find edit that has not been used
                    if edit not in mutants_used:
                        priority_muts.append(edit)
                        mutants_used.append(edit)
                        break
        elif type(edits)==list:
            for m,mutant in enumerate(mut_priority_ls): # Iterate through most clinically-relevant mutations
                if (mutant in set(edits))&(mutant not in mutants_used): # Select a clinically-relevant mutation that has not been used
                    priority_muts.append(mutant)
                    mutants_used.append(mutant)
                    break
            if len(priority_muts)!=e+1: # All clinically-relevant mutations have been used
                for edit in edits: # Find edit that has not been used
                    if edit not in mutants_used:
                        priority_muts.append(edit)
                        mutants_used.append(edit)
                        break
        else: print('Error: pegRNAs_shared["Edits"] is not type string nor list')
    pegRNAs_shared['Priority_mut']=priority_muts

    # Save & return shared sequences library dataframe
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=pegRNAs_shared) 
    return pegRNAs_shared

def priority_edits(pegRNAs: pd.DataFrame | str, pegRNAs_shared: pd.DataFrame | str, df_clinvar: pd.DataFrame | str, 
                   dir:str=None, file:str=None) -> pd.DataFrame:
    ''' 
    priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
    
    Parameters:
    pegRNAs (dataframe | str): pegRNAs library dataframe (or file path)
    pegRNAs_shared (dataframe | str): pegRNAs shared sequences library dataframe (or file path)
    df_clinvar (dataframe | str): ClinVar mutations() dataframe (or file path)
    dir (str, optional): save directory
    file (str, optional): save file

    Dependencies: pandas, ast, & prevalence()
    '''
    # Get pegRNAs shared sequences library and ClinVar mutations() dataframes from file paths if needed
    if type(pegRNAs_shared)==str:
        pegRNAs_shared = io.get(pt=pegRNAs_shared)
    if type(df_clinvar)==str:
        df_clinvar = io.get(pt=df_clinvar)

    # Check if ClinVar mutations dataframe has 'Protein change' column without NaN values
    if df_clinvar["Protein change"].isna().any() == True:
        df_clinvar = mutations(df=df_clinvar, gene_name=pegRNAs_shared['Target_name'][0]) # Get ClinVar mutations dataframe for the gene in pegRNAs shared sequences library

    # Determine priority pegRNAs based on priority mutations from pegRNAs shared sequences library
    pegRNAs_priority = pd.DataFrame()
    for p,priority in enumerate(pegRNAs_shared['Priority_mut']):
        pegRNAs_temp = pegRNAs[pegRNAs['Edit']==priority] # Isolate priority mutations
        pegRNAs_temp = pegRNAs_temp[(pegRNAs_temp['Spacer_sequence']==pegRNAs_shared.iloc[p]['Spacer_sequence'])&(pegRNAs_temp['PBS_sequence']==pegRNAs_shared.iloc[p]['PBS_sequence'])] # Confirm spacer & PBS matches
        pegRNAs_temp.drop_duplicates(subset=['Spacer_sequence'],inplace=True) # Drop redundant pegRNAs (not sure if this is needed)
        pegRNAs_priority = pd.concat([pegRNAs_priority,pegRNAs_temp]).reset_index(drop=True)
        pegRNAs_priority['ClinVar_count'] = [df_clinvar['Protein change'].value_counts()[edit] if edit in df_clinvar['AA_mut'].to_list() else 0 for edit in pegRNAs_priority['Edit']]

    # Save & return pegRNAs priority dataframe
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=pegRNAs_priority) 
    return pegRNAs_priority