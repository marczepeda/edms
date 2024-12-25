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
def mutations(gene_name:str,pt:str,typ:str='tsv'):
    ''' 
    mutations: returns ClinVar mutations dataframe for a given gene.
    
    Parameters:
    gene_name: gene name
    pt: path to ClinVar tsv file

    Dependencies: re, io
    '''
    # Isolate mutations corresponding to gene
    gene = io.get(pt,typ=typ)
    gene = gene[gene['Gene(s)']==gene_name].reset_index(drop=True)
    gene = gene.dropna(subset='Protein change').reset_index(drop=True)

    # Determine mutation positions, AAs before, and AAs after
    befores = []
    nums = []
    afters = []
    aa_muts = []
    for aa_mut in gene['Protein change']:
        num = int(re.findall(r'-?\d+',aa_mut)[0])
        nums.append(num)
        befores.append(aa_mut.split(str(num))[0])
        afters.append(aa_mut.split(str(num))[1])
    gene['AA_position']=nums
    gene['AA_before']=befores
    gene['AA_after']=afters

    return gene

def prevalence(gene: pd.DataFrame):
    ''' 
    prevalence(): returns list of mutations sorted by prevalence on ClinVar
    
    Parameters:
    gene (dataframe): ClinVar mutations dataframe
    
    Dependencies: pandas
'''
    return list(gene['Protein change'].value_counts().keys())

# Prime editing
def priority_muts(pegRNAs: pd.DataFrame, pegRNAs_shared: pd.DataFrame, pt: str):
    ''' 
    priority_muts: returns the shared sequences library dataframe with priority mutations
    
    Parameters:
    pegRNAs (dataframe): pegRNAs library dataframe
    pegRNAs_shared (dataframe): pegRNAs shared sequences library dataframe
    pt (str): path to ClinVar csv file
    
    Dependencies: pandas, ast, prevalence(), & mutations()
    '''
    # Get list of priority mutants from prevalence() & mutations()
    mut_priority_ls = prevalence(gene=mutations(pt=pt))

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

    return pegRNAs_shared

def priority_edits(pegRNAs: pd.DataFrame, pegRNAs_shared: pd.DataFrame, pt: str):
    ''' 
    priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
    
    Parameters:
    pegRNAs (dataframe): pegRNAs library dataframe
    pegRNAs_shared (dataframe): pegRNAs shared sequences library dataframe
    pt (dataframe): path to ClinVar csv file
    gene (dataframe): ClinVar mutations dataframe for the gene from mutations()
    
    Dependencies: pandas, ast, & prevalence()
    '''
    # Get ClinVar mutations dataframe from mutations()
    gene=mutations(pt=pt)

    # Determine priority pegRNAs based on priority mutations from pegRNAs shared sequences library
    pegRNAs_priority = pd.DataFrame()
    for p,priority in enumerate(pegRNAs_shared['Priority_mut']):
        pegRNAs_temp = pegRNAs[pegRNAs['Edit']==priority] # Isolate priority mutations
        pegRNAs_temp = pegRNAs_temp[(pegRNAs_temp['Spacer_sequence']==pegRNAs_shared.iloc[p]['Spacer_sequence'])&(pegRNAs_temp['PBS_sequence']==pegRNAs_shared.iloc[p]['PBS_sequence'])] # Confirm spacer & PBS matches
        pegRNAs_temp.drop_duplicates(subset=['Spacer_sequence'],inplace=True) # Drop redundant pegRNAs (not sure if this is needed)
        pegRNAs_priority = pd.concat([pegRNAs_priority,pegRNAs_temp]).reset_index(drop=True)
        pegRNAs_priority['ClinVar_count'] = [gene['Protein change'].value_counts()[edit] for edit in pegRNAs_priority['Edit']]

    return pegRNAs_priority