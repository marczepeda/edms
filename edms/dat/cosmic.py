''' 
Module: cosmic.py
Author: Marc Zepeda
Created: 2024-09-16
Description: Catalogue Of Somatic Mutations In Cancer

Usage:
[COSMIC database]
- mutations(): returns COSMIC mutations dataframe for a given gene
- prevalence(): returns list of mutations sorted by prevalence on COSMIC

[Prime editing]
- priority_muts: returns the shared sequences library dataframe with priority mutations
- priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
'''
# Import packages
import pandas as pd
import re
import ast
from ..gen import io

# COSMIC database
def mutations(pt:str):
    ''' 
    mutations(): returns COSMIC mutations dataframe for a given gene
    
    Parameters:
    pt (str): path to COSMIC csv file
    
    Dependencies: re & io
    '''
    gene = io.get(pt)
    gene = gene[gene['Type']!='Unknown']

    befores = []
    nums = []
    afters = []
    aa_muts = []
    for aa_mut in gene['AA Mutation']:
        mut = aa_mut.split('.')[1]
        aa_muts.append(mut)
        num = int(re.findall(r'-?\d+',mut)[0])
        nums.append(num)
        befores.append(mut.split(str(num))[0])
        afters.append(mut.split(str(num))[1])
    gene['AA_position']=nums
    gene['AA_before']=befores
    gene['AA_after']=afters
    gene['AA_mut']=aa_muts

    return gene

def prevalence(gene: pd.DataFrame):
    ''' 
    prevalence(): returns list of mutations sorted by prevalence on COSMIC
    
    Parameters:
    gene (dataframe): COSMIC mutations dataframe
    
    Dependencies: pandas
    '''
    return list(gene['AA_mut'].value_counts().keys())

# Prime editing
def priority_muts(pegRNAs: pd.DataFrame, pegRNAs_shared: pd.DataFrame, pt: str):
    ''' 
    priority_muts: returns the shared sequences library dataframe with priority mutations
    
    Parameters:
    pegRNAs (dataframe): pegRNAs library dataframe
    pegRNAs_shared (dataframe): pegRNAs shared sequences library dataframe
    pt (str): path to COSMIC csv file
    
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
    pt (dataframe): path to COSMIC csv file
    gene (dataframe): COSMIC mutations dataframe for the gene from mutations()
    
    Dependencies: pandas, ast, & prevalence()
    '''
    # Get COSMIC mutations dataframe from mutations()
    gene=mutations(pt=pt)
    
    # Determine priority pegRNAs based on priority mutations from pegRNAs shared sequences library
    pegRNAs_priority = pd.DataFrame()
    for p,priority in enumerate(pegRNAs_shared['Priority_mut']):
        pegRNAs_temp = pegRNAs[pegRNAs['Edit']==priority] # Isolate priority mutations
        pegRNAs_temp = pegRNAs_temp[(pegRNAs_temp['Spacer_sequence']==pegRNAs_shared.iloc[p]['Spacer_sequence'])&(pegRNAs_temp['PBS_sequence']==pegRNAs_shared.iloc[p]['PBS_sequence'])] # Confirm spacer & PBS matches
        pegRNAs_temp.drop_duplicates(subset=['Spacer_sequence'],inplace=True) # Drop redundant pegRNAs (not sure if this is needed)
        pegRNAs_priority = pd.concat([pegRNAs_priority,pegRNAs_temp]).reset_index(drop=True)
    pegRNAs_priority['COSMIC_count'] = [gene['AA_mut'].value_counts()[edit] if edit in gene['AA_mut'].to_list() else 0 for edit in pegRNAs_priority['Edit']]

    return pegRNAs_priority