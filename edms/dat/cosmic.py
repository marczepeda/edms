### cosmic.py ###
# Author: Marc Zepeda
# Date: 2024-09-16

# Import packages
import pandas as pd
import re
import ast
from ..gen import io

# COSMIC database methods
''' mutations: Returns COSMIC mutations dataframe for a given gene.
        pt: path to COSMIC csv file
    Dependencies: re,pyMUZ.gen.io
'''
def mutations(pt:str):

    gene = io.get(pt)

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

''' prevalence: Returns list of mutations sorted by prevalence on COSMIC
        gene: COSMIC mutations dataframe
    Dependencies: pandas
'''
def prevalence(gene: pd.DataFrame):
    return list(gene['AA_mut'].value_counts().keys())

# Prime editing methods
''' priority_edits: Returns dataframe of the most clinically-relevant prime edits to prioritize from shared sequences library
        pegRNAs: pegRNAs library dataframe
        pegRNAs_shared: pegRNAs shared sequences library dataframe
        gene: COSMIC mutations dataframe for the gene from mutations()
    Dependencies: pandas,ast,prevalence()
'''
def priority_edits(pegRNAs: pd.DataFrame, pegRNAs_shared: pd.DataFrame, gene: pd.DataFrame):
    
    # Get list of priority mutants from prevalence()
    mut_priority_ls = prevalence(gene=gene)

    # Determine priority mutations for pegRNAs shared sequences library
    priority_muts = []
    for edits in mut_priority_ls['Edits']:
        for mutant in mut_priority_ls:
            if mutant in set(ast.literal_eval(edits)):
                priority_muts.append(mutant)
                break
    pegRNAs_shared['Priority_mut']=priority_muts

    # Determine priority pegRNAs based on priority mutations from pegRNAs shared sequences library
    pegRNAs_priority = pd.DataFrame()
    for p,priority in enumerate(pegRNAs_shared['Priority_mut']):
        pegRNAs_temp = pegRNAs[pegRNAs['Edit']==priority] # Isolate priority mutations
        pegRNAs_temp = pegRNAs_temp[(pegRNAs_temp['Spacer_sequence']==pegRNAs_shared.iloc[p]['Spacer_sequence'])&(pegRNAs_temp['PBS_sequence']==pegRNAs_shared.iloc[p]['PBS_sequence'])] # Confirm spacer & PBS matches
        pegRNAs_temp.drop_duplicates(subset=['Spacer_sequence'],inplace=True) # Drop redundant pegRNAs (not sure if this is needed)
        pegRNAs_priority = pd.concat([pegRNAs_priority,pegRNAs_temp]).reset_index(drop=True)
    pegRNAs_priority['COSMIC_count'] = [gene['AA_mut'].value_counts()[edit] for edit in pegRNAs_priority['Edit']]

    return pegRNAs_priority