### cvar.py ###
# Author: Marc Zepeda
# Date: 2024-09-16

# Import packages
import pandas as pd
import re
import ast
from ..gen import io

# ClinVar database methods
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

# Prime editing methods
def priority_edits(pegRNAs: pd.DataFrame, pegRNAs_shared: pd.DataFrame, gene: pd.DataFrame):
    ''' 
    priority_edits(): returns dataframe of the most clinically-relevant prime edits to prioritize from shared sequences library
    
    Note: I might want to update to resemble cosmic.py (priority_edits & priority_muts)

    Parameters:
    pegRNAs (dataframe): pegRNAs library dataframe
    pegRNAs_shared (dataframe): pegRNAs shared sequences library dataframe
    gene (dataframe): ClinVar mutations dataframe for the gene from mutations()
    
    Dependencies: pandas, ast, & prevalence()
'''
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
    pegRNAs_priority['ClinVar_count'] = [gene['Protein change'].value_counts()[edit] for edit in pegRNAs_priority['Edit']]

    return pegRNAs_priority