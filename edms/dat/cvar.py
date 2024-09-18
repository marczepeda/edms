### cvar.py ###
# Author: Marc Zepeda
# Date: 2024-09-16

# Import packages
import pandas as pd
import re
import ast
from ..gen import io

# ClinVar database methods
''' mutations: Returns ClinVar mutations dataframe for a given gene.
        gene_name: gene name
        pt: path to ClinVar tsv file
    Dependencies: re,pyMUZ.gen.io
'''
def mutations(gene_name:str,pt:str,typ:str='tsv'):

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

''' prevalence: Returns list of mutations sorted by prevalence on ClinVar
        gene: ClinVar mutations dataframe
    Dependencies: pandas
'''
def prevalence(gene: pd.DataFrame):
    return list(gene['Protein change'].value_counts().keys())

# Prime editing methods
''' priority_edits: Returns dataframe of the most clinically-relevant prime edits to prioritize from shared sequences library
        pegRNAs: pegRNAs library dataframe
        pegRNAs_shared: pegRNAs shared sequences library dataframe
        gene: ClinVar mutations dataframe for the gene from mutations()
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
    pegRNAs_priority['ClinVar_count'] = [gene['Protein change'].value_counts()[edit] for edit in pegRNAs_priority['Edit']]

    return pegRNAs_priority