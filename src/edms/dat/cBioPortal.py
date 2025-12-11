'''
Module: cBioPortal.py
Author: Marc Zepeda
Created: 2024-09-16
Description: cBioPortal mutation data

Usage:
- mutations(): Process cBioPortal mutation data for a single protein-coding gene
'''
# Import packages
import pandas as pd
import re
import os

from ..gen import io

def mutations(df: pd.DataFrame | str, wt: str, config: bool=True, dir: str=None, file: str=None) -> pd.DataFrame:
    ''' 
    mutations(): Process cBioPortal mutation data for a single protein-coding gene. Retrieve cBioPortal mutation data from https://www.cbioportal.org/ (GENIE Cohort v18.0-public).

    Parameters:
    df (dataframe or str): DataFrame or path to cBioPortal mutation data TSV file
    wt (str): Wild-type amino acid sequence
    config (bool, optional): Save to configuration directory (Default: True)
    dir (str, optional): Directory to save output file (Default: None)
    file (str, optional): Output filename (Default: None)

    Dependencies: pandas, re, io
    '''
    if isinstance(df, str): # Get dataframe from path if needed
        df = io.get(df)

    # Determine number of patients with each mutation and associated cancer types (sorted by descending frequency)
    df.dropna(subset=['Gene','Cancer Type','Protein Change','Mutation Type'], inplace=True)
    df = df[df['Mutation Type']!='fusion'] # Exclude fusions
    df_cts = df[['Gene','Protein Change','Mutation Type']].value_counts().reset_index(name='counts')
    df_cts['Cancer Types'] = df_cts['Protein Change'].map(df.groupby('Protein Change')['Cancer Type'].apply(lambda x: ', '.join(x.value_counts().index)))

    # Iterate through genes
    dc_df_cts = dict()
    for gene in df_cts['Gene'].unique():
        df_cts_gene = df_cts[df_cts['Gene']==gene]

        # Modify cBioPortal protein change notation to EDMS edit notation
        edit_change = []
        for protein_change,mutation_type in t.zip_cols(df=df_cts_gene, cols=['Protein Change', 'Mutation Type']):
            
            if 'In_Frame_Del' == mutation_type: # In-frame deletion 
                nums = re.findall(r'\d+', protein_change)

                if 'delins' in protein_change and len(nums) == 2: # Amino acid deletion + insertion (EDMS notation change)
                    num1 = int(nums[0])
                    num2 = int(nums[1])
                    delins_i = protein_change.index('delins')
                    edit_change.append(f'{wt[num1-1:num2+1]}{num1}{protein_change[delins_i+6:]}{wt[num2]}')
                    if edit_change[-1][0] != protein_change[0]: # Sanity check
                        raise ValueError(f'In_Frame_Del (delins): cBioPortal protein change notation ({protein_change}) mismatch with WT protein {wt}\nProtein Change: {protein_change}\nEdit Change: {edit_change[-1]}')

                elif len(nums) == 1 and protein_change.endswith('del'): # Single amino acid deletion (EDMS notation change)
                    num = int(nums[0])
                    edit_change.append(f'{wt[num-1]}{wt[num]}{num}{wt[num]}')
                    if edit_change[-1][0] != protein_change[0]: # Sanity check
                        raise ValueError(f'In_Frame_Del (single): cBioPortal protein change notation ({protein_change}) mismatch with WT protein {wt}\nProtein Change: {protein_change}\nEdit Change: {edit_change[-1]}')
                
                elif len(nums) == 2 and protein_change.endswith('del'): # Multiple amino acid deletion (EDMS notation change)
                    num1 = int(nums[0])
                    num2 = int(nums[1])
                    edit_change.append(f'{wt[num1-1:num2+1]}{num1}{wt[num2]}')
                    if edit_change[-1][0] != protein_change[0]: # Sanity check
                        raise ValueError(f'In_Frame_Del (multiple): cBioPortal protein change notation ({protein_change}) mismatch with WT protein {wt}\nProtein Change: {protein_change}\nEdit Change: {edit_change[-1]}')

                else: # Unknown (keep cBioPortal notation)
                    edit_change.append(protein_change)

            elif 'In_Frame_Ins' == mutation_type: # In-frame insertion
                nums = re.findall(r'\d+', protein_change)

                if 'delins' in protein_change and len(nums) == 1: # Amino acid deletion + insertion (EDMS notation change)
                    num = int(nums[0])
                    delins_i = protein_change.index('delins')
                    edit_change.append(f'{wt[num-1]}{num}{protein_change[delins_i+6:]}')
                    if edit_change[-1][0] != protein_change[0]: # Sanity check
                        raise ValueError(f'In_Frame_Ins (delins): cBioPortal protein change notation ({protein_change}) mismatch with WT protein {wt}\nProtein Change: {protein_change}\nEdit Change: {edit_change[-1]}')

                elif len(nums) == 1 and protein_change.endswith('dup'): # Single amino acid duplication (EDMS notation change)
                    num = int(nums[0])
                    edit_change.append(f'{wt[num-1]}{num}{wt[num-1]}{wt[num-1]}')
                    if edit_change[-1][0] != protein_change[0]: # Sanity check
                        raise ValueError(f'In_Frame_Ins (single duplication): cBioPortal protein change notation ({protein_change}) mismatch with WT protein {wt}\nProtein Change: {protein_change}\nEdit Change: {edit_change[-1]}')

                elif len(nums) == 2 and protein_change.endswith('dup'): # Multiple amino acid duplication (EDMS notation change)
                    num1 = int(nums[0])
                    num2 = int(nums[1])
                    edit_change.append(f'{wt[num2-1]}{num2}{wt[num2-1]}{wt[num1-1:num2]}')
                    if edit_change[-1][0] != protein_change[protein_change.find(str(num2))-1]: # Sanity check
                        raise ValueError(f'In_Frame_Ins (multiple duplication): cBioPortal protein change notation ({protein_change}) mismatch with WT protein {wt}\nProtein Change: {protein_change}\nEdit Change: {edit_change[-1]}')

                elif len(nums) == 2 and 'ins' in protein_change: # Amino acid insertion (EDMS notation change)
                    num1 = int(nums[0])
                    num2 = int(nums[1])
                    ins_i = protein_change.index('ins')
                    edit_change.append(f'{wt[num1-1]}{num1}{wt[num1-1]}{protein_change[ins_i+3:]}')
                    if num2 != num1 + 1: # Sanity check
                        raise ValueError(f'In_Frame_Ins (insertion): cBioPortal protein change notation ({protein_change}) has non-consecutive amino acid positions\nProtein Change: {protein_change}\nEdit Change: {edit_change[-1]}')
                    if edit_change[-1][0] != protein_change[0]: # Sanity check
                        raise ValueError(f'In_Frame_Ins (insertion): cBioPortal protein change notation ({protein_change}) mismatch with WT protein {wt}\nProtein Change: {protein_change}\nEdit Change: {edit_change[-1]}')
                
                else: # Unknown (keep cBioPortal notation)
                    edit_change.append(protein_change)

            else: # Other mutation types (keep cBioPortal notation, which may already be equivalent to EDMS notation)
                edit_change.append(protein_change)

        df_cts_gene['Edit Change'] = edit_change

        # Per-gene dataframe 
        if config==True: # saved to config directory
            io.save(dir=os.path.expanduser("~/.config/edms/cBioPortal_mutations"), file=f'{gene}.csv',obj=df_cts_gene)
        dc_df_cts[gene] = df_cts_gene # stored in dictionary
    
    # Concatenate all genes; save & return
    df_cts = pd.concat(dc_df_cts.values(), ignore_index=True)
    if dir is not None and file is not None:
        io.save(df=df_cts, dir=dir, file=file)
    return df_cts