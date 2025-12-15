''' 
Module: cosmic.py
Author: Marc Zepeda
Created: 2024-09-16
Description: Catalogue Of Somatic Mutations In Cancer

Usage:
[COSMIC database]
- mutations(): returns COSMIC mutations dataframe for a given gene
- prevalence(): returns list of mutations sorted by prevalence on COSMIC
- cds_group(): plot COSMIC mutations histogram with CDS regions highlighted in different colors

[Prime editing]
- priority_muts(): returns the shared sequences library dataframe with priority mutations
- priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library

[Base & prime editing accessible mutations]
- editor_mutations(): returns and plots editor accessible COSMIC mutations
'''
# Import packages
import pandas as pd
import re
import ast
from ..gen import io
from ..gen import tidy as t
from ..gen import plot as p

# COSMIC database
def mutations(df: pd.DataFrame | str, dir:str=None, file:str=None) -> pd.DataFrame:
    ''' 
    mutations(): returns COSMIC mutations dataframe for a given gene
    
    Parameters:
    df (dataframe | str): COSMIC dataframe (or file path)
    dir (str, optional): save directory
    file (str, optional): save file
    
    Dependencies: re & io
    '''
    # Isolate mutations corresponding to gene
    if type(df)==str: # Get COSMIC dataframe from file path if needed
        df = io.get(pt=df)
    df = df[df['Type']!='Unknown']

    befores = []
    nums = []
    afters = []
    aa_muts = []
    for aa_mut in df['AA Mutation']:
        mut = aa_mut.split('.')[1]
        aa_muts.append(mut)
        num = int(re.findall(r'-?\d+',mut)[0])
        nums.append(num)
        befores.append(mut.split(str(num))[0])
        afters.append(mut.split(str(num))[1])
    df['AA_position']=nums
    df['AA_before']=befores
    df['AA_after']=afters
    df['AA_mut']=aa_muts

    # Save & return dataframe
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=df) 
    return df

def prevalence(df: pd.DataFrame) -> list:
    ''' 
    prevalence(): returns list of mutations sorted by prevalence on COSMIC
    
    Parameters:
    df (dataframe): COSMIC mutations dataframe
    
    Dependencies: pandas
    '''
    return list(df['AA_mut'].value_counts().keys())


def cds_group(df_cosmic: pd.DataFrame | str, df_cds: pd.DataFrame | str, out_dir: str=None, **plot_kwargs):
    '''
    cds_group(): plot COSMIC mutations histogram with CDS regions highlighted in different colors

    Parameters:
    df_cosmic (str): COSMIC mutations() dataframe (or file path)
    cds_pt (str): CDS dataframe (or file path) with required columns <gene,CDS,start,end>
    out_dir (str, optional): path to output directory
    **plot_kwargs: histogram keyword arguments

    Dependencies: io,plot,os,pandas
    '''
    # Get COSMIC mutations() and CDS datarames from file paths if needed
    if type(df_cosmic)==str: 
        df_cosmic = io.get(pt=df_cosmic)
    if type(df_cds)==str: 
        df_cds = io.get(pt=df_cds)

    # Drop silent COSMIC mutations
    df_cosmic = df_cosmic[(df_cosmic['Type']!='Substitution - coding silent')].reset_index(drop=True)

    # Group COSMIC mutations by CDS region
    regions = []
    for pos in df_cosmic['AA_position']:
        for i,(region,start,end) in enumerate(t.zip_cols(df=df_cds,cols=['CDS','start','end'])):
            if (pos>=start)&(pos<=end): 
                regions.append(region)
                break
            if i==len(df_cds['CDS'])-1: # check for errors
                regions.append(0)
    df_cosmic['CDS']=regions

    # Plot histogram
    p.dist(typ='hist',df=df_cosmic,x='AA_position',bins=df_cds.iloc[-1]['end'],cols='CDS',edgecol=None,
           title='COSMIC Mutations',x_axis='Position (AA)',y_axis=f"{df_cds.iloc[0]['gene']}",
           dir=out_dir,file=f"{df_cds.iloc[0]['gene']}_CDS_group.pdf",
           **plot_kwargs)

# Prime editing
def priority_muts(pegRNAs_shared: pd.DataFrame, df_cosmic: str, dir:str=None, file:str=None) -> pd.DataFrame:
    ''' 
    priority_muts: returns the shared sequences library dataframe with priority mutations
    
    Parameters:
    pegRNAs (dataframe): pegRNAs library dataframe
    pegRNAs_shared (dataframe): pegRNAs shared sequences library dataframe
    df_cosmic (dataframe | str): COSMIC dataframe (or file path)
    dir (str, optional): save directory
    file (str, optional): save file
    
    Dependencies: pandas, ast, prevalence(), & mutations()
    '''
    # Get pegRNAs shared sequences library and ClinVar mutation() dataframes from file paths if needed
    if type(pegRNAs_shared)==str:
        pegRNAs_shared = io.get(pt=pegRNAs_shared)
    if type(df_cosmic)==str:
        df_cosmic = io.get(pt=df_cosmic)

    # Check if COSMIC mutations() dataframe has been processed
    if 'AA_mut' not in df_cosmic.columns:
        df_cosmic = mutations(df=df_cosmic)
    
    # Get list of priority mutants
    mut_priority_ls = prevalence(df=df_cosmic)

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

def priority_edits(pegRNAs: pd.DataFrame | str, pegRNAs_shared: pd.DataFrame | str, df_cosmic: pd.DataFrame | str, 
                   dir:str=None, file:str=None) -> pd.DataFrame:
    ''' 
    priority_edits(): returns a dataframe with the most clinically-relevant prime edits to prioritize from the shared sequences library
    
    Parameters:
    pegRNAs (dataframe | str): pegRNAs library dataframe (or file path)
    pegRNAs_shared (dataframe | str): pegRNAs shared sequences library dataframe (or file path)
    df_cosmic (dataframe | str): COSMIC mutations() dataframe (or file path)
    dir (str, optional): save directory
    file (str, optional): save file

    Dependencies: pandas, ast, & prevalence()
    '''
    # Get pegRNAs shared sequences library and COSMIC mutation() dataframes from file paths if needed
    if type(pegRNAs_shared)==str:
        pegRNAs_shared = io.get(pt=pegRNAs_shared)
    if type(df_cosmic)==str:
        df_cosmic = io.get(pt=df_cosmic)
    
    # Check if COSMIC mutations() dataframe has been processed
    if 'AA_mut' not in df_cosmic.columns:
        df_cosmic = mutations(df=df_cosmic)
    
    # Determine priority pegRNAs based on priority mutations from pegRNAs shared sequences library
    pegRNAs_priority = pd.DataFrame()
    for p,priority in enumerate(pegRNAs_shared['Priority_mut']):
        pegRNAs_temp = pegRNAs[pegRNAs['Edit']==priority] # Isolate priority mutations
        pegRNAs_temp = pegRNAs_temp[(pegRNAs_temp['Spacer_sequence']==pegRNAs_shared.iloc[p]['Spacer_sequence'])&(pegRNAs_temp['PBS_sequence']==pegRNAs_shared.iloc[p]['PBS_sequence'])] # Confirm spacer & PBS matches
        pegRNAs_temp.drop_duplicates(subset=['Spacer_sequence'],inplace=True) # Drop redundant pegRNAs (not sure if this is needed)
        pegRNAs_priority = pd.concat([pegRNAs_priority,pegRNAs_temp]).reset_index(drop=True)
    pegRNAs_priority['COSMIC_count'] = [df_cosmic['AA_mut'].value_counts()[edit] if edit in df_cosmic['AA_mut'].to_list() else 0 for edit in pegRNAs_priority['Edit']]
    
    # Save & return pegRNAs priority dataframe
    if dir is not None and file is not None:
        io.save(dir=dir,file=file,obj=pegRNAs_priority) 
    return pegRNAs_priority

# Base & prime editing accessible mutations
def editor_mutations(df_cosmic: pd.DataFrame | str, df_bescan: pd.DataFrame | str, out_dir: str=None, **plot_kwargs):
    '''
    editor_mutations(): returns and plots editor accessible COSMIC mutations

    Parameters:
    df_cosmic (dataframe | str): COSMIC mutations() dataframe (or file path)
    df_bescan (dataframe | str): BESCAN sgRNA library dataframe (or file path)
    out_dir (str, optional): path to output directory
    **plot_kwargs (optional): Plot keyword arguments

    Dependencies: pandas, plot
    '''
    # Get COSMIC mutation() and BESCAN sgRNA library dataframes from file paths if needed
    if type(df_cosmic)==str:
        df_cosmic = io.get(pt=df_cosmic)
    if type(df_bescan)==str:
        df_bescan = io.get(pt=df_bescan)

    # Get BE (ABE & CBE) mutations
    abe = list(df_bescan['AtoG_mutations'].value_counts().keys())
    cbe = list(df_bescan['CtoT_mutations'].value_counts().keys())

    # Retain all unique BE mutations
    abe_set = set()
    for a in abe:
        abe_set.update(re.split(r'[;/]',a))
    cbe_set = set()
    for a in cbe:
        cbe_set.update(re.split(r'[;/]',a))

    # Combine ABE and CBE mutations
    if '' in abe_set: abe_set.remove('')
    if '' in cbe_set: cbe_set.remove('')
    be_set = abe_set.union(cbe_set)

    # Sort mutations numerically
    abe_list = sorted(list(abe_set), key=lambda x: int(re.search(r'\d+', x).group()))
    cbe_list = sorted(list(cbe_set), key=lambda x: int(re.search(r'\d+', x).group()))
    be_list = sorted(list(be_set), key=lambda x: int(re.search(r'\d+', x).group()))

    # Save BE mutations
    df_bescan2 = pd.DataFrame({'gene': [df_bescan.iloc[0]['gene']],
                               'ABE': [abe_list],
                               'CBE': [cbe_list],
                               'BE': [be_list]})

    # Isolate COSMIC BE mutations
    cosmic_mut_list = list(df_cosmic['AA_mut'])
    abe_cosmic_list = [mut for mut in df_bescan2.iloc[0]['ABE'] if mut in cosmic_mut_list]
    cbe_cosmic_list = [mut for mut in df_bescan2.iloc[0]['CBE'] if mut in cosmic_mut_list]
    be_cosmic_list = [mut for mut in df_bescan2.iloc[0]['BE'] if mut in cosmic_mut_list]

    # Save COSMIC BE mutations
    df_bescan2['ABE_COMSIC']= [abe_cosmic_list]
    df_bescan2['CBE_COMSIC']= [cbe_cosmic_list]
    df_bescan2['BE_COMSIC']= [be_cosmic_list]
    if out_dir is not None:
        io.save(dir=out_dir,
                file=f"{df_bescan.iloc[0]['gene']}_BE_COSMIC.csv",
                obj=df_bescan2)


    # Quantify BE and PE COSMIC mutations 
    cosmic_nodup = df_cosmic.drop_duplicates(subset='AA_mut')
    cosmic_nodup_change = cosmic_nodup[cosmic_nodup['Type']!='Substitution - coding silent']
    cosmic_nodup_change_BE = len(cosmic_nodup_change[[aa_mut in be_cosmic_list for aa_mut in cosmic_nodup_change['AA_mut']]])

    cosmic_nodup_change_type = cosmic_nodup_change['Type'].value_counts()
    cosmic_nodup_change_PE_sub = cosmic_nodup_change_type['Substitution - Missense']+cosmic_nodup_change_type['Substitution - Nonsense']-cosmic_nodup_change_BE
    cosmic_nodup_change_PE_indel = cosmic_nodup_change_type['Insertion - Frameshift']+cosmic_nodup_change_type['Insertion - In frame']+cosmic_nodup_change_type['Deletion - Frameshift']+cosmic_nodup_change_type['Deletion - In frame']

    cosmic_nodup_change_type_complex = len(cosmic_nodup_change['Type'])-cosmic_nodup_change_BE-cosmic_nodup_change_PE_sub-cosmic_nodup_change_PE_indel

    # Save and plot editor accessibly COSMIC mutations 
    df_editor_cosmic = pd.DataFrame({'mutation': ['Substitution (ABE/CBE)','Substitution (PE)','Indel (PE)','Complex'],
                                     'COSMIC': [cosmic_nodup_change_BE,cosmic_nodup_change_PE_sub,cosmic_nodup_change_PE_indel,cosmic_nodup_change_type_complex]})
    df_editor_cosmic['fraction']=df_editor_cosmic['COSMIC']/sum(df_editor_cosmic['COSMIC'])
    df_editor_cosmic['gene']=df_bescan.iloc[0]['gene']
    
    if out_dir is not None:
        io.save(dir=out_dir,
                file=f"{df_bescan.iloc[0]['gene']}_COSMIC_types.csv",
                obj=cosmic_nodup_change_type.reset_index(drop=False))
        io.save(dir=out_dir,
                file=f"{df_bescan.iloc[0]['gene']}_editor_COSMIC.csv",
                obj=df_editor_cosmic)

        p.stack(df=df_editor_cosmic,x='gene',y='fraction',cols='mutation',
                cols_ord=['Substitution (ABE/CBE)','Substitution (PE)','Indel (PE)','Complex'],
                x_ticks_rot=0,x_ticks_ha='center',x_axis='COSMIC',
                title=df_bescan.iloc[0]['gene'],
                figsize=(1,5),dir=out_dir,
                file=f"{df_bescan.iloc[0]['gene']}_editor_COSMIC.pdf",
                **plot_kwargs)