### fastq.py ###
# Author: Marc Zepeda
# Date: 2024-08-05

# Import packages
from Bio.Seq import Seq
from Bio import SeqIO
import gzip
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ..gen import io
from ..gen import tidy as t
from ..gen import plot as p

# Input/Output methods
''' parse_input: Returns flank5, wt, and flank3 from PrimeDesign input file
        pt: path to PrimeDesign input file
    Dependencies: io
'''
def parse_input(pt: str):
    target_sequence = io.get(pt=pt).iloc[0]['target_sequence']
    flank5 = target_sequence.split('(')[0][-10:]
    wt = target_sequence.split('(')[1].split(')')[0]
    flank3 = target_sequence.split(')')[1][:10]
    return flank5,wt,flank3

''' revcom_fastqs: Write reverse complement of fastqs to a new directory
        in_dir: directory with fastq files
        out_dir: new directory with reverse complement fastq files
    Dependencies: Bio.SeqIO, gzip, os, Bio.Seq.Seq
'''
def revcom_fastqs(in_dir: str, out_dir: str):

    os.makedirs(out_dir, exist_ok=True) # Ensure the output directory exists

    for filename in os.listdir(in_dir): # Find all .fastq.gz & .fastq files in the input directory

        if filename.endswith(".fastq.gz"):
            input_fastq_gz = os.path.join(in_dir, filename)
            output_fastq_gz = os.path.join(out_dir, filename)
            print(f"Processing {filename}...")
            with gzip.open(input_fastq_gz, "rt") as infile, gzip.open(output_fastq_gz, "wt") as outfile:
                for record in SeqIO.parse(infile, "fastq"):
                    reverse_complement_seq = record.seq.reverse_complement() # Compute the reverse complement of the sequence
                    reverse_complement_record = record[:] # Create a new record with the reverse complement sequence
                    reverse_complement_record.seq = reverse_complement_seq # Write the new record to the output file
                    SeqIO.write(reverse_complement_record, outfile, "fastq")
            print(f"Saved reverse complement to {output_fastq_gz}")
        
        elif filename.endswith(".fastq"):
            input_fastq = os.path.join(in_dir, filename)
            output_fastq = os.path.join(out_dir, filename)
            print(f"Processing {filename}...")
            with open(input_fastq, "r") as infile, open(output_fastq, "w") as outfile:
                for record in SeqIO.parse(infile, "fastq"):
                    reverse_complement_seq = record.seq.reverse_complement() # Compute the reverse complement of the sequence
                    reverse_complement_record = record[:] # Create a new record with the reverse complement sequence
                    reverse_complement_record.seq = reverse_complement_seq # Write the new record to the output file
                    SeqIO.write(reverse_complement_record, outfile, "fastq")
            print(f"Saved reverse complement to {output_fastq_gz}")

''' get_fastqs: Get fastq files from directory and store records in dataframes in a dictionary
        dir: directory with fastq files
        suf: file suffix (.fastq.gz or .fastq)
        quality: phred quality score threshold
        compressed: zipped file(s)?
    Dependencies: Bio.SeqIO, gzip, os, pandas, Bio.Seq.Seq
'''
def get_fastqs(dir: str,suf='.fastq.gz',quality=0,compressed=True):
    
    # Obtain fastq files
    files = os.listdir(dir)
    fastq_files = [file for file in files if suf in file]
    
    # Make fastqs dictionary
    fastqs = dict()
    for fastq_file in fastq_files:
        if compressed==True: # Compressed fastq files
            ids=[]
            seqs=[]
            phred_scores=[]
            with gzip.open(os.path.join(dir,fastq_file), "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    if all(score >= quality for score in record.letter_annotations['phred_quality']):
                        ids.append(record.id)
                        seqs.append(record.seq)
                        phred_scores.append(record.letter_annotations['phred_quality'])
        else: # Uncompressed fastq files
            with open(os.path.join(dir,fastq_file), "r") as handle:
                ids=[]
                seqs=[]
                phred_scores=[]
                for record in SeqIO.parse(handle, "fastq"):
                    if all(score >= quality for score in record.letter_annotations['phred_quality']):
                        ids.append(record.id)
                        seqs.append(record.seq)
                        phred_scores.append(record.letter_annotations['phred_quality'])
        fastqs[fastq_file[:-len(suf)]]=pd.DataFrame({'id':ids, # Add dataframe to dictionary 
                                                     'seq':seqs,
                                                     'phred_scores':phred_scores})
        print(f'{fastq_file[:-len(suf)]}:\t{len(fastqs[fastq_file[:-len(suf)]])} reads')
    return fastqs

# Determine editing outcomes & distribution methods
''' filter_fastqs: Gets DNA and AA sequence for records within flanks
        fastqs: dictionary from get_fastqs
        pt: path to PrimeDesign input file (Required unless flank5 and flank3 are provided)
        flank5: top strand flanking sequence 5' (Required unless pt is provided)
        flank3: top strand flanking sequence 3' (Required unless pt is provided)
        quality: phred quality score threshold within flanks
    Dependencies: get_fastqs(), pandas, Bio.Seq.Seq
'''
def filter_fastqs(fastqs: dict, pt='', flank5='', flank3='', quality=0):

    # Obtain flank5 and flank3 from pt or check that flank5 and flank3 have been provided
    if pt!='': (flank5,wt,flank3) = parse_input(pt)
    elif (flank5=='')&(flank3==''):
        raise ValueError('pt or flank3 and flank5 must be provided.')
    
    # Remove fastq records that do not have flanks
    fastqs_1=dict()
    for file,fastq in fastqs.items():
        missing_flank = []
        for i,seq in enumerate(fastq['seq']):
            if (seq.find(flank5)==-1)|(seq.find(flank3)==-1): missing_flank.append(i)
        fastqs_1[file] = fastq.drop(missing_flank).reset_index(drop=True)
     
    # Obtain nucleotide and AA sequences within flanks; remove fastq records with phred scores within flanks
    for file,fastq in fastqs_1.items():
        nuc=[]
        prot=[]
        low_phred=[]
        for i,seq in enumerate(fastq['seq']):
            nuc.append(seq[seq.find(flank5)+len(flank5):seq.find(flank3)])
            prot.append(Seq.translate(seq[seq.find(flank5)+len(flank5):seq.find(flank3)]))
            if quality>min(fastq.iloc[i]['phred_scores'][seq.find(flank5)+len(flank5):seq.find(flank3)]): low_phred.append(i)
        fastqs_1[file]['nuc']=nuc
        fastqs_1[file]['prot']=prot
        txt = f'{file}:\t{len(fastqs_1[file])} reads\t=>\t'
        fastqs_1[file] = fastq.drop(low_phred).reset_index(drop=True)
        print(txt + f'{len(fastqs_1[file])} reads')
    
    return fastqs_1

''' genotype: Assign genotypes to sequence records
        fastqs: dictionary from filter_fastqs
        res: First AA number
        pt: path to PrimeDesign input file (Required unless wt is provided)
        wt: Expected wildtype nucleotide sequence (in frame AA; required unless pt is provided)
    Dependencies: get_fastqs(),filter_fastqs(), pandas, Bio.Seq.Seq
'''
def genotype(fastqs: dict, res: int, pt='', wt=''):

    # Obtain wt from pt or check that wt have been provided
    if pt!='': (flank5,wt,flank3) = parse_input(pt)
    elif wt=='':
        raise ValueError('pt or wt must be provided.')

    for file,fastq in fastqs.items():
        edit=[]
        for i in range(len(fastq['prot'])):
            if len(wt)!=len(fastq.iloc[i]['nuc']): edit.append('Indel')
            elif Seq.translate(Seq(wt))==fastq.iloc[i]['prot']: edit.append('WT')
            else:
                e = []
                for j, (a, b) in enumerate(zip(Seq.translate(Seq(wt)), fastq.iloc[i]['prot'])):
                    if a != b: e.append(a+str(j+res)+b)
                if len(e)>1: edit.append(", ".join(e))
                elif len(e)==1: edit.append(e[0])
                else: edit.append('Unknown Edit')
        fastqs[file]['edit']=edit
    return fastqs

''' outcomes: Returns edit count & fraction per sample (tidy format)
        fastqs: dictionary from genotype
    Dependencies: get_fastqs(),filter_fastqs(),genotype(),pandas
'''
def outcomes(fastqs: dict):
    df = pd.DataFrame()
    for file,fastq in fastqs.items():
        temp=pd.DataFrame({'sample':[file]*len(fastq['edit'].value_counts()),
                           'edit':list(fastq['edit'].value_counts().keys()),
                           'count':fastq['edit'].value_counts(),
                           'fraction':fastq['edit'].value_counts()/len(fastq['edit'])})
        df=pd.concat([df,temp]).reset_index(drop=True)
    return df

# Generate DMS grid methods
''' aa_props: Dictionary of AA properties with citations (Generated by ChatGPT)
        Hydrophobicity: https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html
        Weight: from Bio.Data import IUPACData (protein_weights)
        Polarity: https://web.expasy.org/protscale/pscale/PolarityGrantham.html
'''
aa_props = {
    'E': {'name': 'Glutamic acid', 'hydrophobicity': -3.5, 'weight': 147.1, 'polarity': 12.3, 'charge': 'negative'},
    'D': {'name': 'Aspartic acid', 'hydrophobicity': -3.5, 'weight': 133.1, 'polarity': 13.0, 'charge': 'negative'},
    'R': {'name': 'Arginine', 'hydrophobicity': -4.5, 'weight': 174.2, 'polarity': 10.5, 'charge': 'positive'},
    'H': {'name': 'Histidine', 'hydrophobicity': -3.2, 'weight': 155.2, 'polarity': 10.4, 'charge': 'positive'},
    'K': {'name': 'Lysine', 'hydrophobicity': -3.9, 'weight': 146.2, 'polarity': 11.3, 'charge': 'positive'},
    'F': {'name': 'Phenylalanine', 'hydrophobicity': 2.8, 'weight': 165.2, 'polarity': 5.2, 'charge': 'neutral'},
    'Y': {'name': 'Tyrosine', 'hydrophobicity': -1.3, 'weight': 181.2, 'polarity': 6.2, 'charge': 'neutral'},
    'W': {'name': 'Tryptophan', 'hydrophobicity': -0.9, 'weight': 204.2, 'polarity': 5.4, 'charge': 'neutral'},
    'S': {'name': 'Serine', 'hydrophobicity': -0.8, 'weight': 105.1, 'polarity': 9.2, 'charge': 'neutral'},
    'Q': {'name': 'Glutamine', 'hydrophobicity': -3.5, 'weight': 146.2, 'polarity': 10.5, 'charge': 'neutral'},
    'T': {'name': 'Threonine', 'hydrophobicity': -0.7, 'weight': 119.1, 'polarity': 8.6, 'charge': 'neutral'},
    'N': {'name': 'Asparagine', 'hydrophobicity': -3.5, 'weight': 132.1, 'polarity': 11.6, 'charge': 'neutral'},
    'C': {'name': 'Cysteine', 'hydrophobicity': 2.5, 'weight': 121.2, 'polarity': 5.5, 'charge': 'neutral'},
    'P': {'name': 'Proline', 'hydrophobicity': -1.6, 'weight': 115.1, 'polarity': 8.0, 'charge': 'neutral'},
    'A': {'name': 'Alanine', 'hydrophobicity': 1.8, 'weight': 89.1, 'polarity': 8.1, 'charge': 'neutral'},
    'G': {'name': 'Glycine', 'hydrophobicity': -0.4, 'weight': 75.1, 'polarity': 9.0, 'charge': 'neutral'},
    'M': {'name': 'Methionine', 'hydrophobicity': 1.9, 'weight': 149.2, 'polarity': 5.7, 'charge': 'neutral'},
    'V': {'name': 'Valine', 'hydrophobicity': 4.2, 'weight': 117.1, 'polarity': 5.9, 'charge': 'neutral'},
    'I': {'name': 'Isoleucine', 'hydrophobicity': 4.5, 'weight': 131.2, 'polarity': 5.2, 'charge': 'neutral'},
    'L': {'name': 'Leucine', 'hydrophobicity': 3.8, 'weight': 131.2, 'polarity': 4.9, 'charge': 'neutral'},
    '*': {'name': 'Stop', 'hydrophobicity': None, 'weight': None, 'polarity': None, 'charge': None}
}

''' edit_1: Split edit column to before, after, and amino acid number
        df: Fastq outcomes
        col: Edit column name
    Dependencies: pandas
'''
def edit_1(df: pd.DataFrame,col='edit'):
    df_1 = df[(df[col].str.contains(',')==False)&(df[col]!='WT')&(df[col]!='Indel')] # Isolate single AA changes
    df_1['before']=df_1[col].str[0] # Split edit information
    df_1['after']=df_1[col].str[-1]
    df_1['number']=df_1[col].str[1:-1].astype(int)
    return df_1.reset_index(drop=True)

''' dms_tidy: Loads DMS grid data in tidy format split by condition
        df: Fastq outcomes
        cond: Condition column name for splicing fastq outcomes
        wt: Expected wildtype nucleotide sequence (in frame AA)
        res: First AA number
    Dependencies: Bio.Seq.Seq,pandas,numpy,tidy,get_fastqs(),filter_fastqs(),genotype(),edit_1(),aa_props
'''
def dms_tidy(df: pd.DataFrame, cond: str, wt:str, res: str):
    
    wt_prot = Seq(wt).translate(table=1) # Obtain WT protein sequence
    wt_nums = np.arange(res,res+len(wt_prot))
    dc=t.split(edit_1(df),cond) # Isolate single aa change fastq outcomes and split by condition
    
    dc2=dict() # Fill with DMS grid data in tidy format split by condition
    for key,df_cond in dc.items():
        
        wt_fastq = df[(df['edit']=='WT')&(df[cond]==key)] # Obtain WT fastq outcome
        df_cond_DMS=pd.DataFrame(columns=wt_fastq.columns) # Fill with DMS grid data in tidy format
        
        for num in wt_nums: # Iterate through WT protein sequence
            vals=dict() # Create dictionary with all amino acid changes for a given residue
            
            # Add metadata that is the same for all genotypes
            meta = [x for x in df_cond.columns if x not in ['edit','count','fraction','before','after','number']]
            for m in meta: 
                vals[m]=[wt_fastq[m].to_list()[0]]*len(list(aa_props.keys()))
            
            # Create all amino acid changes
            vals['before']=[wt_prot[num-142]]*len(list(aa_props.keys()))
            vals['number']=[num]*len(list(aa_props.keys()))
            vals['after']=list(aa_props.keys())
            vals['edit']=[vals['before'][i]+str(num)+vals['after'][i] for i in range(len(vals['after']))]

            # Fill in counts and fractions for amino acid changes, WT, and none
            counts=[]
            fractions=[]
            num_mut = df_cond[df_cond['number']==num]
            for a in vals['after']:
                if a == wt_prot[num-142]: # Wild type
                    counts.extend(wt_fastq['count'].to_list())
                    fractions.extend(wt_fastq['fraction'].to_list())
                elif a in num_mut['after'].to_list(): # Amino acid change present
                    counts.extend(num_mut[num_mut['after']==a]['count'].to_list())
                    fractions.extend(num_mut[num_mut['after']==a]['fraction'].to_list())
                else: # Amino acid change absent
                    counts.append(0)
                    fractions.append(0)
            vals['count']=counts
            vals['fraction']=fractions
            
            df_cond_DMS = pd.concat([df_cond_DMS,pd.DataFrame(vals)]).reset_index(drop=True) # Append residue DMS data
        df_cond_DMS['number']=df_cond_DMS['number'].astype(int) # Set number as type int
        df_cond_DMS['count']=df_cond_DMS['count'].astype(float) # Set count as type float for plotting
        dc2[key]=df_cond_DMS # Append condition DMS data
    return dc2

# Plot editing outcomes distribution & DMS grid methods
''' subscript: Returns dataframe with subscripts to tick labels
        df: Dataframe
        tick: New tick label
        tick_sub: Previous numeric tick label that will become a subscript
'''
def subscript(df: pd.DataFrame,tick='before',tick_sub='number'):
    ticks = []
    labels = []
    for (t,ts) in set(zip(df[tick],df[tick_sub])):
        ticks.append(ts)
        labels.append('$\\mathrm{'+t+'_{'+str(ts)+'}}$')
    return pd.DataFrame({'tick':ticks,'label':labels}).sort_values(by='tick').reset_index(drop=True)

''' cat: Creates category dependent graphs.
        typ: plot type (bar, box, violin, swarm, strip, point, count, bar_swarm, box_swarm, violin_swarm)
        df: tidy dataframe
        x: x-axis column
        y: y-axis column
    Dependencies: os, matplotlib, seaborn, plot.py
'''
def cat(typ: str,df: pd.DataFrame,x: str,y: str,errorbar=None,cols=None,cols_ord=None,cutoff=0.01,cols_exclude=None,
        file=None,dir=None,palette_or_cmap='colorblind',edgecol='black',lw=1,
        figsize=(10,6),title='',title_size=18,title_weight='bold',
        x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,1),x_ticks_rot=0,xticks=[],
        y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,1),y_ticks_rot=0,yticks=[],
        legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0), 
        **kwargs):
    
    # Omit data smaller than cutoff or excluded
    df_cut=df[df[y]>=cutoff]
    df_other=df[df[y]<cutoff]
    for sample in list(df_other['sample'].value_counts().keys()):
        df_temp = df_other[df_other['sample']==sample]
        df_temp['fraction']=sum(df_temp['fraction'])
        df_temp['edit']='Other'
        df_cut = pd.concat([df_cut,df_temp.iloc[0].to_frame().T])
    
    # Omit excluded
    if type(cols_exclude)==list: 
        for exclude in cols_exclude: df_cut=df_cut[df_cut[cols]!=exclude]
    elif type(cols_exclude)==str: df_cut=df_cut[df_cut[cols]!=cols_exclude]

    # Sort data by genotype position
    if cols_ord==None:
        genotypes = list(df_cut[cols].value_counts().keys())
        positions = list()
        for geno in genotypes:
            numbers = re.findall(r'\d+\.?\d*', geno)
            if len(numbers)==0: positions.append(100000) # Places WT and Indel at the end
            else: positions.append(sum([int(n) for n in numbers])/len(numbers))
        assign = pd.DataFrame({'positions':positions,
                               'genotypes':genotypes})
        cols_ord = list(assign.sort_values(by='positions')['genotypes'])

    p.cat(typ=typ,df=df_cut,x=x,y=y,errorbar=errorbar,cols=cols,cols_ord=cols_ord,cols_exclude=None,
          file=file,dir=dir,palette_or_cmap=palette_or_cmap,edgecol=edgecol,lw=lw,
          figsize=figsize,title=title,title_size=title_size,title_weight=title_weight,
          x_axis=x_axis,x_axis_size=x_axis_size,x_axis_weight=x_axis_weight,x_axis_scale=x_axis_scale,x_axis_dims=x_axis_dims,x_ticks_rot=x_ticks_rot,xticks=xticks,
          y_axis=y_axis,y_axis_size=y_axis_size,y_axis_weight=y_axis_weight,y_axis_scale=y_axis_scale,y_axis_dims=y_axis_dims,y_ticks_rot=y_ticks_rot,yticks=yticks,
          legend_title=legend_title,legend_title_size=legend_title_size,legend_size=legend_size,legend_bbox_to_anchor=legend_bbox_to_anchor,legend_loc=legend_loc,legend_items=legend_items, 
          **kwargs)

''' scat: Creates scatter plot related graphs.
        typ: plot type (scat, line, line_scat)
        df: tidy dataframe
        x: x-axis column
        y: y-axis column
    Dependencies: os, matplotlib, seaborn, plot.py
'''
def scat(typ: str,df: pd.DataFrame,x: str,y: str,cols=None,cols_ord=None,stys=None,cutoff=0.01,cols_exclude=None,
         file=None,dir=None,palette_or_cmap='colorblind',edgecol='black',
         figsize=(10,6),title='',title_size=18,title_weight='bold',
         x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,100),x_ticks_rot=0,xticks=[],
         y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,100),y_ticks_rot=0,yticks=[],
         legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0),
         **kwargs):

    # Omit data smaller than cutoff or excluded
    df_cut=df[df[y]>=cutoff]
    df_other=df[df[y]<cutoff]
    for sample in list(df_other['sample'].value_counts().keys()):
        df_temp = df_other[df_other['sample']==sample]
        df_temp['fraction']=sum(df_temp['fraction'])
        df_temp['edit']='Other'
        df_cut = pd.concat([df_cut,df_temp.iloc[0].to_frame().T])

    # Omit excluded
    if type(cols_exclude)==list: 
        for exclude in cols_exclude: df_cut=df_cut[df_cut[cols]!=exclude]
    elif type(cols_exclude)==str: df_cut=df_cut[df_cut[cols]!=cols_exclude]

    # Sort data by genotype position
    if cols_ord==None:
        genotypes = list(df_cut[cols].value_counts().keys())
        positions = list()
        for geno in genotypes:
            numbers = re.findall(r'\d+\.?\d*', geno)
            if len(numbers)==0: positions.append(100000) # Places WT and Indel at the end
            else: positions.append(sum([int(n) for n in numbers])/len(numbers))
        assign = pd.DataFrame({'positions':positions,
                               'genotypes':genotypes})
        cols_ord = list(assign.sort_values(by='positions')['genotypes'])

    p.scat(typ=typ,df=df_cut,x=x,y=y,cols=cols,cols_ord=cols_ord,cols_exclude=None,
           file=file,dir=dir,palette_or_cmap=palette_or_cmap,edgecol=edgecol,
           figsize=figsize,title=title,title_size=title_size,title_weight=title_weight,
           x_axis=x_axis,x_axis_size=x_axis_size,x_axis_weight=x_axis_weight,x_axis_scale=x_axis_scale,x_axis_dims=x_axis_dims,x_ticks_rot=x_ticks_rot,xticks=xticks,
           y_axis=y_axis,y_axis_size=y_axis_size,y_axis_weight=y_axis_weight,y_axis_scale=y_axis_scale,y_axis_dims=y_axis_dims,y_ticks_rot=y_ticks_rot,yticks=yticks,
           legend_title=legend_title,legend_title_size=legend_title_size,legend_size=legend_size,legend_bbox_to_anchor=legend_bbox_to_anchor,legend_loc=legend_loc,legend_items=legend_items, 
           **kwargs)

''' stack: Creates stacked bar plot
        df: dataframe from outcomes
        x: x-axis column
        y: y-axis column
        cols: colors column
        cutoff: y-axis values needs be greater than (ex: 1%)
    Dependencies: plot.py,re,os,pandas,numpy,matplotlib.pyplot
'''
def stack(df: pd.DataFrame,x='sample',y='fraction',cols='edit',cutoff=0.01,cols_ord=[],x_ord=[],
          file=None,dir=None,cmap='Set2',
          title='Editing Outcomes',title_size=18,title_weight='bold',
          figsize=(10,6),x_axis='',x_axis_size=12,x_axis_weight='bold',x_ticks_rot=45,x_ticks_ha='right',
          y_axis='',y_axis_size=12,y_axis_weight='bold',y_ticks_rot=0,
          legend_title='',legend_title_size=12,legend_size=12,
          legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_ncol=1,**kwargs):
    
    # Omit smaller than cutoff and convert it to other
    df_cut=df[df[y]>=cutoff]
    df_other=df[df[y]<cutoff]
    for sample in list(df_other['sample'].value_counts().keys()):
        df_temp = df_other[df_other['sample']==sample]
        df_temp['fraction']=sum(df_temp['fraction'])
        df_temp['edit']='Other'
        df_cut = pd.concat([df_cut,df_temp.iloc[0].to_frame().T])

    # Sort pivot table columns by genotype position
    if cols_ord==[]:
        genotypes = list(df_cut[cols].value_counts().keys())
        positions = list()
        for geno in genotypes:
            numbers = re.findall(r'\d+\.?\d*', geno)
            if len(numbers)==0: positions.append(100000) # Places WT and Indel at the end
            else: positions.append(sum([int(n) for n in numbers])/len(numbers))
        assign = pd.DataFrame({'positions':positions,
                               'genotypes':genotypes})
        cols_ord = list(assign.sort_values(by='positions')['genotypes'])
    
    # Make stacked barplot
    p.stack(df=df_cut,x=x,y=y,cols=cols,cutoff=0,cols_ord=cols_ord,x_ord=x_ord,
            file=file,dir=dir,cmap=cmap,
            title=title,title_size=title_size,title_weight=title_weight,
            figsize=figsize,x_axis=x_axis,x_axis_size=x_axis_size,x_axis_weight=x_axis_weight,x_ticks_rot=x_ticks_rot,x_ticks_ha=x_ticks_ha,
            y_axis=y_axis,y_axis_size=y_axis_size,y_axis_weight=y_axis_weight,y_ticks_rot=y_ticks_rot,
            legend_title=legend_title,legend_title_size=legend_title_size,legend_size=legend_size,
            legend_bbox_to_anchor=legend_bbox_to_anchor,legend_loc=legend_loc,legend_ncol=legend_ncol,**kwargs)

''' dms_grid: Creates amino acide by residue number heatmaps from tidy-formatted DMS data
        dc: Tidy-formatted DMS dictionary
        x: x-axis (AA residues number column)
        y: y-axis (AA change column)
        vals: Values column
    Dependencies: matplotlib,seaborn,aa_props
'''
def dms_grid(dc: dict,x='number',y='after',vals='count',
             file=None,dir=None,edgecol='black',lw=1,annot=False,cmap="bone_r",
             title='',title_size=12,title_weight='bold',
             x_axis='',x_axis_size=12,x_axis_weight='bold',x_ticks_rot=45,
             y_axis='',y_axis_size=12,y_axis_weight='bold',y_ticks_rot=0,
             **kwargs):
    
    # Make DMS grids
    dc2={key:pd.pivot(df_cond,columns=x,index=y,values=vals).astype(float).reindex(list(aa_props.keys())) for key,df_cond in dc.items()}

    # Create a single figure with multiple heatmap subplots
    fig, axes = plt.subplots(nrows=len(list(dc2.keys())),ncols=1,figsize=(20,7*len(list(dc2.keys()))),sharex=False,sharey=True)
    for (ax, key) in zip(axes, list(dc2.keys())):
        sns.heatmap(dc2[key],annot=annot,cmap=cmap,ax=ax,linecolor=edgecol,linewidths=lw,cbar=True)
        if len(list(dc2.keys()))>1: ax.set_title(key,fontsize=title_size,fontweight=title_weight)  # Add title to subplot
        else: ax.set_title(title,fontsize=title_size,fontweight=title_weight)
        if x_axis=='': ax.set_xlabel(p.re_un_cap(x),fontsize=x_axis_size,fontweight=x_axis_weight) # Add x axis label
        else: ax.set_xlabel(x_axis,fontsize=x_axis_size,fontweight=x_axis_weight)
        if y_axis=='': ax.set_ylabel(p.re_un_cap(y),fontsize=y_axis_size,fontweight=y_axis_weight) # Add y axis label
        else: ax.set_ylabel(y_axis,fontsize=y_axis_size,fontweight=y_axis_weight)
        ax.set_xticklabels(subscript(dc[key])['label'].to_list()) # Change x ticks to have subscript format
        plt.setp(ax.get_xticklabels(), rotation=x_ticks_rot, va='center', ha="right",rotation_mode="anchor") # Format x ticks
        plt.setp(ax.get_yticklabels(), rotation=y_ticks_rot, va='center', ha="right",rotation_mode="anchor") # Format y ticks
        ax.set_facecolor('white')  # Set background to transparent

    # Save & show fig
    if file is not None and dir is not None:
        if not os.path.exists(dir):
            os.mkdir(dir)
        plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight')
    plt.show()
