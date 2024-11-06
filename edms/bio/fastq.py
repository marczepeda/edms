### fastq.py ###
# Author: Marc Zepeda
# Date: 2024-08-05

# Import packages
from Bio.Seq import Seq
from Bio import SeqIO
import gzip
import os
import re
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from ..gen import io
from ..gen import tidy as t
from ..gen import plot as p

# Get rid of warnings
import warnings
warnings.filterwarnings("ignore")

# Input/Output methods
def parse_input(pt: str):
    ''' 
    parse_input(): returns flank5, wt, and flank3 from PrimeDesign input file
    
    Parameters:
    pt (str): path to PrimeDesign input file
    
    Dependencies: io
    '''
    target_sequence = io.get(pt=pt).iloc[0]['target_sequence']
    flank5 = target_sequence.split('(')[0][-10:]
    wt = target_sequence.split('(')[1].split(')')[0]
    flank3 = target_sequence.split(')')[1][:10]
    return flank5,wt,flank3

def revcom_fastqs(in_dir: str, out_dir: str):
    ''' 
    revcom_fastqs(): write reverse complement of fastqs to a new directory
    
    Parameters:
    in_dir (str): directory with fastq files
    out_dir (str): new directory with reverse complement fastq files
    
    Dependencies: Bio.SeqIO, gzip, os, & Bio.Seq.Seq
    '''
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

def comb_fastqs(in_dir: str, out_dir: str, out_file: str):
    ''' 
    comb_fastqs(): Combine fastqs to a new directory

    Combines one or more (un)compressed fastqs files into a single (un)compressed fastq file.

    Parameters:
    in_dir (str): directory with fastq files
    out_dir (str): new directory with combined fastq file
    out_file (str): Name of output fastq file (Needs .fastq or .fastq.gz suffix)
    
    Dependencies: gzip & os
    '''
    os.makedirs(out_dir, exist_ok=True) # Ensure the output directory exists

    if in_file.endswith(".fastq.gz"):
        with gzip.open(os.path.join(out_dir,out_file), 'wt') as out:
            for in_file in os.listdir(in_dir): # Find all .fastq.gz & .fastq files in the input directory
                print(f"Processing {in_file}...")
                if in_file.endswith(".fastq.gz"):
                    with gzip.open(os.path.join(in_dir,in_file), 'rt') as f:
                        for line in f:
                            out.write(line)
                
                elif in_file.endswith(".fastq"):
                    with open(os.path.join(in_dir,in_file), 'r') as f:
                        for line in f:
                            out.write(line)
    
    elif in_file.endswith(".fastq"):
        with open(os.path.join(out_dir,out_file), 'wt') as out:
            for in_file in os.listdir(in_dir): # Find all .fastq.gz & .fastq files in the input directory
                print(f"Processing {in_file}...")
                if in_file.endswith(".fastq.gz"):
                    with gzip.open(os.path.join(in_dir,in_file), 'rt') as f:
                        for line in f:
                            out.write(line)
                
                elif in_file.endswith(".fastq"):
                    with open(os.path.join(in_dir,in_file), 'r') as f:
                        for line in f:
                            out.write(line)

    else: print('out_file needs .fastq or .fastq.gz suffix')

def trim_filter(record,qall:int,qavg:int,qtrim:int,qmask:int,alls:int,avgs:int,trims:int,masks:int):
    ''' 
    trim_filter(): trim and filter fastq sequence based on quality scores
    
    Parameters:
    record: Bio.SeqIO fastq record
    qall (int): phred quality score threshold for all bases for a read to not be discarded
    qtrim (int): phred quality score threshold for trimming reads on both ends
    qavg (int): average phred quality score threshold for a read to not be discarded
    qmask (int): phred quality score threshold for base to not be masked to N
    alls (int): count of records that were dropped due to qall threshold
    avgs (int): count of records that were dropped due to qavg threshold
    trims (int): count of records that were trimmed due to qtrim threshold
    masks (int): count of records that had bases masked due to qmask threshold
    
    Dependencies: Bio.SeqIO, gzip, os, pandas, & Bio.Seq.Seq
    '''
    if all(score >= qall for score in record.letter_annotations['phred_quality']): # All threshold
        if np.mean(record.letter_annotations['phred_quality']) >= qavg: # Avg threshold
            
            quality_scores = record.letter_annotations['phred_quality'] # Set 5' & 3' trim indexes to the start and end
            trim_5 = 0 
            trim_3 = len(quality_scores)
            sequence = record.seq
            
            for i in range(len(quality_scores)): # Find 5' trim
                if quality_scores[i] >= qtrim: break
                trim_5 = i
            for i in reversed(range(len(quality_scores))): # Find 3' trim
                if quality_scores[i] >= qtrim: break
                trim_3 = i
            if (trim_5!=0)|(trim_3!=len(quality_scores)): trims += 1 # Trimmed read

            sequence = sequence[trim_5:trim_3] # Trim the sequence and quality scores
            quality_scores = quality_scores[trim_5:trim_3]

            
            bases = list(sequence) # Mask bases with 'N' threshold
            if masks !=0:
                for i, qual in enumerate(quality_scores):
                    if qual < qmask: bases[i] = 'N'
            sequenceN = Seq('').join(bases) # Update the sequence with the modified version
            if Seq('N') in sequenceN: masks += 1

            return record.id,sequence,sequenceN,quality_scores,alls,avgs,trims,masks
    
        else: return None,None,None,None,alls,avgs+1,trims,masks # Avg threshold not met
    else: return None,None,None,None,alls+1,avgs,trims,masks # All threshold not met

def get_fastqs(dir: str,suf='.fastq.gz',qall=10,qavg=30,qtrim=30,qmask=0,compressed=True,save=True):
    ''' 
    get_fastqs(): get fastq files from directory and store records in dataframes in a dictionary
    
    Parameters:
    dir (str): directory with fastq files
    suf (str, optional): file suffix (.fastq.gz or .fastq)
    qall (int, optional): phred quality score threshold for all bases for a read to not be discarded (Q = -log(err))
    qtrim (int, optional): phred quality score threshold for trimming reads on both ends (Q = -log(err))
    qavg (int, optional): average phred quality score threshold for a read to not be discarded (Q = -log(err))
    qmask (int, optional): phred quality score threshold for base to not be masked to N (Q = -log(err))
    compressed (str, optional): zipped file(s)? (Default: True)
    save (bool, optional): save reads statistics file to local directory (DefaultL True)

    Dependencies: Bio.SeqIO, gzip, os, pandas, Bio.Seq.Seq, & trim_filter()
    '''
    # Obtain fastq files
    files = os.listdir(dir)
    fastq_files = [file for file in files if suf in file]
    
    # Make fastqs dictionary
    fastqs = dict()
    if save == True: out = pd.DataFrame()
    for fastq_file in fastq_files:
        reads = 0 # Keep track of reads & outcomes
        alls = 0
        avgs = 0
        trims = 0
        masks = 0

        if (compressed==True)|(suf=='.fastq.gz'): # Compressed fastq files
            ids=[]
            seqs=[]
            seqsN=[]
            phred_scores=[]
            with gzip.open(os.path.join(dir,fastq_file), "rt") as handle:

                for r,record in enumerate(SeqIO.parse(handle, "fastq")): # Parse reads
                    reads=r+1
                    record_id,record_seq,record_seqN,record_scores,alls,avgs,trims,masks = trim_filter(record,qall,qavg,qtrim,qmask,alls,avgs,trims,masks) # Trim & filter
                    if record_id is not None: # Save id, sequence, & quality scores
                        ids.append(record_id) 
                        seqs.append(record_seq)
                        seqsN.append(record_seqN)
                        phred_scores.append(record_scores)

        else: # Uncompressed fastq files
            ids=[]
            seqs=[]
            seqsN=[]
            phred_scores=[]
            with open(os.path.join(dir,fastq_file), "r") as handle:

                for r,record in enumerate(SeqIO.parse(handle, "fastq")): # Parse reads
                    reads=r+1
                    record_id,record_seq,record_seqN,record_scores,alls,avgs,trims,masks = trim_filter(record,qall,qavg,qtrim,qmask,alls,avgs,trims,masks) # Trim & filter
                    if record_id is not None: # Save id, sequence, sequence masked, & quality scores
                        ids.append(record_id) 
                        seqs.append(record_seq)
                        seqsN.append(record_seqN)
                        phred_scores.append(record_scores)
        
        fastqs[fastq_file[:-len(suf)]]=pd.DataFrame({'id':ids, # Add dataframe to dictionary 
                                                     'seq':seqs,
                                                     'seqN':seqsN,
                                                     'phred_scores':phred_scores})
        print(f'{fastq_file[:-len(suf)]}:\t{reads} reads\t=>\t{len(fastqs[fastq_file[:-len(suf)]])} reads (alls = {alls} & avgs = {avgs});\t{trims} trimmed reads;\t{masks} masked reads')
        if save==True: out = pd.concat([out,
                                        pd.DataFrame({'file': [fastq_file[:-len(suf)]],
                                                      'reads': [reads],
                                                      'reads_filtered': [len(fastqs[fastq_file[:-len(suf)]])],
                                                      'reads_dropped_all': [alls],
                                                      'reads_dropped_avg': [avgs],
                                                      'reads_trimmed': [trims],
                                                      'reads_masked': [masks]})])
    
    if save==True: io.save(dir='.',file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_get_fastqs.csv',obj=out)
    return fastqs

# Determine editing outcomes methods
def region(fastqs: dict, pt='', flank5='', flank3='', save=True, masks=False):
    ''' 
    region(): gets DNA and AA sequence for records within flanks
    
    Parameters:
    fastqs (dict): dictionary from get_fastqs
    pt (str, optional 1): path to PrimeDesign input file (Required unless flank5 and flank3 are provided)
    flank5 (str, optional 2): top strand flanking sequence 5' (Required unless pt is provided)
    flank3 (str, optional 2): top strand flanking sequence 3' (Required unless pt is provided)
    save (bool, optional): save reads statistics file to local directory (Default: True)
    masks (bool, optional): include masked sequence and translsation (Default: False)
    
    Dependencies: pandas & Bio.Seq.Seq
    '''
    # Obtain flank5 and flank3 from pt or check that flank5 and flank3 have been provided
    if pt!='': (flank5,wt,flank3) = parse_input(pt)
    elif (flank5=='')&(flank3==''):
        raise ValueError('pt or flank3 and flank5 must be provided.')
    
    # Check flank lengths
    if (len(flank5)<9)|(len(flank3)<9): print('Warning: flank5 or flank3 less than 9.')

    # Remove fastq records that do not have flanks
    fastqs_1=dict()
    for file,fastq in fastqs.items():
        missing_flank = []
        for i,seq in enumerate(fastq['seq']):
            if (seq.find(flank5)==-1)|(seq.find(flank3)==-1): 
                missing_flank.append(i)

        fastqs_1[file] = fastq.drop(missing_flank).reset_index(drop=True)
     
    # Obtain nucleotide and AA sequences within flanks; remove fastq records with phred scores within flanks
    if save == True: out = pd.DataFrame()
    for file,fastq in fastqs_1.items():
        nuc=[]
        prot=[]
        if masks==True:
            nucN=[]
            protN=[]
        
        for i,seq in enumerate(fastq['seq']):
            nuc.append(seq[seq.find(flank5)+len(flank5):seq.find(flank3)])
            prot.append(Seq.translate(seq[seq.find(flank5)+len(flank5):seq.find(flank3)]))
            if masks==True:
                nucN.append(fastq.iloc[i]['seqN'][seq.find(flank5)+len(flank5):seq.find(flank3)])
                protN.append(Seq.translate(fastq.iloc[i]['seqN'][seq.find(flank5)+len(flank5):seq.find(flank3)]))
        
        fastqs_1[file]['nuc']=nuc
        fastqs_1[file]['prot']=prot
        if masks==True:
            fastqs_1[file]['nucN']=nuc
            fastqs_1[file]['protN']=protN
        
        print(f'{file}:\t{len(fastqs[file])} reads\t=>\t{len(fastqs_1[file])} reads')
        if save==True: out = pd.concat([out,
                                        pd.DataFrame({'file': [file],
                                                      'reads_filtered': [len(fastqs[file])],
                                                      'reads_w_flanks': [len(fastqs_1[file])]})])
    
    if save==True: io.save(dir='.',file=f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_region.csv',obj=out)
    
    return fastqs_1

def genotype(fastqs: dict, res: int, pt='', wt='', masks=False, keepX=False):
    ''' 
    genotype(): assign genotypes to sequence records
    
    Parameters:
    fastqs (dict): dictionary from filter_fastqs
    res (int): first AA number
    pt (str, optional 1): path to PrimeDesign input file (Required unless wt is provided)
    wt (str, optional 2): expected wildtype nucleotide sequence (in frame AA; required unless pt is provided)
    masks (bool, optional): include masked sequence and translsation (Default: False)
    keepX (bool, optional): keep unknown translation (i.e., X) due to sequencing error (Default: False) 
    
    Dependencies: pandas & Bio.Seq.Seq

    Note: Need to add single indels eventually
    '''
    # Obtain wt from pt or check that wt have been provided
    if pt!='': (flank5,wt,flank3) = parse_input(pt)
    elif wt=='':
        raise ValueError('pt or wt must be provided.')

    for file,fastq in fastqs.items():
        edit=[]
        if masks==True: editN=[]
        for i in range(len(fastq['prot'])):
            if len(wt)!=len(fastq.iloc[i]['nuc']): # Add single indels here
                edit.append('Indel')
                if masks==True: editN.append('Indel')
            elif Seq.translate(Seq(wt))==fastq.iloc[i]['prot']: 
                edit.append('WT')
                if masks==True: editN.append('WT')
            else:
                e = []
                if masks==True: eN = []
                
                for j, (a, b) in enumerate(zip(Seq.translate(Seq(wt)), fastq.iloc[i]['prot'])): # Find edits from sequence
                    if a != b: e.append(a+str(j+res)+b)
                if len(e)>1: edit.append(", ".join(e))
                elif len(e)==1: edit.append(e[0])
                else: edit.append('Unknown Edit')

                if masks==True:    
                    for j, (a, b) in enumerate(zip(Seq.translate(Seq(wt)), fastq.iloc[i]['protN'])): # Find edits from masked sequence
                        if (a != b)&(str(b)!='X')&(keepX==False): eN.append(a+str(j+res)+b)
                        elif (a != b)&(keepX==True): eN.append(a+str(j+res)+b)
                    if len(eN)>1: editN.append(", ".join(eN))
                    elif len(eN)==1: editN.append(eN[0])
                    else: editN.append('Masked Edit')

        fastqs[file]['edit']=edit
        if masks==True: fastqs[file]['editN']=editN
        print(f'{file}:\t{len(fastqs[file])} reads')
    
    return fastqs

def outcomes(fastqs: dict, edit='edit'):
    ''' 
    outcomes(): returns edit count & fraction per sample (tidy format)

    Parameters 
    fastqs (dict): dictionary from genotype
    edit (str, optional): edit column name
    
    Dependencies: pandas
    '''
    df = pd.DataFrame()
    for file,fastq in fastqs.items():
        temp=pd.DataFrame({'sample':[file]*len(fastq[edit].value_counts()),
                           edit:list(fastq[edit].value_counts().keys()),
                           'count':fastq[edit].value_counts(),
                           'fraction':fastq[edit].value_counts()/len(fastq[edit])})
        df=pd.concat([df,temp]).reset_index(drop=True)
    return df

# Supporting methods for plots
''' aa_props: dictionary of AA properties with citations (Generated by ChatGPT)
    
    Sources:
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

def edit_1(df: pd.DataFrame,col='edit'):
    ''' 
    edit_1(): split edit column to before, after, and amino acid number
    
    Parameters:
    df (dataframe): fastq outcomes dataframe
    col (str, optional): edit column name
    
    Dependencies: pandas
    '''
    df_1 = df[(df[col].str.contains(',')==False)&(df[col]!='WT')&(df[col]!='Indel')] # Isolate single AA changes
    df_1['before']=df_1[col].str[0] # Split edit information
    df_1['after']=df_1[col].str[-1]
    df_1['number']=df_1[col].str[1:-1].astype(int)
    return df_1.reset_index(drop=True)

def dms_tidy(df: pd.DataFrame, cond: str, wt:str, res: int):
    ''' 
    dms_tidy(): loads DMS grid data in tidy format split by condition
    
    Parameters:
    df (dataframe): fastq outcomes dataframe
    cond (str): Condition column name for splicing fastq outcomes
    wt (str): Expected wildtype nucleotide sequence (in frame AA)
    res (int): First AA number
    
    Dependencies: Bio.Seq.Seq, pandas, numpy, tidy, edit_1(), & aa_props
    '''
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
            vals['before']=[wt_prot[num-res]]*len(list(aa_props.keys()))
            vals['number']=[num]*len(list(aa_props.keys()))
            vals['after']=list(aa_props.keys())
            vals['edit']=[vals['before'][i]+str(num)+vals['after'][i] for i in range(len(vals['after']))]

            # Fill in counts and fractions for amino acid changes, WT, and none
            counts=[]
            fractions=[]
            num_mut = df_cond[df_cond['number']==num]
            for a in vals['after']:
                if a == wt_prot[num-res]: # Wild type
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

def vol_tidy(df: pd.DataFrame, cond: str, wt:str, res: int, psuedocount=1):
    ''' 
    vol_tidy(): loads DMS grid data in tidy format split by condition for volcano plot; uses psuedocounts
    
    Parameters:
    df (dataframe): fastq outcomes dataframe
    cond (str): Condition column name for splicing fastq outcomes
    wt (str): Expected wildtype nucleotide sequence (in frame AA)
    res (int): First AA number
    psuedocount (int, optional): psuedocount to avoid log(0) (Default: 1)
    
    Dependencies: Bio.Seq.Seq, pandas, numpy, tidy, edit_1(), & aa_props
    '''
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
            vals['before']=[wt_prot[num-res]]*len(list(aa_props.keys()))
            vals['number']=[num]*len(list(aa_props.keys()))
            vals['after']=list(aa_props.keys())
            vals['edit']=[vals['before'][i]+str(num)+vals['after'][i] for i in range(len(vals['after']))]

            # Fill in counts (+ psuedocount) for amino acid changes, WT, and none
            counts=[]
            num_mut = df_cond[df_cond['number']==num]
            for a in vals['after']:
                if a == wt_prot[num-res]: counts.append(wt_fastq['count'].to_list()[0]+psuedocount) # Wild type
                elif a in num_mut['after'].to_list(): counts.extend(num_mut[num_mut['after']==a]['count'].to_list()[0]+psuedocount) # Amino acid change present
                else: counts.append(psuedocount) # Amino acid change absent
            vals['count']=counts

            df_cond_DMS = pd.concat([df_cond_DMS,pd.DataFrame(vals)]).reset_index(drop=True) # Append residue DMS data
        
        df_cond_DMS['number']=df_cond_DMS['number'].astype(int) # Set number as type int
        df_cond_DMS['count']=df_cond_DMS['count'].astype(float) # Set count as type float for plotting

        dc2[key]=df_cond_DMS # Append condition DMS data
    
    dc3=dict() # Compare counts for all conditions
    for num in wt_nums: # Iterate through WT protein sequence
        
        for key,df_cond in dc2.items():

            vals=dict() # Create dictionary with all amino acid changes for a given residue
                
            # Add metadata that is the same for all genotypes
            meta = [x for x in df_cond.columns if x not in ['edit','count','fraction','before','after','number']]
            for m in meta: 
                vals[m]=[wt_fastq[m].to_list()[0]]*len(list(aa_props.keys()))

    return dc2

def subscript(df: pd.DataFrame,tick='before',tick_sub='number'):
    ''' 
    subscript(): returns dataframe with subscripts to tick labels
    
    Parameters:
    df (dataframe): dataframe
    tick (str, optional): new tick label column name
    tick_sub (str, optional): previous numeric tick label that will become a subscript

    Dependencies: pandas
    '''
    ticks = []
    labels = []
    for (t,ts) in set(zip(df[tick],df[tick_sub])):
        ticks.append(ts)
        labels.append('$\\mathrm{'+t+'_{'+str(ts)+'}}$')
    return pd.DataFrame({'tick':ticks,'label':labels}).sort_values(by='tick').reset_index(drop=True)

# Plotting mthods
def scat(typ: str,df: pd.DataFrame,x: str,y: str,cols=None,cols_ord=None,stys=None,cutoff=0.01,cols_exclude=None,
         file=None,dir=None,palette_or_cmap='colorblind',edgecol='black',
         figsize=(10,6),title='',title_size=18,title_weight='bold',
         x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,100),x_ticks_rot=0,xticks=[],
         y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,100),y_ticks_rot=0,yticks=[],
         legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0),
         **kwargs):
    ''' 
    scat(): creates scatter plot related graphs.

    Parameters:
    typ (str): plot type (scat, line, line_scat)
    df (dataframe): pandas dataframe
    x (str): x-axis column name
    y (str): y-axis column name
    cols (str, optional): color column name
    cols_ord (list, optional): color column values order
    stys (str, optional): styles column name
    cols_exclude (list, optional): color column values exclude
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn color palette or matplotlib color map
    edgecol (str, optional): point edge color
    figsize (tuple, optional): figure size
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_scale (str, optional): x-axis scale linear, log, etc.
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_ticks_rot (int, optional): x-axis ticks rotation
    xticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    yticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    
    Dependencies: os, matplotlib, seaborn, & plot
    '''
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

def cat(typ: str,df: pd.DataFrame,x: str,y: str,errorbar=None,cols=None,cols_ord=None,cutoff=0.01,cols_exclude=None,
        file=None,dir=None,palette_or_cmap='colorblind',edgecol='black',lw=1,
        figsize=(10,6),title='',title_size=18,title_weight='bold',
        x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,1),x_ticks_rot=0,xticks=[],
        y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,1),y_ticks_rot=0,yticks=[],
        legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0), 
        **kwargs):
    ''' 
    cat: creates category dependent graphs.
    
    Parameters:
    typ (str): plot type (bar, box, violin, swarm, strip, point, count, bar_swarm, box_swarm, violin_swarm)
    df (dataframe): pandas dataframe
    x (str, optional): x-axis column name
    y (str, optional): y-axis column name
    cols (str, optional): color column name
    cols_ord (list, optional): color column values order
    cols_exclude (list, optional): color column values exclude
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn color palette or matplotlib color map
    edgecol (str, optional): point edge color
    lw (int, optional): line width
    errorbar (str, optional): error bar type (sd)
    errwid (int, optional): error bar line width
    errcap (int, optional): error bar cap line width
    figsize (tuple, optional): figure size
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_scale (str, optional): x-axis scale linear, log, etc.
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_ticks_rot (int, optional): x-axis ticks rotation
    xticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    yticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    
    Dependencies: os, matplotlib, seaborn, & plot
    '''
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

def stack(df: pd.DataFrame,x='sample',y='fraction',cols='edit',cutoff=0.01,cols_ord=[],x_ord=[],
          file=None,dir=None,cmap='Set2',
          title='Editing Outcomes',title_size=18,title_weight='bold',
          figsize=(10,6),x_axis='',x_axis_size=12,x_axis_weight='bold',x_ticks_rot=45,x_ticks_ha='right',
          y_axis='',y_axis_size=12,y_axis_weight='bold',y_ticks_rot=0,
          legend_title='',legend_title_size=12,legend_size=12,
          legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_ncol=1,**kwargs):
    ''' 
    stack(): creates stacked bar plot

    Parameters:
    df (dataframe): pandas dataframe
    x (str, optional): x-axis column name
    y (str, optional): y-axis column name
    cols (str, optional): color column name
    cutoff (float, optional): y-axis values needs be greater than (e.g. 1%)
    cols_ord (list, optional): color column values order
    cols_exclude (list, optional): color column values exclude
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    cmap (str, optional): matplotlib color map
    errcap (int, optional): error bar cap line width
    figsize (tuple, optional): figure size
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_ha (str, optional): x-axis ticks horizontal alignment
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_ticks_rot (int, optional): y-axis ticks rotation
    yticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    
    Dependencies: re, os, pandas, numpy, matplotlib.pyplot & plot
    '''
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

def dms_grid(dc: dict,x='number',y='after',vals='count',
             file=None,dir=None,edgecol='black',lw=1,annot=False,cmap="bone_r",
             title='',title_size=12,title_weight='bold',
             x_axis='',x_axis_size=12,x_axis_weight='bold',x_ticks_rot=45,
             y_axis='',y_axis_size=12,y_axis_weight='bold',y_ticks_rot=0,
             **kwargs):
    ''' 
    dms_grid(): creates amino acide by residue number heatmaps from tidy-formatted DMS data
    
    Parameters:
    dc (dict): Tidy-formatted DMS dictionary
    x (str, optional): x-axis column name (AA residues number column)
    y (str, optional): y-axis column name (AA change column)
    vals (str, optional): values column name
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    edgecol (str, optional): point edge color
    lw (int, optional): line width
    annot (bool, optional): annotate values
    cmap (str, optional): matplotlib color map
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_ticks_rot (int, optional): x-axis ticks rotation
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_ticks_rot (int, optional): y-axis ticks rotation
    
    Dependencies: matplotlib, seaborn, pandas, & aa_props
    '''
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
        plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight', format=f'{file.split(".")[-1]}')
    plt.show()

def vol(df: pd.DataFrame,x: str,y: str, stys=None,size=None,
        file=None,dir=None,palette_or_cmap='colorblind',edgecol='black',
        figsize=(10,6),title='',title_size=18,title_weight='bold',
        x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,0),x_ticks_rot=0,xticks=[],
        y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,0),y_ticks_rot=0,yticks=[],
        legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0),legend_ncol=1,
        **kwargs):
    ''' 
    vol(): creates volcano plot
    
    Work in progress...

    Parameters:
    df (dataframe): pandas dataframe
    x (str): x-axis column name
    y (str): y-axis column name
    cols (str, optional): color column name
    stys (str, optional): styles column name
    size (str, optional): size column name
    cols_exclude (list, optional): color column values exclude
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn color palette or matplotlib color map
    edgecol (str, optional): point edge color
    figsize (tuple, optional): figure size
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_scale (str, optional): x-axis scale linear, log, etc.
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_ticks_rot (int, optional): x-axis ticks rotation
    xticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    yticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    
    Dependencies: os, matplotlib, seaborn, pandas, & edit_1()
    '''
    # Strings with subscripts
    log2 = 'log\u2082'
    log10 = 'log\u2081\u2080'
    
    # Log transform data
    df[f'{log2}({x})'] = [-np.log10(xval)/np.log10(2) for xval in df[x]]
    df[f'-{log10}({y})'] = [-np.log10(yval) for yval in df[y]]
    
    # Organize data by significance
    signif = []
    for (log2FC,log10P) in zip(df[f'{log2}({x})'],df[f'-{log10}({y})']):
        if (np.abs(log2FC)>1)&(log10P>-np.log10(0.05)): signif.append(f'{log2}FC & p-value')
        elif (np.abs(log2FC)<=1)&(log10P>-np.log10(0.05)): signif.append('p-value')
        elif (np.abs(log2FC)>1)&(log10P<=-np.log10(0.05)): signif.append(f'{log2}FC')
        else: signif.append('NS')
    df['Significance']=signif
    signif_order = ['NS',f'{log2}FC','p-value',f'{log2}FC & p-value']

    # Organize data by conservation (changed from)
    basic = ['R','K', 'H']
    acidic = ['D','E']
    polar = ['S', 'T', 'N', 'Q', 'Y', 'C']
    nonpolar = ['A','V','L','I','M','F','W','P','G']
    change = []

    df = edit_1(df)
    for (before,after) in zip(df['before'],df['after']):
        if (before in basic)&(after not in basic): change.append('basic')
        elif (before in acidic)&(after not in acidic): change.append('acidic')
        elif (before in polar)&(after not in polar): change.append('polar')
        elif (before in nonpolar)&(after not in nonpolar): change.append('nonpolar')
        else: change.append('Conserved')
    df['Change'] = change

    sty_order = ['Conserved','Basic','Acidic','Polar','Nonpolar']
    mark_order = ['D','^','v','<','>']

    # Organize data by abundance
    ''' size & sizes: 
        size: fold change from values
        sizes=(20, 200): Specifies the minimum and maximum marker sizes.
    '''
    sizes=(1,100)

    # Set dimensions
    if x_axis_dims==(0,0): x_axis_dims=(min(df[f'{log2}({x})']),max(df[f'{log2}({x})']))
    if y_axis_dims==(0,0): y_axis_dims=(0,max(df[f'-{log10}({y})']))

    # Generate figure
    fig, ax = plt.subplots(figsize=(5,5))
    
    # with significance boundraries
    plt.vlines(x=-1, ymin=y_axis_dims[0], ymax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
    plt.vlines(x=1, ymin=y_axis_dims[0], ymax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
    plt.hlines(y=-np.log10(0.05), xmin=x_axis_dims[0], xmax=x_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
    
    # with data
    sns.scatterplot(data=df, x=f'{log2}({x})', y=f'-{log10}({y})', 
                    hue='Significance', hue_order=signif_order, 
                    edgecolor=edgecol, palette='YlOrRd',
                    style_order=sty_order,markers=mark_order,
                    sizes=sizes,
                    ax=ax, **kwargs)
    
    # with labels
    df_signif = df[df['Significance']==f'{log2}FC & p-value']
    adjust_text([plt.text(x=df_signif.iloc[i][f'{log2}({x})'], 
                          y=df_signif.iloc[i][f'-{log10}({y})'],
                          s=edit) for i,edit in enumerate(df_signif['edit'])])

    
    # Move legend to the right of the graph
    ax.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
              bbox_to_anchor=legend_bbox_to_anchor,loc=legend_loc,ncol=legend_ncol)