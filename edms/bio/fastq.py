### fastq.py ###
# Author: Marc Zepeda
# Date: 2024-08-05

# Import packages
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Align import PairwiseAligner
import gzip
import os
import re
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from collections import Counter
from pathlib import Path
from scipy.stats import ttest_ind

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
    io.mkdir(out_dir) # Ensure the output directory exists

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

def unzip_fastqs(in_dir: str, out_dir: str):
    ''' 
    unzip_fastqs(): Unzip gzipped fastqs and write to a new directory.

    Parameters:
    in_dir (str): directory with compresesd fastq files
    out_dir (str): new directory with uncompressed fastq files
    
    Dependencies: gzip & os
    '''
    io.mkdir(out_dir) # Ensure the output directory exists

    for in_file in os.listdir(in_dir): # Find all .fastq.gz & .fastq files in the input directory
        print(f"Processing {in_file}...")
        if in_file.endswith(".fastq.gz"):
            with gzip.open(os.path.join(in_dir,in_file), 'rt') as handle:
                with open(os.path.join(out_dir,in_file.split('.fastq.gz')[0]+'.fastq'), 'wt') as out:
                    for line in handle:
                        out.write(line)

def comb_fastqs(in_dir: str, out_dir: str, out_file: str):
    ''' 
    comb_fastqs(): Combines one or more (un)compressed fastqs files into a single (un)compressed fastq file.

    Parameters:
    in_dir (str): directory with fastq files
    out_dir (str): new directory with combined fastq file
    out_file (str): Name of output fastq file (Needs .fastq or .fastq.gz suffix)
    
    Dependencies: gzip & os
    '''
    io.mkdir(out_dir) # Ensure the output directory exists

    if out_file.endswith(".fastq.gz"):
        with gzip.open(os.path.join(out_dir,out_file), 'wt') as out:
            for in_file in os.listdir(in_dir): # Find all .fastq.gz & .fastq files in the input directory
                print(f"Processing {in_file}...")
                if in_file.endswith(".fastq.gz"):
                    with gzip.open(os.path.join(in_dir,in_file), 'rt') as handle:
                        for line in handle:
                            out.write(line)
                
                elif in_file.endswith(".fastq"):
                    with open(os.path.join(in_dir,in_file), 'r') as handle:
                        for line in handle:
                            out.write(line)
    
    elif out_file.endswith(".fastq"):
        with open(os.path.join(out_dir,out_file), 'wt') as out:
            for in_file in os.listdir(in_dir): # Find all .fastq.gz & .fastq files in the input directory
                print(f"Processing {in_file}...")
                if in_file.endswith(".fastq.gz"):
                    with gzip.open(os.path.join(in_dir,in_file), 'rt') as handle:
                        for line in handle:
                            out.write(line)
                
                elif in_file.endswith(".fastq"):
                    with open(os.path.join(in_dir,in_file), 'r') as f:
                        for line in handle:
                            out.write(line)

    else: print('out_file needs .fastq or .fastq.gz suffix')

# Quantify epegRNA abundance methods
def count_spacers(sample_sheet: str, annotated_lib: str, fastq_dir: str, KEY_INTERVAL=(10,80), 
                  KEY_FLANK5='CGAAACACC', KEY_FLANK3='GTTTAAGA', spacer_col='Spacer_sequence', 
                  dont_trim_G=False, out_dir='', out_file='library_count_spacers.csv', save=True, 
                  return_df=True, save_files=True, plot_out_type='pdf'):
    
    """[Summary]
    Given a set of epegRNA sequences and a FASTQ file, count the reads in the
    FASTQ, assign the reads to epegRNAs, and export the counts to a csv file `out_counts`. All
    sgRNA sequences not found in the reference file (non-perfect matches) are
    written to a separate csv file `out_nc`. 
    count_reads works by reading through a .fastq file and searching for the guide in between
    the subsequences for KEY_FLANK5 and KEY_FLANK3 within the KEY_INTERVAL. These counts are then recorded.

    Parameters
    ----------
    sample_sheet : str or path
        REQUIRED COLS: 'fastq_file', 'counts_file', 'noncounts_file', 'stats_file'
        a sheet with information on sequence id, 
        fastq_R1_file, fastq_R2_file (string or path to the FASTQ file to be processed), 
        out_counts (string or path for the output csv file with perfect sgRNA matches ex: 'counts.csv'),
        out_nc (string or path for the output csv file with non-perfect sgRNA matches ex: 'noncounts.csv'), 
        out_stats (string or path for the output txt file with the read counting statistics ex: 'stats.txt'), 
        condition names, and condition categories
    annotated_lib : str or path
        String or path to the reference file. annotated_lib must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    fastq_dir : str or path, defaults to ''
        String or path to the directory where all fastq files are found. 
    KEY_INTERVAL : tuple, default (10,80)
        Tuple of (KEY_START, KEY_END) that defines the KEY_REGION. Denotes the
        substring within the read to search for the KEY.
    KEY_FLANK5 : str, default 'CGAAACACC'
        Sequence that is expected upstream of the spacer sequence. The default
        is the end of the hU6 promoter.
    KEY_FLANK3 : str, default 'GTTTGAGA'
        Sequence that is expected downstream of the spacer sequence. The
        default is the start of the sgRNA scaffold sequence.
    spacer_col : str, default 'Spacer_sequence'
        Spacer sequence column name in annotated library
    dont_trim_G : bool, default False
        Whether to trim the first G from 21-nt sgRNA sequences to make them 20-nt.
    out_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_file : str or path, defaults to 'counts_library.csv'
        Name of output dataframe with guides and counts. 
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save : bool, default True
        Whether or not to save the resulting dataframe
    save_files : bool, default True
        Whether or not to save individual counts, noncounts, and stats files
    plot_out_type : str, optional, defaults to 'pdf'
        file type of figure output
    
    Dependencies: os,pandas,Path,gzip,matplolib,Counter,numpy,io
    """
    io.mkdir(out_dir) # Make output directory if it does not exist

    sample_filepath = Path(sample_sheet)
    sample_df = pd.read_csv(sample_filepath)
    for colname in ['fastq_file', 'counts_file', 'noncounts_file', 'stats_file', 'condition']: 
        if colname not in sample_df.columns.tolist():
            raise Exception(f"annotated_lib is missing column: {colname}")
    samples = [list(a) for a in zip(sample_df.fastq_file, sample_df.counts_file, 
                                    sample_df.noncounts_file, sample_df.stats_file, 
                                    sample_df.condition)]

    # STEP 1A: OPEN INPUT FILES FOR PROCESSING, CHECK FOR REQUIRED FORMATTING
    # look for 'sgRNA_seq' column, raise Exception if missing
    annotated_lib = Path(annotated_lib)
    df_ref = pd.read_csv(annotated_lib, header=0) # explicit header = first row
    if spacer_col not in df_ref.columns.tolist():
        raise Exception(f'annotated_lib is missing column: {spacer_col}')
    df_ref[spacer_col] = df_ref[spacer_col].str.upper() 
    path = Path.cwd()

    for fastq, counts, nc, stats, cond in samples: 
        # fastq file of reads and paths to all output files, imported from sample_sheet
        in_fastq = Path(fastq_dir) / fastq
        out_counts, out_nc, out_stats = Path(out_dir) / counts, Path(out_dir) / nc, Path(out_dir) / stats        
        # try opening input FASTQ, raise Exception if not possible
        handle = gzip.open(in_fastq, 'rt') if str(in_fastq).endswith('.gz') else open(in_fastq, 'r')

        # STEP 1B: SET UP VARIABLES FOR SCRIPT
        # make dictionary to hold sgRNA counts - sgRNA_seq, count as k,v
        dict_p = {sgRNA:0 for sgRNA in df_ref[spacer_col]}
        list_np = [] # placeholder list for non-perfect matches
        # reads count of: total, perfect match, non perfect match, no key found, not 20bps
        num_reads, num_p_matches, num_np_matches, num_nokey, num_badlength = 0, 0, 0, 0, 0
        KEY_START, KEY_END = KEY_INTERVAL[0], KEY_INTERVAL[1] # set the key interval

        # STEP 2: PROCESS FASTQ FILE READS AND ADD COUNTS TO DICT
        while True: # contains the seq and Qscore etc.
            read = handle.readline()
            if not read: # end of file
                break
            elif read.startswith("@"): # if line is a read header
                read = handle.readline() # next line is the actual sequence
            else:
                continue
            num_reads += 1
            read_sequence = str.upper(str(read))
            key_region = read_sequence[KEY_START:KEY_END]
            key_index = key_region.find(KEY_FLANK5)
            key_rev_index = key_region.rfind(KEY_FLANK3)
            if key_index < 0 or key_rev_index <= key_index: # if keys not found
                num_nokey += 1
                continue
            start_index = key_index + KEY_START + len(KEY_FLANK5)
            end_index = key_rev_index + KEY_START
            guide = read_sequence[start_index:end_index]
            if not dont_trim_G:
                if guide.startswith('G') and len(guide) == 21:
                    guide = guide[1:]
            if len(guide) != 20:
                num_badlength += 1
                continue
            if guide in dict_p:
                dict_p[guide] += 1
                num_p_matches += 1
            else:
                num_np_matches += 1
                list_np.append(guide)
        handle.close()

        # STEP 3: SORT DICTIONARIES AND GENERATE OUTPUT FILES
        # sort perf matches (A-Z) with guides, counts as k,v and output to csv
        df_perfects = pd.DataFrame(data=dict_p.items(), columns=[spacer_col, cond])
        if save_files:
            df_perfects.sort_values(by=cond, ascending=False, inplace=True)
            df_perfects.to_csv(out_counts, index=False)

            weights = np.ones_like(df_perfects[cond]) / len(df_perfects[cond])
            # histogram for df_ref[cond] column
            plt.hist(df_perfects[cond], weights=weights, bins=(len(df_perfects[cond])//5)+1)
            outpath = path / out_dir
            out = stats.split('.')[0] + '_histogram.' + plot_out_type
            plt.title(f"Distributions of sgRNA in {fastq}")
            plt.xlabel('Count of sgRNA')
            plt.ylabel('Proportion of sgRNA')
            plt.savefig(outpath / out, format=plot_out_type)
            plt.clf()

        # add matching counts to dataframe
        df_ref = pd.merge(df_ref, df_perfects, on=spacer_col, how='outer')
        df_ref[cond] = df_ref[cond].fillna(0)

        # now sort non-perfect matches by frequency and output to csv
        dict_np = Counter(list_np) # use Counter to tally up np matches
        nc_name = nc.split("/")[-1]
        df_ncmatches = pd.DataFrame(data=dict_np.items(), columns=[spacer_col, nc_name])
        if save_files:
            df_ncmatches.sort_values(by=nc_name, ascending=False, inplace=True)
            df_ncmatches.to_csv(out_nc, index=False)
        # calculate the read coverage (reads processed / sgRNAs in library)
        num_guides = df_ref[spacer_col].shape[0]
    
        # STEP 4: CALCULATE STATS AND GENERATE STAT OUTPUT FILE
        # percentage of guides that matched perfectly
        pct_p_match = round(num_p_matches/float(num_p_matches + num_np_matches) * 100, 1)
        # percentage of undetected guides (no read counts)
        vals_p = np.fromiter(dict_p.values(), dtype=int)
        guides_no_reads = np.count_nonzero(vals_p==0)
        pct_no_reads = round(guides_no_reads/float(len(dict_p.values())) * 100, 1)
        # skew ratio of top 10% to bottom 10% of guide counts
        top_10 = np.percentile(list(dict_p.values()), 90)
        bottom_10 = np.percentile(list(dict_p.values()), 10)
        if top_10 != 0 and bottom_10 != 0:
            skew_ratio = top_10/bottom_10
        else:
            skew_ratio = 'Not enough perfect matches to determine skew ratio'
        # calculate the read coverage (reads processed / sgRNAs in library)
        coverage = round(num_reads / num_guides, 1)
        # calculate the number of unmapped reads (num_nokey / total_reads)
        pct_unmapped = round((num_nokey / num_reads) * 100, 2)
        # write analysis statistics to statfile
        if save_files:
            with open(out_stats, 'w') as statfile:
                statfile.write('Number of reads processed: ' + str(num_reads) + '\n')
                statfile.write('Number of reads where key was not found: ' + str(num_nokey) + '\n')
                statfile.write('Number of reads where length was not 20bp: ' + str(num_badlength) + '\n')
                statfile.write('Number of perfect guide matches: ' + str(num_p_matches) + '\n')
                statfile.write('Number of nonperfect guide matches: ' + str(num_np_matches) + '\n')
                statfile.write('Number of undetected guides: ' + str(guides_no_reads) + '\n')
                statfile.write('Percentage of unmapped reads (key not found): ' + str(pct_unmapped) + '\n') #
                statfile.write('Percentage of guides that matched perfectly: ' + str(pct_p_match) + '\n') #
                statfile.write('Percentage of undetected guides: ' + str(pct_no_reads) + '\n') #
                statfile.write('Skew ratio of top 10% to bottom 10%: ' + str(skew_ratio) + '\n') #
                statfile.write('Read coverage: ' + str(coverage))
                statfile.close()
                print(str(in_fastq), 'processed')

    plt.close()
    # export files and return dataframes if necessary
    if save: 
        outpath = path / out_dir
        Path.mkdir(outpath, exist_ok=True)
        df_ref.to_csv(outpath / out_file, index=False)
        print('count_reads outputed to', str(outpath / out_file))
    print('Count reads completed')
    if return_df:
        return df_ref

def count_spacers_pbs(sample_sheet: str, annotated_lib: str, fastq_R1_dir: str, fastq_R2_dir: str,
                             KEY_INTERVAL=(10,80), KEY_FLANK5='CGAAACACC', KEY_FLANK3='GTTTAAGA', 
                             spacer_col='Spacer_sequence', dont_trim_G=False, KEY_FLANK5_REV='AACCGCG', 
                             pbs_col='PBS_sequence', pbs_len_col='PBS_length', linker_col='Linker_sequence', out_dir='', 
                             out_file='library_count_spacers_pbs.csv', save=True, return_df=True, 
                             save_files=True, plot_out_type='pdf'):
    
    """[Summary]
    Given a set of epegRNA sequences and a FASTQ file, count the reads in the
    FASTQ, assign the reads to epegRNAs, and export the counts to a csv file `out_counts`. All
    sgRNA sequences not found in the reference file (non-perfect matches) are
    written to a separate csv file `out_nc`. 
    count_reads works by reading through a .fastq file and searching for the guide in between
    the subsequences for KEY_FLANK5 and KEY_FLANK3 within the KEY_INTERVAL. These counts are then recorded.

    Parameters
    ----------
    sample_sheet : str or path
        REQUIRED COLS: 'fastq_R1_file', 'fastq_R2_file', 'counts_file', 'noncounts_file', 'stats_file'
        a sheet with information on sequence id, 
        fastq_R1_file, fastq_R2_file (string or path to the FASTQ file to be processed), 
        out_counts (string or path for the output csv file with perfect sgRNA matches ex: 'counts.csv'),
        out_nc (string or path for the output csv file with non-perfect sgRNA matches ex: 'noncounts.csv'), 
        out_stats (string or path for the output txt file with the read counting statistics ex: 'stats.txt'), 
        condition names, and condition categories
    annotated_lib : str or path
        String or path to the reference file. annotated_lib must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    fastq_R1_dir : str or path, defaults to ''
        String or path to the directory where all fastq R1 files are found. 
    fastq_R2_dir : str or path, defaults to ''
        String or path to the directory where all fastq R2 files are found. 
    KEY_INTERVAL : tuple, default (10,80)
        Tuple of (KEY_START, KEY_END) that defines the KEY_REGION. Denotes the
        substring within the read to search for the KEY.
    KEY_FLANK5 : str, default 'CGAAACACC'
        Sequence that is expected upstream of the spacer sequence. The default
        is the end of the hU6 promoter.
    KEY_FLANK3 : str, default 'GTTTGAGA'
        Sequence that is expected downstream of the spacer sequence. The
        default is the start of the sgRNA scaffold sequence.
    spacer_col : str, default 'Spacer_sequence'
        Spacer sequence column name in annotated library
    dont_trim_G : bool, default False
        Whether to trim the first G from 21-nt sgRNA sequences to make them 20-nt.
    KEY_FLANK5_REV : str, default 'AACCGCG'
        Sequence that is expected downstream of the extebsuib sequence. The default
        is the start of the tevoPreQ1 motif.
    pbs_col : str, default 'PBS_sequence'
        PBS sequence column name in annotated library
    pbs_len_col : str, default 'PBS_length'
        PBS length column name in annotated library
    linker_col : str, default 'Linker_sequence'
        Linker sequence column name in annotated library
    out_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_file : str or path, defaults to 'counts_library.csv'
        Name of output dataframe with guides and counts. 
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save : bool, default True
        Whether or not to save the resulting dataframe
    save_files : bool, default True
        Whether or not to save individual counts, noncounts, and stats files
    plot_out_type : str, optional, defaults to 'pdf'
        file type of figure output
    
    Dependencies: os,pandas,Path,gzip,matplolib,Counter,numpy,io
    """
    io.mkdir(out_dir) # Make output directory if it does not exist

    sample_filepath = Path(sample_sheet)
    sample_df = pd.read_csv(sample_filepath)
    for colname in ['fastq_R1_file', 'fastq_R2_file', 'counts_file', 'noncounts_file', 'stats_file', 'condition']: 
        if colname not in sample_df.columns.tolist():
            raise Exception(f"annotated_lib is missing column: {colname}")
    samples = [list(a) for a in zip(sample_df.fastq_R1_file, sample_df.fastq_R2_file,
                                    sample_df.counts_file, sample_df.noncounts_file, 
                                    sample_df.stats_file, sample_df.condition)]

    # STEP 1A: OPEN INPUT FILES FOR PROCESSING, CHECK FOR REQUIRED FORMATTING
    # look for 'Spacer_squence', 'PBS_squence', & 'Linker_squence' columns, raise Exception if missing
    annotated_lib = Path(annotated_lib)
    df_ref = pd.read_csv(annotated_lib, header=0) # explicit header = first row
    if spacer_col not in df_ref.columns.tolist():
        raise Exception(f'annotated_lib is missing column: {spacer_col}')
    if pbs_col not in df_ref.columns.tolist():
        raise Exception(f'annotated_lib is missing column: {pbs_col}')
    if linker_col not in df_ref.columns.tolist():
        raise Exception(f'annotated_lib is missing column: {linker_col}')
    df_ref[spacer_col] = df_ref[spacer_col].str.upper()
    spacers = set(df_ref[spacer_col])
    df_ref[pbs_col] = df_ref[pbs_col].str.upper() 
    df_ref[linker_col] = df_ref[linker_col].str.upper() 
    path = Path.cwd()

    for fastq_R1, fastq_R2, counts, nc, stats, cond in samples: 
        # fastq file of reads and paths to all output files, imported from sample_sheet
        in_fastq_R1 = Path(fastq_R1_dir) / fastq_R1
        in_fastq_R2 = Path(fastq_R2_dir) / fastq_R2
        out_counts, out_nc, out_stats = Path(out_dir) / counts, Path(out_dir) / nc, Path(out_dir) / stats        
        # try opening input FASTQ, raise Exception if not possible
        reader1 = gzip.open(in_fastq_R1, 'rt') if str(in_fastq_R1).endswith('.gz') else open(in_fastq_R1, 'r')
        reader2 = gzip.open(in_fastq_R2, 'rt') if str(in_fastq_R2).endswith('.gz') else open(in_fastq_R2, 'r')

        # STEP 1B: SET UP VARIABLES FOR SCRIPT
        # make dictionary to hold epegRNA counts - spacer_pbs_linker, count as k,v
        dict_p = {f'{spacer}_{pbs}':0 for (spacer,pbs) in zip(df_ref[spacer_col],df_ref[pbs_col])}
        list_np = [] # placeholder list for non-perfect matches
        # reads count of: total, perfect match, non perfect match, no key found, not 20bps
        num_reads_R1, num_p_matches_R1, num_np_matches_R1, num_nokey_R1, num_badlength_R1 = 0, 0, 0, 0, 0
        num_reads_R2, num_p_matches_R2, num_np_matches_R2, num_nokey_R2, num_badlength_R2 = 0, 0, 0, 0, 0
        KEY_START, KEY_END = KEY_INTERVAL[0], KEY_INTERVAL[1] # set the key interval

        # STEP 2A: PROCESS FASTQ R1 FILE READS AND ADD COUNTS TO DICT
        while True: # contains the seq and Qscore etc.
            read1 = reader1.readline()
            read2 = reader2.readline()
            if not read1: # end of file
                break
            elif read1.startswith("@"): # if line is a read header
                read1 = reader1.readline() # next line is the actual sequence
                read2 = reader2.readline() # next line is the actual sequence
            else:
                continue
            num_reads_R1 += 1
            read_sequence = str.upper(str(read1))
            key_region = read_sequence[KEY_START:KEY_END]
            key_index = key_region.find(KEY_FLANK5)
            key_rev_index = key_region.rfind(KEY_FLANK3)
            if key_index < 0 or key_rev_index <= key_index: # if keys not found
                num_nokey_R1 += 1
                continue
            start_index = key_index + KEY_START + len(KEY_FLANK5)
            end_index = key_rev_index + KEY_START
            guide = read_sequence[start_index:end_index]
            if not dont_trim_G:
                if guide.startswith('G') and len(guide) == 21:
                    guide = guide[1:]
            if len(guide) != 20:
                num_badlength_R1 += 1
                continue
            if guide in spacers:
                num_p_matches_R1 += 1
            else:
                num_np_matches_R1 += 1
                list_np.append(f'{guide}_spacer')
                continue
            
            # STEP 2B: PROCESS FASTQ R2 FILE READS AND ADD COUNTS TO DICT
            num_reads_R2 += 1
            read2_sequence = str.upper(str(read2))
            read2_key_index = read2_sequence.find(KEY_FLANK5_REV)
            if read2_key_index < 0: # if keys not found
                num_nokey_R2 += 1
                continue
            start_index = read2_key_index + len(KEY_FLANK5_REV)
            read2_region = read2_sequence[start_index:]
            df_ref_guide = df_ref[df_ref[spacer_col]==guide]
            df_ref_guide[f'RC_{pbs_col}'] = [str(Seq(pbs).reverse_complement()) for pbs in df_ref_guide[pbs_col]]
            df_ref_guide.sort_values(by=pbs_len_col,ascending=False,inplace=True)
            for i,rc_pbs in enumerate(df_ref_guide[f'RC_{pbs_col}']):
                if rc_pbs in read2_region:
                    dict_p[f'{guide}_{df_ref_guide.iloc[i][pbs_col]}'] += 1
                    num_p_matches_R2 += 1
                    break
                elif i+1>=len(df_ref_guide[f'RC_{pbs_col}']):
                    num_np_matches_R2 += 1
                    list_np.append(f'{guide}_pbs')
                    break
        reader1.close()
        reader2.close()
            

        # STEP 3: SORT DICTIONARIES AND GENERATE OUTPUT FILES
        # sort perf matches (A-Z) with epegRNAs, counts as k,v and output to csv
        df_perfects = pd.DataFrame(data=dict_p.items(), columns=[f'{spacer_col}_{pbs_col}', cond])
        spacers_perfects = []
        pbs_perfects = []
        for seqs in df_perfects[f'{spacer_col}_{pbs_col}']:
            spacers_perfects.append(seqs.split('_')[0])
            pbs_perfects.append(seqs.split('_')[1])
        df_perfects[spacer_col] = spacers_perfects
        df_perfects[pbs_col] = pbs_perfects
        df_perfects = df_perfects[[spacer_col,pbs_col,cond]]
        if save_files:
            df_perfects.sort_values(by=cond, ascending=False, inplace=True)
            df_perfects.to_csv(out_counts, index=False)

            weights = np.ones_like(df_perfects[cond]) / len(df_perfects[cond])
            # histogram for df_ref[cond] column
            plt.hist(df_perfects[cond], weights=weights, bins=(len(df_perfects[cond])//5)+1)
            outpath = path / out_dir
            out = stats.split('.')[0] + '_histogram.' + plot_out_type
            plt.title(f"Distributions of epegRNAs in\n{fastq_R1} & {fastq_R2}")
            plt.xlabel('Count of epegRNA')
            plt.ylabel('Proportion of epegRNA')
            plt.savefig(outpath / out, format=plot_out_type)
            plt.clf()

        # add matching counts to dataframe
        df_ref = pd.merge(df_ref, df_perfects, on=[spacer_col,pbs_col], how='outer')
        df_ref[cond] = df_ref[cond].fillna(0)

        # now sort non-perfect matches by frequency and output to csv
        dict_np = Counter(list_np) # use Counter to tally up np matches
        nc_name = nc.split("/")[-1]
        df_ncmatches = pd.DataFrame(data=dict_np.items(), columns=[f'{spacer_col}_mismatch', nc_name])
        spacers_ncs = []
        mismatches = []
        for spacer_mismatch in df_ncmatches[f'{spacer_col}_mismatch']:
            spacers_ncs.append(spacer_mismatch.split('_')[0])
            mismatches.append('_'.join(spacer_mismatch.split('_')[1:]))
        df_ncmatches[spacer_col] = spacers_ncs
        df_ncmatches['mismatch'] = mismatches
        if save_files:
            df_ncmatches.sort_values(by=nc_name, ascending=False, inplace=True)
            df_ncmatches.to_csv(out_nc, index=False)
        # calculate the read coverage (reads processed / sgRNAs in library)
        num_guides = df_ref[spacer_col].shape[0]
    
        # STEP 4: CALCULATE STATS AND GENERATE STAT OUTPUT FILE
        # percentage of guides that matched perfectly
        pct_p_match_R1 = round(num_p_matches_R1/float(num_p_matches_R1 + num_np_matches_R1) * 100, 1)
        pct_p_match_R2 = round(num_p_matches_R2/float(num_p_matches_R2 + num_np_matches_R2) * 100, 1)
        # percentage of undetected guides (no read counts)
        vals_p = np.fromiter(dict_p.values(), dtype=int)
        guides_no_reads = np.count_nonzero(vals_p==0)
        pct_no_reads = round(guides_no_reads/float(len(dict_p.values())) * 100, 1)
        # skew ratio of top 10% to bottom 10% of guide counts
        top_10 = np.percentile(list(dict_p.values()), 90)
        bottom_10 = np.percentile(list(dict_p.values()), 10)
        if top_10 != 0 and bottom_10 != 0:
            skew_ratio = top_10/bottom_10
        else:
            skew_ratio = 'Not enough perfect matches to determine skew ratio'
        # calculate the read coverage (reads processed / sgRNAs in library)
        coverage = round(num_reads_R1 / num_guides, 1)
        # calculate the number of unmapped reads (num_nokey / total_reads)
        pct_unmapped_R1 = round((num_nokey_R1 / num_reads_R1) * 100, 2)
        pct_unmapped_R2 = round((num_nokey_R2 / num_reads_R2) * 100, 2)
        # write analysis statistics to statfile
        if save_files:
            with open(out_stats, 'w') as statfile:
                statfile.write('Number of R1 reads processed: ' + str(num_reads_R1) + '\n')
                statfile.write('Number of R2 reads processed: ' + str(num_reads_R2) + '\n')
                statfile.write('Number of R1 reads where key was not found: ' + str(num_nokey_R1) + '\n')
                statfile.write('Number of R2 reads where key was not found: ' + str(num_nokey_R2) + '\n')
                statfile.write('Number of R1 reads where length was not 20bp: ' + str(num_badlength_R1) + '\n')
                statfile.write('Number of R1 perfect epegRNA matches: ' + str(num_p_matches_R1) + '\n')
                statfile.write('Number of R2 perfect epegRNA matches: ' + str(num_p_matches_R2) + '\n')
                statfile.write('Number of R1 nonperfect epegRNA matches: ' + str(num_np_matches_R1) + '\n')
                statfile.write('Number of R2 nonperfect epegRNA matches: ' + str(num_np_matches_R2) + '\n')
                statfile.write('Number of undetected epegRNAs: ' + str(guides_no_reads) + '\n')
                statfile.write('Percentage of R1 unmapped reads (key not found): ' + str(pct_unmapped_R1) + '\n') #
                statfile.write('Percentage of R2 unmapped reads (key not found): ' + str(pct_unmapped_R2) + '\n') #
                statfile.write('Percentage of epegRNA that matched perfectly in R1: ' + str(pct_p_match_R1) + '\n') #
                statfile.write('Percentage of epegRNA that matched perfectly in R2: ' + str(pct_p_match_R2) + '\n') 
                statfile.write('Percentage of undetected epegRNAs: ' + str(pct_no_reads) + '\n') #
                statfile.write('Skew ratio of top 10% to bottom 10%: ' + str(skew_ratio) + '\n') #
                statfile.write('Read coverage: ' + str(coverage))
                statfile.close()
                print(str(in_fastq_R1), 'processed')
                print(str(in_fastq_R2), 'processed')

    plt.close()
    # export files and return dataframes if necessary
    if save: 
        outpath = path / out_dir
        Path.mkdir(outpath, exist_ok=True)
        df_ref.to_csv(outpath / out_file, index=False)
        print('count_spacers_pbs outputed to', str(outpath / out_file))
    print('Count reads completed')
    if return_df:
        return df_ref

def count_spacers_pbs_linkers(sample_sheet: str, annotated_lib: str, fastq_R1_dir: str, fastq_R2_dir: str,
                             KEY_INTERVAL=(10,80), KEY_FLANK5='CGAAACACC', KEY_FLANK3='GTTTAAGA', 
                             spacer_col='Spacer_sequence', dont_trim_G=False, KEY_FLANK5_REV='AACCGCG', 
                             pbs_col='PBS_sequence', pbs_len_col='PBS_length', linker_col='Linker_sequence', out_dir='', 
                             out_file='library_count_spacers_pbs_linkers.csv', save=True, return_df=True, 
                             save_files=True, plot_out_type='pdf'):
    
    """[Summary]
    Given a set of epegRNA sequences and a FASTQ file, count the reads in the
    FASTQ, assign the reads to epegRNAs, and export the counts to a csv file `out_counts`. All
    sgRNA sequences not found in the reference file (non-perfect matches) are
    written to a separate csv file `out_nc`. 
    count_reads works by reading through a .fastq file and searching for the guide in between
    the subsequences for KEY_FLANK5 and KEY_FLANK3 within the KEY_INTERVAL. These counts are then recorded.

    Parameters
    ----------
    sample_sheet : str or path
        REQUIRED COLS: 'fastq_R1_file', 'fastq_R2_file', 'counts_file', 'noncounts_file', 'stats_file'
        a sheet with information on sequence id, 
        fastq_R1_file, fastq_R2_file (string or path to the FASTQ file to be processed), 
        out_counts (string or path for the output csv file with perfect sgRNA matches ex: 'counts.csv'),
        out_nc (string or path for the output csv file with non-perfect sgRNA matches ex: 'noncounts.csv'), 
        out_stats (string or path for the output txt file with the read counting statistics ex: 'stats.txt'), 
        condition names, and condition categories
    annotated_lib : str or path
        String or path to the reference file. annotated_lib must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    fastq_R1_dir : str or path, defaults to ''
        String or path to the directory where all fastq R1 files are found. 
    fastq_R2_dir : str or path, defaults to ''
        String or path to the directory where all fastq R2 files are found. 
    KEY_INTERVAL : tuple, default (10,80)
        Tuple of (KEY_START, KEY_END) that defines the KEY_REGION. Denotes the
        substring within the read to search for the KEY.
    KEY_FLANK5 : str, default 'CGAAACACC'
        Sequence that is expected upstream of the spacer sequence. The default
        is the end of the hU6 promoter.
    KEY_FLANK3 : str, default 'GTTTGAGA'
        Sequence that is expected downstream of the spacer sequence. The
        default is the start of the sgRNA scaffold sequence.
    spacer_col : str, default 'Spacer_sequence'
        Spacer sequence column name in annotated library
    dont_trim_G : bool, default False
        Whether to trim the first G from 21-nt sgRNA sequences to make them 20-nt.
    KEY_FLANK5_REV : str, default 'AACCGCG'
        Sequence that is expected downstream of the extebsuib sequence. The default
        is the start of the tevoPreQ1 motif.
    pbs_col : str, default 'PBS_sequence'
        PBS sequence column name in annotated library
    pbs_len_col : str, default 'PBS_length'
        PBS length column name in annotated library
    linker_col : str, default 'Linker_sequence'
        Linker sequence column name in annotated library
    out_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_file : str or path, defaults to 'counts_library.csv'
        Name of output dataframe with guides and counts. 
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save : bool, default True
        Whether or not to save the resulting dataframe
    save_files : bool, default True
        Whether or not to save individual counts, noncounts, and stats files
    plot_out_type : str, optional, defaults to 'pdf'
        file type of figure output

    Dependencies: os,pandas,Path,gzip,matplolib,Counter,numpy,io
    """
    io.mkdir(out_dir) # Make output directory if it does not exist

    sample_filepath = Path(sample_sheet)
    sample_df = pd.read_csv(sample_filepath)
    for colname in ['fastq_R1_file', 'fastq_R2_file', 'counts_file', 'noncounts_file', 'stats_file', 'condition']: 
        if colname not in sample_df.columns.tolist():
            raise Exception(f"annotated_lib is missing column: {colname}")
    samples = [list(a) for a in zip(sample_df.fastq_R1_file, sample_df.fastq_R2_file,
                                    sample_df.counts_file, sample_df.noncounts_file, 
                                    sample_df.stats_file, sample_df.condition)]

    # STEP 1A: OPEN INPUT FILES FOR PROCESSING, CHECK FOR REQUIRED FORMATTING
    # look for 'Spacer_squence', 'PBS_squence', & 'Linker_squence' columns, raise Exception if missing
    annotated_lib = Path(annotated_lib)
    df_ref = pd.read_csv(annotated_lib, header=0) # explicit header = first row
    if spacer_col not in df_ref.columns.tolist():
        raise Exception(f'annotated_lib is missing column: {spacer_col}')
    if pbs_col not in df_ref.columns.tolist():
        raise Exception(f'annotated_lib is missing column: {pbs_col}')
    if linker_col not in df_ref.columns.tolist():
        raise Exception(f'annotated_lib is missing column: {linker_col}')
    df_ref[spacer_col] = df_ref[spacer_col].str.upper()
    spacers = set(df_ref[spacer_col])
    df_ref[pbs_col] = df_ref[pbs_col].str.upper() 
    df_ref[linker_col] = df_ref[linker_col].str.upper() 
    path = Path.cwd()

    for fastq_R1, fastq_R2, counts, nc, stats, cond in samples: 
        # fastq file of reads and paths to all output files, imported from sample_sheet
        in_fastq_R1 = Path(fastq_R1_dir) / fastq_R1
        in_fastq_R2 = Path(fastq_R2_dir) / fastq_R2
        out_counts, out_nc, out_stats = Path(out_dir) / counts, Path(out_dir) / nc, Path(out_dir) / stats        
        # try opening input FASTQ, raise Exception if not possible
        reader1 = gzip.open(in_fastq_R1, 'rt') if str(in_fastq_R1).endswith('.gz') else open(in_fastq_R1, 'r')
        reader2 = gzip.open(in_fastq_R2, 'rt') if str(in_fastq_R2).endswith('.gz') else open(in_fastq_R2, 'r')
        
        # STEP 1B: SET UP VARIABLES FOR SCRIPT
        # make dictionary to hold epegRNA counts - spacer_pbs_linker, count as k,v
        dict_p = {f'{spacer}_{pbs}_{linker}':0 for (spacer,pbs,linker) in zip(df_ref[spacer_col],df_ref[pbs_col],df_ref[linker_col])}
        list_np = [] # placeholder list for non-perfect matches
        # reads count of: total, perfect match, non perfect match, no key found, not 20bps
        num_reads_R1, num_p_matches_R1, num_np_matches_R1, num_nokey_R1, num_badlength_R1 = 0, 0, 0, 0, 0
        num_reads_R2, num_p_matches_R2, num_np_matches_R2, num_nokey_R2, num_badlength_R2 = 0, 0, 0, 0, 0
        KEY_START, KEY_END = KEY_INTERVAL[0], KEY_INTERVAL[1] # set the key interval

        # STEP 2A: PROCESS FASTQ R1 FILE READS AND ADD COUNTS TO DICT
        while True: # contains the seq and Qscore etc.
            read1 = reader1.readline()
            read2 = reader2.readline()
            if not read1: # end of file
                break
            elif read1.startswith("@"): # if line is a read header
                read1 = reader1.readline() # next line is the actual sequence
                read2 = reader2.readline()
            else:
                continue
            num_reads_R1 += 1
            read_sequence = str.upper(str(read1))
            key_region = read_sequence[KEY_START:KEY_END]
            key_index = key_region.find(KEY_FLANK5)
            key_rev_index = key_region.rfind(KEY_FLANK3)
            if key_index < 0 or key_rev_index <= key_index: # if keys not found
                num_nokey_R1 += 1
                continue
            start_index = key_index + KEY_START + len(KEY_FLANK5)
            end_index = key_rev_index + KEY_START
            guide = read_sequence[start_index:end_index]
            if not dont_trim_G:
                if guide.startswith('G') and len(guide) == 21:
                    guide = guide[1:]
            if len(guide) != 20:
                num_badlength_R1 += 1
                continue
            if guide in spacers:
                num_p_matches_R1 += 1
            else:
                num_np_matches_R1 += 1
                list_np.append(f'{guide}_spacer')
                continue
            
            # STEP 2B: PROCESS FASTQ R2 FILE READS AND ADD COUNTS TO DICT
            num_reads_R2 += 1
            read2_sequence = str.upper(str(read2))
            read2_key_index = read2_sequence.find(KEY_FLANK5_REV)
            if read2_key_index < 0: # if keys not found
                num_nokey_R2 += 1
                continue
            start_index = read2_key_index + len(KEY_FLANK5_REV)
            read2_region = read2_sequence[start_index:]
            df_ref_guide = df_ref[df_ref[spacer_col]==guide]
            df_ref_guide['RC_PBS_Linker'] = [str(Seq('').join([pbs,linker]).reverse_complement()) 
                                for (pbs,linker) in zip(df_ref_guide[pbs_col],df_ref_guide[linker_col])]
            df_ref_guide.sort_values(by=pbs_len_col,ascending=False,inplace=True)
            for i,rc_pbs_linker in enumerate(df_ref_guide['RC_PBS_Linker']):
                if rc_pbs_linker in read2_region:
                    dict_p[f'{guide}_{df_ref_guide.iloc[i][pbs_col]}_{df_ref_guide.iloc[i][linker_col]}'] += 1
                    num_p_matches_R2 += 1
                    break
                elif i+1>=len(df_ref_guide['RC_PBS_Linker']):
                    num_np_matches_R2 += 1
                    list_np.append(f'{guide}_pbs_linker')
                    break
        reader1.close()
        reader2.close()  

        # STEP 3: SORT DICTIONARIES AND GENERATE OUTPUT FILES
        # sort perf matches (A-Z) with epegRNAs, counts as k,v and output to csv
        df_perfects = pd.DataFrame(data=dict_p.items(), columns=[f'{spacer_col}_{pbs_col}_{linker_col}', cond])
        spacers_perfects = []
        pbs_perfects = []
        linker_perfects = []
        for seqs in df_perfects[f'{spacer_col}_{pbs_col}_{linker_col}']:
            spacers_perfects.append(seqs.split('_')[0])
            pbs_perfects.append(seqs.split('_')[1])
            linker_perfects.append(seqs.split('_')[2])
        df_perfects[spacer_col] = spacers_perfects
        df_perfects[pbs_col] = pbs_perfects
        df_perfects[linker_col] = linker_perfects
        df_perfects = df_perfects[[spacer_col,pbs_col,linker_col,cond]]
        if save_files:
            df_perfects.sort_values(by=cond, ascending=False, inplace=True)
            df_perfects.to_csv(out_counts, index=False)

            weights = np.ones_like(df_perfects[cond]) / len(df_perfects[cond])
            # histogram for df_ref[cond] column
            plt.hist(df_perfects[cond], weights=weights, bins=(len(df_perfects[cond])//5)+1)
            outpath = path / out_dir
            out = stats.split('.')[0] + '_histogram.' + plot_out_type
            plt.title(f"Distributions of epegRNAs in\n{fastq_R1} & {fastq_R2}")
            plt.xlabel('Count of epegRNA')
            plt.ylabel('Proportion of epegRNA')
            plt.savefig(outpath / out, format=plot_out_type)
            plt.clf()
        
        # add matching counts to dataframe
        df_ref = pd.merge(df_ref, df_perfects, on=[spacer_col,pbs_col,linker_col], how='outer')
        df_ref[cond] = df_ref[cond].fillna(0)

        # now sort non-perfect matches by frequency and output to csv
        dict_np = Counter(list_np) # use Counter to tally up np matches
        nc_name = nc.split("/")[-1]
        df_ncmatches = pd.DataFrame(data=dict_np.items(), columns=[f'{spacer_col}_mismatch', nc_name])
        spacers_ncs = []
        mismatches = []
        for spacer_mismatch in df_ncmatches[f'{spacer_col}_mismatch']:
            spacers_ncs.append(spacer_mismatch.split('_')[0])
            mismatches.append('_'.join(spacer_mismatch.split('_')[1:]))
        df_ncmatches[spacer_col] = spacers_ncs
        df_ncmatches['mismatch'] = mismatches
        if save_files:
            df_ncmatches.sort_values(by=nc_name, ascending=False, inplace=True)
            df_ncmatches.to_csv(out_nc, index=False)
        # calculate the read coverage (reads processed / sgRNAs in library)
        num_guides = df_ref[spacer_col].shape[0]
    
        # STEP 4: CALCULATE STATS AND GENERATE STAT OUTPUT FILE
        # percentage of guides that matched perfectly
        pct_p_match_R1 = round(num_p_matches_R1/float(num_p_matches_R1 + num_np_matches_R1) * 100, 1)
        pct_p_match_R2 = round(num_p_matches_R2/float(num_p_matches_R2 + num_np_matches_R2) * 100, 1)
        # percentage of undetected guides (no read counts)
        vals_p = np.fromiter(dict_p.values(), dtype=int)
        guides_no_reads = np.count_nonzero(vals_p==0)
        pct_no_reads = round(guides_no_reads/float(len(dict_p.values())) * 100, 1)
        # skew ratio of top 10% to bottom 10% of guide counts
        top_10 = np.percentile(list(dict_p.values()), 90)
        bottom_10 = np.percentile(list(dict_p.values()), 10)
        if top_10 != 0 and bottom_10 != 0:
            skew_ratio = top_10/bottom_10
        else:
            skew_ratio = 'Not enough perfect matches to determine skew ratio'
        # calculate the read coverage (reads processed / sgRNAs in library)
        coverage = round(num_reads_R1 / num_guides, 1)
        # calculate the number of unmapped reads (num_nokey / total_reads)
        pct_unmapped_R1 = round((num_nokey_R1 / num_reads_R1) * 100, 2)
        pct_unmapped_R2 = round((num_nokey_R2 / num_reads_R2) * 100, 2)
        # write analysis statistics to statfile
        if save_files:
            with open(out_stats, 'w') as statfile:
                statfile.write('Number of R1 reads processed: ' + str(num_reads_R1) + '\n')
                statfile.write('Number of R2 reads processed: ' + str(num_reads_R2) + '\n')
                statfile.write('Number of R1 reads where key was not found: ' + str(num_nokey_R1) + '\n')
                statfile.write('Number of R2 reads where key was not found: ' + str(num_nokey_R2) + '\n')
                statfile.write('Number of R1 reads where length was not 20bp: ' + str(num_badlength_R1) + '\n')
                statfile.write('Number of R1 perfect epegRNA matches: ' + str(num_p_matches_R1) + '\n')
                statfile.write('Number of R2 perfect epegRNA matches: ' + str(num_p_matches_R2) + '\n')
                statfile.write('Number of R1 nonperfect epegRNA matches: ' + str(num_np_matches_R1) + '\n')
                statfile.write('Number of R2 nonperfect epegRNA matches: ' + str(num_np_matches_R2) + '\n')
                statfile.write('Number of undetected epegRNAs: ' + str(guides_no_reads) + '\n')
                statfile.write('Percentage of R1 unmapped reads (key not found): ' + str(pct_unmapped_R1) + '\n') #
                statfile.write('Percentage of R2 unmapped reads (key not found): ' + str(pct_unmapped_R2) + '\n') #
                statfile.write('Percentage of epegRNA that matched perfectly in R1: ' + str(pct_p_match_R1) + '\n') #
                statfile.write('Percentage of epegRNA that matched perfectly in R2: ' + str(pct_p_match_R2) + '\n') 
                statfile.write('Percentage of undetected epegRNAs: ' + str(pct_no_reads) + '\n') #
                statfile.write('Skew ratio of top 10% to bottom 10%: ' + str(skew_ratio) + '\n') #
                statfile.write('Read coverage: ' + str(coverage))
                statfile.close()
                print(str(in_fastq_R1), 'processed')
                print(str(in_fastq_R2), 'processed')

    plt.close()
    # export files and return dataframes if necessary
    if save: 
        outpath = path / out_dir
        Path.mkdir(outpath, exist_ok=True)
        df_ref.to_csv(outpath / out_file, index=False)
        print('count_spacers_pbs_linkers outputed to', str(outpath / out_file))
    print('Count reads completed')
    if return_df:
        return df_ref

def count_region(fastq_dir: str, end_i: int, start_i:int=0, 
                 out_dir:str=None, out_file:str=None):
    ''' 
    count_region(): returns dataframe with sequence region abundance for every fastq file in a directory

    Parameters:
    fastq_dir (str): path to fastq directory
    end_i (int): region end index (zero-indexed)
    start_i (int, optional): region start index (Default: 0, zero-indexed)
    out_dir (str, optional): path to save directory (Default: None)
    out_file (str, optional): save file name (Default: None)

    Dependencies: pandas,gzip,os,tidy
    '''
    sequences_dc = dict()
    for fastq_file in os.listdir(fastq_dir): # Find all .fastq.gz & .fastq files in the fastq directory
        print(f"Processing {fastq_file}...") # Keep track of sequence regions & reads
        regions = [] 
        reads = 0 

        if fastq_file.endswith(".fastq.gz"): # Compressed fastq
            with gzip.open(os.path.join(fastq_dir,fastq_file), 'rt') as handle:
                for r,record in enumerate(SeqIO.parse(handle, "fastq")): # Parse reads
                    reads=r+1
                    regions.append(''.join(record.seq[start_i:end_i])) # Obtain region
                 
        elif fastq_file.endswith(".fastq"): # Uncompressed fastq
            with open(os.path.join(fastq_dir,fastq_file), 'r') as handle:
                for r,record in enumerate(SeqIO.parse(handle, "fastq")): # Parse reads    
                    reads=r+1
                    regions.append(''.join(record.seq[start_i:end_i])) # Obtain region

        print(f'Completed {reads} reads') # Create dataframe with abundance of regions
        sequences_dc[fastq_file]=pd.Series(regions).value_counts().to_frame().reset_index().rename(columns={'index': 'region'})
        sequences_dc[fastq_file]['fraction']=sequences_dc[fastq_file]['count']/reads
    
    df = t.join(dc=sequences_dc,col='fastq_file') # Join fastq files into single dataframe
    df = df[['fastq_file','region','count','fraction']] # Change column order

    if out_dir is not None and out_file is not None: # Save dataframe (optional)
        io.mkdir(out_dir) # Make output directory if it does not exist
        io.save(dir=out_dir,file=out_file,obj=df)

    return df

def count_alignments(annotated_lib: str, align_col: str, id_col: str, fastq_dir: str, 
                     out_dir: str, fastq_suf='.fastq.gz', match_score=1, mismatch_score=-4,
                     align_max:int=None,plot_suf='.pdf', show=False, **plot_kwargs):
    ''' 
    count_alignments(): get fastq files from directory and store records in dataframes in a dictionary
    
    Parameters:
    annotated_lib (str): path to the annotated library reference file
    align_col (str): align column name in annotated library reference file
    id_col (str): id column name in annotated library reference file
    fastq_dir (str): directory with fastq files
    out_dir (str): directory for output files
    fastq_suf (str, optional): file suffix (.fastq.gz or .fastq)
    match_score (int, optional): match score for pairwise alignment (Default: 1)
    mismatch_score (int, optional): mismatch score for pairwise alignment (Default: -4)
    align_max (int, optional): max alignments per fastq file to save compute (Default: None)
    plot_suf (str, optional): plot type suffix with '.' (Default: '.pdf')
    show (bool, optional): show plots (Default: False)
    **plot_kwargs (optional): plot key word arguments

    Dependencies: Bio.SeqIO, gzip, os, pandas, Bio.Seq.Seq, Bio.PairwiseAligner, & trim_filter()
    '''
    # Intialize the aligner
    aligner = PairwiseAligner()
    aligner.mode = 'global'  # Use 'local' for local alignment
    aligner.match_score = match_score  # Score for a match
    aligner.mismatch_score = mismatch_score/2  # Penalty for a mismatch; applied to both strands
    aligner.open_gap_score = mismatch_score/2  # Penalty for opening a gap; applied to both strands
    aligner.extend_gap_score = mismatch_score/2  # Penalty for extending a gap; applied to both strands

    # Obtain reference file & check for alignment column
    annotatated_lib_name = annotated_lib.split('/')[-1].split('.')[0]
    df_ref = io.get(annotated_lib)
    if align_col not in df_ref.columns.tolist():
        raise Exception(f'{annotated_lib} is missing column: {align_col}') 
    if id_col not in df_ref.columns.tolist():
        raise Exception(f'{annotated_lib} is missing column: {id_col}')
    df_ref[align_col] = df_ref[align_col].str.upper() 

    # Obtain fastq files
    files = os.listdir(fastq_dir)
    fastq_files = [file for file in files if fastq_suf in file]
    
    # Make fastqs dictionary
    fastqs = dict()
    for fastq_file in fastq_files:
        
        # Get reads
        if fastq_suf=='.fastq.gz': # Compressed fastq files
            with gzip.open(os.path.join(fastq_dir,fastq_file), "rt") as handle:
                seqs=[record.seq for record in SeqIO.parse(handle, "fastq")] # Parse reads
        else: # Uncompressed fastq files
            with open(os.path.join(fastq_dir,fastq_file), "r") as handle:    
                seqs=[record.seq for record in SeqIO.parse(handle, "fastq")] # Parse reads
        fastq_name = fastq_file[:-len(fastq_suf)]
        print(f'{fastq_name}:\t{len(seqs)} reads')
        
        # Perform alignments
        print('Perform alignments')
        dc_alignments = {ref:0 for ref in df_ref[align_col]}
        dc_alignments_mismatch_num = {ref:0 for ref in df_ref[align_col]}
        dc_alignments_mismatch_pos = {ref:[] for ref in df_ref[align_col]}
        s=0
        for seq in seqs: # Iterate though sequences
            s+=1
            if align_max is not None:
                if s>align_max: break
            print(f'{s} out of {len(seqs)}')
            seq_alignments_scores = []
            seq_alignments_aligned = []
            for ref in df_ref[align_col]: # Iterate though reference sequences
                seq_alignment = aligner.align(ref, seq[0:len(ref)]) # trim ngs sequence to reference sequence & align
                seq_alignments_scores.append(seq_alignment[0].score) # Save highest alignment score
                seq_alignments_aligned.append(seq_alignment[0].aligned[0]) # Save alignment matches

            # Isolate maximum score alignment
            i = seq_alignments_scores.index(max(seq_alignments_scores))
            ref_i = df_ref.iloc[i][align_col]
            aligned_i = seq_alignments_aligned[i]
            dc_alignments[df_ref.iloc[i][align_col]] = dc_alignments[ref_i]+1

            # Find & quantify mismatches (Change zero-indexed to one-indexed)
            mismatch_pos = []
            if len(aligned_i) == 1: 
                (a1,b1) = aligned_i[0]
                if (a1==0)&(b1==len(ref_i)-1): mismatch_pos.extend([])
                elif a1==0: mismatch_pos.extend([k+1 for k in range(b1+1,len(ref_i))])
                elif b1==len(ref_i)-1: mismatch_pos.extend([k+1 for k in range(0,a1-1)])
                else: mismatch_pos.extend([j+1 for j in range(0,a1-1)] + [k+1 for k in range(b1+1,len(ref_i))])
            else:
                for j in range(len(aligned_i)-1):
                    (a1,b1) = aligned_i[j]
                    (a2,b2) = aligned_i[j+1]
                    if (j==0)&(a1!=0): mismatch_pos.extend([k+1 for k in range(0,a1-1)])
                    if (j==len(aligned_i)-2)&(b2!=len(ref_i)-1): mismatch_pos.extend([k+1 for k in range(b2+1,len(ref_i))])
                    mismatch_pos.extend([k+1 for k in range(b1+1,a2-1)])
            dc_alignments_mismatch_num[ref_i] = dc_alignments_mismatch_num[ref_i] + len(mismatch_pos)
            dc_alignments_mismatch_pos[ref_i] = dc_alignments_mismatch_pos[ref_i] + mismatch_pos

        # Calculate mismatch position fraction of alignments
        dc_alignments_mismatch_pos_fraction = dict()
        for (ref,mismatch_pos) in dc_alignments_mismatch_pos.items():
            dc_alignments_mismatch_pos_fraction[ref] = {pos:mismatch_pos.count(pos) for pos in range(1,len(ref)+1)}

        # Merge alignment dictionaries into a fastq dataframe
        print('Merge alignment dictionaries into a fastq dataframe')
        df_alignments = pd.DataFrame(dc_alignments.items(),columns=[align_col,'alignments'])
        df_alignments_mismatch_num = pd.DataFrame(dc_alignments_mismatch_num.items(),columns=[align_col,'mismatch_num'])
        df_alignments_mismatch_pos = pd.DataFrame(dc_alignments_mismatch_pos.items(),columns=[align_col,'mismatch_pos']) 
        df_fastq = pd.merge(left=df_ref,right=df_alignments,on=align_col)
        df_fastq = pd.merge(left=df_fastq,right=df_alignments_mismatch_num,on=align_col)
        df_fastq = pd.merge(left=df_fastq,right=df_alignments_mismatch_pos,on=align_col)
        
        # Calculate mismatch num & position per alignment
        print('Calculate mismatch num & position per alignment')
        mismatch_num_per_alignment_ls = []
        mismatch_pos_per_alignment_ls = []
        for (ref,mismatch_pos,mismatch_num,alignments) in t.zip_cols(df=df_fastq,cols=[align_col,'mismatch_pos','mismatch_num','alignments']):
            if alignments==0:
                mismatch_num_per_alignment_ls.append(0)
                mismatch_pos_per_alignment_ls.append({pos:0 for pos in range(1,len(ref)+1)})
            else:
                mismatch_num_per_alignment_ls.append(mismatch_num/alignments)
                mismatch_pos_per_alignment_ls.append({pos:mismatch_pos.count(pos)/alignments for pos in range(1,len(ref)+1)})
        df_fastq['mismatch_num_per_alignment'] = mismatch_num_per_alignment_ls
        df_fastq['mismatch_pos_per_alignment'] = mismatch_pos_per_alignment_ls
        
        # Save & append fastq dataframe to fastq dictionary
        print('Save & append fastq dataframe to fastq dictionary')
        io.save(dir=out_dir,file=f'alignment_{annotatated_lib_name}_{fastq_name}.csv',obj=df_fastq)
        fastqs[fastq_name]=df_fastq

        # Plot mismatch position per alignment
        print('Plot mismatch position per alignment')
        
        out_dir_fastq_name = os.path.join(out_dir,fastq_name)
        df_fastq_plot = pd.DataFrame()
        for align,id,mismatch_pos_per_alignment in t.zip_cols(df=df_fastq,cols=[align_col,id_col,'mismatch_pos_per_alignment']):
            df_fastq_plot_align = pd.DataFrame({align_col:[align]*len(mismatch_pos_per_alignment), # Obtain individual alignments
                                                id_col:[id]*len(mismatch_pos_per_alignment),
                                                'mismatch_pos':list(mismatch_pos_per_alignment.keys()),
                                                'mismatch_pos_per_alignment':list(mismatch_pos_per_alignment.values())})
            
            p.scat(typ='line',df=df_fastq_plot_align,x='mismatch_pos',y='mismatch_pos_per_alignment', # Plot mismatches for each alignment
                   title=f'{fastq_name} {id}',x_axis='Alignment Position',y_axis='Mismatches/Alignment',
                   dir=out_dir_fastq_name,file=f'{id.replace(".","_")}{plot_suf}',
                   show=show,**plot_kwargs)
            
            df_fastq_plot = pd.concat(objs=[df_fastq_plot,df_fastq_plot_align]).reset_index(drop=True) # Group alignment mismatches

        p.scat(typ='line',df=df_fastq_plot,x='mismatch_pos',y='mismatch_pos_per_alignment',cols=id_col, # Plot mismatches for each alignment
               title=f'{fastq_name}',x_axis='Alignment Position',y_axis='Mismatches/Alignment',
               dir=out_dir,file=f'{fastq_name}{plot_suf}',legend_ncol=4,
               show=show,**plot_kwargs)
    
    return fastqs

def plot_alignments(fastq_alignments: dict, align_col: str, id_col: str,
                     out_dir: str, plot_suf='.pdf', show=False, **plot_kwargs):
    ''' 
    plot_alignments(): plot fastq alignments dictionary output from count_alignments()
    
    Parameters:
    fastq_alignments (dict): fastq alignments dictionary output from count_alignments()
    align_col (str): align column name in annotated library reference file
    id_col (str): id column name in annotated library reference file
    fastq_dir (str): directory with fastq files
    out_dir (str): directory for output files
    plot_suf (str, optional): plot type suffix with '.' (Default: '.pdf')
    show (bool, optional): show plots (Default: False)
    **plot_kwargs (optional): plot key word arguments

    Dependencies: Bio.SeqIO, gzip, os, pandas, Bio.Seq.Seq, Bio.PairwiseAligner, & trim_filter()
    '''
    for fastq_name,df_fastq in fastq_alignments.items():
        
        # Plot mismatch position per alignment
        print('Plot mismatch position per alignment')
        
        out_dir_fastq_name = os.path.join(out_dir,fastq_name)
        df_fastq_plot = pd.DataFrame()
        for align,id,mismatch_pos_per_alignment in t.zip_cols(df=df_fastq,cols=[align_col,id_col,'mismatch_pos_per_alignment']):
            df_fastq_plot_align = pd.DataFrame({align_col:[align]*len(mismatch_pos_per_alignment), # Obtain individual alignments
                                                id_col:[id]*len(mismatch_pos_per_alignment),
                                                'mismatch_pos':list(mismatch_pos_per_alignment.keys()),
                                                'mismatch_pos_per_alignment':list(mismatch_pos_per_alignment.values())})
            
            p.scat(typ='line',df=df_fastq_plot_align,x='mismatch_pos',y='mismatch_pos_per_alignment', # Plot mismatches for each alignment
                   title=f'{fastq_name} {id}',x_axis='Alignment Position',y_axis='Mismatches/Alignment',
                   dir=out_dir_fastq_name,file=f'{id.replace(".","_")}{plot_suf}',
                   show=show,**plot_kwargs)
            
            df_fastq_plot = pd.concat(objs=[df_fastq_plot,df_fastq_plot_align]).reset_index(drop=True) # Group alignment mismatches

        p.scat(typ='line',df=df_fastq_plot,x='mismatch_pos',y='mismatch_pos_per_alignment',cols=id_col, # Plot mismatches for each alignment
               title=f'{fastq_name}',x_axis='Alignment Position',y_axis='Mismatches/Alignment',
               dir=out_dir,file=f'{fastq_name}{plot_suf}',
               show=show,**plot_kwargs)
# Quantify edit outcome methods
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
            
            if qtrim!=0: # Save compute time if trim is not desired
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
            if masks !=0: # Save compute time if mask is not desired
                for i, qual in enumerate(quality_scores):
                    if qual < qmask: bases[i] = 'N'
            sequenceN = Seq('').join(bases) # Update the sequence with the modified version
            if Seq('N') in sequenceN: masks += 1

            return record.id,sequence,sequenceN,quality_scores,alls,avgs,trims,masks
    
        else: return None,None,None,None,alls,avgs+1,trims,masks # Avg threshold not met
    else: return None,None,None,None,alls+1,avgs,trims,masks # All threshold not met

def get_fastqs(dir: str,suf='.fastq.gz',qall=10,qavg=30,qtrim=30,qmask=0,save=True):
    ''' 
    get_fastqs(): get fastq files from directory and store records in dataframes in a dictionary
    
    Parameters:
    dir (str): directory with fastq files
    suf (str, optional): file suffix (.fastq.gz or .fastq)
    qall (int, optional): phred quality score threshold for all bases for a read to not be discarded (Q = -log(err))
    qtrim (int, optional): phred quality score threshold for trimming reads on both ends (Q = -log(err))
    qavg (int, optional): average phred quality score threshold for a read to not be discarded (Q = -log(err))
    qmask (int, optional): phred quality score threshold for base to not be masked to N (Q = -log(err))
    save (bool, optional): save reads statistics file to local directory (Default: True)

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

        if suf=='.fastq.gz': # Compressed fastq files
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

    Parameters:
    fastqs (dict): dictionary from genotype
    edit (str, optional): edit column name (Default: edit)
    
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

def outcomes_desired(df: pd.DataFrame, desired_edits: list | str, sample_col='sample',
                     edit_col='edit', count_col='count',fraction_col='fraction'):
    ''' 
    outcomes_desired: groups desired edit count & fraction per sample (tidy format)

    Parameters:
    df (DataFrame): dataframe with edit count & fraction per sample (tidy format)
    desired_edits (list or str): list of desired edits (list of str) or desired edits column name (str)
    sample_col (str, optional): sample column name (Default: sample)
    edit_col (str, optional): edit column name (Default: edit)
    count_col (str, optional): count column name (Default: count)
    fraction_col (str, optional): fraction column name (Default: fraction)

    Dependencies: pandas
    '''
    if isinstance(desired_edits, list): desired_edits_col = None # Determine if desired edits is a list or str
    elif isinstance(desired_edits, str): desired_edits_col = desired_edits
    else: TypeError(f'desired_edits = {desired_edits} was not a list or str.')

    df_desired = pd.DataFrame()
    for sample in df[sample_col].value_counts().keys(): # Iterate through samples
        df_sample = df[df[sample_col]==sample].reset_index(drop=True)

        if desired_edits_col: 
            desired_edits = df_sample.iloc[0][desired_edits_col] # Get desired edits list for each sample if the column name was provided
            if isinstance(desired_edits, str): desired_edits = [desired_edits]
        
        i_desired = [] # Store desired edit & corresponding counts & fractions
        count_desired = []
        fraction_desired = []

        for i,(edit,count,fraction) in enumerate(t.zip_cols(df=df_sample,cols=[edit_col,count_col,fraction_col])):
            
            if ', ' in edit: # Search for desired edit within multiple edit outcomes
                edits = edit.split(', ')
                for edit in edits:
                    if edit in desired_edits:
                        i_desired.append(i)
                        count_desired.append(count)
                        fraction_desired.append(fraction)
                        break

            else: # Search for desired within single edit outcome
                if edit in desired_edits:
                    i_desired.append(i)
                    count_desired.append(count)
                    fraction_desired.append(fraction)

        df_sample = df_sample.drop(index=i_desired) # Remove desired edits & combine into 'Desired' edit
        other_cols = [col for col in df_sample.columns if col not in [edit_col,count_col,fraction_col]]
        df_sample_desired = df_sample.iloc[0][other_cols].to_frame().T.reset_index(drop=True)
        df_sample_desired[edit_col] = ['Desired']
        df_sample_desired[count_col] = [sum(count_desired)]
        df_sample_desired[fraction_col] = [sum(fraction_desired)]
        df_sample = pd.concat(objs=[df_sample,df_sample_desired]).reset_index(drop=True)
        df_desired = pd.concat(objs=[df_desired,df_sample]).reset_index(drop=True)

    return df_desired
        
# Supporting methods for DMS plots
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

def dms_cond(df: pd.DataFrame, cond: str, wt:str, res: int, sample='sample', edit='edit', psuedocount=0):
    ''' 
    dms_cond(): returns DMS grid data in tidy format grouped by condition
    
    Parameters:
    df (dataframe): fastq outcomes dataframe
    cond (str): Condition column name for grouping fastq outcomes dataframe
    wt (str): Expected wildtype nucleotide sequence (in frame AA)
    res (int): First AA number
    sample (str, optional): Sample column name for splicing fastq outcomes dataframe (Default: 'sample')
    edit (str, optional): Edit column name within fastq outcomes dataframe (Default: 'edit')
    psuedocount (int, optional): psuedocount to avoid log(0) & /0 (Default: 0)
    
    Dependencies: Bio.Seq.Seq, pandas, numpy, tidy, edit_1(), & aa_props
    '''
    wt_prot = Seq(wt).translate(table=1) # Obtain WT protein sequence
    wt_nums = np.arange(res,res+len(wt_prot))
    print('Isolate single aa change fastq outcomes')
    dc=t.split(edit_1(df),sample) # Isolate single aa change fastq outcomes and split by sample
    
    print('Fill with DMS grid data for each sample:')
    dc2=dict() # Fill with DMS grid data in tidy format split by sample
    for key_sample,df_sample in dc.items():
        print(key_sample)
        wt_fastq = df[(df['edit']=='WT')&(df[sample]==key_sample)] # Obtain WT fastq outcome
        df_sample_DMS=pd.DataFrame(columns=wt_fastq.columns) # Fill with DMS grid data in tidy format
        
        for num in wt_nums: # Iterate through WT protein sequence
            vals=dict() # Create dictionary with all amino acid changes for a given residue
            
            # Add metadata that is the same for all genotypes
            meta = [x for x in df_sample.columns if x not in [edit,'count','fraction','before','after','number']]
            for m in meta: 
                vals[m]=[wt_fastq[m].to_list()[0]]*len(list(aa_props.keys()))
            
            # Create all amino acid changes
            vals['before']=[wt_prot[num-res]]*len(list(aa_props.keys()))
            vals['number']=[num]*len(list(aa_props.keys()))
            vals['after']=list(aa_props.keys())
            vals[edit]=[vals['before'][i]+str(num)+vals['after'][i] for i in range(len(vals['after']))]

            # Fill in counts (+ psuedocount) for amino acid changes, WT, and none
            counts=[]
            num_mut = df_sample[df_sample['number']==num]
            for a in vals['after']:
                if a == wt_prot[num-res]: counts.append(wt_fastq['count'].to_list()[0]+psuedocount) # Wild type
                elif a in num_mut['after'].to_list(): counts.append(num_mut[num_mut['after']==a]['count'].to_list()[0]+psuedocount) # Amino acid change present
                else: counts.append(psuedocount) # Amino acid change absent
            vals['count']=counts
            sum_counts = sum(vals['count'])
            vals['fraction']=[count/sum_counts for count in vals['count']]

            df_sample_DMS = pd.concat([df_sample_DMS,pd.DataFrame(vals)]).reset_index(drop=True) # Append residue DMS data
        
        df_sample_DMS['number']=df_sample_DMS['number'].astype(int) # Set number as type int
        df_sample_DMS['count']=df_sample_DMS['count'].astype(int) # Set count as type int for plotting

        df_sample_DMS[sample] = [key_sample]*df_sample_DMS.shape[0]
        dc2[key_sample]=df_sample_DMS # Append sample DMS data

    print('Group samples by condition:')
    dc3=t.split(t.join(dc2,sample),cond) # Join samples back into 1 dataframe & split by condition
    df_cond_stat = pd.DataFrame()
    for key_cond,df_cond in dc3.items(): # Iterate through conditions
        print(key_cond)
        edit_ls = []
        fraction_avg_ls = []
        fraction_ls = []
        count_avg_ls = []
        before_ls = []
        after_ls = []
        number_ls = []
        for e in df_cond[edit]: # iterate through edits
            df_cond_edit = df_cond[df_cond[edit]==e]
            edit_ls.append(e)
            fraction_avg_ls.append(sum(df_cond_edit['fraction'])/len(df_cond_edit['fraction']))
            fraction_ls.append(df_cond_edit['fraction'].tolist())
            count_avg_ls.append(sum(df_cond_edit['count'])/len(df_cond_edit['count']))
            before_ls.append(df_cond_edit.iloc[0]['before'])
            after_ls.append(df_cond_edit.iloc[0]['after'])
            number_ls.append(df_cond_edit.iloc[0]['number'])
        df_cond_stat = pd.concat([df_cond_stat,
                                  pd.DataFrame({'edit':edit_ls,
                                                'before':before_ls,
                                                'after':after_ls,
                                                'number':number_ls,
                                                'fraction_ls':fraction_ls,
                                                'fraction_avg':fraction_avg_ls,
                                                'count_avg':count_avg_ls,
                                                cond:[key_cond]*len(number_ls)})])
    return df_cond_stat.drop_duplicates(subset=['edit','Description']).reset_index(drop=True)

def dms_comp(df: pd.DataFrame, cond: str, cond_comp: str, wt:str, res: int, sample='sample', edit='edit', psuedocount=1):
    ''' 
    dms_comp(): returns comparison DMS grid dataframe in tidy format split by condition
    
    Parameters:
    df (dataframe): fastq outcomes dataframe
    cond (str): Condition column name for grouping fastq outcomes dataframe
    cond_comp (str): Condition for comparison group
    wt (str): Expected wildtype nucleotide sequence (in frame AA)
    res (int): First AA number
    sample (str, optional): Sample column name for splicing fastq outcomes dataframe (Default: 'sample')
    edit (str, optional): Edit column name within fastq outcomes dataframe (Default: 'edit')
    psuedocount (int, optional): psuedocount to avoid log(0) & /0 (Default: 1)
    
    Dependencies: Bio.Seq.Seq, pandas, numpy, tidy, edit_1(), dms_cond(), & aa_props
    '''
    df_cond_stat = dms_cond(df,cond,wt,res,sample,edit,psuedocount) # Execute dms_cond()

    # Fold change & p-value relative comparison group
    print(f'Compute FC & pval relative to {cond_comp}:')
    df_stat = pd.DataFrame()
    df_comp = df_cond_stat[df_cond_stat[cond]==cond_comp] # Isolate comparison group
    df_other = df_cond_stat[df_cond_stat[cond]!=cond_comp] # From other groups
    for e in set(df_other[edit].tolist()): # iterate through edits
        print(f'{e}')
        df_other_edit = df_other[df_other[edit]==e]
        df_comp_edit = df_comp[df_comp[edit]==e]
        df_other_edit['fraction_avg_compare'] = [df_comp_edit.iloc[0]['fraction_avg']]*df_other_edit.shape[0]
        df_other_edit['count_avg_compare'] = [df_comp_edit.iloc[0]['count_avg']]*df_other_edit.shape[0]
        df_other_edit['FC'] = df_other_edit['fraction_avg']/df_comp_edit.iloc[0]['fraction_avg']
        ttests = [ttest_ind(other_fraction_ls,df_comp_edit.iloc[0]['fraction_ls']) 
                                 for other_fraction_ls in df_other_edit['fraction_ls']]
        df_other_edit['pval'] = [ttest[1] for ttest in ttests]
        df_other_edit['tstat'] = [ttest[0] for ttest in ttests]
        df_stat = pd.concat([df_stat,df_other_edit])
    df_stat['compare'] = [cond_comp]*df_stat.shape[0]
    return df_stat[[edit,'before','after','number','FC','pval','tstat','fraction_avg','fraction_avg_compare','count_avg','count_avg_compare',cond,'compare']].sort_values(by=['number','after']).reset_index(drop=True)

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

# Plot methods
def scat(typ: str,df: pd.DataFrame,x: str,y: str,cols=None,cols_ord=None,stys=None,cutoff=0.01,cols_exclude=None,
         file=None,dir=None,palette_or_cmap='colorblind',edgecol='black',
         figsize=(10,6),title='',title_size=18,title_weight='bold',
         x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,100),x_ticks_rot=0,xticks=[],
         y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,100),y_ticks_rot=0,yticks=[],
         legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0),show=True,
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
    show (bool, optional): show plot (Default: True)
    
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
           legend_title=legend_title,legend_title_size=legend_title_size,legend_size=legend_size,legend_bbox_to_anchor=legend_bbox_to_anchor,legend_loc=legend_loc,legend_items=legend_items,show=show, 
           **kwargs)

def cat(typ: str,df: pd.DataFrame,x: str,y: str,errorbar=None,cols=None,cols_ord=None,cutoff=0.01,cols_exclude=None,
        file=None,dir=None,palette_or_cmap='colorblind',edgecol='black',lw=1,
        figsize=(10,6),title='',title_size=18,title_weight='bold',
        x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_scale='linear',x_axis_dims=(0,1),x_ticks_rot=0,xticks=[],
        y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_scale='linear',y_axis_dims=(0,1),y_ticks_rot=0,yticks=[],
        legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_items=(0,0),show=True,
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
    show (bool, optional): show plot (Default: True)
    
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
          legend_title=legend_title,legend_title_size=legend_title_size,legend_size=legend_size,legend_bbox_to_anchor=legend_bbox_to_anchor,legend_loc=legend_loc,legend_items=legend_items,show=show, 
          **kwargs)

def stack(df: pd.DataFrame,x='sample',y='fraction',cols='edit',cutoff=0.01,cols_ord=[],x_ord=[],
          file=None,dir=None,cmap='Set2',
          title='Editing Outcomes',title_size=18,title_weight='bold',
          figsize=(10,6),x_axis='',x_axis_size=12,x_axis_weight='bold',x_ticks_rot=45,x_ticks_ha='right',
          y_axis='',y_axis_size=12,y_axis_weight='bold',y_ticks_rot=0,
          legend_title='',legend_title_size=12,legend_size=12,
          legend_bbox_to_anchor=(1,1),legend_loc='upper left',legend_ncol=1,show=True,**kwargs):
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
    show (bool, optional): show plot (Default: True)
    
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
            legend_bbox_to_anchor=legend_bbox_to_anchor,legend_loc=legend_loc,legend_ncol=legend_ncol,show=show,**kwargs)

def heat(df: pd.DataFrame, cond: str,x='number',y='after',vals='fraction_avg',vals_dims:tuple=None,
         file=None,dir=None,edgecol='black',lw=1,annot=False,cmap="bone_r",sq=True,cbar=True,
         title='',title_size=12,title_weight='bold',figsize=(20,7),
         x_axis='',x_axis_size=12,x_axis_weight='bold',x_ticks_rot=45,
         y_axis='',y_axis_size=12,y_axis_weight='bold',y_ticks_rot=0,
         show=True,**kwargs):
    ''' 
    heat(): creates heatmap
    
    Parameters:
    df (dataframe): tidy-formatted DMS dataframe (dms_cond() or dms_comp())
    x (str, optional): x-axis column name (AA residues number column)
    y (str, optional): y-axis column name (AA change column)
    vals (str, optional): values column name
    vals_dims (tuple, optional): vals minimum and maximum formatted (vmin, vmax; Default: None)
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    edgecol (str, optional): point edge color
    lw (int, optional): line width
    annot (bool, optional): annotate values
    cmap (str, optional): matplotlib color map
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    figsize (tuple, optional): figure size per subplot
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_ticks_rot (int, optional): x-axis ticks rotation
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_ticks_rot (int, optional): y-axis ticks rotation
    show (bool, optional): show plot (Default: True)
    
    Dependencies: matplotlib, seaborn, pandas, & aa_props
    '''
    # Find min and max values in the dataset for normalization
    if vals_dims is None:
        vmin = df[vals].values.min()
        vmax = df[vals].values.max()
    else:
        vmin = vals_dims[0]
        vmax = vals_dims[1]

    # Make DMS grids
    print('Make DMS grids')
    dc=t.split(df,cond) # Split by condition
    dc2={key:pd.pivot(df_cond,columns=x,index=y,values=vals).astype(float).reindex(list(aa_props.keys())) 
         for key,df_cond in dc.items()} # Generate pivot tables
    
    # Create a single figure with multiple heatmap subplots
    print('Create a single figure with multiple heatmap subplots')
    fig, axes = plt.subplots(nrows=len(list(dc2.keys())),ncols=1,figsize=(figsize[0],figsize[1]*len(list(dc2.keys()))),sharex=False,sharey=True)
    if isinstance(axes, np.ndarray)==False: axes = np.array([axes]) # Make axes iterable if there is only 1 heatmap
    for (ax, key) in zip(axes, list(dc2.keys())):
        print(f'{key}')
        sns.heatmap(dc2[key],annot=annot,cmap=cmap,ax=ax,linecolor=edgecol,linewidths=lw,cbar=cbar,square=sq,vmin=vmin,vmax=vmax, **kwargs)
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
        io.mkdir(dir) # Make output directory if it does not exist
        plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight', format=f'{file.split(".")[-1]}')
    if show: plt.show()

def vol(df: pd.DataFrame,x: str,y: str,size:str=None,size_dims:tuple=None,include_wt=False,
        file=None,dir=None,palette_or_cmap='YlOrRd',edgecol='black',
        figsize=(10,6),title='',title_size=18,title_weight='bold',
        x_axis='',x_axis_size=12,x_axis_weight='bold',x_axis_dims=(0,0),x_ticks_rot=0,xticks=[],
        y_axis='',y_axis_size=12,y_axis_weight='bold',y_axis_dims=(0,0),y_ticks_rot=0,yticks=[],
        legend_title='',legend_title_size=12,legend_size=9,legend_bbox_to_anchor=(1,1),legend_loc='upper left',
        legend_items=(0,0),legend_ncol=1,display_size=True,display_labels=True,return_df=True,show=True,
        **kwargs):
    ''' 
    vol(): creates volcano plot
    
    Parameters:
    df (dataframe): pandas dataframe
    x (str): x-axis column name
    y (str): y-axis column name
    cols (str, optional): color column name
    size (str, optional): size column name
    size_dims (tuple, optional): (minimum,maximum) values in size column (Default: None)
    include_wt (bool, optional): include wildtype (Default: False)
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
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_ticks_rot (int, optional): x-axis ticks rotation
    xticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    yticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    display_size (bool, optional): display size on plot (Default: True)
    display_labels (bool, optional): display labels for significant values (Default: True)
    return_df (bool, optional): return dataframe (Default: True)
    show (bool, optional): show plot (Default: True)
    
    Dependencies: os, matplotlib, seaborn, pandas, & edit_1()
    '''
    # Strings with subscripts
    log2 = 'log\u2082'
    log10 = 'log\u2081\u2080'
    
    # Log transform data
    df[f'{log2}({x})'] = [np.log10(xval)/np.log10(2) for xval in df[x]]
    df[f'-{log10}({y})'] = [-np.log10(yval) for yval in df[y]]
    
    # Organize data by significance
    signif = []
    for (log2FC,log10P) in zip(df[f'{log2}({x})'],df[f'-{log10}({y})']):
        if (np.abs(log2FC)>1)&(log10P>-np.log10(0.05)): signif.append('FC & p-value')
        elif (np.abs(log2FC)<=1)&(log10P>-np.log10(0.05)): signif.append('p-value')
        elif (np.abs(log2FC)>1)&(log10P<=-np.log10(0.05)): signif.append('FC')
        else: signif.append('NS')
    df['Significance']=signif
    signif_order = ['NS','FC','p-value','FC & p-value']

    # Organize data by conservation (changed from)
    basic = ['R','K', 'H']
    acidic = ['D','E']
    polar = ['S', 'T', 'N', 'Q', 'Y', 'C']
    nonpolar = ['A','V','L','I','M','F','W','P','G']
    change = []

    df = edit_1(df)
    for (before,after) in zip(df['before'],df['after']):
        if (before in basic)&(after not in basic): change.append('Basic')
        elif (before in acidic)&(after not in acidic): change.append('Acidic')
        elif (before in polar)&(after not in polar): change.append('Polar')
        elif (before in nonpolar)&(after not in nonpolar): change.append('Nonpolar')
        else: change.append('Conserved')
    df['Change'] = change

    sty_order = ['Conserved','Basic','Acidic','Polar','Nonpolar']
    mark_order = ['D','^','v','<','>']

    # Remove wildtype
    if include_wt==False:
        wt_i = [i for i,(before,after) in enumerate(t.zip_cols(df=df,cols=['before','after'])) if before == after]
        df = df.drop(wt_i,axis=0).reset_index(drop=True)

    # Organize data by abundance
    sizes=(1,100)
    if size_dims is not None: df = df[(df[size]>=size_dims[0])&(df[size]<=size_dims[1])]

    # Set dimensions
    if x_axis_dims==(0,0): x_axis_dims=(min(df[f'{log2}({x})']),max(df[f'{log2}({x})']))
    if y_axis_dims==(0,0): y_axis_dims=(0,max(df[f'-{log10}({y})']))

    # Generate figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # with significance boundraries
    plt.vlines(x=-1, ymin=y_axis_dims[0], ymax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
    plt.vlines(x=1, ymin=y_axis_dims[0], ymax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
    plt.hlines(y=-np.log10(0.05), xmin=x_axis_dims[0], xmax=x_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
    
    # with data
    if display_size==False: size=None
    sns.scatterplot(data=df, x=f'{log2}({x})', y=f'-{log10}({y})', 
                    hue='Significance', hue_order=signif_order, 
                    edgecolor=edgecol, palette=palette_or_cmap, style='Change',
                    style_order=sty_order,markers=mark_order,size=size,
                    sizes=sizes,
                    ax=ax, **kwargs)
    
    # with labels
    if display_labels:
        df_signif = df[df['Significance']=='FC & p-value']
        adjust_text([plt.text(x=df_signif.iloc[i][f'{log2}({x})'], 
                              y=df_signif.iloc[i][f'-{log10}({y})'],
                              s=edit) for i,edit in enumerate(df_signif['edit'])])

    # Set title
    if title=='' and file is not None: title=p.re_un_cap(".".join(file.split(".")[:-1]))
    plt.title(title, fontsize=title_size, fontweight=title_weight)
    
    # Set x axis
    if x_axis=='': x_axis=f'{log2}({x})'
    plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight)
    if xticks==[]: 
        if (x_ticks_rot==0)|(x_ticks_rot==90): plt.xticks(rotation=x_ticks_rot,ha='center')
        else: plt.xticks(rotation=x_ticks_rot,ha='right')
    else: 
        if (x_ticks_rot==0)|(x_ticks_rot==90): plt.xticks(ticks=xticks,labels=xticks,rotation=x_ticks_rot, ha='center')
        else: plt.xticks(ticks=xticks,labels=xticks,rotation=x_ticks_rot,ha='right')

    # Set y axis
    if y_axis=='': y_axis=f'-{log10}({y})'
    plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight)

    if yticks==[]: plt.yticks(rotation=y_ticks_rot)
    else: plt.yticks(ticks=yticks,labels=yticks,rotation=y_ticks_rot)

    # Move legend to the right of the graph
    if legend_items==(0,0): ax.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
                                        bbox_to_anchor=legend_bbox_to_anchor,loc=legend_loc,ncol=legend_ncol)
    else: 
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
                  bbox_to_anchor=legend_bbox_to_anchor,loc=legend_loc,ncol=legend_ncol, # Move right of the graph
                  handles=handles[legend_items[0]:legend_items[1]],labels=labels[legend_items[0]:legend_items[1]]) # Only retains specified labels

    # Save & show fig; return dataframe
    if file is not None and dir is not None:
        io.mkdir(dir) # Make output directory if it does not exist
        plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight', format=f'{file.split(".")[-1]}')
    if show: plt.show()
    if return_df: return df     