''' 
Module: cli.py
Author: Marc Zepeda
Created: 2025-4-01
Description: Command Line Interaction

Usage:
[Common commands for fastq files]
- access(): make all files and subdirectories accessible on Harvard FASRC
- smaller_fastq(): create new subdirectory containing fastqs with the # of reads limited
'''

# Import packages
import subprocess
import os

from . import io

# Common commands for fastq files
def access(pt: str):
    ''' 
    access(): make all files and subdirectories accessible on Harvard FASRC
    
    Parameters:
    pt (str): path to parent directory

    Dependencies: subprocess
    '''
    # Run command in the directory
    command = 'chmod g+r . ; chmod g+rwxs -R . ; chmod g+x .'
    print(f"terminal:\ncd {pt}\n{command}")
    result = subprocess.run(f"{command}", shell=True, cwd=pt, capture_output=True, text=True)
    
    # Print output
    if result.stdout: print(f"output:\n{result.stdout}")
    if result.stderr: print(f"errors:\n{result.stderr}")

def smaller_fastq(pt: str, reads: int, suf: str='.fastq.gz'):
    '''
    smaller_fastq(): create new subdirectory containing fastqs with the # of reads limited

    Parameters:
    pt (str): path to parent directory containing fastq files
    reads (int): maximum # of reads per fastq file
    suf (str, optional): fastq file suffix (Default: .fastq.gz)

    Dependencies: subprocess,os
    '''
    # Get fastq files
    files = os.listdir(pt)
    fastq_files = [file for file in files if suf in file]

    # Make output directory
    out_dir = os.path.join(pt,f'{reads}_reads')
    io.mkdir(out_dir)

    # Run commands in the directory
    print(f"terminal:\ncd {pt}")
    if suf=='.fastq.gz': # gzipped fastq files
        for fastq_file in fastq_files: # Iterate through fastqs

            # Run command
            command = f'gunzip -c {fastq_file} | head -n {4*reads} > {out_dir}/{fastq_file[:-3]}'
            print(f"{command}")
            result = subprocess.run(f"{command}", shell=True, cwd=pt, capture_output=True, text=True)
            
            # Print output/errors
            if result.stdout: print(f"output:\n{result.stdout}")
            if result.stderr: print(f"errors:\n{result.stderr}")

    else: # unzipped fastq files
        for fastq_file in fastq_files: # Iterate through fastqs

            command = f'head -n {4*reads} {fastq_file} > {out_dir}/{fastq_file}'
            print(f"{command}")
            result = subprocess.run(f"{command}", shell=True, cwd=pt, capture_output=True, text=True)
            
            # Print output/errors
            if result.stdout: print(f"output:\n{result.stdout}")
            if result.stderr: print(f"errors:\n{result.stderr}")