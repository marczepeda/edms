''' 
Module: cli.py
Author: Marc Zepeda
Created: 2025-4-01
Description: Command Line Interaction

Usage:
[Common commands for fastq files]
- access(): make all files and subdirectories accessible on Harvard FASRC
- smaller_fastq(): create new subdirectory containing fastqs with the # of reads limited

[Environment variable management]
- detect_shell(): Attempt to detect the user's shell from the SHELL environment variable (i.e., 'bash' or 'zsh').
- create_export_var(): create a persistent export variable by adding it to the user's shell config.
- view_export_vars(): View the current export variables in the user's shell config.
'''

# Import packages
import subprocess
import os
from pathlib import Path
from typing import Literal

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

# Environment variable management
home = Path.home()
shell_configs = {
    'bash': home / '.bashrc',
    'zsh': home / '.zshrc'
}

def detect_shell():
    """
    detect_shell(): Attempt to detect the user's shell from the SHELL environment variable (i.e., 'bash' or 'zsh').
    """
    shell_path = os.environ.get("SHELL", "")
    shell_name = Path(shell_path).name.lower()
    if shell_name in {"bash", "zsh"}:
        return shell_name
    raise ValueError(f"Unsupported or undetected shell: {shell_path}")            

def create_export_var(name: str, pt: str, shell: Literal['bash','zsh']=None):
    """
    create_export_var(): create a persistent export variable by adding it to the user's shell config.
    
    Parameters:
    name (str): The name of the environment variable (e.g. "MYPROJ").
    path (str): The full path the variable should point to.
    shell (str, optional): The shell config to use ('bash' or 'zsh'). If None, auto-detects the shell.
    """
    shell = shell or detect_shell()

    config_file = shell_configs.get(shell)
    if config_file is None:
        raise ValueError(f"Unsupported shell: {shell}")

    # Clean path and environment variable name
    pt = os.path.abspath(os.path.expanduser(pt))
    name = name.upper()

    export_line = f'export {name}="{pt}"\n'

    # Prevent duplicate entries
    if config_file.exists():
        with open(config_file, 'r') as f:
            if export_line.strip() in f.read():
                print(f"Environment variable {name} already exists in {config_file}")
                return

    # Append the new variable
    with open(config_file, 'a') as f:
        f.write(f'\n# Added by script\n{export_line}')

    print(f"Added environment variable: {name}={pt}")
    print(f"Run `source {config_file}` or restart your terminal to apply it.")

def view_export_vars(shell: Literal['bash','zsh']=None):
    """
    view_export_vars(): View the current export variables in the user's shell config.
    
    Parameters:
    shell (str, optional): The shell config to check ('bash' or 'zsh'). If None, auto-detects the shell.
    """
    shell = shell or detect_shell()

    config_file = shell_configs.get(shell)
    if not config_file or not config_file.exists():
        print(f"Shell config file not found for {shell}")
        return

    print(f"Exported environment variables in {config_file}:\n")

    with open(config_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("export "):
                print(line)