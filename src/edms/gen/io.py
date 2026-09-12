'''
Module: io.py
Author: Marc Zepeda
Created: 2024-05-20
Description: Input/Output

Usage:
[Parsing Python literals]
- df_try_parse(): apply try_parse() to dataframe columns and return dataframe
- recursive_json_decode(): recursively decode JSON strings in a nested structure

[Input]
- get(): returns pandas dataframe from a file
- get_dir(): returns python dictionary of dataframe from files within a directory

[Output]
- save(): save .csv file to a specified output file path from obj
- save_dir(): save .csv files to a specified output directory from dictionary of objs

[Input/Output]
- excel_csvs(): exports excel file to .csv files in specified directory
- df_to_dc_txt(): returns pandas DataFrame as a printed text that resembles a Python dictionary
- dc_txt_to_df(): returns a pandas DataFrame from text that resembles a Python dictionary
- in_subs(): moves all files with a given suffix into subfolders named after the files (excluding the suffix).
- out_subs(): delete subdirectories and move their files to the parent directory
- create_sh(): creates a shell script with SLURM job submission parameters for Harvard FASRC cluster.
- create_pipeline(): generates SLURM shell scripts and a submit.sh driver from a pipeline .csv
- combine(): Combine text files matching provided suffixes into a single output file, inserting a header with the original filename before each file's content.
- split_R1_R2(): split paired reads into new R1 and R2 subdirectories at the parent directory
- basespace(): reorganize an Illumina BaseSpace download into a parsable format

[Directory Methods]
- relative_paths(): returns relative paths for all files in a directory including subfolders
- sorted_file_names: returns sorted file names in a directory with the specified suffix
'''

# Import packages
import pandas as pd
import os
import ast
import csv
import shutil
import datetime
import re
import shlex
import json
import importlib.resources
from typing import Literal, Iterable
from pathlib import Path 


from ..utils import try_parse, mkdir, check_outpath
from ..config import get_info
from ..gen import tidy as t

# Parsing Python literals
def df_try_parse(df: pd.DataFrame) -> pd.DataFrame:
    '''
    df_try_parse(): apply try_parse() to dataframe columns and return dataframe

    Parameters: 
    df (dataframe): dataframe with columns to try_parse()

    Dependencies: utils.try_parse()
    '''
    # Apply the parsing function to all columns
    for col in df.columns:
        df[col] = df[col].apply(try_parse)
    return df

def recursive_json_decode(obj):
    '''
    recursive_json_decode(): recursively decode JSON strings in a nested structure

    Parameters:
    obj: object of any type (looking for dict, list, or str)
    
    Dependencies: json
    '''
    if isinstance(obj, dict):
        return {k: recursive_json_decode(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [recursive_json_decode(item) for item in obj]
    elif isinstance(obj, str):
        try:
            decoded = json.loads(obj)
            return recursive_json_decode(decoded)
        except (json.JSONDecodeError, TypeError):
            return obj
    return obj

# Input
def get(pt: str, literal_eval:bool=False, **kwargs) -> pd.DataFrame | dict[pd.DataFrame]:
    ''' 
    get(): returns pandas dataframe from a file
    
    Parameters:    
    pt (str): file path
    literal_evals (bool, optional): convert Python literals encoded as strings (Default: False)
        (1) automatically detects and parses columns containing Python literals (e.g., dict, list, set, tuple) encoded as strings
        (2) recursively evaluates nested structures
    **kwargs: pandas.read_csv() parameters
    
    Dependencies: pandas,ast,utils[try_parse(),recursive_parse()],df_try_parse()
    '''
    suf = pt.split('.')[-1]
    if suf=='csv': 
        if literal_eval: df = df_try_parse(pd.read_csv(filepath_or_buffer=pt,sep=',',**kwargs))
        else: df = pd.read_csv(filepath_or_buffer=pt,sep=',',**kwargs)
        print(f"CSV file: {pt}\nColumns: {', '.join([col for col in df.columns])}")
        return df

    elif suf=='tsv': 
        if literal_eval: df = df_try_parse(pd.read_csv(filepath_or_buffer=pt,sep='\t',**kwargs))
        else: df = pd.read_csv(filepath_or_buffer=pt,sep='\t',**kwargs)
        print(f"TSV file: {pt}\nColumns: {', '.join([col for col in df.columns])}")
        return df

    elif suf=='xlsx': 
        if literal_eval: dc = {sheet_name: df_try_parse(pd.read_excel(pt,sheet_name,**kwargs))
                                for sheet_name in pd.ExcelFile(pt).sheet_names}
        else: dc = {sheet_name: pd.read_excel(pt,sheet_name,**kwargs)
                                for sheet_name in pd.ExcelFile(pt).sheet_names}
        print(f"Excel file: {pt}\nKeys: {', '.join([key for key in dc.keys()])}")
        for key,df in dc.items():
            print(f"Key: {key}\nColumns: {', '.join([col for col in df.columns])}")
        return dc

    elif suf=='html': 
        if literal_eval: return df_try_parse(pd.read_html(pt,**kwargs))
        else: return pd.read_html(pt,**kwargs)
    else: 
        if literal_eval: return df_try_parse(pd.read_csv(filepath_or_buffer=pt,**kwargs))
        else: return pd.read_csv(filepath_or_buffer=pt,**kwargs)
    
def get_dir(dir: str, suf: str='.csv', literal_eval: bool=False, **kwargs) -> dict[pd.DataFrame]:
    ''' 
    get_dir(): returns python dictionary of dataframe from files within a directory
    
    Parameters:
    dir (str): directory path with files
    suf (str): file type suffix
    literal_evals (bool, optional): convert Python literals encoded as strings (Default: False)
        (1) automatically detects and parses columns containing Python literals (e.g., dict, list, set, tuple) encoded as strings
        (2) recursively evaluates nested structures
    **kwargs: pandas.read_csv() parameters
    
    Dependencies: pandas
    '''
    files = [file for file in os.listdir(dir) if file[-len(suf):]==suf]
    dc = {file[:-len(suf)]:get(os.path.join(dir,file),literal_eval,**kwargs) for file in files}
    print(f"Directory: {dir}\nKeys: {', '.join([key for key in dc.keys()])}")
    for key,df in dc.items():
        print(f"Key: {key}\nColumns: {', '.join([col for col in df.columns])}")
    return dc

# Output
def save(obj, file: str, cols: list=[], id: bool=False, sort: bool=True, **kwargs):
    ''' 
    save(): save .csv file to a specified output file path from obj
    
    Parameters:
    obj: dataframe, series, set, or list
    file (str, required): output file path; the output directory is created if needed
    cols (str, list, optional): isolate dataframe column(s)
    id (bool, optional): include dataframe index (False)
    sort (bool, optional): sort set, list, or series before saving (True)
    
    Dependencies: pandas, os, csv & utils.check_outpath()
    '''
    pt = check_outpath(file=file) # Check save path & make output directory
    if pt is None:
        print(f"Warning: Invalid file path; file not saved.\nFile: {file}")
        return
    name = os.path.basename(pt) # File name used for suffix checks & Excel sheet names

    if type(obj)==pd.DataFrame:
        for col in cols: # Check if each element in the list is a string
            if not isinstance(col, str):
                raise ValueError("All elements in the list must be strings.")
        if cols!=[]: obj = obj[cols]
        if name.split('.')[-1]=='tsv': obj.to_csv(pt, index=id, sep='\t', **kwargs)
        elif name.split('.')[-1]=='xlsx': 
            with pd.ExcelWriter(pt) as writer: 
                obj.to_excel(writer,sheet_name='.'.join(name.split('.')[:-1]),index=id) # Dataframe per sheet
        else: obj.to_csv(pt, index=id, **kwargs)
    elif type(obj)==set or type(obj)==list or type(obj)==pd.Series:
        if sort==True: obj2 = sorted(list(obj))
        else: obj2=list(obj)
        with open(pt, 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file, dialect='excel') # Create a CSV writer object
            csv_writer.writerow(obj2) # Write each row of the list to the CSV file
    elif (type(obj)==dict)&(name.split('.')[-1]=='xlsx'):
        for col in cols: # Check if each element in the list is a string
            if not isinstance(col, str):
                raise ValueError("All elements in the list must be strings.")
        with pd.ExcelWriter(pt) as writer:
            if cols!=[]: obj = obj[cols]
            for key,df in obj.items(): 
                if cols!=[]: df = df[cols]
                df.to_excel(writer,sheet_name=key,index=id) # Dataframe per sheet
    else: raise ValueError(f'save() does not work for {type(obj)} objects with {name.split(".")[-1]} files.')

def save_dir(dir: str, dc: dict, suf: str='.csv', **kwargs):
    ''' 
    save_dir(): save .csv files to a specified output directory from dictionary of objs
    
    Parameters:
    dir (str): output directory path
    dc (dict): dictionary of objects (files)
    suf (str, optional): file name suffix (Default: .csv)

    Dependencies: pandas, os, csv, & save()
    '''
    for key,val in dc.items(): save(obj=val, file=os.path.join(dir,f"{key}{suf}"),**kwargs)

# Input/Output
def excel_csvs(pt: str, dir: str='', **kwargs):
    ''' 
    excel_csvs(): exports excel file to .csv files in specified directory    
    
    Parameters:
    pt (str): excel file path
    dir (str, optional): output directory path (Default: same directory as excel file name)
    
    Dependencies: pandas, os, & mkdir
    '''
    if dir=='': dir = '.'.join(pt.split('.')[:-1]) # Get the directory where the Excel file is located
    mkdir(dir) # Make output directory if it does not exist
    for sheet_name in pd.ExcelFile(pt).sheet_names: # Loop through each sheet in the Excel file
        df = pd.read_excel(pd.ExcelFile(pt),sheet_name,**kwargs) # Read the sheet into a DataFrame
        df.to_csv(os.path.join(dir,f"{sheet_name}.csv"),index=False) # Save the DataFrame to a CSV file

def df_to_dc_txt(df: pd.DataFrame) -> str:
    ''' 
    df_to_dc_txt(): returns pandas DataFrame as a printed text that resembles a Python dictionary
    
    Parameters:
    df (dataframe): pandas dataframe
    
    Dependencies: pandas
    '''
    dict_text = "{\n"
    for index, row in df.iterrows():
        dict_text += f"  {index}: {{\n"
        for col in df.columns:
            value = row[col]
            if isinstance(value, str):
                value = f"'{value}'"
            dict_text += f"    '{col}': {value},\n"
        dict_text = dict_text.rstrip(",\n") + "\n  },\n"  # Remove trailing comma for last key-value pair
    dict_text = dict_text.rstrip(",\n") + "\n}"  # Close the main dictionary
    print(dict_text)
    return dict_text

def dc_txt_to_df(dc_txt: str, transpose: bool=True) -> str:
    ''' 
    dc_txt_to_df(): returns a pandas DataFrame from text that resembles a Python dictionary
    
    Parameters:
    dc_txt (str): text that resembles a Python dictionary
    transpose (bool, optional): transpose dataframe (True)
    
    Dependencies: pandas & ast
    '''
    if transpose==True: return pd.DataFrame(ast.literal_eval(dc_txt)).T
    else: return pd.DataFrame(ast.literal_eval(dc_txt))

def in_subs(dir: str, suf: str, 
            group_by: Literal["basename", "prefix"] = "basename",
            prefix_sep: str | None = "_",
            prefix_len: int | None = None,): 
    '''
    in_subs: moves all files with a given suffix into subdirectory named after the files (excluding the suffix).

    Parameters:
    dir (str): Path to the directory containing the files.
    suf (str): File suffix (e.g., '.txt', '.csv') to filter files.
    group_by (Literal, optional): How to group files into subdirectories (Options: "basename" [Default] or "prefix").
        - "basename": Use the full filename without the suffix as the subdirectory name.
        - "prefix": Use a prefix of the filename (before a specified separator or of a specified length) as the subdirectory name.
    prefix_sep (str, optional): Delimiter used to extract prefix (Default: "_"; e.g., sample_001.txt → sample). Ignored if prefix_len is provided.
    prefix_len (int, optional): Number of characters to use as prefix. Overrides prefix_sep if provided.

    Dependences: os, shutil
    '''
    if not os.path.isdir(dir):
        raise ValueError(f"{dir} is not a valid directory.")

    if group_by == "prefix" and prefix_sep is None and prefix_len is None:
        raise ValueError("When group_by='prefix', either prefix_sep or prefix_len must be provided.")

    for filename in os.listdir(dir):
        file_path = os.path.join(dir, filename)

        if not (os.path.isfile(file_path) and filename.endswith(suf)):
            continue

        stem = filename[:-len(suf)]

        if group_by == "basename":
            subdir_name = stem

        elif group_by == "prefix":
            if prefix_len is not None:
                subdir_name = stem[:prefix_len]
            else:
                subdir_name = stem.split(prefix_sep, 1)[0]

        else:
            raise ValueError(f"Unknown group_by mode: {group_by}")

        subdir = os.path.join(dir, subdir_name)
        os.makedirs(subdir, exist_ok=True)

        shutil.move(file_path, os.path.join(subdir, filename))

def out_subs(dir: str):
    """
    out_subs(): Delete subdirectories and move their files to the parent directory.

    Parameters:
    dir (str): Path to the directory containing the files.
    """
    if not os.path.isdir(dir):
        raise ValueError(f"{dir} is not a valid directory.")

    parent_dir = dir

    for root, dirs, files in os.walk(parent_dir, topdown=False):
        for file in files:
            file_path = os.path.join(root, file)
            target_path = os.path.join(parent_dir, file)

            # Already at the top level - nothing to move, and it must not
            # be treated as a "conflict" with itself.
            if os.path.abspath(file_path) == os.path.abspath(target_path):
                continue

            # Resolve name conflicts by appending a counter
            base, ext = os.path.splitext(file)
            counter = 1
            while os.path.exists(target_path):
                target_path = os.path.join(parent_dir, f"{base}_{counter}{ext}")
                counter += 1

            shutil.move(file_path, target_path)

        # Remove empty directories
        for subdir in dirs:
            dir_path = os.path.join(root, subdir)
            if os.path.isdir(dir_path) and not os.listdir(dir_path):
                os.rmdir(dir_path)

def create_sh(file: str,
              cores: int = 1, partition: str='serial_requeue', time: str = '0-00:10', mem: int = 1000, email: str=None,
              python: str = 'python', env: str = 'edms',
              cmd: str = None, log_prefix: str = None):
    '''
    create_sh(): creates a shell script with SLURM job submission parameters for Harvard FASRC cluster.

    Parameters:
    file (str): Output shell script path (e.g., "./scripts/script.sh"); the output directory is created if needed.
    cores (int, optional): The number of cores to request for the job (Default: 1).
    partition (str, optional): The partition to submit the job to (Default: 'serial_requeue').
    time (str, optional): The maximum runtime for the job in D-HH:MM format (Default: '0-00:10').
    mem (int, optional): The amount of memory to request for the job in MB (Default: 1000).
    email (str, optional): The email address (Default: None = get_info('Harvard')['email']).
    python (str, optional): The python module to load (Default: 'python').
    env (str, optional): The conda environment to activate (Default: 'edms').
    cmd (str, optional): Command to run in the body of the script (Default: None = run the
        same-named python script, "python {file stem}.py").
    log_prefix (str, optional): Prefix for the .out/.err SLURM log files (Default: None =
        "{timestamp}_{file stem}").

    Dependencies: os, datetime, utils.check_outpath(), config.get_info()
    '''
    # Check if the file is valid
    if not str(file).endswith('.sh'):
        raise ValueError("File must end with '.sh'")

    # Resolve the output path and make the output directory
    pt = check_outpath(file=file)
    name = os.path.basename(pt)

    # Get email from config if not provided
    if email is None:
        try:
            email = get_info('Harvard')['email']
        except Exception as e:
            print(f"Error retrieving email from config: {e}")
            return

    stem = ".".join(name.split(".")[:-1])

    # Body of the script: an explicit command, or the same-named python script
    if cmd is None:
        body = f'python {stem}.py \t# Run python script'
    else:
        body = cmd

    # Prefix for the SLURM .out/.err logs
    if log_prefix is None:
        log_prefix = f'{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_{stem}'

    # Create .sh file
    try:
        with open(pt, 'w') as file_obj:
            file_obj.write(f'''#!/bin/bash
#
#SBATCH -n {cores} \t# Number of cores
#SBATCH -N 1 \t# Ensure that all cores are on one machine
#SBATCH -t {time} \t# Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p {partition} \t# Partition to submit to
#SBATCH --mem={mem} \t# Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o {log_prefix}_%j.out \t# File to which STDOUT will be written, %j inserts jobid
#SBATCH -e {log_prefix}_%j.err \t# File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-type=ALL \t# Email
#SBATCH --mail-user={email} \t# Email

echo -e "File: {name}\\nTime: {time}\\nMemory: {mem} MB" \t# Print job parameters
module load {python} \t# Load python module
mamba activate {env} \t# Activate conda environment
export PYTHONUNBUFFERED=1 \t# Ensure prints from python script are written to .out file
{body}
''')
    except Exception as e:
        print(f"An error occurred while creating the shell script: {e}")

# Annotated starting template shipped in edms.resources
PIPELINE_EXAMPLE = 'pipeline_example.csv'

# Reserved pipeline .csv columns; every other column becomes a command line flag
PIPELINE_RESERVED = {
    'step', 'sample', 'command', 'scatter', 'depends_on', 'script',
    'sbatch_cores', 'sbatch_partition', 'sbatch_time', 'sbatch_mem',
    'sbatch_email', 'sbatch_python', 'sbatch_env',
}

# Only these spellings are treated as booleans, so that numeric arguments
# such as "-n 0" or "--align_ckpt 1" are never mistaken for a store_true flag.
PIPELINE_TRUE = {'true', 'yes'}
PIPELINE_FALSE = {'false', 'no'}

def _pipeline_cell(value) -> str:
    '''
    _pipeline_cell(): returns a stripped string for a pipeline .csv cell ('' when blank).

    Parameters:
    value: cell value from the pipeline dataframe
    '''
    if value is None:
        return ''
    return str(value).strip()

def _pipeline_flag(col: str, value: str) -> list:
    '''
    _pipeline_flag(): returns command line tokens for one (column, value) pair.

    Column names starting with '-' are used verbatim ('-q'); all others are
    prefixed ('fastq_dir' -> '--fastq_dir'). 'TRUE'/'yes' renders the flag on its
    own (store_true); 'FALSE'/'no' omits it entirely.

    Parameters:
    col (str): pipeline .csv column name
    value (str): pipeline .csv cell value
    '''
    flag = col if col.startswith('-') else f'--{col}'
    low = value.lower()
    if low in PIPELINE_FALSE:
        return []
    if low in PIPELINE_TRUE:
        return [flag]
    return [flag, value]

def _pipeline_var(name: str) -> str:
    '''
    _pipeline_var(): returns a name usable as a bash variable.

    Parameters:
    name (str): step name
    '''
    var = re.sub(r'\W', '_', name)
    if not var or var[0].isdigit():
        var = f's_{var}'
    return var

def _pipeline_split(value: str) -> list:
    '''
    _pipeline_split(): returns a list from a ';' or ',' delimited cell.

    Parameters:
    value (str): pipeline .csv cell value
    '''
    return [v.strip() for v in re.split(r'[;,]', value) if v.strip()]

def create_pipeline(pt: str = None, dir: str = '.', submit_file: str = 'submit.sh',
                    cores: int = 1, partition: str = 'serial_requeue', time: str = '0-00:10',
                    mem: int = 1000, email: str = None, python: str = 'python', env: str = 'edms',
                    prog: str = 'edms', dry_run: bool = False, example: bool = False):
    '''
    create_pipeline(): generates one SLURM shell script per job, plus a submit.sh that
    chains the steps with afterok dependencies, from a pipeline .csv.

    The .csv holds two kinds of rows:
    - Step rows ('sample' blank): one pipeline stage. 'scatter' fans the stage out into one
      script per sample, substituting {sample} into any value.
    - Override rows ('step' and 'sample' both filled): replaces individual fields for that
      one sample of that step (e.g. give one sample more memory).

    Reserved columns:
    - step: stage name; also the script name for un-scattered stages (required)
    - sample: blank on a step row; the sample name on an override row
    - command: edms subcommand, e.g. 'fastq trim' (required on step rows)
    - scatter: ';' or ',' separated sample names, or '@file' to read one name per line
    - depends_on: ';' or ',' separated names of earlier steps to wait on (afterok)
    - script: script filename override; may contain {sample}
    - sbatch_cores / sbatch_partition / sbatch_time / sbatch_mem / sbatch_email /
      sbatch_python / sbatch_env: per-step SLURM and environment settings

    Every other column becomes a flag on the generated command line. Blank cells are
    skipped, so one wide .csv can cover steps that take different arguments. Columns
    whose name starts with '#' are notes for the reader and are ignored, as are rows
    whose 'step' starts with '#'.

    Parameters:
    pt (str): path to the pipeline .csv (or .tsv) file (required unless example=True).
    dir (str, optional): directory to write the scripts into (Default: '.').
    submit_file (str, optional): name of the submission driver script (Default: 'submit.sh').
    cores (int, optional): default number of cores (Default: 1).
    partition (str, optional): default partition (Default: 'serial_requeue').
    time (str, optional): default runtime in D-HH:MM (Default: '0-00:10').
    mem (int, optional): default memory in MB (Default: 1000).
    email (str, optional): default email (Default: None = get_info('Harvard')['email']).
    python (str, optional): default python module to load (Default: 'python').
    env (str, optional): default conda environment (Default: 'edms').
    prog (str, optional): program invoked by each generated command (Default: 'edms').
    example (bool, optional): copy the annotated example pipeline .csv into 'dir' and return,
        instead of generating anything (Default: False). Use it as a starting template.
    dry_run (bool, optional): print what would be generated without writing any files
        (Default: False). The .csv is still fully parsed and validated, so a dry run
        surfaces bad flags, bad dependencies, and colliding script names.

    Dependencies: os, re, shlex, pandas, datetime, utils.mkdir(), config.get_info(), create_sh()
    '''
    # Write the annotated starting template and stop
    if example:
        mkdir(dir)
        out_pt = os.path.join(dir, PIPELINE_EXAMPLE)
        with importlib.resources.files('edms.resources').joinpath(PIPELINE_EXAMPLE).open('r', encoding='utf-8') as fh:
            text = fh.read()
        with open(out_pt, 'w', encoding='utf-8') as fh:
            fh.write(text)
        print(f'Wrote example pipeline spec to {os.path.abspath(out_pt)}\n'
              f'Edit it, then preview with: edms io pipeline -i {PIPELINE_EXAMPLE} --dry_run')
        return out_pt

    if pt is None:
        raise ValueError("pt is required (or pass example=True to write the example .csv)")

    # Read the spec verbatim: dtype=str keeps 15000 from becoming 15000.0 and
    # keep_default_na=False keeps blank cells as '' rather than NaN
    sep = '\t' if str(pt).endswith(('.tsv', '.tab')) else ','
    df = pd.read_csv(pt, dtype=str, keep_default_na=False, sep=sep)
    df.columns = [str(c).strip() for c in df.columns]

    if 'step' not in df.columns:
        raise ValueError("Pipeline file must have a 'step' column")
    if 'command' not in df.columns:
        raise ValueError("Pipeline file must have a 'command' column")

    # Get email from config if not provided
    if email is None:
        try:
            email = get_info('Harvard')['email']
        except Exception as e:
            print(f"Error retrieving email from config: {e}")
            return

    if not dry_run:
        mkdir(dir)

    # '#' columns are notes for the reader; they never become flags
    flag_cols = [c for c in df.columns if c not in PIPELINE_RESERVED and not c.startswith('#')]

    # Separate step rows from per-sample override rows
    steps = []
    steps_by_name = {}
    overrides = []
    for i, row in df.iterrows():
        line = i + 2  # +1 for the header, +1 for 1-based line numbers
        vals = {c: _pipeline_cell(row[c]) for c in df.columns}
        step = vals['step']
        if not step or step.startswith('#'):
            continue  # blank spacer row, or a commented-out row
        sample = vals.get('sample', '')
        if sample:  # override row
            overrides.append((step, sample,
                              {c: v for c, v in vals.items() if v and c not in ('step', 'sample')}, line))
            continue
        if step in steps_by_name:
            raise ValueError(f"Duplicate step '{step}' on line {line}; "
                             f"fill in 'sample' to make it a per-sample override row")
        if not vals['command']:
            raise ValueError(f"Step '{step}' on line {line} is missing a 'command'")
        entry = {'name': step, 'vals': vals, 'line': line}
        steps.append(entry)
        steps_by_name[step] = entry

    if not steps:
        raise ValueError(f"No step rows found in {pt}")

    # Resolve each step's samples
    for entry in steps:
        scatter = entry['vals'].get('scatter', '')
        if scatter.startswith('@'):
            spt = scatter[1:]
            if not os.path.isabs(spt):
                spt = os.path.join(os.path.dirname(os.path.abspath(pt)), spt)
            with open(spt) as fh:
                samples = [ln.strip() for ln in fh if ln.strip() and not ln.lstrip().startswith('#')]
            if not samples:
                raise ValueError(f"Step '{entry['name']}' scatters over '{spt}', which is empty")
        elif scatter:
            samples = _pipeline_split(scatter)
        else:
            samples = [None]
        entry['samples'] = samples

    # Validate overrides point at a real step and a real sample of that step
    for (step, sample, _, line) in overrides:
        if step not in steps_by_name:
            raise ValueError(f"Override on line {line} names step '{step}', which has no step row")
        if sample not in steps_by_name[step]['samples']:
            raise ValueError(f"Override on line {line} names sample '{sample}', "
                             f"which step '{step}' does not scatter over")

    # Validate dependencies point at an earlier step
    seen = set()
    for entry in steps:
        for dep in _pipeline_split(entry['vals'].get('depends_on', '')):
            if dep not in steps_by_name:
                raise ValueError(f"Step '{entry['name']}' depends on '{dep}', which has no step row")
            if dep not in seen:
                raise ValueError(f"Step '{entry['name']}' depends on '{dep}', "
                                 f"which is not defined earlier in {os.path.basename(pt)}")
        seen.add(entry['name'])

    # Build every job (nothing is written until the whole spec has been validated)
    written = []
    scripts_by_step = {}
    jobs = []
    for entry in steps:
        scripts = []
        for sample in entry['samples']:
            vals = dict(entry['vals'])
            if sample is not None:
                for (ostep, osample, odict, _) in overrides:
                    if ostep == entry['name'] and osample == sample:
                        vals.update(odict)

            def sub(v, sample=sample):
                return v.replace('{sample}', sample) if sample is not None else v

            # Script name: explicit 'script' column, else the sample, else the step
            if vals.get('script'):
                stem = sub(vals['script'])
            elif sample is not None:
                stem = sample
            else:
                stem = entry['name']
            if stem.endswith('.sh'):
                stem = stem[:-3]
            if stem in written:
                raise ValueError(f"Two jobs both want to write '{stem}.sh'; add a 'script' column "
                                 f"(e.g. '{entry['name']}_{{sample}}') to keep the names distinct")

            # Command line: the subcommand, then one flag per continued line
            chunks = []
            for col in flag_cols:
                value = vals.get(col, '')
                if not value:
                    continue
                tokens = _pipeline_flag(col, sub(value))
                if tokens:
                    chunks.append(' '.join(shlex.quote(t) for t in tokens))
            base = f"{prog} {' '.join(vals['command'].split())}"
            cmd = base + (' \\\n' + ' \\\n'.join(f'    {c}' for c in chunks) if chunks else '')

            jobs.append(dict(step=entry['name'], file=f'{stem}.sh',
                             cores=int(vals.get('sbatch_cores') or cores),
                             partition=vals.get('sbatch_partition') or partition,
                             time=vals.get('sbatch_time') or time,
                             mem=int(vals.get('sbatch_mem') or mem),
                             email=vals.get('sbatch_email') or email,
                             python=vals.get('sbatch_python') or python,
                             env=vals.get('sbatch_env') or env,
                             cmd=cmd, log_prefix=stem))
            written.append(stem)
            scripts.append(f'{stem}.sh')
        scripts_by_step[entry['name']] = scripts

    # Build submit.sh
    lines = [
        '#!/bin/bash',
        f'# {submit_file}: generated by `edms io pipeline` on '
        f'{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}',
        f'# Source: {os.path.basename(pt)}',
        '#',
        '# Submits every generated script to SLURM, chaining steps with afterok',
        '# dependencies so each stage waits for the previous one to finish cleanly.',
        'set -euo pipefail',
        'cd "$(dirname "${BASH_SOURCE[0]}")"',
        '',
        f'echo "Submitting {len(written)} job(s) across {len(steps)} step(s)"',
        '',
    ]
    for entry in steps:
        name = entry['name']
        var = _pipeline_var(name)
        deps = _pipeline_split(entry['vals'].get('depends_on', ''))
        if deps:
            dep_expr = ':'.join(f'${{{_pipeline_var(d)}_DEP}}' for d in deps)
            dep_arg = f' --dependency=afterok:{dep_expr}'
            header = f'# ---------- {name} (after {", ".join(deps)}) ----------'
        else:
            dep_arg = ''
            header = f'# ---------- {name} ----------'
        lines.append(header)
        lines.append(f'{var}_JIDS=()')
        for script in scripts_by_step[name]:
            lines.append(f'_jid=$(sbatch --parsable{dep_arg} {script})')
            # --parsable can return "jobid;cluster" on federated clusters
            lines.append(f'_jid=${{_jid%%;*}}')
            lines.append(f'{var}_JIDS+=("$_jid")')
            lines.append(f'echo "  {script} -> $_jid"')
        lines.append(f'{var}_DEP=$(IFS=:; echo "${{{var}_JIDS[*]}}")')
        lines.append('')
    lines.append('echo')
    lines.append('echo "All jobs submitted. Track them with: squeue -u $USER"')
    lines.append('')
    submit_text = '\n'.join(lines)

    # Dry run: report exactly what would be generated, without touching the filesystem
    if dry_run:
        print(f'DRY RUN: nothing written. {len(jobs)} job script(s) plus {submit_file} '
              f'would be created in {os.path.abspath(dir)}\n')
        for entry in steps:
            deps = _pipeline_split(entry['vals'].get('depends_on', ''))
            after = f' (after {", ".join(deps)})' if deps else ''
            print(f'=== {entry["name"]}{after} ===')
            for job in [j for j in jobs if j['step'] == entry['name']]:
                print(f'  {job["file"]}  [{job["cores"]} core(s), {job["partition"]}, '
                      f'{job["time"]}, {job["mem"]} MB, module {job["python"]}, env {job["env"]}]')
                for ln in job['cmd'].splitlines():
                    print(f'      {ln}')
            print()
        print(f'=== {submit_file} ===')
        for entry in steps:
            deps = _pipeline_split(entry['vals'].get('depends_on', ''))
            n = len(scripts_by_step[entry['name']])
            after = f'afterok: {", ".join(deps)}' if deps else 'no dependency'
            print(f'  {entry["name"]:<12} {n:>3} job(s)   {after}')
        print(f'\nRe-run without --dry_run to write these files.')
        return written

    # Write every job script, then the submission driver
    for job in jobs:
        create_sh(file=os.path.join(dir, job['file']), cores=job['cores'], partition=job['partition'],
                  time=job['time'], mem=job['mem'], email=job['email'],
                  python=job['python'], env=job['env'],
                  cmd=job['cmd'], log_prefix=job['log_prefix'])

    submit_pt = os.path.join(dir, submit_file)
    with open(submit_pt, 'w') as fh:
        fh.write(submit_text)
    os.chmod(submit_pt, 0o755)

    print(f'Wrote {len(written)} job script(s) to {os.path.abspath(dir)}:')
    for entry in steps:
        print(f'  {entry["name"]}: {", ".join(scripts_by_step[entry["name"]])}')
    print(f'Wrote {submit_file}; run it on the cluster with: bash {submit_file}')

    return written

def combine(
    in_dir: str | Path,
    out_dir: str | Path,
    out_file: str,
    suffixes: Iterable[str] = [".txt", ".log", ".out", ".err"],
    recursive: bool = False,
    full_path: bool = False,
    encoding: str = "utf-8",
) -> Path:
    """
    combine(): Combine text files matching provided suffixes into a single output file, inserting a header with the original filename before each file's content.

    Parameters:
    in_dir (str): Directory to search for input files.
    out_dir (str): Directory to write the combined file.
    out_file (str): Output filename (text file).
    suffixes (Iterable[str]): Iterable of suffixes to match (e.g. [".txt", ".log", ".out", ".err"]).
    recursive (bool): If True, search subdirectories recursively.
    full_path (bool): If True, include full path in header; otherwise only filename.  
    encoding (str): Text encoding to use when reading/writing files.

    Returns: Path to the combined output file.
    """
    in_dir = Path(in_dir)
    mkdir(out_dir) # Has to be before out_path definition
    out_path = Path(out_dir) / out_file
    
    suffixes = tuple(suffixes)
    if not suffixes:
        raise ValueError("suffixes must contain at least one suffix")

    # Collect matching files
    iterator = in_dir.rglob("*") if recursive else in_dir.glob("*")
    files = [
        p for p in iterator
        if p.is_file() and p.name.endswith(suffixes)
    ]

    if not files:
        raise FileNotFoundError(
            f"No files found in {in_dir} matching suffixes={list(suffixes)}"
        )

    # Natural sort (relative path keeps directory structure ordering sensible)
    files.sort(key=lambda p: t.natural_key(str(p.relative_to(in_dir))))

    # Write combined output
    with open(out_path, "w", encoding=encoding, newline="") as out:
        for src in files:
            header_name = str(src) if full_path else src.name

            out.write(f"### SOURCE_FILE: {header_name}\n")

            with open(src, "r", encoding=encoding, errors="replace") as inp:
                shutil.copyfileobj(inp, out)

            # Ensure clean separation between files
            out.write("\n")

    return out_path

def split_R1_R2(dir: str):
    '''
    split_R1_R2(): split paired reads into new R1 and R2 subdirectories at the parent directory

    Parameters:
    dir (str): path to parent directory

    Depedencies: os, shutil, utils.mkdir()
    '''
    r1_dir = os.path.join(dir, 'R1')
    r2_dir = os.path.join(dir, 'R2')

    # Create directories if they don't exist
    mkdir(r1_dir)
    mkdir(r2_dir)

    # Move files based on naming pattern
    for fname in os.listdir(dir):
        if '_R1_' in fname:
            shutil.move(os.path.join(dir, fname), os.path.join(r1_dir, fname))
        elif '_R2_' in fname:
            shutil.move(os.path.join(dir, fname), os.path.join(r2_dir, fname))

    print(f"Moved paired reads into {r1_dir} and {r2_dir}")

def basespace(dir: str='.', suf: str='.fastq.gz', exclude: str | Iterable[str] | None='Undetermined',
              prefix_sep: str='_', dry_run: bool=False) -> list[str]:
    '''
    basespace(): reorganize an Illumina BaseSpace download into a parsable format.

    Inspects every 1st-layer folder of an Illumina BaseSpace download and, for each one that
    holds sample fastqs, runs the standard reorganization in place:
        edms io out
        edms io split_R1_R2
        cd R1; edms io in -s .fastq.gz -g prefix -p _; edms fastq comb -r
        cd ../R2; edms io in -s .fastq.gz -g prefix -p _; edms fastq comb -r

    A 1st-layer folder is only processed if it contains at least one fastq that is not an
    undetermined read; folders holding only undetermined reads (e.g., "MUZ368_600_cycle_1-*")
    or no fastqs at all (e.g., "ICA_Workflows_*") are skipped.

    Parameters:
    dir (str, optional): Path to the BaseSpace download directory (Default: '.').
    suf (str, optional): Fastq file suffix used to find & group reads (Default: '.fastq.gz').
    exclude (str | Iterable[str] | None, optional): Filename prefix, or list of prefixes, marking
        reads that don't count as samples (Default: 'Undetermined'; case-insensitive). Pass None or
        an empty list to let every fastq count.
    prefix_sep (str, optional): Delimiter splitting the sample name from the rest of the fastq
        filename (Default: '_'; e.g., MUZ350-201_S1_L001_R1_001.fastq.gz → MUZ350-201).
    dry_run (bool, optional): Print the commands that would be run for each folder without
        moving any files (Default: False).

    Returns: list of the 1st-layer folder names that were (or would be) processed.

    Dependencies: os, pathlib, tidy.natural_key(), out_subs(), split_R1_R2(), in_subs(),
        bio.fastq.comb_fastqs()
    '''
    from ..bio import fastq as fq # deferred b/c bio.fastq imports gen.io

    if not os.path.isdir(dir):
        raise ValueError(f"{dir} is not a valid directory.")

    # Accept 1 prefix or several; an empty exclude lets every fastq count as a sample
    if exclude is None:
        exclude = []
    elif isinstance(exclude, str):
        exclude = [exclude]
    exclude = [str(e) for e in exclude if str(e)]
    excluded = tuple(e.lower() for e in exclude) # matched case-insensitively
    excluded_txt = ' or '.join(f"'{e}*'" for e in exclude) # reported as the user wrote it

    subdirs = sorted((entry for entry in os.listdir(dir)
                      if os.path.isdir(os.path.join(dir, entry))), key=t.natural_key)

    processed = []
    for entry in subdirs:
        sub = os.path.join(dir, entry)

        # A run folder qualifies if it holds >=1 fastq that isn't an excluded read
        samples = [pt for pt in Path(sub).rglob(f'*{suf}')
                   if pt.is_file() and not pt.name.lower().startswith(excluded)]
        if not samples:
            print(f"\nSkipped {entry}: no '{suf}' files" +
                  (f" outside of {excluded_txt}" if excluded else ""))
            continue

        processed.append(entry)
        print(f"\nProcessing {entry} ({len(samples)} '{suf}' files)...")

        if dry_run: # Show the equivalent commands instead of moving anything
            print(f"  cd {sub}")
            print( "  edms io out")
            print( "  edms io split_R1_R2")
            for i,read in enumerate(['R1','R2']):
                print(f"  cd {read}" if i==0 else f"  cd ../{read}")
                print(f"  edms io in -s {suf} -g prefix -p '{prefix_sep}'")
                print( "  edms fastq comb -r")
            continue

        out_subs(dir=sub) # Flatten the BCLConvert/sample/report subdirectories
        split_R1_R2(dir=sub) # Separate the paired reads into R1 and R2

        for read in ['R1','R2']:
            read_dir = os.path.join(sub, read)
            in_subs(dir=read_dir, suf=suf, group_by='prefix', prefix_sep=prefix_sep) # 1 subdirectory per sample
            fq.comb_fastqs(in_dir=read_dir, out_dir=os.path.join(read_dir, 'combine_fastqs'),
                           recursive=True) # 1 combined fastq per sample (across lanes)

    if not processed:
        print(f"\nNo folders in {os.path.abspath(dir)} contained sample fastqs.")
    else:
        print(f"\n{'Would process' if dry_run else 'Processed'} {len(processed)} folder(s): {', '.join(processed)}")

    return processed

# Directory Methods
def relative_paths(root_dir: str) -> list[str]:
    ''' 
    relative_paths(): returns relative paths for all files in a directory including subfolders
    
    Parameters:
    root_dir (str): root directory path or relative path

    Dependencies: os
    '''
    relative_paths = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            # Get the relative path of the file
            relative_path = os.path.relpath(os.path.join(dirpath, filename), root_dir)
            relative_paths.append(relative_path)
    return relative_paths

def sorted_file_names(dir: str, suf: str='.csv') -> list[str]:
    '''
    sorted_file_names: returns sorted file names in a directory with the specified suffix

    dir (str): directory path or relative path
    suf (str): suffix to parse file names

    Dependencies: os
    '''
    return sorted([file for file in os.listdir(dir) if file[-len(suf):]==suf])