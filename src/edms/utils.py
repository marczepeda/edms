'''
Module: utils.py
Author: Marc Zepeda
Created: 2025-05-05
Description: Ulity

Usage:
[Computation]
- memory(): report current memory
- timer(): report elapsed time
- memory_timer(): report current memory & elapsed time

[Package .csv files]
- load_resource_csv(): load .csv file from resources

[Execute bash script]    
- run_bundled_script(): Run a bundled script from the package script resources.

[Parsing Python literals]
- df_try_parse(): apply try_parse() to dataframe columns and return dataframe
- recursive_json_decode(): recursively decode JSON strings in a nested structure

[Make directories]
- mkdir(): make directory if it does not exist (including parent directories)
- check_outpath(): resolve an output file path and ensure its directory exists

[Filter dataframes]
- apply_filters(): apply pandas query strings and drop_duplicates to a dataframe, reporting row counts

[Supporting argument methods]
- parse_tuple_int(arg): Parse a string argument into a tuple of integers
- parse_tuple_float(arg): Parse a string argument into a tuple of floats
'''
# Import packages
import os
import psutil
import time
import importlib.resources
import pandas as pd
import argparse
import ast
import shutil
import stat
import subprocess 
import sys
import tempfile
from pathlib import Path

# Computation
process = psutil.Process(os.getpid())

def memory(task: str='unspecified') -> tuple:
    '''
    memory(): report current memory

    Parameters:
    task (str, optional): reporting memory for... (Default: unspecified)

    Dependencies: psutil, os
    '''
    mem = process.memory_info().rss / 1e6  # MB
    print(f"{task}:\tMemory: {mem:.2f} MB")
    return task,mem

def timer(task: str='unspecified', reset: bool=False):
    '''
    timer(): report elapsed time

    Parameters:
    task (str, optional): reporting time for... (Default: unspecified)
    reset (bool, optional): reset timer (Default: False)
    
    Dependencies: time
    '''
    now = time.perf_counter() # s
    
    if reset and hasattr(timer, "last_time"): # Reset/start timer
        delattr(timer, "last_time")
        timer.last_time = now
        return
    elif hasattr(timer, "last_time"):
        elapsed = now - timer.last_time # s
        print(f"{task}:\tTime: {elapsed:.2f} s")
        timer.last_time = now
        return task,elapsed
    else: # Start timer
        timer.last_time = now

def memory_timer(task: str='unspecified', reset: bool=False):
    '''
    memory_timer(): report current memory & elapsed time

    Parameters:
    task (str, optional): reporting memory/time for... (Default: unspecified)
    reset (bool, optional): reset timer (Default: False)

    Dependencies: psutil, os, time
    '''
    mem = process.memory_info().rss / 1e6  # MB
    now = time.perf_counter() # s

    if reset and hasattr(timer, "last_time"): # Reset/start timer
        delattr(timer, "last_time")
        timer.last_time = now
        return
    elif hasattr(timer, "last_time"):
        elapsed = now - timer.last_time # s
        print(f"{task}:\tMemory: {mem:.2f} MB\tTime: {elapsed:.2f} s")
        timer.last_time = now
        return task,mem,elapsed
    else: # Start timer
        timer.last_time = now

# Package .csv files
def load_resource_csv(filename: str):
    '''
    laod_resource_csv(): load .csv file from resources
    
    Parameters:
    filename (str): name of .csv file in resources folder

    Dependencies: importlib.resources, pandas
    '''
    with importlib.resources.files("edms.resources").joinpath(filename).open("r", encoding="utf-8") as f:
        return pd.read_csv(f)
    
# Execute bash script
def run_bundled_script(package: str="edms.scripts", relpath: str='autocomplete.sh', args: list[str]=sys.argv[1:]) -> int:
    '''
    run_bundled_script(): Run a bundled script from the package script resources.
    
    Parameters:
    package (str): The package where the resource is located (Default: 'edms.scripts').
    relpath (str): The relative path to the script within the package (Default: 'autocomplete.sh').
    args (list[str]): List of arguments to pass to the script (Default: sys.argv[1:]).
    '''
    # Locate the resource inside the installed wheel
    src = importlib.resources.files(package).joinpath(relpath)
    if not src.is_file():
        print(f"Error: resource {relpath} not found in {package}", file=sys.stderr)
        return 1

    # Copy to a temp file to ensure exec perms on all platforms
    with tempfile.TemporaryDirectory() as td:
        dst = os.path.join(td, os.path.basename(relpath))
        shutil.copyfile(src, dst)
        os.chmod(dst, os.stat(dst).st_mode | stat.S_IXUSR)
        # Run with inherited stdio
        proc = subprocess.run([dst, *args], check=False)
        return proc.returncode

# Parsing Python literals
def try_parse(value):
    """
    try_parse(): try to parse a string as a Python literal; if not possible, return the original value
    
    Parameters:
    value: value of any type (looking for strings)

    Dependencies: ast,recusive_parse()
    """
    if isinstance(value, str):  # Only attempt parsing for strings
        try:
            parsed_value = ast.literal_eval(value)
            # If parsed_value is a dictionary, list, etc., recursively evaluate its contents
            if isinstance(parsed_value, (dict, list, set, tuple)):
                return recursive_parse(parsed_value)
            return parsed_value
        except (ValueError, SyntaxError):
            # Return the value as-is if it can't be parsed
            return value
    return value

def recursive_parse(data):
    """
    recursive_parse(): recursively parse nested structures
    
    Parameters:
    data: data of any type (looking for dict, list, set, or tuple)

    Dependencies: ast,try_parse()
    """
    if isinstance(data, dict):
        return {k: try_parse(v) for k, v in data.items()}
    elif isinstance(data, list):
        return [try_parse(item) for item in data]
    elif isinstance(data, set):
        return {try_parse(item) for item in data}
    elif isinstance(data, tuple):
        return tuple(try_parse(item) for item in data)
    return data  # Return the data as-is if it doesn't match any known structure

# Make directories
def mkdir(dir: str | Path, sep: str='/'):
    '''
    mkdir(): make directory if it does not exist (including parent directories)

    Parameters:
    dir (str): directory path
    sep (str): seperator directory path

    Dependencies: os
    '''
    dirs = os.path.abspath(str(dir)).split(sep)
    for i in range(len(dirs)):
        check_dir = sep.join(dirs[:i+1])
        if (os.path.exists(check_dir)==False)&(i!=0):
            os.mkdir(check_dir)
            print(f'Created {check_dir}')

def check_outpath(file: str | Path | None) -> str | None:
    '''
    check_outpath(): resolve an output file path and ensure its directory exists

    Parameters:
    file (str | Path | None): output file path (if None, no output will be saved)

    Dependencies: os & mkdir()
    '''
    if file is None or isinstance(file, (str, Path))==False:
        return None

    file = os.path.abspath(str(file))
    mkdir(os.path.dirname(file))  # Ensure output directory exists

    return file

# Filter dataframes
def apply_filters(df: pd.DataFrame, query: list[str] | str = None, drop_duplicates: list[str] | str = None,
                  label: str = 'rows', allow_empty: bool = False) -> pd.DataFrame:
    '''
    apply_filters(): apply pandas query strings and drop_duplicates to a dataframe, reporting row counts

    Queries are applied in the order given with pandas\' python engine, so string methods
    (i.e., "Reference_sequence.str.contains(\'ACGTTCAAG\')") and chained comparisons
    (i.e., "difference_count >= 1 and difference_count <= 3") both work. drop_duplicates is
    applied last, keeping the first occurrence.

    Parameters:
    df (dataframe): dataframe to filter
    query (list[str] | str, optional): pandas query string(s) applied in order (Default: None)
    drop_duplicates (list[str] | str, optional): column name(s) to drop duplicates on (Default: None)
    label (str, optional): noun used in the printed row counts (Default: \'rows\')
    allow_empty (bool, optional): return an empty dataframe instead of raising (Default: False)

    Dependencies: pandas
    '''
    if query is None: query = []
    elif isinstance(query, str): query = [query]
    if drop_duplicates is None: drop_duplicates = []
    elif isinstance(drop_duplicates, str): drop_duplicates = [drop_duplicates]

    if not query and not drop_duplicates: # Nothing to do; leave the dataframe untouched
        return df

    print(f'Filtering: {len(df)} {label}')
    for q in query:
        before = len(df)
        try:
            df = df.query(q, engine='python')
        except Exception as e:
            raise ValueError(f"Filter query failed: {q!r}\n"
                             f"{type(e).__name__}: {e}\n"
                             f"Available columns: {list(df.columns)}") from e
        print(f'  --query {q!r}: {before} -> {len(df)} {label}')
        if len(df) == 0 and not allow_empty:
            raise ValueError(f"Filter query {q!r} removed every row; nothing left to work with.")

    if drop_duplicates:
        missing = [col for col in drop_duplicates if col not in df.columns]
        if missing:
            raise ValueError(f"drop_duplicates column(s) {missing} not found. "
                             f"Available columns: {list(df.columns)}")
        before = len(df)
        df = df.drop_duplicates(subset=drop_duplicates)
        print(f'  --drop_duplicates {drop_duplicates}: {before} -> {len(df)} {label}')
        if len(df) == 0 and not allow_empty:
            raise ValueError(f"drop_duplicates on {drop_duplicates} removed every row.")

    return df.reset_index(drop=True)

# Supporting argument methods
def parse_tuple_int(arg):
    '''
    parse_tuple_int(arg): Parse a string argument into a tuple of integers
    '''
    try:
        return tuple(map(int, arg.split(',')))
    except:
        raise argparse.ArgumentTypeError(f"{arg} must be formatted as 'int,int,...'")
    
def parse_tuple_float(arg):
    '''
    parse_tuple_float(arg): Parse a string argument into a tuple of floats
    '''
    try:
        return tuple(map(float, arg.split(',')))
    except:
        raise argparse.ArgumentTypeError(f"{arg} must be formatted as 'float,float,...'")