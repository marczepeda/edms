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
'''
# Import packages
import os
import psutil
import time
import importlib.resources
import pandas as pd

import shutil
import stat
import subprocess 
import sys
import tempfile

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
