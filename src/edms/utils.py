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