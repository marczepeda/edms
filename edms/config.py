''' 
Module: config.py
Author: Marc Zepeda
Created: 2024-12-25
Description: Configuration

Usage:
[Supporting]
- str_dc(): Check that all dictionary components are strings

[CONFIG_FILE]
- load_config(): Load configuration from the file
- save_config(): Save configuration to the file

[Information]
- get_info(): Retrieve information based on id
- set_info(): Store a information based on id  
'''

# Import Packages
import json
import os
import edms.gen.io as io

# Supporting
def str_dc(dc):
    '''
    str_dc(): Check that all dictionary components are strings

    Parameters:
    dc: recursive object from dictionary
    '''
    if isinstance(dc,dict):
        for key,val in dc.items():
            if isinstance(key,str)==False: TypeError(f"{key} was not a string")
            str_dc(val)
    elif isinstance(dc, str): return
    else: TypeError(f"{dc} was not a string")

# CONFIG_FILE
dir = os.path.expanduser("~/.config/edms") # Make directory for configuration file
io.mkdir(dir)
CONFIG_FILE = os.path.join(dir,".config.json") # Define the path for the configuration file

def load_config():
    """
    load_config(): Load configuration from the file
    
    Dependencies: os,json
    """
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, "r") as file:
            return json.load(file)
    return {}  # Return empty dict if config file doesn't exist

def save_config(config):
    """
    save_config(): Save configuration to the file
    
    Dependencies: json
    """
    with open(CONFIG_FILE, "w") as file:
        json.dump(config, file, indent=4)

# Information
def get_info(id: str):
    """
    get_info(): Retrieve information based on id
    
    Parameters:
    id (str): identifier from configuration file

    Dependencies: load_config()
    """
    # Get information
    config = load_config()
    info = config.get(id, None)

    # Report on information status
    if info: print(f"Got {id}: {info}")
    else: print(f"No information found for {id}")
    return info

def set_info(id: str, info: str | dict):
    """
    set_info(): Store a information based on id
    
    Parameters:
    id (str): identifier for configuration file
    info (str | dict): information for configuration file

    Dependencies: load_config(),save_config()
    """
    # Check that info dictionaries only contain strings
    if isinstance(info,dict):
        for key,val in info.items():
            if isinstance(key,str)==False: TypeError(f"{key} was not a string")
            str_dc(val)

    # (Over)write information to configuration file
    config = load_config()
    config[id] = info
    save_config(config)
    print(f"Set {id}: {info}")