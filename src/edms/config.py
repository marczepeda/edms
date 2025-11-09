''' 
Module: config.py
Author: Marc Zepeda
Created: 2024-12-25
Description: Configuration

Usage:
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
from .utils import try_parse, mkdir

# CONFIG_FILE
dir = os.path.expanduser("~/.config/edms") # Make directory for configuration file
mkdir(dir)
CONFIG_FILE = os.path.join(dir,".config.json") # Define the path for the configuration file

def load_config():
    """
    load_config(): Load configuration from the file
    
    Dependencies: os, json
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
def get_info(id: str=None):
    """
    get_info(): Retrieve information based on id
    
    Parameters:
    id (str, optional): identifier from configuration file (Default: None)

    Dependencies: load_config()
    """
    # Get information
    config = load_config()
    info = config.get(id, None)

    # Report on information status
    if info: print(f"Got {id}: {info}")
    else: print(f"Configuration File:\n{json.dumps(config, indent=4)}")
    return info

def set_info(id: str, info: str | dict | list | set | tuple):
    """
    set_info(): Store a information based on id
    
    Parameters:
    id (str): identifier for configuration file
    info (str | dict | list | set | tuple): information for configuration file

    Dependencies: try_parse(),load_config(),save_config()
    """
    # Parse info
    info = try_parse(info)

    # (Over)write information to configuration file
    config = load_config()
    config[id] = info
    save_config(config)
    print(f"Set {id}: {info}")

def del_info(id: str):
    """
    del_info(): Delete information based on id
    
    Parameters:
    id (str): identifier for configuration file

    Dependencies: load_config(),save_config()
    """
    # (Over)write information to configuration file
    config = load_config()
    config.pop(id, None)
    save_config(config)
    print(f"Deleted {id}")