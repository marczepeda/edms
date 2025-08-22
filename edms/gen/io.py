'''
Module: io.py
Author: Marc Zepeda
Created: 2024-05-20
Description: Input/Output

Usage:
[Parsing Python literals]
- try_parse(): try to parse a string as a Python literal; if not possible, return the original value
- recursive_parse(): recursively parse nested structures
- df_try_parse(): apply try_parse() to dataframe columns and return dataframe
- recursive_json_decode(): recursively decode JSON strings in a nested structure

[Input]
- get(): returns pandas dataframe from a file
- get_dir(): returns python dictionary of dataframe from files within a directory

[Output]
- mkdir(): make directory if it does not exist (including parent directories)
- save(): save .csv file to a specified output directory from obj
- save_dir(): save .csv files to a specified output directory from dictionary of objs

[Input/Output]
- excel_csvs(): exports excel file to .csv files in specified directory
- df_to_dc_txt(): returns pandas DataFrame as a printed text that resembles a Python dictionary
- dc_txt_to_df(): returns a pandas DataFrame from text that resembles a Python dictionary
- in_subs(): moves all files with a given suffix into subfolders named after the files (excluding the suffix).
- out_subs(): delete subdirectories and move their files to the parent directory
- create_sh(): creates a shell script with SLURM job submission parameters for Harvard FASRC cluster.
- split_R1_R2(): split paired reads into new R1 and R2 subdirectories at the parent directory

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
import json

from .. import config

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

def df_try_parse(df: pd.DataFrame) -> pd.DataFrame:
    '''
    df_try_parse(): apply try_parse() to dataframe columns and return dataframe

    Parameters: 
    df (dataframe): dataframe with columns to try_parse()

    Dependencies: try_parse()
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
    
    Dependencies: pandas,ast,try_parse(),recursive_parse(),df_try_parse()
    '''
    suf = pt.split('.')[-1]
    if suf=='csv': 
        if literal_eval: return df_try_parse(pd.read_csv(filepath_or_buffer=pt,sep=',',**kwargs))
        else: return pd.read_csv(filepath_or_buffer=pt,sep=',',**kwargs)
    elif suf=='tsv': 
        if literal_eval: return df_try_parse(pd.read_csv(filepath_or_buffer=pt,sep='\t',**kwargs))
        else: return pd.read_csv(filepath_or_buffer=pt,sep='\t',**kwargs)
    elif suf=='xlsx': 
        if literal_eval: 
            dc = {sheet_name:df_try_parse(pd.read_excel(pt,sheet_name,**kwargs))
                                 for sheet_name in pd.ExcelFile(pt).sheet_names}
        else: 
            dc = {sheet_name:pd.read_excel(pt,sheet_name,**kwargs)
                                 for sheet_name in pd.ExcelFile(pt).sheet_names}
        print(f"Excel file: {pt}\nKeys: {', '.join([key for key in dc.keys()])}")
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
    return dc

# Output
def mkdir(dir: str, sep: str='/'):
    '''
    mkdir(): make directory if it does not exist (including parent directories)

    Parameters:
    dir (str): directory path
    sep (str): seperator directory path

    Dependencies: os
    '''
    dirs = dir.split(sep)
    for i in range(len(dirs)):
        check_dir = sep.join(dirs[:i+1])
        if (os.path.exists(check_dir)==False)&(i!=0):
            os.mkdir(check_dir)
            print(f'Created {check_dir}')


def save(dir: str, file: str, obj, cols: list=[], id: bool=False, sort: bool=True, **kwargs):
    ''' 
    save(): save .csv file to a specified output directory from obj
    
    Parameters:
    dir (str): output directory path
    file (str): file name
    obj: dataframe, series, set, or list
    cols (str, list, optional): isolate dataframe column(s)
    id (bool, optional): include dataframe index (False)
    sort (bool, optional): sort set, list, or series before saving (True)
    
    Dependencies: pandas, os, csv & mkdir()
    '''
    mkdir(dir) # Make output directory if it does not exist

    if type(obj)==pd.DataFrame:
        for col in cols: # Check if each element in the list is a string
            if not isinstance(col, str):
                raise ValueError("All elements in the list must be strings.")
        if cols!=[]: obj = obj[cols]
        if file.split('.')[-1]=='tsv': obj.to_csv(os.path.join(dir,file), index=id, sep='\t', **kwargs)
        elif file.split('.')[-1]=='xlsx': 
            with pd.ExcelWriter(os.path.join(dir,file)) as writer: 
                obj.to_excel(writer,sheet_name='.'.join(file.split('.')[:-1]),index=id) # Dataframe per sheet
        else: obj.to_csv(os.path.join(dir,file), index=id, **kwargs)
    elif type(obj)==set or type(obj)==list or type(obj)==pd.Series:
        if sort==True: obj2 = sorted(list(obj))
        else: obj2=list(obj)
        with open(os.path.join(dir,file), 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file, dialect='excel') # Create a CSV writer object
            csv_writer.writerow(obj2) # Write each row of the list to the CSV file
    elif (type(obj)==dict)&(file.split('.')[-1]=='xlsx'):
        for col in cols: # Check if each element in the list is a string
            if not isinstance(col, str):
                raise ValueError("All elements in the list must be strings.")
        with pd.ExcelWriter(os.path.join(dir,file)) as writer:
            if cols!=[]: obj = obj[cols]
            for key,df in obj.items(): 
                if cols!=[]: df = df[cols]
                df.to_excel(writer,sheet_name=key,index=id) # Dataframe per sheet
    else: raise ValueError(f'save() does not work for {type(obj)} objects with {file.split(".")[-1]} files.')

def save_dir(dir: str, dc: dict, suf: str='.csv', **kwargs):
    ''' 
    save_dir(): save .csv files to a specified output directory from dictionary of objs
    
    Parameters:
    dir (str): output directory path
    dc (dict): dictionary of objects (files)
    suf (str, optional): file name suffix (Default: .csv)

    Dependencies: pandas, os, csv, & save()
    '''
    for key,val in dc.items(): save(dir=dir,file=str(key)+suf,obj=val,**kwargs)

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

def in_subs(dir: str, suf: str): 
    '''
    in_subs: moves all files with a given suffix into subdirectory named after the files (excluding the suffix).

    Parameters:
    dir (str): Path to the directory containing the files.
    suf (str): File suffix (e.g., '.txt', '.csv') to filter files.

    Dependences: os, shutil
    '''
    if not os.path.isdir(dir):
        raise ValueError(f"{dir} is not a valid directory.")

    for filename in os.listdir(dir):
        file_path = os.path.join(dir, filename)
        
        if os.path.isfile(file_path) and filename.endswith(suf):
            base_name = filename[:-len(suf)]
            subdir = os.path.join(dir, base_name)

            # Create subfolder if it doesn't exist
            os.makedirs(subdir, exist_ok=True) 

            # Move the file into the subfolder
            new_path = os.path.join(subdir, filename)
            shutil.move(file_path, new_path)

def out_subs(dir: str):
    """
    out_subs(): Delete subdirectories and move their files to the parent directory.

    Parameters:
    dir (str): Path to the directory containing the files.
    """
    if not os.path.isdir(dir):
        raise ValueError(f"{dir} is not a valid directory.")

    for root, dirs, files in os.walk(dir, topdown=False):
        for file in files:
            file_path = os.path.join(root, file)
            target_path = os.path.join(dir, file)

            # Resolve name conflicts by appending a counter
            base, ext = os.path.splitext(file)
            counter = 1
            while os.path.exists(target_path):
                target_path = os.path.join(dir, f"{base}_{counter}{ext}")
                counter += 1

            shutil.move(file_path, target_path)

        # Remove empty directories
        for dir in dirs:
            dir_path = os.path.join(root, dir)
            if os.path.isdir(dir_path) and not os.listdir(dir_path):
                os.rmdir(dir_path)

def create_sh(dir: str, file: str, 
              cores: int = 1, partition: str='serial_requeue', time: str = '0-00:10', mem: int = 1000, email: str=None,
              python: str = 'python/3.12.5-fasrc01', env: str = 'edms'):
    '''
    create_sh(): creates a shell script with SLURM job submission parameters for Harvard FASRC cluster.

    Parameters:
    dir (str): The directory where the shell script will be created.
    file (str): The name of the shell script file to create (e.g., "script.sh").
    cores (int, optional): The number of cores to request for the job (Default: 1).
    partition (str, optional): The partition to submit the job to (Default: 'serial_requeue').
    time (str, optional): The maximum runtime for the job in D-HH:MM format (Default: '0-00:10').
    mem (int, optional): The amount of memory to request for the job in MB (Default: 1000).
    email (str, optional): The email address (Default: None = get_info('Harvard')).
    python (str, optional): The python module to load (Default: 'python/3.12.5-fasrc01').
    env (str, optional): The conda environment to activate (Default: 'edms').

    Dependencies: os, datetime
    '''
    # Check if the directory exists
    mkdir(dir)

    # Check if the file is valid
    if not file.endswith('.sh'):
        raise ValueError("File must end with '.sh'")
    
    # Get email from config if not provided
    if email is None:
        try:
            email = config.get_info('Harvard')['email']
        except Exception as e:
            print(f"Error retrieving email from config: {e}")
            return

    # Create .sh file
    try:
        with open(os.path.join(dir,file), 'w') as file_obj:
            file_obj.write(f'''#!/bin/bash
# 
#SBATCH -n {cores} \t# Number of cores
#SBATCH -N 1 \t# Ensure that all cores are on one machine
#SBATCH -t {time} \t# Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p {partition} \t# Partition to submit to
#SBATCH --mem={mem} \t# Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o {datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_{".".join(file.split(".")[:-1])}_%j.out \t# File to which STDOUT will be written, %j inserts jobid
#SBATCH -e {datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}_{".".join(file.split(".")[:-1])}_%j.err \t# File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-type=ALL \t# Email     
#SBATCH --mail-user={email} \t# Email

echo -e "File: {file}\\nTime: {time}\\nMemory: {mem/1000:.3f} GB" \t# Print job parameters
module load {python} \t# Load python module
mamba activate {env} \t# Activate conda environment
export PYTHONUNBUFFERED=1 \t# Ensure prints from python script are written to .out file
python {".".join(file.split(".")[:-1])}.py \t# Run python script
''')
    except Exception as e:
        print(f"An error occurred while creating the shell script: {e}")

def split_R1_R2(dir: str):
    '''
    split_R1_R2(): split paired reads into new R1 and R2 subdirectories at the parent directory

    Parameters:
    dir (str): path to parent directory

    Depedencies: os, shutil, mkdir()
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