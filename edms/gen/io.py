### io.py ###
# Author: Marc Zepeda
# Date: 2024-05-20

# Import packages
import pandas as pd
import os
import ast
import csv
import subprocess

# Parsing Python literals methods
def try_parse(value):
    """
    try_parse(): try to parse a string as a Python literal; if not possible, return the original value.
    
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
    recursive_parse(): recursively parse nested structures.
    
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

def df_try_parse(df: pd.DataFrame):
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

# Input methods
def get(pt: str,literal_eval=False,**kwargs):
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
        else: pd.read_csv(filepath_or_buffer=pt,sep=',',**kwargs)
    elif suf=='tsv': 
        if literal_eval: return df_try_parse(pd.read_csv(filepath_or_buffer=pt,sep='\t',**kwargs))
        else: pd.read_csv(filepath_or_buffer=pt,sep='\t',**kwargs)
    elif suf=='xlsx': 
        if literal_eval: return {sheet_name:df_try_parse(pd.read_excel(pt,sheet_name,**kwargs))
                                 for sheet_name in pd.ExcelFile(pt).sheet_names}
        else: {sheet_name:pd.read_excel(pt,sheet_name,**kwargs)
                          for sheet_name in pd.ExcelFile(pt).sheet_names}
    elif suf=='html': 
        if literal_eval: return df_try_parse(pd.read_html(pt,**kwargs))
        else: pd.read_html(pt,**kwargs)
    else: 
        if literal_eval: return df_try_parse(pd.read_csv(filepath_or_buffer=pt,**kwargs))
        else: pd.read_csv(filepath_or_buffer=pt,**kwargs)
    
def get_dir(dir: str,suf='.csv',literal_eval=False,**kwargs):
    ''' 
    get_dir(): returns python dictionary of dataframe from files within a directory.
    
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
    return {file[:-len(suf)]:get(os.path.join(dir,file),literal_eval,**kwargs) for file in files}

# Output methods
def mkdir(dir: str, sep='/'):
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


def save(dir: str, file: str, obj, cols=[], id=False, sort=True, **kwargs):
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

def save_dir(dir: str, file_suffix: str, dc: dict, **kwargs):
    ''' 
    save_dir(): save .csv files to a specified output directory from dictionary of objs
    
    Parameters:
    dir (str): output directory path
    file_suffix (str): file name suffix
    dc (dict): dictionary of objects (files)

    Dependencies: pandas, os, csv, & save()
    '''
    for key,val in dc.items(): save(dir=dir,file=key+file_suffix,obj=val,**kwargs)

# Input/Output Methods
def excel_csvs(pt: str,dir='',**kwargs):
    ''' 
    excel_csvs(): exports excel file to .csv files in specified directory.    
    
    Parameters:
    pt (str): excel file path
    dir (str, optional): output directory path (same directory as excel file)
    
    Dependencies: pandas, os, & mkdir
    '''
    if dir=='': dir = '.'.join(pt.split('.')[:-1]) # Get the directory where the Excel file is located
    mkdir(dir) # Make output directory if it does not exist
    for sheet_name in pd.ExcelFile(pt).sheet_names: # Loop through each sheet in the Excel file
        df = pd.read_excel(pd.ExcelFile(pt),sheet_name,**kwargs) # Read the sheet into a DataFrame
        df.to_csv(os.path.join(dir,f"{sheet_name}.csv"),index=False) # Save the DataFrame to a CSV file

def df_to_dc_txt(df: pd.DataFrame):
    ''' 
    df_to_dc_txt(): returns pandas DataFrame as a printed text that resembles a Python dictionary.
    
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

def dc_txt_to_df(dc_txt: str, transpose=True):
    ''' 
    dc_txt_to_df(): returns a pandas DataFrame from text that resembles a Python dictionary.
    
    Parameters:
    dc_txt (str): text that resembles a Python dictionary
    transpose (bool, optional): transpose dataframe (True)
    
    Dependencies: pandas & ast
    '''
    if transpose==True: return pd.DataFrame(ast.literal_eval(dc_txt)).T
    else: return pd.DataFrame(ast.literal_eval(dc_txt))

def echo(message: str):
    """
    echo(): write message to slurm output file using echo

    Parameters:
    message (str): text that will be written to slurm output file.

    Dependencies: subprocess
    """
    try:
        # Use subprocess to run `echo` and append the message to stdout
        subprocess.run(['echo', message], check=True)
    except Exception as e:
        print(f"Failed to write message: {e}")

# Directory Methods
def print_relative_paths(root_dir: str):
    ''' 
    print_relative_paths(): prints relative paths for all files in a directory including subfolders
    
    Parameters:
    root_dir (str): root directory path or relative path

    Dependencies: os
    '''
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            # Get the relative path of the file
            relative_path = os.path.relpath(os.path.join(dirpath, filename), root_dir)
            print(relative_path)