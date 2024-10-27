### io.py ###
# Author: Marc Zepeda
# Date: 2024-05-20

# Import packages
import pandas as pd
import os
import ast
import csv

# Input methods
def get(pt: str,typ='',**kwargs):
    ''' 
    get(): returns pandas dataframe from file.

    Parameters:    
    pt (str): file path
    typ (str, optional): file type (can be determined from suffix)
    
    Dependencies: pandas
    '''
    suf = pt.split('.')[-1]
    if (typ=='csv')|(suf=='csv'): return pd.read_csv(filepath_or_buffer=pt,sep=',',**kwargs)
    elif (typ=='tsv')|(suf=='tsv'): return pd.read_csv(filepath_or_buffer=pt,sep='\t',**kwargs)
    elif (typ=='excel')|(suf=='xlsx'): return {sheet_name:pd.read_excel(pt,sheet_name,**kwargs) for sheet_name in pd.ExcelFile(pt).sheet_names}
    elif (typ=='html')|(suf=='html'): return pd.read_html(pt,**kwargs)
    elif typ=='': return pd.read_csv(filepath_or_buffer=pt,**kwargs)
    else: print("Error: Unknown typ specified")

def get_dir(dir: str,suf='.csv',**kwargs):
    ''' 
    get_dir(): returns python dictionary of dataframe from files within a directory.
    
    Parameters:
    dir (str): directory path with files
    suf (str): file type suffix
    
    Dependencies: pandas
    '''
    files = os.listdir(dir)
    csv_files = [file for file in files if file[-4:]==suf]
    csvs = dict()
    for csv_file in csv_files:
        csvs[csv_file[:-len(suf)]] = get(pt=os.path.join(dir,csv_file),**kwargs)
    return csvs

# Output methods
def save(dir: str, file: str, obj, cols=[], id=False, sort=True, **kwargs):
    ''' 
    save(): save .csv file to a specified output directory from obj
    
    Parameters
    dir (str): output directory path
    file (str): file name
    obj: dataframe, series, set, or list
    cols (str, list, optional): isolate dataframe column(s)
    id (bool, optional): include dataframe index (False)
    sort (bool, optional): sort set, list, or series before saving (True)
    
    Dependencies: pandas, os, & csv
    '''
    if not os.path.exists(dir): # Make output directory if it does not exist
        os.mkdir(dir)

    if type(obj)==pd.DataFrame:
        for col in cols: # Check if each element in the list is a string
            if not isinstance(col, str):
                raise ValueError("All elements in the list must be strings.")
        if cols!=[]: obj = obj[cols]
        obj.to_csv(os.path.join(dir,file), index=id,**kwargs)
    elif type(obj)==set or type(obj)==list or type(obj)==pd.Series:
        if sort==True: obj2 = sorted(list(obj))
        else: obj2=list(obj)
        with open(os.path.join(dir,file), 'w', newline='') as csv_file:
            csv_writer = csv.writer(csv_file, dialect='excel') # Create a CSV writer object
            csv_writer.writerow(obj2) # Write each row of the list to the CSV file

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
    
    Dependencies: pandas & os
    '''
    if dir=='': dir = '.'.join(pt.split('.')[:-1]) # Get the directory where the Excel file is located
    if not os.path.exists(dir): os.mkdir(dir) # Make directory if it does not already exist
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