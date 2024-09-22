### io.py ###
# Author: Marc Zepeda
# Date: 2024-05-20

# Import packages
import pandas as pd
import os
import ast
import csv

# Input methods
''' get: Returns pandas dataframe from csv file.
        pt: file path
        dir: output directory path
    Dependencies: pandas
'''
def get(pt: str,typ='',**kwargs):
    suf = pt.split('.')[-1]
    if (typ=='csv')|(suf=='csv'): return pd.read_csv(filepath_or_buffer=pt,sep=',',**kwargs)
    elif (typ=='tsv')|(suf=='tsv'): return pd.read_csv(filepath_or_buffer=pt,sep='\t',**kwargs)
    elif (typ=='excel')|(suf=='xlsx'): return {sheet_name:pd.read_excel(pt,sheet_name,**kwargs) for sheet_name in pd.ExcelFile(pt).sheet_names}
    elif (typ=='html')|(suf=='html'): return pd.read_html(pt,**kwargs)
    elif typ=='': return pd.read_csv(filepath_or_buffer=pt,**kwargs)
    else: print("Error: Unknown typ specified")

''' get_dir: Returns python dictionary of dataframe from csv files within a directory.
        dir: output directory path
        suf: file type suffix
    Dependencies: pandas
'''
def get_dir(dir: str,suf='.csv',**kwargs):
    files = os.listdir(dir)
    csv_files = [file for file in files if file[-4:]==suf]
    csvs = dict()
    for csv_file in csv_files:
        csvs[csv_file[:-len(suf)]] = get(pt=os.path.join(dir,csv_file),**kwargs)
    return csvs

# Output methods
''' save: Save .csv file to a specified output directory from obj
        dir: output directory path
        file: file name
        obj: Dataframe, set, list, series
    Dependencies: pandas,os,csv
'''
def save(dir: str, file: str, obj, cols=[], id=False, sort=True, **kwargs):
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

''' save_dir: Save .csv files to a specified output directory from dictionary of objs
        dir: output directory path
        file_suffix: file name
        dc: dictionary of objs
    Dependencies: pandas,os,csv, save()
'''
def save_dir(dir: str, file_suffix: str, dc: dict, **kwargs):
    for key,val in dc.items(): save(dir=dir,file=key+file_suffix,obj=val,**kwargs)

# Input/Output Methods
''' excel_csvs: Exports excel file to .csv files in specified directory.
        pt: file path
        dir: output directory path
    Dependencies: pandas,os
'''
def excel_csvs(pt: str,dir='',**kwargs):
    if dir=='': dir = '.'.join(pt.split('.')[:-1]) # Get the directory where the Excel file is located
    if not os.path.exists(dir): os.mkdir(dir) # Make directory if it does not already exist
    for sheet_name in pd.ExcelFile(pt).sheet_names: # Loop through each sheet in the Excel file
        df = pd.read_excel(pd.ExcelFile(pt),sheet_name,**kwargs) # Read the sheet into a DataFrame
        df.to_csv(os.path.join(dir,f"{sheet_name}.csv"),index=False) # Save the DataFrame to a CSV file

''' df_to_dc_txt: Returns pandas DataFrame as a printed text that resembles a Python dictionary.
        df: pandas dataframe
    Dependencies: pandas
'''
def df_to_dc_txt(df: pd.DataFrame):
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

''' dc_txt_to_df: Returns a pandas DataFrame from text that resembles a Python dictionary.
        dc_txt: text that resembles a Python dictionary
        transpose: transpose dataframe (Default: True)
    Dependencies: pandas,ast
'''
def dc_txt_to_df(dc_txt: str, transpose=True):
    if transpose==True: return pd.DataFrame(ast.literal_eval(dc_txt)).T
    else: return pd.DataFrame(ast.literal_eval(dc_txt))

# Directory Methods
''' print_relative paths: Prints relative paths for all files in a directory including subfolders
        root_dir: root directory path or relative path
'''
def print_relative_paths(root_dir):
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            # Get the relative path of the file
            relative_path = os.path.relpath(os.path.join(dirpath, filename), root_dir)
            print(relative_path)