### tidy.py ###
# Author: Marc Zepeda
# Date: 2024-08-03

# Import methods
import pandas as pd
import re

# Methods for dictionary containing dataframes.
def split_by(series, by=', '):
    ''' 
    split_by(): splits elements of list, set, or series by specified seperator
    
    Parameters
    series: list, set or series
    by (str, optional): seperator
    '''
    split_elements = []
    for element in series: 
        if isinstance(element, str): split_elements.extend(element.split(by))
    return split_elements

def isolate(dc: dict, col: str, get, get_col='', get_col_split_by='', want=True, exact=True):
    ''' 
    isolate(): isolate rows in dataframes based specified value(s)
    
    Parameters:
    dc (dict): dictionary
    col (str): df column name
    get: value, set, list, dictionary of dataframes
    get_col (str, optional): dataframe column name with get value
    get_col_split_by (str, optional): get value seperator
    want (bool, optional): do you want the value(s)?
    exact (bool, optional): exact value or contains (Default: exact)
    
    Dependencies: re, pandas, & split_by()
    '''
    if want==True: 
        if get is None: return {key:df[df[col].isnull()==True].reset_index(drop=True) for key,df in dc.items()}
        elif type(get)==set or type(get)==list or type(get)==pd.Series: 
            if exact==True: 
                if get_col_split_by=='': return {key:df[df[col].isin(set(get))==True].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].isin(set(split_by(get,by=get_col_split_by)))==True].reset_index(drop=True) for key,df in dc.items()}
            else: 
                if get_col_split_by=='': return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in get))==True].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in split_by(get,by=get_col_split_by)))==True].reset_index(drop=True) for key,df in dc.items()}
        elif type(get)==dict: 
            if exact==True: 
                if get_col_split_by=='': return {key:df[df[col].isin(set(get[key][get_col]))==True].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].isin(set(split_by(get[key][get_col],by=get_col_split_by)))==True].reset_index(drop=True) for key,df in dc.items()}
            else: 
                if get_col_split_by=='': return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in get[key][get_col]),case=False, na=False)==True].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in split_by(get[key][get_col],by=get_col_split_by)),case=False, na=False)==True].reset_index(drop=True) for key,df in dc.items()}
        else: return {key:df[df[col]==get].reset_index(drop=True) for key,df in dc.items()}
    else: 
        if get is None: return {key:df[df[col].isnull()==False].reset_index(drop=True) for key,df in dc.items()}
        elif type(get)==set or type(get)==list or type(get)==pd.Series: 
            if exact==True: 
                if get_col_split_by=='': return {key:df[df[col].isin(set(get))==False].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].isin(set(split_by(get,by=get_col_split_by)))==False].reset_index(drop=True) for key,df in dc.items()}
            else: 
                if get_col_split_by=='': return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in get))==False].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in split_by(get,by=get_col_split_by)))==False].reset_index(drop=True) for key,df in dc.items()}
        elif type(get)==dict: 
            if exact==True: 
                if get_col_split_by=='': return {key:df[df[col].isin(set(get[key][get_col]))==True].reset_index(drop=False) for key,df in dc.items()}
                else: return {key:df[df[col].isin(set(split_by(get[key][get_col],by=get_col_split_by)))==True].reset_index(drop=False) for key,df in dc.items()}
            else: 
                if get_col_split_by=='': return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in get[key][get_col]),case=False, na=False)==False].reset_index(drop=True) for key,df in dc.items()}
                else: return {key:df[df[col].str.contains('|'.join(re.escape(sub) for sub in split_by(get[key][get_col],by=get_col_split_by)),case=False, na=False)==False].reset_index(drop=True) for key,df in dc.items()}
        else: return {key:df[df[col]!=get].reset_index(drop=True) for key,df in dc.items()}

def modify(dc: dict, col: str, val, axis=1, **kwargs):
    ''' 
    modify(): Returns dictionary containing dataframes new or updated column with specified value(s) or function
    
    Parameters:
    dc (dict): dictionary
    col (str): new/old column name
    val: column value, list, or function (e.g., new_val=lambda df: df['AA Mutation'].split('.')[1])
    axis (int, optional): function is applied to column (1) or row (0)
    
    Dependencies: pandas
'''
    dc2=dict()
    for key,df in dc.items():
        if callable(val): df2 = df.assign(**{col: df.apply(val, axis=axis, **kwargs)})
        else: df2 = df.assign(**{col: val})
        dc2[key]=df2
    return dc2

def melt(dc: dict,id_vars,**kwargs):
    ''' 
    melt(): returns dictionary containing tidy dataframes
    
    Parameters:
    dc: dictionary of dataframes
    id_vars: metadata columns
    
    Dependencies: pandas
    '''
    dc2=dict()
    for key,df in dc.items(): dc2[key]=pd.melt(frame=df,id_vars=id_vars,**kwargs)
    return dc2

def join(dc: dict, col='key'):
    ''' 
    join(): returns a single dataframe from a dictionary of dataframes
    
    Parameters:
    dc (dict): dictionary of dataframes
    col (str, optional): name for keys column
    
    Dependencies: pandas
'''
    df = pd.DataFrame()
    for key,val in dc.items():
        val[col]=key
        df=pd.concat([df,val]).reset_index(drop=True)
    return df

def split(df: pd.DataFrame, key: str):
    ''' 
    split(): returns from a dictionary of dataframes from a single dataframe
    
    Parameters:
    df: dataframe
    key: column for spliting dataframe
    
    Dependencies: pandas
'''
    return {k:df[df[key]==k] for k in list(df[key].value_counts().keys())} 

def merge(data: pd.DataFrame, meta: pd.DataFrame, id, cols: list):
    ''' 
    merge(): adds metadata columns to data dataframe using metadata dataframe
    
    Parameters:
    data (dataframe): data dataframe
    meta (dataframe): metadata dataframe
    id: id(s) column name(s) [str: both, list: data & meta]
    cols (list): list of column names in metadata dataframe
    
    Dependencies: pandas
    '''
    if type(id)==str:
        for c in cols: 
            id_c = dict(zip(meta[id],meta[c]))
            data[c] = [id_c[i] for i in data[id]]
    elif (type(id)==list)&(len(id)==2):
        for c in cols: 
            id_c = dict(zip(meta[id[1]],meta[c]))
            data[c] = [id_c[i] for i in data[id[0]]]
    else: print("Error: id needs to be string or list of 2 strings")
    return data

# Methods for interconverting dictionaries and lists
def dc_to_ls(dc: dict,sep='.'):
    ''' 
    dc_to_ls(): convert a dictionary containing several subdictionaries into a list with all the key value relationships stored as individual values
    
    Parameters:
    dc (dict): dictionary
    sep (str, optional): seperator for subdictionaries for values in the list
    
    Dependencies: 
    '''
    ls = [] # Initialize final list
    
    def recursive_items(dc, sep='', parent_key=''): # Recursive processing submethod
        for key, value in dc.items():
            new_key = f"{parent_key}{sep}{key}" if parent_key else key # Construct the key path (optional, depending on whether you want the full key path in tuples)
            if isinstance(value, dict): recursive_items(value, sep, new_key) # Recursively process the sub-dictionary
            else: ls.append(f"{new_key}{sep}{value}") # Store the key-value pair as a tuple
    
    recursive_items(dc,sep) # Initialized recursive processing
    return ls

def ls_to_dc(ls: list, sep='.'):
    ''' 
    ls_to_dc(): convert a dictionary containing several subdictionaries into a list with all the key value relationships stored as individual values
    
    Parameters:
    ls (list): list
    sep (str, optional): seperator for subdictionaries for values in the list
    
    Dependencies:
    '''

    dc = {} # Initialize final dict

    for item in ls: # Interate through values in the list
        key, value = item.rsplit(sep, 1) # Split final key value relationship by the seperator
        parts = key.split(sep)  # Split the key by the separator
        d = dc
        for part in parts[:-1]:
            d = d.setdefault(part, {})
        d[parts[-1]] = value.strip()  # Assign the value, strip any leading/trailing whitespace

    return dc

# Dataframe methods
def reorder_cols(df: pd.DataFrame, cols: list, keep=True):
    ''' 
    reorder_cols(): returns dataframe with columns reorganized 
    
    Parameters:
    df (dataframe): pandas dataframe
    cols (list): list of column names prioritized in order
    keep (bool, optional): keep columns not listed (Default: True)

    Dependencies: pandas
    '''
    if keep==True: cols.extend([c for c in list(df.columns) if c not in cols]) # Append remaining columns
    return df[cols]