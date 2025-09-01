'''
Module: tidy.py
Author: Marc Zepeda
Created: 2024-08-03
Description:

Usage:
[Dataframe methods]
- reorder_cols(): returns dataframe with columns reorganized 
- zip_cols(): returns zip(dataframe[cols[0]],dataframe[cols[1]],...) for tuple loops
- missing_cols(): returns values from a list if they are not dataframe columns
- merge(): adds metadata columns to data dataframe using metadata dataframe
- shared_group(): group rows with shared values in a column, consolidate unique values from other columns into lists, & append to the original dataframe.
- unique_tuples(): returns a list of unique tuples of value(s) from a dataframe column(s) in order

[Dictionary methods]
- filter_kwargs(): filter **kwargs by specified keywords
- comb_dcs(): Combine a list of dictionaries with the same keys into a single dictionary where keys map to a list of values.

[Methods for dictionary containing dataframes]
- split_by(): splits elements of list, set, or series by specified seperator
- isolate(): isolate rows in dataframes based specified value(s)
- modify(): Returns dictionary containing dataframes new or updated column with specified value(s) or function
- join(): returns a single dataframe from a dictionary of dataframes
- split(): returns from a dictionary of dataframes from a single dataframe

[Methods for interconverting dictionaries and lists]
- dc_to_ls(): convert a dictionary containing several subdictionaries into a list with all the key value relationships stored as individual values
- ls_to_dc(): convert a list with all the key value relationships stored as individual values into a dictionary containing several subdictionaries

[String methods]
- find_all(): Find all indexes of a substring in a string
- color(): return colored text
'''
# Import packages
import pandas as pd
import re

# Dataframe methods
def reorder_cols(df: pd.DataFrame, cols: list, keep: bool=True) -> pd.DataFrame:
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

def zip_cols(df: pd.DataFrame, cols: list) -> zip:
    ''' 
    zip_cols(): returns zip(dataframe[cols[0]],dataframe[cols[1]],...) for tuple loops
    
    Parameters:
    df (dataframe): pandas dataframe
    cols (list): list of column names

    Dependencies: pandas
    ''' 
    return zip(*[df[col] for col in cols])

def missing_cols(df: pd.DataFrame, cols: list) -> list:
    '''
    missing_cols(): returns values from a list if they are not dataframe columns

    Parameters:
    df (dataframe): check this dataframe for col
    cols (list): column name to be checked

    Dependencies: pandas
    '''
    return [col for col in cols if col not in df.columns]

def merge(data: pd.DataFrame, meta: pd.DataFrame, id, cols: list) -> pd.DataFrame:
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

def shared_group(df: pd.DataFrame, shared: str, group: list | str, suffixes = ('','_list')) -> pd.DataFrame:
    """
    shared_group(): group rows with shared values in a column, consolidate unique values from other columns into lists, & append to the original dataframe.

    Parameters:
    df (dataframe): Input pandas DataFrame.
    shared (str): name of column to group by shared values.
    group (list | str): list of other columns that will have unique values grouped that share common values.
    suffixes (tuple, optional): Tuple of suffixes to apply to the columns in the output DataFrame (Default is ('','_list')).
    """
    return pd.merge(left=df, 
                    right=df.groupby(shared)[group].agg(lambda x: list(pd.unique(x))).reset_index(), 
                    on=shared, 
                    how='left', 
                    suffixes=suffixes)

def unique_tuples(df: pd.DataFrame, cols: list | str) -> list:
    ''' 
    unique_tuples(): returns a list of unique tuples of value(s) from a dataframe column(s) in order
    
    Parameters:
    df (dataframe): pandas dataframe
    cols (list | str): column name(s) to get unique pairs from
    
    Dependencies: pandas
    '''
    if type(cols)==str: cols=[cols]
    return list(df[cols].drop_duplicates().itertuples(index=False, name=None))

def vcs_ordered(df: pd.DataFrame, cols: list | str) -> pd.Series:
    ''' 
    vcs_ordered(): returns dataframe.value_counts() in order
    
    Parameters:
    df (dataframe): pandas dataframe
    cols (list | str): column name(s) to get unique values from
    
    Dependencies: pandas
    '''
    if type(cols)==str: cols=[cols]
    vcs = df[cols].value_counts() # Get value counts
    order = unique_tuples(df=df, cols=cols) # Get unique tuples
    return pd.Series([vcs[pair] for pair in order], index=order)

# Dictionary methods
def filter_kwargs(keywords: list, **kwargs) -> dict:
    '''
    filter_kwargs(): filter **kwargs by specified keywords

    Parameters:
    keywords (list): list of keywords to retain
    **kwargs: keyword = arguments to be filtered
    '''
    return {kw:arg for kw,arg in kwargs.items() if kw in keywords and arg is not None}

def comb_dcs(ls: list) -> dict:
    """
    comb_dcs(): Combine a list of dictionaries with the same keys into a single dictionary where keys map to a list of values.

    Parameters:
    ls (list): A list of dictionaries with the same keys.
    
    """
    combined = {}
    for dc in ls:
        for key, value in dc.items():
            combined.setdefault(key, []).append(value)
    return combined

# Methods for dictionary containing dataframes
def split_by(series: list | set | pd.Series, by: str=', ') -> list:
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

def isolate(dc: dict, col: str, get, get_col: str='', get_col_split_by: str='', want: bool=True, exact: bool=True) -> dict[pd.DataFrame]:
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

def modify(dc: dict, col: str, val, axis: int=1, **kwargs) -> dict:
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

def melt(dc: dict, id_vars:list=None, value_vars:list=None, **kwargs) -> dict[pd.DataFrame]:
    ''' 
    melt(): returns dictionary containing tidy dataframes
    
    Parameters:
    dc (dict): dictionary of dataframes
    id_vars (list, optional 1): metadata columns
    value_vars (list, optional 2): data columns
    
    Dependencies: pandas, missing_cols()
    '''
    # Check for id_vars and/or value_vars
    if id_vars is not None and value_vars is not None:
        dc2=dict()
        for key,df in dc.items(): 
            missing_id_vars = missing_cols(df=df,cols=id_vars)
            missing_value_vars = missing_cols(df=df,cols=value_vars)
            if len(missing_id_vars)>0: raise TypeError(f'Error: missing {missing_id_vars} in {key} dataframe columns.')
            if len(missing_value_vars)>0: raise TypeError(f'Error: missing {missing_value_vars} in {key} dataframe columns.')
            dc2[key]=pd.melt(frame=df,id_vars=id_vars,value_vars=value_vars,**kwargs)
    elif id_vars is not None:
        dc2=dict()
        for key,df in dc.items():
            missing_id_vars = missing_cols(df=df,cols=id_vars)
            if len(missing_id_vars)>0: raise TypeError(f'Error: missing {missing_id_vars} in {key} dataframe columns.')
            value_vars = list(df.columns)
            for item in id_vars: value_vars.remove(item)
            dc2[key]=pd.melt(frame=df,id_vars=id_vars,value_vars=value_vars,**kwargs)
    elif value_vars is not None:
        dc2=dict()
        for key,df in dc.items(): 
            missing_value_vars = missing_cols(df=df,cols=value_vars)
            if len(missing_value_vars)>0: raise TypeError(f'Error: missing {missing_value_vars} in {key} dataframe columns.')
            id_vars = list(df.columns)
            for item in value_vars: id_vars.remove(item)
            dc2[key]=pd.melt(frame=df,id_vars=id_vars,value_vars=value_vars,**kwargs)
    else: raise TypeError('Error: specify id_vars and/or value_vars.')
    
    return dc2

def join(dc: dict, col: str='key') -> pd.DataFrame:
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

def split(df: pd.DataFrame, key: str) -> dict[pd.DataFrame]:
    ''' 
    split(): returns from a dictionary of dataframes from a single dataframe
    
    Parameters:
    df (dataframe): dataframe
    key (str): column for spliting dataframe
    
    Dependencies: pandas
    '''
    return {k:df[df[key]==k] for k in list(df[key].value_counts().keys())} 

# Methods for interconverting dictionaries and lists
def dc_to_ls(dc: dict, sep: str='.') -> list:
    ''' 
    dc_to_ls(): convert a dictionary containing several subdictionaries into a list with all the key value relationships stored as individual values
    
    Parameters:
    dc (dict): dictionary
    sep (str, optional): seperator for subdictionaries for values in the list
    '''
    ls = [] # Initialize final list
    
    def recursive_items(dc, sep='', parent_key=''): # Recursive processing submethod
        for key, value in dc.items():
            new_key = f"{parent_key}{sep}{key}" if parent_key else key # Construct the key path (optional, depending on whether you want the full key path in tuples)
            if isinstance(value, dict): recursive_items(value, sep, new_key) # Recursively process the sub-dictionary
            else: ls.append(f"{new_key}{sep}{value}") # Store the key-value pair as a tuple
    
    recursive_items(dc,sep) # Initialized recursive processing
    return ls

def ls_to_dc(ls: list, sep: str='.') -> dict:
    ''' 
    ls_to_dc(): convert a list with all the key value relationships stored as individual values into a dictionary containing several subdictionaries
    
    Parameters:
    ls (list): list
    sep (str, optional): seperator for subdictionaries for values in the list
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

# String methods
def find_all(string: str, substring: str) -> list:
    """
    find_all(): Find all indexes of a substring in a string
    
    Parameters:
    string (str): the string to search within
    substring (str): the substring to search for
    """
    indexes = []
    start = 0
    
    while start < len(string):
        
        start = string.find(substring, start) # Find the next occurrence of the substring
        if start == -1: break # If found, add to the list and move start position forward
        indexes.append(start)
        start += 1  # Move past the current match to find the next one
        
    return indexes

def color(text: str, color: str='green') -> None:
    '''
    color(): return colored text
    
    Parameters:
    text (str): uncolored text
    color (str, optional): text color [red, green, yellow, blue, magenta, cyan, white]
    '''
    colors = {
        'red': '\033[91m',
        'yellow': '\033[93m',
        'green': '\033[92m',
        'cyan': '\033[96m',
        'blue': '\033[94m',
        'magenta': '\033[95m',
        'white': '\033[97m',
        'reset': '\033[0m'
    }
    
    color_code = colors.get(color.lower(), colors['reset']) # Get the color code or default to reset
    reset_code = colors['reset'] # Reset code to clear formatting after the text
    
    return f"{color_code}{text}{reset_code}" # Return the colored text