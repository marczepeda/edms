''' 
Module: plot.py
Author: Marc Zepeda
Created: 2024-08-05
Description: Plot generation

Usage:
[Supporting methods]
- re_un_cap(): replace underscores with spaces and capitalizes each word for a given string
- round_up_pow_10(): rounds up a given number to the nearest power of 10
- round_down_pow_10: rounds down a given number to the nearest power of 10
- log10: returns log10 of maximum value from series or 0
- move_dis_legend(): moves legend for distribution graphs
- extract_pivots(): returns a dictionary of pivot-formatted dataframes from tidy-formatted dataframe
- formatter(): formats, displays, and saves plots
- repeat_palette_cmap(): returns a list of a repeated seaborn palette or matplotlib color map

[Graph methods]
- scat(): creates scatter plot related graphs
- cat(): creates category dependent graphs
- dist(): creates distribution graphs
- heat(): creates heat plot related graphs
- stack(): creates stacked bar plot
- vol(): creates volcano plot

[Color display methods]
- matplotlib_cmaps(): view all matplotlib color maps
- seaborn_palettes(): view all seaborn color palettes
'''

# Import packages
import os
import math
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.ticker import MaxNLocator
import numpy as np
from adjustText import adjust_text
from . import io

# Supporting methods
def re_un_cap(input: str) -> str:
    ''' 
    re_un_cap(): replace underscores with spaces and capitalizes each word for a given string
        
    Parameters:
    input (str): input string
    '''
    # Replace underscores with spaces
    input = input.replace('_', ' ')
    
    # Capitalize first letter of each word
    result = ''
    capitalize_next = True  # first letter too

    for char in input:
        if capitalize_next and char.isalpha():
            result += char.upper()
            capitalize_next = False
        else:
            result += char
            capitalize_next = (char == ' ')
    return result

def round_up_pow_10(number) -> int:
    ''' 
    round_up_pow_10(): rounds up a given number to the nearest power of 10
    
    Parameters:
    number (int or float): input number

    Depedencies: math
    '''
    if number == 0:
        return 0

    exponent = math.ceil(math.log10(abs(number)))
    rounded = math.ceil(number / 10 ** exponent) * 10 ** exponent
    return rounded

def round_down_pow_10(number) -> int:
    ''' 
    round_down_pow_10: rounds down a given number to the nearest power of 10
    
    Parameters:
    number: input number
    
    Dependencies: math
    '''
    
    if number == 0:
        return 0

    exponent = math.floor(math.log10(abs(number)))  # Use floor to round down the exponent
    rounded = math.floor(number / 10 ** exponent) * 10 ** exponent  # Round down the number
    return rounded

def log10(series) -> float:
    ''' 
    log10: returns log10 of maximum value from series or 0
    
    series: series, list, set, or array with values

    Dependencies: numpy
    '''
    return np.log10(np.maximum(series, 1))

def move_dist_legend(ax, legend_loc: str,legend_title_size: int, legend_size: int, legend_bbox_to_anchor: tuple, legend_ncol: tuple):
    ''' 
    move_dis_legend(): moves legend for distribution graphs
    
    Paramemters:
    ax: matplotlib axis
    legend_loc (str): legend location
    legend_title_size (str): legend title font size
    legend_size (str): legend font size
    legend_bbox_to_anchor (tuple): coordinates for bbox anchor
    legend_ncol (tuple): # of columns

    Dependencies: matplotlib.pyplot
    '''
    
    old_legend = ax.legend_
    handles = old_legend.legend_handles
    labels = [t.get_text() for t in old_legend.get_texts()]
    title = old_legend.get_title().get_text()
    ax.legend(handles,labels,loc=legend_loc,bbox_to_anchor=legend_bbox_to_anchor,
              title=title,title_fontsize=legend_title_size,fontsize=legend_size,ncol=legend_ncol)

def extract_pivots(df: pd.DataFrame, x: str, y: str, vars: str='variable', vals: str='value') -> dict[pd.DataFrame]:
    ''' 
    extract_pivots(): returns a dictionary of pivot-formatted dataframes from tidy-formatted dataframe
    
    Parameters:
    df (dataframe): tidy-formatted dataframe
    x (str): x-axis column name
    y (str): y-axis column name
    vars (str, optional): variable column name (variable)
    vals (str, optional): value column name (value)
    
    Dependencies: pandas
    '''
    piv_keys = list(df[vars].value_counts().keys())
    pivots = dict()
    for key in piv_keys:
        pivots[key]=pd.pivot(df[df[vars]==key],index=y,columns=x,values=vals)
    return pivots

def formatter(typ:str,ax,df:pd.DataFrame,x:str,y:str,cols:str,file:str,dir:str,palette_or_cmap:str,
              title:str,title_size:int,title_weight:str,title_font:str,
              x_axis:str,x_axis_size:int,x_axis_weight:str,x_axis_font:str,x_axis_scale:str,x_axis_dims:tuple,x_ticks_rot:int,x_ticks_font:str,x_ticks:list,
              y_axis:str,y_axis_size:int,y_axis_weight:str,y_axis_font:str,y_axis_scale:str,y_axis_dims:tuple,y_ticks_rot:int,y_ticks_font:str,y_ticks:list,
              legend_title:str,legend_title_size:int,legend_size:int,legend_bbox_to_anchor:tuple,legend_loc:str,legend_items:tuple,legend_ncol:int,show:bool,space_capitalize:bool):
    ''' 
    formatter(): formats, displays, and saves plots

    Parameters:
    typ (str): plot type
    ax: matplotlib axis
    df (dataframe): pandas dataframe
    x (str): x-axis column name
    y (str): y-axis column name
    cols (str, optional): color column name
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn color palette or matplotlib color map
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_scale (str, optional): x-axis scale linear, log, etc.
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    show (bool, optional): show plot (Default: True)
    space_capitalize (bool, optional): use re_un_cap() method when applicable (Default: True)

    Dependencies: os, matplotlib, seaborn, io, re_un_cap(), & round_up_pow_10()
    '''
    # Define plot types
    scats = ['scat', 'line', 'line_scat']
    cats = ['bar', 'box', 'violin', 'swarm', 'strip', 'point', 'count', 'bar_strip', 'box_strip', 'violin_strip','bar_swarm', 'box_swarm', 'violin_swarm']
    dists = ['hist', 'kde', 'hist_kde','rid']
    heats = ['ht']
        
    if typ not in heats:
        # Set title
        if title=='' and file is not None: 
            if space_capitalize: title=re_un_cap(".".join(file.split(".")[:-1]))
            else: title=".".join(file.split(".")[:-1])
        plt.title(title, fontsize=title_size, fontweight=title_weight, family=title_font)
        
        # Set x axis
        if x_axis=='': 
            if space_capitalize: x_axis=re_un_cap(x)
            else: x_axis=x
        plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight,fontfamily=x_axis_font)
        if x!='':
            if df[x].apply(lambda row: isinstance(row, (int, float))).all()==True: # Check that x column is numeric
                plt.xscale(x_axis_scale)
                if (x_axis_dims==(0,0))&(x_axis_scale=='log'): plt.xlim(round_down_pow_10(min(df[x])),round_up_pow_10(max(df[x])))
                elif x_axis_dims==(0,0): print('Default x axis dimensions.')
                else: plt.xlim(x_axis_dims[0],x_axis_dims[1])
        if x_ticks==[]: 
            if (x_ticks_rot==0)|(x_ticks_rot==90): plt.xticks(rotation=x_ticks_rot,ha='center',fontfamily=x_ticks_font)
            else: plt.xticks(rotation=x_ticks_rot,ha='right',fontfamily=x_ticks_font)
        else: 
            if (x_ticks_rot==0)|(x_ticks_rot==90): plt.xticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot, ha='center',fontfamily=x_ticks_font)
            else: plt.xticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot,ha='right',fontfamily=x_ticks_font)

        # Set y axis
        if y_axis=='': 
            if space_capitalize: y_axis=re_un_cap(y)
            else: y_axis=y
        plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight,fontfamily=y_axis_font)
        if y!='':
            if df[y].apply(lambda row: isinstance(row, (int, float))).all()==True: # Check that y column is numeric
                plt.yscale(y_axis_scale)
                if (y_axis_dims==(0,0))&(y_axis_scale=='log'): plt.ylim(round_down_pow_10(min(df[y])),round_up_pow_10(max(df[y])))
                elif y_axis_dims==(0,0): print('Default y axis dimensions.')
                else: plt.ylim(y_axis_dims[0],y_axis_dims[1])
        if y_ticks==[]: plt.yticks(rotation=y_ticks_rot,fontfamily=y_ticks_font)
        else: plt.yticks(ticks=y_ticks,labels=y_ticks,rotation=y_ticks_rot,fontfamily=y_ticks_font)

        # Set legend
        if cols is None: print('No legend because cols was not specified.')
        else:
            if legend_title=='': legend_title=cols
            if legend_items==(0,0) and typ not in dists:
                ax.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
                        bbox_to_anchor=legend_bbox_to_anchor,loc=legend_loc,ncol=legend_ncol) # Move legend to the right of the graph
            elif typ not in dists:
                handles, labels = ax.get_legend_handles_labels()
                ax.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
                        bbox_to_anchor=legend_bbox_to_anchor,loc=legend_loc,ncol=legend_ncol, # Move right of the graph
                        handles=handles[legend_items[0]:legend_items[1]],labels=labels[legend_items[0]:legend_items[1]]) # Only retains specified labels
            else: move_dist_legend(ax,legend_loc,legend_title_size,legend_size,legend_bbox_to_anchor,legend_ncol)

    # Save & show fig
    if file is not None and dir is not None:
        io.mkdir(dir) # Make output directory if it does not exist
        plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight', format=f'{file.split(".")[-1]}')
    if show: plt.show()

def repeat_palette_cmap(palette_or_cmap: str, repeats: int):
    '''
    repeat_palette_cmap(): returns a list of a repeated seaborn palette or matplotlib color map

    Parameters:
    palette_or_cmap (str): seaborn palette or matplotlib color map name
    repeats (int): number of color map repeats
    '''
    # Check repeats is a positive integer
    if not isinstance(repeats, int) or repeats <= 0:
        raise ValueError(f"repeats={repeats} must be a positive integer.")
    
    if palette_or_cmap in sns.color_palette(): # Check if cmap is a valid seaborn color palette
        cmap = sns.palettes.SEABORN_PALETTES[palette_or_cmap] # Get the color palette
        return mcolors.ListedColormap(cmap * repeats) # Repeats the color palette
    elif palette_or_cmap in plt.colormaps(): # Check if cmap is a valid matplotlib color map
        cmap = cm.get_cmap(palette_or_cmap) # Get the color map
        return mcolors.ListedColormap([cmap(i) for i in range(cmap.N)] * repeats) # Breaks the color map into a repeated list
    else:
        print(f'{cmap} is not a valid matplotlib color map and did not apply repeat.')
        return cmap

# Graph methods
def scat(typ: str, df: pd.DataFrame | str, x: str, y: str, cols: str=None, cols_ord: list=None, stys: str=None, cols_exclude: list | str=None,
         file: str=None, dir: str=None, palette_or_cmap: str='colorblind', edgecol: str='black',
         figsize: tuple=(10,6), title: str='', title_size: int=18, title_weight: str='bold', title_font: str='Arial',
         x_axis: str='', x_axis_size: int=12, x_axis_weight: str='bold', x_axis_font: str='Arial', x_axis_scale: str='linear', x_axis_dims: tuple=(0,0), x_ticks_rot: int=0, x_ticks_font: str='Arial', x_ticks: list=[],
         y_axis: str='', y_axis_size: int=12, y_axis_weight: str='bold', y_axis_font: str='Arial', y_axis_scale: str='linear', y_axis_dims: tuple=(0,0), y_ticks_rot: int=0, y_ticks_font: str='Arial', y_ticks: list=[],
         legend_title: str='', legend_title_size: int=12, legend_size: int=9, legend_bbox_to_anchor: tuple=(1,1), legend_loc: str='upper left',legend_items: tuple=(0,0), legend_ncol: int=1, show: bool=True, space_capitalize: bool=True, 
         **kwargs):
    ''' 
    scat(): creates scatter plot related graphs

    Parameters:
    typ (str): plot type (scat, line, line_scat)
    df (dataframe | str): pandas dataframe (or file path)
    x (str): x-axis column name
    y (str): y-axis column name
    cols (str, optional): color column name
    cols_ord (list, optional): color column values order
    stys (str, optional): styles column name
    cols_exclude (list | str, optional): color column values exclude
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn color palette or matplotlib color map
    edgecol (str, optional): point edge color
    figsize (tuple, optional): figure size
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    title_font (str, optional): plot title font
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_font (str, optional): x-axis font
    x_axis_scale (str, optional): x-axis scale linear, log, etc.
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_font (str, optional): x-ticks font
    x_ticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_font (str, optional): y-axis font
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks_font (str, optional): y_ticks font
    y_ticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    show (bool, optional): show plot (Default: True)
    space_capitalize (bool, optional): use re_un_cap() method when applicable (Default: True)
    
    Dependencies: os, matplotlib, seaborn, formatter(), re_un_cap(), & round_up_pow_10()
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)

    # Omit excluded data
    if type(cols_exclude)==list: 
        for exclude in cols_exclude: df=df[df[cols]!=exclude]
    elif type(cols_exclude)==str: df=df[df[cols]!=cols_exclude]

    # Set color scheme (Needs to be moved into individual plotting functions)
    color_palettes = ["deep", "muted", "bright", "pastel", "dark", "colorblind", "husl", "hsv", "Paired", "Set1", "Set2", "Set3", "tab10", "tab20"] # List of common Seaborn palettes
    if palette_or_cmap in color_palettes: palette = palette_or_cmap
    elif palette_or_cmap in plt.colormaps(): 
        if cols is not None: # Column specified
            cmap = cm.get_cmap(palette_or_cmap,len(df[cols].value_counts()))
            palette = sns.color_palette([cmap(i) for i in range(cmap.N)])
        else:
            print('Cols not specified. Used seaborn colorblind.')
            palette = 'colorblind'
    else: 
        print('Seaborn color palette or matplotlib color map not specified. Used seaborn colorblind.')
        palette = 'colorblind'

    fig, ax = plt.subplots(figsize=figsize)
    
    if cols is not None and stys is not None:
        if typ=='scat': sns.scatterplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, style=stys, edgecolor=edgecol, palette=palette, ax=ax, **kwargs)
        elif typ=='line': sns.lineplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, style=stys, palette=palette, ax=ax, **kwargs)
        elif typ=='line_scat':
            sns.lineplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, style=stys, palette=palette, ax=ax, **kwargs)  
            sns.scatterplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, style=stys, edgecolor=edgecol, palette=palette, ax=ax, **kwargs)
        else:
            print("Invalid type! scat, line, or line_scat")
            return
    elif cols is not None:
        if typ=='scat': sns.scatterplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, palette=palette, ax=ax, **kwargs)
        elif typ=='line': sns.lineplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, ax=ax, palette=palette, **kwargs)
        elif typ=='line_scat':
            sns.lineplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, palette=palette, ax=ax, **kwargs)  
            sns.scatterplot(data=df, x=x, y=y, hue=cols, hue_order=cols_ord, edgecolor=edgecol, palette=palette, ax=ax, **kwargs)
        else:
            print("Invalid type! scat, line, or line_scat")
            return
    elif stys is not None:
        if typ=='scat': sns.scatterplot(data=df, x=x, y=y, style=stys, edgecolor=edgecol, palette=palette, ax=ax, **kwargs)
        elif typ=='line': sns.lineplot(data=df, x=x, y=y, style=stys, palette=palette, ax=ax, **kwargs)
        elif typ=='line_scat':
            sns.lineplot(data=df, x=x, y=y, style=stys, palette=palette, ax=ax, **kwargs)  
            sns.scatterplot(data=df, x=x, y=y, style=stys, edgecolor=edgecol, palette=palette, ax=ax, **kwargs)
        else:
            print("Invalid type! scat, line, or line_scat")
            return
    else:
        if typ=='scat': sns.scatterplot(data=df, x=x, y=y, edgecolor=edgecol, palette=palette, ax=ax, **kwargs)
        elif typ=='line': sns.lineplot(data=df, x=x, y=y, palette=palette, ax=ax, **kwargs)
        elif typ=='line_scat':
            sns.lineplot(data=df, x=x, y=y, palette=palette, ax=ax, **kwargs)  
            sns.scatterplot(data=df, x=x, y=y, edgecolor=edgecol, palette=palette, ax=ax, **kwargs)
        else:
            print("Invalid type! scat, line, or line_scat")
            return
    
    formatter(typ,ax,df,x,y,cols,file,dir,palette_or_cmap,
              title,title_size,title_weight,title_font,
              x_axis,x_axis_size,x_axis_weight,x_axis_font,x_axis_scale,x_axis_dims,x_ticks_rot,x_ticks_font,x_ticks,
              y_axis,y_axis_size,y_axis_weight,y_axis_font,y_axis_scale,y_axis_dims,y_ticks_rot,y_ticks_font,y_ticks,
              legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol,show,space_capitalize)

def cat(typ:str, df:pd.DataFrame | str, x: str='', y: str='', cols: str=None, cats_ord: list=None, cols_ord: list=None, cols_exclude: list | str=None,
        file: str=None, dir: str=None, palette_or_cmap: str='colorblind', edgecol: str='black', lw: int=1, errorbar: str='sd', errwid: int=1, errcap: float=0.1,
        figsize: tuple=(10,6), title: str='', title_size: int=18, title_weight: str='bold', title_font: str='Arial',
        x_axis: str='', x_axis_size=12, x_axis_weight: str='bold', x_axis_font: str='Arial', x_axis_scale: str='linear', x_axis_dims: tuple=(0,0), x_ticks_rot: int=0, x_ticks_font: str='Arial', x_ticks: list=[],
        y_axis: str='', y_axis_size=12, y_axis_weight: str='bold', y_axis_font: str='Arial', y_axis_scale: str='linear', y_axis_dims: tuple=(0,0), y_ticks_rot: int=0, y_ticks_font: str='Arial', y_ticks: list=[],
        legend_title: str='', legend_title_size: int=12, legend_size: int=9, legend_bbox_to_anchor=(1,1), legend_loc: str='upper left', legend_items: tuple=(0,0), legend_ncol: int=1, show: bool=True, space_capitalize: bool=True, 
        **kwargs):
    ''' 
    cat(): creates category dependent graphs

    Parameters:
    typ (str): plot type (bar, box, violin, swarm, strip, point, count, bar_swarm, box_swarm, violin_swarm)
    df (dataframe | str): pandas dataframe (or file path)
    x (str, optional): x-axis column name
    y (str, optional): y-axis column name
    cols (str, optional): color column name
    cats_ord (list, optional): category column values order (x- or y-axis)
    cols_ord (list, optional): color column values order
    cols_exclude (list | str, optional): color column values exclude
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn color palette or matplotlib color map
    edgecol (str, optional): point edge color
    lw (int, optional): line width
    errorbar (str, optional): error bar type (sd)
    errwid (int, optional): error bar line width
    errcap (int, optional): error bar cap line width
    figsize (tuple, optional): figure size
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    title_font (str, optional): plot title font
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_font (str, optional): x-axis font
    x_axis_scale (str, optional): x-axis scale linear, log, etc.
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_font (str, optional): x-ticks font
    x_ticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_font (str, optional): y-axis font
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_font (str, optional): y_ticks font
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    show (bool, optional): show plot (Default: True)
    space_capitalize (bool, optional): use re_un_cap() method when applicable (Default: True)
    
    Dependencies: os, matplotlib, seaborn, formatter(), re_un_cap(), & round_up_pow_10()
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)
    
    # Omit excluded data
    if type(cols_exclude)==list: 
        for exclude in cols_exclude: df=df[df[cols]!=exclude]
    elif type(cols_exclude)==str: df=df[df[cols]!=cols_exclude]

    # Set color scheme and category column order
    color_palettes = ["deep", "muted", "bright", "pastel", "dark", "colorblind", "husl", "hsv", "Paired", "Set1", "Set2", "Set3", "tab10", "tab20"] # List of common Seaborn palettes
    if palette_or_cmap in color_palettes: palette = palette_or_cmap
    elif palette_or_cmap in plt.colormaps(): 
        if cols is not None: # Column specified
            cmap = cm.get_cmap(palette_or_cmap,len(df[cols].value_counts()))
            palette = sns.color_palette([cmap(i) for i in range(cmap.N)])
        elif (x!='')&(y!=''): # x- and y-axis are specified
            if df[x].apply(lambda row: isinstance(row, str)).all()==True: # Check x column is categorical
                cmap = cm.get_cmap(palette_or_cmap,len(df[x].value_counts()))
                palette = sns.color_palette([cmap(i) for i in range(cmap.N)])
            elif df[y].apply(lambda row: isinstance(row, str)).all()==True: # Check y column is categorical
                cmap = cm.get_cmap(palette_or_cmap,len(df[y].value_counts()))
                palette = sns.color_palette([cmap(i) for i in range(cmap.N)])
        elif (x!='')&(y==''): # x-axis is specified
            if df[x].apply(lambda row: isinstance(row, str)).all()==True: # Check x column is categorical
                cmap = cm.get_cmap(palette_or_cmap,len(df[x].value_counts()))
                palette = sns.color_palette([cmap(i) for i in range(cmap.N)])
        elif (x=='')&(y!=''): # y-axis is specified
            if df[y].apply(lambda row: isinstance(row, str)).all()==True: # Check y column is categorical
                cmap = cm.get_cmap(palette_or_cmap,len(df[y].value_counts()))
                palette = sns.color_palette([cmap(i) for i in range(cmap.N)])
        else: return
    else: 
        print('Seaborn color palette or matplotlib color map not specified. Used seaborn colorblind.')
        palette = 'colorblind'

    fig, ax = plt.subplots(figsize=figsize)

    if cols is not None:

        if typ=='bar': sns.barplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, errcolor=edgecol, errwidth=errwid, capsize=errcap, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='box': sns.boxplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='violin': sns.violinplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='swarm': sns.swarmplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='strip': sns.stripplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='point': sns.pointplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, errwidth=errwid, capsize=errcap, hue=cols, hue_order=cols_ord, palette=palette, ax=ax, **kwargs)
        elif typ=='count': 
            if (x!='')&(y!=''):
                print('Cannot make countplot with both x and y specified.')
                return
            elif x!='': sns.countplot(data=df, x=x, order=cats_ord, hue=cols, hue_order=cols_ord, palette=palette, ax=ax, **kwargs)
            elif y!='': sns.countplot(data=df, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, palette=palette, ax=ax, **kwargs)
            else:
                print('Cannot make countplot without x or y specified.')
                return
        elif typ=='bar_strip':
            sns.barplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, errcolor=edgecol, errwidth=errwid, capsize=errcap, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='box_strip':
            sns.boxplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='violin_strip':
            sns.violinplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='bar_swarm':
            sns.barplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, errcolor=edgecol, errwidth=errwid, capsize=errcap, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='box_swarm':
            sns.boxplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='violin_swarm':
            sns.violinplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        else:
            print('Invalid type! bar, box, violin, swarm, strip, point, count, bar_strip, box_strip, violin_strip, bar_swarm, box_swarm, violin_swarm')
            return

    else: # Cols was not specified
        
        if typ=='bar': sns.barplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, errcolor=edgecol, errwidth=errwid, capsize=errcap, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='box': sns.boxplot(data=df, x=x, y=y, order=cats_ord, linewidth=lw, ax=ax, palette=palette, **kwargs)
        elif typ=='violin': sns.violinplot(data=df, x=x, y=y, order=cats_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='swarm': sns.swarmplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, linewidth=lw, dodge=True, palette=palette, ax=ax, **kwargs)
        elif typ=='strip': sns.stripplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='point': sns.pointplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, errwidth=errwid, capsize=errcap, palette=palette, ax=ax, **kwargs)
        elif typ=='count': 
            if (x!='')&(y!=''):
                print('Cannot make countplot with both x and y specified.')
                return
            elif x!='': sns.countplot(data=df, x=x, order=cats_ord, ax=ax, palette=palette, **kwargs)
            elif y!='': sns.countplot(data=df, y=y, order=cats_ord, ax=ax, palette=palette, **kwargs)
            else:
                print('Cannot make countplot without x or y specified.')
                return
        elif typ=='bar_strip':
            sns.barplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, errcolor=edgecol, errwidth=errwid, capsize=errcap, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='box_strip':
            sns.boxplot(data=df, x=x, y=y, order=cats_ord, linewidth=lw, ax=ax, **kwargs)
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='violin_strip':
            sns.violinplot(data=df, x=x, y=y, order=cats_ord, edgecolor=edgecol, linewidth=lw, ax=ax, palette=palette, **kwargs)
            sns.stripplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, linewidth=lw, ax=ax, palette=palette, **kwargs)
        elif typ=='bar_swarm':
            sns.barplot(data=df, x=x, y=y, order=cats_ord, errorbar=errorbar, errcolor=edgecol, errwidth=errwid, capsize=errcap, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='box_swarm':
            sns.boxplot(data=df, x=x, y=y, order=cats_ord, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        elif typ=='violin_swarm':
            sns.violinplot(data=df, x=x, y=y, order=cats_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            sns.swarmplot(data=df, x=x, y=y, order=cats_ord, color=edgecol, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        else:
            print('Invalid type! bar, box, violin, swarm, strip, point, count, bar_strip, box_strip, violin_strip, bar_swarm, box_swarm, violin_swarm')
            return

    formatter(typ,ax,df,x,y,cols,file,dir,palette_or_cmap,
              title,title_size,title_weight,title_font,
              x_axis,x_axis_size,x_axis_weight,x_axis_font,x_axis_scale,x_axis_dims,x_ticks_rot,x_ticks_font,x_ticks,
              y_axis,y_axis_size,y_axis_weight,y_axis_font,y_axis_scale,y_axis_dims,y_ticks_rot,y_ticks_font,y_ticks,
              legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol,show,space_capitalize)

def dist(typ: str, df: pd.DataFrame | str, x: str, cols: str=None, cols_ord: list=None, cols_exclude: list | str=None, bins: int=40, log10_low: int=0,
        file: str=None, dir: str=None, palette_or_cmap: str='colorblind', edgecol: str='black',lw: int=1, ht: float=1.5, asp: int=5, tp: float=.8, hs: int=0, despine: bool=False,
        figsize=(10,6), title: str='', title_size: int=18, title_weight: str='bold', title_font: str='Arial',
        x_axis: str='', x_axis_size: int=12, x_axis_weight: str='bold', x_axis_font: str='Arial', x_axis_scale: str='linear', x_axis_dims: tuple=(0,0), x_ticks_rot: int=0, x_ticks_font: str='Arial', x_ticks: list=[],
        y_axis: str='', y_axis_size: int=12, y_axis_weight: str='bold', y_axis_font: str='Arial', y_axis_scale: str='linear', y_axis_dims: tuple=(0,0), y_ticks_rot: int=0, y_ticks_font: str='Arial', y_ticks: list=[],
        legend_title: str='', legend_title_size: int=12, legend_size: int=9, legend_bbox_to_anchor: tuple=(1,1), legend_loc: str='upper left', legend_items: tuple=(0,0), legend_ncol: int=1, show: bool=True, space_capitalize: bool=True, 
        **kwargs):
    ''' 
    dist(): creates distribution graphs

    Parameters:
    typ (str): plot type (hist, kde, hist_kde, rid)
    df (dataframe | str): pandas dataframe (or file path)
    x (str): x-axis column name
    cols (str, optional): color column name
    cols_ord (list, optional): color column values order
    cols_exclude (list | str, optional): color column values exclude
    bins (int, optional): # of bins for histogram
    log10_low (int, optional): log scale lower bound
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn color palette or matplotlib color map
    edgecol (str, optional): point edge color
    lw (int, optional): line width
    ht (float, optional): height
    asp (int, optional): aspect
    tp (float, optional): top
    hs (int, optional): hspace
    despine (bool, optional): despine
    figsize (tuple, optional): figure size
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    title_font (str, optional): plot title font
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_font (str, optional): x-axis font
    x_axis_scale (str, optional): x-axis scale linear, log, etc.
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_font (str, optional): x-ticks font
    x_ticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_font (str, optional): y-axis font
    y_axis_scale (str, optional): y-axis scale linear, log, etc.
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks_font (str, optional): y_ticks font
    y_ticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    show (bool, optional): show plot (Default: True)
    space_capitalize (bool, optional): use re_un_cap() method when applicable (Default: True)
    
    Dependencies: os, matplotlib, seaborn, io, formatter(), re_un_cap(), & round_up_pow_10()
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)
    
    # Omit excluded data
    if type(cols_exclude)==list: 
        for exclude in cols_exclude: df=df[df[cols]!=exclude]
    elif type(cols_exclude)==str: df=df[df[cols]!=cols_exclude]

    # Set color scheme (Needs to be moved into individual plotting functions)
    color_palettes = ["deep", "muted", "bright", "pastel", "dark", "colorblind", "husl", "hsv", "Paired", "Set1", "Set2", "Set3", "tab10", "tab20"] # List of common Seaborn palettes
    if palette_or_cmap in color_palettes: palette = palette_or_cmap
    elif palette_or_cmap in plt.colormaps(): 
        if cols is not None: # Column specified
            cmap = cm.get_cmap(palette_or_cmap,len(df[cols].value_counts()))
            palette = sns.color_palette([cmap(i) for i in range(cmap.N)])
        else:
            print('Cols not specified. Used seaborn colorblind.')
            palette = 'colorblind'
    else: 
        print('Seaborn color palette or matplotlib color map not specified. Used seaborn colorblind.')
        palette = 'colorblind'

    if typ=='hist':
        fig, ax = plt.subplots(figsize=figsize)
        if isinstance(bins, int):
            if x_axis_scale=='log': bins = np.logspace(log10(df[x]).min(), log10(df[x]).max(), bins + 1)
            else: bins = np.linspace(df[x].min(), df[x].max(), bins + 1)
        sns.histplot(data=df, x=x, kde=False, bins=bins, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        y=''
        if y_axis=='': y_axis='Count'
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        formatter(typ,ax,df,x,y,cols,file,dir,palette_or_cmap,
                  title,title_size,title_weight,title_font,
                  x_axis,x_axis_size,x_axis_weight,x_axis_font,x_axis_scale,x_axis_dims,x_ticks_rot,x_ticks_font,x_ticks,
                  y_axis,y_axis_size,y_axis_weight,y_axis_font,y_axis_scale,y_axis_dims,y_ticks_rot,y_ticks_font,y_ticks,
                  legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol,show,space_capitalize)
    elif typ=='kde': 
        fig, ax = plt.subplots(figsize=figsize)
        if x_axis_scale=='log':
            df[f'log10({x})']=np.maximum(np.log10(df[x]),log10_low)
            sns.kdeplot(data=df, x=f'log10({x})', hue=cols, hue_order=cols_ord, linewidth=lw, palette=palette, ax=ax, **kwargs)
            x_axis_scale='linear'
            if x_axis=='': x_axis=f'log10({x})'
        else: sns.kdeplot(data=df, x=x, hue=cols, hue_order=cols_ord, linewidth=lw, ax=ax, **kwargs)
        y=''
        if y_axis=='': y_axis='Density'
        formatter(typ,ax,df,x,y,cols,file,dir,palette_or_cmap,
                  title,title_size,title_weight,title_font,
                  x_axis,x_axis_size,x_axis_weight,x_axis_font,x_axis_scale,x_axis_dims,x_ticks_rot,x_ticks_font,x_ticks,
                  y_axis,y_axis_size,y_axis_weight,y_axis_font,y_axis_scale,y_axis_dims,y_ticks_rot,y_ticks_font,y_ticks,
                  legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol,show,space_capitalize)
    elif typ=='hist_kde':
        fig, ax = plt.subplots(figsize=figsize)
        if x_axis_scale=='log':
            df[f'log10({x})']=np.maximum(np.log10(df[x]),log10_low)
            bins = np.logspace(log10(df[x]).min(), log10(df[x]).max(), bins + 1)
            sns.histplot(data=df, x=f'log10({x})', kde=True, bins=bins, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
            x_axis_scale='linear'
            if x_axis=='': x_axis=f'log10({x})'
        else:
            bins = np.linspace(df[x].min(), df[x].max(), bins + 1) 
            sns.histplot(data=df, x=x, kde=True, bins=bins, hue=cols, hue_order=cols_ord, edgecolor=edgecol, linewidth=lw, palette=palette, ax=ax, **kwargs)
        y=''
        if y_axis=='': y_axis='Count'
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        formatter(typ,ax,df,x,y,cols,file,dir,palette_or_cmap,
                  title,title_size,title_weight,title_font,
                  x_axis,x_axis_size,x_axis_weight,x_axis_font,x_axis_scale,x_axis_dims,x_ticks_rot,x_ticks_font,x_ticks,
                  y_axis,y_axis_size,y_axis_weight,y_axis_font,y_axis_scale,y_axis_dims,y_ticks_rot,y_ticks_font,y_ticks,
                  legend_title,legend_title_size,legend_size,legend_bbox_to_anchor,legend_loc,legend_items,legend_ncol,show,space_capitalize)
    elif typ=='rid':
        # Set color scheme
        color_palettes = ["deep", "muted", "bright", "pastel", "dark", "colorblind", "husl", "hsv", "Paired", "Set1", "Set2", "Set3", "tab10", "tab20"] # List of common Seaborn palettes
        if (palette_or_cmap in set(color_palettes))|(palette_or_cmap in set(plt.colormaps())): sns.color_palette(palette_or_cmap)
        else: 
            print('Seaborn color palette or matplotlib color map not specified. Used seaborn colorblind.')
            sns.color_palette('colorblind')
        if x_axis_scale=='log':
            df[f'log10({x})']=np.maximum(np.log10(df[x]),log10_low)
            g = sns.FacetGrid(df, row=cols, hue=cols, col_order=cols_ord, hue_order=cols_ord, height=ht, aspect=asp)
            g.map(sns.kdeplot, f'log10({x})', linewidth=lw, **kwargs)
            if x_axis=='': x_axis=f'log10({x})'
        else:
            g = sns.FacetGrid(df, row=cols, hue=cols, col_order=cols_ord, hue_order=cols_ord, height=ht, aspect=asp)
            g.map(sns.kdeplot, x, linewidth=lw, **kwargs)
            if x_axis=='': x_axis=x
        for ax in g.axes.flatten():
            if x_axis_dims!=(0,0): ax.set_xlim(x_axis_dims[0],x_axis_dims[1]) # This could be an issue with the (0,0) default (Figure out later...)
            ax.set_xlabel(x_axis,fontsize=x_axis_size,fontweight=x_axis_weight,fontfamily=x_axis_font)
        if y_axis=='': y_axis='Density'
        g.set(yticks=y_ticks, ylabel=y_axis) # fontfamily only works on the ax level (Figure out later if I care...)
        g.set_titles("")
        if title=='' and file is not None: 
            if space_capitalize: title=re_un_cap(".".join(file.split(".")[:-1]))
            else: ".".join(file.split(".")[:-1])
        g.figure.suptitle(title, fontsize=title_size, fontweight=title_weight,fontfamily=title_font)
        g.figure.subplots_adjust(top=tp,hspace=hs)
        if despine==False: g.despine(top=False,right=False)
        else: g.despine(left=True)
        if legend_title=='': legend_title=cols
        g.figure.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
                        loc=legend_loc,bbox_to_anchor=legend_bbox_to_anchor)
        if file is not None and dir is not None:
            io.mkdir(dir) # Make output directory if it does not exist
            plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight', format=f'{file.split(".")[-1]}')
        if show: plt.show()
    else:
        print('Invalid type! hist, kde, hist_kde, rid')
        return

def heat(df: pd.DataFrame | str, x: str=None, y: str=None, vars: str=None, vals: str=None, vals_dims:tuple=None,
         file: str=None, dir: str=None, edgecol: str='black', lw: int=1, annot: bool=False, cmap: str="Reds", sq: bool=True, cbar: bool=True,
         title: str='',title_size: int=18, title_weight: str='bold', title_font: str='Arial', figsize: tuple=(10,6),
         x_axis: str='', x_axis_size: int=12, x_axis_weight: str='bold', x_axis_font: str='Arial', x_ticks_rot: int=0, x_ticks_font: str='Arial',
         y_axis: str='', y_axis_size: int=12, y_axis_weight: str='bold', y_axis_font: str='Arial', y_ticks_rot: int=0, y_ticks_font: str='Arial',
         show: bool=True, space_capitalize: bool=True, **kwargs):
    '''
    heat(): creates heat plot related graphs

    Parameters:
    df (dataframe | str): pandas dataframe (or file path)
    x (str, optional): x-axis column name to split tidy-formatted dataframe into a dictionary pivot-formatted dataframes (Default: None)
    y (str, optional): y-axis column name to split tidy-formatted dataframe into a dictionary pivot-formatted dataframes (Default: None)
    vars (str, optional): variable column name to split tidy-formatted dataframe into a dictionary pivot-formatted dataframes (Default: None)
    vals (str, optional): value column name to split tidy-formatted dataframe into a dictionary pivot-formatted dataframes (Default: None)
    vals_dims (tuple, optional): value column minimum and maximum formatted (vmin, vmax; Default: None)
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    edgecol (str, optional): point edge color
    lw (int, optional): line width
    annot (bool, optional): annotate values
    cmap (str, optional): matplotlib color map
    sq (bool, optional): square dimensions (Default: True)
    cbar (bool, optional): show colorbar (Default: True)
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    title_font (str, optional): plot title font
    figsize (tuple, optional): figure size per subplot
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_font (str, optional): x-axis font
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_font (str, optional): x-ticks font
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_font (str, optional): y-axis font
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks_font (str, optional): y-ticks font
    show (bool, optional): show plot (Default: True)
    space_capitalize (bool, optional): use re_un_cap() method when applicable (Default: True)
    
    Dependencies: os, matplotlib, seaborn, formatter(), re_un_cap(), & round_up_pow_10()
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)

    # Determine dataframe type
    if x is None or y is None or vars is None or vals is None: # Pivot-formatted

        # Find min and max values in the dataset for normalization
        if vals_dims is None:
            vmin = df.min().min()
            vmax = df.max().max()
        else:
            vmin = vals_dims[0]
            vmax = vals_dims[1]

        # Create dictionary of pivot-formatted dataframes
        dc = {'Pivot Table': df}
        x = df.columns.name
        y = df.index.name

    else: # Tidy-formatted
        
        # Find min and max values in the dataset for normalization
        if vals_dims is None:
            vmin = df[vals].values.min()
            vmax = df[vals].values.max()
        else:
            vmin = vals_dims[0]
            vmax = vals_dims[1]

        # Create dictionary of pivot-formatted dataframes
        dc = extract_pivots(df=df,x=x,y=y,vars=vars,vals=vals)

    # Create a single figure with multiple heatmap subplots
    fig, axes = plt.subplots(nrows=len(list(dc.keys())),ncols=1,figsize=(figsize[0],figsize[1]*len(list(dc.keys()))),sharex=False,sharey=True)
    if isinstance(axes, np.ndarray)==False: axes = np.array([axes]) # Make axes iterable if there is only 1 heatmap
    for (ax, key) in zip(axes, list(dc.keys())):
        sns.heatmap(dc[key],annot=annot,cmap=cmap,ax=ax,linecolor=edgecol,linewidths=lw,cbar=cbar,square=sq,vmin=vmin,vmax=vmax, **kwargs)
        if len(list(dc.keys()))>1: ax.set_title(key,fontsize=title_size,fontweight=title_weight,fontfamily=title_font)  # Add title to subplot
        else: ax.set_title(title,fontsize=title_size,fontweight=title_weight,fontfamily=title_font)
        if x_axis=='': 
            if space_capitalize: ax.set_xlabel(re_un_cap(x),fontsize=x_axis_size,fontweight=x_axis_weight,fontfamily=x_axis_font) # Add x axis label
            else: ax.set_xlabel(x,fontsize=x_axis_size,fontweight=x_axis_weight,fontfamily=x_axis_font) # Add x axis label
        else: ax.set_xlabel(x_axis,fontsize=x_axis_size,fontweight=x_axis_weight,fontfamily=x_axis_font)
        if y_axis=='': 
            if space_capitalize: ax.set_ylabel(re_un_cap(y),fontsize=y_axis_size,fontweight=y_axis_weight,fontfamily=y_axis_font) # Add y axis label
            else: ax.set_ylabel(y,fontsize=y_axis_size,fontweight=y_axis_weight,fontfamily=y_axis_font) # Add y axis label
        else: ax.set_ylabel(y_axis,fontsize=y_axis_size,fontweight=y_axis_weight,fontfamily=y_axis_font)
        # Format x ticks
        if (x_ticks_rot==0)|(x_ticks_rot==90): plt.setp(ax.get_xticklabels(), rotation=x_ticks_rot, ha="center", va='top', rotation_mode="anchor",fontname=x_ticks_font) 
        else: plt.setp(ax.get_xticklabels(), rotation=x_ticks_rot, ha="right", va='top', rotation_mode="anchor",fontname=x_ticks_font) 
        # Format y ticks
        plt.setp(ax.get_yticklabels(), rotation=y_ticks_rot, va='center', ha="right",rotation_mode="anchor",fontname=y_ticks_font)
        ax.set_facecolor('white')  # Set background to transparent
    
    # Save & show fig
    if file is not None and dir is not None:
        io.mkdir(dir) # Make output directory if it does not exist
        plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight', format=f'{file.split(".")[-1]}')
    if show: plt.show()

def stack(df: pd.DataFrame | str, x: str, y: str, cols: str, cutoff: float=0, cols_ord: list=[], x_ord: list=[],
          file: str=None, dir: str=None, palette_or_cmap: str='tab20', repeats: int=1, errcap: int=4, vertical: bool=True,
          figsize: tuple=(10,6), title: str='', title_size: int=18, title_weight: str='bold', title_font: str='Arial',
          x_axis: str='', x_axis_size: int=12, x_axis_weight: str='bold', x_axis_font: str='Arial', x_ticks_rot: int=0, x_ticks_font: str='Arial',
          y_axis: str='', y_axis_size: int=12, y_axis_weight: str='bold', y_axis_font: str='Arial', y_axis_dims: tuple=(0,0), y_ticks_rot: int=0, y_ticks_font: str='Arial',
          legend_title: str='', legend_title_size: int=12, legend_size: int=12,
          legend_bbox_to_anchor: tuple=(1,1), legend_loc: str='upper left', legend_ncol: int=1, show: bool=True, space_capitalize: bool=True, **kwargs):
    ''' 
    stack(): creates stacked bar plot

    Parameters:
    df (dataframe | str): pandas dataframe (or file path)
    x (str, optional): x-axis column name
    y (str, optional): y-axis column name
    cols (str, optional): color column name
    cutoff (float, optional): y-axis values needs be greater than (e.g. 0)
    cols_ord (list, optional): color column values order
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    palette_or_cmap (str, optional): seaborn palette or matplotlib color map
    repeats (int, optional): number of color palette or map repeats (Default: 1)
    errcap (int, optional): error bar cap line width
    vertical (bool, optional): vertical orientation; otherwise horizontal (Default: True)
    figsize (tuple, optional): figure size
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    title_font (str, optional): plot title font
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_font (str, optional): x-axis font
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_ticks_font (str, optional): x-ticks font
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_font (str, optional): y-axis font
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks_font (str, optional): y-ticks font
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    show (bool, optional): show plot (Default: True)
    space_capitalize (bool, optional): use re_un_cap() method when applicable (Default: True)
    
    Dependencies: re, os, pandas, numpy, matplotlib.pyplot, & io
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)

    # Omit smaller than cutoff and convert it to <cutoff
    df_cut=df[df[y]>=cutoff]
    df_other=df[df[y]<cutoff]
    for x_val in list(df_other[x].value_counts().keys()):
        df_temp = df_other[df_other[x]==x_val]
        df_temp[y]=sum(df_temp[y])
        df_temp[cols]=f'<{cutoff}'
        df_cut = pd.concat([df_cut,df_temp.iloc[:1]])

    # Make pivot table
    df_cut=df[df[y]>=cutoff]
    df_pivot=pd.pivot_table(df_cut, index=x, columns=cols, values=y, aggfunc=np.mean)
    df_pivot_err=pd.pivot_table(df_cut, index=x, columns=cols, values=y, aggfunc=np.std)
    if cols_ord!=[]: df_pivot=df_pivot.reindex(columns=cols_ord)
    if x_ord!=[]: df_pivot=df_pivot.reindex(index=x_ord)

    # Make stacked barplot
    if vertical: # orientation
        df_pivot.plot(kind='bar',yerr=df_pivot_err,capsize=errcap, figsize=figsize,colormap=repeat_palette_cmap(palette_or_cmap,repeats),stacked=True,**kwargs)
        
        # Set x axis
        if x_axis=='': 
            if space_capitalize: x_axis=re_un_cap(x)
            else: x_axis=x
        plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight,fontfamily=x_axis_font)
        if (x_ticks_rot==0)|(x_ticks_rot==90): plt.xticks(rotation=x_ticks_rot,ha='center',fontfamily=x_ticks_font)
        else: plt.xticks(rotation=x_ticks_rot,ha='right',fontfamily=x_ticks_font)
        
        # Set y axis
        if y_axis=='': 
            if space_capitalize: y_axis=re_un_cap(y)
            else: y_axis=y
        plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight,fontfamily=y_axis_font)
        plt.yticks(rotation=y_ticks_rot,fontfamily=y_ticks_font)

        if y_axis_dims==(0,0): print('Default y axis dimensions.')
        else: plt.ylim(y_axis_dims[0],y_axis_dims[1])

    else: # Horizontal orientation
        df_pivot.plot(kind='barh',xerr=df_pivot_err,capsize=errcap, figsize=figsize,colormap=repeat_palette_cmap(palette_or_cmap,repeats),stacked=True,**kwargs)

        # Set y axis
        if x_axis=='': 
            if space_capitalize: x_axis=re_un_cap(x)
            else: x_axis=x
        plt.ylabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight,fontfamily=x_axis_font)
        plt.yticks(rotation=x_ticks_rot,fontfamily=x_ticks_font)
        
        # Set x axis
        if y_axis=='': 
            if space_capitalize: y_axis=re_un_cap(y)
            else: y_axis=y
        plt.xlabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight,fontfamily=y_axis_font)
        if (y_ticks_rot==0)|(y_ticks_rot==90): plt.xticks(rotation=y_ticks_rot,ha='center',fontfamily=y_ticks_font)
        else: plt.xticks(rotation=y_ticks_rot,ha='right',fontfamily=y_ticks_font)

        if y_axis_dims==(0,0): print('Default x axis dimensions.')
        else: plt.xlim(y_axis_dims[0],y_axis_dims[1])
        
    # Set title
    if title=='' and file is not None: title=re_un_cap(".".join(file.split(".")[:-1]))
    plt.title(title, fontsize=title_size, fontweight=title_weight, family=title_font)
    
    # Set legend
    if legend_title=='': 
        if space_capitalize: legend_title=re_un_cap(cols)
        else: legend_title=cols
    plt.legend(title=legend_title, title_fontsize=legend_title_size, fontsize=legend_size, 
               bbox_to_anchor=legend_bbox_to_anchor, loc=legend_loc, ncol=legend_ncol)
    
    # Save & show fig
    if file is not None and dir is not None:
        io.mkdir(dir) # Make output directory if it does not exist
        plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight', format=f'{file.split(".")[-1]}')
    if show: plt.show()

def vol(df: pd.DataFrame | str, x: str, y: str, stys: str=None, size: str=None, size_dims: tuple=None, label: str=None,
        FC_threshold: float=2, pval_threshold: float=0.05, file: str=None, dir: str=None, color: str='lightgray', alpha: float=0.5, edgecol: str='black', vertical: bool=True,
        figsize=(10,6), title: str='', title_size: int=18, title_weight: str='bold', title_font: str='Arial',
        x_axis: str='', x_axis_size: int=12, x_axis_weight: str='bold', x_axis_font: str='Arial', x_axis_dims: tuple=(0,0), x_ticks_rot: int=0, x_ticks_font: str='Arial', x_ticks: list=[],
        y_axis: str='', y_axis_size: int=12, y_axis_weight: str='bold', y_axis_font: str='Arial', y_axis_dims: tuple=(0,0), y_ticks_rot: int=0, y_ticks_font: str='Arial', y_ticks: list=[],
        legend_title: str='', legend_title_size: int=12, legend_size: int=9, legend_bbox_to_anchor: tuple=(1,1), legend_loc: str='upper left',
        legend_items: tuple=(0,0),legend_ncol: int=1 ,display_size: bool=True, display_labels: bool=True, return_df: bool=True, show: bool=True, space_capitalize: bool=True,
        **kwargs) -> pd.DataFrame:
    ''' 
    vol(): creates volcano plot
    
    Parameters:
    df (dataframe | str): pandas dataframe (or file path)
    x (str): x-axis column name (FC)
    y (str): y-axis column name (pval)
    stys (str, optional): style column name
    size (str, optional): size column name
    size_dims (tuple, optional): (minimum,maximum) values in size column (Default: None)
    label (str, optional): label column name
    FC_threshold (float, optional): fold change threshold (Default: 2; log2(2)=1)
    pval_threshold (float, optional): p-value threshold (Default: 0.05; -log10(0.05)=1.3)
    file (str, optional): save plot to filename
    dir (str, optional): save plot to directory
    color (str, optional): matplotlib color for nonsignificant values
    alpha (float, optional): transparency for nonsignificant values (Default: 0.5)
    edgecol (str, optional): point edge color
    vertical (bool, optional): vertical orientation; otherwise horizontal (Default: True)
    figsize (tuple, optional): figure size
    title (str, optional): plot title
    title_size (int, optional): plot title font size
    title_weight (str, optional): plot title bold, italics, etc.
    title_font (str, optional): plot title font
    x_axis (str, optional): x-axis name
    x_axis_size (int, optional): x-axis name font size
    x_axis_weight (str, optional): x-axis name bold, italics, etc.
    x_axis_font (str, optional): x-axis font
    x_axis_dims (tuple, optional): x-axis dimensions (start, end)
    x_ticks_rot (int, optional): x-axis ticks rotation
    x_axis_font (str, optional): x-axis font
    x_ticks (list, optional): x-axis tick values
    y_axis (str, optional): y-axis name
    y_axis_size (int, optional): y-axis name font size
    y_axis_weight (str, optional): y-axis name bold, italics, etc.
    y_axis_font (str, optional): y-axis font
    y_axis_dims (tuple, optional): y-axis dimensions (start, end)
    y_ticks_rot (int, optional): y-axis ticks rotation
    y_ticks_font (str, optional): y_ticks font
    y_ticks (list, optional): y-axis tick values
    legend_title (str, optional): legend title
    legend_title_size (str, optional): legend title font size
    legend_size (str, optional): legend font size
    legend_bbox_to_anchor (tuple, optional): coordinates for bbox anchor
    legend_loc (str): legend location
    legend_ncol (tuple, optional): # of columns
    display_size (bool, optional): display size on plot (Default: True)
    display_labels (bool, optional): display labels for significant values (Default: True)
    return_df (bool, optional): return dataframe (Default: True)
    show (bool, optional): show plot (Default: True)
    space_capitalize (bool, optional): use re_un_cap() method when applicable (Default: True)
    
    Dependencies: os, matplotlib, seaborn, & pandas
    '''
    # Get dataframe from file path if needed
    if type(df)==str:
        df = io.get(pt=df)

    # Strings with subscripts
    log2 = 'log\u2082'
    log10 = 'log\u2081\u2080'
    
    # Log transform data
    df[f'{log2}({x})'] = [np.log10(xval)/np.log10(2) for xval in df[x]]
    df[f'-{log10}({y})'] = [-np.log10(yval) for yval in df[y]]
    
    # Organize data by significance
    signif = []
    for (log2FC,log10P) in zip(df[f'{log2}({x})'],df[f'-{log10}({y})']):
        if (np.abs(log2FC)>=np.log10(FC_threshold)/np.log10(2))&(log10P>=-np.log10(pval_threshold)): signif.append(f'FC & p-value')
        elif (np.abs(log2FC)<np.log10(FC_threshold)/np.log10(2))&(log10P>=-np.log10(pval_threshold)): signif.append('p-value')
        elif (np.abs(log2FC)>=np.log10(FC_threshold)/np.log10(2))&(log10P<-np.log10(pval_threshold)): signif.append('FC')
        else: signif.append('NS')
    df['Significance']=signif
    #signif_order = ['NS','FC','p-value','FC & p-value']

    # Organize data by abundance
    sizes=(1,100)
    if size_dims is not None: df = df[(df[size]>=size_dims[0])&(df[size]<=size_dims[1])]

    # Set dimensions
    if x_axis_dims==(0,0): x_axis_dims=(min(df[f'{log2}({x})']),max(df[f'{log2}({x})']))
    if y_axis_dims==(0,0): y_axis_dims=(0,max(df[f'-{log10}({y})']))

    # Generate figure
    fig, ax = plt.subplots(figsize=figsize)
    
    if vertical: # orientation
        # with significance boundraries
        plt.vlines(x=-np.log10(FC_threshold)/np.log10(2), ymin=y_axis_dims[0], ymax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
        plt.vlines(x=np.log10(FC_threshold)/np.log10(2), ymin=y_axis_dims[0], ymax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
        plt.hlines(y=-np.log10(pval_threshold), xmin=x_axis_dims[0], xmax=x_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
    
        # with data
        if display_size==False: size=None
        sns.scatterplot(data=df[df['Significance']!='FC & p-value'], x=f'{log2}({x})', y=f'-{log10}({y})', 
                        color=color, alpha=alpha,
                        edgecolor=edgecol, style=stys,
                        size=size, sizes=sizes,
                        ax=ax, **kwargs)
        sns.scatterplot(data=df[(df['Significance']=='FC & p-value')&(df[f'{log2}({x})']<0)], 
                        x=f'{log2}({x})', y=f'-{log10}({y})', 
                        hue=f'{log2}({x})',
                        edgecolor=edgecol, palette='Blues_r', style=stys,
                        size=size, sizes=sizes, legend=False,
                        ax=ax, **kwargs)
        sns.scatterplot(data=df[(df['Significance']=='FC & p-value')&(df[f'{log2}({x})']>0)], 
                        x=f'{log2}({x})', y=f'-{log10}({y})', 
                        hue=f'{log2}({x})',
                        edgecolor=edgecol, palette='Reds', style=stys,
                        size=size, sizes=sizes, legend=False,
                        ax=ax, **kwargs)
        
        # with labels
        if display_labels and label is not None:
            df_signif = df[df['Significance']=='FC & p-value']
            adjust_text([plt.text(y=df_signif.iloc[i][f'{log2}({x})'], 
                                  x=df_signif.iloc[i][f'-{log10}({y})'],
                                  s=l) for i,l in enumerate(df_signif[label])])
        
        # Set x axis
        if x_axis=='': x_axis=f'{log2}({x})'
        plt.xlabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight,fontfamily=x_axis_font)
        if x_ticks==[]: 
            if (x_ticks_rot==0)|(x_ticks_rot==90): plt.xticks(rotation=x_ticks_rot,ha='center',fontfamily=x_ticks_font)
            else: plt.xticks(rotation=x_ticks_rot,ha='right',fontfamily=x_ticks_font)
        else: 
            if (x_ticks_rot==0)|(x_ticks_rot==90): plt.xticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot, ha='center',fontfamily=x_ticks_font)
            else: plt.xticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot,ha='right',fontfamily=x_ticks_font)

        # Set y axis
        if y_axis=='': y_axis=f'-{log10}({y})'
        plt.ylabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight,fontfamily=y_axis_font)

        if y_ticks==[]: plt.yticks(rotation=y_ticks_rot,fontfamily=y_ticks_font)
        else: plt.yticks(ticks=y_ticks,labels=y_ticks,rotation=y_ticks_rot,fontfamily=y_ticks_font)

    else: # Horizontal orientation
        # with significance boundraries
        plt.hlines(y=-np.log10(FC_threshold)/np.log10(2), xmin=y_axis_dims[0], xmax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
        plt.hlines(y=np.log10(FC_threshold)/np.log10(2), xmin=y_axis_dims[0], xmax=y_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)
        plt.vlines(x=-np.log10(pval_threshold), ymin=x_axis_dims[0], ymax=x_axis_dims[1], colors='k', linestyles='dashed', linewidth=1)

        # with data
        if display_size==False: size=None
        sns.scatterplot(data=df[df['Significance']!='FC & p-value'], y=f'{log2}({x})', x=f'-{log10}({y})', 
                        color=color, alpha=alpha,
                        edgecolor=edgecol, style=stys,
                        size=size, sizes=sizes,
                        ax=ax, **kwargs)
        sns.scatterplot(data=df[(df['Significance']=='FC & p-value')&(df[f'{log2}({x})']<0)], 
                        y=f'{log2}({x})', x=f'-{log10}({y})', 
                        hue=f'{log2}({x})',
                        edgecolor=edgecol, palette='Blues_r', style=stys,
                        size=size, sizes=sizes, legend=False,
                        ax=ax, **kwargs)
        sns.scatterplot(data=df[(df['Significance']=='FC & p-value')&(df[f'{log2}({x})']>0)], 
                        y=f'{log2}({x})', x=f'-{log10}({y})', 
                        hue=f'{log2}({x})',
                        edgecolor=edgecol, palette='Reds', style=stys,
                        size=size, sizes=sizes, legend=False,
                        ax=ax, **kwargs)
        
        # with labels
        if display_labels and label is not None:
            df_signif = df[df['Significance']=='FC & p-value']
            adjust_text([plt.text(y=df_signif.iloc[i][f'{log2}({x})'], 
                                  x=df_signif.iloc[i][f'-{log10}({y})'],
                                  s=l) for i,l in enumerate(df_signif[label])])
        
        # Set x axis
        if y_axis=='': y_axis=f'-{log10}({y})'
        plt.xlabel(y_axis, fontsize=y_axis_size, fontweight=y_axis_weight,fontfamily=y_axis_font)
        if y_ticks==[]: 
            if (y_ticks_rot==0)|(y_ticks_rot==90): plt.xticks(rotation=y_ticks_rot,ha='center',fontfamily=y_ticks_font)
            else: plt.xticks(rotation=y_ticks_rot,ha='right',fontfamily=y_ticks_font)
        else: 
            if (y_ticks_rot==0)|(y_ticks_rot==90): plt.xticks(ticks=y_ticks,labels=y_ticks,rotation=y_ticks_rot, ha='center',fontfamily=y_ticks_font)
            else: plt.xticks(ticks=y_ticks,labels=y_ticks,rotation=y_ticks_rot,ha='right',fontfamily=y_ticks_font)

        # Set y axis
        if x_axis=='': x_axis=f'{log2}({x})'
        plt.ylabel(x_axis, fontsize=x_axis_size, fontweight=x_axis_weight,fontfamily=x_axis_font)

        if x_ticks==[]: plt.yticks(rotation=x_ticks_rot,fontfamily=x_ticks_font)
        else: plt.yticks(ticks=x_ticks,labels=x_ticks,rotation=x_ticks_rot,fontfamily=x_ticks_font)

    # Set title
    if title=='' and file is not None: 
        if space_capitalize: title=re_un_cap(".".join(file.split(".")[:-1]))
        else: ".".join(file.split(".")[:-1])
    plt.title(title, fontsize=title_size, fontweight=title_weight, family=title_font)

    # Move legend to the right of the graph
    if legend_items==(0,0): ax.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
                                        bbox_to_anchor=legend_bbox_to_anchor,loc=legend_loc,ncol=legend_ncol)
    else: 
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(title=legend_title,title_fontsize=legend_title_size,fontsize=legend_size,
                  bbox_to_anchor=legend_bbox_to_anchor,loc=legend_loc,ncol=legend_ncol, # Move right of the graph
                  handles=handles[legend_items[0]:legend_items[1]],labels=labels[legend_items[0]:legend_items[1]]) # Only retains specified labels

    # Save & show fig; return dataframe
    if file is not None and dir is not None:
        io.mkdir(dir) # Make output directory if it does not exist
        plt.savefig(fname=os.path.join(dir, file), dpi=600, bbox_inches='tight', format=f'{file.split(".")[-1]}')
    if show: plt.show()
    if return_df: return df

# Color display methods
def matplotlib_cmaps():
    ''' 
    matplotlib_cmaps(): view all matplotlib color maps
    
    Dependencies: matplotlib & numpy
    '''
    # Get the list of all available colormaps
    cmaps = plt.colormaps()

    # Create some data to display the colormaps
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))

    # Define how many colormaps to show per row
    n_col = 4
    n_row = len(cmaps) // n_col + 1

    # Create a figure to display the colormaps
    fig, axes = plt.subplots(n_row, n_col, figsize=(12, 15))
    axes = axes.flatten()

    # Loop through all colormaps and display them
    for i, cmap in enumerate(cmaps):
        ax = axes[i]
        ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(cmap))
        ax.set_title(cmap, fontsize=8)
        ax.axis('off')

    # Turn off any unused subplots
    for j in range(i+1, len(axes)):
        axes[j].axis('off')

    # Display plot
    plt.tight_layout()
    plt.show()

def seaborn_palettes():
    ''' 
    seaborn_palettes(): view all seaborn color palettes
    
    Dependencies: matplotlib & seaborn
    '''
    # List of common Seaborn palettes
    palettes = [
        "deep", "muted", "bright", "pastel", "dark", "colorblind",
        "husl", "hsv", "Paired", "Set1", "Set2", "Set3", "tab10", "tab20"
    ]
    
    # Create a figure to display the color palettes
    n_col = 2  # Palettes per row
    n_row = len(palettes) // n_col + 1
    fig, axes = plt.subplots(n_row, n_col, figsize=(10, 8))
    axes = axes.flatten()

    # Loop through the palettes and display them
    for i, palette in enumerate(palettes):
        ax = axes[i]
        colors = sns.color_palette(palette, 10)  # Get the palette with 10 colors
        # Plot the colors as a series of rectangles (like palplot)
        for j, color in enumerate(colors):
            ax.add_patch(plt.Rectangle((j, 0), 1, 1, color=color))
        
        ax.set_xlim(0, len(colors))
        ax.set_ylim(0, 1)
        ax.set_title(palette, fontsize=12)
        ax.axis('off')

    # Turn off any unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    # Display plot
    plt.tight_layout()
    plt.show()